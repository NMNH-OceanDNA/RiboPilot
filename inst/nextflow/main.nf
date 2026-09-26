// RiboPilot WF1 (self-contained): read samples from the project DB via sqlite3,
// preprocess (fastp) + assemble, score the candidate paths, then annotate and
// curate EVERY assembly unit. Assemblies are NOT stitched; each scaffold is
// annotated independently.
//
// The assembler is chosen PER SAMPLE by assemble_opts.assembler: GetOrganelle
// by default, or the baiting/de-novo fallback for samples GetOrganelle could
// not assemble. Both emit the same fastas and both feed the same tail, so there
// is one workflow rather than two. There used to be two, and the copy drifted
// out of step with this file until it no longer ran; see
// dev/design/14_ASSEMBLER_SWAP.md.
//
// Annotation runs with assembly because it is cheap next to
// assembly (a few seconds per scaffold against GetOrganelle's minutes to hours)
// and because the reviewer needs the feature model to decide which scaffolds to
// keep: a rRNA-features column that is empty at the moment of the keep/drop
// decision is worse than useless, since "nothing found" and "not looked at yet"
// render identically.
//
// Re-annotation lives here too. There is no second workflow: units queued after
// review (a trim, an options change, a restored scaffold, an edited sequence)
// enter through query_annotate_units and skip straight to the annotation chain,
// paying for neither assembly nor a second remote BLAST. See
// dev/design/15_WF2_DECOMMISSION.md.

include { preprocess; below_min_depth } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { fallback_assemble } from './modules/fallback.nf'
include { find_rrna; stage_assembly } from './modules/find_rrna.nf'
include { query_samples; query_annotate_units; query_sample_units; score; annotate_tools; reconcile; curate; triage; screen; chimera } from './modules/pipeline_r.nf'
include { blast_genbank } from './modules/blast.nf'
include { blast_ref_fetch; blast_ref_align } from './modules/blast_ref.nf'
include { coverage } from './modules/coverage.nf'

params.publishDir = 'out'
params.minDepth = 0
params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.asmbDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.ref_dir = ''
params.blastdb = ''
params.ncbi_api_key = ''

// A seeds/labels path stored in the database can point into a reference folder
// that has since been deleted. file() on a dead path stages a DANGLING SYMLINK
// rather than failing, and GetOrganelle exits 0 after "seeds not found", so the
// sample would report success having assembled nothing. Blanking the dead path
// sends the sample down the fallback arm instead. A URL is left alone: Nextflow
// fetches it at stage time, and probing it here would be a network call per row.
def live_db(p) {
    (p && (p ==~ '^https?://.*' || file(p).exists())) ? p : ''
}

// The fallback downloads a URL itself, so a dead one fails only its sample.
def is_url(p) { p ==~ '^https?://.*' }
def stage_db(p) { p && !is_url(p) ? file(p) : [] }
def url_db(p) { is_url(p) ? p : '' }

workflow {
    // samples straight from the DB (tab-separated, with header)
    rows = query_samples().splitCsv(header: true, sep: '\t')

    // A sample supplying its own assembly skips preprocess + assembly: its rDNA
    // contigs are found inside what it supplied. It still preprocesses when it
    // named reads, because those reads are what coverage maps back.
    sample_rows = rows.branch { r ->
        user: r.Assembly?.trim() && r.find_rrna?.toString() != '0'
        asis: r.Assembly?.trim()
        reads: true
    }

    pre_in = rows
        .filter { r -> r.R1?.trim() && r.R2?.trim() }
        .map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.pre_cpus, r.pre_mem, r.fastp) }
    preprocess(pre_in)

    // Supplied assemblies: scan for rDNA contigs, or take the file as given when
    // the sample has the search switched off (an ITS-only fragment carries no
    // barrnap anchor and would otherwise be dropped).
    find_rrna(sample_rows.user.map { r ->
        tuple(r.ID, file(params.asmbDir + '/' + r.Assembly), r.asm_cpus, r.asm_mem, r.rrna_tool, r.rrna_kingdom) })
    stage_assembly(sample_rows.asis.map { r -> tuple(r.ID, file(params.asmbDir + '/' + r.Assembly)) })

    // Seeds / labels stay strings here. sqlite3 renders a SQL NULL as '', and
    // assemble_opts.seeds_db is nullable, so file('') would throw for EVERY row
    // -- including fallback rows that never use seeds -- before a task is
    // submitted. The assembler arms below turn them into files, or [] (nothing staged) when blank.
    asm_opts = rows.map { r -> tuple(r.ID, r.asm_cpus, r.asm_mem, r.getOrganelle, live_db(r.seeds_db), live_db(r.labels_db), r.assembler, r.bowtie2, r.spades) }
    // Only samples assembling from reads: a supplied assembly may also name
    // reads (for coverage), and those must not be assembled as well.
    read_ids = sample_rows.reads.map { r -> tuple(r.ID, true) }
    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(read_ids)
        .map { id, reads, _flag -> tuple(id, reads) }
        .join(asm_opts)
    below_min_depth(preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) < (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, after) }
        .join(read_ids)
        .map { id, after, _flag -> tuple(id, after) })

    // Per-sample assembler. The fallback arm is `true` rather than an equality
    // test on purpose: branch DROPS an item matching no condition, silently, so
    // an unrecognised assembler value would vanish from the run instead of
    // assembling. Defaulting the unknown case to the fallback keeps the sample.
    // A sample with no seeds/labels database cannot run GetOrganelle, so it
    // takes the same fallback arm rather than being dropped. live_db() above
    // blanks a path whose file is gone, so those land here too.
    asm_branched = asm_in.branch {
        go: it[7] == 'getorganelle' && it[5] && it[6]
        fb: true
    }
    assemble(asm_branched.go.map { id, reads, cpus, mem, go, seeds, labels, _asm, _bt2, _spades ->
        tuple(id, reads, cpus, mem, go, file(seeds), file(labels)) })
    fallback_assemble(asm_branched.fb.map { id, reads, cpus, mem, _go, seeds, labels, asm, bt2, spades ->
        // a seedless getorganelle sample lands here; give it a real fallback method
        def method = asm == 'getorganelle' ? 'bait' : asm
        def s = method == 'bait' ? seeds : ''
        tuple(id, reads, cpus, mem, method, stage_db(s), stage_db(labels), url_db(s), url_db(labels), bt2, spades) })

    // Load assemblies + score paths (also seeds one annotate unit per scaffold).
    // Both assemblers emit the same fasta contract, so they merge here and
    // everything below is shared.
    score(assemble.out.paths
        .mix(fallback_assemble.out.paths)
        .mix(find_rrna.out.paths)
        .mix(stage_assembly.out.paths))
    // Units queued for annotation only, on samples this run is not assembling.
    // Same tuple shape as the freshly assembled ones, so they mix into one
    // channel and the whole tail below is shared.
    requeued_rows = query_annotate_units().splitCsv(header: true, sep: '\t')
    requeued = requeued_rows.map { r -> tuple(r.ID, r.path, r.scaffold, r.cpus, r.mem) }
    // a scaffold put back after a load-time drop was never searched or screened
    unsearched = requeued_rows.filter { r -> r.searched == '0' }.map { r -> r.ID }.unique()
    // per-scaffold remote GenBank BLAST -> ranked reference candidates + screen,
    // then fetch each candidate reference (sequence + features) from NCBI
    blast_genbank(score.out.mix(unsearched))
    blast_ref_fetch(blast_genbank.out)
    // per-unit coverage (independent branch off score)
    cov_in = score.out.join(preprocess.out.reads.map { id, reads, after -> tuple(id, reads) })
    coverage(cov_in)

    // Annotation fans out per unit straight off score, IN PARALLEL with the
    // remote BLAST. Deliberately not chained after blast_ref_fetch: blast_genbank
    // carries errorStrategy 'ignore', so a sample whose remote search fails
    // emits nothing, and chaining would leave that sample with no annotation at
    // all -- worse than the behaviour being replaced.
    query_sample_units(score.out)
    units = query_sample_units.out
        .flatMap { _id, tsv ->
            tsv.splitCsv(header: true, sep: '\t')
               .collect { r -> tuple(r.ID, r.path, r.scaffold, r.cpus, r.mem) } }

    annotate_tools(units.mix(requeued))
    reconcile(annotate_tools.out)

    // Reference alignment and curation read blast_ref_alignment / blast_screen,
    // so they must follow the reference chain. collect() is a barrier rather than
    // a per-sample join because reconcile emits MANY units per sample while
    // blast_ref_fetch emits one -- a join would silently drop the extras, and a
    // sample whose BLAST failed would drop entirely. ifEmpty keeps the branch
    // alive when no sample produced a reference (offline run, BLAST disabled).
    // The collected list is reduced to one scalar sentinel before combine:
    // combine splices a list operand element-by-element, so carrying the list
    // itself would widen every unit tuple by however many samples were fetched.
    ref_done = blast_ref_fetch.out.collect().map { 'ready' }.ifEmpty('ready')
    aln_in = reconcile.out.combine(ref_done)
        .map { id, upath, uscaf, _ready -> tuple(id, upath, uscaf) }
    blast_ref_align(aln_in)
    curate(blast_ref_align.out)

    // The screen reads only blast_hits, so it hangs off blast_genbank rather
    // than the curate chain. triage joins on it so a sample is triaged once its
    // curation is complete and its screen has been ATTEMPTED. remainder: true is
    // load-bearing: screen carries errorStrategy 'ignore', so an inner join
    // would drop a fully curated sample whenever its screen or its upstream
    // BLAST failed, which is the same trap the ref_done sentinel above avoids.
    screen(blast_genbank.out)

    // Chimera detection needs annotation (reconcile), and reads depth/clips and
    // blast_hits from the db, so it waits on coverage and blast_genbank the same
    // way reference alignment waits on blast_ref_fetch: a scalar barrier, not a
    // per-unit join.
    chim_done = coverage.out.collect().map { 'ready' }.ifEmpty('ready')
        .combine(blast_genbank.out.collect().map { 'ready' }.ifEmpty('ready'))
        .map { 'ready' }
    chimera(reconcile.out.combine(chim_done).map { id, upath, uscaf, _r -> tuple(id, upath, uscaf) })

    // Propose keep/drop/undecided once ALL of a sample's units are curated AND
    // chimera-checked, so sibling concordance and the within-sample depth ratio
    // see the full set. remainder: true keeps a unit whose chimera step was
    // ignored (chimera carries errorStrategy 'ignore').
    unit_done = curate.out.map { id, upath, uscaf -> tuple("${id}.${upath}.${uscaf}", id) }
        .join(chimera.out.map { id, upath, uscaf -> tuple("${id}.${upath}.${uscaf}", true) }, remainder: true)
        .map { _k, id, _c -> id }
    triage(unit_done.map { id -> tuple(id, true) }
        .groupTuple()
        .map { id, _units -> tuple(id, true) }
        .join(screen.out.map { id -> tuple(id, true) }, remainder: true)
        // remainder: true also admits screen-only entries, so a sample with a
        // screen but no curation would reach triage. Requiring the left side
        // keeps triage waiting for curation. Both lines are load-bearing:
        // dropping either one reintroduces a bug the other does not cover.
        .filter { entry -> entry[1] != null }
        .map { id, _a, _b -> id })
}
