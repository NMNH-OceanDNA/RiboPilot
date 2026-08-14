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

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { fallback_assemble } from './modules/fallback.nf'
include { query_samples; query_annotate_units; query_sample_units; score; annotate_tools; reconcile; curate; triage; screen } from './modules/pipeline_r.nf'
include { blast_genbank } from './modules/blast.nf'
include { blast_ref_fetch; blast_ref_align } from './modules/blast_ref.nf'
include { coverage } from './modules/coverage.nf'

params.publishDir = 'out'
params.minDepth = 0
params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.ref_dir = ''
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

workflow {
    // samples straight from the DB (tab-separated, with header)
    rows = query_samples().splitCsv(header: true, sep: '\t')

    pre_in = rows.map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.pre_cpus, r.pre_mem, r.fastp) }
    preprocess(pre_in)

    // Seeds / labels stay strings here. sqlite3 renders a SQL NULL as '', and
    // assemble_opts.seeds_db is nullable, so file('') would throw for EVERY row
    // -- including fallback rows that never use seeds -- before a task is
    // submitted. Only the getOrganelle arm below turns them into files.
    asm_opts = rows.map { r -> tuple(r.ID, r.asm_cpus, r.asm_mem, r.getOrganelle, live_db(r.seeds_db), live_db(r.labels_db), r.assembler) }
    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(asm_opts)

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
    assemble(asm_branched.go.map { id, reads, cpus, mem, go, seeds, labels, _asm ->
        tuple(id, reads, cpus, mem, go, file(seeds), file(labels)) })
    fallback_assemble(asm_branched.fb.map { id, reads, cpus, mem, _go, _seeds, _labels, asm ->
        // a seedless getorganelle sample lands here; give it a real fallback method
        tuple(id, reads, cpus, mem, asm == 'getorganelle' ? 'bait' : asm) })

    // Load assemblies + score paths (also seeds one annotate unit per scaffold).
    // Both assemblers emit the same fasta contract, so they merge here and
    // everything below is shared.
    score(assemble.out.paths.mix(fallback_assemble.out.paths))
    // per-scaffold remote GenBank BLAST -> ranked reference candidates + screen,
    // then fetch each candidate reference (sequence + features) from NCBI
    blast_genbank(score.out)
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
    // Units queued for annotation only, on samples this run is not assembling.
    // Same tuple shape as the freshly assembled ones, so they mix into one
    // channel and the whole tail below is shared.
    requeued = query_annotate_units()
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.ID, r.path, r.scaffold, r.cpus, r.mem) }

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

    // Propose keep/drop/undecided once ALL of a sample's units are curated, so
    // sibling concordance and the within-sample depth ratio see the full set.
    // groupTuple emits one entry per sample when the upstream channel closes.
    triage(curate.out.map { id, upath, uscaf -> tuple(id, "${upath}.${uscaf}") }
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
