// RiboPilot WF1 (self-contained): read samples from the project DB via sqlite3,
// preprocess (fastp) + assemble (GetOrganelle), score the candidate paths, then
// annotate and curate EVERY assembly unit. Assemblies are NOT stitched; each
// scaffold is annotated independently.
//
// Annotation runs here rather than only in WF2 because it is cheap next to
// assembly (a few seconds per scaffold against GetOrganelle's minutes to hours)
// and because the reviewer needs the feature model to decide which scaffolds to
// keep: a rRNA-features column that is empty at the moment of the keep/drop
// decision is worse than useless, since "nothing found" and "not looked at yet"
// render identically. WF2 (main_annotate.nf) remains the RE-annotation entry
// point for trims, annotate_opts changes and un-ignored units.

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { query_samples; query_sample_units; score; annotate_tools; reconcile; curate; triage } from './modules/pipeline_r.nf'
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
params.ncbi_api_key = ''

workflow {
    // samples straight from the DB (tab-separated, with header)
    rows = query_samples().splitCsv(header: true, sep: '\t')

    pre_in = rows.map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.pre_cpus, r.pre_mem, r.fastp) }
    preprocess(pre_in)

    asm_opts = rows.map { r -> tuple(r.ID, r.asm_cpus, r.asm_mem, r.getOrganelle, file(r.seeds_db), file(r.labels_db)) }
    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(asm_opts)
        .map { id, reads, cpus, mem, go, seeds, labels -> tuple(id, reads, cpus, mem, go, seeds, labels) }
    assemble(asm_in)

    // Load assemblies + score paths (also seeds one annotate unit per scaffold).
    score(assemble.out.paths)
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
    annotate_tools(units)
    reconcile(annotate_tools.out)

    // Reference alignment and curation read blast_ref_alignment / blast_screen,
    // so they must follow the reference chain. collect() is a barrier rather than
    // a per-sample join because reconcile emits MANY units per sample while
    // blast_ref_fetch emits one -- a join would silently drop the extras, and a
    // sample whose BLAST failed would drop entirely. ifEmpty keeps the branch
    // alive when no sample produced a reference (offline run, BLAST disabled).
    ref_done = blast_ref_fetch.out.collect().ifEmpty { ['none'] }
    aln_in = reconcile.out.combine(ref_done.map { refs -> tuple(refs) })
        .map { id, upath, uscaf, _refs -> tuple(id, upath, uscaf) }
    blast_ref_align(aln_in)
    curate(blast_ref_align.out)

    // Propose keep/drop/undecided once ALL of a sample's units are curated, so
    // sibling concordance and the within-sample depth ratio see the full set.
    // groupTuple emits one entry per sample when the upstream channel closes.
    triage(curate.out.map { id, upath, uscaf -> tuple(id, "${upath}.${uscaf}") }
        .groupTuple()
        .map { id, _units -> id })
}
