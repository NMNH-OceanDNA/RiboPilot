// RiboPilot WF1 (self-contained): read samples from the project DB via sqlite3,
// preprocess (fastp) + assemble (GetOrganelle), then score -> stitch ->
// annotate -> reconcile as R-script/tool processes that write results back to
// the DB. No nf-sqldb plugin, no R around Nextflow.

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { query_samples; score; stitch; annotate_tools; reconcile; curate } from './modules/pipeline_r.nf'
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

    pre_in = rows.map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.cpus, r.fastp) }
    preprocess(pre_in)

    asm_opts = rows.map { r -> tuple(r.ID, r.cpus, r.getOrganelle, file(r.seeds_db), file(r.labels_db)) }
    asm_in = preprocess.out.reads
        .filter { id, reads, after -> (after as Integer) >= (params.minDepth as Integer) }
        .map { id, reads, after -> tuple(id, reads) }
        .join(asm_opts)
        .map { id, reads, cpus, go, seeds, labels -> tuple(id, reads, cpus, go, seeds, labels) }
    assemble(asm_in)

    // WF1 finalization (self-contained R/tool processes writing the DB)
    score(assemble.out.paths)
    stitch(score.out)
    // remote GenBank BLAST of the stitched sequence -> ranked reference candidates
    blast_genbank(stitch.out.map { id, fasta -> id })
    // fetch each candidate reference (sequence + rRNA/ITS features) from NCBI
    blast_ref_fetch(blast_genbank.out)
    // map trimmed reads back to the stitched sequence -> per-position coverage
    cov_in = stitch.out.join(preprocess.out.reads.map { id, reads, after -> tuple(id, reads) })
    coverage(cov_in)
    annotate_tools(stitch.out)
    reconcile(annotate_tools.out)

    // per-feature sample-vs-reference alignment, once BOTH the reference is
    // fetched and the sample is annotated (reconcile writes the annotations table)
    ref_align_in = blast_ref_fetch.out.map { id -> tuple(id, true) }
        .join(reconcile.out.map { id -> tuple(id, true) })
        .map { id, a, b -> id }
    blast_ref_align(ref_align_in)

    // curate the reconciled model (structural checks always; reference-identity
    // flags whenever blast_ref_align has already written alignments)
    curate(reconcile.out)
}
