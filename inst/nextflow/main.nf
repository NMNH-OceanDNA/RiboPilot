// RiboPilot WF1 (self-contained): read samples from the project DB via sqlite3,
// preprocess (fastp) + assemble (GetOrganelle), then score -> stitch ->
// annotate -> reconcile as R-script/tool processes that write results back to
// the DB. No nf-sqldb plugin, no R around Nextflow.

include { preprocess } from './modules/preprocess.nf'
include { assemble } from './modules/assemble.nf'
include { query_samples; score; stitch; annotate_tools; reconcile } from './modules/pipeline_r.nf'

params.publishDir = 'out'
params.minDepth = 0
params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''

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
    annotate_tools(stitch.out)
    reconcile(annotate_tools.out)
}
