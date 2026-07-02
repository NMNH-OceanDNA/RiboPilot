// RiboPilot WF2 (self-contained): re-annotate samples whose selected path
// changed. Reads queued sample IDs from the DB (sqlite3), then stitch the
// selected path -> barrnap + ITSx -> reconcile, writing results back to the DB.
// Same sqlite3 + R-script-process pattern as WF1.

include { query_annotate_samples; stitch; annotate_tools; reconcile } from './modules/pipeline_r.nf'

params.publishDir = 'out'
params.db = ''
params.projectDir = '.'
params.ribopilot_src = ''
params.taxdump_dir = ''

workflow {
    ids = query_annotate_samples().splitText().map { it.trim() }.filter { it }
    stitch(ids)
    annotate_tools(stitch.out)
    reconcile(annotate_tools.out)
}
