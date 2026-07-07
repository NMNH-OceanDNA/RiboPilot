// RiboPilot WF2 (self-contained): annotate + curate assemblies the user has
// reviewed and locked (assemble_lock = 1), and re-annotate samples whose
// selected path later changed. Reads queued sample IDs from the DB (sqlite3),
// then stitch the selected path -> barrnap + ITSx -> reconcile -> reference
// alignment -> curate, writing results back to the DB. The reference sequences
// were already fetched in WF1, so blast_ref_align just reads the DB. Same
// sqlite3 + R-script-process pattern as WF1.

include { query_annotate_samples; stitch; annotate_tools; reconcile; curate } from './modules/pipeline_r.nf'
include { blast_ref_align } from './modules/blast_ref.nf'

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
    // per-feature sample-vs-reference alignment before curation, so the
    // reference-identity flags are available (opportunistic: no-op if WF1's
    // GenBank BLAST found no reference for this sample)
    blast_ref_align(reconcile.out)
    curate(blast_ref_align.out)
}
