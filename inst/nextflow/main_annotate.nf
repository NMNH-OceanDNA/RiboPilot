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
    rows = query_annotate_samples().splitCsv(header: true, sep: '\t')
    ids = rows.map { r -> r.ID }
    stitch(ids)
    // carry each sample's annotate_opts cpus/memory into annotate_tools
    ann_res = rows.map { r -> tuple(r.ID, r.cpus, r.mem) }
    ann_in = stitch.out.join(ann_res).map { id, fasta, cpus, mem -> tuple(id, fasta, cpus, mem) }
    annotate_tools(ann_in)
    reconcile(annotate_tools.out)
    // per-feature sample-vs-reference alignment before curation, so the
    // reference-identity flags are available (opportunistic: no-op if WF1's
    // GenBank BLAST found no reference for this sample)
    blast_ref_align(reconcile.out)
    curate(blast_ref_align.out)
}
