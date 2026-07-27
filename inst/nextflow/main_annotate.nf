// RiboPilot WF2 (self-contained): RE-annotate + curate assembly UNITS
// (ID, path, scaffold) that are queued for it (annotate_switch IN (1,4)) on
// samples the user has reviewed and locked (assemble_lock = 1).
//
// WF1 now annotates every scaffold as it is assembled, so this is no longer the
// first annotation pass; it is the re-run path for units whose input changed --
// a trimmed assembly, edited annotate_opts, or a scaffold put back after being
// dropped. Each scaffold is annotated INDEPENDENTLY (no stitching): its sequence
// is pulled from the assemblies table -> barrnap + ITSx -> reconcile ->
// reference alignment -> curate, per unit, written back to the DB.

include { query_annotate_samples; annotate_tools; reconcile; curate } from './modules/pipeline_r.nf'
include { blast_ref_align } from './modules/blast_ref.nf'

params.publishDir = 'out'
params.db = ''
params.projectDir = '.'
params.ribopilot_src = ''
params.taxdump_dir = ''

workflow {
    rows = query_annotate_samples().splitCsv(header: true, sep: '\t')
    // one job per (ID, path, scaffold) unit, carrying its compute settings
    units = rows.map { r -> tuple(r.ID, r.path, r.scaffold, r.cpus, r.mem) }
    annotate_tools(units)
    reconcile(annotate_tools.out)
    // per-feature sample-vs-reference alignment (per unit) before curation, so
    // the reference-identity flags are available (no-op if WF1 found no reference)
    blast_ref_align(reconcile.out)
    curate(blast_ref_align.out)
}
