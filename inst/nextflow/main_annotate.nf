// RiboPilot WF2 (self-contained): annotate + curate the assembly UNITS
// (ID, path, scaffold) of samples the user has reviewed and locked
// (assemble_lock = 1). Each non-ignored scaffold is annotated INDEPENDENTLY
// (no stitching): its sequence is pulled from the assemblies table ->
// barrnap + ITSx -> reconcile -> curate, per unit, written back to the DB.
// (Per-scaffold BLAST + reference alignment are added in a later step.)

include { query_annotate_samples; annotate_tools; reconcile; curate } from './modules/pipeline_r.nf'

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
    curate(reconcile.out)
}
