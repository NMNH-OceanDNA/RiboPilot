// RiboPilot recheck: reruns every offline check on the review copies queued on
// Annotation Review (dev/specs/2026-10-08-contig-recheck-design.md, section 7).
// Each copy goes through the same rRNA tool script as the main pipeline's
// annotate_tools (bin/rp_rrna_tools.sh), then rp_recheck.R scores it with the
// main pipeline's check functions in review mode. No reads, no NCBI.
//
// Launched from <project>/.recheck with its own work folder and no -resume
// (nf_recheck_command() in R/nf_command.R), so the main pipeline's run history
// is never touched.

include { query_recheck; recheck_tools; recheck } from './modules/recheck.nf'

params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.blastdb = ''

workflow {
    units = query_recheck().splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.ID, r.path, r.scaffold, r.rev, r.cpus, r.mem, r.mode) }
    recheck_tools(units)
    recheck(recheck_tools.out)
}
