// RiboPilot assembly-fallback workflow: for samples queued with assemble_switch
// = 5 (GetOrganelle failed the QC gate, or a manual cross-check), preprocess +
// reference-guided baiting/de-novo assembly (Option B2 / A), then the standard
// WF1 tail (score -> stitch -> annotate -> reconcile -> curate). Self-contained
// DB I/O, same pattern as main.nf.

include { preprocess } from './modules/preprocess.nf'
include { fallback_assemble } from './modules/fallback.nf'
include { query_fallback_samples; score; stitch; annotate_tools; reconcile; curate } from './modules/pipeline_r.nf'
include { blast_genbank } from './modules/blast.nf'
include { blast_ref_fetch; blast_ref_align } from './modules/blast_ref.nf'
include { coverage } from './modules/coverage.nf'

params.publishDir = 'out'
params.db = ''
params.projectDir = '.'
params.rawDir = ''
params.ribopilot_src = ''
params.taxdump_dir = ''
params.ncbi_api_key = ''
params.method = 'bait'

workflow {
    rows = query_fallback_samples().splitCsv(header: true, sep: '\t')

    pre_in = rows.map { r -> tuple(r.ID, file(params.rawDir + '/' + r.R1), file(params.rawDir + '/' + r.R2), r.cpus, r.fastp) }
    preprocess(pre_in)

    // reference-guided baiting / de-novo assembly on the trimmed reads
    fb_opts = rows.map { r -> tuple(r.ID, r.cpus) }
    fb_in = preprocess.out.reads
        .map { id, reads, after -> tuple(id, reads) }
        .join(fb_opts)
        .map { id, reads, cpus -> tuple(id, reads, cpus, params.method) }
    fallback_assemble(fb_in)

    // WF1 finalization (identical to main.nf)
    score(fallback_assemble.out.paths)
    stitch(score.out)
    blast_genbank(stitch.out.map { id, fasta -> id })
    blast_ref_fetch(blast_genbank.out)
    cov_in = stitch.out.join(preprocess.out.reads.map { id, reads, after -> tuple(id, reads) })
    coverage(cov_in)
    annotate_tools(stitch.out)
    reconcile(annotate_tools.out)
    ref_align_in = blast_ref_fetch.out.map { id -> tuple(id, true) }
        .join(reconcile.out.map { id -> tuple(id, true) })
        .map { id, a, b -> id }
    blast_ref_align(ref_align_in)
    curate(reconcile.out)
}
