// Assembly fallback (Option B2 baiting / Option A de-novo): a drop-in
// replacement for the GetOrganelle `assemble` process. rp_fallback.R baits +
// assembles (SPAdes) + keeps rDNA contigs, emitting the same
// <id>_assembly_<n>.fasta header format so the WF1 tail runs unchanged.
process fallback_assemble {

    tag "${id}"
    cpus { nthreads as int }
    // scales on retry; see the OOM retry policy in nextflow.config
    memory { "${(mem as int) * task.attempt} GB" }
    publishDir { "${launchDir}/${params.publishDir}/${id}/assemble" }, overwrite: true, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(reads), val(nthreads), val(mem), val(method)

    output:
    tuple val(id), path("${id}_assembly_*.fasta"), emit: paths

    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_fallback.R \\
        '${params.db}' '${id}' ${reads[0]} ${reads[1]} '${method}' ${task.cpus}
    """
}
