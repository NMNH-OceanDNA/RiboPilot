process assemble {

    tag "${id}"
    cpus { nthreads as int }
    // scales on retry; see the OOM retry policy in nextflow.config
    memory { "${(mem as int) * task.attempt} GB" }
    publishDir { "${launchDir}/${params.publishDir}/${id}/assemble" }, overwrite: true, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(reads), val(nthreads), val(mem), val(getOrganelle_opts), path(seeds), path(labels)

    output:
    // per-path formatted fastas (headers ">id.path.scaffold topology") + graph
    tuple val(id), path("{${id}_assembly_*.fasta,*.gfa}"), emit: paths
    path "get_org.log.txt", optional: true

    script:
    """
    get_organelle_from_reads.py \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o go_out --overwrite \\
        -s ${seeds} --genes ${labels} \\
        -t ${task.cpus} \\
        ${getOrganelle_opts}

    cp go_out/get_org.log.txt . 2>/dev/null || true
    cp go_out/*.gfa . 2>/dev/null || true

    # order + tag each path fasta -> <id>_assembly_<n>.fasta (bin/ helper)
    format_assembly.sh ${id} go_out
    """
}
