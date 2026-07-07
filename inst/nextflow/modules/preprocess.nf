process preprocess {

    tag "${id}"
    cpus { nthreads as int }
    memory { "${mem} GB" }
    publishDir { "${launchDir}/${params.publishDir}/${id}" }, overwrite: true, pattern: "*_preprocess.*", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(R1), path(R2), val(nthreads), val(mem), val(fastp_opts)

    output:
    // id, trimmed reads, trimmed read count (for the assembly depth gate)
    tuple val(id), path("${id}_preprocess_R{1,2}.fastq.gz"), env('after'), emit: reads
    path "${id}_preprocess.json"

    script:
    """
    fastp --in1 ${R1} --in2 ${R2} \\
        --out1 ${id}_preprocess_R1.fastq.gz --out2 ${id}_preprocess_R2.fastq.gz \\
        --json ${id}_preprocess.json --thread ${task.cpus} ${fastp_opts}
    after=\$(jq '.summary.after_filtering.total_reads' ${id}_preprocess.json)
    """
}
