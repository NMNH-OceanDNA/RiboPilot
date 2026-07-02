// WF1 coverage: map the trimmed reads back to the stitched nrDNA sequence
// (bowtie2) and record per-position read depth (samtools depth -a) so the app
// can draw a coverage track. Best-effort: 'ignore' on failure keeps the run
// going (coverage is informational, not required for annotation).

process coverage {
    cache false
    tag "${id}"
    cpus 4
    errorStrategy 'ignore'

    input:
    tuple val(id), path(stitched), path(reads)

    output:
    val id

    script:
    """
    bowtie2-build --quiet ${stitched} ref
    bowtie2 -x ref -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --no-unal 2> bowtie2.log \
        | samtools sort -@ ${task.cpus} -o aln.bam -
    samtools depth -a aln.bam > depth.txt
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_coverage.R '${params.db}' '${id}' depth.txt
    """
}
