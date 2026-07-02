process annotate {

    tag "${id}"
    cpus { nthreads as int }
    publishDir { "${launchDir}/${params.publishDir}/${id}/annotate" }, overwrite: true, mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(id), path(stitched), val(nthreads)

    output:
    tuple val(id), path("${id}.barrnap.gff"), path("${id}.itsx.positions.txt")

    script:
    """
    # rRNA gene anchors (Rfam eukaryote models under the 'fun' kingdom)
    barrnap --kingdom fun --quiet ${stitched} > ${id}.barrnap.gff

    # ITS spacer boundaries
    ITSx -i ${stitched} -o ${id}.itsx --save_regions none --graphical F \\
        --cpu ${task.cpus} -t all --preserve T || true
    [ -f ${id}.itsx.positions.txt ] || : > ${id}.itsx.positions.txt
    """
}
