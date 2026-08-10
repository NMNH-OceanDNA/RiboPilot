// WF1 coverage (per scaffold): map the trimmed reads back to ALL of the sample's
// non-ignored scaffolds at once (bowtie2, refs named path.scaffold) and record
// per-position read depth, split into one coverage row per assembly unit. The
// scaffold sequences are pulled from the DB so no stitched fasta is needed.
// Best-effort: 'ignore' on failure keeps the run going (coverage is informational).

process coverage {
    cache false
    tag "${id}"
    cpus 4
    memory 8.GB
    errorStrategy 'ignore'

    input:
    tuple val(id), path(reads)

    output:
    val id

    script:
    """
    # build a fasta of this sample's non-ignored scaffolds (headers = path.scaffold)
    sqlite3 -cmd '.timeout 10000' '${params.db}' \
      "SELECT '>'||path||'.'||scaffold||char(10)||sequence FROM assemblies
       WHERE ID='${id}' AND COALESCE(ignore,0)=0;" > scaffolds.fa

    if [ -s scaffolds.fa ]; then
        bowtie2-build --quiet scaffolds.fa ref
        # -q 10 drops multi-mapping reads: the conserved 18S/28S core maps equally
        # well to every scaffold and bowtie2 breaks those ties at random, so
        # unfiltered depths are not comparable between a sample's scaffolds.
        bowtie2 -x ref -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --no-unal 2> bowtie2.log \
            | samtools view -b -q 10 - \
            | samtools sort -@ ${task.cpus} -o aln.bam -
        samtools depth -a aln.bam > depth.txt
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_coverage.R '${params.db}' '${id}' depth.txt
    fi
    """
}
