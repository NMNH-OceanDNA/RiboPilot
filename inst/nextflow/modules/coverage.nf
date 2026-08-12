// WF1 coverage (per scaffold): map the trimmed reads back to the sample's
// non-ignored scaffolds ONE PATH AT A TIME (refs named path.scaffold) and record
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
    # One index PER PATH, not one index holding the whole sample. Paths are
    # competing whole-repeat candidates: alternative resolutions of the SAME
    # molecule, so they are near-identical by construction and a read maps
    # equally well to all of them. Indexed together, the MAPQ filter below then
    # threw away essentially the entire sample -- measured on a 10-path sample,
    # 1.65M read pairs came out as 2 surviving alignments, and 9 of the 10
    # scaffolds got no coverage row at all. Within one path the scaffolds are
    # genuinely distinct contigs, which is the case the filter was written for.
    PATHS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' \
      "SELECT DISTINCT path FROM assemblies
       WHERE ID='${id}' AND COALESCE(ignore,0)=0 ORDER BY path;")

    : > depth.txt
    : > bowtie2.log
    for P in \$PATHS; do
        sqlite3 -cmd '.timeout 10000' '${params.db}' \
          "SELECT '>'||path||'.'||scaffold||char(10)||sequence FROM assemblies
           WHERE ID='${id}' AND path=\$P AND COALESCE(ignore,0)=0;" > path.fa
        [ -s path.fa ] || continue

        bowtie2-build --quiet path.fa ref
        # -q 10 drops multi-mapping reads: the conserved 18S/28S core maps equally
        # well to every scaffold in the path and bowtie2 breaks those ties at
        # random, so unfiltered depths are not comparable between them.
        bowtie2 -x ref -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --no-unal 2> path.bt2.log \
            | samtools view -b -q 10 - \
            | samtools sort -@ ${task.cpus} -o aln.bam -
        # -aa, not -a: -a omits reference sequences that have no alignments at
        # all, so a scaffold the MAPQ filter emptied disappeared from depth.txt
        # and write_coverage_units() wrote no row for it. "Covered by nothing"
        # and "never measured" then rendered identically in the app -- blank.
        samtools depth -aa aln.bam >> depth.txt
        # keep one readable log across the loop, not one file per path
        { echo "=== path \$P ==="; cat path.bt2.log; } >> bowtie2.log
    done

    if [ -s depth.txt ]; then
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_coverage.R '${params.db}' '${id}' depth.txt
    fi
    """
}
