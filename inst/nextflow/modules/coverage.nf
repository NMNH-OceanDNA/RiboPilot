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
    : > clips.txt
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
        # --local: a read that runs past a join is soft-clipped, not dropped,
        # which is what the clip count below measures.
        bowtie2 --local -x ref -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --no-unal 2> path.bt2.log \
            | samtools view -b -q 10 - \
            | samtools sort -@ ${task.cpus} -o aln.bam -
        # -aa, not -a: -a omits reference sequences that have no alignments at
        # all, so a scaffold the MAPQ filter emptied disappeared from depth.txt
        # and write_coverage_units() wrote no row for it. "Covered by nothing"
        # and "never measured" then rendered identically in the app -- blank.
        samtools depth -aa aln.bam >> depth.txt

        # per-position count of reads soft-clipped >= 20 bp at that end. Plain
        # awk (mawk in the image), not gawk: no three-arg match(), CIGAR ops
        # are pulled out with substr()/RSTART/RLENGTH instead.
        samtools view aln.bam | awk -v OFS='\t' '{
            cig = \$6; pos = \$4; ref = \$3
            if (match(cig, /^[0-9]+S/)) {
                len = substr(cig, RSTART, RLENGTH - 1) + 0
                if (len >= 20) print ref, pos
            }
            if (match(cig, /[0-9]+S\$/)) {
                len = substr(cig, RSTART, RLENGTH - 1) + 0
                if (len >= 20) {
                    n = 0; s = cig
                    while (match(s, /^[0-9]+[MDN=X]/)) {
                        n += substr(s, 1, RLENGTH - 1) + 0
                        s = substr(s, RLENGTH + 1)
                    }
                    while (match(s, /[0-9]+[MDN=X]/)) {
                        n += substr(s, RSTART, RLENGTH - 1) + 0
                        s = substr(s, RSTART + RLENGTH)
                    }
                    print ref, pos + n - 1
                }
            }
        }' | sort | uniq -c | awk -v OFS='\t' '{print \$2, \$3, \$1}' >> clips.txt
        # keep one readable log across the loop, not one file per path
        { echo "=== path \$P ==="; cat path.bt2.log; } >> bowtie2.log
    done

    if [ -s depth.txt ]; then
        RIBOPILOT_SRC='${params.ribopilot_src}' rp_coverage.R '${params.db}' '${id}' depth.txt clips.txt
    fi
    """
}
