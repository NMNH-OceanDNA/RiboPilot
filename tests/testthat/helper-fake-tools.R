# Stand-ins for the pipeline's external tools, so the real bin scripts and the
# real main.nf run without Docker. Every tool answers from one real rDNA
# scaffold (fixtures/chimera/clean.fa, 6310 bp): 18S 351-2150, ITS1, 5.8S
# 2651-2807, ITS2, 28S 3991-6310.
local_fake_pipeline_tools <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  sq <- file.path(dir, "scaffold.fa")
  writeLines(c(">scaffold", paste(readLines(test_path("fixtures", "chimera", "clean.fa"))[-1],
    collapse = "")), sq)
  tools <- list(
    fastp = c(
      'while [ $# -gt 0 ]; do case "$1" in',
      '  --in1) i1=$2; shift;; --in2) i2=$2; shift;; --out1) o1=$2; shift;;',
      '  --out2) o2=$2; shift;; --json) js=$2; shift;; esac; shift; done',
      'cp "$i1" "$o1"; cp "$i2" "$o2"',
      paste0('echo \'{"summary":{"before_filtering":{"total_reads":2000},"after_filtering":',
        '{"total_reads":2000,"read1_mean_length":150,"read2_mean_length":150}}}\' > "$js"')),
    get_organelle_from_reads.py = c(
      'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && out=$2; shift; done',
      'mkdir -p "$out"',
      sprintf('cp %s "$out/fake.graph1.1.path_sequence.fasta"', shQuote(sq)),
      'echo "fake GetOrganelle" > "$out/get_org.log.txt"'),
    pybarrnap = c(
      'awk -v OFS="\\t" \'/^>/ {id = substr($1, 2)',
      '  print id, "barrnap:0.9", "rRNA", 351, 2150, "0", "+", ".", "Name=18S_rRNA;product=18S ribosomal RNA"',
      '  print id, "barrnap:0.9", "rRNA", 2651, 2807, "0", "+", ".", "Name=5_8S_rRNA;product=5.8S ribosomal RNA"',
      '  print id, "barrnap:0.9", "rRNA", 3991, 6310, "0", "+", ".", "Name=28S_rRNA;product=28S ribosomal RNA"}\' "${@: -1}"'),
    ITSx = c(
      'while [ $# -gt 0 ]; do case "$1" in -i) in=$2; shift;; -o) out=$2; shift;; esac; shift; done',
      paste0('awk -v OFS="\\t" \'/^>/ {print substr($1, 2), "6310 bp.", "SSU: 351-2150", "ITS1: 2151-2650",',
        ' "5.8S: 2651-2807", "ITS2: 2808-3990", "LSU: 3991-6310", ""}\' "$in" > "$out.positions.txt"'),
      'echo "ITSx run finished at $(date)." > "$out.summary.txt"'),
    blastn = c(
      'case " $* " in',
      paste0('  *" -remote "*) printf "FAKE0001.1\\tUpeneus parvus 18S ribosomal RNA gene, partial sequence;',
        ' internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2,',
        ' complete sequence; and 28S ribosomal RNA gene, partial sequence',
        '\\t99.5\\t100\\t0.0\\t1000\\t11000\\t6310\\t1\\t6310\\t1\\t6310\\t6310\\tq1\\n" ;;'),
      '  *) prev=""; for a in "$@"; do [ "$prev" = "-out" ] && : > "$a"; prev=$a; done ;;',
      'esac'),
    efetch = c(
      'while [ $# -gt 0 ]; do [ "$1" = "-id" ] && ids=$2; shift; done',
      sprintf(paste0('for a in ${ids//,/ }; do echo ">$a Upeneus parvus 18S ribosomal RNA gene, partial sequence";',
        ' tail -n +2 %s; done'), shQuote(sq))),
    makeblastdb = "exit 1",
    `bowtie2-build` = "exit 0",
    bowtie2 = "exit 0",
    samtools = c(
      'case "$1" in',
      '  depth) awk -v OFS="\\t" \'/^>/ {id = substr($1, 2); n = 0; next} {for (i = 1; i <= length($0); i++) print id, ++n, 30}\' path.fa ;;',
      '  sort) while [ $# -gt 0 ]; do [ "$1" = "-o" ] && : > "$2"; shift; done; cat > /dev/null ;;',
      '  view) [ "${@: -1}" = "-" ] && cat > /dev/null ;;',
      'esac',
      'exit 0'))
  for (nm in names(tools)) {
    f <- file.path(dir, nm)
    writeLines(c("#!/usr/bin/env bash", tools[[nm]]), f)
    Sys.chmod(f, "0755")
  }
  list(bin = dir, scaffold = sq)
}
