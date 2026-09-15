# Build the rDNA-filtered test read subsets that new_test_project() downloads
# by default (release tag RP_TEST_READS_TAG in R/init_test_project.R).
#
# Usage (repo root):
#   Rscript dev/make_test_reads.R <full_reads_dir> <out_dir> [project_db] [max_pairs]
#
# full_reads_dir  holds <acc>_R1.fastq.gz / _R2 for the mapping_test.csv samples
#                 (a new_test_project(full_size = TRUE) data/ directory).
# project_db      optional .sqlite of a full-size project that has been through
#                 WF1: each sample's own kept scaffolds are added to the bait so
#                 ITS and flank reads are recruited, not just 18S/28S.
# max_pairs       cap per sample (default 80000).
#
# Needs minimap2 and seqkit on PATH (or MINIMAP2 / SEQKIT env vars pointing at
# the binaries). Publish the outputs with:
#   gh release create <tag> out_dir/*.fastq.gz out_dir/*.sha256 --title "..." --notes "..."
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: make_test_reads.R <full_reads_dir> <out_dir> [project_db] [max_pairs]")
reads_dir <- args[1]; out_dir <- args[2]
project_db <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL
max_pairs <- if (length(args) >= 4) as.integer(args[4]) else 80000L
tool <- function(name) {
  p <- Sys.getenv(toupper(name), Sys.which(name))
  if (!nzchar(p)) stop(name, " not found on PATH; set ", toupper(name), "=")
  p
}
minimap2 <- tool("minimap2"); seqkit <- tool("seqkit")
threads <- max(1L, parallel::detectCores() - 2L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
seeds <- "ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta"
map <- read.csv("inst/test_data/mapping_test.csv", stringsAsFactors = FALSE)
map <- map[!(map$Local %in% TRUE), ]

bait_for <- function(acc) {
  fa <- tempfile(fileext = ".fa")
  file.copy(seeds, fa)
  if (!is.null(project_db)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), project_db); on.exit(DBI::dbDisconnect(con))
    s <- DBI::dbGetQuery(con, "SELECT path, scaffold, sequence FROM assemblies WHERE ID = ? AND COALESCE(ignore, 0) = 0 AND sequence IS NOT NULL", list(acc))
    if (nrow(s)) cat(sprintf(">%s_p%d_s%d\n%s\n", acc, s$path, s$scaffold, s$sequence), file = fa, append = TRUE, sep = "")
    message(sprintf("  bait: %d own scaffold(s) + seeds", nrow(s)))
  }
  fa
}

for (i in seq_len(nrow(map))) {
  acc <- map$ID[i]
  r1 <- file.path(reads_dir, map$R1[i]); r2 <- file.path(reads_dir, map$R2[i])
  o1 <- file.path(out_dir, map$R1[i]); o2 <- file.path(out_dir, map$R2[i])
  if (!file.exists(r1) || !file.exists(r2)) { message(sprintf("[%s] full reads missing, skipped", acc)); next }
  if (file.exists(o1) && file.exists(o2)) { message(sprintf("[%s] already built, skipped", acc)); next }
  message(sprintf("[%s] %s", acc, map$Taxon[i]))
  fa <- bait_for(acc); names <- tempfile(fileext = ".txt")
  # PAF lists mapped reads only: pairs with either mate mapped, capped in file
  # order (shotgun order is random). minimap2 appends /1 /2 when both mates
  # share a name, which SRA-style headers do; strip it or nothing matches.
  cmd <- sprintf("%s -t %d -x sr %s %s %s 2>/dev/null | cut -f1 | sed 's#/[12]$##' | sort -u | head -n %d > %s",
    shQuote(minimap2), threads, shQuote(fa), shQuote(r1), shQuote(r2), max_pairs, shQuote(names))
  system(cmd)  # head closes the pipe early, so the exit status is not informative
  n <- length(readLines(names))
  if (n == 0) stop("mapping produced no reads for ", acc)
  for (io in list(c(r1, o1), c(r2, o2))) {
    if (system2(seqkit, c("grep", "-f", shQuote(names), shQuote(io[1]), "-o", shQuote(io[2])), stderr = FALSE) != 0) stop("seqkit failed for ", acc)
    system(sprintf("cd %s && sha256sum %s > %s.sha256", shQuote(out_dir), shQuote(basename(io[2])), shQuote(basename(io[2]))))
  }
  message(sprintf("  kept %d pairs -> %s (%.1f MB)", n, basename(o1), file.size(o1) / 1e6))
}
