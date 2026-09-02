# Simulate paired reads from the chimera fixtures into inst/test_data.
# Run once from the repo root: Rscript dev/make_chimera_reads.R
set.seed(20260901)
sim_pairs <- function(fa, id, out_dir, depth = 60, rl = 150L, ins = 350, sd = 30, err = 0.01) {
  seq <- paste(readLines(fa)[-1], collapse = "")
  n <- as.integer(depth * nchar(seq) / (2 * rl))
  rc <- function(s) chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
  mut <- function(s) {
    b <- strsplit(s, "")[[1]]; k <- which(runif(length(b)) < err)
    if (length(k)) b[k] <- sapply(b[k], function(x) sample(setdiff(c("A", "C", "G", "T"), x), 1))
    paste(b, collapse = "")
  }
  r1 <- character(n); r2 <- character(n)
  for (i in seq_len(n)) {
    len <- max(rl + 10L, as.integer(rnorm(1, ins, sd)))
    st <- sample.int(nchar(seq) - len + 1L, 1)
    frag <- substr(seq, st, st + len - 1L)
    r1[i] <- mut(substr(frag, 1, rl)); r2[i] <- mut(substr(rc(frag), 1, rl))
  }
  q <- strrep("?", rl)
  w <- function(reads, mate) {
    con <- gzfile(file.path(out_dir, sprintf("%s_R%d.fastq.gz", id, mate)), "w")
    writeLines(as.vector(rbind(sprintf("@%s_%d/%d", id, seq_len(n), mate), reads, "+", q)), con)
    close(con)
  }
  w(r1, 1); w(r2, 2); n
}
fx <- "tests/testthat/fixtures/chimera"
cat("CHIM00 pairs:", sim_pairs(file.path(fx, "clean.fa"), "CHIM00", "inst/test_data"), "\n")
cat("CHIM01 pairs:", sim_pairs(file.path(fx, "host_fungal.fa"), "CHIM01", "inst/test_data"), "\n")
