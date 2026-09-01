# Selecting nuclear rDNA contigs out of a user-supplied assembly.
# See dev/design/20_USER_ASSEMBLIES.md.

gff_line <- function(seqid, name, start, stop) {
  paste(seqid, "barrnap", "rRNA", start, stop, "1e-100", "+", ".",
    paste0("Name=", name, ";product=", name), sep = "\t")
}

test_that("a contig with an 18S anchor is kept and one without is not", {
  gff <- c(
    gff_line("ctg1", "18S_rRNA", 1, 1800),
    gff_line("ctg1", "28S_rRNA", 3000, 6500),
    gff_line("ctg2", "5_8S_rRNA", 100, 250)
  )
  res <- select_rrna_contigs(gff, lengths = c(ctg1 = 12000, ctg2 = 5000))

  expect_equal(res$candidates, "ctg1")
  expect_equal(nrow(res$evidence), 2L)
  expect_equal(res$evidence$selected[res$evidence$contig == "ctg2"], 0L)
  expect_match(res$evidence$reason[res$evidence$contig == "ctg2"], "no 18S or 28S")
  expect_equal(res$evidence$length[res$evidence$contig == "ctg1"], 12000L)
  expect_match(res$evidence$rrna_features[res$evidence$contig == "ctg1"], "18S,28S")
})

test_that("a single supplied rDNA contig passes the scan", {
  res <- select_rrna_contigs(gff_line("rDNA_unit", "28S_rRNA", 500, 4000))
  expect_equal(res$candidates, "rDNA_unit")
  expect_true(is.na(res$evidence$length))
})

test_that("no rRNA hits at all yields no candidates and no evidence", {
  res <- select_rrna_contigs(character(0))
  expect_length(res$candidates, 0L)
  expect_equal(nrow(res$evidence), 0L)

  res <- select_rrna_contigs("# nothing here")
  expect_length(res$candidates, 0L)
  expect_equal(nrow(res$evidence), 0L)
})

test_that("contigs rank by rRNA content and the cap keeps the best", {
  gff <- c(
    gff_line("small", "18S_rRNA", 1, 500),
    gff_line("big", "18S_rRNA", 1, 1800),
    gff_line("big", "28S_rRNA", 2000, 5500),
    gff_line("mid", "18S_rRNA", 1, 1700)
  )
  res <- select_rrna_contigs(gff, max_contigs = 2)

  expect_equal(res$candidates, c("big", "mid"))
  expect_equal(res$evidence$rank, 1:3)
  expect_match(res$evidence$reason[res$evidence$contig == "small"], "beyond the 2 best")
})

test_that("max_contigs of NA or zero means no cap", {
  gff <- c(gff_line("a", "18S_rRNA", 1, 1800), gff_line("b", "28S_rRNA", 1, 3000))
  expect_length(select_rrna_contigs(gff, max_contigs = NA)$candidates, 2L)
  expect_length(select_rrna_contigs(gff, max_contigs = 0)$candidates, 2L)
})
