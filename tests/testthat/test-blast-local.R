# Composing the -db argument. blastn takes a space-separated list, and every
# database in it must carry accession info or the search errors outright.

mk_fake_db <- function(dir, prefix) {
  for (ext in c("nhr", "nin", "nsq")) {
    file.create(file.path(dir, paste0(prefix, ".", ext)))
  }
  file.path(dir, prefix)
}

test_that("a panel alone composes to itself", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  expect_equal(RiboPilot:::local_blast_db_arg(p), p)
})

test_that("decoys present are appended in order", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  dd <- withr::local_tempdir()
  for (n in RiboPilot:::RP_DECOY_DBS) mk_fake_db(dd, n)
  got <- RiboPilot:::local_blast_db_arg(p, dd)
  parts <- strsplit(got, " ", fixed = TRUE)[[1]]
  expect_equal(parts[1], p)
  expect_equal(basename(parts[-1]), RiboPilot:::RP_DECOY_DBS)
})

test_that("a missing decoy is skipped with a warning, not an error", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  dd <- withr::local_tempdir()
  mk_fake_db(dd, "16S_ribosomal_RNA")
  expect_warning(got <- RiboPilot:::local_blast_db_arg(p, dd), "SSU_eukaryote_rRNA")
  parts <- strsplit(got, " ", fixed = TRUE)[[1]]
  expect_length(parts, 2L)
})

test_that("a blank or absent panel yields NULL", {
  expect_null(RiboPilot:::local_blast_db_arg(""))
  expect_null(RiboPilot:::local_blast_db_arg(NULL))
  expect_null(RiboPilot:::local_blast_db_arg(NA_character_))
  d <- withr::local_tempdir()
  expect_null(RiboPilot:::local_blast_db_arg(file.path(d, "absent")))
})

test_that("no decoy directory means panel only, no warning", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  expect_silent(got <- RiboPilot:::local_blast_db_arg(p, NULL))
  expect_equal(got, p)
})
