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
  expected_decoys <- file.path(dd, RiboPilot:::RP_DECOY_DBS)
  expect_equal(parts[-1], expected_decoys)
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

test_that("a stray .notes file without .nin is not a database", {
  d <- withr::local_tempdir()
  p <- file.path(d, "ribopilot_rrna")
  file.create(paste0(p, ".notes"))
  expect_null(RiboPilot:::local_blast_db_arg(p))
})

test_that("a prefix with .nin is detected as a database", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  result <- RiboPilot:::local_blast_db_arg(p)
  expect_equal(result, p)
})

test_that("a panel prefix containing a space warns and returns NULL", {
  d <- withr::local_tempdir()
  p <- file.path(d, "ribopilot rrna")
  file.create(paste0(p, ".nhr"))
  file.create(paste0(p, ".nin"))
  file.create(paste0(p, ".nsq"))
  expect_warning(got <- RiboPilot:::local_blast_db_arg(p), "whitespace")
  expect_null(got)
})

test_that("a decoy path containing a space warns, is skipped, and panel is returned", {
  d <- withr::local_tempdir()
  p <- mk_fake_db(d, "ribopilot_rrna")
  dd <- withr::local_tempdir()
  mk_fake_db(dd, "16S_ribosomal_RNA")
  decoy_with_space <- file.path(dd, "SSU eukaryote_rRNA")
  file.create(paste0(decoy_with_space, ".nhr"))
  file.create(paste0(decoy_with_space, ".nin"))
  file.create(paste0(decoy_with_space, ".nsq"))
  expect_warning(got <- RiboPilot:::local_blast_db_arg(p, dd), "whitespace")
  parts <- strsplit(got, " ", fixed = TRUE)[[1]]
  expect_equal(parts, p)
})
