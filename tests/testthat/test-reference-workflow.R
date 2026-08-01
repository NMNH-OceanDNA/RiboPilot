test_that("the reference workflow ships with both entry points", {
  f <- system.file("nextflow", "main_reference.nf", package = "RiboPilot")
  expect_true(nzchar(f))
  txt <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("workflow COUNT", txt, fixed = TRUE))
  expect_true(grepl("workflow BUILD", txt, fixed = TRUE))
  # the work list is read as a header TSV, matching write_clade_manifest()
  expect_true(grepl("splitCsv(header: true, sep: '\\t')", txt, fixed = TRUE))
})

test_that("the task wrapper detaches stdin from Entrez Direct", {
  f <- system.file("nextflow", "bin", "rp_ref_db.R", package = "RiboPilot")
  expect_true(nzchar(f))
  txt <- paste(readLines(f), collapse = "\n")
  expect_true(grepl("/dev/null", txt, fixed = TRUE))
})
