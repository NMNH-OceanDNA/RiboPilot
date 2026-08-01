test_that("write_clade_manifest writes the four columns Nextflow reads", {
  d <- withr::local_tempdir()
  p <- RiboPilot:::write_clade_manifest(file.path(d, "clades.tsv"),
    data.frame(ref_id = "actinopteri", taxid = 186623L, clade_name = "Actinopteri",
      query = '"Actinopteri"[Organism]', stringsAsFactors = FALSE))
  x <- utils::read.delim(p, stringsAsFactors = FALSE)
  expect_equal(names(x), c("ref_id", "taxid", "clade_name", "query"))
  expect_equal(x$taxid, 186623L)
})

test_that("read_count_files reads every count file and skips junk", {
  d <- withr::local_tempdir()
  writeLines("7898\t11750", file.path(d, "7898.count"))
  writeLines("186623\t9000", file.path(d, "186623.count"))
  writeLines("", file.path(d, "6656.count"))
  x <- RiboPilot:::read_count_files(d)
  expect_equal(nrow(x), 2L)
  expect_equal(x$n_records[x$taxid == 7898L], 11750L)
})

test_that("merge_counts replaces stale rows and keeps the rest", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = c(7898L, 33208L), n_records = c(10L, 500L)),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out <- RiboPilot:::merge_counts(d, data.frame(taxid = 7898L, n_records = 11750L))
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 500L)
  expect_equal(nrow(RiboPilot:::read_counts_cache(d)), 2L)
})
