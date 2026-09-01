test_that("new_db creates the chimera tables and coverage.clips", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  db <- file.path(d, ".sqlite")
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  tabs <- DBI::dbListTables(con)
  expect_true(all(c("chimera", "chimera_windows") %in% tabs))
  expect_true("clips" %in% DBI::dbListFields(con, "coverage"))
  expect_true("segments" %in% DBI::dbListFields(con, "chimera"))
  expect_true("lineage" %in% DBI::dbListFields(con, "chimera_windows"))
})
