# The in-app reference-database step's tables on a fresh project database.

test_that("a fresh project database has the reference tables and stores ref_dir", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  mapping <- data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = file.path(dir, "S1_R1.fq.gz"), R2 = file.path(dir, "S1_R2.fq.gz"))
  file.create(mapping$R1, mapping$R2)
  mapping_fn <- file.path(dir, "mapping.csv")
  utils::write.csv(mapping, mapping_fn, row.names = FALSE)
  refs <- withr::local_tempdir()
  seeds_db <- testthat::test_path("..", "..", "ref_dbs", "getOrganelle", "seeds",
    "actinopterygii_rrna_seeds.fasta")
  labels_db <- testthat::test_path("..", "..", "ref_dbs", "getOrganelle", "labels",
    "actinopterygii_rrna_labels.fasta")
  skip_if_not(file.exists(seeds_db) && file.exists(labels_db))
  new_db(db_path = file.path(dir, ".sqlite"), mapping_fn = mapping_fn, ref_dir = refs,
    seeds_db = seeds_db, labels_db = labels_db)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_true("reference" %in% DBI::dbListTables(con))
  expect_true("ref_id" %in% DBI::dbListFields(con, "samples"))
  expect_equal(DBI::dbGetQuery(con, "SELECT ref_dir FROM reference_config")$ref_dir[1], refs)
})
