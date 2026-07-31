# Schema v19: the in-app reference-database step.

mk_v18_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT,
    lineage_taxids TEXT, taxonomy_lock INT, taxonomy_switch INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT PRIMARY KEY, cpus INT,
    memory INT, getOrganelle TEXT, seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INT,
    max_scaffolds INT, min_assembly_length INT, assembler TEXT, bowtie2 TEXT, spades TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa;Chordata', '33208;7711', 1, 2)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1', 1, 0, 'default', 'default')")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default',4,16,'','s.fa','l.fa',NULL,10,10,1000,'getorganelle','--local','--isolate')")
  DBI::dbExecute(con, "PRAGMA user_version = 18")
  con
}

test_that("migrating a v18 project adds the reference tables and ref_id", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("reference" %in% DBI::dbListTables(con))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]],
    RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("reference" %in% DBI::dbListTables(con))
  expect_true("reference_config" %in% DBI::dbListTables(con))
  expect_true("ref_id" %in% DBI::dbListFields(con, "samples"))
  for (f in c("ref_id", "clade_taxid", "clade_name", "clade_rank", "query", "n_samples",
              "n_records", "reference_switch", "source", "seeds_db", "labels_db",
              "log_file", "pid", "fail_reason", "time_stamp")) {
    expect_true(f %in% DBI::dbListFields(con, "reference"), info = f)
  }
})

test_that("the migration seeds exactly one reference_config row with the default window", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  cfg <- DBI::dbGetQuery(con, "SELECT * FROM reference_config")
  expect_equal(nrow(cfg), 1L)
  expect_true(is.na(cfg$ref_dir[1]))
  expect_equal(cfg$min_records[1], 1000)
  expect_equal(cfg$max_records[1], 25000)
})

test_that("the migration is idempotent and leaves existing rows alone", {
  con <- mk_v18_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  DBI::dbExecute(con, "UPDATE reference_config SET min_records = 250")
  migrate_db(con)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM reference_config")), 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT min_records FROM reference_config")$min_records[1], 250)
  expect_equal(DBI::dbGetQuery(con, "SELECT taxid FROM samples")$taxid[1], 7955)
})

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
