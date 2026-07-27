# Schema v9: retained hit distribution + rank-aligned lineage columns.
# Migrations must be additive and idempotent, and must not disturb existing rows.

# A minimal v8-shaped project DB, as a project created by the previous release
# would look.
mk_v8_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT PRIMARY KEY, taxid INT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT,
    blast_accession TEXT, blast_species TEXT, blast_pident REAL, blast_qcovs REAL, blast_evalue REAL,
    ignore INT, PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT,
    title TEXT, species TEXT, staxid INT, lineage TEXT, pident REAL, qcovs REAL, evalue REAL,
    is_nuclear_rrna INT, time_stamp INT, PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, evalue REAL, time_stamp INT,
    PRIMARY KEY (ID, path, scaffold, rank))")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa; Chordata')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'ACGT','AA1','Danio rerio',99,80,1e-90,0)")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,'AA1','t','Danio rerio',7955,'Metazoa',99,80,1e-90,1,1)")
  DBI::dbExecute(con, "PRAGMA user_version = 8")
  con
}

test_that("migrating v8 to v9 adds the new tables and columns", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("blast_hits" %in% DBI::dbListTables(con))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 9L)
  expect_true("blast_hits" %in% DBI::dbListTables(con))
  expect_true("blast_source" %in% DBI::dbListFields(con, "assemblies"))
  expect_true("lineage_taxids" %in% DBI::dbListFields(con, "samples"))
  for (f in c("lineage_taxids", "rank_distance", "bitscore")) {
    expect_true(f %in% DBI::dbListFields(con, "blast_screen"), info = f)
  }
  expect_true("bitscore" %in% DBI::dbListFields(con, "blast_ref_candidates"))
})

test_that("migration preserves existing rows and leaves new columns NULL", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  s <- DBI::dbGetQuery(con, "SELECT * FROM samples")
  expect_equal(s$ID, "S1")
  expect_equal(s$taxid, 7955L)
  expect_equal(s$lineage, "Metazoa; Chordata")
  expect_true(is.na(s$lineage_taxids))

  a <- DBI::dbGetQuery(con, "SELECT * FROM assemblies")
  expect_equal(a$blast_accession, "AA1")
  expect_equal(a$sequence, "ACGT")
  expect_true(is.na(a$blast_source))   # provenance unknown for pre-v9 rows

  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_screen")$n, 1)
})

test_that("migrate_db is idempotent", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession)
                       VALUES('S1',1,1,1,1,'AA1')")
  expect_silent(migrate_db(con))
  expect_silent(migrate_db(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 9L)
  # a second migrate must not wipe rows written after the first
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_hits")$n, 1)
})

test_that("a fresh new_db() matches a migrated one", {
  dir <- withr::local_tempdir()
  map <- file.path(dir, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", R1 = "a_R1.fq", R2 = "a_R2.fq"),
    map, row.names = FALSE)
  db <- file.path(dir, ".sqlite")
  new_db(db_path = db, mapping_fn = map,
    seeds_db = map, labels_db = map)   # any existing file passes the DB check
  fresh <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(fresh), add = TRUE)

  expect_equal(DBI::dbGetQuery(fresh, "PRAGMA user_version")[[1]], 9L)
  expect_true("blast_hits" %in% DBI::dbListTables(fresh))
  expect_true("blast_source" %in% DBI::dbListFields(fresh, "assemblies"))
  expect_true("lineage_taxids" %in% DBI::dbListFields(fresh, "samples"))

  mig <- mk_v8_db(); on.exit(DBI::dbDisconnect(mig), add = TRUE)
  migrate_db(mig)
  # the columns v9 introduces must exist in both, or code diverges by vintage
  for (tbl in c("blast_screen", "blast_ref_candidates", "blast_hits")) {
    expect_true(all(DBI::dbListFields(mig, tbl) %in% DBI::dbListFields(fresh, tbl)), info = tbl)
  }
})
