# Schema migrations. v9 = retained hit distribution + rank-aligned lineages,
# v10 = median per-position depth.
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
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1', 7955, 'Metazoa; Chordata')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'ACGT','AA1','Danio rerio',99,80,1e-90,0)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',2,0,'default','default')")
  DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,1,'AA1','t','Danio rerio',7955,'Metazoa',99,80,1e-90,1,1)")
  DBI::dbExecute(con, "PRAGMA user_version = 8")
  con
}

test_that("migrating a v8 project adds the v9 tables and columns", {
  con <- mk_v8_db(); on.exit(DBI::dbDisconnect(con))
  expect_false("blast_hits" %in% DBI::dbListTables(con))

  migrate_db(con)

  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
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
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
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

  expect_equal(DBI::dbGetQuery(fresh, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
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

test_that("v10 adds median_depth and back-fills it from the stored tracks", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    min_depth INT, max_depth INT, depths TEXT, time_stamp INT)")
  # a track with a heavy right tail: mean is pulled well above the typical base
  spiky <- paste(c(rep(10, 90), rep(1000, 10)), collapse = ",")
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,1,109,10,1000,?,1)", params = list(spiky))
  DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,2,0,0,0,'',1)")
  DBI::dbExecute(con, "PRAGMA user_version = 9")

  migrate_db(con)
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], RiboPilot:::RIBOPILOT_SCHEMA_VERSION)
  expect_true("median_depth" %in% DBI::dbListFields(con, "coverage"))
  # migration alone leaves it NULL: parsing every track would stall app startup
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT median_depth FROM coverage WHERE scaffold=1")$median_depth))

  expect_equal(backfill_coverage_median(con), 1L)   # the empty track is skipped
  got <- DBI::dbGetQuery(con, "SELECT mean_depth, median_depth FROM coverage WHERE scaffold=1")
  expect_equal(got$median_depth, 10)                # typical base, not 109
  expect_equal(got$mean_depth, 109)
  expect_equal(backfill_coverage_median(con), 0L)   # idempotent
})

test_that("write_coverage stores both statistics", {
  skip_if_not_installed("withr")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    median_depth REAL, min_depth INT, max_depth INT, depths TEXT, time_stamp INT)")
  f <- withr::local_tempfile()
  writeLines(sprintf("ref\t%d\t%d", 1:100, c(rep(10, 90), rep(1000, 10))), f)
  RiboPilot::write_coverage(con, "S1", 1, 1, f)
  d <- DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, max_depth FROM coverage")
  expect_equal(d$median_depth, 10)
  expect_equal(d$mean_depth, 109)
  expect_equal(d$max_depth, 1000L)
})

test_that("v12 adds triage_opts, its seed row, and the rule columns", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  expect_true("triage_opts" %in% DBI::dbListTables(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM triage_opts")$n, 1L)

  # every rule key has a column, and the seed carries the shipped defaults
  r <- DBI::dbGetQuery(con, "SELECT * FROM triage_opts WHERE triage_opts = 'default'")
  for (k in RiboPilot:::triage_rule_keys()) {
    expect_true(paste0("rule_", k) %in% names(r))
  }
  expect_equal(r$rule_no_feat_non_rdna, "drop")
  expect_equal(r$rule_odd_one_out, "flag")

  expect_true("triage_opts" %in% DBI::dbListFields(con, "assemble"))
  expect_true("triage_rule" %in% DBI::dbListFields(con, "assemblies"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 12L)
})
