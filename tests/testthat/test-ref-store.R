# The shared store is read by many projects and, on a cluster, by many users.
# Reading must not require write access.

mk_ro_db <- function(dir) {
  p <- file.path(dir, "ro.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "CREATE TABLE t (a TEXT)")
  DBI::dbExecute(con, "INSERT INTO t VALUES ('x')")
  DBI::dbDisconnect(con)
  Sys.chmod(p, "0444")
  p
}

test_that("open_sqlite_ro reads a write-protected database", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT a FROM t")$a, "x")
})

test_that("open_sqlite_ro refuses to write", {
  d <- withr::local_tempdir()
  p <- mk_ro_db(d)
  con <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(con))
  expect_error(DBI::dbExecute(con, "INSERT INTO t VALUES ('y')"))
})

test_that("open_sqlite_ro errors on a missing file rather than creating one", {
  d <- withr::local_tempdir()
  p <- file.path(d, "absent.sqlite")
  expect_error(RiboPilot:::open_sqlite_ro(p))
  expect_false(file.exists(p))
})

test_that("open_sqlite_ro refuses writes even when the file is writable", {
  d <- withr::local_tempdir()
  p <- file.path(d, "rw.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "CREATE TABLE t (a TEXT)")
  DBI::dbExecute(con, "INSERT INTO t VALUES ('x')")
  DBI::dbDisconnect(con)
  ro <- RiboPilot:::open_sqlite_ro(p)
  on.exit(DBI::dbDisconnect(ro))
  expect_equal(DBI::dbGetQuery(ro, "SELECT a FROM t")$a, "x")
  expect_error(DBI::dbExecute(ro, "INSERT INTO t VALUES ('y')"))
})

mk_store_dir <- function() withr::local_tempdir(.local_envir = parent.frame())

seq_fixture <- function() data.frame(
  accession = c("A.1", "B.1"), species = c("Genus a", "Genus b"),
  title = c("Genus a 18S", "Genus b 28S"), sequence = c("ACGT", "TTTT"),
  length = c(4L, 4L), stringsAsFactors = FALSE)

ann_fixture <- function() data.frame(
  accession = c("A.1", "A.1", "B.1"), element = c("18S", "ITS1", "28S"),
  type = c("rRNA", "spacer", "rRNA"), start = c(1L, 3L, 1L), stop = c(2L, 4L, 4L),
  strand = c("+", "+", "-"), stringsAsFactors = FALSE)

test_that("a store round-trips sequences and annotations", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  n <- RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(),
    kingdom = "euk", taxon = "all", tool = "pybarrnap")
  expect_equal(as.integer(n), 2L)
  got <- RiboPilot:::ref_store_get(st, c("A.1", "B.1"))
  expect_equal(sort(got$seqs$accession), c("A.1", "B.1"))
  expect_equal(nrow(got$ann), 3L)
  expect_equal(sort(got$ann$element[got$ann$accession == "A.1"]), c("18S", "ITS1"))
})

test_that("ref_store_get returns only what was asked for", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  got <- RiboPilot:::ref_store_get(st, "A.1")
  expect_equal(got$seqs$accession, "A.1")
  expect_true(all(got$ann$accession == "A.1"))
})

test_that("ref_store_get on an unknown accession returns empty frames, not an error", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  got <- RiboPilot:::ref_store_get(st, "ZZ.9")
  expect_equal(nrow(got$seqs), 0L)
  expect_equal(nrow(got$ann), 0L)
})

test_that("re-putting an accession replaces rather than duplicating", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_sequences")$n, 2L)
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_annotations")$n, 3L)
})

test_that("the store records how its coordinates were produced", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "M", "pybarrnap")
  r <- DBI::dbGetQuery(st, "SELECT kingdom, taxon, tool FROM ref_sequences WHERE accession='A.1'")
  expect_equal(r$kingdom, "euk")
  expect_equal(r$taxon, "M")
  expect_equal(r$tool, "pybarrnap")
})

test_that("open_ref_store returns NULL for a store that does not exist", {
  d <- mk_store_dir()
  expect_null(RiboPilot:::open_ref_store(d, write = FALSE))
})

test_that("an existing store opens read-only without a write lock", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)
  Sys.chmod(RiboPilot:::ref_store_path(d), "0444")
  ro <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(ro))
  expect_equal(nrow(RiboPilot:::ref_store_get(ro, "A.1")$seqs), 1L)
})

mk_project_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT,
    rank INT, accession TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences(accession TEXT PRIMARY KEY,
    species TEXT, title TEXT, taxid TEXT, lineage TEXT, sequence TEXT, length INT,
    time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations(accession TEXT, element TEXT,
    type TEXT, start INT, stop INT, strand TEXT, time_stamp INT,
    PRIMARY KEY (accession, element))")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,1,'A.1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,2,'ZZ.9')")
  con
}

test_that("priming copies a stored accession into the project and skips unknown ones", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::prime_from_ref_store(con, d)
  expect_equal(as.integer(n), 1L)
  s <- DBI::dbGetQuery(con, "SELECT accession FROM blast_ref_sequences")
  expect_equal(s$accession, "A.1")
  a <- DBI::dbGetQuery(con, "SELECT element FROM blast_ref_annotations ORDER BY element")
  expect_equal(a$element, c("18S", "ITS1"))
})

test_that("priming is a no-op when no store exists", {
  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  d <- mk_store_dir()
  expect_equal(as.integer(RiboPilot:::prime_from_ref_store(con, d)), 0L)
})

test_that("priming does not overwrite an accession the project already has", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence) VALUES ('A.1','KEEPME')")
  RiboPilot:::prime_from_ref_store(con, d)
  expect_equal(DBI::dbGetQuery(con, "SELECT sequence FROM blast_ref_sequences WHERE accession='A.1'")$sequence,
    "KEEPME")
})
