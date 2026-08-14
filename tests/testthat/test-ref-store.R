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

test_that("a corrupt store is a no-op rather than an error", {
  d <- mk_store_dir()
  p <- RiboPilot:::ref_store_path(d)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  writeBin(charToRaw("this is not a sqlite file"), p)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  n <- NULL
  expect_warning(n <- RiboPilot:::prime_from_ref_store(con, d), "reference store")
  expect_equal(as.integer(n), 0L)
})

test_that("pushing writes fetched records into the store", {
  d <- mk_store_dir()
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 Genus a 18S ribosomal RNA", "ACGT",
               ">B.1 Genus b 28S ribosomal RNA", "TTTT"), fa)
  ann <- data.frame(seqid = c("A.1", "B.1"), element = c("18S", "28S"),
    type = c("rRNA", "rRNA"), start = c(1L, 1L), stop = c(4L, 4L),
    strand = c("+", "-"), stringsAsFactors = FALSE)
  n <- RiboPilot:::push_to_ref_store(d, fa, ann, "euk", "all", "pybarrnap")
  expect_equal(as.integer(n), 2L)

  st <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(st))
  got <- RiboPilot:::ref_store_get(st, c("A.1", "B.1"))
  expect_equal(sort(got$seqs$accession), c("A.1", "B.1"))
  expect_equal(got$seqs$species[got$seqs$accession == "A.1"], "Genus a")
  expect_equal(nrow(got$ann), 2L)
})

test_that("pushing to an unwritable reference folder is not an error", {
  d <- mk_store_dir()
  Sys.chmod(d, "0555")
  on.exit(Sys.chmod(d, "0755"), add = TRUE)
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "ACGT"), fa)
  expect_equal(as.integer(RiboPilot:::push_to_ref_store(d, fa, NULL, "euk", "all", "pybarrnap")), 0L)
})

test_that("a per-accession insert collision is skipped, not fatal to the batch", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "euk", "all", "pybarrnap")
  DBI::dbDisconnect(st)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  # Simulate the race: another process primes A.1 concurrently, between this
  # process's SELECT of candidates and this process's own INSERT for A.1.
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence) VALUES ('A.1','RACED')")
  seq_row <- data.frame(accession = "A.1", species = "Genus a", title = "Genus a 18S",
    sequence = "ACGT", length = 4L, stringsAsFactors = FALSE)
  ann_rows <- ann_fixture()[ann_fixture()$accession == "A.1", , drop = FALSE]
  ts <- as.integer(Sys.time())
  ok1 <- NULL
  expect_warning(ok1 <- RiboPilot:::prime_copy_accession(con, "A.1", seq_row, ann_rows, ts), "A.1")
  expect_false(ok1)
  expect_equal(DBI::dbGetQuery(con, "SELECT sequence FROM blast_ref_sequences WHERE accession='A.1'")$sequence,
    "RACED")

  # the collision must not wedge the connection for the next accession
  seq_row2 <- data.frame(accession = "B.1", species = "Genus b", title = "Genus b 28S",
    sequence = "TTTT", length = 4L, stringsAsFactors = FALSE)
  ok2 <- RiboPilot:::prime_copy_accession(con, "B.1", seq_row2, NULL, ts)
  expect_true(ok2)
  expect_equal(DBI::dbGetQuery(con, "SELECT accession FROM blast_ref_sequences ORDER BY accession")$accession,
    c("A.1", "B.1"))
})
