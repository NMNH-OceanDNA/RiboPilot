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
  n <- RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), tool = "pybarrnap")
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
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
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

test_that("ref_store_get returns a row written under the current constants, and not one written under different values", {
  # Pins the self-healing property: if REF_ANNOT_KINGDOM/REF_ANNOT_TAXON are
  # ever changed, rows written under the old values must stop matching so
  # they get re-fetched and re-annotated rather than silently reused.
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))

  cur <- seq_fixture()[1, , drop = FALSE]
  cur$accession <- "CUR.1"
  RiboPilot:::ref_store_put(st, cur, NULL, "pybarrnap")

  # a row as if written under different kingdom/taxon values, bypassing
  # ref_store_put, which can now only ever write the current constants
  DBI::dbExecute(st, "INSERT INTO ref_sequences
    (accession, species, title, sequence, length, kingdom, taxon, tool, ribopilot_version, stored_utc)
    VALUES ('OLD.1', 'Genus x', 'Genus x 18S', 'ACGT', 4, 'bac', 'all', 'pybarrnap', '0.0.0', '2020-01-01T00:00:00Z')")

  got <- RiboPilot:::ref_store_get(st, c("CUR.1", "OLD.1"))
  expect_equal(got$seqs$accession, "CUR.1")
})

test_that("re-putting an accession replaces rather than duplicating", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_sequences")$n, 2L)
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_annotations")$n, 3L)
})

test_that("the store records how its coordinates were produced", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
  r <- DBI::dbGetQuery(st, "SELECT kingdom, taxon, tool FROM ref_sequences WHERE accession='A.1'")
  expect_equal(r$kingdom, RiboPilot:::REF_ANNOT_KINGDOM)
  expect_equal(r$taxon, RiboPilot:::REF_ANNOT_TAXON)
  expect_equal(r$tool, "pybarrnap")
})

test_that("open_ref_store returns NULL for a store that does not exist", {
  d <- mk_store_dir()
  expect_null(RiboPilot:::open_ref_store(d, write = FALSE))
})

test_that("an existing store opens read-only without a write lock", {
  # Covers a read-only DIRECTORY, not just a read-only file: the recommended
  # HPC deployment is a read-only directory, and WAL (the pre-fix default on a
  # local filesystem) needs to create a -shm file even just to read, which a
  # read-only directory refuses.
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
  DBI::dbDisconnect(st)
  p <- RiboPilot:::ref_store_path(d)
  store_dir <- dirname(p)
  Sys.chmod(p, "0444")
  Sys.chmod(store_dir, "0555")
  # restore before withr's tempdir cleanup runs, or a still-read-only directory
  # cannot be cleaned up; after = FALSE runs this ahead of that cleanup
  on.exit({
    Sys.chmod(store_dir, "0755")
    Sys.chmod(p, "0644")
  }, add = TRUE, after = FALSE)
  ro <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(ro), add = TRUE)
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
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
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
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
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
  n <- RiboPilot:::push_to_ref_store(d, fa, ann, "pybarrnap")
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
  expect_equal(as.integer(RiboPilot:::push_to_ref_store(d, fa, NULL, "pybarrnap")), 0L)
})

test_that("push_to_ref_store fails closed on a blank or NA ref_dir", {
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "ACGT"), fa)

  expect_equal(as.integer(RiboPilot:::push_to_ref_store("", fa, NULL, "pybarrnap")), 0L)
  expect_false(dir.exists("/refstore"))

  expect_equal(as.integer(RiboPilot:::push_to_ref_store(NA_character_, fa, NULL, "pybarrnap")), 0L)
  expect_false(dir.exists(file.path(wd, "NA")))
})

test_that("prime_from_ref_store fails closed on a blank or NA ref_dir", {
  wd <- withr::local_tempdir()
  withr::local_dir(wd)

  con <- mk_project_db()
  on.exit(DBI::dbDisconnect(con))
  expect_equal(as.integer(RiboPilot:::prime_from_ref_store(con, "")), 0L)

  # plant a real store at the bogus relative path ref_store_path(NA) resolves
  # to, so only a fixed guard stands between this and a coincidental read
  decoy <- RiboPilot:::open_ref_store("NA", write = TRUE)
  RiboPilot:::ref_store_put(decoy, seq_fixture(), ann_fixture(), "pybarrnap")
  DBI::dbDisconnect(decoy)

  expect_equal(as.integer(RiboPilot:::prime_from_ref_store(con, NA_character_)), 0L)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT accession FROM blast_ref_sequences")), 0L)
})

test_that("a stored row carries honest pipeline provenance", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "ribopilot-annotate")
  r <- DBI::dbGetQuery(st, "SELECT tool, ribopilot_version FROM ref_sequences WHERE accession='A.1'")
  expect_equal(r$tool, "ribopilot-annotate")
  expect_true(nzchar(r$ribopilot_version))
})

test_that("a per-accession insert collision is skipped, not fatal to the batch", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  RiboPilot:::ref_store_put(st, seq_fixture(), ann_fixture(), "pybarrnap")
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

test_that("ann$tool round-trips into ref_annotations", {
  d <- mk_store_dir()
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 Genus a 18S ribosomal RNA", "ACGT"), fa)
  ann <- data.frame(seqid = "A.1", element = "18S", type = "rRNA",
    start = 1L, stop = 4L, strand = "+", tool = "barrnap", stringsAsFactors = FALSE)
  n <- RiboPilot:::push_to_ref_store(d, fa, ann, "ribopilot-annotate")
  expect_equal(as.integer(n), 1L)

  st <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(st))
  r <- DBI::dbGetQuery(st, "SELECT tool FROM ref_annotations WHERE accession='A.1' AND element='18S'")
  expect_equal(r$tool, "barrnap")
})

test_that("ref_dir_ok does not error on a zero-length input", {
  expect_false(RiboPilot:::ref_dir_ok(character(0)))
  expect_false(RiboPilot:::ref_dir_ok(NULL))
})

test_that("pushing skips a record with an empty sequence, like write_blast_refs does", {
  d <- mk_store_dir()
  fa <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 Genus a 18S ribosomal RNA", "ACGT",
               ">EMPTY.1 nothing here", ""), fa)
  n <- RiboPilot:::push_to_ref_store(d, fa, NULL, "pybarrnap")
  expect_equal(as.integer(n), 1L)

  st <- RiboPilot:::open_ref_store(d, write = FALSE)
  on.exit(DBI::dbDisconnect(st))
  got <- RiboPilot:::ref_store_get(st, c("A.1", "EMPTY.1"))
  expect_equal(got$seqs$accession, "A.1")
})

test_that("a mid-loop failure in ref_store_put returns the partial count, not zero", {
  d <- mk_store_dir()
  st <- RiboPilot:::open_ref_store(d, write = TRUE)
  on.exit(DBI::dbDisconnect(st))
  bad <- seq_fixture()
  bad$accession[2] <- NA_character_  # violates NOT NULL, forcing the second row to fail
  n <- RiboPilot:::ref_store_put(st, bad, ann_fixture(), "pybarrnap")
  expect_equal(as.integer(n), 1L)
  expect_equal(DBI::dbGetQuery(st, "SELECT COUNT(*) n FROM ref_sequences")$n, 1L)
})
