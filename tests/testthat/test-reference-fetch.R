# WF1 reference fetch: what efetch hands us is not always a record.

mk_ref_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences (accession TEXT PRIMARY KEY,
    species TEXT, title TEXT, taxid TEXT, lineage TEXT, sequence TEXT, length INTEGER,
    time_stamp INTEGER)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations (accession TEXT, element TEXT,
    type TEXT, start INTEGER, stop INTEGER, strand TEXT, time_stamp INTEGER,
    PRIMARY KEY (accession, element))")
  con
}

test_that("a headerless or truncated file is rejected", {
  d <- withr::local_tempdir()
  p <- file.path(d, "bare.fa")
  writeLines(c("", "ACGTACGTACGTACGT"), p)
  expect_false(RiboPilot:::is_fasta_file(p))
})

test_that("blast_ref_fetch validates the efetch response before storing it", {
  nf <- readLines(system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot"))
  txt <- paste(nf, collapse = "\n")
  expect_true(grepl("<ERROR>", txt, fixed = TRUE))
  expect_true(grepl("<eFetchResult>", txt, fixed = TRUE))
  expect_true(grepl("F a i l e d", txt, fixed = TRUE))
  expect_true(grepl("^Error", txt, fixed = TRUE))
  # a failed chunk is logged rather than aborting the whole batch
  expect_true(any(grepl("efetch failed for chunk", nf, fixed = TRUE)))
  # the sqlite timeout stays on every sqlite3 call
  sq <- grep("sqlite3", nf, value = TRUE)
  expect_true(all(grepl(".timeout 10000", sq, fixed = TRUE)))
})

test_that("write_blast_refs stores every record in a multi-fasta", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(
    ">ACC1.1 Danio rerio 18S ribosomal RNA gene",
    paste(rep("ACGT", 50), collapse = ""),
    ">ACC2.1 Danio rerio 28S ribosomal RNA gene",
    paste(rep("TGCA", 60), collapse = "")), f)
  # seqid order deliberately does NOT match the fasta record order (ACC1.1
  # then ACC2.1): write_blast_refs must join ann by seqid, not positionally
  ann <- data.frame(
    seqid = c("ACC2.1", "ACC1.1"), element = c("28S", "18S"),
    type = c("rRNA", "rRNA"), start = c(1L, 1L), stop = c(240L, 200L),
    length = c(240L, 200L), strand = c("+", "+"),
    model_cov = c(100, 100), tool = c("barrnap", "barrnap"),
    stringsAsFactors = FALSE)
  n <- RiboPilot:::write_blast_refs(con, f, ann)
  expect_equal(as.integer(n), 2L)
  s <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_sequences ORDER BY accession")
  expect_equal(s$accession, c("ACC1.1", "ACC2.1"))
  expect_equal(nchar(s$sequence[1]), 200L)
  expect_equal(s$species[1], "Danio rerio")
  a <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_annotations ORDER BY accession")
  expect_equal(a$element, c("18S", "28S"))
})

test_that("write_blast_refs does not concatenate records into one sequence", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "AAAA", ">B.1 y", "CCCC"), f)
  RiboPilot:::write_blast_refs(con, f, NULL)
  s <- DBI::dbGetQuery(con, "SELECT accession, sequence FROM blast_ref_sequences ORDER BY accession")
  expect_equal(s$sequence, c("AAAA", "CCCC"))
})

test_that("write_blast_refs replaces rather than duplicating on a re-run", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", "AAAA"), f)
  RiboPilot:::write_blast_refs(con, f, NULL)
  RiboPilot:::write_blast_refs(con, f, NULL)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_sequences")$n, 1L)
})

test_that("a PDB chain is stored under the accession BLAST reports", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">pdb|4V6X|A5 Chain A5, 28S ribosomal RNA", strrep("ACGT", 100),
    ">MN000001.1 Danio rerio 28S ribosomal RNA", strrep("TGCA", 100)), f)
  RiboPilot:::normalize_pdb_deflines(f)
  expect_equal(readLines(f)[1], ">4V6X_A5 Chain A5, 28S ribosomal RNA")
  RiboPilot:::write_blast_refs(con, f, NULL)
  expect_setequal(DBI::dbGetQuery(con, "SELECT accession FROM blast_ref_sequences")$accession,
    c("4V6X_A5", "MN000001.1"))
  bin <- readLines(system.file("nextflow", "bin", "rp_ref_store.R", package = "RiboPilot"))
  expect_lt(grep("normalize_pdb_deflines", bin)[1], grep("annotate_records", bin)[1])
})

test_that("stop_if_nothing_stored fires only when input existed but nothing landed", {
  con <- mk_ref_db()
  on.exit(DBI::dbDisconnect(con))
  # two headers, no sequence lines: both records fail write_blast_refs' nzchar() check
  f <- withr::local_tempfile(fileext = ".fa")
  writeLines(c(">A.1 x", ">B.1 y"), f)

  n <- RiboPilot:::write_blast_refs(con, f, NULL)
  expect_equal(as.integer(n), 0L)

  expect_error(RiboPilot:::stop_if_nothing_stored(f, n), "stored nothing")
  # storing at least one record clears the guard
  expect_equal(RiboPilot:::stop_if_nothing_stored(f, 1L), 1L)
  # an empty fasta (no input records) is not this guard's job either
  empty <- withr::local_tempfile(fileext = ".fa")
  file.create(empty)
  expect_equal(RiboPilot:::stop_if_nothing_stored(empty, 0L), 0L)
})

test_that("the fetch module asks for fasta only, not feature tables", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- readLines(p, warn = FALSE)
  expect_false(any(grepl("-format ft", txt, fixed = TRUE)))
  expect_true(any(grepl("-format fasta", txt, fixed = TRUE)))
})

test_that("the fetch module batches accessions rather than looping one at a time", {
  p <- system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast_ref.nf"
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_true(grepl("rp_ref_store.R", txt, fixed = TRUE))
  expect_false(grepl("rp_ref_fetch.R", txt, fixed = TRUE))
})
