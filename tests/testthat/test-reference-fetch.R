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

# 30 bases of leading spacer carrying IUPAC codes, then a 20-base "18S" at 31..50.
ref_seq_parts <- list(
  lead = "ACGTRYSWKMBDHVACGTNRYACGTACGTAC",
  s18  = "GGGGTTTTCCCCAAAAGGGG"
)

write_ref_fixture <- function(dir, acc = "MZ000001") {
  seq <- paste0(ref_seq_parts$lead, ref_seq_parts$s18)
  fa <- file.path(dir, paste0(acc, ".fa"))
  writeLines(c(paste0(">", acc, ".1 Fake species small subunit ribosomal RNA gene"),
    substr(seq, 1, 25), substr(seq, 26, nchar(seq))), fa)
  ft <- file.path(dir, paste0(acc, ".ft"))
  writeLines(c(paste0(">Feature ", acc, ".1"),
    paste(nchar(ref_seq_parts$lead) + 1, nchar(seq), "rRNA", sep = "\t"),
    paste("", "", "", "product", "18S ribosomal RNA", sep = "\t")), ft)
  list(fa = fa, ft = ft, seq = seq, acc = acc)
}

test_that("IUPAC ambiguity codes are kept, so the sequence keeps its length", {
  d <- withr::local_tempdir()
  f <- write_ref_fixture(d)
  seq <- RiboPilot:::read_fasta_seq(f$fa)
  expect_equal(seq, f$seq)
  expect_equal(nchar(seq), nchar(f$seq))
  expect_true(grepl("R", seq, fixed = TRUE))
})

test_that("whitespace and digits are still stripped", {
  d <- withr::local_tempdir()
  p <- file.path(d, "wrapped.fa")
  writeLines(c(">X", "   1 acgt rysw", "  11 kmbd hv"), p)
  expect_equal(RiboPilot:::read_fasta_seq(p), "acgtryswkmbdhv")
})

test_that("the stored sequence and the parsed coordinates land on the same subsequence", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  f <- write_ref_fixture(d)

  n <- RiboPilot:::write_blast_ref(con, f$acc, f$fa, f$ft)
  expect_equal(n, nchar(f$seq))

  stored <- DBI::dbGetQuery(con, "SELECT sequence, length FROM blast_ref_sequences")
  ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM blast_ref_annotations")
  expect_equal(stored$length, nchar(f$seq))
  expect_equal(ann$element, "18S")
  # this is the substr() that blast_ref_align.R applies downstream
  expect_equal(substr(stored$sequence, ann$start, ann$stop), ref_seq_parts$s18)
})

test_that("an XML error body is rejected instead of parsed into a short sequence", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  p <- file.path(d, "err.fa")
  writeLines(c("<?xml version=\"1.0\"?>", "<eFetchResult>",
    "<ERROR>Unable to obtain query #1</ERROR>", "</eFetchResult>"), p)

  expect_equal(RiboPilot:::read_fasta_seq(p), "")
  expect_equal(RiboPilot:::write_blast_ref(con, "MZ000002", p, file.path(d, "none.ft")), 0L)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM blast_ref_sequences")), 0L)
})

test_that("a headerless or truncated file is rejected", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  p <- file.path(d, "bare.fa")
  writeLines(c("", "ACGTACGTACGTACGT"), p)
  expect_false(RiboPilot:::is_fasta_file(p))
  expect_equal(RiboPilot:::write_blast_ref(con, "MZ000003", p, file.path(d, "none.ft")), 0L)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM blast_ref_sequences")), 0L)
})

test_that("an empty or missing result is still handled", {
  con <- mk_ref_db(); on.exit(DBI::dbDisconnect(con))
  d <- withr::local_tempdir()
  empty <- file.path(d, "empty.fa"); file.create(empty)
  expect_equal(RiboPilot:::read_fasta_seq(empty), "")
  expect_equal(RiboPilot:::read_fasta_seq(file.path(d, "gone.fa")), "")
  expect_equal(RiboPilot:::write_blast_ref(con, "MZ000004", empty, file.path(d, "gone.ft")), 0L)
  # a header with no bases is not a record either
  hdr <- file.path(d, "hdr.fa"); writeLines(">MZ000005.1 Fake species", hdr)
  expect_equal(RiboPilot:::write_blast_ref(con, "MZ000005", hdr, file.path(d, "gone.ft")), 0L)
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM blast_ref_sequences")), 0L)
})

test_that("blast_ref_fetch validates the efetch response before storing it", {
  nf <- readLines(system.file("nextflow", "modules", "blast_ref.nf", package = "RiboPilot"))
  txt <- paste(nf, collapse = "\n")
  expect_true(grepl("<ERROR>", txt, fixed = TRUE))
  expect_true(grepl("<eFetchResult>", txt, fixed = TRUE))
  expect_true(grepl("F a i l e d", txt, fixed = TRUE))
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
  ann <- data.frame(
    seqid = c("ACC1.1", "ACC2.1"), element = c("18S", "28S"),
    type = c("rRNA", "rRNA"), start = c(1L, 1L), stop = c(200L, 240L),
    length = c(200L, 240L), strand = c("+", "+"),
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
