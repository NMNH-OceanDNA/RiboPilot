# compute_ref_alignment() stores per-element identity and coverage, not the
# aligned strings themselves.

mk_aln_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_alignment (ID TEXT, path INT, scaffold INT,
    accession TEXT, element TEXT, pid REAL, aln_cov REAL, aln_len INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
    element TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT,
    rank INT, accession TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_sequences(accession TEXT, sequence TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_annotations(accession TEXT, element TEXT,
    start INT, stop INT, strand TEXT)")

  s <- paste(rep("ACGTACGTAC", 20), collapse = "")   # 200 bp
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,'18S',1,200)")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates VALUES('S1',1,1,1,'ACC1')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences VALUES('ACC1',?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations VALUES('ACC1','18S',1,200,'+')")
  con
}

test_that("compute_ref_alignment writes identity and coverage but no aligned strings", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 1L)
  r <- DBI::dbGetQuery(con, "SELECT * FROM blast_ref_alignment")
  expect_equal(nrow(r), 1L)
  expect_equal(r$element[1], "18S")
  expect_equal(r$pid[1], 100)
  expect_equal(r$aln_len[1], 200L)
  expect_false("sample_aln" %in% names(r))
})

test_that("a degenerate element alignment is not stored", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  # replace the reference with a short unrelated spacer so the ends-free
  # alignment covers almost none of the sample element
  DBI::dbExecute(con, "UPDATE blast_ref_sequences SET sequence = ? WHERE accession = 'ACC1'",
    params = list(paste(rep("TTTTTTTTTT", 3), collapse = "")))
  DBI::dbExecute(con, "UPDATE blast_ref_annotations SET stop = 30 WHERE accession = 'ACC1'")
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM blast_ref_alignment")$n, 0L)
})

test_that("a good alignment is still stored", {
  con <- mk_aln_db()
  on.exit(DBI::dbDisconnect(con))
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)
  expect_equal(as.integer(n), 1L)
})

aln_seed <- function(con) {
  s <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
    VALUES ('S1', 1, 1, 200, ?)", params = list(s))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length)
    VALUES ('S1', 1, 1, '18S', 1, 200, 200)")
  DBI::dbExecute(con, "INSERT INTO blast_ref_candidates (ID, path, scaffold, rank, accession, species)
    VALUES ('S1', 1, 1, 1, 'ACC1', 'Danio rerio')")
  DBI::dbExecute(con, "INSERT INTO blast_ref_sequences (accession, sequence) VALUES ('ACC1', ?)",
    params = list(s))
  DBI::dbExecute(con, "INSERT INTO blast_ref_annotations (accession, element, start, stop, strand)
    VALUES ('ACC1', '18S', 1, 200, '+')")
  s
}

# The copy: one base changed every 10 bp over the first 100, and 18S cut to 1..100.
aln_edit_copy <- function(con, s) {
  b <- strsplit(s, "")[[1]]
  i <- seq(10, 100, by = 10)
  b[i] <- chartr("ACGT", "CATG", b[i])
  m <- paste(b, collapse = "")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?", params = list(m))
  DBI::dbExecute(con, "UPDATE review_annotations SET stop = 100, length = 100")
  m
}

test_that("compute_ref_alignment(review = TRUE) aligns the copy and writes only its mirror", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  s <- aln_seed(con)
  expect_equal(as.integer(RiboPilot:::compute_ref_alignment(con, "S1", 1, 1)), 1L)
  aln_edit_copy(con, s)
  n <- RiboPilot:::compute_ref_alignment(con, "S1", 1, 1, review = TRUE)
  expect_equal(as.integer(n), 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT element, pid, aln_len FROM review_ref_alignment"),
    data.frame(element = "18S", pid = 90, aln_len = 100L))
  expect_equal(DBI::dbGetQuery(con, "SELECT element, pid, aln_len FROM blast_ref_alignment"),
    data.frame(element = "18S", pid = 100, aln_len = 200L))
})

test_that("the live alignment reads the copy when there is one, else the original", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  s <- aln_seed(con)
  sample_row <- function() {
    m <- RiboPilot:::compute_feature_msa(con, "S1", 1L, 1L, "18S")
    gsub("-", "", as.character(m)[["sample|18S"]])
  }
  expect_equal(sample_row(), s)
  m <- aln_edit_copy(con, s)
  expect_equal(sample_row(), substr(m, 1, 100))
})
