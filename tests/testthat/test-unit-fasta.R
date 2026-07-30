# FASTA text for the Assemble modal's "Copy FASTA" button.

mk_fasta_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT)")
  # Taxon is what the user supplied; sci_name is what taxonomy resolution assigned
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio NCBI')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?)",
    params = list(strrep("ACGT", 40)))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,2,'ACGTACGTAA')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',2,1,'TTTT')")
  con
}

test_that("the header carries the sample, the user's taxon and the unit", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  txt <- RiboPilot:::unit_fasta_text(con, "S1", data.frame(path = 1L, scaffold = 2L))
  expect_equal(txt[1], ">S1 Danio rerio p1.s2")
  expect_equal(txt[2], "ACGTACGTAA")
})

test_that("the assigned scientific name is not used as the taxon", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  txt <- RiboPilot:::unit_fasta_text(con, "S1", data.frame(path = 1L, scaffold = 2L))
  expect_false(any(grepl("NCBI", txt, fixed = TRUE)))
})

test_that("several units come back in one blob, in the order given", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  units <- data.frame(path = c(2L, 1L), scaffold = c(1L, 2L))
  txt <- RiboPilot:::unit_fasta_text(con, "S1", units)
  expect_equal(grep("^>", txt, value = TRUE), c(">S1 Danio rerio p2.s1", ">S1 Danio rerio p1.s2"))
})

test_that("long sequences are wrapped, short ones are not split", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  txt <- RiboPilot:::unit_fasta_text(con, "S1", data.frame(path = 1L, scaffold = 1L))
  seq_lines <- txt[-1]
  expect_true(all(nchar(seq_lines) <= 60))
  expect_equal(paste(seq_lines, collapse = ""), strrep("ACGT", 40))
})

test_that("a unit with no sequence is skipped rather than emitted empty", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',3,1,NULL)")
  txt <- RiboPilot:::unit_fasta_text(con, "S1",
    data.frame(path = c(3L, 1L), scaffold = c(1L, 2L)))
  expect_equal(grep("^>", txt, value = TRUE), ">S1 Danio rerio p1.s2")
})

test_that("no units gives no text", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::unit_fasta_text(con, "S1", data.frame()), character(0))
})

test_that("a missing taxon leaves the header well formed", {
  con <- mk_fasta_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE samples SET Taxon = NULL")
  txt <- RiboPilot:::unit_fasta_text(con, "S1", data.frame(path = 1L, scaffold = 2L))
  expect_equal(txt[1], ">S1 p1.s2")
})
