# Helpers the schema step 4 upgrade uses to write a trimmed contig's original
# back: the coverage summary (median_depth is the column the contig table and
# triage read) and the chimera rows, moved by a shift or put back from the JSON
# an older trim stored in assembly_backup.chimera.

mk_trim_db <- function(dir, seq, depths) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT,
                       length INT, sequence TEXT, edited INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT,
                       element TEXT, start INT, stop INT, length INT, edited INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT,
                       mean_depth REAL, median_depth REAL, min_depth INT, max_depth INT,
                       depths TEXT, clips TEXT, time_stamp INT)")
  DBI::dbExecute(con, RiboPilot:::CHIMERA_DDL)
  DBI::dbExecute(con, RiboPilot:::CHIMERA_WINDOWS_DDL)
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,?,?,0)",
    params = list(nchar(seq), seq))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, mean_depth, median_depth, min_depth,
                       max_depth, depths, time_stamp) VALUES('S1',1,1,?,?,?,?,?,0)",
    params = list(mean(depths), stats::median(depths), min(depths), max(depths),
      paste(depths, collapse = ",")))
  con
}

cov_row <- function(con) {
  DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, min_depth, max_depth, depths
                          FROM coverage WHERE ID='S1' AND path=1 AND scaffold=1")
}

# the unit's chimera rows as JSON, in the form an older trim stored them
chimera_json <- function(con) {
  ch <- DBI::dbGetQuery(con, "SELECT * FROM chimera")
  w <- DBI::dbGetQuery(con, "SELECT * FROM chimera_windows")
  as.character(jsonlite::toJSON(list(chimera = ch, windows = w), dataframe = "rows", na = "null", digits = NA))
}

test_that("write_coverage_summary rewrites median_depth with the track and skips an empty one", {
  dir <- withr::local_tempdir()
  dep <- c(rep(2, 50), rep(20, 50))
  con <- mk_trim_db(dir, strrep("ACGT", 25), dep)
  on.exit(DBI::dbDisconnect(con))

  expect_true(RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep[51:100]))
  cv <- cov_row(con)
  expect_equal(cv$median_depth, 20)
  expect_equal(cv$mean_depth, 20)
  expect_equal(cv$min_depth, 20L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 50L)

  RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep)
  cv <- cov_row(con)
  expect_equal(cv$median_depth, stats::median(dep))
  expect_equal(cv$mean_depth, mean(dep))
  expect_equal(cv$min_depth, 2L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 100L)

  expect_false(RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, numeric(0)))
  expect_equal(cov_row(con), cv)
})

test_that("chimera_shift moves junctions and windows both ways, and chimera_restore puts the JSON back", {
  dir <- withr::local_tempdir()
  seq <- paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")
  con <- mk_trim_db(dir, seq, rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1',1,1,'severe','30,70')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  ch <- function() DBI::dbGetQuery(con, "SELECT junctions, segments FROM chimera")
  win <- function() DBI::dbGetQuery(con, "SELECT start, stop FROM chimera_windows ORDER BY start")
  snap <- chimera_json(con)

  # keep 41..100: what was cut goes, the rest moves 40 bp left
  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 40L, 60L)
  expect_equal(ch()$junctions, "30")
  expect_equal(win()$start, 11L)
  expect_equal(jsonlite::fromJSON(ch()$segments)$start, 11L)

  # the reverse shift puts the kept rows back in full-length coordinates
  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, -40L, 100L)
  expect_equal(ch()$junctions, "70")
  expect_equal(win()$start, 51L)

  RiboPilot:::chimera_restore(con, "S1", 1L, 1L, snap)
  expect_equal(ch()$junctions, "30,70")
  expect_equal(win()$start, c(1L, 51L))
})

test_that("shifting away every junction clears the chimera label, and chimera_restore brings it back", {
  dir <- withr::local_tempdir()
  seq <- paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")
  con <- mk_trim_db(dir, seq, rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, junctions, reason)
    VALUES ('S1',1,1,0.9,'severe','30','1-40 C1; 41-90 C2; junction near 30')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  ch <- function() DBI::dbGetQuery(con, "SELECT label, score, junctions, reason FROM chimera")
  snap <- chimera_json(con)

  RiboPilot:::chimera_shift(con, "S1", 1L, 1L, 44L, 56L)
  expect_equal(ch()$label, "clean")
  expect_equal(ch()$score, 0)
  expect_match(ch()$reason, "junctions trimmed away", fixed = TRUE)
  expect_false(grepl("41-90", ch()$reason, fixed = TRUE))

  RiboPilot:::chimera_restore(con, "S1", 1L, 1L, snap)
  expect_equal(ch()$label, "severe")
  expect_equal(ch()$junctions, "30")
  expect_equal(nrow(DBI::dbGetQuery(con, "SELECT * FROM chimera_windows")), 2L)
})
