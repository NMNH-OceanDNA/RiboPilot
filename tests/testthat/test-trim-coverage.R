# Trimming and restoring a unit rewrite its coverage summary. median_depth is
# the column the scaffold table and triage actually read, and the depth track
# can be shorter than the sequence (a truncated `samtools depth -a`, or a
# re-assembly coverage was never re-run on).

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

test_that("a trim rewrites median_depth, not just mean/min/max", {
  dir <- withr::local_tempdir()
  dep <- c(rep(2, 50), rep(20, 50))
  con <- mk_trim_db(dir, strrep("ACGT", 25), dep)
  on.exit(DBI::dbDisconnect(con))

  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, from = 51, to = 100)

  cv <- cov_row(con)
  expect_equal(cv$median_depth, 20)
  expect_equal(cv$mean_depth, 20)
  expect_equal(cv$min_depth, 20L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 50L)
})

test_that("a restore puts the original median_depth back", {
  dir <- withr::local_tempdir()
  dep <- c(rep(2, 50), rep(20, 50))
  con <- mk_trim_db(dir, strrep("ACGT", 25), dep)
  on.exit(DBI::dbDisconnect(con))

  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, from = 51, to = 100)
  expect_equal(cov_row(con)$median_depth, 20)

  RiboPilot:::restore_assembly(con, "S1", 1L, 1L)

  cv <- cov_row(con)
  expect_equal(cv$median_depth, stats::median(dep))
  expect_equal(cv$mean_depth, mean(dep))
  expect_equal(cv$min_depth, 2L)
  expect_equal(cv$max_depth, 20L)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 100L)
})

test_that("a depth track shorter than the trim start leaves the coverage alone", {
  dir <- withr::local_tempdir()
  dep <- rep(7, 40)                       # 40 positions of depth for 100 bases
  con <- mk_trim_db(dir, strrep("ACGT", 25), dep)
  on.exit(DBI::dbDisconnect(con))
  before <- cov_row(con)

  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, from = 60, to = 100)

  cv <- cov_row(con)
  expect_false(anyNA(unlist(cv)))
  expect_equal(cv$mean_depth, before$mean_depth)
  expect_equal(cv$median_depth, before$median_depth)
  expect_equal(cv$depths, before$depths)
})

test_that("a trim overlapping the end of a short depth track keeps the covered part", {
  dir <- withr::local_tempdir()
  con <- mk_trim_db(dir, strrep("ACGT", 25), c(rep(1, 20), rep(9, 20)))
  on.exit(DBI::dbDisconnect(con))

  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, from = 21, to = 100)

  cv <- cov_row(con)
  expect_equal(cv$median_depth, 9)
  expect_equal(length(strsplit(cv$depths, ",")[[1]]), 20L)
})

test_that("a trim moves the chimera junctions and windows, and a restore moves them back", {
  dir <- withr::local_tempdir()
  seq <- paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")
  con <- mk_trim_db(dir, seq, rep(20, 100))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE coverage SET clips = '0,1,2'")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1',1,1,'severe','30,70')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage, lineage_taxids)
    VALUES ('S1',1,1,1,40,99,'K;P;C1;o;f;g;s','1;2;3;;;;'), ('S1',1,1,51,90,99,'K;P;C2;o;f;g;s','1;2;4;;;;')")
  ch <- function() DBI::dbGetQuery(con, "SELECT junctions, segments FROM chimera")
  win <- function() DBI::dbGetQuery(con, "SELECT start, stop FROM chimera_windows ORDER BY start")

  RiboPilot:::trim_assembly(con, "S1", 1L, 1L, from = 41, to = 100)
  expect_equal(ch()$junctions, "30")
  expect_equal(win()$start, 11L)
  expect_equal(jsonlite::fromJSON(ch()$segments)$start, 11L)
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT clips FROM coverage")$clips))

  RiboPilot:::restore_assembly(con, "S1", 1L, 1L)
  expect_equal(ch()$junctions, "70")
  expect_equal(win()$start, 51L)
})
