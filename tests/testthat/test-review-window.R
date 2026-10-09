# Coverage window and mapped BLAST hits of a review copy (SPEC S4.5).

# Kept contigs S1.1.<scaffold> of `len` bp, done in the main pipeline.
rw_units <- function(con, scaffolds, len) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in scaffolds) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', ?, ?, 0)", params = list(s, len, strrep("ACGTA", len / 5)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, '')", params = list(s))
  }
}

# Make the contig's copy and set its kept window, in original coordinates.
rw_window <- function(con, scaffold, from, to) {
  RiboPilot:::review_copy(con, "S1", 1L, scaffold)
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = ?, trim_end = ?, length = ?
    WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(from, to, to - from + 1L, scaffold))
}

rw_track <- function(x) if (is.null(x)) NA_character_ else paste(x, collapse = ",")

test_that("review_coverage() cuts the original tracks to the copy's window", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:4, 10L)
  cov <- function(s, depths, clips) DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips)
    VALUES ('S1', 1, ?, ?, ?)", params = list(s, rw_track(depths), rw_track(clips)))
  cov(1, 11:20, 0:9)
  cov(2, 31:35, 40:49)
  cov(3, 51:60, NULL)
  rc <- function(s) RiboPilot:::review_coverage(con, "S1", 1L, s)

  rw_window(con, 2, 3L, 7L)
  rw_window(con, 3, 2L, 9L)
  # no copy: the whole original tracks
  expect_equal(rc(1), list(depths = as.numeric(11:20), clips = as.numeric(0:9)))
  # no coverage row
  expect_equal(rc(4), list(depths = NULL, clips = NULL))

  rw_window(con, 1, 3L, 7L)
  expect_equal(rc(1), list(depths = as.numeric(13:17), clips = as.numeric(2:6)))
  # a depth track shorter than the contig is cut, never padded
  expect_equal(rc(2), list(depths = as.numeric(33:35), clips = as.numeric(42:46)))
  # clips lost by an old trim
  expect_equal(rc(3), list(depths = as.numeric(52:59), clips = NULL))
  # a window that starts past the end of a track
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = 6, trim_end = 9 WHERE scaffold = 2")
  expect_equal(rc(2), list(depths = NULL, clips = as.numeric(45:48)))
  rw_window(con, 4, 1L, 5L)
  expect_equal(rc(4), list(depths = NULL, clips = NULL))

  expect_equal(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE scaffold = 1"),
    data.frame(depths = rw_track(11:20), clips = rw_track(0:9)))
})

test_that("review_coverage() reads a coverage table without a clips column", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:"); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE coverage (ID TEXT, path INT, scaffold INT, mean_depth REAL, depths TEXT)")
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  DBI::dbExecute(con, "INSERT INTO coverage VALUES ('S1', 1, 1, 3, '1,2,3,4,5')")
  expect_equal(RiboPilot:::review_coverage(con, "S1", 1L, 1L), list(depths = as.numeric(1:5), clips = NULL))
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, trim_start, trim_end, rev)
    VALUES ('S1', 1, 1, 2, 4, 1)")
  expect_equal(RiboPilot:::review_coverage(con, "S1", 1L, 1L), list(depths = as.numeric(2:4), clips = NULL))
})

test_that("review_coverage_summary() rounds as write_coverage_summary() stores", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  dep <- c(2, 9, 4, 4, 1, 3)
  expect_identical(RiboPilot:::review_coverage_summary(dep),
    list(mean_depth = 3.83, median_depth = 3.5, min_depth = 1L, max_depth = 9L))
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold) VALUES ('S1', 1, 1)")
  RiboPilot:::write_coverage_summary(con, "S1", 1L, 1L, dep)
  expect_equal(as.list(DBI::dbGetQuery(con, "SELECT mean_depth, median_depth, min_depth, max_depth FROM coverage")),
    RiboPilot:::review_coverage_summary(dep))
  expect_null(RiboPilot:::review_coverage_summary(NULL))
  expect_null(RiboPilot:::review_coverage_summary(numeric(0)))
})

# HSPs of contig S1.1.1 (300 bp). qcovs is stored as BLAST computes it on the
# whole contig: the union of a subject's HSPs over 300, rounded. Subject 7 has
# no query coordinates (an old hits file).
RW_HITS <- data.frame(
  hit_rank = c(1L, 1L, 2L, 3L, 3L, 3L, 4L, 4L, 5L, 6L, 6L, 7L),
  hsp = c(1L, 2L, 1L, 1L, 2L, 3L, 1L, 2L, 1L, 1L, 2L, 1L),
  qstart = c(1L, 61L, 1L, 241L, 251L, 51L, 230L, 200L, 250L, 40L, 30L, NA),
  qend = c(300L, 90L, 50L, 300L, 280L, 65L, 260L, 250L, 300L, 60L, 51L, NA),
  sstart = c(1L, 500L, 1L, 1L, 30L, 1000L, 9L, 77L, 1L, 5L, 2L, NA),
  send = c(300L, 470L, 50L, 60L, 59L, 1014L, 39L, 127L, 51L, 25L, 23L, NA),
  aln_len = c(300L, 31L, 50L, 60L, 30L, 15L, 31L, 51L, 51L, 21L, 22L, 99L),
  qcovs = c(100, 100, 17, 25, 25, 25, 20, 20, 17, 10, 10, 0))

rw_hits <- function(con) {
  h <- cbind(data.frame(ID = "S1", path = 1L, scaffold = 1L), RW_HITS,
    accession = paste0("H", RW_HITS$hit_rank), species = "Danio rerio", hit_class = "rdna",
    pident = 90 + seq_len(nrow(RW_HITS)) / 10, bitscore = 1000 - seq_len(nrow(RW_HITS)),
    evalue = 1e-30, slen = 2000L, time_stamp = 9L)
  DBI::dbAppendTable(con, "blast_hits", h)
  # another contig's hit, outside the window of its own copy
  DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession, pident, qcovs,
      aln_len, qstart, qend, sstart, send) VALUES ('S1', 1, 2, 1, 1, 'X1', 99, 33, 100, 1, 100, 1, 100)")
}

rw_orig_hits <- function(con, scaffold) {
  DBI::dbGetQuery(con, sprintf("SELECT %s FROM blast_hits WHERE ID = 'S1' AND path = 1 AND scaffold = ?
    ORDER BY hit_rank, hsp", paste(RiboPilot:::ddl_cols(RiboPilot:::BLAST_HITS_DDL), collapse = ", ")),
    params = list(scaffold))
}

test_that("review_map_hits() refuses a contig without a copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1, 300L)
  rw_hits(con)
  expect_error(RiboPilot:::review_map_hits(con, "S1", 1L, 1L), "No review copy for contig S1.1.1.", fixed = TRUE)
})

test_that("review_map_hits() on an untrimmed copy returns the original hits", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:2, 300L)
  rw_hits(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  rw_window(con, 2, 101L, 300L)
  orig <- rw_orig_hits(con, 1)
  expect_equal(RiboPilot:::review_map_hits(con, "S1", 1L, 1L), orig[!is.na(orig$qstart), ], ignore_attr = TRUE)
})

test_that("review_map_hits() clips and shifts hits to the kept window", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rw_units(con, 1:3, 300L)
  rw_hits(con)
  rw_window(con, 2, 101L, 300L)
  rw_window(con, 1, 51L, 250L)
  rw_window(con, 3, 1L, 300L)
  orig <- rw_orig_hits(con, 1)
  m <- RiboPilot:::review_map_hits(con, "S1", 1L, 1L)

  expect_named(m, DBI::dbListFields(con, "review_blast_hits"))
  # subjects 2 and 7 fall outside; the others keep their stored rank and hsp
  expect_equal(m[, c("hit_rank", "hsp", "qstart", "qend", "sstart", "send", "aln_len", "qcovs")],
    data.frame(hit_rank = c(1L, 1L, 3L, 3L, 4L, 4L, 5L, 6L, 6L), hsp = c(1L, 2L, 1L, 3L, 1L, 2L, 1L, 1L, 2L),
      qstart = c(1L, 11L, 191L, 1L, 180L, 150L, 200L, 1L, 1L), qend = c(200L, 40L, 200L, 15L, 200L, 200L, 200L, 10L, 1L),
      sstart = c(NA, 500L, NA, 1000L, NA, 77L, NA, NA, NA), send = c(NA, 470L, NA, 1014L, NA, 127L, NA, NA, NA),
      aln_len = c(200L, 31L, 10L, 15L, 21L, 51L, 1L, 10L, 1L), qcovs = c(100, 100, 13, 13, 26, 26, 1, 5, 5)))
  # everything else as stored
  kept <- match(paste(m$hit_rank, m$hsp), paste(orig$hit_rank, orig$hsp))
  rest <- setdiff(names(m), c("qstart", "qend", "sstart", "send", "aln_len", "qcovs"))
  expect_equal(m[, rest], orig[kept, rest], ignore_attr = TRUE)

  # appendable to the mirror as the recheck writes it
  DBI::dbExecute(con, "DELETE FROM review_blast_hits WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  DBI::dbAppendTable(con, "review_blast_hits", m)
  expect_equal(DBI::dbGetQuery(con, sprintf("SELECT %s FROM review_blast_hits WHERE scaffold = 1 ORDER BY hit_rank, hsp",
    paste(names(m), collapse = ", "))), m)

  # nothing overlaps the window, or the original has no hits
  for (s in 2:3) {
    z <- RiboPilot:::review_map_hits(con, "S1", 1L, s)
    expect_equal(nrow(z), 0L)
    expect_named(z, names(m))
  }
  expect_equal(rw_orig_hits(con, 1), orig)
})
