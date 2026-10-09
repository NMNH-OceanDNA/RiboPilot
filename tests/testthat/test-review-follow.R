# The main pipeline's reruns and option changes follow a contig's review copy (SPEC S5.3).

# Five annotated, done 100 bp contigs of S1, each with one feature and a
# current screen row.
rf_units <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:5) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence)
      VALUES ('S1', 1, ?, 'linear', 100, ?)", params = list(s, strrep("ACGT", 25)))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        annotate_opts, curate_opts, n_features, missing, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, 'default', 'default', 1, '', '')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand, status)
      VALUES ('S1', 1, ?, '18S', 1, 90, '+', 'ok')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version)
      VALUES ('S1', 1, ?, 'on_target', ?)", params = list(s, RiboPilot:::SCREEN_VERSION))
  }
}

test_that("review_follow_original drops a State copy, invalidates an edited one, and spares a locked one", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rf_units(con)
  for (s in 2:5) RiboPilot:::review_copy(con, "S1", 1L, s)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE scaffold IN (3, 4)")
  DBI::dbExecute(con, "UPDATE review_annotations SET start = 11 WHERE scaffold IN (3, 4)")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold IN (4, 5)")
  follow <- function(s) RiboPilot:::review_follow_original(con, "S1", 1L, s)
  n_rows <- function(s) vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = ?", t),
    params = list(s))$n, integer(1))
  snap <- function(s) list(n_rows(s),
    DBI::dbGetQuery(con, "SELECT edited, check_switch, CAST(rev AS TEXT) AS rev, time_stamp
      FROM review_units WHERE scaffold = ?", params = list(s)),
    DBI::dbGetQuery(con, "SELECT start, status FROM review_annotations WHERE scaffold = ?",
      params = list(s)))
  rev3 <- DBI::dbGetQuery(con, "SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE scaffold = 3")$rev
  before <- lapply(4:5, snap)
  expect_equal(unname(n_rows(3)), c(1L, 1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 0L))

  expect_invisible(follow(1))
  expect_equal(follow(1), "none")
  expect_equal(sum(n_rows(1)), 0L)

  expect_equal(follow(2), "dropped")
  expect_equal(sum(n_rows(2)), 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy FROM review_units_v WHERE scaffold = 2")$has_copy, 0L)

  expect_equal(follow(3), "invalidated")
  expect_equal(DBI::dbGetQuery(con, "SELECT edited, check_switch, rev = CAST(? AS INTEGER) + 1 AS bumped
    FROM review_units WHERE scaffold = 3", params = list(rev3)),
    data.frame(edited = 1L, check_switch = 1L, bumped = 1L))
  a3 <- DBI::dbGetQuery(con, "SELECT start, status FROM review_annotations WHERE scaffold = 3")
  expect_equal(a3$start, 11L)
  expect_true(is.na(a3$status))
  expect_equal(n_rows(3)[["review_screen"]], 0L)

  expect_equal(follow(4), "locked")
  expect_equal(follow(5), "locked")
  expect_equal(lapply(4:5, snap), before)

  # inside an open transaction it nests as a savepoint
  DBI::dbBegin(con)
  expect_equal(follow(3), "invalidated")
  DBI::dbRollback(con)
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("a contig without a copy is answered from review_units alone", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REVIEW_UNITS_DDL)
  expect_equal(RiboPilot:::review_follow_original(con, "S1", 1L, 1L), "none")
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("each per-contig pipeline script follows the copy of the contig it wrote", {
  bin <- system.file("nextflow", "bin", package = "RiboPilot")
  for (f in c("rp_reconcile.R", "rp_blast.R", "rp_searched.R", "rp_curate.R", "rp_chimera.R",
              "rp_screen.R")) {
    txt <- readLines(file.path(bin, f), warn = FALSE)
    expect_equal(sum(grepl("RiboPilot:::review_follow_original(con, a[2], ", txt, fixed = TRUE)), 1L,
      info = f)
  }
  chim <- readLines(file.path(bin, "rp_chimera.R"), warn = FALSE)
  expect_true(any(grepl("if (!is.null(r)) RiboPilot:::review_follow_original", chim, fixed = TRUE)))
})

test_that("rp_screen.R follows only the contigs whose screen it rewrote, and rp_curate.R its contig", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not_installed("pkgload")
  src <- RiboPilot:::ribopilot_src_default()
  skip_if(!nzchar(src), "not a source checkout")
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rf_units(con)
  for (s in 1:3) RiboPilot:::review_copy(con, "S1", 1L, s)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE scaffold = 3")
  # only contigs 2 and 3 lack a current screen row, so the screen rewrites only theirs
  DBI::dbExecute(con, "UPDATE screen SET screen_version = NULL WHERE scaffold IN (2, 3)")
  run <- function(script, args = character(0)) processx::run("Rscript",
    c(file.path(system.file("nextflow", "bin", package = "RiboPilot"), script),
      DBI::dbGetInfo(con)$dbname, "S1", args),
    env = c("current", RIBOPILOT_SRC = src), error_on_status = FALSE)
  copies <- function() DBI::dbGetQuery(con,
    "SELECT scaffold, check_switch FROM review_units ORDER BY scaffold")

  res <- run("rp_screen.R")
  expect_equal(res$status, 0L, info = res$stderr)
  expect_match(res$stdout, "screen: S1 - 5 unit(s)", fixed = TRUE)
  expect_equal(copies(), data.frame(scaffold = c(1L, 3L), check_switch = c(2L, 1L)))

  res <- run("rp_curate.R", c("1", "1"))
  expect_equal(res$status, 0L, info = res$stderr)
  expect_equal(copies(), data.frame(scaffold = 3L, check_switch = 1L))
})
