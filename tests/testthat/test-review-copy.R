# Review copy lifecycle: first touch, invalidate, drop (SPEC S5.1-S5.4).

RC_SEQ <- strrep("ACGTTGCA", 25)

# Contigs S1.1.1 and S1.1.2 (200 bp), each with rows in every original a copy
# reads. Only columns a version 1 database has, so the upgrade test can use it.
rc_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:2) {
    p <- list(s)
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore,
        blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
      VALUES ('S1', 1, ?, 'linear', 200, ?, 0, 'A1', 'Danio rerio', 99.5, 88, 1e-50, 'core_nt', 7)",
      params = list(s, RC_SEQ))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, time_stamp)
      VALUES ('S1', 1, ?, 2, 0, 3, '5.8S,ITS2', 1, 2, 'truncated', ?, 9)",
      params = list(s, c("itsx_x", NA)[s]))
    for (f in list(list("18S", 1L, 60L, 1L, NA), list("ITS1", 61L, 100L, NA, NA),
                   list("28S", 101L, 200L, NA, 1L))) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length,
          strand, model_coverage, partial_start, partial_stop, tool, status, flags, edited, notes, time_stamp)
        VALUES ('S1', 1, ?, ?, 'rRNA', ?, ?, ?, '+', 95.5, ?, ?, 'barrnap', 'ok', 'f', 0, 'n', 9)",
        params = list(s, f[[1]], f[[2]], f[[3]], f[[3]] - f[[2]] + 1L, f[[4]], f[[5]]))
    }
    for (r in 1:2) {
      DBI::dbExecute(con, "INSERT INTO blast_hits (ID, path, scaffold, hit_rank, hsp, accession, species,
          pident, qcovs, evalue, bitscore, aln_len, qstart, qend, sstart, send, slen, time_stamp)
        VALUES ('S1', 1, ?, ?, 1, ?, 'Danio rerio', 99, 90, 1e-40, 300, 180, 11, 190, 1, 180, 1800, 9)",
        params = list(s, r, paste0("H", r)))
    }
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, accession, species, pident, qcovs,
        time_stamp) VALUES ('S1', 1, ?, 'H1', 'Danio rerio', 99, 90, 9)", params = p)
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version,
        variant_class, variant_group, variant_identity, time_stamp)
      VALUES ('S1', 1, ?, 'on_target', 3, 'major', 'g1', 99.1, 9)", params = p)
    DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, reason, time_stamp)
      VALUES ('S1', 1, ?, 0.1, 'clean', 'none', 9)", params = p)
    for (w in c(1L, 101L)) {
      DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, accession, pident)
        VALUES ('S1', 1, ?, ?, ?, 'H1', 98)", params = list(s, w, w + 99L))
    }
    DBI::dbExecute(con, "INSERT INTO blast_ref_alignment (ID, path, scaffold, accession, element, pid,
        aln_cov, aln_len, time_stamp) VALUES ('S1', 1, ?, 'R1', '18S', 98, 100, 60, 9)", params = p)
  }
}

# One row in each mirror, as a discarded recheck or an edit would leave them.
RC_NT <- length(RiboPilot:::REVIEW_TABLES)

rc_mirror_rows <- function(con, s) {
  p <- list(s)
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, ?, '18S')",
    params = p)
  DBI::dbExecute(con, "INSERT INTO review_blast_hits (ID, path, scaffold, hit_rank) VALUES ('S1', 1, ?, 9)",
    params = p)
  for (t in c("review_blast_screen", "review_screen", "review_chimera")) {
    DBI::dbExecute(con, sprintf("INSERT INTO %s (ID, path, scaffold) VALUES ('S1', 1, ?)", t), params = p)
  }
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start) VALUES ('S1', 1, ?, 7)",
    params = p)
  DBI::dbExecute(con, "INSERT INTO review_ref_alignment (ID, path, scaffold, accession, element)
    VALUES ('S1', 1, ?, 'R9', '28S')", params = p)
  DBI::dbExecute(con, "INSERT INTO review_feature_blast (ID, path, scaffold, element) VALUES ('S1', 1, ?, 'ITS1')",
    params = p)
  DBI::dbExecute(con, "INSERT INTO review_feature_blast_hits (ID, path, scaffold, element, hit_rank)
    VALUES ('S1', 1, ?, 'ITS1', 1)", params = p)
}

# Rows of a table, by column name, in a stable order.
rc_rows <- function(con, tbl, cols = DBI::dbListFields(con, tbl), path = 1L, scaffold = NULL) {
  ord <- paste(cols, collapse = ", ")
  if (is.null(scaffold)) {
    return(DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", ord, tbl, ord)))
  }
  DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s WHERE ID = 'S1' AND path = ? AND scaffold = ? ORDER BY %s",
    ord, tbl, ord), params = list(as.integer(path), as.integer(scaffold)))
}

rc_counts <- function(con, scaffold, path = 1L) {
  vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = ? AND scaffold = ?", t),
    params = list(as.integer(path), as.integer(scaffold)))$n, integer(1))
}

rc_originals <- function(con) {
  lapply(stats::setNames(nm = c("assemblies", "annotate", names(RiboPilot:::REVIEW_MIRRORS))),
    function(t) rc_rows(con, t))
}

rc_unit <- function(con, scaffold) {
  DBI::dbGetQuery(con, "SELECT sequence, length, trim_start, trim_end, edited, check_switch, check_note,
      CAST(rev AS TEXT) AS rev, checked_at, n_features, missing, chimeric, n_flags, curate_flags,
      itsx_warnings, itsx_origin, itsx_scores, blast_accession, blast_species, blast_pident, blast_qcovs,
      blast_evalue, blast_source, time_stamp
    FROM review_units WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(as.integer(scaffold)))
}

rc_now_ms <- function() as.numeric(Sys.time()) * 1000

# Copy, invalidate, and drop on a seeded database (new_db or an upgraded one).
rc_lifecycle <- function(con) {
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M', itsx_scores = 'M:50'")
  orig <- rc_originals(con)
  mirrors <- RiboPilot:::REVIEW_MIRRORS

  # first touch: the unit's originals, the view's check state, a fresh rev
  expect_false(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  t0 <- rc_now_ms()
  expect_true(RiboPilot:::review_copy(con, "S1", 1L, 1L))
  t1 <- rc_now_ms()
  expect_true(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  u <- rc_unit(con, 1)
  expect_equal(u[, c("sequence", "length", "trim_start", "trim_end", "edited", "check_switch",
      "check_note", "checked_at", "n_features", "missing", "chimeric", "n_flags", "curate_flags",
      "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession", "blast_species", "blast_pident",
      "blast_qcovs", "blast_evalue", "blast_source")],
    data.frame(sequence = RC_SEQ, length = 200L, trim_start = 1L, trim_end = 200L, edited = 0L,
      check_switch = 2L, check_note = NA_character_, checked_at = NA_integer_, n_features = 3L,
      missing = "5.8S,ITS2", chimeric = 1L, n_flags = 2L, curate_flags = "truncated",
      itsx_warnings = "itsx_x", itsx_origin = "M", itsx_scores = "M:50", blast_accession = "A1",
      blast_species = "Danio rerio", blast_pident = 99.5, blast_qcovs = 88, blast_evalue = 1e-50,
      blast_source = "core_nt"))
  expect_match(u$rev, "^[0-9]{13}$")
  expect_gte(as.numeric(u$rev), floor(t0))
  expect_lte(as.numeric(u$rev), ceiling(t1))
  expect_gte(u$time_stamp, floor(t0 / 1000))
  expect_lte(u$time_stamp, ceiling(t1 / 1000))
  for (o in names(mirrors)) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, mirrors[[o]], cols, scaffold = 1), rc_rows(con, o, cols, scaffold = 1), info = o)
  }
  expect_equal(nrow(rc_rows(con, "review_annotations", scaffold = 1)), 3L)
  expect_true(all(is.na(rc_rows(con, "review_annotations", c("partial_start_user", "partial_stop_user")))))
  expect_equal(unname(rc_counts(con, 2)), rep(0L, RC_NT))
  expect_identical(rc_originals(con), orig)

  # a second touch changes nothing
  expect_false(RiboPilot:::review_copy(con, "S1", 1L, 1L))
  expect_equal(rc_unit(con, 1), u)

  # the user's partial wins over curation's on the review layer
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 0 WHERE scaffold = 1 AND element = '18S'")
  eff <- DBI::dbGetQuery(con, paste("SELECT element,", RiboPilot:::REVIEW_EFF_PARTIALS,
    "FROM review_annotations_v WHERE scaffold = 1 ORDER BY start"))
  expect_equal(eff$partial_start, c(0L, NA, NA))
  expect_equal(eff$partial_stop, c(NA, NA, 1L))

  # rows a discarded recheck left for a contig without a copy are replaced;
  # done without a stored ITSx value starts as not run
  rc_mirror_rows(con, 2)
  expect_true(RiboPilot:::review_copy(con, "S1", 1L, 2L))
  for (o in names(mirrors)) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, mirrors[[o]], cols, scaffold = 2), rc_rows(con, o, cols, scaffold = 2), info = o)
  }
  expect_equal(rc_unit(con, 2)$check_switch, 1L)
  expect_true(is.na(rc_unit(con, 2)$itsx_warnings))

  # invalidate after a finished recheck and an edit that removed ITS1
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, check_note = 'old', checked_at = 5,
    time_stamp = 1 WHERE scaffold = 1")
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE scaffold = 1 AND element = 'ITS1'")
  other <- lapply(stats::setNames(nm = RiboPilot:::REVIEW_TABLES[-1]), function(t) rc_rows(con, t, scaffold = 2))
  other_unit <- rc_unit(con, 2)
  rev0 <- u$rev
  t0 <- rc_now_ms()
  rev1 <- RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  t1 <- rc_now_ms()
  v <- rc_unit(con, 1)
  expect_equal(rev1, v$rev)
  expect_equal(DBI::dbGetQuery(con, "SELECT CAST(rev - 1 AS TEXT) AS r FROM review_units WHERE scaffold = 1")$r,
    rev0)
  expect_equal(v[, c("check_switch", "check_note", "checked_at", "n_features", "missing", "n_flags",
      "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores")],
    data.frame(check_switch = 1L, check_note = NA_character_, checked_at = NA_integer_, n_features = 2L,
      missing = "5.8S,ITS1,ITS2", n_flags = NA_integer_, curate_flags = NA_character_,
      itsx_warnings = NA_character_, itsx_origin = NA_character_, itsx_scores = NA_character_))
  # kept as copied: sequence, window, the ITSx chimera flag, and the BLAST summary
  keep <- c("sequence", "length", "trim_start", "trim_end", "edited", "chimeric", "blast_accession",
    "blast_species", "blast_pident", "blast_qcovs", "blast_evalue", "blast_source")
  expect_equal(v[, keep], u[, keep])
  expect_gte(v$time_stamp, floor(t0 / 1000))
  expect_lte(v$time_stamp, ceiling(t1 / 1000))
  ra <- rc_rows(con, "review_annotations", c("element", "start", "stop", "model_coverage", "tool", "status",
    "flags", "partial_start", "partial_stop", "partial_start_user"), scaffold = 1)
  expect_equal(ra$element, c("18S", "28S"))
  expect_equal(ra$start, c(1L, 101L))
  expect_equal(ra$model_coverage, c(95.5, 95.5))
  expect_equal(ra$tool, c("barrnap", "barrnap"))
  expect_true(all(is.na(ra[, c("status", "flags", "partial_start", "partial_stop")])))
  expect_equal(ra$partial_start_user, c(0L, NA))
  expect_equal(unname(rc_counts(con, 1)), c(1L, 2L, rep(0L, RC_NT - 2L)))
  # the other contig's copy and every original are untouched
  expect_identical(lapply(stats::setNames(nm = RiboPilot:::REVIEW_TABLES[-1]), function(t)
    rc_rows(con, t, scaffold = 2)), other)
  expect_equal(rc_unit(con, 2), other_unit)
  expect_identical(rc_originals(con), orig)
  rev2 <- RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT CAST(rev - 2 AS TEXT) AS r FROM review_units
    WHERE scaffold = 1 AND rev = CAST(? AS INTEGER)", params = list(rev2))$r, rev0)

  # drop: the contig reads its originals again; the other copy stays
  RiboPilot:::review_drop(con, "S1", 1L, 1L)
  expect_false(RiboPilot:::review_has_copy(con, "S1", 1L, 1L))
  expect_equal(unname(rc_counts(con, 1)), rep(0L, RC_NT))
  expect_equal(unname(rc_counts(con, 2)), c(1L, 3L, 2L, 1L, 1L, 1L, 2L, 1L, 0L, 0L))
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1"), data.frame(has_copy = 0L, length = 200L, check_switch = 2L))
  expect_identical(rc_originals(con), orig)
  expect_no_error(RiboPilot:::review_drop(con, "S1", 1L, 1L))
}

test_that("with_savepoint() commits, rolls back, and nests inside an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  ws <- RiboPilot:::with_savepoint
  DBI::dbExecute(con, "CREATE TABLE t (x INTEGER)")
  xs <- function() DBI::dbGetQuery(con, "SELECT x FROM t ORDER BY x")$x

  expect_equal(ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (1)")
    "value"
  }), "value")
  expect_error(ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (2)")
    stop("boom")
  }), "boom")
  expect_equal(xs(), 1L)
  expect_false(RSQLite::sqliteIsTransacting(con))

  # inside DBI::dbWithTransaction: a failed block undoes only its own writes
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (3)")
    expect_error(ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (4)")
      stop("inner")
    }), "inner")
    ws(con, ws(con, DBI::dbExecute(con, "INSERT INTO t VALUES (5)")))
    expect_true(RSQLite::sqliteIsTransacting(con))
  })
  expect_equal(xs(), c(1L, 3L, 5L))

  # a nested failure caught inside an outer block
  ws(con, {
    DBI::dbExecute(con, "INSERT INTO t VALUES (6)")
    try(ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (7)")
      stop("nested")
    }), silent = TRUE)
  })
  expect_equal(xs(), c(1L, 3L, 5L, 6L))

  # leaving the block early rolls back like an error
  f <- function() {
    ws(con, {
      DBI::dbExecute(con, "INSERT INTO t VALUES (8)")
      return("early")
    })
    "late"
  }
  expect_equal(f(), "early")
  expect_equal(xs(), c(1L, 3L, 5L, 6L))
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("with_savepoint() takes the write lock before its first statement", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con), add = TRUE)
  other <- DBI::dbConnect(RSQLite::SQLite(), con@dbname)
  on.exit(DBI::dbDisconnect(other), add = TRUE)
  DBI::dbExecute(other, "PRAGMA busy_timeout = 0")
  DBI::dbExecute(con, "CREATE TABLE t (x INTEGER)")
  got <- RiboPilot:::with_savepoint(con, {
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM t")
    tryCatch({
      DBI::dbExecute(other, "INSERT INTO t VALUES (1)")
      "written"
    }, error = conditionMessage)
  })
  expect_match(got, "locked")
  expect_equal(DBI::dbExecute(other, "INSERT INTO t VALUES (1)"), 1)
})

test_that("review_new_rev() is the time in milliseconds with every digit", {
  withr::local_options(digits = 3)
  t0 <- rc_now_ms()
  r <- RiboPilot:::review_new_rev()
  t1 <- rc_now_ms()
  expect_type(r, "character")
  expect_match(r, "^[0-9]{13}$")
  expect_gte(as.numeric(r), floor(t0))
  expect_lte(as.numeric(r), ceiling(t1))
})

test_that("review_counts() counts the copy's features as reconcile_one does", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 1L),
    list(n_features = 0L, missing = "18S,5.8S,28S,ITS1,ITS2"))
  for (e in c("28S", "ITS1", "18S", "5.8S", "ITS2")) {
    DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, 1, ?)",
      params = list(e))
  }
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element) VALUES ('S1', 1, 2, '18S')")
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 1L), list(n_features = 5L, missing = ""))
  expect_equal(RiboPilot:::review_counts(con, "S1", 1L, 2L),
    list(n_features = 1L, missing = "5.8S,28S,ITS1,ITS2"))
})

test_that("copy, invalidate, and drop on a new project", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  rc_lifecycle(con)
})

test_that("copy, invalidate, and drop on a version 1 database holding data, upgraded through every step", {
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  rc_seed(con)
  DBI::dbDisconnect(con)

  suppressMessages(RiboPilot::backwards_compatibility(dir))
  con <- RiboPilot:::open_project_db(db)
  on.exit(DBI::dbDisconnect(con))
  expect_equal(max(DBI::dbGetQuery(con, "SELECT schema_version FROM schema_history")$schema_version),
    RiboPilot:::current_schema_version())
  rc_lifecycle(con)
})

test_that("review_copy() copies by column name on a hand-built database", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length,
      strand, status, flags, partial_start, edited, time_stamp)
    VALUES ('S1', 2, 1, '18S', 'rRNA', 5, 40, 36, '-', 'ok', 'f', 1, 0, 9)")
  DBI::dbExecute(con, "INSERT INTO blast_screen (ID, path, scaffold, accession, species, pident, time_stamp)
    VALUES ('S1', 2, 1, 'H1', 'Danio rerio', 97.5, 9)")
  RiboPilot:::review_copy(con, "S1", 2L, 1L)
  for (o in c("annotations", "blast_screen")) {
    cols <- DBI::dbListFields(con, o)
    expect_equal(rc_rows(con, RiboPilot:::REVIEW_MIRRORS[[o]], cols, path = 2, scaffold = 1),
      rc_rows(con, o, cols, path = 2, scaffold = 1), info = o)
  }
  expect_equal(DBI::dbGetQuery(con, "SELECT path, scaffold, sequence, length, trim_end, check_switch,
      n_features, missing, blast_accession, blast_qcovs FROM review_units"),
    data.frame(path = 2L, scaffold = 1L, sequence = "ACGT", length = 5980L, trim_end = 5980L,
      check_switch = 1L, n_features = 3L, missing = "", blast_accession = "A1", blast_qcovs = 88))

  RiboPilot:::review_invalidate(con, "S1", 2L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT n_features, missing FROM review_units"),
    data.frame(n_features = 1L, missing = "5.8S,28S,ITS1,ITS2"))
  expect_equal(unname(rc_counts(con, 1, path = 2)), c(1L, 1L, rep(0L, RC_NT - 2L)))
  RiboPilot:::review_drop(con, "S1", 2L, 1L)
  expect_equal(unname(rc_counts(con, 1, path = 2)), rep(0L, RC_NT))
})

test_that("review_copy() refuses a contig with no annotate or assemblies row and writes nothing", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  DBI::dbExecute(con, "DELETE FROM annotate WHERE scaffold = 2")
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 2L), "No sequence for contig S1.1.2.", fixed = TRUE)
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 9L), "No sequence for contig S1.1.9.", fixed = TRUE)
  expect_equal(unname(rc_counts(con, 2)), rep(0L, RC_NT))
  expect_error(RiboPilot:::review_invalidate(con, "S1", 1L, 1L), "No review copy for contig S1.1.1.",
    fixed = TRUE)
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("a failure part way through review_copy() leaves no copy", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  DBI::dbExecute(con, "CREATE TRIGGER rc_boom BEFORE INSERT ON review_ref_alignment
    BEGIN SELECT RAISE(ABORT, 'rc_boom'); END")
  expect_error(RiboPilot:::review_copy(con, "S1", 1L, 1L), "rc_boom")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, RC_NT))
  expect_false(RSQLite::sqliteIsTransacting(con))
})

test_that("review_copy, review_invalidate, and review_drop join an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  rc_seed(con)
  expect_error(DBI::dbWithTransaction(con, {
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    stop("outer failed")
  }), "outer failed")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, RC_NT))

  DBI::dbWithTransaction(con, {
    RiboPilot:::review_copy(con, "S1", 1L, 1L)
    RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
    RiboPilot:::review_copy(con, "S1", 1L, 2L)
    RiboPilot:::review_drop(con, "S1", 1L, 2L)
  })
  expect_equal(rc_unit(con, 1)$check_switch, 1L)
  expect_equal(unname(rc_counts(con, 1)), c(1L, 3L, rep(0L, RC_NT - 2L)))
  expect_equal(unname(rc_counts(con, 2)), rep(0L, RC_NT))
})

test_that("a re-assembly deletes the copies of changed and vanished contigs", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  out <- withr::local_tempdir()
  sq <- lapply(1:4, function(k) withr::with_seed(k, paste(sample(c("A", "C", "G", "T"), 1200, TRUE),
    collapse = "")))
  fa <- function(x) writeLines(unlist(lapply(names(x), function(k)
    c(sprintf(">S1.1.%s linear c%s", k, k), x[[k]]))), file.path(out, "S1_assembly_1.fasta"))
  fa(list("1" = sq[[1]], "2" = sq[[2]], "3" = sq[[3]]))
  RiboPilot:::load_assemblies(con, out, "S1")
  for (s in 1:3) {
    RiboPilot:::review_copy(con, "S1", 1L, s)
    rc_mirror_rows(con, s)
  }
  expect_equal(unname(rc_counts(con, 1)), rep(1L, RC_NT))

  # scaffold 1 comes back with another sequence, 2 unchanged, 3 gone
  fa(list("1" = sq[[4]], "2" = sq[[2]]))
  RiboPilot:::load_assemblies(con, out, "S1")
  expect_equal(unname(rc_counts(con, 1)), rep(0L, RC_NT))
  expect_equal(unname(rc_counts(con, 2)), rep(1L, RC_NT))
  expect_equal(unname(rc_counts(con, 3)), rep(0L, RC_NT))

  RiboPilot:::reset_unit_evidence(con, "S1", 1L, 2L)
  expect_equal(unname(rc_counts(con, 2)), rep(0L, RC_NT))
})
