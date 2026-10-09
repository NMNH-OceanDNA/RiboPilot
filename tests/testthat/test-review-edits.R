# Annotation Review edits on a contig's review copy: trim, boundary nudge,
# 5'/3' choice, Reset, and State (SPEC S5.1, S5.2).

RE_SEQ <- withr::with_seed(4, replicate(2, paste(sample(c("A", "C", "G", "T"), 100, TRUE), collapse = "")))

# S1.1.1 has its features on the plus strand, S1.1.2 on the minus strand.
RE_FEATS <- list(
  list(element = c("18S", "ITS1", "5.8S", "ITS2", "28S"), start = c(11L, 41L, 51L, 61L, 71L),
    stop = c(40L, 50L, 60L, 70L, 95L), strand = "+"),
  list(element = c("28S", "ITS2", "5.8S", "ITS1", "18S"), start = c(6L, 31L, 41L, 51L, 61L),
    stop = c(30L, 40L, 50L, 60L, 90L), strand = "-"))

# Two done 100 bp contigs; curation marked every feature's 3' end partial.
re_seed <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:2) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', 100, ?, 0)", params = list(s, RE_SEQ[s]))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features,
        missing, n_flags, curate_flags, itsx_warnings)
      VALUES ('S1', 1, ?, 2, 0, 5, '', 1, 'truncated', '')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths) VALUES ('S1', 1, ?, ?)",
      params = list(s, paste(1:100, collapse = ",")))
    f <- RE_FEATS[[s]]
    for (i in 1:5) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, length, strand,
          partial_stop, status, edited)
        VALUES ('S1', 1, ?, ?, ?, ?, ?, ?, 1, 'ok', 0)",
        params = list(s, f$element[i], f$start[i], f$stop[i], f$stop[i] - f$start[i] + 1L, f$strand))
    }
  }
}

re_ints <- function(d) {
  k <- setdiff(names(d), c("element", "sequence", "missing"))
  d[k] <- lapply(d[k], as.integer)
  d
}

re_feats <- function(con, s) {
  re_ints(DBI::dbGetQuery(con, "SELECT element, start, stop, length, edited, partial_start_user,
      partial_stop_user FROM review_annotations WHERE ID = 'S1' AND path = 1 AND scaffold = ? ORDER BY start",
    params = list(as.integer(s))))
}

re_unit <- function(con, s) {
  re_ints(DBI::dbGetQuery(con, "SELECT sequence, length, trim_start, trim_end, edited, check_switch,
      n_features, missing, n_flags FROM review_units WHERE ID = 'S1' AND path = 1 AND scaffold = ?",
    params = list(as.integer(s))))
}

re_rev <- function(con, s, minus = 0L) {
  DBI::dbGetQuery(con, "SELECT CAST(rev - ? AS TEXT) AS rev FROM review_units
    WHERE ID = 'S1' AND path = 1 AND scaffold = ?", params = list(as.integer(minus), as.integer(s)))$rev
}

re_rows <- function(con, s) {
  sum(vapply(RiboPilot:::REVIEW_TABLES, function(t) DBI::dbGetQuery(con, sprintf(
    "SELECT COUNT(*) AS n FROM %s WHERE ID = 'S1' AND path = 1 AND scaffold = ?", t),
    params = list(as.integer(s)))$n, integer(1)))
}

re_originals <- function(con) {
  lapply(stats::setNames(nm = c("assemblies", "annotate", "annotations", "coverage")), function(t) {
    cols <- paste(DBI::dbListFields(con, t), collapse = ", ")
    DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", cols, t, cols))
  })
}

re_feat_df <- function(element, start, stop, edited, ps, pe) {
  data.frame(element = element, start = start, stop = stop, length = stop - start + 1L,
    edited = edited, partial_start_user = ps, partial_stop_user = pe)
}

test_that("a trim cuts the copy from the original window and marks cut ends on the plus strand", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 80), 60L)
  expect_equal(re_unit(con, 1), data.frame(sequence = substr(RE_SEQ[1], 21, 80), length = 60L,
    trim_start = 21L, trim_end = 80L, edited = 1L, check_switch = 1L, n_features = 5L, missing = "",
    n_flags = NA_))
  # 18S is cut on the left (its 5' end); 28S on the right (its 3' end)
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 60L), rep(1L, 5), c(1L, NA_, NA_, NA_, NA_), c(NA_, NA_, NA_, NA_, 1L)))
  # curation's partial ends are cleared until the recheck; the cut ends stay
  eff <- DBI::dbGetQuery(con, paste("SELECT", RiboPilot:::REVIEW_EFF_PARTIALS, "FROM review_annotations_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1 ORDER BY start"))
  expect_equal(as.integer(eff$partial_start), c(1L, NA_, NA_, NA_, NA_))
  expect_equal(as.integer(eff$partial_stop), c(NA_, NA_, NA_, NA_, 1L))
  rev1 <- re_rev(con, 1)

  # a second trim works in the copy's coordinates: 5..50 of the copy is 25..70
  # of the original; ITS2 now ends exactly at the right end, and 28S is gone
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 1L, 5, 50), 46L)
  expect_equal(re_unit(con, 1), data.frame(sequence = substr(RE_SEQ[1], 25, 70), length = 46L,
    trim_start = 25L, trim_end = 70L, edited = 1L, check_switch = 1L, n_features = 4L, missing = "28S",
    n_flags = NA_))
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element[1:4], c(1L, 17L, 27L, 37L),
    c(16L, 26L, 36L, 46L), rep(1L, 4), c(1L, NA_, NA_, NA_), rep(NA_, 4)))
  expect_equal(re_rev(con, 1, 1L), rev1)
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, orig_length, length, trim_start, trim_end
      FROM review_units_v WHERE ID = 'S1' AND path = 1 AND scaffold = 1"),
    data.frame(has_copy = 1L, orig_length = 100L, length = 46L, trim_start = 25L, trim_end = 70L))

  # the other contig and every original are untouched
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("on the minus strand a left cut marks the 3' end and a right cut the 5' end", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  # 28S is cut on the left; 18S ends exactly at the new right end, so it is not partial
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 2L, 11, 90), 80L)
  expect_equal(re_feats(con, 2), re_feat_df(RE_FEATS[[2]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 80L), rep(1L, 5), rep(NA_, 5), c(1L, NA_, NA_, NA_, NA_)))

  # a 5'/3' choice the user made stays on an uncut feature; a cut end overrides it
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start_user = 0
    WHERE scaffold = 2 AND element IN ('ITS1', '18S')")
  expect_equal(RiboPilot:::review_trim(con, "S1", 1L, 2L, 1, 70), 70L)
  expect_equal(re_feats(con, 2), re_feat_df(RE_FEATS[[2]]$element, c(1L, 21L, 31L, 41L, 51L),
    c(20L, 30L, 40L, 50L, 70L), rep(1L, 5), c(NA_, NA_, NA_, 0L, 1L), c(1L, NA_, NA_, NA_, NA_)))
  expect_equal(re_unit(con, 2)[, c("sequence", "length", "trim_start", "trim_end")],
    data.frame(sequence = substr(RE_SEQ[2], 11, 80), length = 70L, trim_start = 11L, trim_end = 80L))
  expect_equal(re_rows(con, 1), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("a trim clamps to the copy's sequence and refuses bad ranges and missing sequences", {
  dir <- mk_app_project()
  app <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(app))
  # stored length 6120, sequence 'ACGT': the sequence decides
  expect_equal(RiboPilot:::review_trim(app, "S1", 1L, 1L, 2, 6000), 3L)
  expect_equal(RiboPilot:::review_trim(app, "S1", 1L, 1L, -5, 2), 2L)
  expect_equal(DBI::dbGetQuery(app, "SELECT sequence, length, trim_start, trim_end FROM review_units"),
    data.frame(sequence = "CG", length = 2L, trim_start = 2L, trim_end = 3L))

  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con), add = TRUE)
  re_seed(con)
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, ignore) VALUES ('S1', 1, 3, 100, 0)")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, n_features)
    VALUES ('S1', 1, 3, 2, 0)")
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 9L, 1, 10), "No sequence for contig S1.1.9.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 3L, 1, 10), "No sequence for contig S1.1.3.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 1L, 50, 40), "Invalid trim range.", fixed = TRUE)
  expect_error(RiboPilot:::review_trim(con, "S1", 1L, 1L, 150, 200), "Invalid trim range.", fixed = TRUE)
  expect_equal(vapply(c(1, 3, 9), function(s) re_rows(con, s), integer(1)), c(0L, 0L, 0L))
})

test_that("a boundary nudge writes the copy and marks its checks not run", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_
  f0 <- re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop, rep(0L, 5), rep(NA_, 5),
    rep(NA_, 5))

  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40)
  f <- f0; f[1, c("start", "length", "edited")] <- c(5L, 36L, 1L)
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_unit(con, 1), data.frame(sequence = RE_SEQ[1], length = 100L, trim_start = 1L,
    trim_end = 100L, edited = 1L, check_switch = 1L, n_features = 5L, missing = "", n_flags = NA_))
  rev1 <- re_rev(con, 1)

  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "28S", 71, 99)
  f[5, c("stop", "length", "edited")] <- c(99L, 29L, 1L)
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1, 1L), rev1)

  expect_error(RiboPilot:::review_set_bounds(con, "S1", 1L, 2L, "XX", 1, 5), "No feature XX on contig S1.1.2.",
    fixed = TRUE)
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("a 5'/3' choice writes the user column and keeps the checks", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  NA_ <- NA_integer_

  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "stop", 0)
  rev1 <- re_rev(con, 1)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "start", 1)
  expect_equal(re_feats(con, 1), re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop,
    c(1L, 0L, 0L, 0L, 1L), c(NA_, NA_, NA_, NA_, 1L), c(0L, NA_, NA_, NA_, NA_)))
  expect_equal(re_unit(con, 1)[, c("edited", "check_switch", "n_flags")],
    data.frame(edited = 1L, check_switch = 2L, n_flags = 1L))
  expect_equal(re_rev(con, 1), rev1)
  # the user's choice wins over curation's, which stays until the next recheck
  eff <- DBI::dbGetQuery(con, paste("SELECT", RiboPilot:::REVIEW_EFF_PARTIALS, "FROM review_annotations_v
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1 ORDER BY start"))
  expect_equal(as.integer(eff$partial_start), c(NA_, NA_, NA_, NA_, 1L))
  expect_equal(as.integer(eff$partial_stop), c(0L, 1L, 1L, 1L, 1L))

  expect_error(RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "middle", 1), "should be one of")
  expect_error(RiboPilot:::review_set_partial(con, "S1", 1L, 2L, "XX", "start", 1),
    "No feature XX on contig S1.1.2.", fixed = TRUE)
  expect_equal(re_rows(con, 2), 0L)
  expect_identical(re_originals(con), orig)
})

test_that("Reset writes the session's start back and invalidates only when a boundary moved", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  NA_ <- NA_integer_
  backup <- function() DBI::dbGetQuery(con, "SELECT element, start, stop, partial_start_user, partial_stop_user
    FROM review_annotations_v WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  f0 <- re_feat_df(RE_FEATS[[1]]$element, RE_FEATS[[1]]$start, RE_FEATS[[1]]$stop, rep(0L, 5), rep(NA_, 5),
    rep(NA_, 5))

  # Edit ticked on a contig without a copy, nothing changed: nothing to write
  b0 <- backup()
  expect_false(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b0))
  expect_equal(re_rows(con, 1), 0L)

  # a nudge and a toggle, then Reset
  RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 1)
  rev1 <- re_rev(con, 1)
  expect_true(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b0))
  f <- f0; f$edited[c(1, 5)] <- 1L
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1, 1L), rev1)
  expect_equal(re_unit(con, 1)$check_switch, 1L)

  # a session with only toggles; the backup may name a feature the copy no longer has
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, n_flags = 4 WHERE scaffold = 1")
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 1)
  b1 <- rbind(backup(), data.frame(element = "XX", start = 1L, stop = 2L, partial_start_user = NA,
    partial_stop_user = NA))
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "start", 0)
  RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "28S", "stop", 0)
  rev2 <- re_rev(con, 1)
  expect_false(RiboPilot:::review_reset_session(con, "S1", 1L, 1L, b1))
  f$partial_stop_user[5] <- 1L
  expect_equal(re_feats(con, 1), f)
  expect_equal(re_rev(con, 1), rev2)
  expect_equal(re_unit(con, 1)[, c("check_switch", "n_flags")], data.frame(check_switch = 2L, n_flags = 4L))
})

test_that("State queues or holds copies and refuses locked and waiting contigs", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  # 1 done; 2 locked; 3 running and 6 queued in the main pipeline; 4 never
  # annotated; 5 queued in the main pipeline but already has a copy
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold = 2")
  for (s in 3:6) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES ('S1', 1, ?, 4, 'ACGT', 0)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features)
      VALUES ('S1', 1, ?, ?, 0, ?)", params = list(s, c(4L, 2L, 1L, 1L)[s - 2], c(0L, NA, 0L, 0L)[s - 2]))
  }
  RiboPilot:::review_copy(con, "S1", 1L, 5L)
  rev5 <- re_rev(con, 5)
  u <- data.frame(ID = "S1", path = 1L, scaffold = 1:6)
  wait <- "press Update on the Assemble step first"
  sw <- function() DBI::dbGetQuery(con, "SELECT scaffold, check_switch FROM review_units ORDER BY scaffold")

  t0 <- as.integer(Sys.time())
  expect_equal(RiboPilot:::review_set_state(con, u, 1), c(NA, "unlock first", wait, wait, NA, wait))
  expect_equal(sw(), data.frame(scaffold = c(1L, 5L), check_switch = c(1L, 1L)))
  expect_equal(re_rev(con, 5, 1L), rev5)
  rev1 <- re_rev(con, 1); rev5 <- re_rev(con, 5)

  expect_equal(RiboPilot:::review_set_state(con, u, 0), c(NA, "unlock first", NA, NA, NA, NA))
  expect_equal(sw(), data.frame(scaffold = c(1L, 3:6), check_switch = rep(0L, 5)))
  expect_equal(c(re_rev(con, 1), re_rev(con, 5)), c(rev1, rev5))
  expect_true(all(DBI::dbGetQuery(con, "SELECT time_stamp FROM review_units")$time_stamp >= t0))
  expect_equal(re_rows(con, 2), 0L)

  # a locked contig's copy is left alone too
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE scaffold = 1")
  expect_equal(RiboPilot:::review_set_state(con, u[1, ], 1), "unlock first")
  expect_equal(re_unit(con, 1)$check_switch, 0L)
  expect_error(RiboPilot:::review_set_state(con, u, 2), "1 (Ready to run) or 0 (On hold)", fixed = TRUE)
})

test_that("each edit is one atomic write that also joins an open transaction", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  re_seed(con)
  orig <- re_originals(con)
  one <- data.frame(ID = "S1", path = 1L, scaffold = 1L)
  edits <- list(
    function() RiboPilot:::review_trim(con, "S1", 1L, 1L, 21, 80),
    function() RiboPilot:::review_set_bounds(con, "S1", 1L, 1L, "18S", 5, 40),
    function() RiboPilot:::review_set_partial(con, "S1", 1L, 1L, "18S", "start", 1),
    function() RiboPilot:::review_set_state(con, one, 1),
    function() RiboPilot:::review_set_state(con, one, 0))
  DBI::dbExecute(con, "CREATE TRIGGER re_boom BEFORE UPDATE ON review_units
    BEGIN SELECT RAISE(ABORT, 're_boom'); END")
  for (e in edits) {
    expect_error(e(), "re_boom")
    expect_equal(re_rows(con, 1), 0L)
  }
  DBI::dbExecute(con, "DROP TRIGGER re_boom")
  expect_false(RSQLite::sqliteIsTransacting(con))

  expect_error(DBI::dbWithTransaction(con, {
    for (e in edits) e()
    stop("outer failed")
  }), "outer failed")
  expect_equal(re_rows(con, 1), 0L)
  DBI::dbWithTransaction(con, {
    for (e in edits) e()
    RiboPilot:::review_reset_session(con, "S1", 1L, 1L, data.frame(element = "18S", start = 11L, stop = 40L,
      partial_start_user = NA, partial_stop_user = NA))
  })
  expect_equal(re_unit(con, 1)[, c("length", "check_switch")], data.frame(length = 60L, check_switch = 1L))
  expect_equal(re_feats(con, 1)$start[1:2], c(11L, 21L))
  expect_identical(re_originals(con), orig)
})
