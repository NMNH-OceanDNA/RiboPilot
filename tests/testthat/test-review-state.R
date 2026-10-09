# Check state on Annotation Review (SPEC S6.1, S6.3, S8).

# Contigs 1-10 of S1 (assembly locked) cover every state, contig 11 is dropped,
# and S2's assembly is not locked. Copies hold only the columns read here.
rst_units <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'S1'")
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon) VALUES ('S2', 'Danio rerio')")
  DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_lock) VALUES ('S2', 0)")
  orig <- data.frame(
    ID = c(rep("S1", 11), "S2"), scaffold = c(1:11, 1L),
    annotate_switch = c(1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 4L, 2L, 1L),
    annotate_lock = c(rep(0L, 9), 1L, 0L, 0L),
    n_features = c(2L, NA, rep(2L, 10)),
    itsx_warnings = c(NA, NA, NA, "", NA, "", "", "", "", NA, "", NA),
    ignore = c(rep(0L, 10), 1L, 0L))
  for (i in seq_len(nrow(orig))) {
    o <- orig[i, ]
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES (?, 1, ?, 4, 'ACGT', ?)", params = list(o$ID, o$scaffold, o$ignore))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
      n_features, itsx_warnings) VALUES (?, 1, ?, ?, ?, ?, ?)",
      params = list(o$ID, o$scaffold, o$annotate_switch, o$annotate_lock, o$n_features, o$itsx_warnings))
  }
  # copies: 5 failed (its original is queued again; the copy wins), 6 running,
  # 7 on hold, 8 with no state stored, 9 done
  cp <- data.frame(scaffold = 5:9, check_switch = c(3L, 4L, 0L, NA, 2L),
    check_note = c("boom", NA, NA, NA, NA))
  for (i in seq_len(nrow(cp))) {
    DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, trim_start,
      trim_end, edited, check_switch, check_note, rev, n_features, itsx_warnings)
      VALUES ('S1', 1, ?, 'ACGT', 4, 1, 4, 1, ?, ?, CAST(? AS INTEGER), 2, '')",
      params = list(cp$scaffold[i], cp$check_switch[i], cp$check_note[i], "1760000000000"))
  }
}

test_that("review_unit_state names every state, and waiting matches the recheck's waiting count", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rst_units(con)
  st <- RiboPilot:::review_unit_state(con)
  expect_equal(names(st), c("ID", "path", "scaffold", "has_copy", "check_switch", "annotate_switch",
    "annotate_lock", "n_features", "check_note", "state"))
  expect_equal(st$ID, rep("S1", 10))
  expect_equal(st$scaffold, 1:10)
  expect_equal(st$has_copy, c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 0L))
  expect_equal(st$state, c("waiting", "waiting", "not run", "done", "failed", "running", "on hold",
    "not run", "done", "waiting"))
  expect_equal(st$check_note, c(NA, NA, NA, NA, "boom", NA, NA, NA, NA, NA))
  expect_equal(RiboPilot:::review_state_text(st$state), c("waiting for Assemble Update",
    "waiting for Assemble Update", "not run", "done", "failed", "running", "on hold", "not run",
    "done", "waiting for Assemble Update"))
  # the Update alert counts the same contigs, less the locked one
  expect_equal(RiboPilot:::nf_recheck_waiting(dirname(con@dbname)), 2L)

  # given units: their order, any lock or drop state, and a unit that does not exist
  u <- data.frame(ID = c("S1", "S2", "S1", "SX"), path = 1L, scaffold = c(9L, 1L, 11L, 1L))
  st <- RiboPilot:::review_unit_state(con, u)
  expect_equal(st[c("ID", "path", "scaffold")], u)
  expect_equal(st$state, c("done", "waiting", "done", "not run"))
  expect_equal(nrow(RiboPilot:::review_unit_state(con, u[0, ])), 0L)
})

test_that("unit_checks on the review layer reads the copies and shows the state of unfinished checks", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotate")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence)
      VALUES ('S1', 1, ?, 4, 'ACGT')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, n_features, itsx_warnings)
      VALUES ('S1', 1, ?, ?, 2, '')", params = list(s, if (s == 4) 1L else 2L))
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, ?, 'on_target')",
      params = list(s))
  }
  copy <- function(s, sw, flags, itsx) DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold,
      sequence, length, trim_start, trim_end, edited, check_switch, rev, n_features, curate_flags,
      itsx_warnings) VALUES ('S1', 1, ?, 'ACG', 3, 1, 3, 1, ?, CAST(? AS INTEGER), 2, ?, ?)",
    params = list(s, sw, "1760000000000", flags, itsx))
  # 2: a done copy the screen calls off target, with a contaminant flag, a
  # flagged feature, and an ITSx check that could not decide
  copy(2L, 2L, "contaminant", "not_run(ITSx did not finish)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 2, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, status)
    VALUES ('S1', 1, 2, '28S', 'flag')")
  # 3: an edited copy not rechecked yet, under the severe chimera a discarded recheck left
  copy(3L, 1L, "", NA)
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label) VALUES ('S1', 1, 3, 'severe')")

  u <- RiboPilot:::unit_checks(con, "unit", review = TRUE)
  expect_equal(u$scaffold, 1:4)
  expect_equal(u$checks, c("ok", "off target, contaminant, ITSx not checked, feature flags", "not run",
    "waiting for Assemble Update"))
  expect_equal(u$checks_tone, c("success", "danger", "info", "info"))

  # the Assemble step reads the originals, which no copy changed
  o <- RiboPilot:::unit_checks(con, "unit")
  expect_equal(o$checks, rep("ok", 4))
  expect_equal(o$checks_tone, rep("success", 4))
})
