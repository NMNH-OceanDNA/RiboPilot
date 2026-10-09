# Edit mode and State on Annotation Review write the contig's review copy
# (SPEC S5.1, S5.2, D6, D7); the original tables never change.

# mk_app_project() with every contig 600 bp; p1.s1 holds 18S, ITS1, and 5.8S on
# the + strand. Every contig is done (annotated, ITSx stored) and unlocked.
mk_handler_project <- function(env = parent.frame()) {
  dir <- mk_app_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 600", params = list(strrep("ACGT", 150)))
  for (f in list(c("18S", 10, 200), c("ITS1", 201, 260), c("5.8S", 261, 420))) {
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, strand, status, flags)
      VALUES ('S1', 1, 1, ?, ?, ?, ?, '+', 'ok', '')",
      params = list(f[1], RiboPilot:::rrna_element_type(f[1]), as.integer(f[2]), as.integer(f[3])))
  }
  fixture_mark_done(con)
  dir
}

hq <- function(sql, params = NULL) RiboPilot:::app_query(sql, params)
hx <- function(sql) {
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, sql)
}

test_that("a nudge edits the review copy and queues its checks; marks and Save make no copy and run no curation", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  curated <- 0L
  local_mocked_bindings(curate_one = function(...) curated <<- curated + 1L)
  copies <- function() hq("SELECT COUNT(*) AS n FROM review_units")$n

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(mark_reviewed = 1, mark_partial = 1, annotate_notes = "check the 18S start")
    session$setInputs(tbl__reactable__selected = 1:4, reviewed_btn = 1, partial_btn = 1)
    session$setInputs(modal_hidden = 1)
    expect_equal(copies(), 0)

    expect_equal(active_feature(), "18S")
    session$setInputs(feat_edit = TRUE)
    r0 <- refresh()
    session$setInputs(nt_5_out = 1)
    expect_equal(hq("SELECT start FROM annotations WHERE element = '18S'")$start, 10L)
    expect_equal(hq("SELECT start, stop, length, edited, status FROM review_annotations WHERE element = '18S'"),
      data.frame(start = 9L, stop = 200L, length = 192L, edited = 1L, status = NA_character_))
    expect_equal(hq("SELECT path, scaffold, edited, check_switch FROM review_units"),
      data.frame(path = 1L, scaffold = 1L, edited = 1L, check_switch = 1L))
    expect_equal(feature_rows()$start, c(9L, 201L, 261L))
    expect_equal(nudge_ctx()$start, 9L)

    # the table reloads once the window closes, not on every nudge
    expect_equal(refresh(), r0)
    session$setInputs(modal_hidden = 2)
    expect_equal(refresh(), r0 + 1)

    session$setInputs(feat_save = 1)
    expect_equal(curated, 0L)
    expect_equal(hq("SELECT check_switch FROM review_units")$check_switch, 1L)
  })
})

test_that("Reset and Discard put back the session's bounds and 5'/3' choices", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  # curation's verdict, which Reset must not turn into a user choice
  hx("UPDATE annotations SET partial_stop = 1 WHERE element = '18S'")
  feat <- function() hq("SELECT start, stop, partial_start_user, partial_stop_user FROM review_annotations WHERE element = '18S'")
  unit <- function() hq("SELECT check_switch, CAST(rev AS TEXT) AS rev FROM review_units")
  raised <- function(rev) hq("SELECT rev - CAST(? AS INTEGER) AS d FROM review_units", list(rev))$d

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(feat_edit = TRUE)
    session$setInputs(toggle_partial_start = 1)
    expect_equal(feat(), data.frame(start = 10L, stop = 200L, partial_start_user = 1L, partial_stop_user = NA_integer_))
    u0 <- unit()
    # a 5'/3' choice leaves the checks as they were
    expect_equal(u0$check_switch, 2L)
    session$setInputs(feat_reset = 1)
    expect_equal(feat(), data.frame(start = 10L, stop = 200L, partial_start_user = NA_integer_,
      partial_stop_user = NA_integer_))
    expect_equal(unit(), u0)

    session$setInputs(nt_5_out = 1)
    expect_equal(feat()$start, 9L)
    expect_equal(raised(u0$rev), 1)
    hx("UPDATE review_units SET check_switch = 2")
    session$setInputs(feat_reset = 2)
    expect_equal(feat()$start, 10L)
    expect_equal(unit()$check_switch, 1L)
    expect_equal(raised(u0$rev), 2)

    # leaving Edit with unsaved nudges and choosing Discard does the same
    session$setInputs(feat_edit = FALSE)
    session$setInputs(feat_edit = TRUE)
    session$setInputs(nt_3_in = 1)
    expect_equal(feat()$stop, 199L)
    session$setInputs(feat_edit = FALSE)
    session$setInputs(unsaved_confirm = FALSE)
    expect_equal(feat()$stop, 200L)
    expect_equal(raised(u0$rev), 4)
  })
  expect_equal(hq("SELECT start, stop, partial_start, partial_stop FROM annotations WHERE element = '18S'"),
    data.frame(start = 10L, stop = 200L, partial_start = NA_integer_, partial_stop = 1L))
})

test_that("the details window reads the copy: features, effective 5'/3' ends, sequence, and length", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  # the copy keeps 181..600 of the original: 18S is gone and 5.8S runs to the copy's end
  DBI::dbExecute(con, "UPDATE review_units SET sequence = ?, length = 420, trim_start = 181, trim_end = 600",
    params = list(strrep("T", 420)))
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE element = '18S'")
  DBI::dbExecute(con, "UPDATE review_annotations SET start = start - 180, stop = stop - 180")
  DBI::dbExecute(con, "UPDATE review_annotations SET stop = 420 WHERE element = '5.8S'")
  # curation said 5' partial; the user said no
  DBI::dbExecute(con, "UPDATE review_annotations SET partial_start = 1, partial_start_user = 0 WHERE element = 'ITS1'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    expect_equal(active_feature(), "ITS1")
    f <- feature_rows()
    expect_equal(f$element, c("ITS1", "5.8S"))
    expect_equal(f$start, c(21L, 81L))
    expect_equal(f$partial_start, c(0L, 0L))
    expect_equal(f$partial_stop, c(0L, 1L))
    expect_equal(output$feat_edit_active, "Active feature: ITS1")
    active_feature("5.8S")
    session$flushReact()
    expect_equal(output$feat_edit_active, "Active feature: 5.8S (partial 3')")
    cx <- nudge_ctx()
    expect_equal(c(cx$start, cx$stop), c(81L, 420L))
    expect_equal(cx$seq, strrep("T", 420))
  })
})

test_that("State on Annotation Review queues or holds the recheck on copies and refuses locked and waiting contigs", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE path = 2 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE path = 2 AND scaffold = 2")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  said <- character(0)
  local_mocked_bindings(mp_alert = function(title, text = NULL, ...) said <<- c(said, paste(title, text)))
  copies <- function() hq("SELECT path, scaffold, check_switch FROM review_units ORDER BY path, scaffold")
  rev_of <- function(p, s) hq("SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE path = ? AND scaffold = ?",
    list(p, s))$rev

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    row <- function(p, s) { d <- samples_data(); which(d$path == p & d$scaffold == s) }
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(2, 1), row(2, 2)))
    session$setInputs(annot_state_btn = 1)
    session$setInputs(new_state = "1", update_state = 1)
    expect_equal(copies(), data.frame(path = 1L, scaffold = 1L, check_switch = 1L))
    expect_equal(said, paste("Left 2 contigs unchanged",
      "S1 p2.s1 (unlock first), S1 p2.s2 (press Update on the Assemble step first)."))

    # On hold keeps rev; Ready to run again raises it
    r1 <- rev_of(1, 1)
    session$setInputs(state_row = as.character(row(1, 1)))
    session$setInputs(new_state = "0", update_state = 2)
    expect_equal(copies()$check_switch, 0L)
    expect_equal(rev_of(1, 1), r1)
    session$setInputs(state_row = as.character(row(1, 1)))
    session$setInputs(new_state = "1", update_state = 3)
    expect_equal(copies()$check_switch, 1L)
    expect_equal(hq("SELECT rev - CAST(? AS INTEGER) AS d FROM review_units", list(r1))$d, 1)

    # On hold on a contig without a copy makes one, held
    session$setInputs(state_row = as.character(row(1, 2)))
    session$setInputs(new_state = "0", update_state = 4)
    expect_equal(copies(), data.frame(path = c(1L, 1L), scaffold = c(1L, 2L), check_switch = c(1L, 0L)))
  })
  # the main pipeline's queue is not touched
  expect_equal(hq("SELECT annotate_switch FROM annotate ORDER BY path, scaffold")$annotate_switch, c(2L, 2L, 2L, 1L))
})

test_that("a contig with no copy and no ITSx result shows On hold, not Ready to run", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  hx("UPDATE annotate SET itsx_warnings = NULL WHERE path = 1 AND scaffold = 2")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    sw <- function(p, s) d$annotate_switch[d$path == p & d$scaffold == s]
    expect_equal(sw(1, 2), 0L)
    expect_equal(sw(1, 1), 2L)
  })
})

test_that("a nudge that clamps to the same bounds makes no copy", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  hx("UPDATE annotations SET start = 1 WHERE element = '18S'")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(feat_edit = TRUE)
    session$setInputs(nt_5_out = 1)
    expect_equal(hq("SELECT COUNT(*) AS n FROM review_units")$n, 0)
    expect_false(isTRUE(edit_dirty()))
  })
})

test_that("Edit and Trim refuse a contig waiting for the Assemble step's Update", {
  skip_if_not_installed("withr")
  dir <- mk_handler_project()
  withr::local_dir(dir)
  hx("UPDATE annotate SET annotate_switch = 1 WHERE path = 1 AND scaffold = 1")
  said <- character(0)
  local_mocked_bindings(rp_toast = function(text, ...) said <<- c(said, text),
    rp_alert = function(text, ...) said <<- c(said, text))
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(feat_edit = TRUE)
    expect_false(isTRUE(edit_active()))
    session$setInputs(trim_btn = 1)
    expect_equal(said, rep("Press Update on the Assemble step first.", 2))
    expect_equal(hq("SELECT COUNT(*) AS n FROM review_units")$n, 0)
  })
})
