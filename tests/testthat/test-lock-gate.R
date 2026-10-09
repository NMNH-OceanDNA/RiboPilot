# Lock on Annotation Review waits for finished checks (S6.2); Exclude & Lock does not.

lock_gate_units <- function(con) {
  seq <- paste(rep("ACGT", 50), collapse = "")
  DBI::dbExecute(con, "DELETE FROM annotate")
  # scaffold: annotate_switch, n_features, itsx_warnings, problematic
  rows <- list(c(1, 2, 3, "", "no"), c(2, 1, 3, NA, "no"), c(3, 2, 3, NA, "no"), c(4, 2, 3, "", "no"),
    c(5, 2, 3, "", "no"), c(6, 2, 3, "", "no"), c(7, 2, 3, "", "no"), c(8, 2, 3, "", "no"),
    c(9, 1, 3, NA, "yes"), c(10, 2, NA, NA, "no"))
  for (r in rows) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
      VALUES ('S1', 1, ?, 'linear', 200, ?, 0)", params = list(as.integer(r[1]), seq))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features,
      itsx_warnings, problematic) VALUES ('S1', 1, ?, ?, 0, ?, ?, ?)",
      params = list(as.integer(r[1]), as.integer(r[2]), as.integer(r[3]), r[4], r[5]))
  }
  for (s in 4:8) RiboPilot:::review_copy(con, "S1", 1L, s)
  for (x in list(c(4, 1), c(5, 4), c(6, 3), c(7, 0), c(8, 2))) {
    DBI::dbExecute(con, "UPDATE review_units SET check_switch = ? WHERE ID = 'S1' AND path = 1 AND scaffold = ?",
      params = list(as.integer(x[2]), as.integer(x[1])))
  }
}

test_that("review_lock_gate names the next step for each state and passes done and excluded contigs", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  lock_gate_units(con)
  u <- data.frame(ID = "S1", path = 1L, scaffold = c(8L, 1:7, 9:10))
  expect_identical(RiboPilot:::review_lock_gate(con, u), c(
    NA,                                                  # copy, done
    NA,                                                  # no copy, done
    "Press Update on the Assemble step first.",          # no copy, queued for the main pipeline
    "Set State to Ready to run, then press Update.",     # no copy, ITSx never checked
    "Press Update first.",                               # copy, not run
    "Wait for the running Update to finish.",            # copy, running
    "Set State to Ready to run, then press Update.",     # copy, failed
    "Set State to Ready to run, then press Update.",     # copy, on hold
    NA,                                                  # excluded, though queued
    "Press Update on the Assemble step first."))         # no copy, never annotated
  expect_identical(RiboPilot:::review_lock_gate(con, u[0, ]), character(0))
})

test_that("Lock on Annotation Review refuses contigs whose checks are not done", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  said <- character(0)
  local_mocked_bindings(mp_alert = function(title, text = NULL, ...) said <<- c(said, paste(title, text)))
  lk <- function() RiboPilot:::app_query(
    "SELECT COALESCE(annotate_lock, 0) AS v FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")$v
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- function(p, s) which(d$path == p & d$scaffold == s)[1]

    # every contig is waiting for the main pipeline: toolbar, row icon, and details window refuse
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(1, 2)), annot_lock_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_match(said[length(said)], "Checks not finished", fixed = TRUE)
    expect_match(said[length(said)], "S1 p1.s1, S1 p1.s2: Press Update on the Assemble step first.", fixed = TRUE)
    session$setInputs(lock_row = as.character(row(1, 1)))
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    n <- length(said)
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(lock_and_close = 1)
    session$flushReact()
    session$setInputs(modal_lock = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_equal(length(said), n + 2L)

    # one done contig with one not done: nothing is locked
    set("UPDATE annotate SET annotate_switch = 2, itsx_warnings = '' WHERE path = 1 AND scaffold = 1")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(tbl__reactable__selected = c(row(1, 1), row(1, 2)), annot_lock_btn = 2)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    expect_match(said[length(said)], "S1 p1.s2: Press Update on the Assemble step first.", fixed = TRUE)
    expect_false(grepl("p1.s1", said[length(said)], fixed = TRUE))

    # the done contig alone locks, from the toolbar and from the details window
    session$setInputs(tbl__reactable__selected = row(1, 1), annot_lock_btn = 3)
    session$flushReact()
    expect_equal(lk(), c(1L, 0L, 0L, 0L))
    set("UPDATE annotate SET annotate_lock = 0")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(modal_lock = 2)
    session$flushReact()
    expect_equal(lk(), c(1L, 0L, 0L, 0L))

    # Exclude & Lock is exempt, and so is a plain Lock on an excluded contig
    set("UPDATE annotate SET annotate_lock = 0")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(tbl__reactable__selected = row(1, 2), problematic_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 1L, 0L, 0L))
    set("UPDATE annotate SET problematic = 'yes' WHERE path = 2 AND scaffold = 1")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(lock_row = as.character(row(2, 1)))
    session$flushReact()
    expect_equal(lk(), c(0L, 1L, 1L, 0L))
  })
})

test_that("Lock on Annotation Review asks about flags found on the review copy", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  fixture_mark_done(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 2, itsx_warnings = '', curate_flags = 'contaminant'
    WHERE ID = 'S1' AND path = 1 AND scaffold = 1")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  local_mocked_bindings(mp_alert = function(...) NULL)
  lk <- function() RiboPilot:::app_query(
    "SELECT COALESCE(annotate_lock, 0) AS v FROM annotate WHERE ID = 'S1' ORDER BY path, scaffold")$v
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- function(p, s) which(d$path == p & d$scaffold == s)[1]
    session$setInputs(tbl__reactable__selected = row(1, 1), annot_lock_btn = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
    session$setInputs(lock_confirm = TRUE)
    expect_equal(lk(), c(1L, 0L, 0L, 0L))

    # the details window asks too
    set("UPDATE annotate SET annotate_lock = 0")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(modal_lock = 1)
    session$flushReact()
    expect_equal(lk(), c(0L, 0L, 0L, 0L))
  })
})
