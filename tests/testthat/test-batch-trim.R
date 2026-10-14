# Batch Trim and Restore from the Assemble and Annotation Review toolbars.

FIVE <- list(c("18S", 51, 200), c("ITS1", 201, 260), c("5.8S", 261, 320),
  c("ITS2", 321, 400), c("28S", 401, 550))

# mk_app_project()'s four units, 600 bp each: p1.s1 and p2.s2 have all five
# features, p1.s2 lacks ITS2 and 28S, p2.s1 is annotation-locked.
mk_batch_project <- function(env = parent.frame()) {
  dir <- mk_app_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 600", params = list(withr::with_seed(1, paste(sample(c("A", "C", "G", "T"), 600, TRUE), collapse = ""))))
  for (u in list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))) {
    feats <- if (all(u == c(1, 2))) FIVE[1:3] else FIVE
    for (f in feats) DBI::dbExecute(con,
      "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand) VALUES ('S1', ?, ?, ?, ?, ?, '+')",
      params = list(u[1], u[2], f[1], as.integer(f[2]), as.integer(f[3])))
  }
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE path = 2 AND scaffold = 1")
  fixture_mark_done(con)
  dir
}

units4 <- data.frame(ID = "S1", path = c(1L, 1L, 2L, 2L), scaffold = c(1L, 2L, 1L, 2L))
len <- function(con) DBI::dbGetQuery(con, "SELECT length FROM review_units_v ORDER BY path, scaffold")$length

test_that("only unlocked units with all five features are eligible, with reasons", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  p <- RiboPilot:::batch_trim_plan(con, units4)
  expect_equal(p$reason, c(NA, "missing ITS2, 28S", "annotation locked", NA))
  expect_equal(p$from, c(51L, NA, NA, 51L))
  expect_equal(p$to, c(550L, NA, NA, 550L))
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, flank = 20)$to, c(570L, NA, NA, 570L))
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)$reason,
    c("no edits to restore", "no edits to restore", "annotation locked", "no edits to restore"))
})

test_that("contigs waiting for the Assemble step's Update are refused; copies are trimmed", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 4")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 4")
  w <- "Press Update on the Assemble step first."
  expect_equal(RiboPilot:::batch_trim_plan(con, units4)$reason, c(w, w, "annotation locked", w))
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  expect_equal(RiboPilot:::batch_trim_plan(con, units4)$reason[1], NA_character_)
  expect_false("running" %in% names(formals(RiboPilot:::batch_trim_plan)))
})

test_that("batch trim cuts review copies to gene ends or a flank, and restore deletes them", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ok <- function(p) p[is.na(p$reason), , drop = FALSE]
  orig <- function() list(
    DBI::dbGetQuery(con, "SELECT path, scaffold, length, sequence FROM assemblies ORDER BY path, scaffold"),
    DBI::dbGetQuery(con, "SELECT path, scaffold, element, start, stop, partial_start, partial_stop
      FROM annotations ORDER BY path, scaffold, start"))
  before <- orig()

  err <- RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(err, c(NA_character_, NA_character_))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT path, scaffold, trim_start, trim_end, edited, check_switch
    FROM review_units ORDER BY path, scaffold"),
    data.frame(path = c(1L, 2L), scaffold = c(1L, 2L), trim_start = 51L, trim_end = 550L,
      edited = 1L, check_switch = 1L))
  # a cut exactly at the gene ends leaves no partial end
  expect_false(any(DBI::dbGetQuery(con,
    "SELECT COALESCE(partial_start_user, 0) + COALESCE(partial_stop_user, 0) n FROM review_annotations")$n > 0))
  expect_equal(orig(), before)

  # trimmed again with a 20 bp flank: the copy already sits at its gene ends,
  # so the flank is clipped to the contig and nothing changes
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, flank = 20)$reason[1], "nothing to trim")

  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c(NA, "no edits to restore", "annotation locked", NA))
  expect_equal(RiboPilot:::batch_trim_apply(con, ok(p), restore = TRUE), c(NA_character_, NA_character_))
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM review_units")$n, 0L)

  # flank mode, then gene ends: the second trim plans from the copy, in copy coordinates
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, flank = 20)))
  expect_equal(len(con), c(540L, 600L, 600L, 540L))
  p <- RiboPilot:::batch_trim_plan(con, units4)
  expect_equal(c(p$from[1], p$to[1]), c(21L, 520L))
  RiboPilot:::batch_trim_apply(con, ok(p))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT trim_start, trim_end FROM review_units WHERE path = 1 AND scaffold = 1"),
    data.frame(trim_start = 51L, trim_end = 550L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT start, stop FROM review_annotations_v WHERE path = 1 AND scaffold = 1 AND element = '18S'"),
    data.frame(start = 51L, stop = 200L))
  expect_equal(orig(), before)
})

test_that("one failing unit does not stop the rest", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  p <- RiboPilot:::batch_trim_plan(con, units4)
  p <- p[is.na(p$reason), , drop = FALSE]
  p$ID[1] <- "nope"
  err <- RiboPilot:::batch_trim_apply(con, p)
  expect_match(err[1], "No sequence")
  expect_true(is.na(err[2]))
  expect_equal(len(con)[4], 500L)
})

test_that("the Trim and Restore buttons act on eligible ticked units only", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  withr::local_dir(dir)
  toasts <- character(0)
  local_mocked_bindings(mp_toast = function(text, ...) toasts <<- c(toasts, text))
  lens <- function() RiboPilot:::app_query("SELECT length FROM review_units_v ORDER BY path, scaffold")$length
  tip <- RiboPilot:::MP_TRIM_TIP()
  expect_false(grepl("running run", tip[["trim"]], fixed = TRUE))
  expect_match(tip[["restore"]], "boundary edits", fixed = TRUE)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = seq_len(nrow(samples_data())), annot_trim_btn = 1)
    html <- output$batch_trim_summary$html
    expect_match(html, "2 contigs are trimmed. Checks show \"not run\" until you run Update in the Annotation Review step.", fixed = TRUE)
    expect_false(grepl("re-curated", html, fixed = TRUE))
    expect_match(html, "S1 p1 s2: missing ITS2, 28S", fixed = TRUE)
    expect_match(html, "S1 p2 s1: annotation locked", fixed = TRUE)
    expect_match(html, "max-height:200px", fixed = TRUE)
    session$setInputs(batch_trim_mode = "flank", batch_trim_flank = 20, batch_trim_go = 1)
    expect_equal(lens(), c(540L, 600L, 600L, 540L))
    expect_equal(utils::tail(toasts, 1), "2 contigs trimmed, 2 skipped.")

    session$setInputs(annot_restore_btn = 1)
    expect_match(output$batch_trim_summary$html, "S1 p1 s2: no edits to restore", fixed = TRUE)
    expect_match(output$batch_trim_summary$html, "boundary edits, and 5'/3' choices are deleted", fixed = TRUE)
    session$setInputs(batch_trim_go = 2)
    expect_equal(lens(), rep(600L, 4))

    # only the ticked unit
    d <- samples_data()
    row <- which(d$path == 2 & d$scaffold == 2)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1)
    session$setInputs(batch_trim_mode = "annotated", batch_trim_go = 3)
    expect_equal(lens(), c(600L, 600L, 600L, 500L))
    expect_equal(utils::tail(toasts, 1), "1 contig trimmed.")
    session$setInputs(annot_restore_btn = 1, batch_trim_go = 4)
    expect_equal(lens(), rep(600L, 4))
    expect_equal(utils::tail(toasts, 1), "1 contig restored.")
  })
})

test_that("only the trimmed end of a copy stops being an assembly end", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE annotations SET start = 1 WHERE path = 1 AND scaffold = 1 AND element = '18S'")
  expect_equal(RiboPilot:::unit_trimmed(con, "S1", 1L, 1L), c(FALSE, FALSE))
  RiboPilot:::review_trim(con, "S1", 1L, 1L, 1, 550)
  tr <- RiboPilot:::unit_trimmed(con, "S1", 1L, 1L)
  expect_equal(tr, c(FALSE, TRUE))
  d <- DBI::dbGetQuery(con, paste("SELECT element, start, stop, strand,", RiboPilot:::REVIEW_EFF_PARTIALS,
    "FROM review_annotations_v WHERE path = 1 AND scaffold = 1 ORDER BY start"))
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p5[1])
  expect_false(e$p3[5])
  # on the minus strand the untouched left end is the 3' end
  d$strand <- "-"
  e <- RiboPilot:::partial_ends(d, 550, tr)
  expect_true(e$p3[1])
  expect_false(e$p5[5])
  RiboPilot:::batch_trim_apply(con, units4[1, ], restore = TRUE)
  expect_equal(RiboPilot:::unit_trimmed(con, "S1", 1L, 1L), c(FALSE, FALSE))
})

test_that("batch Restore offers every contig with a copy, not only trimmed ones", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  # an untrimmed copy, as State makes one
  RiboPilot:::review_copy(con, "S1", 2L, 2L)
  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c("no edits to restore", "no edits to restore", "annotation locked", NA))
  expect_equal(RiboPilot:::batch_trim_apply(con, p[4, ], restore = TRUE), NA_character_)
  expect_false(RiboPilot:::review_has_copy(con, "S1", 2L, 2L))
})

test_that("ETS and IGS do not widen the gene-ends window", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1', 1, 1, 'ETS', 1, 50, '+'), ('S1', 1, 1, 'IGS', 551, 600, '+')")
  p <- RiboPilot:::batch_trim_plan(con, units4[1, ])
  expect_equal(c(p$from, p$to), c(51L, 550L))
})

test_that("a batch trim refreshes an open unit window, and flank mode needs a flank", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  withr::local_dir(dir)
  local_mocked_bindings(mp_toast = function(text, ...) NULL)
  db <- function() RiboPilot:::app_query(
    "SELECT element, start, stop FROM review_annotations_v WHERE path = 2 AND scaffold = 2 ORDER BY start")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- which(d$path == 2 & d$scaffold == 2)
    session$setInputs(ann_details = row)
    expect_equal(feature_rows()$start[1], 51L)
    expect_equal(unname(trim_span()), c(51, 550, 600))
    # an Edit session left open when the window closes
    session$setInputs(feat_edit = TRUE)
    expect_false(is.null(edit_backup()))
    session$setInputs(annotate_close = 1)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1)

    session$setInputs(batch_trim_mode = "flank", batch_trim_flank = NA)
    expect_match(output$batch_trim_summary$html, "Enter a flank of at least 1 bp.", fixed = TRUE)
    session$setInputs(batch_trim_go = 1)
    expect_equal(db()$start[1], 51L)
    expect_false(is.null(edit_backup()))

    # the trim ends that session, so a later Reset cannot write pre-trim coordinates
    session$setInputs(batch_trim_mode = "annotated", batch_trim_go = 2)
    expect_equal(db()$start[1], 1L)
    expect_null(edit_backup())
    expect_false(edit_active())
    session$setInputs(ann_details = row)
    expect_equal(unname(trim_span()), c(1, 500, 500))
  })
})

test_that("the viewer is re-sent after a trim, not only when the unit changes", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  withr::local_dir(dir)
  local_mocked_bindings(mp_toast = function(text, ...) NULL)
  drawn <- 0L
  local_mocked_bindings(seqview_send = function(...) { drawn <<- drawn + 1L; invisible() })

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- which(d$path == 1 & d$scaffold == 1)
    session$setInputs(ann_details = row)
    session$flushReact()
    n <- drawn
    session$setInputs(annotate_close = 1)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1,
      batch_trim_mode = "annotated", batch_trim_go = 1)
    session$flushReact()
    expect_gt(drawn, n)
  })
})
