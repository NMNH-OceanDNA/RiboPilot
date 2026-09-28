# Batch Trim and Restore from the Assemble and Annotation Review toolbars.

FIVE <- list(c("18S", 51, 200), c("ITS1", 201, 260), c("5.8S", 261, 320),
  c("ITS2", 321, 400), c("28S", 401, 550))

# mk_app_project()'s four units, 600 bp each: p1.s1 and p2.s2 have all five
# features, p1.s2 lacks ITS2 and 28S, p2.s1 is annotation-locked.
mk_batch_project <- function(env = parent.frame()) {
  dir <- mk_app_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  for (sql in c("ALTER TABLE assemblies ADD COLUMN edited INT",
    "ALTER TABLE annotations ADD COLUMN length INT", "ALTER TABLE annotations ADD COLUMN edited INT",
    "ALTER TABLE annotations ADD COLUMN strand TEXT", "ALTER TABLE coverage ADD COLUMN clips TEXT")) {
    DBI::dbExecute(con, sql)
  }
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 600", params = list(strrep("ACGT", 150)))
  for (u in list(c(1, 1), c(1, 2), c(2, 1), c(2, 2))) {
    feats <- if (all(u == c(1, 2))) FIVE[1:3] else FIVE
    for (f in feats) DBI::dbExecute(con,
      "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand) VALUES ('S1', ?, ?, ?, ?, ?, '+')",
      params = list(u[1], u[2], f[1], as.integer(f[2]), as.integer(f[3])))
  }
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1 WHERE path = 2 AND scaffold = 1")
  dir
}

units4 <- data.frame(ID = "S1", path = c(1L, 1L, 2L, 2L), scaffold = c(1L, 2L, 1L, 2L))
len <- function(con) DBI::dbGetQuery(con, "SELECT length FROM assemblies ORDER BY path, scaffold")$length

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
    c("not trimmed", "not trimmed", "annotation locked", "not trimmed"))
})

test_that("batch trim cuts to gene ends or a flank, sets no partial ends, and restores", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ok <- function(p) p[is.na(p$reason), , drop = FALSE]

  err <- RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(err, c(NA_character_, NA_character_))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assembly_backup")$n, 2L)
  expect_false(any(DBI::dbGetQuery(con,
    "SELECT COALESCE(partial_start, 0) + COALESCE(partial_stop, 0) n FROM annotations")$n > 0))

  # trimmed again with a 20 bp flank: the unit already sits at its gene ends,
  # so the flank is clipped to the scaffold and nothing changes
  expect_equal(RiboPilot:::batch_trim_plan(con, units4, flank = 20)$reason[1], "nothing to trim")

  p <- RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)
  expect_equal(p$reason, c(NA, "not trimmed", "annotation locked", NA))
  RiboPilot:::batch_trim_apply(con, ok(p), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))

  # flank mode, then gene ends: restore still returns the original
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, flank = 20)))
  expect_equal(len(con), c(540L, 600L, 600L, 540L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4)))
  expect_equal(len(con), c(500L, 600L, 600L, 500L))
  RiboPilot:::batch_trim_apply(con, ok(RiboPilot:::batch_trim_plan(con, units4, restore = TRUE)), restore = TRUE)
  expect_equal(len(con), rep(600L, 4))
  expect_equal(DBI::dbGetQuery(con,
    "SELECT start, stop FROM annotations WHERE path = 1 AND scaffold = 1 AND element = '18S'"),
    data.frame(start = 51L, stop = 200L))
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

test_that("the Trim and Restore buttons on both tables act on eligible units only", {
  skip_if_not_installed("withr")
  dir <- mk_batch_project()
  withr::local_dir(dir)
  toasts <- character(0)
  local_mocked_bindings(mp_toast = function(text, ...) toasts <<- c(toasts, text))
  lens <- function() RiboPilot:::app_query("SELECT length FROM assemblies ORDER BY path, scaffold")$length

  shiny::testServer(RiboPilot:::app_server, {
    # Assemble: one ticked sample covers all four of its scaffolds
    session$setInputs(mode = "Assemble")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1L, asmb_trim_btn = 1)
    html <- output$batch_trim_summary$html
    expect_match(html, "2 units are trimmed", fixed = TRUE)
    expect_match(html, "S1 p1 s2: missing ITS2, 28S", fixed = TRUE)
    expect_match(html, "S1 p2 s1: annotation locked", fixed = TRUE)
    expect_match(html, "max-height:200px", fixed = TRUE)
    session$setInputs(batch_trim_mode = "flank", batch_trim_flank = 20, batch_trim_go = 1)
    expect_equal(lens(), c(540L, 600L, 600L, 540L))
    expect_equal(utils::tail(toasts, 1), "2 units trimmed, 2 skipped.")

    session$setInputs(asmb_restore_btn = 1)
    expect_match(output$batch_trim_summary$html, "S1 p1 s2: not trimmed", fixed = TRUE)
    session$setInputs(batch_trim_go = 2)
    expect_equal(lens(), rep(600L, 4))

    # Annotation Review: only the ticked unit
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    row <- which(d$path == 2 & d$scaffold == 2)
    session$setInputs(tbl__reactable__selected = row, annot_trim_btn = 1)
    session$setInputs(batch_trim_mode = "annotated", batch_trim_go = 3)
    expect_equal(lens(), c(600L, 600L, 600L, 500L))
    expect_equal(utils::tail(toasts, 1), "1 unit trimmed.")
    session$setInputs(annot_restore_btn = 1, batch_trim_go = 4)
    expect_equal(lens(), rep(600L, 4))
    expect_equal(utils::tail(toasts, 1), "1 unit restored.")
  })
})
