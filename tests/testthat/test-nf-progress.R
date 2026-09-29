test_that("a redrawn process row overwrites rather than repeats", {
  s <- RiboPilot:::nf_progress_empty()
  s <- RiboPilot:::nf_progress_parse(s, c("N E X T F L O W  ~  version 25.10.6",
    "executor >  local (1)", "[a1/b2c3] preprocess (S1) | 0 of 1"))
  s <- RiboPilot:::nf_progress_parse(s, "[a1/b2c3] preprocess (S1) | 1 of 1")
  txt <- RiboPilot:::nf_progress_text(s)
  expect_equal(length(s$proc), 1L)
  expect_equal(lengths(regmatches(txt, gregexpr("preprocess", txt, fixed = TRUE))), 1L)
  expect_false(grepl("0 of 1", txt, fixed = TRUE))
  expect_true(grepl("1 of 1", txt, fixed = TRUE))
  expect_true(grepl("executor >  local (1)", txt, fixed = TRUE))
})

test_that("the newest executor line wins", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("executor >  local (1)", "executor >  local (4)"))
  expect_equal(s$executor, "executor >  local (4)")
})

test_that("banner lines are header and post-process lines are footer", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("N E X T F L O W", "[a1/b2c3] count_clade (Metazoa) | 1 of 1", "Completed at: now"))
  expect_equal(s$header, "N E X T F L O W")
  expect_equal(s$footer, "Completed at: now")
  expect_equal(names(s$proc), "count_clade")
})

test_that("sections render in header, executor, process, footer order", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(), c(
    "N E X T F L O W  ~  version 25.10.6",
    "Launching `main.nf`",
    "executor >  local (2)",
    "[a1/b2c3] preprocess (S1) | 1 of 1",
    "[d4/e5f6] build_ref_db (actinopteri) | 1 of 1",
    "Completed at: now",
    "Duration: 1m"
  ))
  expect_equal(RiboPilot:::nf_progress_text(s), paste(
    "N E X T F L O W  ~  version 25.10.6",
    "Launching `main.nf`",
    "executor >  local (2)",
    "[a1/b2c3] preprocess (S1) | 1 of 1",
    "[d4/e5f6] build_ref_db (actinopteri) | 1 of 1",
    "Completed at: now",
    "Duration: 1m",
    sep = "\n"))
})

test_that("ANSI escapes are stripped but the line survives", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    "\033[32m[a1/b2c3] build_ref_db (actinopteri) | 1 of 1\033[0m")
  txt <- RiboPilot:::nf_progress_text(s)
  expect_false(grepl("\033", txt, fixed = TRUE))
  expect_equal(names(s$proc), "build_ref_db")
  expect_equal(txt, "[a1/b2c3] build_ref_db (actinopteri) | 1 of 1")
})

test_that("a cached process row is recognised and keyed like any other", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    "[a1/b2c3] build_ref_db (actinopteri) | 0 of 1")
  s <- RiboPilot:::nf_progress_parse(s,
    "[-        ] build_ref_db (actinopteri) | 1 of 1, cached: 1")
  expect_equal(names(s$proc), "build_ref_db")
  expect_equal(RiboPilot:::nf_progress_text(s),
    "[-        ] build_ref_db (actinopteri) | 1 of 1, cached: 1")
  expect_equal(s$footer, character(0))
})

test_that("blank lines are skipped rather than banked as header or footer", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("", "   ", "N E X T F L O W", "",
      "[a1/b2c3] build_ref_db (x) | 1 of 1", "", "  ", "Completed at: now", ""))
  expect_equal(s$header, "N E X T F L O W")
  expect_equal(s$footer, "Completed at: now")
  expect_equal(RiboPilot:::nf_progress_text(s), paste(
    "N E X T F L O W",
    "[a1/b2c3] build_ref_db (x) | 1 of 1",
    "Completed at: now",
    sep = "\n"))
})

test_that("an empty state renders as an empty string", {
  expect_equal(RiboPilot:::nf_progress_text(RiboPilot:::nf_progress_empty()), "")
})

test_that("the footer keeps only the most recent lines", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("[a1/b2c3] build_ref_db (x) | 1 of 1", paste0("note ", 1:50)))
  txt <- RiboPilot:::nf_progress_text(s, footer_n = 40)
  expect_equal(length(s$footer), 50L)
  expect_false(grepl("note 10\n", txt, fixed = TRUE))
  expect_true(grepl("note 11\n", txt, fixed = TRUE))
  expect_true(grepl("note 50", txt, fixed = TRUE))
})

test_that("the footer stops growing once it is past the cap", {
  n <- RiboPilot:::NF_FOOTER_MAX
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("[a1/b2c3] build_ref_db (x) | 1 of 1", paste0("note ", seq_len(n + 500))))
  expect_equal(length(s$footer), n)
  expect_equal(utils::tail(s$footer, 1), paste0("note ", n + 500))
  # a long run appends line by line, not in one block
  s2 <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    "[a1/b2c3] build_ref_db (x) | 1 of 1")
  for (i in seq_len(n + 10)) s2 <- RiboPilot:::nf_progress_parse(s2, paste0("warn ", i))
  expect_equal(length(s2$footer), n)
  expect_true(grepl(paste0("warn ", n + 10), RiboPilot:::nf_progress_text(s2), fixed = TRUE))
})

test_that("the hidden-process line is kept once and cleared by a later redraw", {
  more <- "Plus 1 more processes waiting for tasks..."
  row <- "[a1/b2c3] preprocess (S1) | 0 of 1"
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    rep(c("executor >  local (1)", row, more), 3))
  txt <- RiboPilot:::nf_progress_text(s)
  expect_equal(lengths(regmatches(txt, gregexpr("Plus 1 more", txt, fixed = TRUE))), 1L)
  expect_length(s$footer, 0L)
  s <- RiboPilot:::nf_progress_parse(s, c("executor >  local (2)", row))
  expect_false(grepl("Plus 1 more", RiboPilot:::nf_progress_text(s), fixed = TRUE))
})
