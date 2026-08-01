test_that("a redrawn process row overwrites rather than repeats", {
  s <- RiboPilot:::nf_progress_empty()
  s <- RiboPilot:::nf_progress_parse(s, c("N E X T F L O W  ~  version 25.10.6",
    "executor >  local (1)", "[a1/b2c3] preprocess (S1) | 0 of 1"))
  s <- RiboPilot:::nf_progress_parse(s, "[a1/b2c3] preprocess (S1) | 1 of 1")
  txt <- RiboPilot:::nf_progress_text(s)
  expect_equal(length(gregexpr("preprocess", txt)[[1]]), 1L)
  expect_true(grepl("1 of 1", txt, fixed = TRUE))
  expect_true(grepl("executor >  local (1)", txt, fixed = TRUE))
})

test_that("banner lines are header and post-process lines are footer", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("N E X T F L O W", "[a1/b2c3] count_clade (Metazoa) | 1 of 1", "Completed at: now"))
  expect_true(grepl("N E X T F L O W", RiboPilot:::nf_progress_text(s), fixed = TRUE))
  expect_true(grepl("Completed at: now", RiboPilot:::nf_progress_text(s), fixed = TRUE))
})

test_that("ANSI escapes are stripped", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    "\033[32m[a1/b2c3] build_ref_db (actinopteri) | 1 of 1\033[0m")
  expect_false(grepl("\033", RiboPilot:::nf_progress_text(s), fixed = TRUE))
})

test_that("an empty state renders as an empty string", {
  expect_equal(RiboPilot:::nf_progress_text(RiboPilot:::nf_progress_empty()), "")
})

test_that("the footer keeps only the most recent lines", {
  s <- RiboPilot:::nf_progress_parse(RiboPilot:::nf_progress_empty(),
    c("[a1/b2c3] build_ref_db (x) | 1 of 1", paste0("note ", 1:50)))
  txt <- RiboPilot:::nf_progress_text(s, footer_n = 40)
  expect_false(grepl("note 10\n", txt, fixed = TRUE))
  expect_true(grepl("note 50", txt, fixed = TRUE))
})
