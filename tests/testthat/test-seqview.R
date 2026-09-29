test_that("the viewer payload carries genes, depth, and clipped-read share", {
  f <- data.frame(element = c("18S", "ITS1"), start = c(1, 11), stop = c(10, 18),
    strand = c("+", "-"), partial_start = c(1, 0), partial_stop = c(0, NA))
  p <- RiboPilot:::seqview_payload("acgtacgtacgtacgtacgt", f, depths = c(rep(20, 19), 5),
    clips = c(10, rep(0, 19)), selected = "ITS1")
  expect_equal(p$len, 20L)
  expect_equal(p$seq, "ACGTACGTACGTACGTACGT")
  expect_equal(p$selected, "ITS1")
  expect_equal(vapply(p$features, `[[`, "", "row"), c("18S", "ITS1"))
  expect_true(p$features[[1]]$partial5)
  expect_false(p$features[[2]]$partial3)
  expect_equal(p$features[[2]]$dir, "-")
  expect_equal(p$features[[1]]$color, "#0072B2")
  expect_equal(p$clip[1], 0.5)
  expect_true(is.na(p$clip[20]))
  expect_equal(p$cutoff, RiboPilot:::CHIMERA_CLIP_FRAC)
})

test_that("the viewer payload leaves out tracks it has no data for", {
  f <- data.frame(element = character(0), start = numeric(0), stop = numeric(0), strand = character(0))
  p <- RiboPilot:::seqview_payload("ACGT", f)
  expect_null(p$depth); expect_null(p$clip); expect_null(p$selected)
  expect_length(p$features, 0L)
})

test_that("both details windows show the viewer instead of the coverage plot", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  html <- paste(as.character(RiboPilot:::annotate_modal("S1", 1L, 1L)), collapse = " ")
  expect_match(html, "seqv_ann-canvas", fixed = TRUE)
  expect_false(grepl("coverage_plot", html, fixed = TRUE))
})
