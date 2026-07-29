# coverage.depths is a per-position blob, several kb per scaffold. It is
# downsampled server-side before it reaches a table cell.

test_that("downsampling preserves the ends and bounds the length", {
  x <- as.numeric(1:1000)
  y <- RiboPilot:::downsample_depths(x, n = 100)
  expect_lte(length(y), 100L)
  expect_equal(y[1], 1)
  expect_equal(y[length(y)], 1000)
})

test_that("a short track is returned unchanged", {
  x <- c(3, 4, 5)
  expect_equal(RiboPilot:::downsample_depths(x, n = 100), x)
})

test_that("an empty or unparseable blob yields NA rather than erroring", {
  expect_true(is.na(RiboPilot:::depth_track_string(NA_character_)))
  expect_true(is.na(RiboPilot:::depth_track_string("")))
  expect_true(is.na(RiboPilot:::depth_track_string("not,a,number")))
})

test_that("a comma-separated blob round-trips into a bounded track string", {
  blob <- paste(1:500, collapse = ",")
  s <- RiboPilot:::depth_track_string(blob, n = 50)
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  expect_lte(length(parts), 50L)
  expect_equal(as.numeric(parts[1]), 1)
})

test_that("the sparkline renderer draws an SVG and never evaluates its input", {
  js <- as.character(RiboPilot:::rt_depth_sparkline())
  expect_match(js, "svg", fixed = TRUE)
  expect_match(js, "polyline", fixed = TRUE)
  expect_false(grepl("eval(", js, fixed = TRUE))
  expect_false(grepl("innerHTML", js, fixed = TRUE))
})
