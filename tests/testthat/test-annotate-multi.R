# annotate_records() maps tool output back to the right record and the right
# coordinate frame. Tool output is injected rather than invoking real binaries.

test_that("annotate_records demultiplexes barrnap and ITSx output by record", {
  bar <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("18S", "18S"),
    start = c(10L, 20L), stop = c(1810L, 1920L),
    strand = c("+", "+"), model_cov = c(100, 100), stringsAsFactors = FALSE)
  its <- data.frame(
    seqid = c("ACC1", "ACC2"), element = c("ITS1", "ITS2"),
    start = c(1811L, 1921L), stop = c(2000L, 2100L), stringsAsFactors = FALSE)
  out <- RiboPilot:::reconcile_records(bar, its)
  expect_equal(sort(unique(out$seqid)), c("ACC1", "ACC2"))
  # the load-bearing property is that a record's coordinates come from its own
  # rows, not from whichever row sorted first
  expect_equal(out$start[out$seqid == "ACC1" & out$element == "18S"], 10L)
  expect_equal(out$start[out$seqid == "ACC2" & out$element == "18S"], 20L)
  expect_true("ITS1" %in% out$element[out$seqid == "ACC1"])
  expect_true("ITS2" %in% out$element[out$seqid == "ACC2"])
})

test_that("minus-strand ITSx coordinates are mapped back to the input frame", {
  # a 1000 bp record; ITSx searched the reverse complement and reported 100-200
  # there, which is 801-901 in the input frame
  flipped <- data.frame(seqid = "ACC1", element = "5.8S",
    start = 100L, stop = 200L, stringsAsFactors = FALSE)
  back <- RiboPilot:::flip_itsx_coords(flipped, c(ACC1 = 1000L))
  expect_equal(back$start, 801L)
  expect_equal(back$stop, 901L)
})

test_that("flip_itsx_coords leaves records it has no length for untouched", {
  d <- data.frame(seqid = "ACCX", element = "5.8S", start = 5L, stop = 10L,
    stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::flip_itsx_coords(d, c(ACC1 = 1000L)), d)
})

test_that("records over the ITSx ceiling are excluded from the ITSx input only", {
  lens <- c(SMALL = 5000L, HUGE = 200000L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "SMALL")
})
