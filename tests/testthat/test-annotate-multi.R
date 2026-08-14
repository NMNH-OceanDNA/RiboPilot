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
  back <- RiboPilot:::flip_coords(flipped, c(ACC1 = 1000L))
  expect_equal(back$start, 801L)
  expect_equal(back$stop, 901L)
})

test_that("flip_coords leaves records it has no length for untouched", {
  d <- data.frame(seqid = "ACCX", element = "5.8S", start = 5L, stop = 10L,
    stringsAsFactors = FALSE)
  expect_equal(RiboPilot:::flip_coords(d, c(ACC1 = 1000L)), d)
})

test_that("records over the ITSx ceiling are excluded from the ITSx input only", {
  lens <- c(SMALL = 5000L, HUGE = 200000L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "SMALL")
})

test_that("itsx_eligible includes the exact ceiling and excludes one bp over", {
  max_bp <- RiboPilot:::ITSX_MAX_BP
  lens <- c(AT = max_bp, OVER = max_bp + 1L)
  expect_equal(RiboPilot:::itsx_eligible(lens), "AT")
})

test_that("reconciling a minus-strand record in the revcomp frame keeps all five elements", {
  # ACC1 is 5000 bp, stored in non-sense (minus-strand) orientation. barrnap
  # reports anchors in the INPUT frame, descending (28S low ... 18S high, all
  # strand "-"). ITSx was run on the reverse complement and reports its
  # regions in THAT frame (ascending, canonical), unflipped -- this is what
  # annotate_records must do: flip bar into the revcomp frame, reconcile
  # there against ITSx's native output, then flip the finished model back.
  L <- 5000L
  lens <- c(ACC1 = L)
  bar <- data.frame(
    seqid = "ACC1", element = c("28S", "5.8S", "18S"),
    start = c(1L, 3651L, 4001L), stop = c(3450L, 3800L, 5000L),
    strand = "-", model_cov = 100, stringsAsFactors = FALSE)
  its <- data.frame(
    seqid = "ACC1", element = c("ITS1", "5.8S", "ITS2"),
    start = c(1001L, 1201L, 1351L), stop = c(1200L, 1350L, 1550L),
    stringsAsFactors = FALSE)

  bar_rc <- RiboPilot:::flip_coords(bar, lens)
  out_rc <- RiboPilot:::reconcile_records(bar_rc, its)
  out <- RiboPilot:::flip_coords(out_rc, lens)

  expect_setequal(out$element, c("18S", "ITS1", "5.8S", "ITS2", "28S"))
  expect_true(all(out$start >= 1L & out$stop <= L))
  expect_true(all(out$stop >= out$start))
})
