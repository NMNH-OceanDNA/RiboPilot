# Intragenomic variants: are two on-target units the same molecule assembled
# twice, or two genuinely different rDNA copies from one individual. Restricting
# the comparison to on_target units is what makes any difference between them
# intragenomic BY CONSTRUCTION (doc 12 section 1), so no threshold has to
# separate a paralog from a contaminant. That separation already happened.
#
# Only the ANCHORS are compared. The spacers are what varies in both cases and
# so cannot discriminate between them.

ann_row <- function(element, start, stop, strand = "+") {
  data.frame(element = element, start = start, stop = stop, strand = strand,
    stringsAsFactors = FALSE)
}
anns <- function(...) do.call(rbind, list(...))

test_that("the anchor string is 18S + 5.8S + 28S in element order", {
  # 18S(1-10) ITS1(11-15) 5.8S(16-23) ITS2(24-27) 28S(28-39)
  seq <- paste0(strrep("A", 10), strrep("C", 5), strrep("G", 8), strrep("T", 4),
    strrep("A", 12))
  a <- anns(ann_row("18S", 1, 10), ann_row("ITS1", 11, 15),
            ann_row("5.8S", 16, 23), ann_row("ITS2", 24, 27),
            ann_row("28S", 28, 39))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8), strrep("A", 12)))
})

test_that("element order is fixed, not coordinate order", {
  # a circularly permuted assembly puts 28S first; the same molecule must not
  # produce a different anchor string and read as a paralog of itself
  seq <- paste0(strrep("G", 8), strrep("A", 10))
  a <- anns(ann_row("28S", 1, 8), ann_row("18S", 9, 18))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8)))
})

test_that("a unit with two of three anchors is still comparable", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10), ann_row("28S", 11, 18))
  expect_equal(nchar(RiboPilot:::screen_anchor_string(seq, a)), 18L)
})

test_that("a unit with no anchors at all is not comparable", {
  seq <- strrep("A", 30)
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, anns(ann_row("ITS1", 1, 10)))))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, data.frame())))
  expect_true(is.na(RiboPilot:::screen_anchor_string(NA_character_,
    anns(ann_row("18S", 1, 10)))))
})

test_that("an anchor with unusable coordinates is skipped, not fatal", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10), ann_row("5.8S", NA, NA), ann_row("28S", 11, 18))
  expect_equal(nchar(RiboPilot:::screen_anchor_string(seq, a)), 18L)
})

test_that("an all-minus unit is reverse-complemented rather than compared antisense", {
  # orient_sequence() normally leaves a unit on '+', so this is defensive: a unit
  # that escaped orientation must not be compared against an oriented one and
  # read as a paralog on the strength of being antisense
  fwd <- paste0(strrep("A", 10), strrep("G", 8))
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd)))
  a <- anns(ann_row("28S", 1, 8, "-"), ann_row("18S", 9, 18, "-"))
  expect_equal(RiboPilot:::screen_anchor_string(rc, a), fwd)
})

test_that("a mixed-strand unit is not comparable", {
  # concatenating a sense and an antisense span produces a string that means
  # nothing; a mixed-strand unit is a chimera signal, not a paralog
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10, "+"), ann_row("28S", 11, 18, "-"))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, a)))
})
