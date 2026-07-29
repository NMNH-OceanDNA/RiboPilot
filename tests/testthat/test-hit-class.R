# Three-state classification of a GenBank defline. This layer answers "is this
# record nuclear rDNA", NOT "is this the host": a fungal or human rDNA record is
# correctly rdna, and catching it as a contaminant is the taxonomy channel's job.

fixtures <- read.csv(testthat::test_path("fixtures", "hit_titles.csv"),
  stringsAsFactors = FALSE, colClasses = "character")

test_that("every fixture title classifies as labelled", {
  got <- RiboPilot:::classify_hit_title(fixtures$title)
  bad <- which(got != fixtures$expected)
  expect_equal(got, fixtures$expected,
    info = paste("mismatched:", paste(fixtures$title[bad], collapse = " | ")))
})

test_that("the corpus actually covers all three classes and both sources", {
  # a corpus that drifted to one class would make the test above vacuous
  expect_setequal(unique(fixtures$expected), RiboPilot:::HIT_CLASSES)
  expect_setequal(unique(fixtures$source), c("real", "authored"))
  expect_gte(sum(fixtures$expected == "uninformative" & fixtures$source == "real"), 26L)
})

test_that("the organellar veto beats an rDNA phrase in the same title", {
  expect_equal(RiboPilot:::classify_hit_title(
    "Zea mays chloroplast large subunit ribosomal RNA gene, complete sequence"), "other")
  expect_equal(RiboPilot:::classify_hit_title(
    "Danio rerio mitochondrial 18S ribosomal RNA gene"), "other")
})

test_that("a prokaryotic record does not match the broadened rRNA gene pattern", {
  expect_equal(RiboPilot:::classify_hit_title(
    "Escherichia coli 16S ribosomal RNA gene, partial sequence"), "other")
  # but a record naming both keeps its eukaryotic reading
  expect_equal(RiboPilot:::classify_hit_title(
    "Synthetic construct 16S and 18S ribosomal RNA genes"), "rdna")
})

test_that("an assembly title naming an rDNA locus is rdna, not uninformative", {
  # precedence: the title says what the locus is, so the container is irrelevant
  expect_equal(RiboPilot:::classify_hit_title(
    "Danio rerio chromosome 5, 18S ribosomal RNA gene, complete sequence"), "rdna")
})

test_that("NA and empty titles are uninformative rather than an error", {
  expect_equal(RiboPilot:::classify_hit_title(NA_character_), "uninformative")
  expect_equal(RiboPilot:::classify_hit_title(""), "uninformative")
  expect_equal(RiboPilot:::classify_hit_title(character(0)), character(0))
})

test_that("consensus reads the informative hits and ignores the rest", {
  expect_equal(RiboPilot:::consensus_hit_class(
    c("rdna", "rdna", "rdna"), c(900, 890, 880)), 1L)
  expect_equal(RiboPilot:::consensus_hit_class(
    c("other", "other", "other"), c(900, 890, 880)), 0L)
  # the reference project's case: every top hit is a chromosome record
  expect_true(is.na(RiboPilot:::consensus_hit_class(
    c("uninformative", "uninformative", "uninformative"), c(900, 890, 880))))
})

test_that("one mislabelled top hit cannot flip the unit", {
  # the outlier holds the TOP bitscore, which is exactly the case hit #1 lost to
  cls <- c("other", "rdna", "rdna", "rdna")
  bs <- c(1000, 990, 985, 980)
  expect_equal(RiboPilot:::consensus_hit_class(cls, bs), 1L)
})

test_that("hits far below the top bitscore do not vote", {
  # 3 rdna near the top, a pile of weak other hits that must not outweigh them
  cls <- c(rep("rdna", 3), rep("other", 20))
  bs <- c(1000, 995, 990, rep(200, 20))
  expect_equal(RiboPilot:::consensus_hit_class(cls, bs), 1L)
})

test_that("a tie resolves to rdna", {
  # a false 0 costs a spurious triage flag; a false 1 only costs its absence
  expect_equal(RiboPilot:::consensus_hit_class(c("rdna", "other"), c(900, 900)), 1L)
})

test_that("fewer than three hits still vote, and empty input is NA", {
  expect_equal(RiboPilot:::consensus_hit_class("rdna", 900), 1L)
  expect_equal(RiboPilot:::consensus_hit_class(c("other", "rdna"), c(900, 100)), 0L)
  expect_true(is.na(RiboPilot:::consensus_hit_class(character(0), numeric(0))))
  expect_true(is.na(RiboPilot:::consensus_hit_class(NA_character_, NA_real_)))
})
