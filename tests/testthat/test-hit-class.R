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
