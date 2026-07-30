# The validation gate of doc 08 section 5, made real. A labelled corpus the
# screen's verdicts are MEASURED against rather than asserted about.
#
# Offline and deterministic by construction: this file reads two committed
# fixtures and nothing else. It never runs blastn, opens a socket, or touches a
# project DB. The capture that produced the fixtures is dev/make_screen_fixture.R,
# which is dev-only and absent from the built package.

fx_hits <- read.csv(gzfile(testthat::test_path("fixtures", "screen_hits.csv.gz")),
  stringsAsFactors = FALSE)
fx_labels <- read.csv(testthat::test_path("fixtures", "screen_labels.csv"),
  stringsAsFactors = FALSE)

fx_units <- function() {
  stats::setNames(lapply(fx_labels$unit, function(u) list(
    hits = fx_hits[fx_hits$unit == u, , drop = FALSE],
    label = fx_labels[fx_labels$unit == u, , drop = FALSE])), fx_labels$unit)
}

# score one fixture unit exactly as the pipeline would
fx_score <- function(x) {
  RiboPilot:::screen_channel_taxonomy(x$hits, x$label$sample_taxids)
}

fx_where <- function(u, f) u[vapply(u, f, logical(1))]

test_that("the fixture is internally consistent", {
  # a label with no hits, or hits with no label, is a corrupt fixture and must
  # fail loudly rather than be silently skipped
  expect_setequal(unique(fx_hits$unit), fx_labels$unit)
  expect_equal(anyDuplicated(fx_labels$unit), 0L)
  expect_false(anyNA(fx_labels$sample_taxids))
  expect_true(all(nzchar(fx_labels$sample_taxids)))
  expect_true(all(fx_labels$expected_verdict %in% RiboPilot:::SCREEN_VERDICTS))
  ec <- fx_labels$expected_class[!is.na(fx_labels$expected_class)]
  expect_true(all(ec %in% RiboPilot:::SCREEN_CLASSES))
  # every unit carries at least one hit row, or it scores no_reference for a
  # reason that has nothing to do with the screen
  expect_true(all(table(fx_hits$unit) > 0))
  expect_true(all(fx_labels$unit %in% names(table(fx_hits$unit))))
})

test_that("the clean-host corpus yields ZERO false positives", {
  # doc 08 section 5: a false positive on an under-represented taxon is the
  # failure that makes the tool useless. This is the assertion the gate is about.
  host <- fx_where(fx_units(), function(x) identical(x$label$category, "host"))
  expect_gt(length(host), 0)
  verdicts <- vapply(host, function(x) fx_score(x)$verdict, character(1))
  expect_equal(names(verdicts)[verdicts == "off_target"], character(0))
})

test_that("every host unit gets its expected verdict", {
  host <- fx_where(fx_units(), function(x) identical(x$label$category, "host"))
  for (nm in names(host)) {
    expect_equal(fx_score(host[[nm]])$verdict,
      host[[nm]]$label$expected_verdict, info = nm)
  }
})
