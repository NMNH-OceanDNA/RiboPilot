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

# --- the spike-in half ---------------------------------------------------------

fx_spikes <- function() fx_where(fx_units(), function(x) identical(x$label$source, "spike"))

test_that("every spike-in unit gets its expected verdict and class", {
  # per unit, never in aggregate: an aggregate floor is what lets one whole
  # category silently stop working
  spikes <- fx_spikes()
  expect_gt(length(spikes), 0)
  for (nm in names(spikes)) {
    r <- fx_score(spikes[[nm]])
    lab <- spikes[[nm]]$label
    expect_equal(r$verdict, lab$expected_verdict, info = nm)
    if (!is.na(lab$expected_class)) expect_equal(r$class, lab$expected_class, info = nm)
  }
})

test_that("an uncertain verdict is never counted as catching a contaminant", {
  # counting uncertain as a catch would let the 0.6 gate be tuned to look good
  # by declining to answer
  catchable <- fx_where(fx_spikes(), function(x)
    identical(x$label$expected_verdict, "off_target"))
  expect_gt(length(catchable), 0)
  for (nm in names(catchable)) {
    expect_equal(fx_score(catchable[[nm]])$verdict, "off_target", info = nm)
  }
})

test_that("each cross-kingdom spike is caught with its own subclass", {
  u <- fx_units()
  for (pair in list(c("fungal_on_fish", "fungal"),
                    c("algal_on_fish", "algal_plant"),
                    c("arthropod_on_fish", "metazoan_other"))) {
    r <- fx_score(u[[pair[1]]])
    expect_equal(r$verdict, "off_target", info = pair[1])
    expect_equal(r$class, pair[2], info = pair[1])
  }
})

test_that("a bacterial contaminant is caught, not reported as unjudged", {
  # doc 10 classifies a 16S title with no eukaryotic marker as `other` BY DESIGN.
  # That class is real evidence about what the scaffold is, so it votes: the unit
  # used to have zero rows left and came back no_reference, i.e. UNJUDGED.
  x <- fx_units()[["bacterial_on_fish"]]
  r <- fx_score(x)
  expect_true(all(x$hits$hit_class == "other"))
  expect_equal(r$verdict, "off_target")
  expect_equal(r$host_concordant, 0L)
  expect_equal(r$hit_class, "other")
  expect_match(r$reason, "non-rDNA records", fixed = TRUE)
  # NCBI's 2024 reclassification puts a bacterial kingdom (Pseudomonadati) in the
  # kingdom slot rather than taxid 2, so both generations must classify.
  expect_equal(r$class, "bacterial")
})

test_that("the same human distribution is missed on a fish host, caught on an invertebrate", {
  # holds the contaminant constant and varies ONLY the host, so the screen's
  # boundary is measured rather than argued. Homo sapiens is Chordata, and so is
  # a ray-finned fish; concordance is tested at phylum rank or deeper.
  u <- fx_units()
  expect_true(all(c("human_on_fish", "human_on_invert") %in% names(u)))
  same <- u[["human_on_fish"]]; diff <- u[["human_on_invert"]]
  expect_equal(nrow(same$hits), nrow(diff$hits))       # identical evidence
  expect_equal(same$hits$bitscore, diff$hits$bitscore)
  expect_false(identical(same$label$sample_taxids, diff$label$sample_taxids))
  expect_equal(fx_score(same)$verdict, "on_target")    # different host,
  expect_equal(fx_score(diff)$verdict, "off_target")   # different answer
})

test_that("the confusion matrix is reported", {
  u <- fx_units()
  got <- vapply(u, function(x) fx_score(x)$verdict, character(1))
  want <- vapply(u, function(x) x$label$expected_verdict, character(1))
  cat("\n--- screen validation ---\n")
  print(table(expected = want, got = got))
  n_host <- sum(vapply(u, function(x) identical(x$label$category, "host"), logical(1)))
  fp <- sum(want == "on_target" & got == "off_target")
  catch <- want == "off_target"
  cat(sprintf("host units %d | false positives %d | recall %d/%d\n",
    n_host, fp, sum(catch & got == "off_target"), sum(catch)))
  expect_equal(fp, 0L)
  expect_equal(sum(catch & got == "off_target"), sum(catch))
})

test_that("a resolvable but weakly supported contaminant is uncertain, not a catch", {
  # the 0.6 gate is otherwise unexercised: every captured unit sits either well
  # above it or at confidence 0. A contaminant the screen declines to call must
  # NOT be counted as caught, or the gate could be tuned to look good by
  # answering less often.
  x <- fx_units()[["fungal_weak_on_fish"]]
  r <- fx_score(x)
  expect_equal(r$verdict, "uncertain")
  expect_lt(r$confidence, 0.6)
  expect_equal(r$host_concordant, 0L)   # the evidence still points off-target
  expect_false(is.na(r$lca_taxid))      # and is still recorded
  # and it is excluded from the recall numerator
  catchable <- fx_where(fx_spikes(), function(y)
    identical(y$label$expected_verdict, "off_target"))
  expect_false("fungal_weak_on_fish" %in% names(catchable))
})
