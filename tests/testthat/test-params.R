# Clade ruleset dispatch, per-element identity floors, and lineage mapping.

test_that("rrna_params returns every known ruleset and errors on unknown", {
  for (r in RiboPilot:::RRNA_RULESETS) {
    p <- RiboPilot:::rrna_params(r)
    expect_equal(p$name, r)
    expect_true(is.list(p$lengths))
    expect_true(all(c("18S", "5.8S", "28S") %in% names(p$lengths)))
  }
  expect_error(RiboPilot:::rrna_params("octocorali"), "Unknown ruleset")
  expect_error(RiboPilot:::rrna_params("actinopterygii "), "Unknown ruleset")
  # NULL / empty fall back to generic rather than erroring
  expect_equal(RiboPilot:::rrna_params(NULL)$name, "metazoa_generic")
  expect_equal(RiboPilot:::rrna_params("")$name, "metazoa_generic")
})

test_that("invert 18S range is tighter than generic", {
  gen <- RiboPilot:::rrna_params("metazoa_generic")$lengths[["18S"]]
  inv <- RiboPilot:::rrna_params("octocorallia")$lengths[["18S"]]
  expect_true(inv[1] >= gen[1] && inv[2] <= gen[2])
})

test_that("min_anchor_pid is per-element and ordered 28S < 5.8S <= 18S", {
  p <- RiboPilot:::rrna_params("metazoa_generic")
  f18 <- RiboPilot:::anchor_pid_floor(p$min_anchor_pid, "18S")
  f58 <- RiboPilot:::anchor_pid_floor(p$min_anchor_pid, "5.8S")
  f28 <- RiboPilot:::anchor_pid_floor(p$min_anchor_pid, "28S")
  expect_true(f28 < f58)
  expect_true(f58 <= f18)
  # scalar back-compat
  expect_equal(RiboPilot:::anchor_pid_floor(85, "18S"), 85)
  expect_true(is.na(RiboPilot:::anchor_pid_floor(NULL, "18S")))
})

test_that("ruleset_for_lineage maps clades, incl. 2022 octocoral orders", {
  rl <- RiboPilot:::ruleset_for_lineage
  expect_equal(rl("Metazoa; Echinodermata; Asteroidea; Valvatida"), "asteroidea")
  expect_equal(rl("Metazoa; Cnidaria; Scyphozoa; Semaeostomeae"), "scyphozoa")
  expect_equal(rl("Metazoa; Cnidaria; Anthozoa; Octocorallia; Malacalcyonacea"), "octocorallia")
  # lineage with NO Octocorallia token, only the 2022 order name
  expect_equal(rl("Metazoa; Cnidaria; Anthozoa; Scleralcyonacea; Keratoisididae"), "octocorallia")
  expect_equal(rl("Metazoa; Cnidaria; Anthozoa; Alcyonacea"), "octocorallia")
  expect_equal(rl("Metazoa; Chordata; Actinopteri; Gadus"), "actinopterygii")
  expect_equal(rl("Metazoa; Mollusca; Gastropoda"), "metazoa_generic")
  expect_equal(rl(NA_character_), "metazoa_generic")
  expect_equal(rl(""), "metazoa_generic")
})
