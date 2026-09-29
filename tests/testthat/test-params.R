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
  expect_true(is.na(RiboPilot:::anchor_pid_floor(NULL, "18S")))
  expect_true(is.na(RiboPilot:::anchor_pid_floor(c("18S" = 90), "ITS1")))
})

test_that("every ruleset carries a screen sublist and its weights still sum to 1", {
  for (rs in RiboPilot:::RRNA_RULESETS) {
    p <- RiboPilot:::rrna_params(rs)
    expect_true(is.list(p$screen), info = rs)
    for (k in c("support", "resolution", "separation", "min_confidence",
                "min_isomer_anchor_pid")) {
      expect_true(!is.null(p$screen[[k]]), info = paste(rs, k))
    }
    expect_equal(sum(unlist(p$screen[c("support", "resolution", "separation")])), 1,
      info = rs)
  }
})

test_that("relocating the screen tunables changed none of their values", {
  # doc 11 section 5: the weights are priors. Building the sublist must not
  # silently retune them, or every stored verdict shifts for a reason nobody chose.
  p <- RiboPilot:::rrna_params("metazoa_generic")$screen
  expect_equal(p$support, 0.30)
  expect_equal(p$resolution, 0.30)
  expect_equal(p$separation, 0.40)
  expect_equal(p$min_confidence, 0.6)
  expect_equal(p$min_isomer_anchor_pid, 99.5)
})

test_that("a clade ruleset can retune the screen without touching the default", {
  base <- RiboPilot:::rrna_params("metazoa_generic")$screen
  for (rs in RiboPilot:::RRNA_RULESETS) {
    expect_equal(RiboPilot:::rrna_params(rs)$screen, base, info = rs)
  }
  # the mechanism exists even though no clade uses it yet: modifyList on the
  # sublist must not drop the untouched keys
  tuned <- utils::modifyList(base, list(min_isomer_anchor_pid = 99.9))
  expect_equal(tuned$support, 0.30)
  expect_equal(tuned$min_isomer_anchor_pid, 99.9)
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
