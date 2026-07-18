# Small pure helpers: alignment coverage, N-run detection, DB journal choice,
# constants coherence.

test_that("alignment_coverage is aligned-sample-bases / element length", {
  expect_equal(RiboPilot:::alignment_coverage("ACGT", 4), 1)
  expect_equal(RiboPilot:::alignment_coverage("AC--GT", 6), 4 / 6)   # gaps not counted
  expect_lte(RiboPilot:::alignment_coverage("ACGTACGT", 4), 1)        # capped at 1
  expect_true(is.na(RiboPilot:::alignment_coverage("ACGT", 0)))
})

test_that("element constants are internally consistent", {
  expect_true(all(RiboPilot:::RRNA_GENES %in% RiboPilot:::RRNA_ELEMENT_ORDER))
  expect_true(all(RiboPilot:::RRNA_SPACERS %in% RiboPilot:::RRNA_ELEMENT_ORDER))
  expect_equal(RiboPilot:::rrna_element_type("18S"), "rRNA")
  expect_equal(RiboPilot:::rrna_element_type("ITS1"), "misc_RNA")
  expect_equal(RiboPilot:::rrna_element_product("ITS2"), "internal transcribed spacer 2")
})

test_that("sqlite journal mode is overridable", {
  old <- Sys.getenv("RIBOPILOT_SQLITE_JOURNAL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("RIBOPILOT_SQLITE_JOURNAL") else
    Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = old))
  Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = "delete")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/x.sqlite"), "delete")
  Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = "wal")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/x.sqlite"), "wal")
})
