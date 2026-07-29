# Nextflow version classification + NXF_VER pin selection (pure helpers).

test_that("nf_version_status classifies against the supported range", {
  f <- RiboPilot:::nf_version_status
  expect_equal(f(NA_character_), "missing")
  expect_equal(f("23.10.0"), "too_old")
  expect_equal(f("24.04.4"), "too_old")
  expect_equal(f("24.10.0"), "supported")   # lower bound inclusive
  expect_equal(f("25.04.8"), "supported")
  expect_equal(f("25.10.6"), "supported")
  expect_equal(f("26.0.0"), "too_new")      # broken-from inclusive
  expect_equal(f("26.1.0"), "too_new")
})

test_that("nf_pin_for pins in-range to itself and out-of-range to the max", {
  f <- RiboPilot:::nf_pin_for
  expect_true(is.na(f(NA_character_)))
  # supported -> preserve the original string, including a leading zero
  expect_equal(f("25.04.8"), "25.04.8")
  expect_equal(f("24.10.0"), "24.10.0")
  # out of range -> latest verified-good
  expect_equal(f("23.10.0"), RiboPilot:::NF_MAX_SUPPORTED)
  expect_equal(f("26.1.0"), RiboPilot:::NF_MAX_SUPPORTED)
})

test_that("nxf_ver_prefix is empty when Nextflow is missing", {
  # nf_pin_version(<bogus binary>) -> NA -> empty prefix, so nothing is pinned.
  expect_equal(RiboPilot:::nxf_ver_prefix("definitely-not-a-real-nextflow-bin"), "")
})

# WF1 pipeline wiring: the screen must run per sample and triage must wait for
# it, so a sample is triaged only once both its curation and its screen exist.

main_nf <- readLines(system.file("nextflow", "main.nf", package = "RiboPilot"))
pipeline_r_nf <- readLines(system.file("nextflow", "modules", "pipeline_r.nf", package = "RiboPilot"))

test_that("main.nf includes screen from pipeline_r.nf", {
  include_line <- grep("^include \\{.*\\} from './modules/pipeline_r.nf'$", main_nf, value = TRUE)
  expect_length(include_line, 1)
  expect_match(include_line, "\\bscreen\\b")
})

test_that("screen process exists in pipeline_r.nf with errorStrategy 'ignore'", {
  start <- grep("^process screen \\{", pipeline_r_nf)
  expect_length(start, 1)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]
  expect_true(any(grepl("errorStrategy 'ignore'", body)))
  expect_true(any(grepl("cache false", body)))
})

test_that("triage's input joins on screen.out", {
  start <- grep("^\\s*triage\\(", main_nf)
  expect_length(start, 1)
  end <- start - 1 + which(grepl("^\\s*\\}\\s*$", main_nf[start:length(main_nf)]) |
    grepl("\\)\\)$", main_nf[start:length(main_nf)]))[1]
  block <- paste(main_nf[start:end], collapse = "\n")
  expect_match(block, "\\.join\\(screen\\.out")
})
