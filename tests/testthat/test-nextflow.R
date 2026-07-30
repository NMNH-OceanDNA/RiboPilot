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
  # screen carries errorStrategy 'ignore', so an inner join would silently drop
  # a fully curated sample whenever its screen (or its upstream BLAST) failed.
  # remainder: true is what keeps that sample in triage; guard it specifically,
  # not just the presence of a join.
  expect_match(block, "\\.join\\(screen\\.out\\.map \\{ id -> tuple\\(id, true\\) \\}, remainder: true\\)")
})

# The assembly fallback is a per-sample assembler inside WF1, not a second
# workflow. It used to be a second workflow, and the copy drifted out of step
# with main.nf until it no longer ran. These assertions read the workflow files
# as text because the defect that motivated them was a channel arity mismatch,
# which no R-level test can see.

test_that("query_samples selects the assembler and declares it in the header", {
  start <- grep("^process query_samples \\{", pipeline_r_nf)
  expect_length(start, 1)
  end <- start - 1 + which(pipeline_r_nf[start:length(pipeline_r_nf)] == "}")[1]
  body <- pipeline_r_nf[start:end]

  hdr <- grep("printf 'ID", body, value = TRUE)
  expect_length(hdr, 1)
  expect_match(hdr, "assembler")
  expect_true(any(grepl("AS assembler", body)))

  # The header is hand-written because sqlite3 -header emits nothing for zero
  # rows, so it and the select list must agree on column count or splitCsv
  # misreads every field.
  n_hdr <- length(strsplit(sub(".*printf '", "", sub("\\\\n' > .*", "", hdr)), "\\\\t")[[1]])
  expect_equal(n_hdr, 12L)
})

test_that("main.nf branches on the assembler and merges both arms into score", {
  br <- grep("asm_branched = asm_in\\.branch", main_nf)
  expect_length(br, 1)
  block <- paste(main_nf[br:(br + 4)], collapse = "\n")
  expect_match(block, "go: it\\[7\\] == 'getorganelle'")
  # `fb: true`, not an equality test: branch silently DROPS an item matching no
  # condition, so an unrecognised assembler would vanish from the run.
  expect_match(block, "fb: true")

  expect_true(any(grepl("^\\s*assemble\\(asm_branched\\.go", main_nf)))
  expect_true(any(grepl("^\\s*fallback_assemble\\(asm_branched\\.fb", main_nf)))
  expect_true(any(grepl("score\\(assemble\\.out\\.paths\\.mix\\(fallback_assemble\\.out\\.paths\\)\\)", main_nf)))
})

test_that("main.nf includes fallback_assemble and still passes seeds to GetOrganelle", {
  expect_true(any(grepl("^include \\{ fallback_assemble \\} from './modules/fallback.nf'$", main_nf)))
  go_call <- grep("assemble\\(asm_branched\\.go", main_nf)
  expect_length(go_call, 1)
  block <- paste(main_nf[go_call:(go_call + 1)], collapse = " ")
  expect_match(block, "tuple\\(id, reads, cpus, mem, go, seeds, labels\\)")
})

test_that("there is no second assembly workflow left to drift", {
  nf_dir <- system.file("nextflow", package = "RiboPilot")
  wf <- setdiff(list.files(nf_dir, pattern = "\\.nf$"), "nextflow.config")
  expect_setequal(wf, c("main.nf", "main_annotate.nf"))
  expect_false(file.exists(file.path(nf_dir, "main_fallback.nf")))

  # query_fallback_samples went with the workflow that called it.
  expect_false(any(grepl("query_fallback_samples", pipeline_r_nf)))

  # The stitch machinery survives pending doc 08 open question 2, but nothing
  # reachable calls it any more: no workflow includes or invokes stitch.
  annotate_nf <- readLines(system.file("nextflow", "main_annotate.nf", package = "RiboPilot"))
  for (nm in list(main = main_nf, annotate = annotate_nf)) {
    expect_false(any(grepl("\\bstitch\\b", nm[!grepl("^\\s*//", nm)])))
  }
})
