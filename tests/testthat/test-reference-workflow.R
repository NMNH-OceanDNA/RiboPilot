nf_txt <- paste(readLines(
  system.file("nextflow", "main_reference.nf", package = "RiboPilot")), collapse = "\n")
wrapper_file <- system.file("nextflow", "bin", "rp_ref_db.R", package = "RiboPilot")

test_that("the reference workflow ships with both entry points", {
  expect_true(nzchar(system.file("nextflow", "main_reference.nf", package = "RiboPilot")))
  expect_true(grepl("workflow COUNT", nf_txt, fixed = TRUE))
  expect_true(grepl("workflow BUILD", nf_txt, fixed = TRUE))
  # the work list is read as a header TSV, matching write_clade_manifest()
  expect_true(grepl("splitCsv(header: true, sep: '\\t')", nf_txt, fixed = TRUE))
})

test_that("the workflow declares exactly the params the launcher passes", {
  proj <- withr::local_tempdir()
  args <- nf_reference_command(proj, "COUNT", file.path(proj, "clades.tsv"),
    file.path(proj, "ref"), file.path(proj, "out"), api_key = "KEY", refseq = TRUE)
  passed <- sub("^--", "", grep("^--", args, value = TRUE))
  declared <- sub("^params\\.([A-Za-z_]+).*$", "\\1",
    grep("^params\\.", strsplit(nf_txt, "\n")[[1]], value = TRUE))
  expect_setequal(passed, declared)
})

test_that("the workflow publishes where reference_db_paths() looks", {
  p <- reference_db_paths("/refdir", "REFID")
  seeds_rel <- file.path(basename(dirname(p$seeds)),
    sub("^REFID_rrna", "*", basename(p$seeds)))
  labels_rel <- file.path(basename(dirname(p$labels)),
    sub("^REFID_rrna", "*", basename(p$labels)))
  expect_true(grepl(paste0("publishDir params.refDir, mode: 'copy', pattern: '",
    seeds_rel, "'"), nf_txt, fixed = TRUE))
  expect_true(grepl(paste0("publishDir params.refDir, mode: 'copy', pattern: '",
    labels_rel, "'"), nf_txt, fixed = TRUE))
  expect_true(grepl(paste0("path \"", seeds_rel, "\""), nf_txt, fixed = TRUE))
  expect_true(grepl(paste0("path \"", labels_rel, "\""), nf_txt, fixed = TRUE))
  # the _rrna suffix rule belongs to reference_db_paths(); nothing may restate it
  expect_false(grepl("_rrna", nf_txt, fixed = TRUE))
  expect_true(any(grepl("build_reference_db",
    readLines(wrapper_file), fixed = TRUE)))
})

test_that("the task wrapper detaches stdin from Entrez Direct", {
  expect_true(nzchar(wrapper_file))
  # every invocation, in the wrapper and in the builder it delegates to
  for (src in c(wrapper_file, test_path("..", "..", "R", "build_reference_db.R"))) {
    txt <- readLines(src)
    calls <- sum(grepl("system2(", txt, fixed = TRUE))
    detached <- sum(grepl("stdin = \"/dev/null\"", txt, fixed = TRUE))
    expect_gt(calls, 0)
    expect_equal(detached, calls)
  }
})

test_that("the wrapper reaches package internals the way a container will", {
  txt <- readLines(wrapper_file)
  # rp_env.R falls back to library(RiboPilot), which exports neither of these
  for (fn in c("app_sys", "parse_record_count")) {
    expect_false(fn %in% getNamespaceExports("RiboPilot"))
    hits <- grep(paste0("\\b", fn, "\\("), txt, value = TRUE)
    expect_gt(length(hits), 0)
    expect_true(all(grepl(paste0("RiboPilot:::", fn, "("), hits, fixed = TRUE)))
  }
})
