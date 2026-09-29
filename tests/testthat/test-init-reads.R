# A read file named in the mapping file but absent on disk used to surface only
# per-task, deep inside a work/ directory, as an opaque fastp failure.

# Writes a mapping file plus the read files it names, and returns the paths.
mk_reads_project <- function(ids = "S1", present = TRUE, R2 = NULL) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  data_dir <- file.path(dir, "data")
  dir.create(data_dir)
  mapping <- data.frame(
    ID = ids,
    Taxon = "Danio rerio",
    R1 = paste0(ids, "_R1.fq.gz"),
    R2 = if (is.null(R2)) paste0(ids, "_R2.fq.gz") else R2,
    stringsAsFactors = FALSE
  )
  keep <- rep_len(present, length(ids))
  for (i in which(keep)) {
    reads <- c(mapping$R1[i], mapping$R2[i])
    reads <- reads[!is.na(reads) & nzchar(reads)]
    file.create(file.path(data_dir, reads))
  }
  map_fn <- file.path(dir, "mapping.csv")
  utils::write.csv(mapping, map_fn, row.names = FALSE)
  list(dir = dir, data_dir = data_dir, mapping_fn = map_fn, db = file.path(dir, ".sqlite"))
}

# any existing file passes the seeds/labels check
mk_reads_db <- function(p, ...) {
  new_db(db_path = p$db, mapping_fn = p$mapping_fn, data_path = p$data_dir,
    seeds_db = p$mapping_fn, labels_db = p$mapping_fn, ...)
}

test_that("new_db accepts a mapping whose read files are all present", {
  p <- mk_reads_project(c("S1", "S2"))
  expect_no_error(mk_reads_db(p))
})

test_that("a missing R1 errors and names the sample and the path", {
  p <- mk_reads_project("S1")
  file.remove(file.path(p$data_dir, "S1_R1.fq.gz"))
  expect_error(mk_reads_db(p), "S1: .*S1_R1\\.fq\\.gz")
  expect_error(mk_reads_db(p), "not found")
})

test_that("a missing R2 errors and names the sample and the path", {
  p <- mk_reads_project("S1")
  file.remove(file.path(p$data_dir, "S1_R2.fq.gz"))
  expect_error(mk_reads_db(p), "S1: .*S1_R2\\.fq\\.gz")
})

test_that("several missing files are reported in one error", {
  p <- mk_reads_project(c("S1", "S2", "S3"), present = c(TRUE, FALSE, FALSE))
  err <- expect_error(mk_reads_db(p))
  expect_match(conditionMessage(err), "S2_R1\\.fq\\.gz")
  expect_match(conditionMessage(err), "S2_R2\\.fq\\.gz")
  expect_match(conditionMessage(err), "S3_R1\\.fq\\.gz")
  expect_false(grepl("S1_R1\\.fq\\.gz", conditionMessage(err)))
})

test_that("a long list of missing files is truncated", {
  ids <- sprintf("S%02d", 1:20)
  p <- mk_reads_project(ids, present = FALSE)
  err <- expect_error(mk_reads_db(p))
  msg <- conditionMessage(err)
  expect_match(msg, "\\.\\.\\. and 30 more")
  # header + 10 entries + the truncation line
  expect_equal(length(strsplit(msg, "\n")[[1]]), 12L)
})

test_that("an empty or NA R2 reports unsupported single-end input", {
  for (r2 in list("", NA_character_)) {
    p <- mk_reads_project("S1", R2 = r2)
    err <- expect_error(mk_reads_db(p))
    expect_match(conditionMessage(err), "Single-end reads are not supported")
    expect_match(conditionMessage(err), "S1")
  }
})

test_that("check_reads = FALSE skips the check", {
  p <- mk_reads_project("S1", present = FALSE)
  expect_no_error(mk_reads_db(p, check_reads = FALSE))
})

test_that("new_project checks the reads and writes a private .config", {
  skip_on_os("windows")
  p <- mk_reads_project("S1")
  proj <- withr::local_tempdir()

  expect_error(
    suppressWarnings(suppressMessages(new_project(path = proj, mapping_fn = p$mapping_fn,
      data_path = p$dir, executor = "conda", Rproj = FALSE,
      custom_seeds_db = p$mapping_fn, custom_labels_db = p$mapping_fn,
      ncbi_api_key = "secret"))),
    "not found"
  )

  expect_warning(suppressMessages(new_project(path = proj, mapping_fn = p$mapping_fn,
    data_path = p$data_dir, executor = "conda", Rproj = FALSE, force = TRUE,
    custom_seeds_db = p$mapping_fn, custom_labels_db = p$mapping_fn,
    ncbi_api_key = "secret")), "not supported")

  cfg <- file.path(proj, ".config")
  expect_true(file.exists(cfg))
  expect_equal(as.character(file.info(cfg)$mode), "600")
  # capped at this machine, so an oversized request cannot end the run
  lim <- grep("resourceLimits", readLines(cfg), value = TRUE)
  expect_match(lim, sprintf("cpus: %d, memory: %d.GB", ps::ps_cpu_count(),
    as.integer(floor(ps::ps_system_memory()$total / 1024^3))), fixed = TRUE)
})

test_that("new_project runs outside RStudio with the default Rproj = TRUE", {
  skip_if(requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  p <- mk_reads_project("S1")
  proj <- withr::local_tempdir()
  expect_no_error(suppressWarnings(suppressMessages(new_project(path = proj,
    mapping_fn = p$mapping_fn, data_path = p$data_dir, executor = "conda",
    custom_seeds_db = p$mapping_fn, custom_labels_db = p$mapping_fn))))
  expect_true(file.exists(file.path(proj, ".config")))
})
