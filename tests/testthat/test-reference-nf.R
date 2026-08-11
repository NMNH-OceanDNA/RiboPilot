test_that("write_clade_manifest writes the four columns Nextflow reads", {
  d <- withr::local_tempdir()
  p <- RiboPilot:::write_clade_manifest(file.path(d, "clades.tsv"),
    data.frame(ref_id = "actinopteri", taxid = 186623L, clade_name = "Actinopteri",
      query = '"Actinopteri"[Organism]', stringsAsFactors = FALSE))
  x <- utils::read.delim(p, stringsAsFactors = FALSE)
  expect_equal(names(x), c("ref_id", "taxid", "clade_name", "query"))
  expect_equal(x$taxid, 186623L)
  lines <- readLines(p)
  expect_equal(lines[1], "ref_id\ttaxid\tclade_name\tquery")
  expect_equal(lines[2], 'actinopteri\t186623\tActinopteri\t"Actinopteri"[Organism]')
})

test_that("read_count_files reads every count file and skips junk", {
  d <- withr::local_tempdir()
  writeLines("7898\t11750", file.path(d, "7898.count"))
  writeLines("186623\t9000", file.path(d, "186623.count"))
  writeLines("", file.path(d, "6656.count"))
  x <- RiboPilot:::read_count_files(d)
  expect_equal(nrow(x), 2L)
  expect_equal(x$n_records[x$taxid == 7898L], 11750L)
})

test_that("merge_counts replaces stale rows and keeps the rest", {
  d <- withr::local_tempdir()
  utils::write.table(data.frame(taxid = c(7898L, 33208L), n_records = c(10L, 500L)),
    file.path(d, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  out <- RiboPilot:::merge_counts(d, data.frame(taxid = 7898L, n_records = 11750L))
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 500L)
  expect_equal(nrow(RiboPilot:::read_counts_cache(d)), 2L)
})

test_that("nf_reference_command names the entry point and every param", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  a <- RiboPilot:::nf_reference_command(d, "COUNT", file.path(d, "clades.tsv"),
    "/refs", file.path(d, "out"))
  expect_equal(a[1], "run")
  expect_true("-entry" %in% a)
  expect_equal(a[which(a == "-entry") + 1], "COUNT")
  expect_equal(a[which(a == "--refDir") + 1], "/refs")
  expect_true("--clades" %in% a && "--outDir" %in% a)
})

test_that("refseq and the api key only appear when asked for", {
  d <- withr::local_tempdir()
  writeLines("process { executor = 'local' }", file.path(d, ".config"))
  plain <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o")
  expect_false("--refseq" %in% plain)
  expect_false("--ncbi_api_key" %in% plain)
  keyed <- RiboPilot:::nf_reference_command(d, "BUILD", "c.tsv", "/refs", "o",
    api_key = "abc", refseq = TRUE)
  expect_equal(keyed[which(keyed == "--ncbi_api_key") + 1], "abc")
  expect_equal(keyed[which(keyed == "--refseq") + 1], "true")
})

test_that("refDir is never bind-mounted, so a folder that does not exist yet cannot fail the run", {
  d <- withr::local_tempdir()
  ref_dir <- file.path(d, "refs")   # deliberately not created
  for (engine in c("docker", "singularity")) {
    writeLines(RiboPilot:::container_engine_block(engine), file.path(d, ".config"))
    a <- RiboPilot:::nf_reference_command(d, "COUNT", "c.tsv", ref_dir, "o")
    expect_false("-process.containerOptions" %in% a)
    expect_equal(a[which(a == "--refDir") + 1], ref_dir)
  }
})

test_that("failed_task_errors reports the stderr of a task the run ignored", {
  d <- withr::local_tempdir()
  work <- file.path(d, "work_ref")
  task <- file.path(work, "aa", "474b0c1122")
  dir.create(task, recursive = TRUE)
  writeLines(c("WARNING: skipping mount of /nope",
    "FATAL:   container creation failed: mount source /nope doesn't exist"),
    file.path(task, ".command.err"))
  log <- file.path(d, "count.log")
  writeLines(c("executor >  sge (1)",
    "[aa/474b0c] NOTE: Process `COUNT:count_clade (Mullidae)` terminated with an error exit status (255) -- Error is ignored"),
    log)
  errs <- RiboPilot:::failed_task_errors(log, work)
  expect_length(errs, 1L)
  expect_true(grepl("container creation failed", errs))
})

test_that("failed_task_errors is quiet when nothing failed", {
  d <- withr::local_tempdir()
  log <- file.path(d, "count.log")
  writeLines(c("executor >  local (1)", "[97/e0f118] COUNT:count_clade (X) | 1 of 1"), log)
  expect_length(RiboPilot:::failed_task_errors(log, file.path(d, "work_ref")), 0L)
  expect_length(RiboPilot:::failed_task_errors(file.path(d, "missing.log"), d), 0L)
})

test_that("run_clade_counts skips cached clades and merges what the run produced", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  utils::write.table(data.frame(taxid = 33208L, n_records = 900000L),
    file.path(ref, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  seen <- new.env()
  runner <- function(project, args, log = NULL) {
    seen$args <- args
    out <- args[which(args == "--outDir") + 1]
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines("7898\t11750", file.path(out, "7898.count"))
    0L
  }
  clades <- data.frame(ref_id = c("metazoa", "actinopterygii"), taxid = c(33208L, 7898L),
    clade_name = c("Metazoa", "Actinopterygii"),
    query = c('"Metazoa"[Organism]', '"Actinopterygii"[Organism]'),
    stringsAsFactors = FALSE)
  out <- RiboPilot:::run_clade_counts(d, ref, clades, runner = runner)
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
  expect_equal(out$n_records[out$taxid == 33208L], 900000L)
  # only the uncached clade was sent to the run
  man <- utils::read.delim(seen$args[which(seen$args == "--clades") + 1], stringsAsFactors = FALSE)
  expect_equal(man$taxid, 7898L)
})

test_that("a clade cached as NA is counted again, not skipped forever", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  # the pre-Nextflow counter wrote NA on a failed fetch; that is not an answer
  writeLines(c("taxid\tn_records", "7898\tNA"), file.path(ref, "counts.tsv"))
  runner <- function(project, args, log = NULL) {
    out <- args[which(args == "--outDir") + 1]
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines("7898\t11750", file.path(out, "7898.count"))
    0L
  }
  out <- RiboPilot:::run_clade_counts(d, ref,
    data.frame(ref_id = "a", taxid = 7898L, clade_name = "Actinopterygii",
      query = "q", stringsAsFactors = FALSE), runner = runner)
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
})

test_that("run_clade_counts writes a log the caller can quote", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  runner <- function(project, args, log = NULL) {
    RiboPilot:::append_run_log(log, "ERROR ~ Unknown config attribute")
    1L
  }
  out <- NULL
  suppressWarnings(out <- RiboPilot:::run_clade_counts(d, ref,
    data.frame(ref_id = "a", taxid = 7898L, clade_name = "Actinopterygii",
      query = "q", stringsAsFactors = FALSE), runner = runner))
  log <- attr(out, "log")
  expect_true(file.exists(log))
  expect_true(grepl("Unknown config attribute", RiboPilot:::log_tail(log, 3), fixed = TRUE))
})

test_that("run_clade_counts launches nothing when every clade is cached", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  utils::write.table(data.frame(taxid = 7898L, n_records = 11750L),
    file.path(ref, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  runner <- function(project, args, log = NULL) stop("must not launch")
  out <- RiboPilot:::run_clade_counts(d, ref,
    data.frame(ref_id = "a", taxid = 7898L, clade_name = "Actinopterygii",
      query = "q", stringsAsFactors = FALSE), runner = runner)
  expect_equal(out$n_records[out$taxid == 7898L], 11750L)
})

test_that("a failed run still merges whatever counts landed", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  runner <- function(project, args, log = NULL) {
    out <- args[which(args == "--outDir") + 1]
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    writeLines("7898\t11750", file.path(out, "7898.count"))
    1L
  }
  out <- NULL
  expect_warning(
    out <- RiboPilot:::run_clade_counts(d, ref,
      data.frame(ref_id = c("a", "b"), taxid = c(7898L, 6656L),
        clade_name = c("Actinopterygii", "Arthropoda"), query = c("q1", "q2"),
        stringsAsFactors = FALSE), runner = runner),
    "status")
  expect_equal(nrow(out), 1L)
})

test_that("run_clade_counts rejects an unvalidated project before deleting anything", {
  d <- withr::local_tempdir()
  ref <- file.path(d, "refs"); dir.create(ref)
  clades <- data.frame(ref_id = "a", taxid = 7898L, clade_name = "Actinopterygii",
    query = "q", stringsAsFactors = FALSE)
  runner <- function(project, args, log = NULL) stop("must not launch")

  expect_error(RiboPilot:::run_clade_counts("", ref, clades, runner = runner))
  expect_error(RiboPilot:::run_clade_counts(NA_character_, ref, clades, runner = runner))
  expect_error(RiboPilot:::run_clade_counts(file.path(d, "nope"), ref, clades, runner = runner))

  # nothing under the temp dir was touched
  expect_false(dir.exists(file.path(d, ".reference")))
  expect_false(file.exists("/.reference"))
})

test_that("the project's config key wins over the environment variable", {
  d <- withr::local_tempdir()
  withr::local_envvar(c(NCBI_API_KEY = "from_env"))
  expect_equal(RiboPilot:::project_ncbi_key(d), "from_env")  # no config yet
  writeLines(c("params {", "  ncbi_api_key = 'from_project'", "}"), file.path(d, ".config"))
  expect_equal(RiboPilot:::project_ncbi_key(d), "from_project")
})

test_that("an unset or placeholder config key falls back to the environment", {
  d <- withr::local_tempdir()
  withr::local_envvar(c(NCBI_API_KEY = "from_env"))
  writeLines(c("params {", "  ncbi_api_key = ''", "}"), file.path(d, ".config"))
  expect_equal(RiboPilot:::project_ncbi_key(d), "from_env")
  writeLines(c("params {", "  ncbi_api_key = 'NA'", "}"), file.path(d, ".config"))
  expect_equal(RiboPilot:::project_ncbi_key(d), "from_env")
})
