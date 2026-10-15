# Hydra session setup: the PATH fix that makes qsub visible to Nextflow.

test_that("hydra_setup warns and changes nothing off the cluster", {
  testthat::local_mocked_bindings(is_hydra_cluster = function() FALSE, .package = "RiboPilot")
  before <- Sys.getenv("PATH")
  expect_warning(res <- RiboPilot::hydra_setup(), "Not running on the NMNH Hydra")
  expect_false(res)
  expect_equal(Sys.getenv("PATH"), before)
})

test_that("hydra_setup puts the grid engine and Java at the front of PATH", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    check_nextflow_version = function(...) invisible("supported"),
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  withr::local_envvar(c(PATH = "/usr/bin:/bin", HOME = "/home/someone"))
  expect_message(res <- RiboPilot::hydra_setup(), "Hydra environment applied")
  expect_true(res)
  p <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  expect_equal(p[1:3], c("/home/someone/bin", RiboPilot:::HYDRA_BIN_DIRS))
  expect_true(all(c("/usr/bin", "/bin") %in% p))
})

test_that("a directory already on PATH is moved to the front, not duplicated", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    check_nextflow_version = function(...) invisible("supported"),
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  uge <- RiboPilot:::HYDRA_BIN_DIRS[2]
  withr::local_envvar(c(PATH = paste("/usr/bin", uge, sep = ":"), HOME = "/home/someone"))
  suppressMessages(RiboPilot::hydra_setup())
  p <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  expect_equal(sum(p == uge), 1L)
  expect_true(match(uge, p) < match("/usr/bin", p))
})

test_that("the app fixes PATH when a Hydra session has no scheduler on it", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() NULL,
    check_nextflow_version = function(...) invisible("supported"),
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  withr::local_envvar(c(PATH = "/usr/bin:/bin", HOME = "/home/someone"))
  expect_message(res <- RiboPilot:::ensure_hydra_setup(), "Hydra environment applied")
  expect_true(res)
  expect_true(RiboPilot:::HYDRA_BIN_DIRS[2] %in%
    strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]])
})

test_that("library() applies the Hydra PATH fix before probing Nextflow", {
  calls <- character(0)
  testthat::local_mocked_bindings(
    ensure_hydra_setup = function() { calls <<- c(calls, "hydra"); invisible(TRUE) },
    nf_installed_version = function(...) { calls <<- c(calls, "probe"); "25.04.0" },
    .package = "RiboPilot")
  rm(list = ls(RiboPilot:::.nf_probe_cache), envir = RiboPilot:::.nf_probe_cache)
  withr::defer(rm(list = ls(RiboPilot:::.nf_probe_cache), envir = RiboPilot:::.nf_probe_cache))
  suppressPackageStartupMessages(RiboPilot:::.onAttach("lib", "RiboPilot"))
  expect_equal(calls, c("hydra", "probe"))
})

test_that("nothing happens once the scheduler and nextflow are visible", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() list(cmd = "qsub", name = "SGE/PBS"),
    nf_on_path = function(...) TRUE,
    .package = "RiboPilot")
  before <- Sys.getenv("PATH")
  expect_silent(res <- RiboPilot:::ensure_hydra_setup())
  expect_false(res)
  expect_equal(Sys.getenv("PATH"), before)
})

test_that("the PATH fix still runs when qsub is visible but nextflow is not", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() list(cmd = "qsub", name = "SGE/PBS"),
    nf_on_path = function(...) FALSE,
    check_nextflow_version = function(...) invisible("supported"),
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  withr::local_envvar(c(PATH = "/usr/bin:/bin", HOME = "/home/someone"))
  expect_message(res <- RiboPilot:::ensure_hydra_setup(), "Hydra environment applied")
  expect_true(res)
  expect_true("/home/someone/bin" %in%
    strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]])
})

test_that("the Hydra head job goes to the workflow-manager queue", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  s <- RiboPilot:::job_script("/scratch/me/proj", "nextflow run x", NULL, "ribopilot")
  expect_match(s, "^#!/bin/sh")
  expect_match(s, "#\\$ -q lTWFM\\.sq", fixed = FALSE)
  expect_match(s, "#\\$ -l wfmq")
  expect_match(s, "#\\$ -pe mthread 1")
  expect_match(s, "mres=24G,h_data=24G,h_vmem=64G", fixed = TRUE)
  expect_match(s, "#\\$ -S /bin/sh")
  # absolute log path: the submitting session's cwd is not the project dir
  expect_match(s, "#\\$ -o /scratch/me/proj/\\.runs/logs/ribopilot\\.\\$JOB_ID\\.log")
  expect_match(s, "nextflow run x", fixed = TRUE)
})

test_that("off Hydra the generic SGE header is used", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  s <- RiboPilot:::job_script("/tmp/proj", "nextflow run x",
    list(cmd = "qsub", name = "SGE/PBS"), "ribopilot")
  expect_false(grepl("lTWFM.sq", s, fixed = TRUE))
  # no profile, so no parallel environment is guessed
  expect_false(grepl("-pe", s, fixed = TRUE))
  expect_match(s, "#$ -l h_rt=48:00:00", fixed = TRUE)
  expect_match(s, "#\\$ -S /bin/bash")
})

test_that("job scripts report how nextflow exited", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  s <- RiboPilot:::job_script("/scratch/me/proj", "nextflow run x", NULL, "ribopilot")
  # nextflow's own summary never reaches the job log under -ansi-log false, so a
  # finished run must still leave a definitive last line
  expect_match(s, "nextflow exited with status", fixed = TRUE)
  expect_match(s, "exit $status", fixed = TRUE)
  # -cwd would resolve to whatever directory the app was in, not the project
  expect_false(grepl("-cwd", s, fixed = TRUE))
})

test_that("generic job scripts write their log inside the project", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    .package = "RiboPilot")
  for (s in list(list(cmd = "qsub"), list(cmd = "sbatch"), list(cmd = "bsub"))) {
    txt <- RiboPilot:::job_script("/tmp/proj", "nextflow run x", s, "ribopilot")
    # logs go to the project's logs/ dir, not its root
    expect_match(txt, "-o /tmp/proj/.runs/logs/ribopilot", fixed = TRUE)
    expect_match(txt, "nextflow exited with status", fixed = TRUE)
  }
  expect_false(grepl("-cwd",
    RiboPilot:::job_script("/tmp/proj", "nextflow run x", list(cmd = "qsub"), "ribopilot"),
    fixed = TRUE))
})

test_that("one submission's script and both its logs share a timestamped name", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    nf_command = function(...) c("run", "main.nf"),
    .package = "RiboPilot")
  base <- RiboPilot:::job_basename("pipeline", as.POSIXct("2026-08-12 09:45:03", tz = "UTC"))
  expect_equal(base, "pipeline_2026-08-12_09-45-03")

  txt <- RiboPilot:::build_submit_script(dir, sched = list(cmd = "qsub"), base = base)
  # the scheduler log, the Nextflow log, and the job name all carry the one base
  expect_match(txt, file.path(dir, ".runs", "logs", paste0(base, ".$JOB_ID.log")), fixed = TRUE)
  expect_match(txt, file.path(dir, ".runs", "nextflow", paste0(base, ".nextflow.log")), fixed = TRUE)
  expect_match(txt, paste0("#$ -N ", base), fixed = TRUE)
  # the directories exist by the time the job runs: a scheduler will not create a
  # missing -o directory, it fails the job before anything can say why
  expect_true(dir.exists(file.path(dir, ".runs", "logs")))
  expect_true(dir.exists(file.path(dir, ".runs", "nextflow")))
})

test_that("a second submission does not overwrite the first one's log paths", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    nf_command = function(...) c("run", "main.nf"),
    .package = "RiboPilot")
  t1 <- as.POSIXct("2026-08-12 09:45:03", tz = "UTC")
  a <- RiboPilot:::build_submit_script(dir, sched = list(cmd = "qsub"),
    base = RiboPilot:::job_basename("pipeline", t1))
  b <- RiboPilot:::build_submit_script(dir, sched = list(cmd = "qsub"),
    base = RiboPilot:::job_basename("pipeline", t1 + 60))
  expect_false(identical(a, b))
})

test_that("the head job follows the project profile", {
  skip_if_not_installed("withr")
  testthat::local_mocked_bindings(is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_, .package = "RiboPilot")
  prof <- withr::local_tempdir(); proj <- withr::local_tempdir()
  job <- "ribopilot_2026-09-24_10-00-00"
  hdr <- function(s) {
    suppressMessages(generate_config(s, scheduler = s, queue = "q1", account = "a1",
      container_cache = "/c", profile_dir = prof, overwrite = TRUE))
    file.copy(file.path(prof, paste0("config.", s)), file.path(proj, ".config"), overwrite = TRUE)
    RiboPilot:::job_script(proj, "nextflow run x", NULL, job)
  }
  s <- hdr("slurm")
  expect_match(s, "#SBATCH -p q1", fixed = TRUE)
  expect_match(s, "#SBATCH --account=a1", fixed = TRUE)
  expect_match(s, paste0(job, "_%j.log"), fixed = TRUE)
  expect_false(grepl("%%", s, fixed = TRUE))
  s <- hdr("sge")
  expect_match(s, "#$ -pe mthread 2", fixed = TRUE)
  expect_match(s, "#$ -q q1", fixed = TRUE)
  expect_match(s, "#$ -V", fixed = TRUE)
  s <- hdr("pbs")
  expect_match(s, "#PBS -l select=1:ncpus=2:mem=8gb", fixed = TRUE)
  expect_match(s, "#PBS -l walltime=48:00:00", fixed = TRUE)
  expect_match(s, "#PBS -A a1", fixed = TRUE)
  expect_match(s, "#PBS -N ribopilot_2026-", fixed = TRUE)
  expect_false(grepl("#$", s, fixed = TRUE))
  s <- hdr("lsf")
  expect_match(s, "#BSUB -q q1", fixed = TRUE)
  expect_match(s, "#BSUB -W 48:00", fixed = TRUE)
  expect_match(s, paste0(job, ".%J.log"), fixed = TRUE)
  expect_false(grepl("%%", s, fixed = TRUE))
  for (x in c("slurm", "sge", "pbs", "lsf")) {
    expect_match(hdr(x), "export NXF_OPTS='-Xms512m -Xmx4g'", fixed = TRUE, info = x)
  }
})

test_that("a closure clusterOptions is left out of the head job", {
  skip_if_not_installed("withr")
  testthat::local_mocked_bindings(is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_, .package = "RiboPilot")
  proj <- withr::local_tempdir()
  writeLines(c("process {", "  executor = 'slurm'",
    "  clusterOptions = { \"--mem=${task.memory}\" }", "}"), file.path(proj, ".config"))
  s <- RiboPilot:::job_script(proj, "nextflow run x", NULL, "ribopilot")
  expect_false(grepl("{", s, fixed = TRUE))
})

test_that("the Hydra profile never asks SGE for h_rss", {
  skip_on_cran()
  skip_if_not_installed("withr")
  skip_if(Sys.which("nextflow") == "", "nextflow not installed")
  pin <- RiboPilot:::nf_pin_version()
  skip_if(is.na(pin), "no supported Nextflow")
  dir <- withr::local_tempdir()
  withr::local_dir(dir)
  dir.create("fake")
  writeLines(c("#!/bin/sh", "echo 'Your job 4242 (\"nf\") has been submitted'"), "fake/qsub")
  writeLines(c("#!/bin/sh", "true"), "fake/qstat")
  Sys.chmod(c("fake/qsub", "fake/qstat"), "755")
  writeLines(c(
    "process find_rrna_scan {", "  cpus { nthreads as int }",
    "  memory { \"${(mem as int) * task.attempt} GB\" }",
    "  input: tuple val(id), val(nthreads), val(mem)", "  script: 'true'", "}",
    "process find_rrna_pick {", "  cpus 1", "  memory { \"${4 * task.attempt} GB\" }",
    "  input: val(id)", "  script: 'true'", "}",
    "process stage_assembly {", "  cpus 1", "  memory 2.GB",
    "  input: val(id)", "  script: 'true'", "}",
    "workflow find_rrna {", "  take: x", "  main:", "  find_rrna_scan(x)",
    "  find_rrna_pick(x.map { it[0] })", "}",
    "workflow {", "  find_rrna(Channel.of(tuple('S1', 4, 12)))",
    "  stage_assembly(Channel.of('S1'))", "}"), "main.nf")
  cfg <- RiboPilot:::fill_config(readLines(RiboPilot:::app_sys("config.NMNH_Hydra")),
    list(CONTAINER_ID = "", RAW_DIR = "NA", ASMB_DIR = dir, MIN_DEPTH = "0",
      NCBI_API_KEY = "", EXTRA_MOUNTS = ""))
  writeLines(c(cfg, "singularity.enabled = false", "executor.pollInterval = '1 sec'"),
    "hydra.config")
  p <- processx::process$new("nextflow", c("-q", "run", "main.nf", "-c", "hydra.config",
      "-work-dir", "work", "--projectDir", dir),
    env = c("current", PATH = paste(file.path(dir, "fake"), Sys.getenv("PATH"), sep = ":"),
      NXF_VER = pin, NXF_DISABLE_CHECK_LATEST = "true"),
    stdout = NULL, stderr = NULL, cleanup_tree = TRUE)
  withr::defer(p$kill_tree())
  runs <- character(0)
  t0 <- Sys.time()
  while (length(runs) < 3 && difftime(Sys.time(), t0, units = "secs") < 120) {
    Sys.sleep(1)
    runs <- list.files("work", pattern = "^\\.command\\.run$", recursive = TRUE,
      all.files = TRUE, full.names = TRUE)
  }
  expect_length(runs, 3)
  hdr <- unlist(lapply(runs, function(f) grep("^#\\$ -l ", readLines(f), value = TRUE)))
  expect_false(any(grepl("h_rss|mem_free", hdr)))
  expect_true(any(grepl("mres=12G,h_data=3G", hdr, fixed = TRUE)))
  expect_equal(sum(grepl("mres=4G,h_data=4G", hdr, fixed = TRUE)), 2L)
})


test_that("a recheck submission's script and logs share one base name", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() FALSE,
    nf_pin_version = function(...) NA_character_,
    nf_recheck_command = function(path) c("run", "main_recheck.nf"),
    .package = "RiboPilot")
  base <- RiboPilot:::job_basename("recheck", as.POSIXct("2026-10-08 09:45:03", tz = "UTC"))
  expect_equal(base, "recheck_2026-10-08_09-45-03")
  txt <- RiboPilot:::build_recheck_submit_script(dir, sched = list(cmd = "qsub"), base = base)
  expect_match(txt, file.path(dir, ".runs", "logs", paste0(base, ".$JOB_ID.log")), fixed = TRUE)
  expect_match(txt, file.path(dir, ".runs", "nextflow", paste0(base, ".nextflow.log")), fixed = TRUE)
  expect_match(txt, paste0("#$ -N ", base), fixed = TRUE)
  expect_true(dir.exists(file.path(dir, ".runs", "logs")))
  expect_true(dir.exists(file.path(dir, ".runs", "nextflow")))
})
