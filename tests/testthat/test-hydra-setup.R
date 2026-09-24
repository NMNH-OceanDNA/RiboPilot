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

test_that("nothing happens once the scheduler is visible", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() list(cmd = "qsub", name = "SGE/PBS"),
    .package = "RiboPilot")
  before <- Sys.getenv("PATH")
  expect_silent(res <- RiboPilot:::ensure_hydra_setup())
  expect_false(res)
  expect_equal(Sys.getenv("PATH"), before)
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
  expect_match(s, "#\\$ -o /scratch/me/proj/logs/ribopilot\\.\\$JOB_ID\\.log")
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
    expect_match(txt, "-o /tmp/proj/logs/ribopilot", fixed = TRUE)
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
  base <- RiboPilot:::job_basename("ribopilot", as.POSIXct("2026-08-12 09:45:03", tz = "UTC"))
  expect_equal(base, "ribopilot_2026-08-12_09-45-03")

  txt <- RiboPilot:::build_submit_script(dir, sched = list(cmd = "qsub"), base = base)
  # the scheduler log, the Nextflow log, and the job name all carry the one base
  expect_match(txt, file.path(dir, "logs", paste0(base, ".$JOB_ID.log")), fixed = TRUE)
  expect_match(txt, file.path(dir, "nextflow_logs", paste0(base, ".nextflow.log")), fixed = TRUE)
  expect_match(txt, paste0("#$ -N ", base), fixed = TRUE)
  # the directories exist by the time the job runs: a scheduler will not create a
  # missing -o directory, it fails the job before anything can say why
  expect_true(dir.exists(file.path(dir, "logs")))
  expect_true(dir.exists(file.path(dir, "nextflow_logs")))
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
    base = RiboPilot:::job_basename("ribopilot", t1))
  b <- RiboPilot:::build_submit_script(dir, sched = list(cmd = "qsub"),
    base = RiboPilot:::job_basename("ribopilot", t1 + 60))
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
