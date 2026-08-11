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
  expect_match(s, "#\\$ -o /scratch/me/proj/ribopilot\\.\\$JOB_ID\\.log")
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
  # Hydra has no "mpi" parallel environment, and neither do most SGE sites
  expect_match(s, "#\\$ -pe mthread 2")
  expect_match(s, "#\\$ -S /bin/bash")
})
