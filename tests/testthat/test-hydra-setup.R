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

test_that("the app warns when a Hydra session has no scheduler on PATH", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() NULL,
    .package = "RiboPilot")
  expect_message(res <- RiboPilot:::warn_if_hydra_unconfigured(), "hydra_setup")
  expect_true(res)
})

test_that("no warning once the scheduler is visible", {
  testthat::local_mocked_bindings(
    is_hydra_cluster = function() TRUE,
    detect_scheduler = function() list(cmd = "qsub", name = "SGE/PBS"),
    .package = "RiboPilot")
  expect_silent(res <- RiboPilot:::warn_if_hydra_unconfigured())
  expect_false(res)
})
