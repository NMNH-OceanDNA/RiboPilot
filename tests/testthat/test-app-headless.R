# Which sessions get the cluster job script instead of a local "Run from App".

test_that("not opening a browser does not turn on cluster mode", {
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "", PATH = ""))
  expect_false(RiboPilot:::detect_headless())
})

test_that("RStudio Server is headless even though it is interactive", {
  # the case this exists for: an RStudio Server session on a cluster login node
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "server", PATH = ""))
  expect_true(RiboPilot:::detect_headless())
})

test_that("a batch scheduler on PATH means a cluster session", {
  bin <- withr::local_tempdir()
  writeLines("#!/bin/sh\n", file.path(bin, "sbatch"))
  Sys.chmod(file.path(bin, "sbatch"), "0755")
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "", PATH = bin))
  expect_true(RiboPilot:::detect_headless())
})

test_that("a plain desktop session is not headless", {
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "desktop", PATH = withr::local_tempdir()))
  expect_false(RiboPilot:::detect_headless())
})

test_that("the headless run modal offers the job script instead of a local run", {
  html <- as.character(RiboPilot:::run_modal(3, headless = TRUE,
    sched = list(cmd = "qsub", name = "SGE/PBS"), script = "#!/bin/bash\necho hi"))
  expect_true(grepl("submit_job", html, fixed = TRUE))
  expect_true(grepl("submit_script_text", html, fixed = TRUE))
  expect_false(grepl("nf_start", html, fixed = TRUE))
  expect_true(grepl(".runs/nextflow/", html, fixed = TRUE))
  ref_html <- as.character(RiboPilot:::reference_run_modal(1, "/refs", headless = TRUE,
    sched = list(cmd = "sbatch", name = "Slurm"), script = "x"))
  expect_true(grepl(".runs/nextflow/", ref_html, fixed = TRUE))
})

test_that("a tunneled launch prints the local address and keeps Run from App", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir(); stamped_sqlite(dir)
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = ""))
  withr::local_options(list(RiboPilot.headless = NULL, RiboPilot.local_browser = NULL))
  testthat::local_mocked_bindings(check_nextflow_version = function(...) invisible("supported"),
    ensure_hydra_setup = function() invisible(FALSE), detect_scheduler = function() NULL,
    .package = "RiboPilot")
  msgs <- testthat::capture_messages(run_app(dir, port = 4321, launch.browser = FALSE))
  expect_true(any(grepl("Open http://localhost:4321", msgs, fixed = TRUE)))
  expect_false(any(grepl("Cluster mode", msgs, fixed = TRUE)))
  expect_false(getOption("RiboPilot.headless"))
  msgs <- testthat::capture_messages(run_app(dir, port = 4321, launch.browser = FALSE, headless = TRUE))
  expect_true(any(grepl("Pass headless = FALSE", msgs, fixed = TRUE)))
})

test_that("a folder with no project names new_project()", {
  expect_error(run_app(withr::local_tempdir()), "new_project()", fixed = TRUE)
})


test_that("the recheck window prefixes every id and has no Resume box", {
  has_id <- function(html, id) grepl(sprintf('id="%s"', id), html, fixed = TRUE)
  local_ids <- c("nf_cmd_text", "nf_progress", "nf_prog_all", "nf_start", "nf_stop",
    "nf_close", "rp_gears")
  html <- as.character(RiboPilot:::run_modal(2, pfx = "rc_",
    title = "RiboPilot recheck - 2 contig(s) to update", resume = FALSE))
  for (id in local_ids) {
    expect_true(has_id(html, paste0("rc_", id)), info = id)
    expect_false(has_id(html, id), info = id)
  }
  expect_true(grepl("RiboPilot recheck - 2 contig(s) to update", html, fixed = TRUE))
  expect_false(grepl("nf_resume", html, fixed = TRUE))
  hl <- as.character(RiboPilot:::run_modal(2, headless = TRUE,
    sched = list(cmd = "qsub", name = "SGE/PBS"), script = "x", pfx = "rc_", resume = FALSE))
  for (id in c("submit_script_text", "submit_job", "save_script", "nf_close", "nf_cmd_text")) {
    expect_true(has_id(hl, paste0("rc_", id)), info = id)
    expect_false(has_id(hl, id), info = id)
  }
  # the Assemble step's window keeps its ids, title, and Resume box
  def <- as.character(RiboPilot:::run_modal(3))
  for (id in c(local_ids, "nf_resume")) expect_true(has_id(def, id), info = id)
  expect_true(grepl("RiboPilot pipeline - 3 sample(s) or contig(s) to update", def, fixed = TRUE))
})
