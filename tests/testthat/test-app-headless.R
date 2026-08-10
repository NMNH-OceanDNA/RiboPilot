# Which sessions get the cluster job script instead of a local "Run from App".

test_that("a session with no browser of its own is headless", {
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "", PATH = ""))
  expect_true(RiboPilot:::detect_headless(FALSE))
})

test_that("RStudio Server is headless even though it is interactive", {
  # the case this exists for: an RStudio Server session on a cluster login node
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "server", PATH = ""))
  expect_true(RiboPilot:::detect_headless(TRUE))
})

test_that("a batch scheduler on PATH means a cluster session", {
  bin <- withr::local_tempdir()
  writeLines("#!/bin/sh\n", file.path(bin, "sbatch"))
  Sys.chmod(file.path(bin, "sbatch"), "0755")
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "", PATH = bin))
  expect_true(RiboPilot:::detect_headless(TRUE))
})

test_that("a plain desktop session is not headless", {
  withr::local_envvar(c(RSTUDIO_PROGRAM_MODE = "desktop", PATH = withr::local_tempdir()))
  expect_false(RiboPilot:::detect_headless(TRUE))
})

test_that("the headless run modal offers the job script instead of a local run", {
  html <- as.character(RiboPilot:::run_modal(3, headless = TRUE,
    sched = list(cmd = "qsub", name = "SGE/PBS"), script = "#!/bin/bash\necho hi"))
  expect_true(grepl("submit_job", html, fixed = TRUE))
  expect_true(grepl("submit_script_text", html, fixed = TRUE))
  expect_false(grepl("nf_start", html, fixed = TRUE))
})
