# Work-dir browser parsing and the file-manager shell-out.

mk_task_log <- function(dir, lines) {
  writeLines(lines, file.path(dir, ".nextflow.log"))
  dir
}

th_line <- function(name, workdir, exit = NULL) {
  paste0("Jan-01 00:00:00.000 [Task monitor] DEBUG n.p.TaskPollingMonitor - ",
    "Task completed > TaskHandler[id: 1; name: ", name, "; status: COMPLETED; ",
    if (!is.null(exit)) paste0("exit: ", exit, "; ") else "",
    "error: -; workDir: ", workdir, "]")
}

test_that("a task with no parseable exit code reports an unknown status", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  mk_task_log(dir, c(
    th_line("preprocess (S1)", file.path(dir, "wd1"), exit = "0"),
    th_line("assemble (S1)", file.path(dir, "wd2"), exit = "137"),
    th_line("annotate (S1)", file.path(dir, "wd3"))
  ))
  d <- RiboPilot:::find_workdirs(dir, "S1")
  expect_equal(d$process, c("preprocess", "assemble", "annotate"))
  expect_equal(d$status, c("success", "failed", "unknown"))
  expect_false(any(is.na(d$status)))
})

test_that("a work dir containing an 's' is parsed whole", {
  skip_if_not_installed("withr")
  # TRE reads \s inside a bracket expression as a literal backslash and 's', so
  # the old [^]\s] cut every path at its first "s". Paths then collided and the
  # de-duplication collapsed the rows down to one.
  dir <- withr::local_tempdir()
  wd <- file.path(dir, "scaffolds", "assemble_step")
  mk_task_log(dir, c(
    th_line("assemble (S1)", wd, exit = "0"),
    th_line("screen (S1)", file.path(dir, "screen_step"), exit = "0")))
  d <- RiboPilot:::find_workdirs(dir, "S1")
  expect_equal(nrow(d), 2)
  expect_equal(d$workdir[1], wd)
  expect_false(any(grepl("^/tmp/[^/]*$", d$workdir)))
})

test_that("open_path ignores an absent, empty or NA path", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  bin <- withr::local_tempdir()
  marker <- file.path(bin, "opened")
  opener <- if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open"
  writeLines(c("#!/bin/sh", paste0("echo \"$1\" >> ", shQuote(marker))),
    file.path(bin, opener))
  Sys.chmod(file.path(bin, opener), "0755")
  withr::local_path(bin, action = "prefix")

  RiboPilot:::open_path(NULL)
  RiboPilot:::open_path("")
  RiboPilot:::open_path(NA_character_)
  RiboPilot:::open_path(character(0))

  RiboPilot:::open_path("/tmp/ribopilot-open-path")
  for (i in 1:40) { if (file.exists(marker)) break; Sys.sleep(0.05) }
  expect_true(file.exists(marker))
  expect_equal(readLines(marker), "/tmp/ribopilot-open-path")
})
