# Saved run reports under .runs/, and the Run Reports window.

rr_log <- function(end = "Sep-08 13:00:57.110 [main] DEBUG nextflow.script.ScriptRunner - > Execution complete -- Goodbye") c(
  "Sep-08 12:51:07.607 [main] INFO  nextflow.cli.CmdRun - N E X T F L O W  ~  version 25.10.6",
  "Sep-08 12:51:09.134 [main] DEBUG nextflow.Session - Run name: suspicious_boltzmann",
  "Sep-08 12:51:13.321 [Task submitter] INFO  nextflow.Session - [a7/8fd59c] Submitted process > preprocess (S1)",
  "Sep-08 12:51:13.330 [Task submitter] INFO  nextflow.Session - [27/614924] Submitted process > annotate (S2.1.3)",
  "Sep-08 12:51:14.293 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > GridTaskHandler[jobId: 9; id: 11; name: annotate (S2.1.3); status: COMPLETED; exit: 0; error: -; workDir: /proj/work/27/614924]",
  "Sep-08 13:00:57.091 [main] DEBUG n.trace.WorkflowStatsObserver - Workflow completed > WorkflowStats[succeededCount=2; failedCount=0; ignoredCount=0; cachedCount=1; pendingCount=0; submittedCount=0; runningCount=0; retriesCount=0; abortedCount=0; ]",
  end)

put_log <- function(d, base, lines = rr_log()) {
  f <- file.path(RiboPilot:::run_dir(d, "nextflow"), paste0(base, ".nextflow.log"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, f)
  f
}

test_that("a run's base name carries its workflow and start time", {
  t <- as.POSIXct("2026-09-08 12:51:00")
  b <- RiboPilot:::job_basename("reference", t)
  expect_equal(RiboPilot:::run_workflow(b), "reference")
  expect_equal(RiboPilot:::run_base_time(b), t)
  expect_equal(RiboPilot:::run_log_base(paste0("/x/", b, ".nextflow.log")), b)
})

test_that("task tags map to samples, unit tags included", {
  expect_equal(RiboPilot:::tag_sample(c("S1", "S2.1.3", "S2.4")), c("S1", "S2", "S2"))
  expect_equal(RiboPilot:::tag_sample("A.B.1.2", ids = c("A", "A.B")), "A.B")
})

test_that("log stamps get the year of the log, not a future one", {
  ref <- as.POSIXct("2026-01-02 10:00:00")
  t <- RiboPilot:::nf_log_time(c("Dec-31 23:00:00.000", "Jan-02 09:00:00.000"), ref)
  expect_equal(format(t, "%Y"), c("2025", "2026"))
  expect_equal(RiboPilot:::fmt_duration(3723), "1h 2m 3s")
  expect_equal(RiboPilot:::fmt_duration(5), "5s")
})

test_that("finished runs get a saved report; unfinished ones are listed without one", {
  d <- withr::local_tempdir()
  put_log(d, "pipeline_2026-09-08_12-51-00")
  put_log(d, "pipeline_2026-09-09_08-00-00", rr_log(end = NULL))
  expect_equal(RiboPilot:::sync_run_reports(d), "pipeline_2026-09-08_12-51-00")
  expect_equal(RiboPilot:::sync_run_reports(d), character(0))
  txt <- readLines(RiboPilot:::run_report_paths(d, "pipeline_2026-09-08_12-51-00")[["txt"]])
  expect_true(any(grepl("Result:    Finished, exit status 0", txt, fixed = TRUE)))
  expect_true(any(grepl("Samples:   2", txt, fixed = TRUE)))

  runs <- RiboPilot:::list_run_reports(d)
  expect_equal(runs$base, c("pipeline_2026-09-09_08-00-00", "pipeline_2026-09-08_12-51-00"))
  expect_equal(runs$result, c("unfinished", "finished"))
  expect_equal(runs$n_samples[2], 2L)
  expect_equal(runs$cached[2], 1L)
  expect_true(is.na(runs$report[1]))
})

test_that("a job's report waits for its exit trailer and uses its status", {
  d <- withr::local_tempdir()
  base <- "reference_2026-09-08_12-51-00"
  put_log(d, base)
  job <- file.path(RiboPilot:::run_dir(d, "jobs"), paste0(base, ".sh"))
  dir.create(dirname(job), recursive = TRUE)
  writeLines(RiboPilot:::NF_EXIT_TRAILER, job)
  expect_equal(RiboPilot:::sync_run_reports(d), character(0))
  dir.create(RiboPilot:::run_dir(d, "logs"))
  writeLines("[RiboPilot] nextflow exited with status 1",
    file.path(RiboPilot:::run_dir(d, "logs"), paste0(base, "_42.log")))
  expect_equal(RiboPilot:::sync_run_reports(d), base)
  r <- RiboPilot:::list_run_reports(d)
  expect_equal(r$result, "failed")
  expect_equal(r$exit_status, 1L)
  expect_equal(r$launch, "job")
  expect_equal(r$workflow, "reference")
})

test_that("grid task handlers show up in the work dir browser", {
  d <- withr::local_tempdir()
  put_log(d, "pipeline_2026-09-08_12-51-00")
  w <- RiboPilot:::find_workdirs(d, "S2")
  expect_equal(w$workdir, "/proj/work/27/614924")
})

test_that("result cells carry the label", {
  h <- RiboPilot:::run_result_html(c("finished", "failed", NA))
  expect_match(h[1], "Finished")
  expect_match(h[2], "mp-pill-danger")
  expect_equal(h[3], "")
})

test_that("the Run Reports window writes, lists, and shows a finished run", {
  d <- withr::local_tempdir()
  base <- "pipeline_2026-09-08_12-51-00"
  put_log(d, base)
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    session$setInputs(rr_open = 1)
    expect_true(file.exists(RiboPilot:::run_report_paths(d, base)[["txt"]]))
    expect_equal(h$runs()$base, base)
    expect_match(output$rr_tbl, "New")
    session$setInputs(rr_tbl__reactable__selected = 1)
    expect_match(output$rr_detail$html, "RiboPilot run report")
  })
})

test_that("the Run Reports window copes with a project that has no runs", {
  d <- withr::local_tempdir()
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    expect_no_error(session$setInputs(refresh = 1))
    session$setInputs(rr_open = 1)
    expect_equal(nrow(h$runs()), 0L)
  })
})


test_that("Run Reports names a recheck run as such", {
  d <- withr::local_tempdir()
  put_log(d, "recheck_2026-09-08_12-51-00")
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    session$setInputs(rr_open = 1)
    expect_equal(h$runs()$workflow, "recheck")
    expect_match(output$rr_tbl, "Recheck", fixed = TRUE)
    expect_false(grepl("Pipeline", output$rr_tbl, fixed = TRUE))
  })
})

test_that("an open Run Reports window picks up a run that finishes", {
  d <- withr::local_tempdir()
  base <- "pipeline_2026-09-08_12-51-00"
  put_log(d, base, rr_log(end = NULL))
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    session$setInputs(rr_open = 1)
    session$setInputs(rr_tbl__reactable__selected = 1)
    output$rr_poll
    expect_match(output$rr_detail$html, "No report yet")
    expect_equal(h$runs()$result, "unfinished")
    put_log(d, base)
    session$elapse(15000)
    expect_equal(h$runs()$result, "finished")
    expect_match(output$rr_detail$html, "RiboPilot run report")
    expect_match(output$rr_tbl, "Finished")
  })
})

test_that("an open Run Reports window lists new runs every 15 seconds", {
  d <- withr::local_tempdir()
  shiny::testServer(function(input, output, session) {
    h <- RiboPilot:::run_reports_handlers(input, output, session, function() d,
      refresh = shiny::reactive(input$refresh))
  }, {
    session$setInputs(rr_open = 1)
    output$rr_poll
    expect_equal(nrow(h$runs()), 0L)
    put_log(d, "pipeline_2026-09-08_12-51-00")
    session$elapse(15000)
    expect_equal(nrow(h$runs()), 1L)
    output$rr_poll
    put_log(d, "pipeline_2026-09-09_12-51-00", rr_log(end = NULL))
    session$elapse(15000)
    expect_equal(h$runs()$result, c("unfinished", "finished"))
  })
})

test_that("a job's id is read from its scheduler log name", {
  d <- withr::local_tempdir()
  dir.create(RiboPilot:::run_dir(d, "logs"), recursive = TRUE)
  expect_true(is.na(RiboPilot:::run_job_id(d, "pipeline_2026-09-08_12-51-00")))
  file.create(file.path(RiboPilot:::run_dir(d, "logs"), "pipeline_2026-09-08_12-51-00.77.log"))
  expect_equal(RiboPilot:::run_job_id(d, "pipeline_2026-09-08_12-51-00"), "77")
  file.create(file.path(RiboPilot:::run_dir(d, "logs"), "reference_2026-09-08_12-51-00_78.log"))
  expect_equal(RiboPilot:::run_job_id(d, "reference_2026-09-08_12-51-00"), "78")
  file.create(file.path(RiboPilot:::run_dir(d, "logs"), "recheck_2026-09-08_12-51-00.log"))
  expect_true(is.na(RiboPilot:::run_job_id(d, "recheck_2026-09-08_12-51-00")))
})

test_that("scheduler replies classify a job as alive, gone, or unknown", {
  f <- RiboPilot:::job_alive_from
  expect_true(f("squeue", "123 main  job  R", 0L))
  expect_false(f("squeue", character(0), 0L))
  expect_false(f("squeue", "slurm_load_jobs error: Invalid job id specified", 1L))
  expect_true(f("qstat", c("job_number: 123", "sge_o_home: /home/x"), 0L))
  expect_false(f("qstat", "Following jobs do not exist:", 1L))
  expect_true(is.na(f("qstat", "error: commlib error", 1L)))
  expect_true(f("bjobs", "123 me RUN normal", 0L))
  expect_false(f("bjobs", "Job <123> is not found", 255L))
  expect_false(f("bjobs", "123 me EXIT normal", 0L))
  expect_true(is.na(RiboPilot:::job_alive(NA_character_)))
})

test_that("a job the scheduler no longer has gets a killed report", {
  d <- withr::local_tempdir()
  base <- "pipeline_2026-09-08_12-51-00"
  put_log(d, base, rr_log(end = NULL))
  job <- file.path(RiboPilot:::run_dir(d, "jobs"), paste0(base, ".sh"))
  dir.create(dirname(job), recursive = TRUE)
  writeLines(RiboPilot:::NF_EXIT_TRAILER, job)
  dir.create(RiboPilot:::run_dir(d, "logs"))
  writeLines("job started", file.path(RiboPilot:::run_dir(d, "logs"), paste0(base, ".77.log")))

  testthat::local_mocked_bindings(job_alive = function(id) TRUE, .package = "RiboPilot")
  expect_equal(RiboPilot:::sync_run_reports(d), character(0))
  expect_equal(RiboPilot:::open_job_runs(d), base)

  testthat::local_mocked_bindings(job_alive = function(id) NA, .package = "RiboPilot")
  expect_equal(RiboPilot:::sync_run_reports(d), character(0))

  testthat::local_mocked_bindings(job_alive = function(id) { expect_equal(id, "77"); FALSE },
    .package = "RiboPilot")
  expect_equal(RiboPilot:::sync_run_reports(d), base)
  r <- RiboPilot:::list_run_reports(d)
  expect_equal(r$result, "killed")
  expect_true(is.na(r$exit_status))
  expect_equal(r$launch, "job")
  txt <- readLines(RiboPilot:::run_report_paths(d, base)[["txt"]])
  expect_true(any(grepl("Result:    Killed", txt, fixed = TRUE)))
  expect_equal(RiboPilot:::open_job_runs(d), character(0))
  expect_match(RiboPilot:::run_result_html("killed"), "mp-pill-danger")
})
