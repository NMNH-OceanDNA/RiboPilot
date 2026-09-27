# What the app says when a run ends, and the project's one-run-at-a-time lock.

# Lines copied from real Nextflow 25.10.6 logs (a local docker run of a toy
# pipeline): the formats the report reads.
nf_log_lines <- function(work) c(
  "Sep-24 11:47:16.706 [main] DEBUG nextflow.cli.CmdRun - N E X T F L O W  ~  version 25.10.6",
  "Sep-24 11:47:17.116 [main] DEBUG nextflow.Session - Run name: elegant_wegener",
  sprintf("Sep-24 11:47:18.128 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 3; name: assemble (S3); status: COMPLETED; exit: 127; error: -; workDir: %s/b7/fa6579]", work),
  "Sep-24 11:47:18.138 [TaskFinalizer-1] INFO  nextflow.processor.TaskProcessor - [b7/fa6579] NOTE: Process `assemble (S3)` terminated with an error exit status (127) -- Error is ignored",
  sprintf("Sep-24 11:47:18.191 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 1; name: assemble (S1); status: COMPLETED; exit: 137; error: -; workDir: %s/67/154f91]", work),
  "Sep-24 11:47:18.196 [TaskFinalizer-3] INFO  nextflow.processor.TaskProcessor - [67/154f91] NOTE: Process `assemble (S1)` terminated with an error exit status (137) -- Execution is retried (1)",
  sprintf("Sep-24 11:47:18.510 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 5; name: assemble (S1); status: COMPLETED; exit: 137; error: -; workDir: %s/e8/41c41e]", work),
  "Sep-24 11:47:18.512 [TaskFinalizer-5] INFO  nextflow.processor.TaskProcessor - [e8/41c41e] NOTE: Process `assemble (S1)` terminated with an error exit status (137) -- Error is ignored",
  sprintf("Sep-24 11:47:18.600 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 6; name: assemble (S2); status: COMPLETED; exit: 137; error: -; workDir: %s/aa/000001]", work),
  "Sep-24 11:47:18.610 [TaskFinalizer-6] INFO  nextflow.processor.TaskProcessor - [aa/000001] NOTE: Process `assemble (S2)` terminated with an error exit status (137) -- Execution is retried (1)",
  sprintf("Sep-24 11:47:18.700 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 7; name: assemble (S2); status: COMPLETED; exit: 0; error: -; workDir: %s/aa/000002]", work),
  sprintf("Sep-24 12:09:57.372 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 8; name: score (S4); status: COMPLETED; exit: 0; error: -; workDir: %s/82/c1366e]", work),
  "Sep-24 12:09:57.385 [TaskFinalizer-1] INFO  nextflow.processor.TaskProcessor - [82/c1366e] NOTE: Missing output file(s) `S4.txt` expected by process `score (S4)` -- Error is ignored",
  "Sep-24 11:47:18.817 [main] DEBUG n.trace.WorkflowStatsObserver - Workflow completed > WorkflowStats[succeededCount=1; failedCount=5; ignoredCount=3; cachedCount=0; pendingCount=0; submittedCount=0; runningCount=0; retriesCount=2; abortedCount=0; succeedDuration=498ms]")

test_that("the report names each task whose last attempt failed, with its exit code and cause", {
  d <- withr::local_tempdir()
  work <- file.path(d, "work")
  dir.create(file.path(work, "e8", "41c41e"), recursive = TRUE)
  writeLines(c("SPAdes started", "", "== Error ==  finished abnormally, OS return value: -9"),
    file.path(work, "e8", "41c41e", ".command.err"))
  log <- file.path(d, "run.nextflow.log")
  writeLines(nf_log_lines(work), log)

  f <- RiboPilot:::nf_log_facts(log)
  expect_equal(f$run, "elegant_wegener")
  expect_equal(f$version, "25.10.6")
  expect_equal(f$stats[["ignored"]], 3L)
  # S2 failed once and then succeeded on retry, so it is not a failure
  expect_setequal(f$tasks$item, c("S1", "S3", "S4"))
  s1 <- f$tasks[f$tasks$item == "S1", ]
  expect_equal(s1$process, "assemble")
  expect_equal(s1$exit, "137")
  expect_equal(s1$tries, 2L)
  expect_match(s1$cause, "out of memory")
  expect_match(s1$stderr, "OS return value: -9")
  # a task that exited 0 but lost its output says so, not "exit 0"
  expect_match(f$tasks$cause[f$tasks$item == "S4"], "Missing output")
})

test_that("a version mismatch is the report's lead", {
  d <- withr::local_tempdir()
  work <- file.path(d, "work")
  dir.create(file.path(work, "e8", "41c41e"), recursive = TRUE)
  writeLines(c("SPAdes started", "",
    "Version mismatch: this run was started by RiboPilot 9.9.9, but ..."),
    file.path(work, "e8", "41c41e", ".command.err"))
  log <- file.path(d, "run.nextflow.log")
  writeLines(nf_log_lines(work), log)

  rep <- RiboPilot:::nf_run_report(d, log, 0L, "pipeline", 0)
  expect_match(rep$lead, "^Version mismatch:")
})

test_that("the run error is the abort cause, or else the first ERROR line", {
  d <- withr::local_tempdir()
  a <- file.path(d, "a.log")
  writeLines(c("Sep-24 11:47:37.916 [TaskFinalizer-5] DEBUG nextflow.Session - Session aborted -- Cause: Process `preprocess (S9)` terminated with an error exit status (3)",
    "Sep-24 11:47:37.920 [main] ERROR nextflow.cli.Launcher - something later"), a)
  expect_match(RiboPilot:::nf_log_facts(a)$error, "^Process `preprocess \\(S9\\)`")
  b <- file.path(d, "b.log")
  writeLines("Sep-24 11:47:59.489 [main] ERROR nextflow.cli.Launcher - Unknown config attribute `x`", b)
  expect_equal(RiboPilot:::nf_log_facts(b)$error, "Unknown config attribute `x`")
  expect_true(is.na(RiboPilot:::nf_log_facts(file.path(d, "none.log"))$error))
})

test_that("a task that failed before it started, or wrote bad bytes, still gets a report", {
  d <- withr::local_tempdir()
  writeBin(as.raw(c(0x65, 0x72, 0x72, 0xff, 0x0a)), file.path(d, ".command.err"))
  log <- file.path(d, "run.nextflow.log")
  writeLines(c(nf_log_lines(d)[1:2],
    "Sep-24 11:47:37.900 [TaskFinalizer-1] ERROR nextflow.processor.TaskProcessor - Error executing process > 'preprocess (S9)'",
    sprintf("Sep-24 11:47:38.000 [Task monitor] DEBUG n.processor.TaskPollingMonitor - Task completed > TaskHandler[id: 2; name: assemble (S1); status: COMPLETED; exit: 1; error: -; workDir: %s]", d),
    "Sep-24 11:47:38.100 [TaskFinalizer-2] INFO  nextflow.processor.TaskProcessor - [07/65d217] NOTE: Process `assemble (S1)` terminated with an error exit status (1) -- Error is ignored"), log)
  f <- RiboPilot:::nf_log_facts(log)
  expect_equal(f$tasks$exit[f$tasks$item == "S9"], "-")
  expect_match(f$tasks$stderr[f$tasks$item == "S1"], "^err")
  expect_equal(RiboPilot:::nf_run_report(d, log, 1L, "pipeline", 0)$type, "error")
})

test_that("exit codes map to a likely cause", {
  f <- RiboPilot:::nf_exit_cause
  expect_match(f(137L), "out of memory")
  expect_match(f(140L), "scheduler")
  expect_match(f(143L), "terminated")
  expect_match(f(127L), "not found")
  expect_match(f(247L), "out of memory")
  expect_match(f(NA), "no exit status")
  expect_equal(f(1L), "the tool stopped with an error")
})

test_that("only a clean run is reported green, and every failure carries a copyable report", {
  d <- withr::local_tempdir()
  log <- file.path(d, "run.nextflow.log")
  writeLines(nf_log_lines(file.path(d, "work")), log)
  clean <- file.path(d, "clean.log")
  writeLines(nf_log_lines(file.path(d, "work"))[1:2], clean)

  ok <- RiboPilot:::nf_run_report(d, clean, 0L, "pipeline", 0)
  expect_equal(ok$type, "message")
  expect_equal(ok$text, "")

  bad <- RiboPilot:::nf_run_report(d, log, 0L, "pipeline", 0)
  expect_equal(bad$type, "warning")
  expect_match(bad$text, "RiboPilot:", fixed = TRUE)
  expect_match(bad$text, "Nextflow:  25.10.6", fixed = TRUE)
  expect_match(bad$text, paste("Log:      ", log), fixed = TRUE)
  expect_match(bad$text, "S1 | assemble | 137 | killed", fixed = TRUE)

  expect_equal(RiboPilot:::nf_run_report(d, log, 1L, "pipeline", 0)$type, "error")
  stopped <- RiboPilot:::nf_run_report(d, clean, 1L, "pipeline", 0, how = "stopped")
  expect_equal(stopped$title, "Run stopped")
  expect_match(stopped$lead, "Resume")
})

test_that("the report lists the samples this run left Failed, with their notes", {
  d <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, assemble_notes TEXT,
    find_rrna_notes TEXT, time_stamp INT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES
    ('NEW', 3, 'below min_depth: 812 reads after trimming, both mates counted (minimum 2000)', NULL, 200),
    ('OLD', 3, 'failed assembly', NULL, 50), ('OK', 2, NULL, NULL, 200)")
  DBI::dbDisconnect(con)
  bad <- RiboPilot:::run_failed_samples(d, 100)
  expect_equal(bad$ID, "NEW")
  r <- RiboPilot:::nf_run_report(d, file.path(d, "none.log"), 0L, "pipeline", 100)
  expect_equal(r$type, "warning")
  expect_match(r$lead, "Run finished, but 1 sample(s) failed.", fixed = TRUE)
  expect_match(r$text, "NEW | below min_depth: 812", fixed = TRUE)
  expect_match(r$text, "Too few reads", fixed = TRUE)
})

test_that("each pipeline note gets a plain cause and a next step", {
  h <- RiboPilot:::failure_hint
  expect_match(h("below min_depth: 812 reads after trimming")$why, "Too few reads")
  expect_match(h("3 or more paths exceed max_paths (2)")$next_step, "Max paths")
  expect_match(h("failed assembly")$next_step, "type a new name in Parameter set name", fixed = TRUE)
  expect_match(h("5 scaffolds in a path exceeds max_scaffolds (4)")$next_step, "Max scaffolds")
  expect_match(h("2 scaffold(s) under min_assembly_length (500 bp)")$next_step, "SPAdes, baited")
  expect_match(h("no rRNA genes found in the supplied assembly")$why, "assembly you supplied")
  expect_match(h("rRNA search incomplete: 1 of 3 part(s) of the supplied assembly could not be scanned")$next_step,
    "Ready to run")
  expect_match(h("failed assembly: SPAdes exit 247: finished abnormally, OS return value: -9")$why,
    "memory")
  expect_match(h("failed assembly")$next_step, "SPAdes, baited")
  expect_match(h(NA_character_)$why, "no rDNA scaffolds")
  expect_match(h("failed assembly: no usable seed database to bait reads with")$next_step,
    "Seed database and Label database", fixed = TRUE)
  expect_match(h("failed assembly: no usable label database to pick out rDNA contigs")$why,
    "no seed or label database")
})

test_that("a project lock refuses a second run and clears once the run is gone", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  m <- RiboPilot:::run_lock_take(d, "pipeline")
  p <- processx::process$new("sleep", "30")
  Sys.unsetenv(m)
  held <- RiboPilot:::run_lock_holder(d)
  expect_false(is.null(held))
  expect_match(RiboPilot:::run_lock_message(held, d), "already going")
  p$kill()
  expect_null(RiboPilot:::run_lock_holder(d))
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("a lock taken on another computer is reported, not cleared", {
  d <- withr::local_tempdir()
  write.dcf(data.frame(kind = "pipeline", marker = "RIBOPILOT_RUN_1_1", host = "login2.cluster",
    started = "2026-09-24 10:00"), RiboPilot:::run_lock_path(d))
  held <- RiboPilot:::run_lock_holder(d)
  expect_equal(held$host, "login2.cluster")
  msg <- RiboPilot:::run_lock_message(held, d)
  expect_match(msg, "login2.cluster", fixed = TRUE)
  expect_match(msg, RiboPilot:::run_lock_path(d), fixed = TRUE)
  expect_true(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("stopping a run interrupts it first and kills only a run that ignores that", {
  skip_on_os("windows")
  skip_on_cran()
  p <- processx::process$new("sleep", "30")
  expect_true(RiboPilot:::nf_stop_run(p, "RIBOPILOT_RUN_none", wait = 5))
  expect_false(p$is_alive())
  expect_equal(p$get_exit_status(), -2L)   # SIGINT, not SIGKILL
  stubborn <- processx::process$new("bash", c("-c", "trap '' INT; sleep 30"))
  Sys.sleep(0.5)
  expect_false(RiboPilot:::nf_stop_run(stubborn, "RIBOPILOT_RUN_none", wait = 1))
  expect_false(stubborn$is_alive())
})

test_that("a run started with SIGINT ignored is stopped with SIGTERM, not killed", {
  skip_on_os("windows")
  skip_on_cran()
  deaf <- processx::process$new("bash", c("-c", "trap '' INT; exec sleep 30"))
  Sys.sleep(0.5)
  expect_true(RiboPilot:::nf_stop_run(deaf, "RIBOPILOT_RUN_none", wait = 10, term_wait = 1))
  expect_equal(deaf$get_exit_status(), -15L)
})

test_that("a project's task containers are the nxf ones mounting a folder inside it", {
  d <- withr::local_tempdir()
  testthat::local_mocked_bindings(docker_out = function(args) {
    if (args[1] == "ps") return("c1\nc2\nc3\n")
    paste0(c(sprintf("c1\t/data/refs\t%s/work/ab/cd", d), "c2\t/other/project/work/ab",
      sprintf("c3\t%s2/work/ab", d)), collapse = "\n")
  }, .package = "RiboPilot")
  expect_equal(RiboPilot:::run_containers(d), "c1")
})

test_that("runs pin Nextflow and switch off its update check", {
  testthat::local_mocked_bindings(nf_pin_version = function(...) "25.10.6", .package = "RiboPilot")
  env <- RiboPilot:::nf_env()
  expect_equal(env[["NXF_VER"]], "25.10.6")
  expect_equal(env[["NXF_DISABLE_CHECK_LATEST"]], "true")
  expect_equal(unname(env[1]), "current")
})
