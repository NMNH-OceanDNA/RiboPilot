# The Update window: what it refuses, how it stops a run, and what it says at
# the end. Nextflow is replaced by a stand-in handle; nothing is spawned.

mk_run_project <- function(env = parent.frame()) {
  d <- withr::local_tempdir(.local_envir = env)
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  suppressWarnings(suppressMessages(RiboPilot::new_db(db_path = file.path(d, ".sqlite"),
    mapping_fn = map, check_reads = FALSE, seeds_db = fa, labels_db = fa)))
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1")
  DBI::dbDisconnect(con)
  d
}

# a processx stand-in: alive until finished or interrupted
fake_nf <- function(seen) {
  st <- new.env(); st$alive <- TRUE
  list(get_pid = function() 1L, is_alive = function() st$alive,
    read_output_lines = function() "executor >  local (1)",
    read_error_lines = function() character(0),
    interrupt = function() { seen$interrupted <- TRUE; st$alive <- FALSE; TRUE },
    kill = function() { seen$killed <- TRUE; st$alive <- FALSE },
    wait = function(...) NULL,
    get_exit_status = function() if (isTRUE(seen$interrupted)) 1L else 0L,
    finish = function() st$alive <- FALSE)
}

mock_run <- function(seen, env = parent.frame()) {
  testthat::local_mocked_bindings(
    nf_pin_version = function(...) "25.10.6",
    nf_command = function(...) c("run", "main.nf"),
    nf_spawn = function(cmd, args, wd, env) {
      seen$args <- args
      seen$p <- if (is.null(seen$spawn_error)) fake_nf(seen) else stop(seen$spawn_error, call. = FALSE)
      seen$p
    },
    rp_alert = function(msg, type = "message", title = NULL, session = NULL) {
      seen$alert <- list(msg = paste(as.character(msg), collapse = " "), type = type, title = title)
    },
    .package = "RiboPilot", .env = env)
}

test_that("Update refuses to open the run window when nothing is queued", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 2")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2")
  DBI::dbDisconnect(con)
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
  })
  expect_equal(seen$alert$title, "Nothing to update")
})

test_that("a missing Nextflow is reported, not a crash, and leaves no lock", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); seen$spawn_error <- "Nextflow not found on PATH."
  mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
  })
  expect_equal(seen$alert$type, "error")
  expect_match(seen$alert$msg, "Nextflow not found on PATH")
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("a finished run gets a report, logs to nextflow_logs, and frees the project", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
    expect_true(file.exists(RiboPilot:::run_lock_path(d)))
    seen$p$finish()
    session$elapse(300)
  })
  expect_equal(seen$alert$title, "Run finished")
  expect_equal(dirname(seen$args[which(seen$args == "-log") + 1]), file.path(d, "nextflow_logs"))
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("Stop interrupts Nextflow instead of killing it", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
    session$setInputs(nf_stop = 1)
    session$elapse(300)
  })
  expect_true(isTRUE(seen$interrupted))
  expect_null(seen$killed)
  expect_equal(seen$alert$title, "Run stopped")
})

test_that("closing the browser stops the run gracefully", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
  })
  expect_true(isTRUE(seen$interrupted))
  expect_null(seen$killed)
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
  expect_length(list.files(file.path(d, "nextflow_logs"), "_report\\.txt$"), 1)
})

test_that("the project stays locked while a task container of the run still runs", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  box <- new.env(); box$ids <- "c1"
  testthat::local_mocked_bindings(
    run_containers = function(project, all = FALSE) box$ids,
    run_clear_containers = function(project) seen$cleared <- TRUE,
    .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
    seen$p$finish()
    session$elapse(300)
    expect_true(file.exists(RiboPilot:::run_lock_path(d)))
    box$ids <- character(0)
    session$elapse(300)
  })
  expect_true(isTRUE(seen$cleared))
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("a second window cannot start a run while one holds the project lock", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  d <- mk_run_project()
  withr::local_dir(d)
  m <- RiboPilot:::run_lock_take(d, "pipeline")
  other <- processx::process$new("sleep", "30")
  Sys.unsetenv(m)
  on.exit(other$kill(), add = TRUE)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
  })
  expect_match(seen$alert$msg, "already going")
  expect_null(seen$p)
})

test_that("a failed sample's details say why and what to try next", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 3,
    assemble_notes = '3 or more paths exceed max_paths (2)'")
  DBI::dbDisconnect(con)
  withr::local_dir(d)
  html <- as.character(RiboPilot:::failure_banner("S1"))
  expect_match(html, "This sample failed", fixed = TRUE)
  expect_match(html, "Max paths", fixed = TRUE)
  expect_match(html, "3 or more paths exceed max_paths (2)", fixed = TRUE)
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = 2")
  DBI::dbDisconnect(con)
  expect_null(RiboPilot:::failure_banner("S1"))
})

test_that("a scaffold of a failed sample says BLAST did not run, not that NCBI failed", {
  expect_match(RiboPilot:::no_screen_note(failed = TRUE, blast_on = TRUE, running = FALSE),
    "did not run because the assembly failed")
  expect_match(RiboPilot:::no_screen_note(failed = FALSE, blast_on = TRUE, running = FALSE),
    "did not return")
})

test_that("an LSF job script goes to bsub on stdin, others as an argument", {
  skip_on_os("windows")
  bin <- withr::local_tempdir()
  writeLines(c("#!/bin/sh", 'echo "args=$# first=$(head -n 1)"'), file.path(bin, "bsub"))
  writeLines(c("#!/bin/sh", 'echo "args=$#"'), file.path(bin, "sbatch"))
  Sys.chmod(file.path(bin, c("bsub", "sbatch")), "0755")
  withr::local_path(bin, action = "prefix")
  f <- file.path(bin, "job.sh")
  writeLines(c("#BSUB -q long", "nextflow run"), f)
  expect_equal(RiboPilot:::submit_job("bsub", f), "args=0 first=#BSUB -q long")
  expect_equal(RiboPilot:::submit_job("sbatch", f), "args=1")
})
