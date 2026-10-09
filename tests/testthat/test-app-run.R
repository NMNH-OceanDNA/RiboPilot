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
    nf_recheck_command = function(path) c("run", "main_recheck.nf"),
    nf_spawn = function(cmd, args, wd, env) {
      seen$args <- args
      seen$wd <- wd
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

test_that("a finished run gets a report, logs to .runs/nextflow, and frees the project", {
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
  expect_equal(dirname(seen$args[which(seen$args == "-log") + 1]), file.path(d, ".runs", "nextflow"))
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

test_that("Close reads Stop and Close only while a run is going", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  labels <- character(0)
  testthat::local_mocked_bindings(
    updateActionButton = function(session, inputId, label = NULL, ...) {
      if (inputId == "nf_close") labels <<- c(labels, label)
    }, .package = "shiny")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$setInputs(asmb_run_btn = 1)
    session$setInputs(nf_start = 1)
    session$setInputs(nf_stop = 1)
    session$elapse(300)
  })
  expect_equal(labels, c("Stop and Close", "Close", "Close"))
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
  expect_length(list.files(file.path(d, ".runs", "reports"), "_report\\.txt$"), 1)
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
    "found no match, or it could")
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

test_that("opening a project on an image for another version warns and carries on", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  writeLines(c("process {", "    container = 'macguigand/ribopilot:0.0.0.1'", "}"),
    file.path(dir, ".config"))
  withr::local_dir(dir)
  seen <- character(0)
  local_mocked_bindings(mp_alert = function(title, ...) seen <<- c(seen, title))
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    expect_true(nrow(samples_data()) > 0)
  })
  expect_equal(seen, "Container version does not match RiboPilot")
})


# ---- Update on Annotation Review: the recheck window ----

# S1 with one kept, annotated contig whose review copy is queued for a recheck
mk_recheck_project <- function(env = parent.frame()) {
  d <- mk_run_project(env)
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, n_features = 3")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
    VALUES ('S1', 1, 1, 4, 'ACGT', 0)")
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length,
    trim_start, trim_end, edited, check_switch, rev, time_stamp)
    VALUES ('S1', 1, 1, 'ACGT', 4, 1, 4, 1, 1, 1, 100)")
  DBI::dbDisconnect(con)
  d
}

# what query_recheck does when the run picks the contig up
start_copy <- function(d) {
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 4, time_stamp = ?",
    params = list(as.integer(Sys.time())))
  DBI::dbDisconnect(con)
}

copy_state <- function(d) {
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con, "SELECT check_switch, check_note FROM review_units")
}

test_that("contigs waiting for the Assemble step are unlocked, kept, uncopied, and not done", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "DELETE FROM annotate")
  u <- data.frame(scaffold = 1:7, annotate_switch = c(1L, 2L, 2L, 4L, 1L, 1L, 1L),
    annotate_lock = c(0L, 0L, 0L, NA, 1L, 0L, 0L), n_features = c(NA, 3L, NA, 3L, NA, NA, NA),
    ignore = c(0L, 0L, 0L, 0L, 0L, 1L, 0L))
  for (i in seq_len(nrow(u))) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
      VALUES ('S1', 1, ?, 4, 'ACGT', ?)", params = list(u$scaffold[i], u$ignore[i]))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
      n_features, annotate_opts, curate_opts) VALUES ('S1', 1, ?, ?, ?, ?, 'default', 'default')",
      params = list(u$scaffold[i], u$annotate_switch[i], u$annotate_lock[i], u$n_features[i]))
  }
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, check_switch, rev, time_stamp)
    VALUES ('S1', 1, 7, 1, 1, 100)")
  DBI::dbDisconnect(con)
  # 1 queued, 3 never annotated, 4 running with a NULL lock; not 2 (done), 5 (locked),
  # 6 (dropped), or 7 (has a copy)
  expect_identical(RiboPilot:::nf_recheck_waiting(d), 3L)
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  expect_identical(RiboPilot:::nf_recheck_waiting(d), 0L)
})

test_that("with nothing queued, the recheck says how many contigs wait for the Assemble step", {
  skip_if_not_installed("withr")
  d <- mk_run_project()
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1, assemble_switch = 2")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, length, sequence, ignore)
    VALUES ('S1', 1, 1, 4, 'ACGT', 0)")
  DBI::dbDisconnect(con)
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
  })
  expect_equal(seen$alert$title, "Nothing to update")
  expect_match(seen$alert$msg, "Waiting for the Assemble step's Update: 1 contig.", fixed = TRUE)
  expect_null(seen$p)
})

test_that("the recheck will not open while another run holds the project", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  d <- mk_recheck_project()
  withr::local_dir(d)
  m <- RiboPilot:::run_lock_take(d, "pipeline")
  other <- processx::process$new("sleep", "30")
  Sys.unsetenv(m)
  on.exit(other$kill(), add = TRUE)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
  })
  expect_match(seen$alert$msg, "already going")
  expect_null(seen$p)
})

test_that("Update on Annotation Review opens the recheck window, not the pipeline's", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  testthat::local_mocked_bindings(run_modal = function(n, ...) {
    seen$modal <- list(n = n, ...)
    shiny::modalDialog()
  }, .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    expect_equal(output$rc_nf_cmd_text, "nextflow run main_recheck.nf")
  })
  expect_equal(seen$modal$n, 1L)
  expect_equal(seen$modal$pfx, "rc_")
  expect_false(seen$modal$resume)
  expect_equal(seen$modal$title, "RiboPilot recheck - 1 contig(s) to update")
})

test_that("a recheck runs from .recheck, reports contigs it left running, then queues them again", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_nf_start = 1)
    expect_true(file.exists(RiboPilot:::run_lock_path(d)))
    start_copy(d)
    seen$p$finish()
    session$elapse(300)
  })
  expect_equal(seen$wd, file.path(d, ".recheck"))
  expect_true(dir.exists(file.path(d, ".recheck")))
  expect_true("main_recheck.nf" %in% seen$args)
  lg <- seen$args[which(seen$args == "-log") + 1]
  expect_equal(dirname(lg), file.path(d, ".runs", "nextflow"))
  expect_match(basename(lg), "^recheck_")
  # the report is written before the release, so the contig counts as a failure
  expect_equal(seen$alert$title, "Run finished with failures")
  expect_match(seen$alert$msg, "S1.1.1 | recheck did not finish", fixed = TRUE)
  expect_length(list.files(file.path(d, ".runs", "reports"), "^recheck_.*_report\\.txt$"), 1)
  expect_equal(copy_state(d)$check_switch, 1L)
  expect_equal(copy_state(d)$check_note, "recheck did not finish")
  expect_false(file.exists(RiboPilot:::run_lock_path(d)))
})

test_that("a stopped recheck points back to Annotation Review and queues its contigs again", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_nf_start = 1)
    start_copy(d)
    session$setInputs(rc_nf_stop = 1)
    session$elapse(300)
  })
  expect_true(isTRUE(seen$interrupted))
  expect_equal(seen$alert$title, "Run stopped")
  expect_match(seen$alert$msg, "Press Update on Annotation Review.", fixed = TRUE)
  expect_equal(copy_state(d)$check_switch, 1L)
})

test_that("a headless recheck saves its job script under the name of its logs", {
  skip_if_not_installed("withr")
  d <- mk_recheck_project()
  withr::local_dir(d)
  withr::local_options(list(RiboPilot.headless = TRUE))
  seen <- new.env(); mock_run(seen)
  testthat::local_mocked_bindings(detect_scheduler = function() NULL, .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_submit_script_text = "   ")
    session$setInputs(rc_save_script = 1)
  })
  f <- Sys.glob(file.path(d, ".runs", "jobs", "recheck_*.sh"))
  expect_length(f, 1)
  expect_length(Sys.glob(file.path(d, ".runs", "jobs", "pipeline_*.sh")), 0)
  s <- readLines(f)
  base <- sub("\\.sh$", "", basename(f))
  expect_true(any(grepl(paste0(base, ".nextflow.log"), s, fixed = TRUE)))
  expect_true(paste("cd", shQuote(file.path(d, ".recheck"))) %in% s)
  expect_true(any(grepl("'main_recheck.nf'", s, fixed = TRUE)))
  expect_match(seen$alert$msg, "Saved job script to", fixed = TRUE)
})

test_that("Submit as job hands the recheck script to the scheduler once", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  d <- mk_recheck_project()
  withr::local_dir(d)
  withr::local_options(list(RiboPilot.headless = TRUE))
  bin <- withr::local_tempdir()
  calls <- file.path(bin, "calls.txt")
  writeLines(c("#!/bin/sh", sprintf('echo "$@" >> %s', shQuote(calls)),
    'echo "Your job 12345 has been submitted"'), file.path(bin, "qsub"))
  Sys.chmod(file.path(bin, "qsub"), "0755")
  withr::local_path(bin, action = "prefix")
  skip_if(!identical(RiboPilot:::detect_scheduler()$cmd, "qsub"), "another scheduler on PATH")
  seen <- new.env(); mock_run(seen)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(annot_run_btn = 1)
    session$setInputs(rc_submit_job = 1)
    while (!later::loop_empty()) later::run_now()
    expect_match(seen$alert$msg, "Submitted via qsub", fixed = TRUE)
    session$setInputs(rc_submit_job = 2)
    while (!later::loop_empty()) later::run_now()
    expect_match(seen$alert$msg, "already submitted", fixed = TRUE)
  })
  expect_length(readLines(calls), 1L)
  expect_match(readLines(calls), "/\\.runs/jobs/recheck_[0-9_-]+\\.sh")
})
