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

test_that("a submitted run's tasks are found in nextflow_logs/", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "nextflow_logs"))
  wd <- file.path(dir, "work", "ab", "cdef12")
  writeLines(paste0("Sep-24 10:00:00.000 [Task monitor] DEBUG n.processor.TaskPollingMonitor - ",
    "Task completed > TaskHandler[jobId: 4242; id: 7; name: assemble (S1); status: COMPLETED; ",
    "exit: 140; error: -; workDir: ", wd, " started: 1727172000000; exited: 2026-09-24T10:00:00Z; ]"),
    file.path(dir, "nextflow_logs", "ribopilot_2026-09-24_10-00-00.nextflow.log"))
  d <- RiboPilot:::find_workdirs(dir, "S1")
  expect_equal(d$process, "assemble")
  expect_equal(d$status, "failed")
  expect_equal(d$workdir, wd)
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

# A fake desktop opener on PATH that records what it was asked to open, so a
# test can tell "opened the file manager" from "did not".
mk_opener_probe <- function(env) {
  bin <- withr::local_tempdir(.local_envir = env)
  marker <- file.path(bin, "opened.txt")
  opener <- if (Sys.info()[["sysname"]] == "Darwin") "open" else "xdg-open"
  writeLines(c("#!/bin/sh", paste0("echo \"$1\" >> ", shQuote(marker))),
    file.path(bin, opener))
  Sys.chmod(file.path(bin, opener), "0755")
  withr::local_path(bin, action = "prefix", .local_envir = env)
  marker
}

test_that("a remote session never reaches for a desktop file manager", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  marker <- mk_opener_probe(rlang::current_env())
  # RStudio Server on a login node: xdg-open there reaches no desktop, and
  # RStudio intercepts it into a URL that 404s
  withr::local_options(list(RiboPilot.headless = TRUE))
  testthat::local_mocked_bindings(rstudio_show_dir = function(path) TRUE, .package = "RiboPilot")

  expect_true(RiboPilot:::open_path("/tmp/ribopilot-remote"))
  Sys.sleep(0.2)
  expect_false(file.exists(marker))
})

test_that("a remote session with no Files pane reports that it opened nothing", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  marker <- mk_opener_probe(rlang::current_env())
  withr::local_options(list(RiboPilot.headless = TRUE))
  testthat::local_mocked_bindings(rstudio_show_dir = function(path) FALSE, .package = "RiboPilot")

  # FALSE, not a silent no-op: the caller has to tell the user where the folder
  # is, since a button that appears to do nothing is the bug being fixed
  expect_false(RiboPilot:::open_path("/tmp/ribopilot-remote-nopane"))
  Sys.sleep(0.2)
  expect_false(file.exists(marker))
})

test_that("a desktop session still opens the OS file manager", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  marker <- mk_opener_probe(rlang::current_env())
  withr::local_options(list(RiboPilot.headless = FALSE))

  expect_true(RiboPilot:::open_path("/tmp/ribopilot-desktop"))
  for (i in 1:40) { if (file.exists(marker)) break; Sys.sleep(0.05) }
  expect_true(file.exists(marker))
  expect_equal(readLines(marker), "/tmp/ribopilot-desktop")
})

test_that("a browser on another computer never reaches for a desktop file manager", {
  skip_if_not_installed("withr")
  skip_on_os("windows")
  marker <- mk_opener_probe(rlang::current_env())
  withr::local_options(list(RiboPilot.headless = FALSE, RiboPilot.local_browser = FALSE))
  testthat::local_mocked_bindings(rstudio_show_dir = function(path) FALSE, .package = "RiboPilot")

  expect_false(RiboPilot:::open_path("/tmp/ribopilot-remote-browser"))
  Sys.sleep(0.2)
  expect_false(file.exists(marker))
})

test_that("unit_evidence surfaces the screen reason", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, sci_name TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT, ignore INT,
    triage TEXT, triage_reason TEXT, triage_source TEXT, source_contig TEXT,
    blast_species TEXT, blast_accession TEXT, blast_pident REAL, blast_qcovs REAL, blast_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT, status TEXT, flags TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    median_depth REAL, depths TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL, qcovs REAL,
    bitscore REAL, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE path_scores(ID TEXT, path INT, score REAL)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, itsx_warnings TEXT)")
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, RiboPilot:::CHIMERA_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_SEARCHED_DDL)
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio',NULL)")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID,path,scaffold,length) VALUES('S1',1,1,500)")
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_reason,screen_version)
     VALUES('S1',1,1,'no_reference',
     'the bundled rDNA database returned no usable hit for this scaffold; a search at NCBI (BLAST Options, Search) may identify it',1)")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  d <- RiboPilot:::unit_evidence("S1")
  expect_true("screen_reason" %in% names(d))
  expect_match(d$screen_reason[1], "search at NCBI")
})

test_that("the Files pane helper declines outside RStudio rather than erroring", {
  # this test process is not an RStudio session, which is the same position a
  # tunnelled plain-R session on a login node is in
  expect_false(RiboPilot:::rstudio_show_dir(tempdir()))
})

test_that("unit_evidence lists each flagged feature with its warnings", {
  con <- chimera_test_con()
  db <- DBI::dbGetInfo(con)$dbname
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, status, flags) VALUES
    ('S1',1,1,'18S',1,2,'flag','length,truncated'), ('S1',1,1,'ITS1',3,3,'ok',NULL), ('S1',1,1,'28S',4,4,'flag','overlap')")
  DBI::dbDisconnect(con)
  withr::local_options(list(ribopilot.db = db))
  withr::local_dir(dirname(db))
  d <- RiboPilot:::unit_evidence("S1")
  expect_equal(d$elements, "18S,ITS1,28S")
  expect_equal(d$element_flags, "18S: length, truncated; 28S: overlap")
  expect_match(as.character(RiboPilot:::rt_element_strip()), "element_flags", fixed = TRUE)
})
