# Small pure helpers: alignment coverage, N-run detection, DB journal choice,
# constants coherence.

test_that("alignment_coverage is aligned-sample-bases / element length", {
  expect_equal(RiboPilot:::alignment_coverage("ACGT", 4), 1)
  expect_equal(RiboPilot:::alignment_coverage("AC--GT", 6), 4 / 6)   # gaps not counted
  expect_lte(RiboPilot:::alignment_coverage("ACGTACGT", 4), 1)        # capped at 1
  expect_true(is.na(RiboPilot:::alignment_coverage("ACGT", 0)))
})

test_that("element constants are internally consistent", {
  expect_true(all(RiboPilot:::RRNA_GENES %in% RiboPilot:::RRNA_ELEMENT_ORDER))
  expect_true(all(RiboPilot:::RRNA_SPACERS %in% RiboPilot:::RRNA_ELEMENT_ORDER))
  expect_equal(RiboPilot:::rrna_element_type("18S"), "rRNA")
  expect_equal(RiboPilot:::rrna_element_type("ITS1"), "misc_RNA")
  expect_equal(RiboPilot:::rrna_element_product("ITS2"), "internal transcribed spacer 2")
})

test_that("a connection failure is unknown even when the message contains 403", {
  # port 403 refuses the connection; the old substring test read the port number
  # in the message as an HTTP status and reported the database as missing, which
  # check_reference_db() escalates to an error.
  expect_equal(RiboPilot:::url_status("http://127.0.0.1:403/nope"), "unknown")
  expect_warning(RiboPilot:::check_reference_db("http://127.0.0.1:403/nope.fa", "seeds"),
    "Could not verify")
})

test_that("config_container_engine reads the enabled flag, not the block", {
  f <- RiboPilot:::config_container_engine
  expect_equal(f(c("docker { enabled = false }", "singularity { enabled = true }")), "singularity")
  expect_equal(f(c("singularity { enabled = false }", "docker { enabled = true }")), "docker")
  expect_equal(f(c("apptainer {", "  enabled = true", "}")), "singularity")
  expect_true(is.na(f(c("docker { enabled = false }", "singularity { enabled = false }"))))
  expect_true(is.na(f("process { executor = 'slurm' }")))
})

test_that("app connections carry the pipeline's busy timeout", {
  # The app reads and writes the same file as running Nextflow tasks; a bare
  # dbConnect has a busy timeout of 0 and errors on the first contended write.
  dir <- withr::local_tempdir()
  con <- RiboPilot::open_project_db(file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA busy_timeout")[[1]], 10000)

  withr::local_options(golem.app.prod = FALSE)
  local_mocked_bindings(app_project = function() dir, .package = "RiboPilot")
  acon <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(acon), add = TRUE)
  expect_equal(DBI::dbGetQuery(acon, "PRAGMA busy_timeout")[[1]], 10000)
})

test_that("filesystem type is probed once per directory", {
  old <- Sys.getenv("RIBOPILOT_SQLITE_JOURNAL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("RIBOPILOT_SQLITE_JOURNAL") else
    Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = old))
  Sys.unsetenv("RIBOPILOT_SQLITE_JOURNAL")
  rm(list = ls(RiboPilot:::.fs_cache), envir = RiboPilot:::.fs_cache)
  calls <- 0L
  local_mocked_bindings(
    filesystem_type = function(path) { calls <<- calls + 1L; "network" },
    .package = "RiboPilot")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/rp_fs/a.sqlite"), "delete")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/rp_fs/b.sqlite"), "delete")
  expect_equal(calls, 1L)
})

test_that("sqlite journal mode is overridable", {
  old <- Sys.getenv("RIBOPILOT_SQLITE_JOURNAL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("RIBOPILOT_SQLITE_JOURNAL") else
    Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = old))
  Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = "delete")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/x.sqlite"), "delete")
  Sys.setenv(RIBOPILOT_SQLITE_JOURNAL = "wal")
  expect_equal(RiboPilot:::sqlite_journal_mode("/tmp/x.sqlite"), "wal")
})
