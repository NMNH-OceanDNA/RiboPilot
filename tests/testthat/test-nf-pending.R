# One workflow, one Update button, so one number: samples queued for assembly
# plus units queued for annotation only. It must agree with what the two
# Nextflow queries select, or the modal promises work the run will not do.

mk_pending_project <- function(dir) {
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, assemble_lock INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT,
    annotate_switch INT, annotate_lock INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, ignore INT)")
  DBI::dbDisconnect(con)
  db
}

test_that("nf_pending_count adds assembly work to annotate-only work", {
  dir <- withr::local_tempdir()
  db <- mk_pending_project(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))

  expect_equal(RiboPilot:::nf_pending_count(dir), 0L)

  # a sample queued for assembly
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('A',1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 1L)

  # a second sample, assembled already, with one unit queued for annotation
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('B',2,0)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('B',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('B',1,1,1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)

  # a queued unit on the sample that is about to be assembled is NOT extra work:
  # assembly re-seeds it, and the Nextflow query excludes it for that reason
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('A',1,1,0)")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('A',1,1,1,0)")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)

  # a locked sample is not work
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1 WHERE ID = 'A'")
  expect_equal(RiboPilot:::nf_pending_count(dir), 2L)  # A drops as assembly, returns as a unit
})

test_that("nf_pending_count takes no workflow argument", {
  # There is one workflow; a caller passing "annotate" would silently get a
  # count for something that no longer exists.
  expect_named(formals(RiboPilot:::nf_pending_count), "path")
})

test_that("every nf_command() caller matches its signature", {
  # A real run caught run_assembly() still passing the removed workflow
  # argument, which the suite missed because nothing calls these drivers with a
  # project. Check the calls against the formals instead.
  fns <- list(run_assembly = RiboPilot::run_assembly,
    build_submit_script = RiboPilot:::build_submit_script,
    app_server = RiboPilot:::app_server)
  ok <- names(formals(RiboPilot:::nf_command))
  for (nm in names(fns)) {
    txt <- paste(deparse(fns[[nm]]), collapse = " ")
    calls <- regmatches(txt, gregexpr("nf_command\\([^)]*\\)", txt))[[1]]
    for (cl in calls) {
      args <- strsplit(sub("^nf_command\\(", "", sub("\\)$", "", cl)), ",")[[1]]
      args <- trimws(args)
      named <- grepl("=", args, fixed = TRUE)
      # only the project path may be positional. run_assembly() passed a second
      # positional argument, "assemble", which silently bound to `resume` once
      # the workflow parameter was removed, and the run died at launch.
      expect_lte(sum(!named), 1L, label = paste(nm, cl))
      for (a in sub("\\s*=.*", "", args[named])) {
        expect_true(a %in% ok, info = paste(nm, cl, a))
      }
    }
  }
})
