# A sample's state must read "in progress" (4) while WF1 still has work to do
# for it, and "done" (2) only once nothing is left. Assembly loading is the
# start of that window, triage the end.

mk_state_db <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT,
                       length INT, length_raw INT, sequence TEXT, ignore INT, edited INT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, paths INT, scaffolds INT, length TEXT,
                       topology TEXT, assemble_switch INT, assemble_notes TEXT,
                       assemble_opts TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble_opts(assemble_opts TEXT, max_paths INT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
                       annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT,
                       PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES('default', 10)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',NULL,NULL,NULL,NULL,1,NULL,'default',NULL)")
  con
}

write_state_fasta <- function(dir, id = "S1") {
  ad <- file.path(dir, id, "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste0(">", id, ".1.1 circular scaffold_1--101+"), strrep("ACGT", 200),
    paste0(">", id, ".1.2 linear scaffold_2--103+"), strrep("TGCA", 150)),
    file.path(ad, paste0(id, "_assembly_1.fasta")))
  invisible(ad)
}

sample_state <- function(con, id = "S1") {
  DBI::dbGetQuery(con, "SELECT assemble_switch FROM assemble WHERE ID = ?",
    params = list(id))$assemble_switch
}

test_that("loading an assembly leaves the sample in progress, not done", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  write_state_fasta(dir)

  RiboPilot:::load_assemblies(con, dir, "S1")

  # the units it just seeded are queued for annotation, so WF1 is not finished
  expect_equal(sample_state(con), 4L)
})

test_that("a failed assembly still reads failed", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::load_assemblies(con, dir, "S1")   # no fasta written
  expect_equal(sample_state(con), 3L)
})

test_that("refresh_sample_state settles a sample once no unit is queued", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  write_state_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")

  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2 WHERE scaffold = 1")
  RiboPilot:::refresh_sample_state(con, "S1")
  expect_equal(sample_state(con), 4L)           # scaffold 2 still queued

  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2 WHERE scaffold = 2")
  RiboPilot:::refresh_sample_state(con, "S1")
  expect_equal(sample_state(con), 2L)
})

test_that("units the user parked do not hold the sample open", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  write_state_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2 WHERE scaffold = 1")

  # held (state 0) is a deliberate skip, not pending work
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 0 WHERE scaffold = 2")
  RiboPilot:::refresh_sample_state(con, "S1")
  expect_equal(sample_state(con), 2L)

  # and an ignored scaffold is invisible everywhere else, so it cannot pend here
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE scaffold = 2")
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE scaffold = 2")
  RiboPilot:::refresh_sample_state(con, "S1")
  expect_equal(sample_state(con), 2L)
})

test_that("dropping the last queued unit settles the sample", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  write_state_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2 WHERE scaffold = 1")
  expect_equal(sample_state(con), 4L)

  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE ID='S1' AND path=1 AND scaffold=2")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 1L)
  expect_equal(sample_state(con), 2L)

  # putting it back reopens the sample: it has to be annotated again
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 0 WHERE ID='S1' AND path=1 AND scaffold=2")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 0L)
  expect_equal(sample_state(con), 4L)
})

test_that("refresh_sample_state does not overwrite a state the pipeline did not set", {
  dir <- withr::local_tempdir()
  con <- mk_state_db(dir); on.exit(DBI::dbDisconnect(con))
  write_state_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2")

  for (st in c(0L, 1L, 3L)) {   # held, queued to re-run, failed
    DBI::dbExecute(con, "UPDATE assemble SET assemble_switch = ?", params = list(st))
    RiboPilot:::refresh_sample_state(con, "S1")
    expect_equal(sample_state(con), st, info = paste("state", st))
  }
})

# the triage path reads a wider slice of the schema than assembly loading does
mk_triage_state_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT,
    ignore INT, triage TEXT, triage_reason TEXT, triage_source TEXT, triage_rule TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, n_features INT,
    PRIMARY KEY (ID, path, scaffold))")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT,
    mean_depth REAL, median_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, species TEXT,
    lineage TEXT, lineage_taxids TEXT, rank_distance INT, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE screen(ID TEXT, path INT, scaffold INT, screen_verdict TEXT,
    screen_class TEXT, screen_confidence REAL, lca_taxon TEXT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, sci_name TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_switch INT, triage_opts TEXT)")
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',4,'default')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','1;2')")
  for (s in 1:2) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,?,5000,0,NULL,NULL,NULL,NULL)",
      params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,?,2,0,'default','default','no',1)",
      params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,?,100,100)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(s))
  }
  con
}

test_that("reaching triage marks the sample done", {
  con <- mk_triage_state_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(sample_state(con), 4L)

  RiboPilot:::auto_triage_sample(con, "S1")

  expect_equal(sample_state(con), 2L)
})

test_that("triage settles the state even when it decides nothing", {
  # every rule off is a deliberate no-op that returns early, but the sample has
  # still been through the whole workflow by the time triage runs
  con <- mk_triage_state_db(); on.exit(DBI::dbDisconnect(con))
  for (k in RiboPilot:::triage_rule_keys()) {
    DBI::dbExecute(con, sprintf("UPDATE triage_opts SET rule_%s = 'off'", k))
  }
  off <- RiboPilot:::triage_all_off(
    RiboPilot:::normalize_triage_actions(RiboPilot:::triage_opts_for(con, "S1")))
  skip_if_not(off, "fixture no longer has every rule off")

  expect_equal(RiboPilot:::auto_triage_sample(con, "S1"), 0L)
  expect_equal(sample_state(con), 2L)
})

test_that("a scaffold annotation found nothing on does not hold the sample open", {
  # annotation ran and legitimately returned no features; triage then reviews the
  # unit and leaves it kept. That must not re-queue it, or the sample can never
  # leave "in progress".
  con <- mk_triage_state_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 2")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0 WHERE scaffold = 2")

  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 0L)

  st <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(st$annotate_switch, 2L)
  expect_equal(sample_state(con), 2L)
})

test_that("a unit that was never annotated still queues", {
  con <- mk_triage_state_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 2")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL WHERE scaffold = 2")

  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 0L)

  st <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(st$annotate_switch, 1L)
  expect_equal(sample_state(con), 4L)
})
