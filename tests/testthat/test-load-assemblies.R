# load_assemblies() replaces a sample's rows wholesale on every (re)assembly.
# The curator's keep/drop decisions must survive that, or changing assemble_opts
# (which re-queues the sample) silently un-ignores every dropped scaffold.

mk_asm_db <- function(dir) {
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

# GetOrganelle-style output: one path fasta, two scaffolds.
write_asm_fasta <- function(dir, id = "S1") {
  ad <- file.path(dir, id, "assemble")
  dir.create(ad, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    paste0(">", id, ".1.1 circular scaffold_1--101+,102+"), strrep("ACGT", 200),
    paste0(">", id, ".1.2 linear scaffold_2--103+"), strrep("TGCA", 150)),
    file.path(ad, paste0(id, "_assembly_1.fasta")))
  invisible(ad)
}

test_that("a re-assembly preserves per-scaffold ignore flags", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  write_asm_fasta(dir)

  RiboPilot:::load_assemblies(con, dir, "S1")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM assemblies")$n, 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate")$n, 2)

  # WF1 annotates every unit as it is assembled, so give scaffold 2 a feature
  # model before dropping it -- that is the work a destructive drop would lose.
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,2,'18S')")

  # curator drops scaffold 2
  DBI::dbExecute(con, "UPDATE assemblies SET ignore = 1 WHERE ID='S1' AND path=1 AND scaffold=2")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 1L)

  # HELD, not deleted: the unit is already invisible downstream (everything
  # filters on assemblies.ignore), so dropping its curation buys nothing and
  # makes "drop it, look again, put it back" cost a re-annotation.
  st <- DBI::dbGetQuery(con, "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(nrow(st), 1L)
  expect_equal(st$annotate_switch, 0L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotations WHERE scaffold=2")$n, 1)

  # re-run assembly over the same output
  RiboPilot:::load_assemblies(con, dir, "S1")
  ig <- DBI::dbGetQuery(con, "SELECT scaffold, ignore FROM assemblies WHERE ID='S1' ORDER BY scaffold")
  expect_equal(ig$ignore, c(0L, 1L))
  # a re-assembly must not treat a curator-dropped unit as an orphan
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate")$n, 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotations WHERE scaffold=2")$n, 1)
})

test_that("un-ignoring a held unit does not re-queue an already-annotated one", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  write_asm_fasta(dir)
  RiboPilot:::load_assemblies(con, dir, "S1")
  DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,2,'18S')")

  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 1L)   # drop
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 2L, 0L)   # put back
  st <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=2")
  expect_equal(st$annotate_switch, 2L)   # already annotated -> done, no WF2 round trip

  # a unit with no feature model still queues for annotation
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold=1")
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 1L, 1L)
  RiboPilot:::reconcile_ignore_unit(con, "S1", 1L, 1L, 0L)
  st1 <- DBI::dbGetQuery(con,
    "SELECT annotate_switch FROM annotate WHERE ID='S1' AND path=1 AND scaffold=1")
  expect_equal(st1$annotate_switch, 1L)
})

test_that("assemble.scaffolds records total units, not the largest path's count", {
  dir <- withr::local_tempdir()
  con <- mk_asm_db(dir)
  on.exit(DBI::dbDisconnect(con))
  ad <- write_asm_fasta(dir)
  # a second path with one scaffold: 3 units total, max 2 per path
  writeLines(c(">S1.2.1 circular scaffold_1--104+", strrep("ACGT", 180)),
    file.path(ad, "S1_assembly_2.fasta"))

  RiboPilot:::load_assemblies(con, dir, "S1")
  a <- DBI::dbGetQuery(con, "SELECT paths, scaffolds FROM assemble WHERE ID='S1'")
  expect_equal(a$paths, 2)
  expect_equal(a$scaffolds, 3)
})
