# The scaffold table's clicks are keyed on the row's own uid, not its position.
# Previously the handler re-queried with a fixed ORDER BY and indexed by row
# number, so any change to the display order acted on the WRONG scaffold -- and
# the action on the other end DELETEs that unit's annotation and curation.

mk_app_project <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT, rank TEXT,
    match_type TEXT, lineage TEXT, lineage_taxids TEXT, taxonomy_switch INT, taxonomy_lock INT,
    export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT, assemble_switch INT,
    assemble_opts TEXT, blast_opts TEXT, paths INT, scaffolds INT, length TEXT, assemble_notes TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE preprocess(ID TEXT, pre_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT, length INT,
    sequence TEXT, ignore INT, blast_species TEXT, blast_accession TEXT, blast_pident REAL,
    blast_qcovs REAL, blast_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, problematic TEXT,
    partial TEXT, n_features INT, missing TEXT, n_flags INT, annotate_notes TEXT, export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL,
    qcovs REAL, evalue REAL, bitscore REAL, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, bitscore REAL, evalue REAL)")
  DBI::dbExecute(con, "CREATE TABLE path_scores(ID TEXT, path INT, length INT, depth REAL, anchors TEXT,
    blast_pident REAL, hit_id TEXT, hit_species TEXT, score REAL, rank INT, auto_select INT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio',7955,'species','exact',NULL,NULL,1,1,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1,2,'default','default',2,4,'100',NULL,1)")
  DBI::dbExecute(con, "INSERT INTO preprocess VALUES('S1','default')")
  # four units; the SHORTEST is p2.s2, so a length-sorted view reorders the table
  for (u in list(c(1,1,6120), c(1,2,1442), c(2,1,5980), c(2,2,890))) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',?,?,'linear',?,'ACGT',0,'Danio rerio','A1',99,88,'core_nt')",
      params = list(u[1], u[2], u[3]))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',?,?,1,0,'default','default','no','no','no',3,'',0,NULL,NULL,1)",
      params = list(u[1], u[2]))
    DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',?,?,100)", params = list(u[1], u[2]))
  }
  # Bring the fixture up to the current schema rather than hand-listing every
  # column: run_app() migrates on open, so this is what the app actually sees,
  # and the fixture then tracks future migrations for free.
  RiboPilot::migrate_db(con)
  DBI::dbDisconnect(con)
  dir
}

test_that("ignoring by uid acts on that exact scaffold, whatever the row order", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    expect_equal(nrow(scaffold_rows()), 4)

    # p2.s2 is the LAST row by (path, scaffold) but the FIRST by length. A
    # position-keyed handler would disagree with itself here.
    session$setInputs(scaf_ignore = "S1|2|2")
    session$flushReact()

    ig <- RiboPilot:::app_query(
      "SELECT path, scaffold, ignore FROM assemblies WHERE ID='S1' ORDER BY path, scaffold")
    expect_equal(ig$ignore, c(0L, 0L, 0L, 1L))       # only p2.s2
    # its annotate unit is HELD, not deleted, and only that one is held
    left <- RiboPilot:::app_query(
      "SELECT path, scaffold, annotate_switch FROM annotate WHERE ID='S1' ORDER BY path, scaffold")
    expect_equal(nrow(left), 4)
    held <- left$annotate_switch == 0L
    expect_equal(which(held), 4L)
    expect_true(left$path[held] == 2 && left$scaffold[held] == 2)
  })
})

test_that("a uid for a unit that does not exist is ignored, not misapplied", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_ignore = "S1|9|9")
    session$flushReact()
    ig <- RiboPilot:::app_query("SELECT ignore FROM assemblies WHERE ID='S1'")
    expect_true(all(ig$ignore == 0))     # nothing touched
    session$setInputs(scaf_ignore = "garbage")
    session$flushReact()
    expect_true(all(RiboPilot:::app_query("SELECT ignore FROM assemblies WHERE ID='S1'")$ignore == 0))
  })
})

test_that("the scaffold table actually shows its evidence columns", {
  # defaultColDef sets show = FALSE (unit_evidence carries many internal
  # columns), so every column meant to be seen needs an explicit show = TRUE.
  # Missing one is invisible to a "does it render" check.
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    j <- jsonlite::fromJSON(as.character(output$assemble_scaffolds), simplifyVector = FALSE)
    cols <- j$x$tag$attribs$columns
    vis <- vapply(Filter(function(c) !identical(c$show, FALSE), cols),
      function(c) c$id, character(1))
    for (want in c("unit", "length", "depth", "elements", "screen_species",
                   "shares", "odd_one_out", "include", "hits")) {
      expect_true(want %in% vis, info = paste("column not visible:", want))
    }
    # the uid must travel to the browser (the click key) but stay out of sight
    expect_true("uid" %in% vapply(cols, function(c) c$id, character(1)))
    expect_false("uid" %in% vis)
  })
})

test_that("the evidence table exposes a uid for every unit", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  d <- RiboPilot:::unit_evidence("S1")
  expect_equal(d$uid, c("S1|1|1", "S1|1|2", "S1|2|1", "S1|2|2"))
  expect_equal(d$unit, c("p1.s1", "p1.s2", "p2.s1", "p2.s2"))
  # unknown / empty ID yields no rows (typed frame for a real ID, bare for none)
  expect_equal(nrow(RiboPilot:::unit_evidence("nope")), 0L)
  expect_equal(nrow(RiboPilot:::unit_evidence("")), 0L)
})
