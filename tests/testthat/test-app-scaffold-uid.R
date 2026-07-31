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
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    min_depth INT, max_depth INT, depths TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL,
    qcovs REAL, evalue REAL, bitscore REAL, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, bitscore REAL, evalue REAL)")
  DBI::dbExecute(con, "CREATE TABLE path_scores(ID TEXT, path INT, length INT, depth REAL, anchors TEXT,
    blast_pident REAL, hit_id TEXT, hit_species TEXT, score REAL, rank INT, auto_select INT)")
  DBI::dbExecute(con, "INSERT INTO path_scores (ID,path,score,rank,auto_select) VALUES('S1',1,0.912,1,1)")
  DBI::dbExecute(con, "INSERT INTO path_scores (ID,path,score,rank,auto_select) VALUES('S1',2,0.431,2,0)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio',7955,'species','exact',NULL,NULL,1,1,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1,2,'default','default',2,4,'100',NULL,1)")
  DBI::dbExecute(con, "INSERT INTO preprocess VALUES('S1','default')")
  # four units; the SHORTEST is p2.s2, so a length-sorted view reorders the table
  for (u in list(c(1,1,6120), c(1,2,1442), c(2,1,5980), c(2,2,890))) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',?,?,'linear',?,'ACGT',0,'Danio rerio','A1',99,88,'core_nt')",
      params = list(u[1], u[2], u[3]))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',?,?,1,0,'default','default','no','no','no',3,'',0,NULL,NULL,1)",
      params = list(u[1], u[2]))
    DBI::dbExecute(con, "INSERT INTO coverage (ID,path,scaffold,mean_depth) VALUES('S1',?,?,100)",
      params = list(u[1], u[2]))
  }
  # Bring the fixture up to the current schema rather than hand-listing every
  # column: run_app() migrates on open, so this is what the app actually sees,
  # and the fixture then tracks future migrations for free.
  RiboPilot::migrate_db(con)
  DBI::dbDisconnect(con)
  dir
}

test_that("triaging by uid acts on that exact scaffold, whatever the row order", {
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
    session$setInputs(scaf_triage = "S1|2|2|drop")
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

test_that("a uid for a unit that does not exist is dropped, not misapplied", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_triage = "S1|9|9|drop")
    session$flushReact()
    ig <- RiboPilot:::app_query("SELECT ignore FROM assemblies WHERE ID='S1'")
    expect_true(all(ig$ignore == 0))     # nothing touched
    session$setInputs(scaf_triage = "garbage")
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
                   "shares", "odd_one_out", "triage", "path_score")) {
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

test_that("bulk triage previews before it writes, and applying honours the preview", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  # make every unit droppable: no features, and a non-rDNA top hit
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "DELETE FROM annotations")
  for (u in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID,path,scaffold,species,is_nuclear_rrna) VALUES('S1',?,?,'Homo sapiens',0)",
      params = list(u[1], u[2]))
  }
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()

    # opening the modal must not change anything
    session$setInputs(triage_btn = 1, triage_overwrite_user = FALSE)
    session$flushReact()
    expect_true(all(RiboPilot:::app_query("SELECT COALESCE(ignore,0) i FROM assemblies WHERE ID='S1'")$i == 0))
    prev <- output$triage_preview$html      # renderUI returns list(html, deps)
    expect_true(nzchar(prev))
    expect_match(prev, "Apply would change 4 of 4", fixed = TRUE)
    expect_match(prev, "4 newly dropped", fixed = TRUE)
    expect_match(prev, "Resulting state across 1 sample(s): 4 dropped", fixed = TRUE)
    expect_match(prev, "NO scaffolds kept")  # the outcome that silently loses a sample

    session$setInputs(triage_apply = 1)
    session$flushReact()
    r <- RiboPilot:::app_query(
      "SELECT triage, triage_source, COALESCE(ignore,0) i FROM assemblies WHERE ID='S1'")
    expect_true(all(r$triage == "drop"))
    expect_true(all(r$triage_source == "auto"))
    expect_true(all(r$i == 1L))
  })
})

test_that("the modal previews the stored preset, not the shipped defaults", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  # same droppable evidence, but the project's preset says flag, not drop. No
  # triage_rule_* input is ever set here, standing in for the window before the
  # browser echoes the controls back: the preview must already honour the preset.
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "DELETE FROM annotations")
  for (u in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID,path,scaffold,species,is_nuclear_rrna) VALUES('S1',?,?,'Homo sapiens',0)",
      params = list(u[1], u[2]))
  }
  RiboPilot:::write_triage_opts(con, list(no_feat_non_rdna = "flag", odd_one_out = "off"))
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(triage_btn = 1, triage_overwrite_user = FALSE)
    session$flushReact()

    expect_equal(triage_actions()$no_feat_non_rdna, "flag")
    expect_equal(triage_actions()$odd_one_out, "off")
    prev <- output$triage_preview$html
    expect_match(prev, "0 newly dropped", fixed = TRUE)
    expect_match(prev, "Resulting state across 1 sample(s): 0 dropped", fixed = TRUE)
    expect_match(prev, "No rule is set to Drop", fixed = TRUE)
  })
})

test_that("preview says so when a Drop rule is armed but everything it matches is hand-decided", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  # same droppable evidence as above, but every unit is already a user
  # decision: with overwrite off, auto-triage must skip all of them
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "DELETE FROM annotations")
  for (u in list(c(1,1), c(1,2), c(2,1), c(2,2))) {
    DBI::dbExecute(con, "INSERT INTO blast_screen (ID,path,scaffold,species,is_nuclear_rrna) VALUES('S1',?,?,'Homo sapiens',0)",
      params = list(u[1], u[2]))
    DBI::dbExecute(con,
      "UPDATE assemblies SET triage = 'keep', triage_source = 'user' WHERE ID='S1' AND path=? AND scaffold=?",
      params = list(u[1], u[2]))
  }
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(triage_btn = 1, triage_overwrite_user = FALSE)
    session$flushReact()
    prev <- output$triage_preview$html
    expect_true(nzchar(prev))
    expect_match(prev, "nothing in this selection will be dropped", fixed = TRUE)
  })
})

test_that("turning every rule off says Apply changes nothing, not 'Would drop'", {
  # the reported bug: with every rule off, the preview kept saying "Would drop
  # N" because it reported the END STATE (which includes standing hand-decided
  # drops), not what Apply itself would touch
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  # p1.s1 hand-decided drop, p1.s2 hand-decided keep: standing curator calls
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage='drop', triage_source='user' WHERE ID='S1' AND path=1 AND scaffold=1")
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage='keep', triage_source='user' WHERE ID='S1' AND path=1 AND scaffold=2")
  # p2.s1 and p2.s2 already auto-triaged to keep by an earlier run, matching
  # exactly what every rule off would propose again
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage='keep', triage_source='auto' WHERE ID='S1' AND path=2 AND scaffold=1")
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage='keep', triage_source='auto' WHERE ID='S1' AND path=2 AND scaffold=2")
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(triage_btn = 1, triage_overwrite_user = FALSE)
    off_inputs <- stats::setNames(as.list(rep("off", length(RiboPilot:::triage_rule_keys()))),
      paste0("triage_rule_", RiboPilot:::triage_rule_keys()))
    do.call(session$setInputs, off_inputs)
    session$flushReact()

    expect_true(all(vapply(RiboPilot:::triage_rule_keys(),
      function(k) triage_actions()[[k]] == "off", logical(1))))
    prev <- output$triage_preview$html
    expect_true(nzchar(prev))
    expect_match(prev, "Apply would change nothing", fixed = TRUE)
    expect_false(grepl("Would drop", prev, fixed = TRUE))
  })
})

test_that("the all-off message replaces the no-Drop note instead of also printing it", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(triage_btn = 1, triage_overwrite_user = FALSE)
    off_inputs <- stats::setNames(as.list(rep("off", length(RiboPilot:::triage_rule_keys()))),
      paste0("triage_rule_", RiboPilot:::triage_rule_keys()))
    do.call(session$setInputs, off_inputs)
    session$flushReact()

    prev <- output$triage_preview$html
    expect_match(prev, "Every rule is off, so Apply would change nothing.", fixed = TRUE)
    expect_false(grepl("No rule is set to Drop", prev, fixed = TRUE))
  })
})

test_that("Copy FASTA copies the ticked scaffolds, keyed on the table's own selection", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  copied <- NULL
  # intercept the browser-side write; the assertion is about what we hand it
  testthat::local_mocked_bindings(copy_to_clipboard = function(text) {
    copied <<- text; invisible(text)
  })

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()

    # nothing ticked: warn, do not put an empty blob on the clipboard
    session$setInputs(copy_fasta = 1)
    session$flushReact()
    expect_null(copied)

    d <- scaffold_rows()
    pick <- which(d$path == 2 & d$scaffold == 1)
    session$setInputs(`assemble_scaffolds__reactable__selected` = pick)
    session$setInputs(copy_fasta = 2)
    session$flushReact()

    expect_false(is.null(copied))
    expect_match(copied, ">S1 Danio rerio p2.s1", fixed = TRUE)
    expect_false(grepl("p1.s1", copied, fixed = TRUE))
  })
})

test_that("the scaffold table actually shows its tick-box column", {
  # the table hides columns by default (defaultColDef show = FALSE) and reactable
  # applies that to its built-in selection column too, so the checkboxes were
  # created and then hidden: the button was there with no way to pick a row
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()

    # testServer hands back the widget already rendered to JSON
    w <- jsonlite::fromJSON(output$assemble_scaffolds, simplifyVector = FALSE)
    cols <- w$x$tag$attribs$columns
    sel <- Filter(function(c) identical(c$id, ".selection"), cols)
    expect_length(sel, 1)
    expect_true(isTRUE(sel[[1]]$show))
  })
})

test_that("clicking a row examines that scaffold in the panels below", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()

    session$setInputs(scaf_examine = "S1|2|2")
    session$flushReact()
    expect_equal(sel_hit_unit(), list(path = 2L, scaffold = 2L))
    expect_match(output$unit_screen$html, "Scaffold p2.s2", fixed = TRUE)
    expect_match(output$unit_coverage_label$html, "Scaffold p2.s2", fixed = TRUE)

    # keyed on the row's own uid, so display order cannot misdirect it
    session$setInputs(scaf_examine = "S1|1|2")
    session$flushReact()
    expect_equal(sel_hit_unit(), list(path = 1L, scaffold = 2L))
  })
})

test_that("ticking rows for Copy FASTA does not move the examined scaffold", {
  # the two controls mean different things; that is the whole point of the split
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(scaf_examine = "S1|1|1")
    session$flushReact()

    session$setInputs(`assemble_scaffolds__reactable__selected` = c(3L, 4L))
    session$flushReact()

    expect_equal(sel_hit_unit(), list(path = 1L, scaffold = 1L))
  })
})

test_that("the tick column is headed with the copy icon, and the hits link is gone", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()

    w <- jsonlite::fromJSON(output$assemble_scaffolds, simplifyVector = FALSE)
    cols <- w$x$tag$attribs$columns
    sel <- Filter(function(c) identical(c$id, ".selection"), cols)[[1]]
    expect_match(sel$header, "clipboard")
    # the old per-row link was the only way in and read as an action, not a picker
    expect_length(Filter(function(c) identical(c$id, "hits"), cols), 0)
  })
})

test_that("the scaffold table carries each unit's path score", {
  # the score is computed per candidate path, so every scaffold of a path shows
  # the same number; it replaced the separate "Candidate paths" panel.
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    d <- scaffold_rows()
    expect_equal(d$path_score[d$path == 1], c(0.912, 0.912))
    expect_equal(d$path_score[d$path == 2], c(0.431, 0.431))
  })
})

test_that("the scaffold table does not report topology", {
  # GetOrganelle labels an rDNA result "circular" when it disentangles the graph
  # as a cycle, which a tandem array does by construction. Nuclear rDNA is not a
  # circular molecule, so the column asserted something false about every
  # successfully assembled repeat.
  q <- paste(deparse(RiboPilot:::unit_evidence), collapse = " ")
  expect_false(grepl("a.topology", q, fixed = TRUE))
  srv <- paste(deparse(RiboPilot:::app_server), collapse = " ")
  expect_false(grepl('name = "Topology"', srv, fixed = TRUE))
})

test_that("every visible scaffold column explains itself on hover", {
  # the help lives in SCAFFOLD_COL_HELP; a column added without a tooltip is
  # invisible to a plain "does it render" check.
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    j <- jsonlite::fromJSON(as.character(output$assemble_scaffolds), simplifyVector = FALSE)
    cols <- Filter(function(c) !identical(c$show, FALSE), j$x$tag$attribs$columns)
    for (cd in cols) {
      if (identical(cd$id, ".selection")) next     # its own icon carries the hint
      expect_true(is.character(cd$header) && grepl("title=", cd$header, fixed = TRUE),
        info = paste("no hover help on column:", cd$id))
    }
  })
})
