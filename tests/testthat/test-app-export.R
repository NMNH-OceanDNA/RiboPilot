# The export dialog's header editor: templates are checked live and saved into
# the selected configuration name.

mk_export_project <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT, rank TEXT,
    match_type TEXT, lineage TEXT, lineage_taxids TEXT, taxonomy_switch INT, taxonomy_lock INT,
    export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT, assemble_switch INT,
    assemble_opts TEXT, blast_opts TEXT, paths INT, scaffolds INT, length TEXT, assemble_notes TEXT, find_rrna_notes TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE preprocess(ID TEXT, pre_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT, length INT,
    sequence TEXT, ignore INT, blast_species TEXT, blast_accession TEXT, blast_pident REAL,
    blast_qcovs REAL, blast_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, problematic TEXT,
    partial TEXT, n_features INT, missing TEXT, n_flags INT, annotate_notes TEXT, itsx_warnings TEXT, export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT,
    type TEXT, start INT, stop INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    min_depth INT, max_depth INT, depths TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL,
    qcovs REAL, evalue REAL, bitscore REAL, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, bitscore REAL, evalue REAL)")
  DBI::dbExecute(con, "CREATE TABLE path_scores(ID TEXT, path INT, length INT, depth REAL, anchors TEXT,
    blast_pident REAL, hit_id TEXT, hit_species TEXT, score REAL, rank INT, auto_select INT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio',7955,'species','exact',NULL,NULL,1,1,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1,2,'default','default',1,1,'600',NULL,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO preprocess VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'linear',600,?,0,'Danio rerio','A1',99,88,'core_nt')",
    params = list(paste(rep("ACGT", 150), collapse = "")))
  # annotation-locked, so the unit is exportable
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,2,1,'default','default','yes','no','no',3,'',0,NULL,NULL,NULL,1)")
  for (f in list(c("18S", 10, 200), c("ITS1", 201, 260), c("5.8S", 261, 420))) {
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,1,?,?,?,?)",
      params = list(f[1], RiboPilot:::rrna_element_type(f[1]),
        as.integer(f[2]), as.integer(f[3])))
  }
  DBI::dbExecute(con, "INSERT INTO coverage (ID,path,scaffold,mean_depth) VALUES('S1',1,1,100)")
  fixture_fill_schema(con)
  DBI::dbDisconnect(con)
  dir
}

# the debounced boxes only settle once the mock clock has moved past the delay
settle <- function(session) {
  session$flushReact()
  session$elapse(600)
  session$flushReact()
}

test_that("the header status lines pass the defaults and fail an unknown column", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(export_btn = 1)
    session$setInputs(
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER)
    settle(session)

    expect_equal(hdr_rec_chk()$level, "ok")
    expect_equal(hdr_mrk_chk()$level, "ok")
    expect_match(output$export_hdr_status$html, "circle-check")

    session$setInputs(fasta_header = "{seqid} {not_a_column}")
    settle(session)
    expect_equal(hdr_rec_chk()$level, "error")
    expect_match(hdr_rec_chk()$message, "not_a_column")
    expect_match(output$export_hdr_status$html, "circle-xmark")
  })
})

test_that("Save template stores the current boxes under the typed name", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(export_btn = 1)
    session$setInputs(
      fasta_header = "{seqid} {Taxon} custom",
      fasta_header_marker = "{seqid}_{element}", export_markers = TRUE)
    settle(session)
    session$setInputs(export_template = "mine")
    settle(session)

    # typing a name alone writes nothing
    expect_equal(export_template(), "mine")
    expect_equal(nrow(RiboPilot:::app_query("SELECT * FROM export_opts WHERE export_opts = 'mine'")), 0L)

    session$setInputs(save_template = 1)
    settle(session)
    row <- RiboPilot:::app_query("SELECT * FROM export_opts WHERE export_opts = 'mine'")
    expect_equal(nrow(row), 1L)
    expect_equal(row$fasta_header, "{seqid} {Taxon} custom")
    expect_equal(row$fasta_header_marker, "{seqid}_{element}")
    # the new name is a second configuration, not a rename of the old one
    expect_true(all(c("default", "mine") %in%
      RiboPilot:::app_query("SELECT export_opts FROM export_opts")$export_opts))
  })
})

test_that("editing a box does not save until Save, and invalid text is never saved", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(export_btn = 1)
    session$setInputs(
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER, export_markers = TRUE)
    settle(session)
    session$setInputs(export_template = "batch1")
    session$setInputs(save_template = 1)
    settle(session)

    session$setInputs(fasta_header = "{seqid} edited {Taxon}")
    settle(session)
    saved <- RiboPilot:::app_query("SELECT fasta_header FROM export_opts WHERE export_opts = 'batch1'")
    expect_equal(saved$fasta_header, RiboPilot:::DEFAULT_FASTA_HEADER)

    session$setInputs(save_template = 2)
    settle(session)
    saved <- RiboPilot:::app_query("SELECT fasta_header FROM export_opts WHERE export_opts = 'batch1'")
    expect_equal(saved$fasta_header, "{seqid} edited {Taxon}")

    # a template naming a column that does not exist is never written
    session$setInputs(fasta_header = "{seqid} {nope}")
    settle(session)
    session$setInputs(save_template = 3)
    settle(session)
    saved <- RiboPilot:::app_query("SELECT fasta_header FROM export_opts WHERE export_opts = 'batch1'")
    expect_equal(saved$fasta_header, "{seqid} edited {Taxon}")
  })
})

test_that("exporting over an existing bundle asks before replacing it", {
  # the group folder is rebuilt from scratch, so a silent overwrite would take
  # an earlier bundle with it
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  DBI::dbDisconnect(con)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)

    # nothing there yet: the export just runs
    session$setInputs(export_go = 1)
    settle(session)
    expect_true(file.exists(file.path(dir, "out", "export", "g1", "g1.fasta")))

    # a marker of the earlier bundle, to prove the second run did not delete it
    keep <- file.path(dir, "out", "export", "g1", "keep_me.txt")
    writeLines("mine", keep)
    session$setInputs(export_go = 2)
    settle(session)
    expect_true(file.exists(keep))          # asked, did not act

    # answering the warning does the export, and the folder is rebuilt
    session$setInputs(export_overwrite_go = 1)
    settle(session)
    expect_false(file.exists(keep))
    expect_true(file.exists(file.path(dir, "out", "export", "g1", "g1.fasta")))
  })
})

test_that("export asks before shipping a unit flagged as a contaminant", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1', reviewed = 'no', curate_flags = 'contaminant'")
  DBI::dbDisconnect(con)
  fasta <- file.path(dir, "out", "export", "g1", "g1.fasta")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)

    session$setInputs(export_go = 1)
    settle(session)
    expect_false(file.exists(fasta))          # asked, did not act

    session$setInputs(export_screen_skip = 1)
    settle(session)
    expect_false(file.exists(fasta))          # the only unit was left out

    session$setInputs(export_go = 2)
    settle(session)
    expect_warning(session$setInputs(export_screen_keep = 1), "contaminant")
    settle(session)
    expect_true(file.exists(fasta))
    smry <- utils::read.csv(file.path(dir, "out", "export", "g1", "g1_summary.csv"))
    expect_match(smry$submission_risk, "contaminant")
  })

  # a curator who marked it Reviewed is not asked again
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET reviewed = 'yes'")
  DBI::dbDisconnect(con)
  unlink(file.path(dir, "out"), recursive = TRUE)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    expect_warning(session$setInputs(export_go = 1), "contaminant")
    settle(session)
    expect_true(file.exists(fasta))
  })
})

test_that("detail modal sample navigation steps through the table order", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$setInputs(ann_details = 1)
    session$setInputs(nav_order = list(ID = c("S1", "S1"), path = c("1", "1"), scaffold = c("1", "2")))
    settle(session)
    expect_equal(sel_unit()$scaffold, 1L)
    html <- output$sample_nav$html
    expect_match(html, "chevron-left")
    expect_match(html, "S1 [1/2]", fixed = TRUE)

    session$setInputs(nav_next = 1)
    settle(session)
    expect_equal(sel_unit()$scaffold, 2L)
    expect_match(output$sample_nav$html, "shinyjs-disabled\" id=\"nav_next\"")
    expect_no_match(output$sample_nav$html, "shinyjs-disabled\" id=\"nav_prev\"")

    session$setInputs(nav_pick = "S1|1|1")
    settle(session)
    expect_equal(sel_unit()$scaffold, 1L)
  })
})

test_that("assigning a group across two assemblers offers a split", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon) VALUES('S2','Danio rerio')")
  DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_lock, assemble_opts) VALUES('S2',1,'default_bait')")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES('S2',1,1,'linear',600,?,0)", params = list(paste(rep("ACGT", 150), collapse = "")))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_lock) VALUES('S2',1,1,1)")
  DBI::dbExecute(con, "INSERT INTO assemble_opts (assemble_opts, assembler) VALUES ('default_bait','bait')")
  DBI::dbDisconnect(con)
  groups <- function() RiboPilot:::app_query(
    "SELECT ID, export_group FROM annotate ORDER BY ID")$export_group

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1:2, group_name = "g")
    session$setInputs(assign_group_btn = 1)
    session$setInputs(assign_group_go = 1)
    session$flushReact()
    # mixed assemblers: nothing written until the warning is answered
    expect_true(all(is.na(groups())))

    session$setInputs(group_split = 1)
    session$flushReact()
    expect_equal(groups(), c("g-getorganelle", "g-bait"))

    # the override keeps one name for everything
    session$setInputs(tbl__reactable__selected = 1:2, group_name = "one")
    session$setInputs(assign_group_go = 2)
    session$setInputs(group_keep_one = 1)
    session$flushReact()
    expect_equal(groups(), c("one", "one"))
  })
})

test_that("export asks before shipping a severe chimera, not a minor one", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1', reviewed = 'no'")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, reason) VALUES ('S1',1,1,'minor','depth step at 300')")
  DBI::dbDisconnect(con)
  fasta <- file.path(dir, "out", "export", "g1", "g1.fasta")
  go <- function(n) shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    session$setInputs(export_go = n)
    settle(session)
  })

  go(1)
  expect_true(file.exists(fasta))          # minor: no question asked

  unlink(file.path(dir, "out"), recursive = TRUE)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE chimera SET label = 'severe', reason = '1-250 Actinopteri; 251-600 Eurotiomycetes'")
  DBI::dbDisconnect(con)
  go(1)
  expect_false(file.exists(fasta))         # severe: asked, did not act
})

test_that("a busy database during Assign Group is reported and the session stays open", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  local_mocked_bindings(set_export_group = function(...) stop("database is locked"),
    .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1, group_name = "g")
    session$setInputs(assign_group_go = 1)
    expect_false(session$isClosed())
    expect_true(is.na(RiboPilot:::app_query("SELECT export_group FROM annotate")$export_group))
  })
})
