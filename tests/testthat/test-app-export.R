# The export dialog's header editor: templates are checked live and saved into
# the selected configuration name.

mk_export_project <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT, rank TEXT,
    match_type TEXT, lineage TEXT, lineage_taxids TEXT, taxonomy_switch INT, taxonomy_lock INT,
    export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT, assemble_switch INT,
    assemble_opts TEXT, blast_opts TEXT, paths INT, scaffolds INT, length TEXT, assemble_notes TEXT, find_rrna_notes TEXT, user_notes TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE preprocess(ID TEXT, pre_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT, length INT,
    sequence TEXT, ignore INT, blast_species TEXT, blast_accession TEXT, blast_pident REAL,
    blast_qcovs REAL, blast_source TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, problematic TEXT,
    partial TEXT, n_features INT, missing TEXT, n_flags INT, annotate_notes TEXT, itsx_warnings TEXT, export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT,
    type TEXT, start INT, stop INT, status TEXT, flags TEXT)")
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
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1,2,'default','default',1,1,'600',NULL,NULL,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO preprocess VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1,'linear',600,?,0,'Danio rerio','A1',99,88,'core_nt')",
    params = list(paste(rep("ACGT", 150), collapse = "")))
  # annotation-locked, so the unit is exportable
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,2,1,'default','default','yes','no','no',3,'',0,NULL,'',NULL,1)")
  for (f in list(c("18S", 10, 200), c("ITS1", 201, 260), c("5.8S", 261, 420))) {
    DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, type, start, stop) VALUES('S1',1,1,?,?,?,?)",
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

test_that("Assign Group refuses a name that cannot be a folder", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1, group_name = "Batch 1")
    session$setInputs(assign_group_btn = 1)
    session$setInputs(assign_group_go = 1)
    session$flushReact()
    expect_true(all(is.na(RiboPilot:::app_query("SELECT export_group FROM annotate")$export_group)))
    session$setInputs(group_name = "Batch_1")
    session$setInputs(assign_group_go = 2)
    session$flushReact()
    expect_equal(RiboPilot:::app_query("SELECT export_group FROM annotate")$export_group, "Batch_1")
  })
})

test_that("locking a sample with no taxon asks first, and locks on Lock anyway", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), ".sqlite")
  DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 0")
  DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon, match_type, taxonomy_switch, taxonomy_lock)
    VALUES ('S2', 'Morus sp.', 'ambiguous', 2, 0)")
  DBI::dbDisconnect(con)
  locks <- function() {
    c2 <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
    on.exit(DBI::dbDisconnect(c2))
    DBI::dbGetQuery(c2, "SELECT taxonomy_lock FROM samples ORDER BY ID")$taxonomy_lock
  }
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Taxonomy")
    session$setInputs(tbl__reactable__selected = 1L)
    session$setInputs(tax_lock_btn = 1)
    expect_equal(locks(), c(1L, 0L))       # a resolved sample locks at once
    session$setInputs(tbl__reactable__selected = 2L)
    session$setInputs(tax_lock_btn = 2)
    expect_equal(locks(), c(1L, 0L))       # no taxon: asked, not locked
    session$setInputs(lock_confirm = TRUE)
    expect_equal(locks(), c(1L, 1L))
  })
})

test_that("Lock asks before locking a flagged unit that is not marked Reviewed", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, reviewed = 'no', curate_flags = 'contaminant'")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  q <- function(sql) RiboPilot:::app_query(sql)[[1]]
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    # Assemble: the sample's kept unit is flagged, so Lock asks first
    session$setInputs(mode = "Assemble")
    session$flushReact()
    d <- samples_data()
    expect_equal(names(d)[match("ID", names(d)) + 1], "checks")
    expect_equal(d$checks, "contaminant")
    session$setInputs(tbl__reactable__selected = 1L, asmb_lock_btn = 1)
    expect_equal(q("SELECT assemble_lock FROM assemble"), 0L)
    session$setInputs(lock_confirm = TRUE)
    expect_equal(q("SELECT assemble_lock FROM assemble"), 1L)

    # Annotation Review toolbar: same question, Cancel leaves it unlocked
    session$setInputs(mode = "Annotate")
    session$flushReact()
    d <- samples_data()
    expect_equal(names(d)[match("scaffold", names(d)) + 1], "checks")
    session$setInputs(tbl__reactable__selected = 1L, annot_lock_btn = 1)
    session$setInputs(lock_confirm = FALSE)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 0L)

    # the unit window asks too, and Lock anyway locks and closes
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(lock_and_close = 1)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 0L)
    session$setInputs(lock_confirm = TRUE)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 1L)

    # a Reviewed unit locks at once
    set("UPDATE annotate SET annotate_lock = 0, reviewed = 'yes'")
    session$setInputs(modal_lock = 1)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 1L)
  })
})

test_that("unlocking a sample in Assemble unlocks its units in Annotation Review", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  testthat::local_mocked_bindings(confirmSweetAlert = function(...) NULL, .package = "shinyWidgets")
  q <- function(sql) RiboPilot:::app_query(sql)[[1]]
  expect_true(all(q("SELECT annotate_lock FROM annotate WHERE ID = 'S1'") == 1))
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1L, asmb_lock_btn = 1)
    session$setInputs(lock_confirm = TRUE)
    expect_equal(q("SELECT assemble_lock FROM assemble WHERE ID = 'S1'"), 0L)
    expect_true(all(q("SELECT annotate_lock FROM annotate WHERE ID = 'S1'") == 0))
  })
})

test_that("a row's lock icon asks about a flagged unit every time", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, reviewed = 'no', curate_flags = 'contaminant'")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  q <- function(sql) RiboPilot:::app_query(sql)[[1]]
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble", lock_confirm_skip = TRUE)
    session$flushReact()
    session$setInputs(lock_row = "1")
    expect_equal(q("SELECT assemble_lock FROM assemble"), 0L)
    session$setInputs(lock_confirm = TRUE)
    expect_equal(q("SELECT assemble_lock FROM assemble"), 1L)
    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(lock_row = "1")
    expect_equal(q("SELECT annotate_lock FROM annotate"), 0L)
  })
})

test_that("Lock does not ask about a flagged unit marked Exclude from export", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, reviewed = 'no', problematic = 'yes',
    curate_flags = 'contaminant'")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  q <- function(sql) RiboPilot:::app_query(sql)[[1]]
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1L, asmb_lock_btn = 1)
    expect_equal(q("SELECT assemble_lock FROM assemble"), 1L)

    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1L, annot_lock_btn = 1)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 1L)

    set("UPDATE annotate SET annotate_lock = 0")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(modal_lock = 1)
    expect_equal(q("SELECT annotate_lock FROM annotate"), 1L)
  })
})

test_that("Exclude & Lock excludes and locks in one click, from the toolbar and the unit window", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, problematic = 'no'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  q <- function(sql) RiboPilot:::app_query(sql)
  set <- function(sql) { c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, sql); DBI::dbDisconnect(c2) }

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    session$flushReact()
    session$setInputs(tbl__reactable__selected = 1L, problematic_btn = 1)
    expect_equal(q("SELECT problematic, annotate_lock FROM annotate"),
      data.frame(problematic = "yes", annotate_lock = 1L))
    # a locked unit is left alone
    session$setInputs(problematic_btn = 2)
    expect_equal(q("SELECT problematic FROM annotate")[[1]], "yes")
    # unlocked, the same button includes it again and leaves it unlocked
    set("UPDATE annotate SET annotate_lock = 0")
    refresh(refresh() + 1); session$flushReact()
    session$setInputs(problematic_btn = 3)
    expect_equal(q("SELECT problematic, annotate_lock FROM annotate"),
      data.frame(problematic = "no", annotate_lock = 0L))

    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(mark_problematic = 1)
    expect_equal(q("SELECT problematic, annotate_lock FROM annotate"),
      data.frame(problematic = "yes", annotate_lock = 1L))
  })
})

test_that("an excluded path does not count against one locked path per sample", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_lock, problematic) VALUES ('S1', 2, 1, 1, 'yes')")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0 WHERE path = 1")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  local_mocked_bindings(mp_alert = function(...) NULL)
  u <- data.frame(ID = "S1", path = 1L, scaffold = 1L)
  expect_true(RiboPilot:::need_one_path(u, session = NULL))
  c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE annotate SET problematic = 'no' WHERE path = 2"); DBI::dbDisconnect(c2)
  expect_false(RiboPilot:::need_one_path(u, session = NULL))
})

# the review layer reads these; fixture_review_schema() adds any the fixture lacks
add_trim_cols <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  fixture_review_schema(con)
}

test_that("a locked unit is read only in the annotation window", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  add_trim_cols(dir)
  withr::local_dir(dir)
  snap <- function() RiboPilot:::app_query(
    "SELECT asm.length, an.problematic, an.annotate_notes,
       (SELECT stop FROM annotations WHERE element = '18S') AS stop18,
       (SELECT COUNT(*) FROM review_units) AS n_copies
     FROM assemblies asm JOIN annotate an ON an.ID = asm.ID")

  # the window says so, and offers Close instead of Lock & Close
  html <- paste(as.character(RiboPilot:::annotate_modal("S1", 1L, 1L)), collapse = " ")
  expect_match(html, "This contig is locked", fixed = TRUE)
  expect_false(grepl("lock_and_close", html, fixed = TRUE))

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    before <- snap()
    active_feature("18S")
    session$setInputs(feat_edit = TRUE)
    session$flushReact()
    session$setInputs(nt_3_in = 1)
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    session$setInputs(mark_problematic = 1)
    session$setInputs(annotate_notes = "changed after sign-off")
    session$flushReact()
    expect_equal(snap(), before)
  })
})

test_that("a trim writes the review copy, and a second trim offers the copy's span", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  add_trim_cols(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0")
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = NULL")
  DBI::dbExecute(con, "UPDATE coverage SET depths = ?", params = list(paste(1:600, collapse = ",")))
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label, junctions) VALUES ('S1', 1, 1, 'severe', '300')")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  alerts <- character(0)
  local_mocked_bindings(rp_alert = function(msg, ...) alerts <<- c(alerts, msg))

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$setInputs(trim_mode = "annotated")
    session$flushReact()
    expect_equal(unname(trim_span()), c(10, 420, 600))
    # the junction hint is a check result, shown only once the contig's checks are done
    expect_error(output$trim_chimera_hint)
    c2 <- RiboPilot:::app_con(); fixture_mark_done(c2); DBI::dbDisconnect(c2)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_equal(output$trim_chimera_hint, "suggested junctions: 300")

    active_feature("18S")
    session$setInputs(feat_edit = TRUE)
    expect_false(is.null(edit_backup()))
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    session$flushReact()
    expect_equal(utils::tail(alerts, 1), "Contig trimmed. Its checks show not run until you press Update.")
    # the trim ends the Edit session, so a later Reset cannot write pre-trim coordinates
    expect_null(edit_backup())
    expect_false(edit_active())
    expect_equal(unname(trim_span()), c(1, 411, 411))
    expect_equal(RiboPilot:::app_query("SELECT trim_start, trim_end, check_switch FROM review_units"),
      data.frame(trim_start = 10L, trim_end = 420L, check_switch = 1L))

    # not run on the copy: no hint; once done, the copy's own junctions
    expect_error(output$trim_chimera_hint)
    c2 <- RiboPilot:::app_con()
    DBI::dbExecute(c2, "UPDATE review_units SET check_switch = 2")
    DBI::dbExecute(c2, "INSERT INTO review_chimera (ID, path, scaffold, label, junctions) VALUES ('S1', 1, 1, 'severe', '291')")
    DBI::dbDisconnect(c2)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_equal(output$trim_chimera_hint, "suggested junctions: 291")

    # the trim window plots the copy: the coverage cut to the window, genes in copy coordinates
    g <- RiboPilot:::plot_trim("S1", 1L, 1L, 1, 411)
    expect_equal(range(g$data$pos), c(1L, 411L))
    expect_equal(g$data$depth[1], 10)
    f <- g$layers[[1]]$data
    expect_equal(c(f$start[f$element == "18S"], f$stop[f$element == "18S"]), c(1L, 191L))

    # a 5 bp flank around the trimmed copy keeps every base of 18S
    session$setInputs(trim_mode = "flank", trim_flank = 5)
    session$flushReact()
    expect_equal(trim_preset(), c(1, 411))
    session$setInputs(trim_from = 1, trim_to = 411, trim_apply = 2)
    session$flushReact()
    expect_equal(RiboPilot:::app_query("SELECT start, stop FROM review_annotations_v WHERE element = '18S'"),
      data.frame(start = 1L, stop = 191L))
    # the original stays as the pipeline left it, for the Assemble step
    expect_equal(RiboPilot:::app_query("SELECT start, stop FROM annotations WHERE element = '18S'"),
      data.frame(start = 10L, stop = 200L))
    expect_equal(RiboPilot:::app_query("SELECT length FROM assemblies")$length, 600L)
  })
})

test_that("Restore deletes the copy and keeps the marks; a contig without a copy has nothing to restore", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, annotate_notes = 'keep me', export_group = 'g1'")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  alerts <- character(0)
  local_mocked_bindings(rp_alert = function(msg, ...) alerts <<- c(alerts, msg))
  asked <- character(0)
  local_mocked_bindings(mp_confirm = function(id, title, text, ...) asked <<- c(asked, text))
  copies <- function() RiboPilot:::app_query("SELECT COUNT(*) AS n FROM review_units")$n
  marks <- function() RiboPilot:::app_query(
    "SELECT reviewed, problematic, partial, annotate_notes, export_group, annotate_lock FROM annotate")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(restore_btn = 1)
    expect_equal(utils::tail(alerts, 1), "This contig has no edits to restore.")
    expect_length(asked, 0)

    # marks write annotate only and never make a copy
    session$setInputs(mark_partial = 1)
    expect_equal(copies(), 0)
    session$setInputs(trim_from = 10, trim_to = 420, trim_apply = 1)
    expect_equal(copies(), 1)
    before <- marks()
    expect_equal(before$partial, "yes")

    active_feature("18S")
    session$setInputs(feat_edit = TRUE)
    session$setInputs(restore_btn = 2)
    expect_match(asked[1], "411 bp now and 600 bp originally", fixed = TRUE)
    expect_match(asked[1], "Marks, notes, group, and lock stay.", fixed = TRUE)
    session$setInputs(restore_confirm = TRUE)
    expect_equal(copies(), 0)
    expect_equal(marks(), before)
    expect_null(edit_backup())
    expect_equal(utils::tail(alerts, 1), "Contig restored to its original (600 bp).")
    expect_equal(unname(trim_span()), c(10, 420, 600))
  })
})

test_that("export needs both locks, and the export window names what it leaves out", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export", export_group_sel = "g1")
    session$flushReact()
    expect_equal(nrow(samples_data()), 0)
    expect_equal(nrow(export_units_df()), 0)
    expect_match(output$export_summary$html, "S1 p1.s1 (assembly not locked)", fixed = TRUE)

    con <- RiboPilot:::app_con()
    DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
    DBI::dbExecute(con, "UPDATE annotate SET problematic = 'yes'")
    DBI::dbDisconnect(con)
    session$setInputs(refresh = 1)
    session$flushReact()
    expect_equal(nrow(export_units_df()), 0)
    expect_match(output$export_summary$html, "S1 p1.s1 (excluded from export)", fixed = TRUE)
  })
})

test_that("annotation options on the Assemble tab answer to the annotation lock", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0, annotate_switch = 2")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  said <- character(0)
  local_mocked_bindings(mp_alert = function(title, text = NULL, ...) said <<- c(said, paste(title, text)))
  switch_now <- function() RiboPilot:::app_query("SELECT annotate_switch FROM annotate")$annotate_switch

  # assembly locked, unit unlocked: the cell is a live link and Save re-queues annotation
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$flushReact()
    j <- jsonlite::fromJSON(as.character(output$tbl), simplifyVector = FALSE)
    cd <- Filter(function(c) identical(c$id, "annotate_opts"), j$x$tag$attribs$columns)[[1]]
    expect_false(grepl("assemble_lock", paste(unlist(cd$cell), collapse = " "), fixed = TRUE))
    session$setInputs(set_annotate_opts = 1)
    session$flushReact()
    session$setInputs(no_name = "default", no_edit = FALSE, update_no = 1)
    session$flushReact()
  })
  expect_equal(switch_now(), 1L)

  # unit annotation-locked: refused, naming that lock and the tab to unlock it in
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1, annotate_switch = 2")
  DBI::dbDisconnect(con)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    session$flushReact()
    session$setInputs(set_annotate_opts = 1)
    session$flushReact()
  })
  expect_equal(switch_now(), 2L)
  expect_true(any(grepl("annotation is locked.*Annotation Review", said)))
})

test_that("export_project skips a group with no unit that has both locks", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g2'")
  DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_lock, assemble_switch) VALUES ('S0', 0, 2)")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, ignore) VALUES ('S0', 1, 1, 0)")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_lock, problematic, export_group)
    VALUES ('S0', 1, 1, 1, 'no', 'g1')")
  DBI::dbDisconnect(con)
  expect_equal(suppressMessages(RiboPilot::export_project(dir)), "g2")
})

test_that("the 5' and 3' partial toggles keep with Save and undo with Reset", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  add_trim_cols(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 0")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  p <- function() as.integer(unlist(RiboPilot:::app_query(
    "SELECT partial_start_user, partial_stop_user FROM review_annotations_v WHERE element = '18S'")))

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    active_feature("18S")
    session$setInputs(toggle_partial_start = 1)   # not editing: ignored
    expect_equal(p(), c(NA_integer_, NA_integer_))
    expect_equal(output$feat_edit_active, "Active feature: 18S")

    session$setInputs(feat_edit = TRUE)
    session$setInputs(toggle_partial_start = 2)
    expect_equal(p(), c(1L, NA_integer_))
    expect_equal(output$feat_edit_active, "Active feature: 18S (partial 5')")
    session$setInputs(feat_reset = 1)
    expect_equal(p(), c(NA_integer_, NA_integer_))

    session$setInputs(feat_edit = FALSE)
    session$setInputs(feat_edit = TRUE)
    session$setInputs(toggle_partial_stop = 1)
    session$setInputs(toggle_partial_start = 3)
    session$setInputs(toggle_partial_start = 4)   # toggled back off
    session$setInputs(feat_save = 1)
    expect_equal(p(), c(0L, 1L))

    # a gene at the scaffold end is automatically partial; a click clears it explicitly
    c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE review_annotations SET start = 1, partial_start_user = NULL WHERE element = '18S'")
    DBI::dbDisconnect(c2)
    edit_tick(edit_tick() + 1)
    session$flushReact()
    expect_equal(output$feat_edit_active, "Active feature: 18S (partial 5', 3')")
    session$setInputs(toggle_partial_start = 5)
    expect_equal(p(), c(0L, 1L))

    # a locked unit is read only
    c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE annotate SET annotate_lock = 1"); DBI::dbDisconnect(c2)
    session$setInputs(toggle_partial_stop = 2)
    expect_equal(p(), c(0L, 1L))
  })
  # the choices live on the review copy; the original feature never changes
  expect_equal(RiboPilot:::app_query("SELECT start, partial_start, partial_stop FROM annotations WHERE element = '18S'"),
    data.frame(start = 10L, partial_start = NA_integer_, partial_stop = NA_integer_))
})

test_that("the Export window offers token chips, and choosing fields returns to the window", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  DBI::dbExecute(con, "UPDATE samples SET GEOME_BCID = 'ark:/1/A'")
  DBI::dbExecute(con, "INSERT INTO meta_export_fields VALUES ('geome:combo:lat_lon')")
  DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "S1", source = "GEOME", level = "Event", depth = 0L, ref = "ark:/1/A",
    field = c("decimalLatitude", "decimalLongitude"), value = c("17.5", "-149.8")))
  DBI::dbDisconnect(con)
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    html <- as.character(export_modal_current())
    expect_match(html, "Available metadata", fixed = TRUE)
    expect_match(html, "mp-token-list", fixed = TRUE)
    expect_match(html, "data-target=\"fasta_header\"", fixed = TRUE)
    expect_match(html, "data-insert=\"{Taxon}\"", fixed = TRUE)
    expect_match(html, "data-insert=\"[lat_lon={geome_lat_lon}]\"", fixed = TRUE)
    expect_match(html, "data-totals=\"{&quot;g1&quot;:1}\"", fixed = TRUE)
    expect_match(html, "id=\"token_fields_geome\"", fixed = TRUE)

    session$setInputs(fasta_header = "{seqid} edited", export_markers = FALSE, export_template = "default",
      export_group_sel = "g1")
    session$setInputs(token_fields_geome = 1)
    expect_equal(export_snap()$fasta_header, "{seqid} edited")
    back <- as.character(export_modal_current())
    expect_match(back, "{seqid} edited", fixed = TRUE)
    expect_null(export_snap())
    session$setInputs(token_fields_gbif = 1)
    expect_false(is.null(export_snap()))
    session$setInputs(meta_view_closed = 1)
    expect_null(export_snap())
  })
})

test_that("a conflicting {geome_lat_lon} asks before exporting, and Export anyway proceeds", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  DBI::dbExecute(con, "UPDATE samples SET GEOME_BCID = 'ark:/1/A'")
  DBI::dbExecute(con, "INSERT INTO meta_export_fields VALUES ('geome:combo:lat_lon')")
  add <- function(source, field, value) DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "S1", source = source, level = if (source == "GBIF") "Occurrence" else "Event", depth = 0L, ref = "r", field = field, value = value))
  add("GEOME", c("decimalLatitude", "decimalLongitude"), c("10.5", "20.5"))
  add("GBIF", c("decimalLatitude", "decimalLongitude"), c("-33.1", "151.2"))
  DBI::dbDisconnect(con)
  asked <- character(0)
  local_mocked_bindings(mp_confirm = function(id, title, ...) asked <<- c(asked, title))
  fasta <- file.path(dir, "out", "export", "g1", "g1.fasta")

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = "{seqid} [lat_lon={geome_lat_lon}]",
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    session$setInputs(export_go = 1)
    settle(session)
    expect_equal(asked, "Sample metadata disagrees")
    expect_false(file.exists(fasta))

    session$setInputs(specimen_confirm = TRUE)
    settle(session)
    settle(session)
    expect_true(file.exists(fasta))
  })

  unlink(file.path(dir, "out"), recursive = TRUE)
  asked <- character(0)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1",
      fasta_header = RiboPilot:::DEFAULT_FASTA_HEADER,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    session$setInputs(export_go = 1)
    settle(session)
    expect_length(asked, 0)
    expect_true(file.exists(fasta))
  })
})

test_that("Cancel on the metadata disagreement window brings the Export window back with the edited header", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  DBI::dbExecute(con, "UPDATE samples SET GEOME_BCID = 'ark:/1/A'")
  DBI::dbExecute(con, "INSERT INTO meta_export_fields VALUES ('geome:combo:lat_lon')")
  add <- function(source, field, value) DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "S1", source = source, level = if (source == "GBIF") "Occurrence" else "Event", depth = 0L, ref = "r", field = field, value = value))
  add("GEOME", c("decimalLatitude", "decimalLongitude"), c("10.5", "20.5"))
  add("GBIF", c("decimalLatitude", "decimalLongitude"), c("-33.1", "151.2"))
  DBI::dbDisconnect(con)
  local_mocked_bindings(mp_confirm = function(...) invisible())
  shown <- character(0)
  local_mocked_bindings(showModal = function(ui, ...) shown <<- c(shown, as.character(ui)))
  hdr <- "{seqid} [lat_lon={geome_lat_lon}]"

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1", export_template = "default", fasta_header = hdr,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    session$setInputs(export_go = 1)
    settle(session)
    shown <<- character(0)
    session$setInputs(specimen_confirm = FALSE)
    settle(session)
    expect_length(shown, 1)
    expect_match(shown, "[lat_lon={geome_lat_lon}]", fixed = TRUE)
    expect_false(file.exists(file.path(dir, "out", "export", "g1", "g1.fasta")))
  })
})

test_that("NMNH user: status, cell edit, and Export anyway", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE annotate SET export_group = 'g1'")
  # real projects key samples on ID; the CSV upload upserts on it
  ddl <- DBI::dbGetQuery(con, "SELECT sql FROM sqlite_master WHERE name = 'samples'")$sql
  DBI::dbExecute(con, "ALTER TABLE samples RENAME TO samples_old")
  DBI::dbExecute(con, sub("ID TEXT", "ID TEXT PRIMARY KEY", ddl, fixed = TRUE))
  DBI::dbExecute(con, "INSERT INTO samples SELECT * FROM samples_old")
  DBI::dbDisconnect(con)
  local_mocked_bindings(.gbif_get = function(path) stop("offline"))
  tpl <- RiboPilot:::nmnh_template_add(RiboPilot:::DEFAULT_FASTA_HEADER)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(export_group_sel = "g1", fasta_header = tpl,
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER,
      export_markers = FALSE, export_report = FALSE)
    settle(session)
    expect_false(nmnh_res()$ok)
    expect_match(output$nmnh_status$html, "1 of 1 sample needs attention")
    expect_true(hdr_rec_chk()$ok)
    expect_match(output$nmnh_copy$html, "Copy table")
    session$setInputs(nmnh_report = 1)
    expect_match(as.character(output$nmnh_table), "nmnh_edit")
    csv <- utils::read.csv(output$nmnh_csv)
    expect_equal(names(csv), c("ID", "Taxon", "specimen_voucher", "voucherURI"))

    # Export with problems asks first
    session$setInputs(export_go = 1)
    settle(session)
    expect_false(file.exists(file.path(dir, "out", "export", "g1", "g1.fasta")))

    session$setInputs(nmnh_edit = list(id = "S1", field = "voucher", value = "USNM:FISH:7"))
    settle(session)
    expect_equal(nmnh_res()$nmnh_specimen_voucher, "USNM:FISH:7")
    v <- app_query("SELECT specimen_voucher, specimen_voucher_source FROM nmnh_vouchers WHERE ID = 'S1'")
    expect_equal(v$specimen_voucher, "USNM:FISH:7")
    expect_equal(v$specimen_voucher_source, "entered")
    expect_false("specimen_voucher" %in% names(app_query("SELECT * FROM samples")))
    # a fixed row stays in the report with its check; an invalid entry stays visible
    tbl <- as.character(output$nmnh_table)
    expect_match(tbl, "USNM:FISH:7", fixed = TRUE)
    expect_match(tbl, "circle-check", fixed = TRUE)
    session$setInputs(nmnh_edit = list(id = "S1", field = "uri", value = "not-an-ark"))
    settle(session)
    tbl <- as.character(output$nmnh_table)
    expect_match(tbl, "not-an-ark", fixed = TRUE)
    expect_match(tbl, "circle-xmark", fixed = TRUE)

    f <- file.path(dir, "up.csv")
    utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio", specimen_voucher = "USNM:FISH:8",
      voucherURI = ""), f, row.names = FALSE)
    session$setInputs(nmnh_upload = list(datapath = f))
    settle(session)
    expect_equal(nmnh_res()$nmnh_specimen_voucher, "USNM:FISH:8")
    expect_equal(app_query("SELECT specimen_voucher_source FROM nmnh_vouchers WHERE ID = 'S1'")[[1]], "upload")

    session$setInputs(nmnh_confirm = TRUE)
    settle(session)
    h <- readLines(file.path(dir, "out", "export", "g1", "g1.fasta"))[1]
    expect_match(h, "[specimen_voucher=USNM:FISH:8] Danio", fixed = TRUE)
    expect_false(grepl("voucherURI", h))
  })
})

test_that("turning on marker export after NMNH adds the NMNH fields to the marker header", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  withr::local_dir(dir)
  local_mocked_bindings(.gbif_get = function(path) stop("offline"))
  sent <- list()
  local_mocked_bindings(updateTextAreaInput = function(session, inputId, value = NULL, ...) {
    sent[[inputId]] <<- value
  })
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Export")
    session$setInputs(nmnh_user = TRUE, export_markers = FALSE,
      fasta_header = RiboPilot:::nmnh_template_add(RiboPilot:::DEFAULT_FASTA_HEADER),
      fasta_header_marker = RiboPilot:::DEFAULT_FASTA_HEADER_MARKER)
    session$setInputs(export_markers = TRUE)
    expect_true(RiboPilot:::nmnh_template_on(sent$fasta_header_marker))
  })
})

test_that("the export fixture answers every review view", {
  skip_if_not_installed("withr")
  dir <- mk_export_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  for (v in names(RiboPilot:::REVIEW_VIEWS_DDL)) {
    expect_no_error(DBI::dbGetQuery(con, sprintf("SELECT ID, path, scaffold FROM %s LIMIT 1", v)))
  }
  # locked, annotated, and checked by ITSx: done
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v"),
    data.frame(has_copy = 0L, length = 600L, check_switch = 2L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_annotations_v")$n, 3L)
})
