# The main table at scale: filters, refreshes, and what bulk actions reach.

# Four taxonomy-locked samples for the Assemble table: S1 assembly-locked and
# finished, S2 and S4 unlocked and finished, S3 unlocked and failed.
mk_table_project <- function(env = parent.frame()) {
  dir <- mk_app_project(env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  for (id in c("S2", "S3", "S4")) {
    DBI::dbExecute(con, "INSERT INTO samples (ID, Taxon, taxonomy_switch, taxonomy_lock) VALUES (?, 'Danio rerio', 2, 1)",
      params = list(id))
    DBI::dbExecute(con, "INSERT INTO assemble (ID, assemble_lock, assemble_switch, assemble_opts, blast_opts, time_stamp)
      VALUES (?, 0, ?, 'default', 'default', 1)", params = list(id, if (id == "S3") 3L else 2L))
    DBI::dbExecute(con, "INSERT INTO preprocess VALUES (?, 'default')", params = list(id))
  }
  DBI::dbDisconnect(con)
  dir
}

# every updateReactable() call the server makes, newest last
mock_table_updates <- function(env = parent.frame()) {
  calls <- new.env(); calls$list <- list()
  testthat::local_mocked_bindings(
    updateReactable = function(outputId, data = NULL, selected = NULL, page = NULL, ...) {
      calls$list[[length(calls$list) + 1]] <- list(data = data, selected = selected, page = page)
    }, .package = "reactable", .env = env)
  calls
}

test_that("the Lock and State filters take rows out of the table, not just out of view", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    expect_equal(sort(samples_data()$ID), c("S1", "S2", "S3", "S4"))
    session$setInputs(Assemble_state_filter = "3")
    expect_equal(samples_data()$ID, "S3")
    session$setInputs(Assemble_state_filter = RiboPilot:::STATE_CHOICES, Assemble_lock_filter = "0")
    expect_equal(sort(samples_data()$ID), c("S2", "S3", "S4"))
    expect_equal(nrow(all_rows()), 4)
  })
})

test_that("a refresh keeps the page and re-selects the same samples by ID", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  up <- mock_table_updates()
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    invisible(output$tbl)
    d <- samples_data()
    i <- match(c("S3", "S4"), d$ID)
    session$setInputs(tbl__reactable__selected = i, tbl__reactable__page = 3L)
    # S1 leaves the Assemble table, so every later row moves up one
    c2 <- RiboPilot:::app_con(); DBI::dbExecute(c2, "UPDATE samples SET taxonomy_lock = 0 WHERE ID = 'S1'")
    DBI::dbDisconnect(c2)
    session$setInputs(refresh = 1)
    last <- up$list[[length(up$list)]]
    expect_equal(last$page, 3L)
    expect_equal(sort(last$data$ID[last$selected]), c("S3", "S4"))
  })
})

test_that("bulk actions skip rows the search or a column filter hides, and the selection drops them", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  up <- mock_table_updates()
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    d <- samples_data()
    i <- match(c("S2", "S4"), d$ID)
    session$setInputs(tbl__reactable__selected = i)
    session$setInputs(tbl_visible = list(rows = i[1]))   # the search now hides S4
    expect_equal(selected_ids(), "S2")
    expect_equal(up$list[[length(up$list)]]$selected, i[1])
    session$setInputs(tbl_visible = list(rows = integer(0)))   # nothing matches
    expect_equal(selected_ids(), character(0))
    session$setInputs(tbl_visible = list(rows = i[1]))
    session$setInputs(asmb_lock_btn = 1)
    expect_equal(RiboPilot:::app_query("SELECT ID FROM assemble WHERE assemble_lock = 1 ORDER BY ID")$ID,
      c("S1", "S2"))
  })
})

test_that("bulk triage acts on the unlocked samples selected when it opened", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  seen <- character(0)
  testthat::local_mocked_bindings(auto_triage_sample = function(con, id, ...) {
    seen <<- c(seen, id); 0L
  }, .package = "RiboPilot")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    d <- samples_data()
    session$setInputs(tbl__reactable__selected = match(c("S1", "S2"), d$ID))
    session$setInputs(triage_btn = 1)
    expect_equal(triage_ids(), "S2")                 # S1 is assembly-locked
    session$setInputs(tbl__reactable__selected = integer(0), refresh = 1)
    session$setInputs(triage_apply = 1)
    expect_equal(seen, "S2")
  })
})

test_that("a refresh right after a re-render does not carry the old table's selection", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  up <- mock_table_updates()
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Taxonomy")
    invisible(output$tbl)
    session$setInputs(tbl__reactable__selected = c(1L, 2L), tbl__reactable__page = 1L)
    session$setInputs(mode = "Assemble")
    invisible(output$tbl)
    session$setInputs(refresh = 1)
    last <- up$list[[length(up$list)]]
    expect_length(last$selected, 0)
    expect_null(last$page)
    session$setInputs(tbl__reactable__selected = 2L)
    session$setInputs(refresh = 2)
    last <- up$list[[length(up$list)]]
    expect_equal(last$data$ID[last$selected], samples_data()$ID[2])
  })
})

test_that("a row's lock icon acts on that row alone, with the same prompts as Lock", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  asked <- 0L
  testthat::local_mocked_bindings(confirmSweetAlert = function(...) asked <<- asked + 1L,
    .package = "shinyWidgets")
  lk <- function() RiboPilot:::app_query("SELECT ID FROM assemble WHERE assemble_lock = 1 ORDER BY ID")$ID
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble", Assemble_lock_filter = "0")
    row <- function(id) match(id, samples_data()$ID)
    session$setInputs(tbl__reactable__selected = row(c("S2", "S3")))
    session$setInputs(tbl_visible = list(rows = row("S4")))   # the search hides S2 and S3
    session$setInputs(lock_row = as.character(row("S4")))
    expect_equal(lk(), c("S1", "S4"))
    expect_equal(asked, 0L)
    session$setInputs(Assemble_lock_filter = c("0", "1"), tbl_visible = NULL)
    session$setInputs(lock_row = as.character(row("S4")))     # unlocking asks
    expect_equal(c(asked, length(lk())), c(1L, 2L))
    session$setInputs(lock_confirm = TRUE)
    expect_equal(lk(), "S1")
    session$setInputs(lock_row = as.character(row("S1")))
    session$setInputs(lock_confirm_skip = TRUE, lock_confirm = FALSE)
    expect_equal(c(asked, length(lk())), c(2L, 1L))
    session$setInputs(lock_row = as.character(row("S1")))     # not asked again
    expect_equal(c(asked, length(lk())), c(2L, 0L))
  })
})

test_that("a row's state icon sets the state of that row alone", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    row <- function(id) match(id, samples_data()$ID)
    session$setInputs(tbl__reactable__selected = row(c("S2", "S3")))
    session$setInputs(state_row = as.character(row("S4")))
    session$setInputs(new_state = "0", update_state = 1)
    expect_equal(RiboPilot:::app_query("SELECT ID FROM assemble WHERE assemble_switch = 0")$ID, "S4")
  })
})

test_that("the Metadata column sits after Taxon in every step that lists samples, and its helpers stay hidden", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "UPDATE assemble SET assemble_lock = 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_lock = 1, export_group = 'g1'")
  DBI::dbDisconnect(con)
  shiny::testServer(RiboPilot:::app_server, {
    for (m in c("Taxonomy", "Assemble", "Annotate", "Export")) {
      session$setInputs(mode = m)
      session$flushReact()
      j <- jsonlite::fromJSON(as.character(output$tbl), simplifyVector = FALSE)
      cols <- j$x$tag$attribs$columns
      id <- vapply(cols, function(c) c$id %||% "", "")
      expect_equal(id[match("Taxon", id) + 1], "specimen", info = m)
      sp <- cols[[match("specimen", id)]]
      expect_equal(sp$name, "Metadata", info = m)
      expect_false(isFALSE(sp$show), info = m)
      expect_false(any(vapply(cols[id %in% c("specimen_message", "specimen_icons")],
        function(c) !isFALSE(c$show), TRUE)), info = m)
    }
  })
})

test_that("a shown metadata field is a column in every non-Reference view", {
  skip_if_not_installed("withr")
  dir <- mk_table_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  RiboPilot:::.meta_ensure_tables(con)
  DBI::dbAppendTable(con, "meta_records", data.frame(
    ID = "S1", source = "GBIF", level = "Occurrence", depth = 0L, ref = "r",
    field = "sex", value = "Female"))
  f <- RiboPilot:::meta_view_fields(con)
  RiboPilot:::meta_view_save(con, f, "gbif:raw:Occurrence:sex", FALSE)
  col <- f$col[f$key == "gbif:raw:Occurrence:sex"]
  DBI::dbDisconnect(con)
  shiny::testServer(RiboPilot:::app_server, {
    for (m in c("Taxonomy", "Assemble", "Annotate")) {
      session$setInputs(mode = m)
      d <- all_rows()
      expect_true(col %in% names(d), info = m)
      expect_equal(d[[col]][d$ID == "S1"][1], "Female", info = m)
      if (m != "Annotate") expect_equal(d[[col]][d$ID == "S2"][1], "", info = m)
    }
    session$setInputs(mode = "Reference")
    expect_false(col %in% names(all_rows()))
  })
})

test_that("Annotation Review reads the review layer and shows the state of unfinished checks", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  # p1.s1 done with an ITSx flag; p1.s2 still queued for the main pipeline;
  # p2.s1 a copy whose recheck failed; p2.s2 a rechecked copy
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, itsx_warnings = 'itsx_chimeric'
    WHERE path = 1 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 2, itsx_warnings = '' WHERE path = 2")
  copy <- function(s, len, sw, note, species) DBI::dbExecute(con, "INSERT INTO review_units (ID, path,
      scaffold, sequence, length, trim_start, trim_end, edited, check_switch, check_note, rev,
      n_features, itsx_warnings, blast_species)
    VALUES ('S1', 2, ?, 'ACGT', ?, 1, ?, 1, ?, ?, CAST(? AS INTEGER), 3, '', ?)",
    params = list(s, len, len, sw, note, "1760000000000", species))
  copy(1L, 500L, 3L, "tool crashed", "Copy sp")
  copy(2L, 800L, 2L, NA, "Danio rerio")
  # p2.s2's original screen and chimera results differ from its copy's
  DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 2, 2, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1', 2, 2, 'severe')")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 2, 2, 'on_target')")
  DBI::dbExecute(con, "INSERT INTO feature_blast (ID, path, scaffold, element, status, lineage)
    VALUES ('S1', 1, 1, 'ITS1', 'no_hit', NULL), ('S1', 1, 1, '28S', 'too_small', NULL),
           ('S1', 1, 1, '18S', 'ok', NULL)")
  DBI::dbDisconnect(con)
  withr::local_dir(dir)
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    d <- all_rows()
    d <- d[order(d$path, d$scaffold), ]
    expect_equal(d$feature_blast_warnings,
      c("ITS1: no hit; 28S: no BLAST, feature too small", "", "", ""))
    expect_equal(match("feature_blast_warnings", names(d)), match("chimera_reason", names(d)) + 1L)
    expect_equal(d$length, c(6120L, 1442L, 500L, 800L))
    expect_equal(d$blast_species, c("Danio rerio", "Danio rerio", "Copy sp", "Danio rerio"))
    expect_equal(d$screen_verdict, c(NA, NA, NA, "on_target"))
    expect_true(all(is.na(d$chimera_label)))
    expect_equal(d$annotate_switch, c(2L, 5L, 3L, 2L))
    expect_equal(d$check_state, c(NA, "waiting for Assemble Update", "failed", NA))
    expect_equal(d$check_note, c(NA, NA, "tool crashed", NA))
    expect_equal(d$checks, c("ITSx", "waiting for Assemble Update", "failed", "ok"))
    expect_equal(d$checks_tone, c("warning", "info", "info", "success"))

    session$setInputs(Annotate_state_filter = "5")
    expect_equal(unit_key(samples_data()$ID, samples_data()$path, samples_data()$scaffold),
      unit_key("S1", 1L, 2L))
    j <- jsonlite::fromJSON(as.character(output$tbl), simplifyVector = FALSE)
    cols <- j$x$tag$attribs$columns
    names(cols) <- vapply(cols, function(c) c$id %||% "", "")
    expect_false(cols$check_state$show)
    expect_false(cols$check_note$show)
    for (k in c("n_flags", "itsx_warnings", "genbank_check", "screen_verdict", "chimera_label", "feature_blast_warnings", "blast_species",
                "blast_accession", "blast_pident", "blast_qcovs", "blast_source")) {
      expect_match(cols[[k]]$cell, "row['check_state']", fixed = TRUE, info = k)
      expect_true(isTRUE(cols[[k]]$html), info = k)
    }
    expect_match(cols$annotate_switch$cell, "(cellInfo.row || {})['check_note']", fixed = TRUE)
    expect_match(cols$annotate_switch$cell, "Waiting for Assemble Update", fixed = TRUE)
  })
})

test_that("the Feature BLAST panel lists one row per feature once checks are done", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  on.exit(DBI::dbDisconnect(con))
  fixture_mark_done(con)
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  DBI::dbExecute(con, "INSERT INTO review_feature_blast (ID, path, scaffold, element, length, status, low_conf, species)
    VALUES ('S1',1,1,'28S',80,'ok',1,'Danio rerio'), ('S1',1,1,'18S',1800,'no_hit',0,NULL)")
  DBI::dbExecute(con, "INSERT INTO review_feature_blast_hits (ID, path, scaffold, element, hit_rank, accession, title)
    VALUES ('S1',1,1,'28S',1,'AB000001','Danio rerio 28S')")
  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    j <- jsonlite::fromJSON(as.character(output$feature_blast_tbl), simplifyVector = FALSE)
    expect_equal(unlist(j$x$tag$attribs$data$element), c("18S", "28S"))
    expect_equal(unlist(j$x$tag$attribs$data$status), c("No hit", "OK"))
    expect_match(j$x$tag$attribs$data$note[[2]], "Low confidence")
  })
})

test_that("the GenBank BLAST column has help text and renders like the ITSx one", {
  expect_match(RiboPilot:::SCAFFOLD_COL_HELP$genbank_check, "GenBank", fixed = TRUE)
  expect_true(grepl("'not checked'", RiboPilot:::SCAFFOLD_COL_HELP$genbank_check, fixed = TRUE))
  srv <- paste(deparse(RiboPilot:::app_server), collapse = " ")
  expect_true(grepl('genbank_check = cd("genbank_check"', srv, fixed = TRUE))
  expect_true(grepl('name = "GenBank BLAST"', srv, fixed = TRUE))
})
