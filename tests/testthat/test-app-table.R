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
