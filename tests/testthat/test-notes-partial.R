# Sample notes on the Assemble step, unit notes in the table, and Partial = "no".

test_that("new units start with Partial set to no", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT partial FROM annotate")$partial, "no")
  DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold) VALUES ('S1', 1, 2)")
  expect_equal(DBI::dbGetQuery(con, "SELECT partial FROM annotate WHERE scaffold = 2")$partial, "no")
  expect_true("user_notes" %in% DBI::dbListFields(con, "assemble"))
})

test_that("the Assemble details notes box saves, shows in the table, and clears", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)
  note <- function() RiboPilot:::app_query("SELECT user_notes FROM assemble WHERE ID = 'S1'")$user_notes

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(assemble_user_notes = "check the 28S end")
    expect_equal(note(), "check the 28S end")
    session$setInputs(assemble_user_notes = "")
    expect_equal(note(), "check the 28S end")
    expect_true("user_notes" %in% names(samples_data()))
    session$setInputs(assemble_user_notes_clear = 1)
    expect_true(is.na(note()))
  })
})

test_that("typed notes reach the table once the details window closes", {
  skip_if_not_installed("withr")
  dir <- mk_app_project()
  withr::local_dir(dir)

  shiny::testServer(RiboPilot:::app_server, {
    session$setInputs(mode = "Assemble")
    sel_id("S1")
    session$flushReact()
    session$setInputs(assemble_user_notes = "check the 28S end")
    session$setInputs(modal_hidden = 1)
    expect_equal(samples_data()$user_notes[samples_data()$ID == "S1"], "check the 28S end")

    session$setInputs(mode = "Annotate")
    sel_unit(list(ID = "S1", path = 1L, scaffold = 1L))
    session$flushReact()
    session$setInputs(annotate_notes = "recheck ITS2")
    session$setInputs(modal_hidden = 2)
    d <- samples_data()
    expect_equal(d$annotate_notes[d$ID == "S1" & d$path == 1 & d$scaffold == 1], "recheck ITS2")
  })
})
