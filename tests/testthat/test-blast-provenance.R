# A scaffold with no GenBank screen used to get one message whether BLAST was
# switched off or switched on and never returned, and the seed-file hit that had
# already named a species was not shown at all.

mk_blast_opts_project <- function(run_blast, env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_opts(blast_opts TEXT, run_blast INTEGER)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO blast_opts VALUES('default',?)", params = list(run_blast))
  dir
}

test_that("blast_enabled reads run_blast, and defaults to on when unknown", {
  skip_if_not_installed("withr")
  withr::local_dir(mk_blast_opts_project(1L))
  expect_true(RiboPilot:::blast_enabled("S1"))
  expect_true(RiboPilot:::blast_enabled("no_such_sample"))
  expect_true(RiboPilot:::blast_enabled(NULL))

  withr::local_dir(mk_blast_opts_project(0L))
  expect_false(RiboPilot:::blast_enabled("S1"))
})

test_that("the stored hit is labelled by the search that produced it", {
  r <- data.frame(blast_species = "Danio rerio", blast_accession = "MW829490.1",
    blast_pident = 99.5, blast_source = "seeds", stringsAsFactors = FALSE)
  html <- as.character(RiboPilot:::seed_hit_block(r))
  expect_match(html, "Danio rerio")
  expect_match(html, "MW829490.1")
  expect_match(html, "NOT a GenBank identification")

  # a legacy row has no source recorded; it is a seed hit, so it keeps the caveat
  r$blast_source <- ""
  expect_match(as.character(RiboPilot:::seed_hit_block(r)), "NOT a GenBank identification")

  r$blast_source <- "core_nt"
  html <- as.character(RiboPilot:::seed_hit_block(r))
  expect_match(html, "GenBank search")
  expect_false(grepl("seed file", html, fixed = TRUE))
})

test_that("no stored hit renders nothing", {
  expect_null(RiboPilot:::seed_hit_block(
    data.frame(blast_species = NA_character_, blast_accession = NA_character_,
      blast_pident = NA_real_, blast_source = NA_character_)))
  expect_null(RiboPilot:::seed_hit_block(
    data.frame(blast_species = "", blast_accession = "", blast_pident = NA_real_,
      blast_source = "")))
})
