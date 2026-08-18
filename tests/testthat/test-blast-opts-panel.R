# The BLAST options modal used to say nothing about the panel it searches, offer
# no way to turn the remote search on or off, and present the Entrez query as if
# it applied to both searches. These cover the panel note, the two new toggles,
# and the remote-only labelling.

mk_manifest <- function(dir, version = "2026.08.18", name = "ribopilot_rrna") {
  jsonlite::write_json(list(
    manifest_schema = 1L, db_name = name, db_version = version,
    built_utc = "2026-08-18T15:27:33Z",
    counts = list(records_after_dedup = 598202L, distinct_taxids = 185708L,
      total_bases = 688824954)
  ), file.path(dir, "manifest.json"), auto_unbox = TRUE)
  dir
}

test_that("panel_summary states version, build date, records and taxa", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  m <- RiboPilot:::read_panel_manifest(mk_manifest(d))
  s <- RiboPilot:::panel_summary(m)
  expect_match(s, "2026\\.08\\.18")
  expect_match(s, "2026-08-18")           # build date, not the timestamp
  expect_match(s, "598,202")              # thousands separators, not raw digits
  expect_match(s, "185,708")
  expect_false(grepl("15:27:33", s))
})

test_that("panel_summary survives a manifest with no counts block", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  jsonlite::write_json(list(db_name = "ribopilot_rrna", db_version = "2026.08.18"),
    file.path(d, "manifest.json"), auto_unbox = TRUE)
  s <- RiboPilot:::panel_summary(RiboPilot:::read_panel_manifest(d))
  expect_match(s, "2026\\.08\\.18")
  expect_false(grepl("NA", s))
})

test_that("panel_note_text prefers the panel on disk, then what the last run recorded", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  m <- RiboPilot:::read_panel_manifest(mk_manifest(d))

  # panel visible to the app (conda / local run): full detail
  expect_match(RiboPilot:::panel_note_text(m, "recorded text"), "598,202")

  # panel only inside the container: fall back to what the run wrote down
  expect_match(RiboPilot:::panel_note_text(NULL, "ribopilot_rrna 2026.08.18, 4 records"),
    "4 records")

  # nothing yet: say so plainly rather than implying there is no panel
  txt <- RiboPilot:::panel_note_text(NULL, "")
  expect_match(txt, "first time a local search runs")
  expect_false(grepl("NA", txt))
})

test_that("record_blastdb_provenance stores the summary the modal displays", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  mk_manifest(d)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")

  RiboPilot:::record_blastdb_provenance(con, file.path(d, "ribopilot_rrna"))
  got <- DBI::dbGetQuery(con,
    "SELECT blastdb_version, blastdb_info FROM reference_config WHERE rowid_key = 1")
  expect_equal(got$blastdb_version, "ribopilot_rrna 2026.08.18")
  expect_match(got$blastdb_info, "598,202")
})

test_that("blast_opts carries the remote toggles, remote off and fallback on", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = "S1", Taxon = "Danio rerio",
    R1 = "S1_R1.fq.gz", R2 = "S1_R2.fq.gz"), map, row.names = FALSE)
  db <- file.path(d, ".sqlite")
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  row <- DBI::dbGetQuery(con, "SELECT * FROM blast_opts WHERE blast_opts = 'default'")
  expect_equal(row$remote_blast, 0L)
  expect_equal(row$remote_fallback, 1L)
})

test_that("opts_row writes both remote toggles as integers", {
  input <- list(bo_run = TRUE, bo_entrez = "", bo_extra = "", bo_mts = 50,
    bo_remote = FALSE, bo_fallback = TRUE)
  row <- RiboPilot:::opts_row("bo", "default", input, base = list())
  expect_equal(row$remote_blast, 0L)
  expect_equal(row$remote_fallback, 1L)

  input$bo_remote <- TRUE
  expect_equal(RiboPilot:::opts_row("bo", "default", input, base = list())$remote_blast, 1L)
})

test_that("the modal shows the panel note, both toggles, and marks Entrez remote-only", {
  html <- as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
    RiboPilot:::blast_opts_modal(
    cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
      extra_opts = "-evalue 1e-20", max_target_seqs = 50L,
      remote_blast = 0L, remote_fallback = 1L,
      panel_note = "Bundled BLAST panel ribopilot_rrna 2026.08.18, 598,202 sequences."),
    choices = "default", n = 1)))
  expect_match(html, "598,202")
  expect_match(html, "bo_remote")
  expect_match(html, "bo_fallback")
  expect_match(html, "remote BLAST only", fixed = TRUE)
  # the old title claimed every search went to GenBank
  expect_false(grepl("GenBank BLAST options", html, fixed = TRUE))
})

test_that("blast.nf honours both toggles", {
  p <- system.file("nextflow", "modules", "blast.nf", package = "RiboPilot")
  if (!nzchar(p)) p <- "inst/nextflow/modules/blast.nf"
  skip_if_not(file.exists(p))
  txt <- paste(readLines(p, warn = FALSE), collapse = "\n")

  expect_match(txt, "remote_blast", info = "blast.nf ignores the Remote BLAST toggle")
  expect_match(txt, "remote_fallback", info = "blast.nf ignores the remote fallback toggle")
  # the escalation itself: a local unit with no usable hit retries remotely
  expect_match(txt, "RF\\b", info = "no remote-fallback variable in the script")
})
