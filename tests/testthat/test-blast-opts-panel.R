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
  expect_match(s, "2026\\.08\\.18")     # the version is the harvest date
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

test_that("panel_note prefers the database on disk, then what the last run recorded", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  m <- RiboPilot:::read_panel_manifest(mk_manifest(d))

  # visible to the app (conda / local run): full detail
  expect_match(RiboPilot:::panel_note(m, "recorded text")$text, "598,202")

  # only inside the container: fall back to what the run wrote down
  expect_match(RiboPilot:::panel_note(NULL, "ribopilot_rrna 2026.08.18, 4 records")$text,
    "4 records")

  # nothing yet: say so plainly rather than implying there is no database
  note <- RiboPilot:::panel_note(NULL, "")
  expect_false(note$known)
  expect_match(note$text, "first local search")
  expect_false(grepl("NA", note$text))
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


test_that("the modal shows the panel note, both toggles, and marks Entrez remote-only", {
  html <- as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
    RiboPilot:::blast_opts_modal(
    cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
      extra_opts = "-evalue 1e-20", max_target_seqs = 50L,
      remote_blast = 0L, remote_fallback = 1L,
      panel_note = "ribopilot_rrna 2026.08.18, 598,202 sequences", panel_known = TRUE),
    choices = "default", n = 1)))
  expect_match(html, "598,202")
  expect_match(html, "bo_source")
  expect_match(html, "Entrez query", fixed = TRUE)
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

test_that("panel_summary reports the sequence volume as well as the counts", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  s <- RiboPilot:::panel_summary(RiboPilot:::read_panel_manifest(mk_manifest(d)))
  expect_match(s, "689 Mbp")
})

test_that("config_container_image reads the image out of a project config", {
  lines <- c("docker {", "  enabled = true", "}", "process {",
    "    container = 'ribopilot:0.0.0.9011'", "}")
  expect_equal(RiboPilot:::config_container_image(lines), "ribopilot:0.0.0.9011")
  expect_equal(RiboPilot:::config_container_image(c("process {", "}")), "")
})

test_that("container_version_gap flags an official image for another version only", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  cfg <- function(img) writeLines(c("process {", sprintf("    container = '%s'", img), "}"),
    file.path(d, ".config"))
  expect_null(RiboPilot:::container_version_gap(d))
  now <- paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot"))
  cfg(now)
  expect_null(RiboPilot:::container_version_gap(d))
  cfg("ribopilot:batch4")
  expect_null(RiboPilot:::container_version_gap(d))
  cfg("macguigand/ribopilot:0.0.0.1")
  g <- RiboPilot:::container_version_gap(d)
  expect_equal(g$configured, "macguigand/ribopilot:0.0.0.1")
  expect_equal(g$expected, now)
})

test_that("the note falls back to asking the image, and caches what it learns", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  writeLines(c("docker {", "  enabled = true", "}", "process {",
    "  container = 'ribopilot:test'", "}"), file.path(d, ".config"))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")

  calls <- 0
  probe <- function(...) { calls <<- calls + 1; "ribopilot_rrna 2026.08.18, 42 sequences" }
  note <- RiboPilot:::project_panel_note(con, "", project = d, probe = probe)
  expect_match(note$text, "42 sequences")
  expect_equal(calls, 1)

  # cached in the project, so opening the modal again does not re-run the image
  note2 <- RiboPilot:::project_panel_note(con, "", project = d, probe = probe)
  expect_match(note2$text, "42 sequences")
  expect_equal(calls, 1)
})

test_that("a probe that cannot reach the image leaves the placeholder, not an error", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(d, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key) VALUES (1)")
  note <- RiboPilot:::project_panel_note(con, "", project = d,
    probe = function(...) stop("no docker here"))
  expect_false(note$known)
})

test_that("the user-facing wording says reference database, never panel", {
  html <- as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
    RiboPilot:::blast_opts_modal(
      cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
        extra_opts = "", max_target_seqs = 50L, remote_blast = 0L,
        remote_fallback = 1L, panel_note = RiboPilot:::panel_note(NULL, "")$text,
        panel_known = FALSE),
      choices = "default", n = 1)))
  # help_panel_* are popover element ids from tool_help_icon(), not wording
  visible <- gsub("help_panel_[a-z]+", "", html)
  expect_false(grepl("\\bpanel\\b", visible, ignore.case = TRUE))
  expect_match(html, "Bundled database")
  expect_false(grepl("panel", RiboPilot:::panel_note(NULL, "")$text, ignore.case = TRUE))
})

test_that("the modal's note follows the project, not the session's working dir", {
  skip_if_not_installed("withr")
  proj <- withr::local_tempdir()
  dbdir <- file.path(proj, "db"); dir.create(dbdir)
  mk_manifest(dbdir)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(proj, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, RiboPilot:::REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con, "INSERT INTO reference_config (rowid_key, blastdb) VALUES (1, ?)",
    params = list(file.path(dbdir, "ribopilot_rrna")))

  # the app is normally launched from somewhere else entirely
  elsewhere <- withr::local_tempdir()
  withr::local_dir(elsewhere)
  cur <- RiboPilot:::with_panel_note(
    data.frame(blast_opts = "default", stringsAsFactors = FALSE), con, project = proj)
  expect_match(cur$panel_note, "598,202")
})

test_that("the summary does not print the build date twice", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  # version IS the harvest date, so a build on the same day adds nothing
  s <- RiboPilot:::panel_summary(RiboPilot:::read_panel_manifest(mk_manifest(d)))
  expect_match(s, "2026\\.08\\.18")
  expect_false(grepl("built", s))

  # a build that lagged the harvest is worth stating
  jsonlite::write_json(list(db_name = "ribopilot_rrna", db_version = "2026.08.18",
    built_utc = "2026-09-02T10:00:00Z",
    counts = list(records_after_dedup = 10L)),
    file.path(d, "manifest.json"), auto_unbox = TRUE)
  s2 <- RiboPilot:::panel_summary(RiboPilot:::read_panel_manifest(d))
  expect_match(s2, "built 2026-09-02")
})

test_that("panel_note reports whether the database is known, not just its text", {
  skip_if_not_installed("withr")
  d <- withr::local_tempdir()
  known <- RiboPilot:::panel_note(RiboPilot:::read_panel_manifest(mk_manifest(d)), "")
  expect_true(known$known)
  expect_match(known$text, "598,202")
  # the lead-in belongs to the callout, not the sentence
  expect_false(grepl("Bundled", known$text))

  unknown <- RiboPilot:::panel_note(NULL, "")
  expect_false(unknown$known)
  expect_match(unknown$text, "first local search")
})

test_that("the note renders as a callout, not as a bordered field", {
  box <- as.character(RiboPilot:::panel_note_box(
    list(known = TRUE, text = "ribopilot_rrna 2026.08.18, 42 sequences")))
  # a full 1px border plus a radius is what made it read as a disabled input
  expect_false(grepl("border:1px solid", box, fixed = TRUE))
  expect_match(box, "border-left")
  expect_match(box, "role=\"note\"", fixed = TRUE)
  expect_match(box, "Bundled database")          # lead-in, inside the block
  expect_match(box, "42 sequences")

  # the unknown state differs by wording, not only by colour
  un <- as.character(RiboPilot:::panel_note_box(list(known = FALSE, text = "x")))
  expect_match(un, "Not yet recorded")
})

test_that("Entrez shows whenever a search can reach NCBI, and says so", {
  mk <- function(remote, fallback) {
    as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
      RiboPilot:::blast_opts_modal(
        cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
          extra_opts = "", max_target_seqs = 50L, remote_blast = remote,
          remote_fallback = fallback, panel_note = "x", panel_known = TRUE),
        choices = "default", n = 1)))
  }
  hidden <- "display:none"
  # fallback alone still sends unnamed scaffolds to NCBI with this filter
  expect_false(grepl(hidden, sub(".*bo_entrez_group", "", mk(0L, 1L))))
  expect_false(grepl(hidden, sub(".*bo_entrez_group", "", mk(1L, 0L))))
  # neither: no search reaches NCBI, so the filter is irrelevant
  expect_match(sub("(?s).*id=\"bo_entrez_group\"", "", mk(0L, 0L), perl = TRUE),
    hidden, fixed = TRUE)
  expect_match(mk(0L, 1L), "any search that goes to NCBI")
  expect_false(grepl("remote BLAST only", mk(0L, 1L), fixed = TRUE))
})


test_that("the modal states the metazoan rDNA scope once, not three times", {
  html <- as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
    RiboPilot:::blast_opts_modal(
      cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
        extra_opts = "", max_target_seqs = 50L, remote_blast = 0L,
        remote_fallback = 1L, panel_note = "x", panel_known = TRUE),
      choices = "default", n = 1)))
  expect_lte(length(gregexpr("metazoan", html)[[1]]), 1)
})

test_that("the source choice round-trips through the two stored columns", {
  f <- RiboPilot:::blast_source
  expect_equal(f(0L, 0L), "bundled")
  expect_equal(f(0L, 1L), "escalate")
  expect_equal(f(1L, 0L), "remote")
  # legacy and hand-edited rows: remote wins, because that is what blast.nf does
  expect_equal(f(1L, 1L), "remote")
  # missing columns fall to the shipped default
  expect_equal(f(NA, NA), "escalate")

  g <- RiboPilot:::blast_source_cols
  expect_equal(g("bundled"), list(remote_blast = 0L, remote_fallback = 0L))
  expect_equal(g("escalate"), list(remote_blast = 0L, remote_fallback = 1L))
  # the impossible state is not writable: choosing NCBI only clears the flag
  expect_equal(g("remote"), list(remote_blast = 1L, remote_fallback = 0L))
})

test_that("opts_row writes the columns the radio implies", {
  base <- list(bo_run = TRUE, bo_entrez = "", bo_extra = "", bo_mts = 50)
  row <- RiboPilot:::opts_row("bo", "default", c(base, list(bo_source = "remote")),
    base = list())
  expect_equal(row$remote_blast, 1L)
  expect_equal(row$remote_fallback, 0L)

  row2 <- RiboPilot:::opts_row("bo", "default", c(base, list(bo_source = "bundled")),
    base = list())
  expect_equal(row2$remote_blast, 0L)
  expect_equal(row2$remote_fallback, 0L)
})

test_that("the modal offers one three-way choice, not two checkboxes", {
  mk <- function(remote, fallback) {
    as.character(shiny::withReactiveDomain(shiny::MockShinySession$new(),
      RiboPilot:::blast_opts_modal(
        cur = list(blast_opts = "default", run_blast = 1L, entrez_query = "",
          extra_opts = "", max_target_seqs = 50L, remote_blast = remote,
          remote_fallback = fallback, panel_note = "x", panel_known = TRUE),
        choices = "default", n = 1)))
  }
  html <- mk(0L, 1L)
  expect_match(html, "bo_source")
  expect_false(grepl("bo_remote", html, fixed = TRUE))
  expect_false(grepl("bo_fallback", html, fixed = TRUE))
  # the contradiction the two checkboxes allowed has nowhere to live now
  expect_false(grepl("Not used when every search", html, fixed = TRUE))

  # the stored state selects the matching option
  expect_match(mk(1L, 1L), 'value="remote"[^>]*checked')
  expect_match(mk(0L, 0L), 'value="bundled"[^>]*checked')

  # Entrez is irrelevant only when nothing reaches NCBI
  expect_match(sub("(?s).*id=\"bo_entrez_group\"", "", mk(0L, 0L), perl = TRUE),
    "display:none", fixed = TRUE)
  expect_false(grepl("display:none", sub(".*bo_entrez_group", "", mk(0L, 1L))))
})
