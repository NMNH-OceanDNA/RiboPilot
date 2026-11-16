# Review copies and the views over them (SPEC S4.2, S4.3, S11).

RS_SEQ <- strrep("ACGTACGTAC", 10)

# Four 100 bp contigs of S1 in the original tables. Scaffold 1 is done without
# a stored ITSx value, 2 is done, 3 is running in the main pipeline with no
# stored ITSx value. Scaffold 5 is a dropped contig with no annotate row; the
# review layer does not list it.
rs_originals <- function(con) {
  DBI::dbExecute(con, "DELETE FROM annotate")
  for (s in 1:4) {
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence,
        blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
      VALUES ('S1', 1, ?, 'linear', 100, ?, 'A1', 'Danio rerio', 99.5, 88, 1e-50, 'core_nt', 999)",
      params = list(s, RS_SEQ))
    DBI::dbExecute(con, "INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock,
        n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, time_stamp)
      VALUES ('S1', 1, ?, ?, 0, 2, 'ITS1', 0, 3, 'length', ?, 111)",
      params = list(s, c(2L, 2L, 4L, 2L)[s], c(NA, "", NA, "")[s]))
    for (f in list(c("18S", 1, 50), c("28S", 60, 100))) {
      DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
        VALUES ('S1', 1, ?, ?, ?, ?, '+')", params = list(s, f[1], as.integer(f[2]), as.integer(f[3])))
    }
    DBI::dbExecute(con, "INSERT INTO screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, ?, 'on_target')",
      params = list(s))
    DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, label) VALUES ('S1', 1, ?, 'clean')",
      params = list(s))
    for (w in c(1L, 51L)) {
      DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, ?, ?, ?)",
        params = list(s, w, w + 49L))
    }
  }
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore)
    VALUES ('S1', 1, 5, 'linear', 100, ?, 1)", params = list(RS_SEQ))
}

# The stored ITSx origin of every contig; a copy of scaffold 4 (trimmed to
# 11..90, one feature, no chimera row); and rows a discarded recheck left for
# scaffold 1, which has no copy.
rs_copy <- function(con) {
  DBI::dbExecute(con, "UPDATE annotate SET itsx_origin = 'M', itsx_scores = 'M:50'")
  DBI::dbExecute(con, "INSERT INTO review_units (ID, path, scaffold, sequence, length, trim_start,
      trim_end, edited, check_switch, check_note, rev, checked_at, n_features, missing, chimeric,
      n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores, blast_accession, blast_species,
      blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp)
    VALUES ('S1', 1, 4, ?, 80, 11, 90, 1, 3, 'boom', CAST(? AS INTEGER), 333, 1, '28S', 1, 0, '',
      'itsx_x', 'F', 'F:10', 'B2', 'Copy sp', 90, 70, 1e-10, 'local', 222)",
    params = list(substr(RS_SEQ, 11, 90), "1760000000123"))
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, start, stop, strand,
      partial_start_user) VALUES ('S1', 1, 4, '18S', 1, 40, '+', 1)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 4, 'uncertain')")
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, 4, 1, 50)")
  DBI::dbExecute(con, "INSERT INTO review_annotations (ID, path, scaffold, element, start, stop)
    VALUES ('S1', 1, 1, 'ITS1', 51, 59)")
  DBI::dbExecute(con, "INSERT INTO review_screen (ID, path, scaffold, screen_verdict) VALUES ('S1', 1, 1, 'off_target')")
  DBI::dbExecute(con, "INSERT INTO review_chimera (ID, path, scaffold, label) VALUES ('S1', 1, 1, 'chimeric')")
  DBI::dbExecute(con, "INSERT INTO review_chimera_windows (ID, path, scaffold, start, stop) VALUES ('S1', 1, 1, 7, 8)")
}

rs_check_views <- function(con) {
  u <- DBI::dbGetQuery(con, "SELECT scaffold, has_copy, topology, orig_length, sequence, length,
      trim_start, trim_end, edited, check_switch, check_note, CAST(rev AS TEXT) AS rev, checked_at,
      n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores,
      blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source, time_stamp
    FROM review_units_v ORDER BY scaffold")
  expect_equal(u$scaffold, 1:4)
  expect_equal(u$has_copy, c(0L, 0L, 0L, 1L))
  expect_equal(u$topology, rep("linear", 4))
  expect_equal(u$orig_length, rep(100L, 4))
  expect_equal(u$sequence, c(rep(RS_SEQ, 3), substr(RS_SEQ, 11, 90)))
  expect_equal(u$length, c(100L, 100L, 100L, 80L))
  expect_equal(u$trim_start, c(1L, 1L, 1L, 11L))
  expect_equal(u$trim_end, c(100L, 100L, 100L, 90L))
  expect_equal(u$edited, c(0L, 0L, 0L, 1L))
  # done without a stored ITSx value reads as not run
  expect_equal(u$check_switch, c(1L, 2L, 4L, 3L))
  expect_equal(u$check_note, c(NA, NA, NA, "boom"))
  expect_equal(u$rev, c(NA, NA, NA, "1760000000123"))
  expect_equal(u$checked_at, c(NA, NA, NA, 333L))
  expect_equal(u$n_features, c(2L, 2L, 2L, 1L))
  expect_equal(u$missing, c("ITS1", "ITS1", "ITS1", "28S"))
  expect_equal(u$chimeric, c(0L, 0L, 0L, 1L))
  expect_equal(u$n_flags, c(3L, 3L, 3L, 0L))
  expect_equal(u$curate_flags, c("length", "length", "length", ""))
  expect_equal(u$itsx_warnings, c(NA, "", NA, "itsx_x"))
  expect_equal(u$itsx_origin, c("M", "M", "M", "F"))
  expect_equal(u$itsx_scores, c("M:50", "M:50", "M:50", "F:10"))
  expect_equal(u$blast_accession, c("A1", "A1", "A1", "B2"))
  expect_equal(u$blast_species, c("Danio rerio", "Danio rerio", "Danio rerio", "Copy sp"))
  expect_equal(u$blast_pident, c(99.5, 99.5, 99.5, 90))
  expect_equal(u$blast_qcovs, c(88, 88, 88, 70))
  expect_equal(u$blast_evalue, c(1e-50, 1e-50, 1e-50, 1e-10))
  expect_equal(u$blast_source, c("core_nt", "core_nt", "core_nt", "local"))
  expect_equal(u$time_stamp, c(111L, 111L, 111L, 222L))

  a <- DBI::dbGetQuery(con, "SELECT scaffold, element, start, stop, partial_start_user
    FROM review_annotations_v ORDER BY scaffold, start")
  expect_equal(a$scaffold, c(1L, 1L, 2L, 2L, 3L, 3L, 4L))
  expect_equal(a$element, c(rep(c("18S", "28S"), 3), "18S"))
  expect_equal(a$stop, c(rep(c(50L, 100L), 3), 40L))
  expect_equal(a$partial_start_user, c(rep(NA, 6), 1L))

  s <- DBI::dbGetQuery(con, "SELECT scaffold, screen_verdict FROM review_screen_v ORDER BY scaffold")
  expect_equal(s$scaffold, 1:4)
  expect_equal(s$screen_verdict, c("on_target", "on_target", "on_target", "uncertain"))
  ch <- DBI::dbGetQuery(con, "SELECT scaffold, label FROM review_chimera_v ORDER BY scaffold")
  expect_equal(ch$scaffold, 1:3)
  expect_equal(ch$label, rep("clean", 3))
  w <- DBI::dbGetQuery(con, "SELECT scaffold, start FROM review_chimera_windows_v ORDER BY scaffold, start")
  expect_equal(w$scaffold, c(1L, 1L, 2L, 2L, 3L, 3L, 4L))
  expect_equal(w$start, c(rep(c(1L, 51L), 3), 1L))
}

test_that("new_db() makes each review table like its original, and annotate gains the ITSx columns", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  info <- function(t) DBI::dbGetQuery(con, sprintf("PRAGMA table_info(%s)", t))[,
    c("name", "type", "notnull", "dflt_value", "pk")]
  m <- RiboPilot:::REVIEW_MIRRORS
  expect_equal(names(m), c("annotations", "blast_hits", "blast_screen", "screen", "chimera",
    "chimera_windows", "blast_ref_alignment", "feature_blast", "feature_blast_hits"))
  expect_equal(RiboPilot:::REVIEW_TABLES, c("review_units", unname(m)))
  for (orig in names(m)[-1]) expect_equal(info(m[[orig]]), info(orig), info = orig)
  ra <- info("review_annotations")
  expect_equal(nrow(ra), 20L)
  expect_equal(ra[1:18, ], info("annotations"))
  expect_equal(ra$name[19:20], c("partial_start_user", "partial_stop_user"))
  expect_equal(ra$type[19:20], c("INTEGER", "INTEGER"))

  ru <- info("review_units")
  expect_equal(ru$name, c("ID", "path", "scaffold", "sequence", "length", "trim_start", "trim_end",
    "edited", "check_switch", "check_note", "rev", "checked_at", "n_features", "missing", "chimeric",
    "n_flags", "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession",
    "blast_species", "blast_pident", "blast_qcovs", "blast_evalue", "blast_source", "time_stamp",
    "genbank_check"))
  expect_equal(ru$pk, c(1:3, rep(0L, 25)))
  expect_equal(ru$notnull[ru$name %in% c("ID", "path", "scaffold", "edited", "rev")], rep(1L, 5))
  expect_equal(ru$dflt_value[ru$name == "edited"], "0")

  expect_equal(utils::tail(DBI::dbListFields(con, "annotate"), 4),
    c("export_time_stamp", "itsx_origin", "itsx_scores", "genbank_check"))

  views <- DBI::dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type = 'view'")$name
  expect_setequal(views, c(names(RiboPilot:::REVIEW_VIEWS_DDL), names(RiboPilot:::FEATURE_BLAST_VIEWS_DDL)))
  expect_equal(DBI::dbListFields(con, "review_units_v"), c("ID", "path", "scaffold", "has_copy",
    "topology", "orig_length", "sequence", "length", "trim_start", "trim_end", "edited",
    "check_switch", "check_note", "rev", "checked_at", "n_features", "missing", "chimeric", "n_flags",
    "curate_flags", "itsx_warnings", "itsx_origin", "itsx_scores", "blast_accession", "blast_species",
    "blast_pident", "blast_qcovs", "blast_evalue", "blast_source", "time_stamp", "genbank_check"))
  expect_equal(DBI::dbListFields(con, "review_annotations_v"), ra$name)
  for (t in c("screen", "chimera", "chimera_windows")) {
    expect_equal(DBI::dbListFields(con, paste0("review_", t, "_v")), DBI::dbListFields(con, t), info = t)
  }
})

test_that("review_apply_schema() is idempotent", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  before <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name")
  RiboPilot:::review_apply_schema(con)
  RiboPilot:::review_apply_schema(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name"), before)
})

test_that("the views read a contig's copy when it has one and its originals otherwise", {
  con <- chimera_test_con()
  on.exit(DBI::dbDisconnect(con))
  rs_originals(con)
  rs_copy(con)
  rs_check_views(con)

  # a copy whose features were all trimmed away shows none, not the original ones
  n_feat <- function() DBI::dbGetQuery(con,
    "SELECT COUNT(*) AS n FROM review_annotations_v WHERE scaffold = 4")$n
  DBI::dbExecute(con, "DELETE FROM review_annotations WHERE scaffold = 4")
  expect_equal(n_feat(), 0)

  # without its review_units row the contig reads its originals again
  DBI::dbExecute(con, "DELETE FROM review_units WHERE scaffold = 4")
  expect_equal(n_feat(), 2)
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_verdict FROM review_screen_v WHERE scaffold = 4")$screen_verdict,
    "on_target")
  expect_equal(DBI::dbGetQuery(con, "SELECT has_copy, length, check_switch FROM review_units_v WHERE scaffold = 4"),
    data.frame(has_copy = 0L, length = 100L, check_switch = 2L))
})

test_that("a version 1 database holding data upgrades through step 4 and its views work", {
  dir <- withr::local_tempdir()
  db <- file.path(dir, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES ('created', 1, '0.1.0', 0, NULL)")
  rs_originals(con)
  DBI::dbDisconnect(con)

  suppressMessages(backwards_compatibility(dir))
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  h <- DBI::dbGetQuery(con, "SELECT schema_version, note FROM schema_history ORDER BY rowid")
  expect_equal(h$schema_version[h$note %in% "review copies for edited contigs; ITSx origin and scores"], 5L)
  expect_equal(utils::tail(DBI::dbListFields(con, "annotate"), 4),
    c("export_time_stamp", "itsx_origin", "itsx_scores", "genbank_check"))
  rs_copy(con)
  rs_check_views(con)
})

test_that("the app fixture has every column the review layer reads, and its views answer", {
  dir <- mk_app_project()
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  ref <- chimera_test_con()
  for (t in c("assemblies", "annotate", "annotations", "coverage", "screen", "chimera",
              "chimera_windows", "blast_hits", "blast_screen", "blast_ref_alignment")) {
    expect_true(all(DBI::dbListFields(ref, t) %in% DBI::dbListFields(con, t)), info = t)
  }
  DBI::dbDisconnect(ref)
  for (v in names(RiboPilot:::REVIEW_VIEWS_DDL)) {
    expect_no_error(DBI::dbGetQuery(con, sprintf("SELECT ID, path, scaffold FROM %s LIMIT 1", v)))
  }
  u <- DBI::dbGetQuery(con, "SELECT has_copy, trim_end, check_switch FROM review_units_v ORDER BY path, scaffold")
  expect_equal(u$has_copy, rep(0L, 4))
  expect_equal(u$trim_end, c(6120L, 1442L, 5980L, 890L))
  expect_equal(u$check_switch, rep(1L, 4))

  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = 'itsx_x' WHERE path = 2 AND scaffold = 2")
  fixture_mark_done(con)
  d <- DBI::dbGetQuery(con, "SELECT v.check_switch, an.itsx_warnings FROM review_units_v v
    JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold
    ORDER BY v.path, v.scaffold")
  expect_equal(d$check_switch, rep(2L, 4))
  expect_equal(d$itsx_warnings, c("", "", "", "itsx_x"))

  before <- DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name")
  fixture_review_schema(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT type, name, sql FROM sqlite_master ORDER BY name"), before)
})

test_that("genbank_check follows itsx_warnings through copy, view and invalidate", {
  con <- chimera_test_con(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT OR REPLACE INTO assemblies (ID, path, scaffold, sequence, length) VALUES ('S1',1,1,'ACGT',4)")
  DBI::dbExecute(con, "INSERT OR REPLACE INTO annotate (ID, path, scaffold, genbank_check) VALUES ('S1',1,1,'genbank_low_coverage(ITS 70%)')")
  expect_equal(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units_v")$genbank_check, "genbank_low_coverage(ITS 70%)")
  RiboPilot:::review_copy(con, "S1", 1L, 1L)
  expect_equal(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units")$genbank_check, "genbank_low_coverage(ITS 70%)")
  RiboPilot:::review_invalidate(con, "S1", 1L, 1L)
  expect_true(is.na(DBI::dbGetQuery(con, "SELECT genbank_check FROM review_units_v")$genbank_check))
  expect_equal(RiboPilot:::current_schema_version(), 7L)
})
