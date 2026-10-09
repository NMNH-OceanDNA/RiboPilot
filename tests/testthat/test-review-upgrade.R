# Schema step 4 turns the trims and edits an older RiboPilot wrote into the
# original tables into review copies, then drops assembly_backup (SPEC S11.3-11.5).

RU_SEQ <- withr::with_seed(7, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
RU_OTHER <- withr::with_seed(8, paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = ""))
RU_DEP <- c(rep(5, 100), rep(15, 100))

# A version 4 database: the 0.1.0 schema plus steps 1-3, with sample S1
ru_v4 <- function(dir) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  sql <- paste(readLines(test_path("fixtures", "schema_v1.sql")), collapse = "\n")
  for (s in strsplit(sql, "\n;(\n|$)")[[1]]) DBI::dbExecute(con, s)
  for (st in RiboPilot:::db_migrations()[1:3]) st$fn(con)
  RiboPilot:::record_schema_event(con, "created", 4L)
  DBI::dbExecute(con, "INSERT INTO samples (ID) VALUES ('S1')")
  con
}

# Contig S1.1.<s>: 200 bp, three + strand features, a depth track, and a
# chimera result over three windows
ru_unit <- function(con, s, lock = 0L, sw = 2L, itsx = "itsx_x") {
  x <- function(sql, ...) DBI::dbExecute(con, sql, params = list(...))
  x("INSERT INTO assemblies (ID, path, scaffold, topology, length, sequence, ignore, edited)
     VALUES ('S1', 1, ?, 'linear', 200, ?, 0, 0)", s, RU_SEQ)
  x("INSERT INTO annotate (ID, path, scaffold, annotate_switch, annotate_lock, n_features, missing,
       n_flags, curate_flags, itsx_warnings) VALUES ('S1', 1, ?, ?, ?, 99, 'x', 7, 'stale', ?)",
    s, sw, lock, itsx)
  for (f in list(list("18S", 21L, 80L), list("ITS1", 81L, 120L), list("28S", 121L, 190L))) {
    x("INSERT INTO annotations (ID, path, scaffold, element, type, start, stop, length, strand,
         model_coverage, tool, status, edited)
       VALUES ('S1', 1, ?, ?, 'rRNA', ?, ?, ?, '+', 99, 'barrnap', 'stale', 0)",
      s, f[[1]], f[[2]], f[[3]], f[[3]] - f[[2]] + 1L)
  }
  x("INSERT INTO coverage (ID, path, scaffold, mean_depth, median_depth, min_depth, max_depth, depths)
     VALUES ('S1', 1, ?, 10, 10, 5, 15, ?)", s, paste(RU_DEP, collapse = ","))
  x("INSERT INTO chimera (ID, path, scaffold, score, label, junctions)
     VALUES ('S1', 1, ?, 0.9, 'severe', '30,150')", s)
  for (w in c(1L, 51L, 101L)) {
    x("INSERT INTO chimera_windows (ID, path, scaffold, start, stop, lineage, lineage_taxids, pident)
       VALUES ('S1', 1, ?, ?, ?, 'K;P;C1;o;f;g;s', '1;2;3;;;;', 99)", s, w, w + 49L)
  }
}

# What the old trim_assembly() left after cutting S1.1.<s> to 41..170: the
# backup row (snapshots as JSON), the cut sequence, features shifted with cut
# ends partial and edited = 1, the depth track cut, and the chimera rows shifted
ru_trim <- function(con, s, snap = TRUE) {
  rows <- function(t) DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s WHERE scaffold = ?",
    paste(DBI::dbListFields(con, t), collapse = ", "), t), params = list(s))
  chim <- if (snap) {
    as.character(jsonlite::toJSON(list(chimera = rows("chimera"), windows = rows("chimera_windows")),
      dataframe = "rows", na = "null", digits = NA))
  } else NA_character_
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length, annotations,
      coverage, chimera, time_stamp) VALUES ('S1', 1, ?, ?, 200, ?, ?, ?, 0)",
    params = list(s, RU_SEQ, as.character(jsonlite::toJSON(rows("annotations"), dataframe = "rows",
      na = "null")), paste(RU_DEP, collapse = ","), chim))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = 130, edited = 1 WHERE scaffold = ?",
    params = list(substr(RU_SEQ, 41, 170), s))
  DBI::dbExecute(con, "UPDATE annotations SET start = MAX(start - 40, 1), stop = MIN(stop - 40, 130),
      edited = 1, partial_start = CASE WHEN start < 41 THEN 1 ELSE partial_start END,
      partial_stop = CASE WHEN stop > 170 THEN 1 ELSE partial_stop END
    WHERE scaffold = ?", params = list(s))
  DBI::dbExecute(con, "UPDATE annotations SET length = stop - start + 1 WHERE scaffold = ?", params = list(s))
  RiboPilot:::write_coverage_summary(con, "S1", 1L, s, RU_DEP[41:170])
  RiboPilot:::chimera_shift(con, "S1", 1L, s, 40L, 130L)
}

test_that("step 4 turns old trims and edits into review copies and lists what it could not convert", {
  dir <- withr::local_tempdir()
  con <- ru_v4(dir)
  for (s in 1:12) {
    ru_unit(con, s, lock = as.integer(s %in% c(4, 5, 12)), sw = if (s == 4) 1L else 2L,
      itsx = if (s %in% c(4, 5)) NA_character_ else if (s == 12) "" else "itsx_x")
  }
  # 1 trimmed; 2 trimmed, then restored; 3 nudged, its 5' end set by hand;
  # 4 locked, State set by hand; 5 trimmed and locked, no ITSx value, no chimera
  # snapshot; 6 its sequence is not in the backup; 7 backup without a sequence;
  # 8 snapshots holding a column the tables no longer have; 9 its features
  # cannot be written back; 10 trimmed, then lost its annotate row; 11 fails
  # part way; 12 locked and done
  for (s in c(1, 5:11)) ru_trim(con, s, snap = s != 5)
  DBI::dbExecute(con, "INSERT INTO assembly_backup (ID, path, scaffold, sequence, length)
    VALUES ('S1', 1, 2, ?, 200)", params = list(RU_SEQ))
  DBI::dbExecute(con, "UPDATE annotations SET start = 25, length = 56, edited = 1, partial_start = 1
    WHERE scaffold = 3 AND element = '18S'")
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = ? WHERE scaffold = 6", params = list(RU_OTHER))
  DBI::dbExecute(con, "UPDATE assembly_backup SET sequence = NULL WHERE scaffold = 7")
  DBI::dbExecute(con, "UPDATE assembly_backup SET annotations = replace(annotations, ?, ?),
      chimera = replace(chimera, ?, ?) WHERE scaffold = 8",
    params = list('"element":', '"retired":1,"element":', '"label":', '"retired":1,"label":'))
  DBI::dbExecute(con, "CREATE TRIGGER ru_boom BEFORE INSERT ON annotations WHEN NEW.scaffold = 9
    BEGIN SELECT RAISE(ABORT, 'ru_boom'); END")
  DBI::dbExecute(con, "DELETE FROM annotate WHERE scaffold = 10")
  DBI::dbDisconnect(con)

  real <- RiboPilot:::review_invalidate
  local_mocked_bindings(review_invalidate = function(con, id, path, scaffold) {
    if (scaffold == 11) stop("disk full")
    real(con, id, path, scaffold)
  })
  msgs <- testthat::capture_messages(backwards_compatibility(dir))
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::project_schema_version(file.path(dir, ".sqlite")), 5L)
  expect_false("assembly_backup" %in% DBI::dbListTables(con))

  # the copies hold the trimmed or edited state; unlocked ones wait for Update
  nr <- "not_run(not checked before upgrade)"
  expect_equal(DBI::dbGetQuery(con, "SELECT scaffold, length, trim_start, trim_end, edited, check_switch,
      n_features, missing, n_flags, itsx_warnings FROM review_units ORDER BY scaffold"),
    data.frame(scaffold = c(1L, 3:9), length = c(130L, 200L, 200L, rep(130L, 5)),
      trim_start = c(41L, 1L, 1L, 41L, 1L, 1L, 41L, 1L), trim_end = c(170L, 200L, 200L, 170L, 130L, 130L, 170L, 130L),
      edited = c(1L, 1L, 0L, rep(1L, 5)), check_switch = c(1L, 1L, 2L, 2L, rep(1L, 4)),
      n_features = 3L, missing = "5.8S,ITS2", n_flags = c(NA, NA, 7L, 7L, rep(NA, 4)),
      itsx_warnings = c(NA, NA, nr, nr, rep(NA, 4))))
  expect_equal(DBI::dbGetQuery(con, "SELECT sequence FROM review_units WHERE scaffold IN (1, 5, 8)")$sequence,
    rep(substr(RU_SEQ, 41, 170), 3))
  rf <- function(s) DBI::dbGetQuery(con, "SELECT element, start, stop, partial_start, partial_stop,
      partial_start_user, partial_stop_user, status FROM review_annotations WHERE scaffold = ? ORDER BY start",
    params = list(s))
  feats <- data.frame(element = c("18S", "ITS1", "28S"), start = c(1L, 41L, 81L), stop = c(40L, 80L, 130L),
    partial_start = c(1L, NA, NA), partial_stop = c(NA, NA, 1L), partial_start_user = c(1L, NA, NA),
    partial_stop_user = c(NA, NA, 1L), status = "stale")
  # locked: the copy keeps its check results
  expect_equal(rf(5), feats)
  # unlocked: the cut ends stay as the user's choice, curation's values are cleared
  feats[c("partial_start", "partial_stop", "status")] <- list(NA_integer_, NA_integer_, NA_character_)
  expect_equal(rf(1), feats)
  expect_equal(rf(3)[c("element", "start", "partial_start_user", "partial_start")],
    data.frame(element = c("18S", "ITS1", "28S"), start = c(25L, 81L, 121L),
      partial_start_user = c(1L, NA, NA), partial_start = NA_integer_))
  expect_equal(DBI::dbGetQuery(con, "SELECT scaffold, start FROM review_chimera_windows ORDER BY scaffold, start"),
    data.frame(scaffold = c(4L, 4L, 4L, 5L, 5L), start = c(1L, 51L, 101L, 11L, 61L)))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM review_annotations WHERE scaffold = 11")$n, 0L)

  # trimmed originals are put back and curated again; the others keep their state
  expect_equal(DBI::dbGetQuery(con, "SELECT length, sequence = ? AS full FROM assemblies ORDER BY scaffold",
      params = list(RU_SEQ)),
    data.frame(length = c(rep(200L, 5), 130L, 130L, 200L, 130L, 130L, 130L, 200L),
      full = c(rep(1L, 5), 0L, 0L, 1L, 0L, 0L, 0L, 1L)))
  expect_equal(DBI::dbGetQuery(con, "SELECT start FROM annotations WHERE element = '18S' ORDER BY scaffold")$start,
    c(21L, 21L, 25L, 21L, 21L, 1L, 1L, 21L, 1L, 1L, 1L, 21L))
  expect_equal(DBI::dbGetQuery(con, "SELECT element, start, stop, edited FROM annotations WHERE scaffold = 8
      ORDER BY start"),
    data.frame(element = c("18S", "ITS1", "28S"), start = c(21L, 81L, 121L), stop = c(80L, 120L, 190L),
      edited = 0L))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM annotations
      WHERE scaffold IN (1, 5, 8) AND status = 'stale'")$n, 0L)
  an <- DBI::dbGetQuery(con, "SELECT scaffold, itsx_warnings, curate_flags FROM annotate ORDER BY scaffold")
  back <- an$scaffold %in% c(1, 5, 8)
  expect_true(all(is.na(an$itsx_warnings[back])))
  expect_false(any(an$curate_flags[back] %in% "stale"))
  expect_equal(an$itsx_warnings[an$scaffold %in% c(3, 6, 7, 9, 11)], rep("itsx_x", 5))
  expect_equal(an$curate_flags[!back], rep("stale", 8))
  cv <- DBI::dbGetQuery(con, "SELECT median_depth, depths FROM coverage WHERE scaffold IN (1, 9) ORDER BY scaffold")
  expect_equal(cv$median_depth, c(10, 15))
  expect_equal(lengths(strsplit(cv$depths, ",")), c(200L, 130L))
  ch <- function(s) list(
    DBI::dbGetQuery(con, "SELECT label, junctions FROM chimera WHERE scaffold = ?", params = list(s)),
    DBI::dbGetQuery(con, "SELECT start FROM chimera_windows WHERE scaffold = ? ORDER BY start",
      params = list(s))$start)
  expect_equal(ch(1), list(data.frame(label = "severe", junctions = "30,150"), c(1L, 51L, 101L)))
  # no snapshot, or one that no longer fits: the rows are shifted back
  expect_equal(ch(5), list(data.frame(label = "severe", junctions = "150"), c(51L, 101L)))
  expect_equal(ch(8), ch(5))
  expect_equal(ch(9)[[2]], c(11L, 61L))

  v <- DBI::dbGetQuery(con, "SELECT scaffold, has_copy, orig_length, check_switch FROM review_units_v
    ORDER BY scaffold")
  expect_equal(v, data.frame(scaffold = c(1:9, 11:12), has_copy = c(1L, 0L, rep(1L, 7), 0L, 0L),
    orig_length = c(rep(200L, 5), 130L, 130L, 200L, 130L, 130L, 200L),
    check_switch = c(1L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L)))

  m <- grep("S1.1.", msgs, fixed = TRUE, value = TRUE)
  expect_length(m, 1)
  expect_match(m, "Moved the edits of 7 contig(s)", fixed = TRUE)
  expect_match(m, "S1.1.6 (original not restored: the trimmed sequence is not in the saved original)",
    fixed = TRUE)
  expect_match(m, "S1.1.7 (original not restored: the saved original has no sequence)", fixed = TRUE)
  expect_match(m, "S1.1.8 (chimera result shifted back instead of restored: ", fixed = TRUE)
  expect_match(m, "S1.1.9 (original not restored: ru_boom", fixed = TRUE)
  expect_match(m, "Not converted, left as they were: S1.1.11 (disk full).", fixed = TRUE)
  expect_match(m, "Locked contigs marked done without a recheck: S1.1.4.", fixed = TRUE)
  expect_match(m, "Trimmed contigs left as they are (dropped or never annotated): S1.1.10.", fixed = TRUE)
  expect_no_match(m, "S1\\.1\\.(1|2|3|5|12)\\b", perl = TRUE)
})
