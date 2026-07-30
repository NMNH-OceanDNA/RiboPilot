#!/usr/bin/env Rscript
# Generator for the screen validation fixture (dev/design/13_VALIDATION_FIXTURES.md).
#
# DEV ONLY. Never run by tests and never shipped: dev/ is in .Rbuildignore, so
# this file is absent from the built package. Tests read the committed fixture
# and never reach the network.
#
# Usage:
#   Rscript dev/make_screen_fixture.R host <project_db> [taxdump_date] [--force]
#   Rscript dev/make_screen_fixture.R spikes <project_db> <taxdump_dir> [taxdump_date]
#
# Regenerating erases the provenance of numbers the tests assert against, so
# host mode refuses to clobber an existing fixture without --force.

suppressMessages(pkgload::load_all(quiet = TRUE))

FIX_DIR <- file.path("tests", "testthat", "fixtures")
HITS_F <- file.path(FIX_DIR, "screen_hits.csv.gz")
LAB_F <- file.path(FIX_DIR, "screen_labels.csv")

HIT_COLS <- c("unit", "hit_rank", "hsp", "bitscore", "hit_class",
  "lineage_taxids", "lineage")
LAB_COLS <- c("unit", "category", "expected_verdict", "expected_class", "source",
  "sample_taxids", "accession", "taxdump_date", "note")

#' Append rows to the fixture, creating it when absent
#' @noRd
fixture_append <- function(hits, labels) {
  dir.create(FIX_DIR, recursive = TRUE, showWarnings = FALSE)
  stopifnot(identical(names(hits), HIT_COLS), identical(names(labels), LAB_COLS))
  if (file.exists(HITS_F)) {
    old_h <- utils::read.csv(gzfile(HITS_F), stringsAsFactors = FALSE)
    old_l <- utils::read.csv(LAB_F, stringsAsFactors = FALSE)
    # a re-captured unit REPLACES its old rows; keeping both would silently
    # double its hits and change every score derived from them
    old_h <- old_h[!old_h$unit %in% hits$unit, , drop = FALSE]
    old_l <- old_l[!old_l$unit %in% labels$unit, , drop = FALSE]
    hits <- rbind(old_h, hits); labels <- rbind(old_l, labels)
  }
  con <- gzfile(HITS_F, "w"); utils::write.csv(hits, con, row.names = FALSE); close(con)
  utils::write.csv(labels, LAB_F, row.names = FALSE)
  cat("fixture:", nrow(labels), "units,", nrow(hits), "hit rows\n")
}

#' Capture every on_target unit of a real project as the clean-host corpus
#' @noRd
fixture_capture_host <- function(db_path, taxdump_date = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con))
  u <- DBI::dbGetQuery(con,
    "SELECT s.ID, s.path, s.scaffold, sm.lineage_taxids
       FROM screen s JOIN samples sm ON sm.ID = s.ID
      WHERE s.screen_verdict = 'on_target'
      ORDER BY s.ID, s.path, s.scaffold")
  if (nrow(u) == 0) stop("no on_target units in ", db_path, call. = FALSE)
  hits <- list(); labs <- list()
  for (i in seq_len(nrow(u))) {
    unit <- sprintf("host_%s_%d.%d", u$ID[i], u$path[i], u$scaffold[i])
    h <- DBI::dbGetQuery(con,
      "SELECT hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
         FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?
        ORDER BY hit_rank, hsp",
      params = list(u$ID[i], u$path[i], u$scaffold[i]))
    if (nrow(h) == 0) next   # a unit with no hits scores no_reference, not on_target
    hits[[length(hits) + 1]] <- cbind(unit = unit, h)
    labs[[length(labs) + 1]] <- data.frame(unit = unit, category = "host",
      expected_verdict = "on_target", expected_class = NA_character_,
      source = "real", sample_taxids = u$lineage_taxids[i],
      accession = NA_character_, taxdump_date = taxdump_date,
      note = "clean host scaffold from the reference project",
      stringsAsFactors = FALSE)
  }
  fixture_append(do.call(rbind, hits)[, HIT_COLS], do.call(rbind, labs)[, LAB_COLS])
  invisible(length(labs))
}

a <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% a
a <- a[a != "--force"]
if (length(a) >= 1 && a[1] == "host") {
  if (file.exists(HITS_F) && !force) stop("fixture exists; pass --force", call. = FALSE)
  if (file.exists(HITS_F)) { unlink(HITS_F); unlink(LAB_F) }
  fixture_capture_host(a[2], taxdump_date = if (length(a) >= 3) a[3] else NA_character_)
}
