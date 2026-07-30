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

# The spike mode runs inside ribopilot:test, which has RiboPilot installed but
# not pkgload; the host mode runs from a source checkout. Support both.
if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  suppressMessages(pkgload::load_all(quiet = TRUE))
} else {
  suppressMessages(library(RiboPilot))
}

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

# --- spike-in capture (network; run inside ribopilot:test) ---------------------

# Accessions are RESOLVED at generation time, never hardcoded: baking in an
# accession nobody verified is how a fixture starts lying. The resolved
# accession is written into the label so what was fetched stays on the record.
SPIKES <- data.frame(
  key   = c("fungal", "bacterial", "algal", "arthropod", "human"),
  query = c(
    "Fusarium[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Escherichia coli[Organism] AND 16S ribosomal RNA[Title] AND 1200:1700[SLEN]",
    "Arabidopsis thaliana[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Calanus[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Homo sapiens[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]"),
  stringsAsFactors = FALSE)

#' Resolve one accession for a query
#' @noRd
spike_accession <- function(query) {
  out <- suppressWarnings(system2("sh", c("-c", shQuote(sprintf(
    "esearch -db nuccore -query %s | efetch -format acc | head -1",
    shQuote(query)))), stdout = TRUE, stderr = FALSE))
  acc <- trimws(paste(out, collapse = ""))
  if (!nzchar(acc)) stop("no accession for query: ", query, call. = FALSE)
  acc
}

#' The same remote blastn the pipeline runs (inst/nextflow/modules/blast.nf)
#'
#' Wrapped in `timeout`. `blastn -remote` blocks indefinitely when the NCBI URL
#' API stops answering, and it does: one capture sat 85 minutes at zero CPU
#' waiting on a query that normally takes minutes. Failing fast and retrying
#' beats an unbounded wait nobody can distinguish from progress.
#' @noRd
spike_blast <- function(fasta, out_file, secs = 900L) {
  fmt <- paste("6 saccver stitle pident qcovs evalue staxids bitscore length",
    "qstart qend sstart send")
  ok <- system2("timeout", c(as.character(secs), "blastn",
    "-remote", "-db", "core_nt", "-query", shQuote(fasta),
    "-outfmt", shQuote(fmt), "-max_target_seqs", "50", "-max_hsps", "5",
    "-task", "megablast"), stdout = out_file)
  if (identical(as.integer(ok), 124L)) {
    stop("blastn timed out after ", secs, "s on ", basename(fasta),
      "; NCBI is not answering, retry later", call. = FALSE)
  }
  if (!identical(as.integer(ok), 0L) || !file.exists(out_file) ||
      file.size(out_file) == 0) {
    stop("blastn returned nothing for ", fasta, call. = FALSE)
  }
  invisible(out_file)
}

#' Is this a usable nucleotide FASTA, not an error page or a truncated record?
#'
#' `file.size() > 0` is too weak: eutils answers a dropped connection with an
#' HTML body and HTTP 200, which would be accepted and BLASTed as if it were
#' sequence.
#' @noRd
spike_fasta_ok <- function(fa) {
  if (!file.exists(fa) || file.size(fa) == 0) return(FALSE)
  ln <- tryCatch(readLines(fa, warn = FALSE), error = function(e) character(0))
  if (length(ln) < 2 || !startsWith(ln[1], ">")) return(FALSE)
  seq <- paste(ln[-1], collapse = "")
  nchar(gsub("[^ACGTNacgtn]", "", seq)) >= 500
}

#' Capture each spike-in's hit distribution through the pipeline's own parser
#'
#' A scratch DB is used because parse_blast_candidates() writes into blast_hits
#' and needs the real schema. The scratch sample's lineage does not decide
#' anything: scoring happens in the test, against whatever host lineage the
#' label carries.
#' @noRd
fixture_capture_spikes <- function(db_path, taxdump_dir, taxdump_date = NA_character_) {
  hc <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  fish <- DBI::dbGetQuery(hc,
    "SELECT lineage_taxids t FROM samples WHERE lineage_taxids IS NOT NULL LIMIT 1")$t[1]
  DBI::dbDisconnect(hc)
  if (is.na(fish) || !nzchar(fish)) stop("no host lineage in ", db_path, call. = FALSE)

  # A COPY of the real project DB, not a hand-built one: parse_blast_candidates()
  # writes blast_hits, blast_ref_candidates and blast_screen, and reproducing
  # that schema here would be a second copy of the DDL free to drift from it.
  # The copy is scratch, so the real project is never touched.
  tmp <- tempfile("rp_fix_"); dir.create(tmp)
  scratch <- file.path(tmp, ".sqlite")
  if (!file.copy(db_path, scratch)) stop("could not copy ", db_path, call. = FALSE)
  sc <- DBI::dbConnect(RSQLite::SQLite(), scratch)
  on.exit(DBI::dbDisconnect(sc), add = TRUE)
  # The copy arrives carrying the project's OWN blast_hits, which share the
  # path/scaffold numbers the spikes use. Clear them: a readback that finds real
  # host hits mixed into a spike produces a distribution labelled as a
  # contaminant while actually being the host, which is worse than no fixture.
  DBI::dbExecute(sc, "DELETE FROM blast_hits")
  DBI::dbExecute(sc, "INSERT INTO samples (ID, lineage_taxids) VALUES('SPIKE',?)",
    params = list(fish))

  for (i in seq_len(nrow(SPIKES))) {
    key <- SPIKES$key[i]
    acc <- spike_accession(SPIKES$query[i])
    fa <- file.path(tmp, paste0(key, ".fa"))
    # eutils drops connections intermittently, so retry rather than lose the run
    for (try in 1:3) {
      system2("efetch", c("-db", "nuccore", "-id", acc, "-format", "fasta"), stdout = fa)
      if (spike_fasta_ok(fa)) break
      Sys.sleep(5 * try)
    }
    if (!spike_fasta_ok(fa)) stop("efetch gave no usable fasta for ", acc, call. = FALSE)
    cat("spike", key, "accession", acc, "\n")
    hf <- file.path(tmp, paste0(key, ".hits"))
    spike_blast(fa, hf)
    parse_blast_candidates(sc, "SPIKE", path = 1L, scaffold = i,
      hits_file = hf, taxdump_dir = taxdump_dir)
    h <- DBI::dbGetQuery(sc,
      "SELECT hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
         FROM blast_hits WHERE ID = 'SPIKE' AND path = 1 AND scaffold = ?
        ORDER BY hit_rank, hsp",
      params = list(i))
    if (nrow(h) == 0) stop("no parsed hits for spike ", key, call. = FALSE)
    cat("  parsed", nrow(h), "hit rows\n")
    unit <- paste0(key, "_on_fish")
    # append per spike, not once at the end: remote BLAST is flaky, and losing
    # four completed searches because the fifth dropped its connection is a cost
    # nobody should pay twice
    fixture_append(cbind(unit = unit, h)[, HIT_COLS],
      data.frame(unit = unit, category = key,
        expected_verdict = NA_character_, expected_class = NA_character_,
        source = "spike", sample_taxids = fish, accession = acc,
        taxdump_date = taxdump_date, note = "",
        stringsAsFactors = FALSE)[, LAB_COLS])
  }
  invisible(nrow(SPIKES))
}

if (length(a) >= 1 && a[1] == "spikes") {
  fixture_capture_spikes(a[2], a[3],
    taxdump_date = if (length(a) >= 4) a[4] else NA_character_)
}
