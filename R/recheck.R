# The recheck: reruns every offline check on one contig's review copy
# (dev/specs/2026-10-08-contig-recheck-design.md, section 7.3).

#' blastn settings of the local search in inst/nextflow/modules/blast.nf
#'
#' The recheck repeats that search in R. test-recheck.R checks these values
#' against blast.nf, so the two cannot drift.
#' @noRd
RECHECK_BLAST <- list(
  outfmt = "6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid",
  task = "megablast",
  retry_task = "dc-megablast",
  max_hsps = 5L,
  retry_below_qcovs = 25,
  default_max_target_seqs = 50L
)

#' Local megablast of one sequence, with blast.nf's dc-megablast retry
#'
#' Retries with dc-megablast when megablast finds nothing or its best query
#' coverage is under `RECHECK_BLAST$retry_below_qcovs`, and keeps the retry only
#' when it exits 0 with output. Never searches NCBI.
#' @param seq Query sequence.
#' @param dbarg blastn `-db` value (rp_dbarg.R's panel list), passed as one argument.
#' @param max_target_seqs,extra_opts The sample's BLAST options; `extra_opts`
#'   is split on whitespace, as blast.nf's unquoted `${XO}` is.
#' @param workdir Folder for query.fa and the hits files.
#' @return Path of the hits file, or NULL when megablast failed.
#' @noRd
recheck_blast_local <- function(seq, dbarg, max_target_seqs, extra_opts, workdir,
                                threads = 1L, blastn = "blastn") {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  query <- file.path(workdir, "query.fa")
  writeLines(c(">query", seq), query)
  xo <- strsplit(trimws(extra_opts %|NA|% ""), "[[:space:]]+")[[1]]
  run <- function(task, out) {
    st <- suppressWarnings(system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query),
      "-outfmt", shQuote(RECHECK_BLAST$outfmt), "-max_target_seqs", max_target_seqs,
      "-max_hsps", RECHECK_BLAST$max_hsps, "-task", task, "-num_threads", threads,
      shQuote(xo[nzchar(xo)])), stdout = out, stderr = FALSE))
    identical(as.integer(st), 0L)
  }
  hits <- file.path(workdir, "hits.txt")
  if (!run(RECHECK_BLAST$task, hits)) return(NULL)
  if (file.size(hits) == 0) writeLines("NO_SIGNIFICANT_HITS", hits)
  h <- parse_blast_hits_file(hits)
  if (nrow(h) == 0 || max(h$qcovs, 0, na.rm = TRUE) < RECHECK_BLAST$retry_below_qcovs) {
    dc <- file.path(workdir, "hits_dc.txt")
    if (run(RECHECK_BLAST$retry_task, dc) && isTRUE(file.size(dc) > 0)) file.rename(dc, hits)
  }
  hits
}

#' Fresh BLAST hits for a review copy (SPEC S7.3 step 2)
#'
#' `local` searches the panel again with the sample's options. Without a panel
#' (`dbarg` NULL or "") or when that search fails, it maps the original hits
#' instead, as `map` does. A copy whose original has no hits at all keeps the
#' BLAST summary it copied (a seed-database name, for example). Writes
#' review_blast_hits, review_blast_screen, and review_units.blast_*.
#' @return "local" or "map", whichever was used.
#' @noRd
recheck_hits <- function(con, id, path, scaffold, mode, dbarg = NULL, taxdump_dir = NULL,
                         threads = 1L, workdir = tempfile("recheck_")) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  unit <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  p <- list(id, path, scaffold)
  if (identical(mode, "local") && isTRUE(nzchar(dbarg, keepNA = TRUE))) {
    o <- DBI::dbGetQuery(con, "SELECT b.extra_opts, b.max_target_seqs FROM assemble a
      JOIN blast_opts b ON a.blast_opts = b.blast_opts WHERE a.ID = ?", params = list(id))
    seq <- DBI::dbGetQuery(con, paste("SELECT sequence FROM review_units", unit), params = p)$sequence
    hits <- recheck_blast_local(seq, dbarg,
      o$max_target_seqs[1] %|NA|% RECHECK_BLAST$default_max_target_seqs,
      o$extra_opts[1] %|NA|% "", workdir, threads)
    if (!is.null(hits)) {
      parse_blast_candidates(con, id, path, scaffold, hits, taxdump_dir = taxdump_dir, review = TRUE)
      return("local")
    }
  }
  if (DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM blast_hits", unit), params = p)$n == 0) {
    with_savepoint(con, for (t in c("review_blast_hits", "review_blast_screen")) {
      DBI::dbExecute(con, paste("DELETE FROM", t, unit), params = p)
    })
    return("map")
  }
  m <- review_map_hits(con, id, path, scaffold)
  with_savepoint(con, {
    DBI::dbExecute(con, paste("DELETE FROM review_blast_hits", unit), params = p)
    DBI::dbAppendTable(con, "review_blast_hits", m)
  })
  derive_blast_products(con, id, path, scaffold, review = TRUE)
  "map"
}

#' Refresh the model coverage of a copy's barrnap features (SPEC S7.3 step 5)
#'
#' Each barrnap-sourced feature (tool "barrnap", "barrnap+ITSx", ...) takes the
#' coverage of the barrnap hit with its element and strand that overlaps it
#' most; a feature no hit overlaps gets 0. Other features are left alone.
#' @param gff barrnap GFF of the copy's sequence, or NULL to change nothing.
#' @noRd
recheck_model_coverage <- function(con, id, path, scaffold, gff) {
  if (is.null(gff)) return(invisible())
  p <- list(id, as.integer(path), as.integer(scaffold))
  bar <- if (file.exists(gff)) parse_barrnap(readLines(gff, warn = FALSE)) else data.frame()
  f <- DBI::dbGetQuery(con, "SELECT element, start, stop, strand, tool FROM review_annotations
    WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  f <- f[which(startsWith(f$tool, "barrnap")), , drop = FALSE]
  cov <- vapply(seq_len(nrow(f)), function(i) {
    b <- if (nrow(bar)) bar[bar$element == f$element[i] & bar$strand == f$strand[i], , drop = FALSE] else bar
    ov <- if (nrow(b)) pmin(b$stop, f$stop[i]) - pmax(b$start, f$start[i]) + 1 else 0
    if (max(ov) < 1) 0 else b$model_cov[which.max(ov)]
  }, numeric(1))
  with_savepoint(con, for (i in seq_len(nrow(f))) {
    DBI::dbExecute(con, "UPDATE review_annotations SET model_coverage = ?
      WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(list(cov[i]), p, list(f$element[i])))
  })
  invisible(cov)
}

#' Recheck one contig's review copy (SPEC S7.3)
#'
#' Reruns every offline check on the copy in the order of the spec, then marks
#' it done, but only while the copy still has the `rev` the run was queued with
#' and is still queued or running (1 or 4); an edit, On hold, or Restore during
#' the run discards the result. An error in a step marks the copy failed under
#' the same guard. Opens no transaction around the steps: each check opens its
#' own.
#' @param rev The copy's `rev` as text, as query_recheck emitted it.
#' @param mode "local" (search the panel again) or "map" (map the original hits).
#' @param tools_dir Folder holding rp_rrna_tools.sh's output for the copy.
#' @return "done", "stale", or "failed", invisibly.
#' @noRd
recheck_unit <- function(con, id, path, scaffold, rev, mode, tools_dir, dbarg = NULL,
                         taxdump_dir = NULL, label_db = NULL, threads = 1L,
                         workdir = tempfile("recheck_")) {
  stopifnot(is.character(rev), length(rev) == 1)
  if (!isTRUE(nzchar(dbarg, keepNA = TRUE))) dbarg <- NULL
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  unit <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  guard <- paste(unit, "AND rev = CAST(? AS INTEGER) AND check_switch IN (1, 4)")
  p <- list(id, path, scaffold)
  pg <- c(p, list(as.character(rev)))
  if (DBI::dbGetQuery(con, paste("SELECT COUNT(*) AS n FROM review_units", guard), params = pg)$n == 0) {
    return(invisible("stale"))
  }
  st <- file.path(tools_dir, "unit.tools_status")
  st <- if (file.exists(st)) trimws(readLines(st, warn = FALSE)) else character(0)
  ok <- function(tool) paste0(tool, "=0") %in% st
  err <- tryCatch({
    recheck_hits(con, id, path, scaffold, mode, dbarg = dbarg, taxdump_dir = taxdump_dir,
      threads = threads, workdir = file.path(workdir, "blast"))
    screen_one(con, id, path, scaffold, force = TRUE, review = TRUE)
    cnt <- review_counts(con, id, path, scaffold)
    DBI::dbExecute(con, paste("UPDATE review_units SET n_features = ?, missing = ?", unit),
      params = c(list(cnt$n_features, cnt$missing), p))
    recheck_model_coverage(con, id, path, scaffold,
      if (ok("barrnap")) file.path(tools_dir, "unit.barrnap.gff"))
    if (ok("itsx")) {
      pos <- readLines(file.path(tools_dir, "unit.itsx.positions.txt"), warn = FALSE)
      DBI::dbExecute(con, paste("UPDATE review_units SET chimeric = ?", unit),
        params = c(list(as.integer(itsx_is_chimeric(pos))), p))
    }
    tryCatch(itsx_check_one(con, id, path, scaffold,
      out_prefix = file.path(tools_dir, "unit.itsx_check"), review = TRUE),
      error = function(e) itsx_store_not_run(con, id, path, scaffold,
        paste("ITSx check failed:", conditionMessage(e)), review = TRUE))
    compute_ref_alignment(con, id, path, scaffold, review = TRUE)
    curate_one(con, id, path, scaffold, review = TRUE)
    chimera_unit(con, id, path, scaffold, dbarg = dbarg, taxdump_dir = taxdump_dir,
      label_db = label_db, threads = threads, workdir = file.path(workdir, "chimera"), review = TRUE)
    NULL
  }, error = function(e) conditionMessage(e))
  ts <- as.integer(Sys.time())
  if (is.null(err)) {
    n <- DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 2, checked_at = ?,
      check_note = NULL, time_stamp = ?", guard), params = c(list(ts, ts), pg))
    return(invisible(if (n == 1) "done" else "stale"))
  }
  n <- DBI::dbExecute(con, paste("UPDATE review_units SET check_switch = 3, check_note = ?, time_stamp = ?",
    guard), params = c(list(err, ts), pg))
  invisible(if (n == 1) "failed" else "stale")
}

#' After an in-app recheck run: copies left running go back to not run
#'
#' @return The number of copies released.
#' @noRd
recheck_release_running <- function(con) {
  DBI::dbExecute(con, "UPDATE review_units SET check_switch = 1, check_note = 'recheck did not finish',
    time_stamp = ? WHERE check_switch = 4", params = list(as.integer(Sys.time())))
}
