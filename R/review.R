# Review copies: the per-contig working copy that Annotation Review edits and
# the recheck scores, so the originals keep the main pipeline's results.
# Design: dev/specs/2026-10-08-contig-recheck-design.md (S4, S5).

#' Effective partial ends on the review layer: the user's choice, else curation's
#' @noRd
REVIEW_EFF_PARTIALS <- paste("COALESCE(partial_start_user, partial_start) AS partial_start,",
  "COALESCE(partial_stop_user, partial_stop) AS partial_stop")

#' Run code atomically, inside or outside an open transaction
#'
#' With no transaction open it runs BEGIN IMMEDIATE ... COMMIT, so the write
#' lock is taken before the first read and busy_timeout applies. Inside an open
#' transaction it runs SAVEPOINT ... RELEASE. An error, or leaving code early
#' (return() inside code, an interrupt), rolls back and leaves no transaction or
#' savepoint open; the error propagates. Code inside must not call
#' DBI::dbWithTransaction(). Returns the value of code.
#' @noRd
with_savepoint <- function(con, code, name = "rp_review") {
  outer <- !RSQLite::sqliteIsTransacting(con)
  DBI::dbExecute(con, if (outer) "BEGIN IMMEDIATE" else paste0("SAVEPOINT ", name))
  ok <- FALSE
  on.exit(if (!ok) {
    if (outer) {
      try(DBI::dbExecute(con, "ROLLBACK"), silent = TRUE)
    } else {
      try(DBI::dbExecute(con, paste0("ROLLBACK TO ", name)), silent = TRUE)
      try(DBI::dbExecute(con, paste0("RELEASE ", name)), silent = TRUE)
    }
  })
  res <- force(code)
  DBI::dbExecute(con, if (outer) "COMMIT" else paste0("RELEASE ", name))
  ok <- TRUE
  res
}

#' A new review_units.rev: the current time in milliseconds, as text
#'
#' rev does not fit an R integer; carry it as text and compare it in SQL
#' (rev = CAST(? AS INTEGER)).
#' @noRd
review_new_rev <- function() format(round(as.numeric(Sys.time()) * 1000), scientific = FALSE)

#' @noRd
review_has_copy <- function(con, id, path, scaffold) {
  nrow(DBI::dbGetQuery(con,
    "SELECT 1 AS x FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))) > 0
}

#' Make a contig's review copy from its originals, unless it has one (S5.1)
#'
#' Rows a discarded recheck left in the mirrors are deleted first. The
#' review_units row comes from the review_units_v fallback row, so the copy
#' starts with the check state the view showed. Mirror rows are copied by column
#' name, since column order differs between new_db(), upgraded databases, and
#' hand-built ones.
#' @return invisible TRUE when a copy was made, FALSE when one existed.
#' @noRd
review_copy <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  cols <- paste(setdiff(ddl_cols(REVIEW_UNITS_DDL), c("check_note", "rev", "checked_at", "time_stamp")),
    collapse = ", ")
  invisible(with_savepoint(con, {
    made <- !review_has_copy(con, id, path, scaffold)
    if (made) {
      for (m in REVIEW_MIRRORS) {
        DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", m),
          params = p)
      }
      n <- DBI::dbExecute(con, sprintf(
        "INSERT INTO review_units (%s, rev, time_stamp)
         SELECT %s, CAST(? AS INTEGER), ? FROM review_units_v
         WHERE ID = ? AND path = ? AND scaffold = ? AND has_copy = 0", cols, cols),
        params = c(list(review_new_rev(), as.integer(Sys.time())), p))
      if (n == 0) stop(sprintf("No sequence for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
      for (orig in names(REVIEW_MIRRORS)) {
        mc <- paste(intersect(DBI::dbListFields(con, orig), DBI::dbListFields(con, REVIEW_MIRRORS[[orig]])),
          collapse = ", ")
        DBI::dbExecute(con, sprintf(
          "INSERT INTO %s (%s) SELECT %s FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
          REVIEW_MIRRORS[[orig]], mc, mc, orig), params = p)
      }
    }
    made
  }))
}

#' Feature count and missing elements of a copy, as reconcile_one counts them
#' @noRd
review_counts <- function(con, id, path, scaffold) {
  el <- DBI::dbGetQuery(con,
    "SELECT element FROM review_annotations WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))$element
  list(n_features = length(el),
    missing = paste(setdiff(c(RRNA_GENES, RRNA_SPACERS), el), collapse = ","))
}

#' Mark a copy's checks not run after an edit (S5.2)
#'
#' Raises rev, clears the check results and the curation-owned feature columns,
#' and deletes the copy's check rows. Counts are recomputed, not cleared;
#' chimeric and blast_* stay as copied until the recheck rewrites them.
#' @return invisible new rev, as text.
#' @noRd
review_invalidate <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  invisible(with_savepoint(con, {
    cnt <- review_counts(con, id, path, scaffold)
    n <- DBI::dbExecute(con,
      "UPDATE review_units SET check_switch = 1, rev = rev + 1, check_note = NULL, checked_at = NULL,
         n_features = ?, missing = ?, n_flags = NULL, curate_flags = NULL, itsx_warnings = NULL,
         itsx_origin = NULL, itsx_scores = NULL, time_stamp = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(cnt$n_features, cnt$missing, as.integer(Sys.time())), p))
    if (n == 0) stop(sprintf("No review copy for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
    for (m in REVIEW_MIRRORS[c("blast_hits", "blast_screen", "screen", "chimera", "chimera_windows",
                               "blast_ref_alignment")]) {
      DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", m),
        params = p)
    }
    DBI::dbExecute(con,
      "UPDATE review_annotations SET status = NULL, flags = NULL, partial_start = NULL, partial_stop = NULL
       WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
    DBI::dbGetQuery(con,
      "SELECT CAST(rev AS TEXT) AS rev FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?",
      params = p)$rev
  }))
}

#' Delete a contig's review copy, so it reads its originals again (S5.2 Restore)
#' @noRd
review_drop <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  with_savepoint(con, for (t in REVIEW_TABLES) {
    DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", t), params = p)
  })
  invisible()
}


#' A contig's read coverage on the review layer (S4.5)
#'
#' The original coverage tracks, cut to the copy's kept window when the contig
#' has a copy. A short track is cut, never padded. A track is NULL when the
#' coverage row or the track is missing, or nothing of it lies in the window.
#' Coverage is read with SELECT * because older projects and test fixtures lack
#' the clips column.
#' @return list(depths, clips): numeric vectors in copy coordinates, or NULL.
#' @noRd
review_coverage <- function(con, id, path, scaffold) {
  p <- list(id, as.integer(path), as.integer(scaffold))
  cv <- DBI::dbGetQuery(con, "SELECT * FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  w <- DBI::dbGetQuery(con,
    "SELECT trim_start, trim_end FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  cut <- function(col) {
    if (!nrow(cv) || !col %in% names(cv)) return(NULL)
    v <- parse_depths(cv[[col]][1])
    if (nrow(w)) v <- v[seq_along(v) >= w$trim_start & seq_along(v) <= w$trim_end]
    if (length(v)) v
  }
  list(depths = cut("depths"), clips = cut("clips"))
}

#' Coverage summary of a depth track, rounded as write_coverage_summary() rounds it
#' @return list(mean_depth, median_depth, min_depth, max_depth), or NULL for an
#'   empty track.
#' @noRd
review_coverage_summary <- function(dep) {
  if (!length(dep)) return(NULL)
  list(mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2),
    min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)))
}

#' A copy's BLAST hits mapped from the original's (S4.5)
#'
#' The original HSPs that overlap the copy's kept window, clipped to it and
#' shifted to copy coordinates. A clipped HSP gets aln_len = its clipped span
#' and NA sstart and send; pident and bitscore stay as stored. qcovs is
#' recomputed per subject as the union of its kept spans over the copy length,
#' rounded as BLAST rounds it. hit_rank and hsp keep their stored values, so a
#' subject or HSP clipped away leaves a gap.
#' @return data.frame with the blast_hits columns; zero rows when nothing overlaps.
#' @noRd
review_map_hits <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  w <- DBI::dbGetQuery(con,
    "SELECT trim_start, trim_end FROM review_units WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
  if (!nrow(w)) stop(sprintf("No review copy for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
  h <- DBI::dbGetQuery(con, sprintf(
    "SELECT %s FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY hit_rank, hsp",
    paste(ddl_cols(BLAST_HITS_DDL), collapse = ", ")), params = p)
  h <- h[!is.na(h$qstart) & !is.na(h$qend) & h$qend >= w$trim_start & h$qstart <= w$trim_end, , drop = FALSE]
  qs <- pmax(h$qstart, w$trim_start)
  qe <- pmin(h$qend, w$trim_end)
  cut <- qs != h$qstart | qe != h$qend
  h$aln_len[cut] <- qe[cut] - qs[cut] + 1L
  h$sstart[cut] <- NA
  h$send[cut] <- NA
  h$qstart <- qs - w$trim_start + 1L
  h$qend <- qe - w$trim_start + 1L
  span <- stats::ave(seq_len(nrow(h)), h$hit_rank,
    FUN = function(i) length(unique(unlist(Map(seq, h$qstart[i], h$qend[i])))))
  h$qcovs <- floor(100 * span / (w$trim_end - w$trim_start + 1) + 0.5)
  rownames(h) <- NULL
  h
}
