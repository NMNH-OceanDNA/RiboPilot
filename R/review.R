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
      # an older schema mid-upgrade may lack mirrors a later step adds
      mirrors <- REVIEW_MIRRORS[REVIEW_MIRRORS %in% DBI::dbListTables(con)]
      for (m in mirrors) {
        DBI::dbExecute(con, sprintf("DELETE FROM %s WHERE ID = ? AND path = ? AND scaffold = ?", m),
          params = p)
      }
      n <- DBI::dbExecute(con, sprintf(
        "INSERT INTO review_units (%s, rev, time_stamp)
         SELECT %s, CAST(? AS INTEGER), ? FROM review_units_v
         WHERE ID = ? AND path = ? AND scaffold = ? AND has_copy = 0", cols, cols),
        params = c(list(review_new_rev(), as.integer(Sys.time())), p))
      if (n == 0) stop(sprintf("No sequence for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
      for (orig in names(mirrors)) {
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
         itsx_origin = NULL, itsx_scores = NULL, genbank_check = NULL, time_stamp = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(cnt$n_features, cnt$missing, as.integer(Sys.time())), p))
    if (n == 0) stop(sprintf("No review copy for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
    for (m in intersect(REVIEW_MIRRORS[c("blast_hits", "blast_screen", "screen", "chimera", "chimera_windows",
                               "blast_ref_alignment", "feature_blast", "feature_blast_hits")], DBI::dbListTables(con))) {
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
  with_savepoint(con, for (t in intersect(REVIEW_TABLES, DBI::dbListTables(con))) {
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

#' Trim a contig's review copy to positions from..to of the copy (S5.2)
#'
#' The sequence is cut from the original by the new window. Features are
#' shifted, clamped, or dropped by the window; a cut end gets
#' partial_*_user = 1 in the feature's own orientation. Makes the copy on first
#' touch and marks its checks not run.
#' @return invisible new length.
#' @noRd
review_trim <- function(con, id, path, scaffold, from, to) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  key <- "WHERE ID = ? AND path = ? AND scaffold = ?"
  invisible(with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    r <- DBI::dbGetQuery(con,
      "SELECT r.sequence, r.trim_start, a.sequence AS orig FROM review_units r
       JOIN assemblies a ON a.ID = r.ID AND a.path = r.path AND a.scaffold = r.scaffold
       WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?", params = p)
    if (nrow(r) == 0 || is.na(r$sequence[1]) || !nzchar(r$sequence[1])) {
      stop(sprintf("No sequence for contig %s.%d.%d.", id, path, scaffold), call. = FALSE)
    }
    from <- max(1L, as.integer(from)); to <- min(nchar(r$sequence[1]), as.integer(to))
    if (!isTRUE(from <= to)) stop("Invalid trim range.", call. = FALSE)
    w <- as.integer(r$trim_start[1]) + c(from, to) - 1L
    sq <- substr(r$orig[1], w[1], w[2]); n <- nchar(sq)
    DBI::dbExecute(con, paste("UPDATE review_units SET sequence = ?, length = ?, trim_start = ?, trim_end = ?,
      edited = 1", key), params = c(list(sq, n, w[1], w[2]), p))
    ann <- DBI::dbGetQuery(con, paste("SELECT element, start, stop, strand FROM review_annotations", key),
      params = p)
    fkey <- paste(key, "AND element = ?")
    for (i in seq_len(nrow(ann))) {
      fp <- c(p, list(ann$element[i]))
      ns <- ann$start[i] - from + 1L; nt <- ann$stop[i] - from + 1L
      if (nt < 1 || ns > n) {
        DBI::dbExecute(con, paste("DELETE FROM review_annotations", fkey), params = fp)
        next
      }
      minus <- identical(ann$strand[i], "-")
      cut <- c(if (ns < 1) (if (minus) "partial_stop_user" else "partial_start_user"),
        if (nt > n) (if (minus) "partial_start_user" else "partial_stop_user"))
      for (col in cut) {
        DBI::dbExecute(con, sprintf("UPDATE review_annotations SET %s = 1 %s", col, fkey), params = fp)
      }
      ns <- max(1L, ns); nt <- min(n, nt)
      DBI::dbExecute(con, paste("UPDATE review_annotations SET start = ?, stop = ?, length = ?, edited = 1", fkey),
        params = c(list(ns, nt, nt - ns + 1L), fp))
    }
    review_invalidate(con, id, path, scaffold)
    n
  }))
}

#' Write one feature of a copy and mark feature and copy edited
#' @noRd
review_edit_feature <- function(con, id, path, scaffold, element, set, vals) {
  n <- DBI::dbExecute(con, sprintf("UPDATE review_annotations SET %s, edited = 1
    WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?", set),
    params = c(vals, list(id, path, scaffold, element)))
  if (n == 0) stop(sprintf("No feature %s on contig %s.%d.%d.", element, id, path, scaffold), call. = FALSE)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
}

#' Move a feature's boundaries on the copy (Edit mode nudge); checks go to not run
#' @return invisible new rev, as text.
#' @noRd
review_set_bounds <- function(con, id, path, scaffold, element, start, stop) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  start <- as.integer(start); stop <- as.integer(stop)
  invisible(with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    review_edit_feature(con, id, path, scaffold, element, "start = ?, stop = ?, length = ?",
      list(start, stop, stop - start + 1L))
    review_invalidate(con, id, path, scaffold)
  }))
}

#' Set the user's 5' (start) or 3' (stop) partial choice on a copy's feature
#'
#' No check reads partial ends, so the check state is kept.
#' @noRd
review_set_partial <- function(con, id, path, scaffold, element, end = c("start", "stop"), value) {
  end <- match.arg(end)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  with_savepoint(con, {
    review_copy(con, id, path, scaffold)
    review_edit_feature(con, id, path, scaffold, element, sprintf("partial_%s_user = ?", end),
      list(as.integer(value)))
  })
  invisible()
}

#' Write an Edit session's starting bounds and user partials back (Reset, Discard)
#'
#' backup: data.frame(element, start, stop, partial_start_user,
#' partial_stop_user) read from review_annotations_v when Edit was ticked.
#' Features the copy no longer has are skipped.
#' @return invisible TRUE when a boundary moved back and the checks were
#'   marked not run again; FALSE otherwise, including without a copy.
#' @noRd
review_reset_session <- function(con, id, path, scaffold, backup) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  p <- list(id, path, scaffold)
  invisible(with_savepoint(con, if (!review_has_copy(con, id, path, scaffold)) FALSE else {
    cur <- DBI::dbGetQuery(con,
      "SELECT element, start, stop FROM review_annotations WHERE ID = ? AND path = ? AND scaffold = ?",
      params = p)
    b <- backup[backup$element %in% cur$element, , drop = FALSE]
    m <- match(b$element, cur$element)
    moved <- !identical(as.integer(c(cur$start[m], cur$stop[m])), as.integer(c(b$start, b$stop)))
    for (i in seq_len(nrow(b))) {
      DBI::dbExecute(con,
        "UPDATE review_annotations SET start = ?, stop = ?, length = ?, partial_start_user = ?,
           partial_stop_user = ?
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = c(list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, b$partial_start_user[i],
          b$partial_stop_user[i]), p, list(b$element[i])))
    }
    if (moved) review_invalidate(con, id, path, scaffold)
    moved
  }))
}

#' State on Annotation Review: 1 Ready to run (queue the recheck), 0 On hold
#'
#' Ready to run makes the copy if needed and marks its checks not run; On hold
#' sets check_switch 0 and keeps rev. Locked contigs are refused, and so is
#' Ready to run on a contig without a copy whose original the main pipeline
#' has queued or never annotated.
#' @param units data.frame(ID, path, scaffold).
#' @return character, NA where the state was set, else the reason.
#' @noRd
review_set_state <- function(con, units, state) {
  state <- as.integer(state)
  if (!isTRUE(state %in% 0:1)) stop("state must be 1 (Ready to run) or 0 (On hold).", call. = FALSE)
  vapply(seq_len(nrow(units)), function(i) {
    id <- units$ID[i]; path <- as.integer(units$path[i]); scaffold <- as.integer(units$scaffold[i])
    p <- list(id, path, scaffold)
    with_savepoint(con, {
      a <- DBI::dbGetQuery(con, "SELECT annotate_lock, annotate_switch, n_features FROM annotate
        WHERE ID = ? AND path = ? AND scaffold = ?", params = p)
      if (isTRUE(a$annotate_lock[1] == 1)) {
        "unlock first"
      } else if (state == 1L && !review_has_copy(con, id, path, scaffold) &&
                 (a$annotate_switch[1] %in% c(1, 4) || is.na(a$n_features[1]))) {
        "press Update on the Assemble step first"
      } else {
        review_copy(con, id, path, scaffold)
        if (state == 1L) {
          review_invalidate(con, id, path, scaffold)
        } else {
          DBI::dbExecute(con, "UPDATE review_units SET check_switch = 0, time_stamp = ?
            WHERE ID = ? AND path = ? AND scaffold = ?", params = c(list(as.integer(Sys.time())), p))
        }
        NA_character_
      }
    })
  }, character(1))
}

#' Table a check function reads or writes
#'
#' The name itself when review is FALSE. In review mode the mirrored tables
#' map to their review_* mirror, and assemblies (sequence, length, blast_*) and
#' annotate (result columns only) map to review_units. Shared tables and
#' coverage are never routed, so asking for one is an error.
#' @noRd
review_tbl <- function(name, review = FALSE) {
  if (!review) return(name)
  if (name %in% names(REVIEW_MIRRORS)) return(REVIEW_MIRRORS[[name]])
  if (name %in% c("assemblies", "annotate")) return("review_units")
  stop("No review table for ", name, ".", call. = FALSE)
}

#' Follow a main-pipeline rerun on a contig's review copy (S5.3)
#'
#' The pipeline's per-contig scripts call this after they rewrite a contig's
#' originals, and option changes call it at once. An unlocked copy made only
#' by State (edited 0) is deleted, so the contig reads its new originals; an
#' unlocked edited copy keeps its edits and its checks go back to not run; a
#' locked copy is left alone, since the recheck skips locked contigs.
#' @return invisible "none" (no copy), "locked", "dropped", or "invalidated".
#' @noRd
review_follow_original <- function(con, id, path, scaffold) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  if (!review_has_copy(con, id, path, scaffold)) return(invisible("none"))
  invisible(with_savepoint(con, {
    u <- DBI::dbGetQuery(con,
      "SELECT r.edited, COALESCE(an.annotate_lock, 0) AS locked FROM review_units r
       LEFT JOIN annotate an ON an.ID = r.ID AND an.path = r.path AND an.scaffold = r.scaffold
       WHERE r.ID = ? AND r.path = ? AND r.scaffold = ?", params = list(id, path, scaffold))
    if (nrow(u) == 0) {
      "none"
    } else if (u$locked[1] == 1) {
      "locked"
    } else if (u$edited[1] == 0) {
      review_drop(con, id, path, scaffold)
      "dropped"
    } else {
      review_invalidate(con, id, path, scaffold)
      "invalidated"
    }
  }))
}

#' Check state of contigs, as Annotation Review shows it
#'
#' One row per kept contig of the assembly-locked samples, ordered by ID, path,
#' scaffold; or, with `units` (data.frame(ID, path, scaffold)), one row per
#' unit in its order, matched in R. "waiting": no copy, and the original is
#' queued for the main pipeline or was never annotated (the contigs
#' nf_recheck_waiting() counts). Otherwise check_switch decides: 2 done,
#' 1 or NULL not run, 4 running, 3 failed, 0 on hold.
#' @noRd
review_unit_state <- function(con, units = NULL) {
  d <- DBI::dbGetQuery(con, paste(
    "SELECT v.ID, v.path, v.scaffold, v.has_copy, v.check_switch, an.annotate_switch,
       an.annotate_lock, an.n_features, v.check_note
     FROM review_units_v v
     JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold",
    if (is.null(units)) {
      "JOIN assemblies asm ON asm.ID = v.ID AND asm.path = v.path AND asm.scaffold = v.scaffold
     JOIN assemble a ON a.ID = v.ID
     WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0"
    },
    "ORDER BY v.ID, v.path, v.scaffold"))
  if (!is.null(units)) {
    d <- d[match(unit_key(units$ID, units$path, units$scaffold), unit_key(d$ID, d$path, d$scaffold)), ,
      drop = FALSE]
    d$ID <- as.character(units$ID)
    d$path <- as.integer(units$path)
    d$scaffold <- as.integer(units$scaffold)
    rownames(d) <- NULL
  }
  waiting <- d$has_copy %in% 0 & (d$annotate_switch %in% c(1, 4) | is.na(d$n_features))
  st <- unname(c(`0` = "on hold", `1` = "not run", `2` = "done", `3` = "failed",
    `4` = "running")[as.character(d$check_switch)])
  d$state <- as.character(ifelse(waiting, "waiting", ifelse(is.na(st), "not run", st)))
  d
}

WAITING_MSG <- "Press Update on the Assemble step first."

#' Words for review_unit_state()'s state, as the table cells show them
#' @noRd
review_state_text <- function(state) {
  unname(c(done = "done", `not run` = "not run", running = "running", failed = "failed",
    `on hold` = "on hold", waiting = "waiting for Assemble Update")[state])
}

#' Why each unit cannot be locked on Annotation Review (S6.2)
#'
#' NA where Lock may go ahead: the checks are done, or the contig is excluded
#' (an excluded contig is never exported). Otherwise the next step to take.
#' @noRd
review_lock_gate <- function(con, units) {
  if (!NROW(units)) return(character(0))
  st <- review_unit_state(con, units)
  ex <- DBI::dbGetQuery(con, "SELECT ID, path, scaffold FROM annotate WHERE problematic = 'yes'")
  again <- "Set State to Ready to run, then press Update."
  msg <- unname(c(waiting = WAITING_MSG, `not run` = again,
    running = "Wait for the running Update to finish.", failed = again, `on hold` = again,
    done = NA_character_)[st$state])
  msg[st$state == "not run" & st$has_copy %in% 1] <- "Press Update first."
  msg[unit_key(units$ID, units$path, units$scaffold) %in% unit_key(ex$ID, ex$path, ex$scaffold)] <- NA_character_
  msg
}
