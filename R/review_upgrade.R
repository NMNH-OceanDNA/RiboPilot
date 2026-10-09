# Schema step 4 (version 4 to 5): the trims and edits an older RiboPilot wrote
# into the original tables become review copies, and trimmed originals are put
# back from assembly_backup (SPEC S11.3-11.4).

#' Convert old trims and edits into review copies (schema step 4)
#'
#' Runs inside the upgrade transaction. Each contig converts in its own
#' savepoint, so a contig that fails is rolled back and listed instead of
#' stopping the upgrade. One message lists the converted contigs, the ones with
#' a problem or left as they were, and the locked ones marked done.
#' @noRd
review_upgrade <- function(con) {
  uid <- function(d) paste(d$ID, d$path, d$scaffold, sep = ".")
  u <- DBI::dbGetQuery(con,
    "SELECT a.ID, a.path, a.scaffold, a.sequence AS cur, b.ID IS NOT NULL AS has_bak,
       b.sequence AS bak, b.annotations AS bak_ann, b.coverage AS bak_cov, b.chimera AS bak_chim,
       COALESCE(an.annotate_lock, 0) AS locked
     FROM assemblies a
     JOIN annotate an ON an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold
     LEFT JOIN assembly_backup b ON b.ID = a.ID AND b.path = a.path AND b.scaffold = a.scaffold
     WHERE (b.sequence IS NOT NULL AND a.sequence IS NOT b.sequence)
        OR EXISTS (SELECT 1 FROM annotations f WHERE f.ID = a.ID AND f.path = a.path
                     AND f.scaffold = a.scaffold AND f.edited = 1)
     ORDER BY a.ID, a.path, a.scaffold")
  n <- 0L
  bad <- skip <- done <- character()
  for (i in seq_len(nrow(u))) {
    r <- tryCatch(with_savepoint(con, review_upgrade_unit(con, u[i, ]), name = "rp_upgrade"),
      error = function(e) e)
    if (inherits(r, "error")) {
      skip <- c(skip, sprintf("%s (%s)", uid(u[i, ]), conditionMessage(r)))
    } else {
      n <- n + 1L
      if (length(r)) bad <- c(bad, sprintf("%s (%s)", uid(u[i, ]), paste(r, collapse = "; ")))
    }
  }

  lk <- DBI::dbGetQuery(con,
    "SELECT v.ID, v.path, v.scaffold FROM review_units_v v
     JOIN annotate an ON an.ID = v.ID AND an.path = v.path AND an.scaffold = v.scaffold
     WHERE v.has_copy = 0 AND an.annotate_lock = 1 AND v.check_switch IS NOT 2
     ORDER BY v.ID, v.path, v.scaffold")
  for (i in seq_len(nrow(lk))) {
    x <- lk[i, ]
    r <- tryCatch(with_savepoint(con, {
      review_copy(con, x$ID, x$path, x$scaffold)
      review_upgrade_done(con, x$ID, x$path, x$scaffold)
    }, name = "rp_upgrade"), error = function(e) e)
    if (inherits(r, "error")) {
      skip <- c(skip, sprintf("%s (%s)", uid(x), conditionMessage(r)))
    } else {
      done <- c(done, uid(x))
    }
  }

  left <- uid(DBI::dbGetQuery(con,
    "SELECT b.ID, b.path, b.scaffold FROM assembly_backup b
     JOIN assemblies a ON a.ID = b.ID AND a.path = b.path AND a.scaffold = b.scaffold
     WHERE b.sequence IS NOT NULL AND a.sequence IS NOT b.sequence
       AND NOT EXISTS (SELECT 1 FROM annotate an
         WHERE an.ID = b.ID AND an.path = b.path AND an.scaffold = b.scaffold)
     ORDER BY b.ID, b.path, b.scaffold"))

  txt <- c(
    if (n) sprintf(paste("Moved the edits of %d contig(s) to review copies on Annotation Review; the",
      "Assemble step shows the original contigs. Press Update on Annotation Review to recheck the",
      "unlocked ones."), n),
    if (length(bad)) paste0("Converted with problems: ", paste(bad, collapse = "; "), "."),
    if (length(skip)) paste0("Not converted, left as they were: ", paste(skip, collapse = "; "), "."),
    if (length(done)) paste0("Locked contigs marked done without a recheck: ", paste(done, collapse = ", "), "."),
    if (length(left)) paste0("Trimmed contigs left as they are (dropped or never annotated): ",
      paste(left, collapse = ", "), "."))
  if (length(txt)) message(paste(txt, collapse = "\n"))
  invisible()
}

#' Convert one edited contig of an older project; returns its problems
#'
#' The copy takes the contig's current (trimmed or edited) state, with the
#' partial ends of edited features kept as the user's choice. A trimmed
#' contig's original is written back from its backup row and curated again;
#' when that is not possible the original stays as it is and the copy keeps
#' its full length.
#' @return character, empty when the contig converted cleanly.
#' @noRd
review_upgrade_unit <- function(con, u) {
  id <- u$ID; path <- as.integer(u$path); scaffold <- as.integer(u$scaffold)
  key <- list(id, path, scaffold)
  review_copy(con, id, path, scaffold)
  DBI::dbExecute(con, "UPDATE review_units SET edited = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
    params = key)
  DBI::dbExecute(con,
    "UPDATE review_annotations SET partial_start_user = partial_start, partial_stop_user = partial_stop
     WHERE ID = ? AND path = ? AND scaffold = ? AND edited = 1", params = key)
  probs <- character()
  if (u$has_bak == 1 && !identical(u$cur, u$bak)) {
    off <- if (anyNA(c(u$cur, u$bak)) || !nzchar(u$cur)) -1L else
      as.integer(regexpr(u$cur, u$bak, fixed = TRUE))
    if (is.na(u$bak)) {
      probs <- "original not restored: the saved original has no sequence"
    } else if (off < 1L) {
      probs <- "original not restored: the trimmed sequence is not in the saved original"
    } else {
      r <- tryCatch(with_savepoint(con, review_upgrade_original(con, u, off), name = "rp_upgrade_orig"),
        error = function(e) e)
      probs <- if (inherits(r, "error")) paste("original not restored:", conditionMessage(r)) else c(r,
        tryCatch({
          curate_one(con, id, path, scaffold)
          NULL
        }, error = function(e) paste("curation not rerun on the original:", conditionMessage(e))))
    }
  }
  if (u$locked == 1) review_upgrade_done(con, id, path, scaffold) else review_invalidate(con, id, path, scaffold)
  probs
}

#' Write a trimmed contig's original back from its backup row (schema step 4)
#'
#' off is where the trimmed sequence starts in the saved original; it also sets
#' the copy's window. The stored ITSx value may describe the trimmed sequence,
#' so it is cleared.
#' @return a problem when the saved chimera result could not be put back and
#'   the rows were shifted back instead, else character(0).
#' @noRd
review_upgrade_original <- function(con, u, off) {
  id <- u$ID; path <- as.integer(u$path); scaffold <- as.integer(u$scaffold)
  key <- list(id, path, scaffold)
  ann <- if (!is.na(u$bak_ann)) jsonlite::fromJSON(u$bak_ann)
  if (!is.data.frame(ann) && !identical(ann, list())) stop("the saved features cannot be read", call. = FALSE)
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, length = ? WHERE ID = ? AND path = ? AND scaffold = ?",
    params = c(list(u$bak, nchar(u$bak)), key))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = key)
  if (is.data.frame(ann) && nrow(ann)) {
    ann$ID <- id; ann$path <- path; ann$scaffold <- scaffold
    DBI::dbAppendTable(con, "annotations", ann[intersect(names(ann), DBI::dbListFields(con, "annotations"))])
  }
  write_coverage_summary(con, id, path, scaffold, parse_depths(u$bak_cov))
  DBI::dbExecute(con, "UPDATE annotate SET itsx_warnings = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
    params = key)
  DBI::dbExecute(con, "UPDATE review_units SET trim_start = ?, trim_end = ? WHERE ID = ? AND path = ? AND scaffold = ?",
    params = c(list(off, off + nchar(u$cur) - 1L), key))
  shift_back <- function() if (off > 1L) chimera_shift(con, id, path, scaffold, 1L - off, nchar(u$bak))
  if (is.na(u$bak_chim)) {
    shift_back()
    return(character())
  }
  tryCatch({
    chimera_restore(con, id, path, scaffold, u$bak_chim)
    character()
  }, error = function(e) {
    shift_back()
    paste("chimera result shifted back instead of restored:", conditionMessage(e))
  })
}

#' Mark a converted copy of a locked contig done without a recheck (schema step 4)
#' @noRd
review_upgrade_done <- function(con, id, path, scaffold) {
  cnt <- review_counts(con, id, path, scaffold)
  DBI::dbExecute(con,
    "UPDATE review_units SET check_switch = 2, n_features = ?, missing = ?,
       itsx_warnings = COALESCE(itsx_warnings, 'not_run(not checked before upgrade)'), time_stamp = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(cnt$n_features, cnt$missing, as.integer(Sys.time()), id, as.integer(path),
      as.integer(scaffold)))
}
