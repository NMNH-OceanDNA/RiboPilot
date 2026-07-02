# Project-level taxonomy: resolve every sample's Taxon, let the user confirm ALL
# of them, then freeze the confirmed taxid + lineage into the project database.
#
# Two-step, user-in-the-loop (mirrors how the Shiny app will present a modal):
#   resolve_project_taxonomy(path)  -> writes taxonomy_review.csv for review
#   (user edits confirmed_taxid for any rows they want to change)
#   apply_project_taxonomy(path)    -> freezes confirmed values into `samples`

#' Resolve all sample Taxon strings and write a review table for confirmation
#'
#' Downloads (once) a pinned NCBI taxdump into `<path>/.taxonomy`, resolves every
#' unique `Taxon` in the project, and writes `taxonomy_review.csv` listing EVERY
#' sample's proposal (not only failures) for the user to confirm or correct. Edit
#' the `confirmed_taxid` column for any row, then run [apply_project_taxonomy()].
#'
#' @param path Project directory.
#' @param taxdump_dir Taxdump location (default `<path>/.taxonomy`).
#' @param force_taxdump Re-download the taxdump even if present.
#' @return (invisibly) the review data.frame.
#' @export
resolve_project_taxonomy <- function(path = ".",
                                     taxdump_dir = file.path(path, ".taxonomy"),
                                     force_taxdump = FALSE) {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) {
    stop("No project database at ", db, call. = FALSE)
  }
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  samples <- DBI::dbGetQuery(con, "SELECT ID, Taxon FROM samples")

  ensure_taxdump(taxdump_dir, force = force_taxdump)

  resolved <- resolve_taxonomy(unique(samples$Taxon), taxdump_dir = taxdump_dir)

  review <- data.frame(
    Taxon = resolved$Taxon,
    match_type = resolved$match_type,
    proposed_taxid = resolved$taxid,
    proposed_name = resolved$sci_name,
    rank = resolved$rank,
    lineage = resolved$lineage,
    alternatives = resolved$alternatives,
    confirmed_taxid = resolved$taxid, # user edits this column to override
    stringsAsFactors = FALSE
  )
  out_csv <- file.path(path, "taxonomy_review.csv")
  utils::write.csv(review, out_csv, row.names = FALSE, na = "")

  n_flag <- sum(review$match_type != "exact")
  message(
    "Resolved ", nrow(review), " unique taxa (", n_flag, " fuzzy/unresolved).\n",
    "Review and edit `confirmed_taxid` in:\n  ", out_csv, "\n",
    "Then run apply_project_taxonomy('", path, "') to freeze into the database."
  )
  invisible(review)
}

#' Freeze confirmed taxids + lineages from the review CSV into `samples`
#'
#' Reads `taxonomy_review.csv`, (re)fetches lineages for any `confirmed_taxid`
#' the user changed, and writes `taxid` + `lineage` into the `samples` table
#' joined on `Taxon`. Rows with no confirmed taxid are flagged (left NULL).
#'
#' @param path Project directory.
#' @param taxdump_dir Taxdump location (default `<path>/.taxonomy`).
#' @return (invisibly) the frozen per-Taxon data.frame.
#' @export
apply_project_taxonomy <- function(path = ".",
                                   taxdump_dir = file.path(path, ".taxonomy")) {
  csv <- file.path(path, "taxonomy_review.csv")
  if (!file.exists(csv)) {
    stop("No taxonomy_review.csv in ", path, ". Run resolve_project_taxonomy() first.", call. = FALSE)
  }
  review <- utils::read.csv(csv, stringsAsFactors = FALSE)
  review$confirmed_taxid <- suppressWarnings(as.integer(review$confirmed_taxid))

  # refresh lineage for taxids the user changed away from the proposal
  changed <- !is.na(review$confirmed_taxid) &
    (is.na(review$proposed_taxid) | review$confirmed_taxid != review$proposed_taxid)
  review$final_lineage <- review$lineage
  if (any(changed) && have_taxdump(taxdump_dir)) {
    dd <- if (file.exists(file.path(taxdump_dir, "names.dmp"))) taxdump_dir else extract_taxdump(taxdump_dir)
    if (!identical(dd, taxdump_dir)) on.exit(unlink(dd, recursive = TRUE), add = TRUE)
    info <- tk_taxid_info(review$confirmed_taxid[changed], dd)
    review$final_lineage[changed] <- info$lineage[match(review$confirmed_taxid[changed], info$taxid)]
  }

  db <- file.path(path, ".sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  upd <- DBI::dbSendStatement(con, "UPDATE samples SET taxid = ?, lineage = ? WHERE Taxon = ?")
  for (i in seq_len(nrow(review))) {
    DBI::dbBind(upd, list(
      review$confirmed_taxid[i],
      review$final_lineage[i] %|NA|% NA,
      review$Taxon[i]
    ))
  }
  DBI::dbClearResult(upd)

  n_ok <- sum(!is.na(review$confirmed_taxid))
  message("Froze taxid + lineage for ", n_ok, "/", nrow(review),
    " taxa into samples (", nrow(review) - n_ok, " left unresolved).")
  invisible(review)
}

#' Write scored candidate paths for a sample into the `path_scores` table
#'
#' @param con Open DBI connection to the project database.
#' @param ID Sample ID.
#' @param scored A scored candidate data.frame from [score_assembly()] /
#'   [score_candidates()].
#' @return (invisibly) NULL. Existing rows for `ID` are replaced.
#' @export
write_path_scores <- function(con, ID, scored) {
  if (is.null(scored) || nrow(scored) == 0) {
    return(invisible())
  }
  cols <- c("length", "depth", "anchors", "blast_pident", "hit_id", "hit_taxid",
    "hit_species", "blast_lineage", "s_tax", "s_identity", "s_anchors",
    "s_depth", "s_length", "score", "rank", "auto_select")
  df <- data.frame(ID = ID, path = scored$path, stringsAsFactors = FALSE)
  for (c in cols) df[[c]] <- if (!is.null(scored[[c]])) scored[[c]] else NA
  df$auto_select <- as.integer(df$auto_select)
  df$time_stamp <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM path_scores WHERE ID = ?", params = list(ID))
  DBI::dbAppendTable(con, "path_scores", df)
  invisible()
}
