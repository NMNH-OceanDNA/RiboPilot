# Where built reference databases live, and what is already there.
#
# The folder is deliberately outside any project so one clade database serves
# every project (and, on a shared drive, every user) that needs that clade. The
# manifest is the only index of what a folder holds, so a fresh project can
# offer an existing database instead of downloading it again.

#' Project-level reference settings, with the folder resolved
#'
#' Resolution order: the project setting, then `RIBOPILOT_REF_DIR`, then
#' [ribopilot_ref_dir()].
#' @param con Open project DB connection.
#' @return list(ref_dir, min_records, max_records).
#' @noRd
ref_config <- function(con) {
  d <- tryCatch(DBI::dbGetQuery(con, "SELECT * FROM reference_config WHERE rowid_key = 1"),
    error = function(e) data.frame())
  stored <- if (nrow(d) > 0) d$ref_dir[1] else NA_character_
  dir <- if (!is.na(stored) && nzchar(trimws(stored))) trimws(stored) else ribopilot_ref_dir()
  list(
    ref_dir = path.expand(dir),
    min_records = if (nrow(d) > 0) as.integer(d$min_records[1] %|NA|% 1000L) else 1000L,
    max_records = if (nrow(d) > 0) as.integer(d$max_records[1] %|NA|% 25000L) else 25000L
  )
}

#' Write project-level reference settings (only the arguments supplied)
#' @param con Open project DB connection.
#' @param ref_dir,min_records,max_records New values, or NULL to leave as is.
#' @return (invisibly) the resulting [ref_config()].
#' @noRd
set_ref_config <- function(con, ref_dir = NULL, min_records = NULL, max_records = NULL) {
  if (!is.null(ref_dir)) {
    DBI::dbExecute(con, "UPDATE reference_config SET ref_dir = ? WHERE rowid_key = 1",
      params = list(if (nzchar(trimws(ref_dir))) path.expand(trimws(ref_dir)) else NA_character_))
  }
  if (!is.null(min_records)) {
    DBI::dbExecute(con, "UPDATE reference_config SET min_records = ? WHERE rowid_key = 1",
      params = list(as.integer(min_records)))
  }
  if (!is.null(max_records)) {
    DBI::dbExecute(con, "UPDATE reference_config SET max_records = ? WHERE rowid_key = 1",
      params = list(as.integer(max_records)))
  }
  invisible(ref_config(con))
}

#' Is `dir` usable as a reference folder? (created if missing)
#' @return list(ok, msg).
#' @noRd
validate_ref_dir <- function(dir) {
  dir <- path.expand(trimws(dir %||% ""))
  if (!nzchar(dir)) return(list(ok = FALSE, msg = "Give a folder path."))
  if (file.exists(dir) && !dir.exists(dir)) {
    return(list(ok = FALSE, msg = paste0(dir, " exists but is a file, not a folder.")))
  }
  if (!dir.exists(dir)) {
    ok <- dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) return(list(ok = FALSE, msg = paste0("Could not create ", dir, ".")))
  }
  if (unname(file.access(dir, mode = 2)) != 0) {
    return(list(ok = FALSE, msg = paste0(dir, " is not writable.")))
  }
  list(ok = TRUE, msg = dir)
}

#' Path to a reference folder's manifest
#' @noRd
ref_manifest_path <- function(dir) file.path(path.expand(dir), "manifest.tsv")

#' Manifest columns, in file order
#' @noRd
REF_MANIFEST_COLS <- c("clade", "taxid", "rank", "query", "n_records", "built", "seeds", "labels")

#' Read a reference folder's manifest
#' @return data.frame with [REF_MANIFEST_COLS]; zero rows when there is no manifest.
#' @noRd
read_ref_manifest <- function(dir) {
  empty <- data.frame(clade = character(0), taxid = integer(0), rank = character(0),
    query = character(0), n_records = integer(0), built = character(0),
    seeds = character(0), labels = character(0), stringsAsFactors = FALSE)
  p <- ref_manifest_path(dir)
  if (!file.exists(p)) return(empty)
  d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE, colClasses = "character"),
    error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0 || !all(REF_MANIFEST_COLS %in% names(d))) return(empty)
  d <- d[, REF_MANIFEST_COLS, drop = FALSE]
  d$taxid <- suppressWarnings(as.integer(d$taxid))
  d$n_records <- suppressWarnings(as.integer(d$n_records))
  d
}

#' Record a built database in the folder's manifest (one row per clade, newest wins)
#' @return (invisibly) the manifest path.
#' @noRd
append_ref_manifest <- function(dir, clade, taxid, rank, query, n_records, seeds, labels) {
  dir.create(path.expand(dir), recursive = TRUE, showWarnings = FALSE)
  m <- read_ref_manifest(dir)
  m <- m[m$clade != clade, , drop = FALSE]
  m <- rbind(m, data.frame(
    clade = clade, taxid = as.integer(taxid), rank = rank %||% NA_character_,
    query = query, n_records = as.integer(n_records),
    built = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    seeds = seeds, labels = labels, stringsAsFactors = FALSE))
  p <- ref_manifest_path(dir)
  utils::write.table(m, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  invisible(p)
}
