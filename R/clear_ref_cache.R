# What a reference folder caches: the clade databases and their index, the
# GenBank record counts, and the store of fetched reference sequences. The
# BLAST panel and container images that may sit beside them are left alone.
REF_CACHE_ITEMS <- c("manifest.tsv", "counts.tsv", "seeds", "labels", "refstore")

# Bound here so tests can stand in for the person at the keyboard.
interactive <- NULL
readline <- NULL

#' Cached items present in a reference folder
#' @return Character vector of paths.
#' @noRd
ref_cache_found <- function(ref_dir) {
  p <- file.path(ref_dir, REF_CACHE_ITEMS)
  p[file.exists(p)]
}

#' Clear the shared reference cache
#'
#' Deletes the cached reference data in a reference folder: every clade
#' database (seeds and labels), the list of built databases, the GenBank record
#' counts, and the stored reference sequences. The BLAST panel and container
#' images are not touched.
#'
#' The reference folder is shared. Every RiboPilot project that points at it,
#' yours or another user's, loses these files and must build its reference
#' databases again in the Reference step. A run in progress in any of those
#' projects can fail. Nothing can be restored afterwards.
#'
#' To give one project its own references instead, leave the shared folder
#' alone and point that project at a different folder: open the settings in the
#' app's Reference step, or pass `ref_dir` to [new_project()].
#'
#' The function lists what it found, then asks you to type `clear`. It refuses
#' to run outside an interactive R session, so a script cannot clear the cache
#' unattended.
#'
#' @param ref_dir Reference folder to clear (default: the `RIBOPILOT_REF_DIR`
#'   environment variable, else the per-user default).
#' @return (invisibly) the paths deleted; none when cancelled.
#' @export
clear_ref_cache <- function(ref_dir = ribopilot_ref_dir()) {
  if (!ref_dir_ok(ref_dir)) stop("Give a reference folder path.", call. = FALSE)
  ref_dir <- path.expand(trimws(ref_dir))
  found <- ref_cache_found(ref_dir)
  if (length(found) == 0) {
    message("Nothing to clear in ", ref_dir, ".")
    return(invisible(character(0)))
  }
  if (unname(file.access(ref_dir, mode = 2)) != 0) {
    stop(ref_dir, " is not writable, so its cache cannot be cleared.", call. = FALSE)
  }
  files <- c(found[!dir.exists(found)], list.files(found[dir.exists(found)],
    recursive = TRUE, full.names = TRUE))
  n_db <- nrow(read_ref_manifest(ref_dir))
  message(
    "WARNING: this deletes the shared RiboPilot reference cache.\n\n",
    "  Folder:  ", ref_dir, "\n",
    "  Holds:   ", n_db, " reference database(s), ", length(files), " file(s), ",
    format(structure(sum(file.size(files), na.rm = TRUE), class = "object_size"), units = "auto"), "\n",
    "  Deletes: ", paste(basename(found), collapse = ", "), "\n\n",
    "This affects EVERY RiboPilot project that uses this folder, not just one,\n",
    "including other users' projects on a shared drive. Each will have to build\n",
    "its reference databases again, and a run in progress may fail.\n",
    "This cannot be undone.\n\n",
    "Only want fresh references for one project? Do not clear this folder.\n",
    "Point that project at a new, empty reference folder instead: in the app,\n",
    "Reference step > settings; or new_project(ref_dir = \"<new folder>\").\n")
  if (!interactive()) {
    stop("clear_ref_cache() must be confirmed by a person. Run it in an interactive R session.",
      call. = FALSE)
  }
  if (!identical(trimws(readline("Type clear to delete, or anything else to cancel: ")), "clear")) {
    message("Cancelled. Nothing was deleted.")
    return(invisible(character(0)))
  }
  unlink(found, recursive = TRUE)
  left <- found[file.exists(found)]
  if (length(left)) warning("Could not delete: ", paste(left, collapse = ", "), call. = FALSE)
  message("Cleared ", ref_dir, ".")
  invisible(setdiff(found, left))
}
