# Local BLAST against the shipped rDNA panel, composed with NCBI decoy sets so a
# non-metazoan contaminant is identified rather than merely unmatched.

#' NCBI prebuilt databases searched alongside the panel.
#'
#' The panel is metazoan nuclear rDNA, so anything else returns nothing and the
#' screen loses the class it would have reported. These restore it. They are
#' useless as a reference panel (the eukaryote sets are years stale) but a decoy
#' identifier only has to resolve to a kingdom.
#' @noRd
RP_DECOY_DBS <- c("16S_ribosomal_RNA", "SSU_eukaryote_rRNA", "LSU_eukaryote_rRNA")

#' TRUE when a BLAST database prefix has its volume files on disk
#' @noRd
blast_db_exists <- function(prefix) {
  if (!(length(prefix) == 1 && !is.na(prefix) && nzchar(prefix))) {
    return(FALSE)
  }
  index_file <- Sys.glob(paste0(prefix, c(".nin", ".nal")))
  length(index_file) > 0
}

#' Compose the -db argument for a local search
#'
#' Every database in a multi-database list must be built with -parse_seqids;
#' mixing one that carries accession info with one that does not fails the
#' search outright. blastn cannot express a whitespace-containing prefix in a
#' multi-database list, so whitespace in the panel is rejected; whitespace in a
#' decoy is skipped with a warning.
#'
#' @param panel Path prefix of the RiboPilot panel.
#' @param decoy_dir Directory holding the NCBI decoy databases, or NULL.
#' @return A single space-separated string, or NULL when the panel is unusable.
#' @noRd
local_blast_db_arg <- function(panel, decoy_dir = NULL) {
  if (!blast_db_exists(panel)) return(NULL)
  if (grepl("\\s", panel)) {
    warning("panel prefix contains whitespace, unusable: ", panel, call. = FALSE)
    return(NULL)
  }
  out <- panel
  if (!is.null(decoy_dir) && length(decoy_dir) == 1 && !is.na(decoy_dir) &&
      nzchar(decoy_dir)) {
    for (n in RP_DECOY_DBS) {
      p <- file.path(decoy_dir, n)
      if (grepl("\\s", p)) {
        warning("decoy path contains whitespace, skipping: ", p, call. = FALSE)
        next
      }
      if (blast_db_exists(p)) {
        out <- c(out, p)
      } else {
        warning("decoy database not found, skipping: ", p, call. = FALSE)
      }
    }
  }
  paste(out, collapse = " ")
}

#' Files a panel must have before blastn is pointed at it
#'
#' blast_db_exists() is deliberately looser (it gates a user-supplied prefix).
#' A directory holding only .nin passes that and then segfaults blastn, so
#' discovery, which chooses on the user's behalf, demands the whole set.
#' @noRd
PANEL_REQUIRED_EXT <- c("nin", "nsq", "nhr", "ndb")

#' TRUE when every required volume file is present and non-empty
#' @noRd
panel_complete <- function(prefix) {
  f <- paste0(prefix, ".", PANEL_REQUIRED_EXT)
  all(file.exists(f)) && all(file.info(f)$size > 0)
}

#' Newest usable panel installed under `db_dir`, or NULL
#'
#' Container-side only. The result must never be written to
#' reference_config.blastdb: an in-image path there would make the mount
#' builder bind a host directory over the panel and hide it.
#' @noRd
discover_panel <- function(db_dir = ribopilot_db_dir()) {
  if (length(db_dir) != 1 || is.na(db_dir) || !nzchar(db_dir)) return(NULL)
  db_dir <- path.expand(db_dir)
  if (!dir.exists(db_dir) || grepl("\\s", db_dir)) return(NULL)
  subs <- list.dirs(db_dir, full.names = FALSE, recursive = FALSE)
  subs <- subs[nzchar(subs) & !startsWith(subs, ".")]
  if (!length(subs)) return(NULL)
  v <- suppressWarnings(numeric_version(subs, strict = FALSE))
  keep <- !is.na(v)
  if (!any(keep)) return(NULL)
  subs <- subs[keep]
  subs <- subs[order(v[keep], decreasing = TRUE)]
  for (s in subs) {
    p <- file.path(db_dir, s, "ribopilot_rrna")
    if (panel_complete(p)) return(p)
  }
  NULL
}
