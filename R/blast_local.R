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
    all_nin <- Sys.glob(file.path(decoy_dir, "*.nin"))
    all_nal <- Sys.glob(file.path(decoy_dir, "*.nal"))
    all_dbs <- c(all_nin, all_nal)
    all_prefixes <- unique(sub("\\.n(in|al)$", "", all_dbs))
    expected_prefixes <- file.path(decoy_dir, RP_DECOY_DBS)
    unknown_prefixes <- setdiff(all_prefixes, expected_prefixes)
    has_whitespace_decoy <- any(grepl("\\s", unknown_prefixes))
    if (has_whitespace_decoy) {
      for (p in unknown_prefixes) {
        if (grepl("\\s", p)) {
          warning("decoy path contains whitespace, skipping: ", p, call. = FALSE)
        }
      }
      return(paste(out, collapse = " "))
    }
    for (n in RP_DECOY_DBS) {
      p <- file.path(decoy_dir, n)
      if (grepl("\\s", p)) {
        warning("decoy path contains whitespace, skipping: ", p, call. = FALSE)
        next
      }
      if (blast_db_exists(p)) {
        out <- c(out, p)
      } else {
        warning("decoy database not found, skipping: ", n, call. = FALSE)
      }
    }
  }
  paste(out, collapse = " ")
}
