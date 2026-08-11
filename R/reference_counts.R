# How many nuclear rRNA records GenBank holds for a clade.
#
# This is the evidence the clade suggestion walks on: a clade nobody has
# sequenced makes a useless seed database, and one covering half of Metazoa
# makes an unspecific one. The count comes from the same query the builder will
# run, so the number shown is the number that will be downloaded.

#' GenBank organism term for a clade
#'
#' Prefers the numeric taxid when one is supplied: `txid<ID>[Organism:exp]`
#' resolves against a stable NCBI node, so a clade renamed since the project's
#' pinned taxdump still returns its records. Falls back to the name for callers
#' that only have one.
#' @param name Clade name.
#' @param taxid Optional NCBI taxid.
#' @noRd
clade_query <- function(name, taxid = NULL) {
  taxid <- suppressWarnings(as.integer(taxid %||% NA_integer_))
  if (length(taxid) == 1 && !is.na(taxid)) {
    return(sprintf("txid%d[Organism:exp]", taxid))
  }
  sprintf("\"%s\"[Organism]", trimws(name))
}

#' Pull the record count out of the builder script's output
#' @param lines Character vector of script stdout.
#' @return integer, or NA when the script did not report a count.
#' @noRd
parse_record_count <- function(lines) {
  hit <- grep("^Found [0-9]+ matching", lines, value = TRUE)
  if (length(hit) == 0) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Found ([0-9]+) .*$", "\\1", hit[1])))
}

#' Path to a reference folder's record-count cache
#' @noRd
ref_counts_path <- function(dir) file.path(path.expand(dir), "counts.tsv")

#' Read a reference folder's record-count cache
#' @noRd
read_counts_cache <- function(dir) {
  empty <- data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE)
  p <- ref_counts_path(dir)
  if (!file.exists(p)) return(empty)
  d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d) || !all(c("taxid", "n_records") %in% names(d))) return(empty)
  d$taxid <- suppressWarnings(as.integer(d$taxid))
  d$n_records <- suppressWarnings(as.integer(d$n_records))
  d[, c("taxid", "n_records"), drop = FALSE]
}

#' Cache-only count lookup: NA for anything not counted yet
#' @noRd
cache_count_fn <- function(dir) {
  function(taxid, name) {
    cache <- read_counts_cache(dir)
    i <- match(suppressWarnings(as.integer(taxid)), cache$taxid)
    if (is.na(i)) NA_integer_ else cache$n_records[i]
  }
}
