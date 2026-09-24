# WF1 reference fetch: store the NCBI records blast_ref.nf fetched in batch, so the
# app can align each sample feature against its reference.

#' Split fasta deflines ">ACC.V description" into accession and title
#' @param full Character vector of defline strings (Biostrings names()).
#' @return list(acc, title)
#' @noRd
parse_defline_ids <- function(full) {
  list(acc = sub("\\s.*$", "", full), title = trimws(sub("^\\S+\\s*", "", full)))
}

#' TRUE if the file's first non-blank line is a FASTA header
#' @noRd
is_fasta_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) return(FALSE)
  ls <- readLines(path, warn = FALSE)
  ls <- ls[nzchar(trimws(ls))]
  length(ls) > 0 && grepl("^>", ls[1])
}

#' Store every record of a multi-record fasta as a reference
#'
#' Sequence comes from the fasta; element coordinates come from `ann` rather than
#' from an NCBI feature table. An accession with no rows in `ann` is still stored,
#' with no annotations.
#'
#' @param con Open project DB connection.
#' @param fasta_file Multi-record fasta, deflines `>ACC.V description`.
#' @param ann Annotation frame carrying seqid, element, type, start, stop, strand,
#'   or NULL.
#' @return (invisibly) the number of accessions stored.
#' @export
write_blast_refs <- function(con, fasta_file, ann = NULL) {
  if (!is_fasta_file(fasta_file)) return(invisible(0L))
  seqs <- Biostrings::readDNAStringSet(fasta_file)
  if (length(seqs) == 0) return(invisible(0L))
  ts <- as.integer(Sys.time())
  ids <- parse_defline_ids(names(seqs))
  acc <- ids$acc
  title <- ids$title
  n <- 0L
  for (i in seq_along(seqs)) {
    s <- as.character(seqs[[i]])
    if (!nzchar(s)) next
    a <- ann
    rows <- if (!is.null(a) && nrow(a) > 0) a[a$seqid == acc[i], , drop = FALSE] else NULL
    DBI::dbWithTransaction(con, {
      DBI::dbExecute(con, "DELETE FROM blast_ref_sequences WHERE accession = ?",
        params = list(acc[i]))
      DBI::dbAppendTable(con, "blast_ref_sequences", data.frame(
        accession = acc[i], species = extract_species(title[i]), title = title[i],
        taxid = NA_character_, lineage = NA_character_, sequence = s,
        length = nchar(s), time_stamp = ts, stringsAsFactors = FALSE))
      DBI::dbExecute(con, "DELETE FROM blast_ref_annotations WHERE accession = ?",
        params = list(acc[i]))
      if (!is.null(rows) && nrow(rows) > 0) {
        DBI::dbAppendTable(con, "blast_ref_annotations", data.frame(
          accession = acc[i], element = rows$element, type = rows$type,
          start = rows$start, stop = rows$stop, strand = rows$strand,
          time_stamp = ts, stringsAsFactors = FALSE))
      }
    })
    n <- n + 1L
  }
  invisible(n)
}

#' Fail loudly when a fetched batch stored nothing
#'
#' Mirrors the old per-accession guard (`WANTED > 0 && STORED == 0`): a fasta
#' that held at least one input record but stored none means NCBI returned
#' bytes that then failed to store, which is worth stopping on even though
#' `refs.fa` itself was non-empty.
#'
#' @param fasta_file Path to the fasta passed to [write_blast_refs()]; named
#'   in the error and re-read here to count input records.
#' @param n_stored Number of records [write_blast_refs()] actually stored.
#' @return (invisibly) `n_stored`, when the check passes.
#' @export
stop_if_nothing_stored <- function(fasta_file, n_stored) {
  n_in <- if (file.exists(fasta_file)) {
    sum(grepl("^>", readLines(fasta_file, warn = FALSE)))
  } else {
    0L
  }
  if (n_in > 0 && n_stored == 0) {
    stop("stored nothing from ", fasta_file, " (", n_in, " input record(s))", call. = FALSE)
  }
  invisible(n_stored)
}
