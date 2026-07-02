# WF1 reference alignment: per-feature nucleotide alignment of a sample's
# rRNA/ITS elements against its top BLAST reference, for the app's alignment
# view. Pure R (Biostrings/pwalign); no external aligner needed.

#' Global pairwise nucleotide alignment wrapper (pwalign, fallback Biostrings)
#' @noRd
pairwise_nt <- function(a, b) {
  fn <- if (requireNamespace("pwalign", quietly = TRUE)) {
    pwalign::pairwiseAlignment
  } else {
    getExportedValue("Biostrings", "pairwiseAlignment")
  }
  fn(Biostrings::DNAString(a), Biostrings::DNAString(b), type = "global")
}

#' Compute per-feature sample-vs-reference alignments for one sample
#'
#' Aligns each rRNA/ITS element that both the sample and its top reference
#' annotate, writing aligned strings + percent identity to `blast_ref_alignment`.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return (invisibly) number of element alignments written.
#' @export
compute_ref_alignment <- function(con, id) {
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM blast_ref_alignment WHERE ID = ?", params = list(id))

  acc <- app_query_con(con, "SELECT reference_accession FROM stitched WHERE ID = ?", id)$reference_accession
  if (length(acc) == 0 || is.na(acc[1]) || !nzchar(acc[1])) return(invisible(0L))
  acc <- acc[1]

  s_seq <- app_query_con(con, "SELECT sequence FROM stitched WHERE ID = ?", id)$sequence
  r_seq <- app_query_con(con, "SELECT sequence FROM blast_ref_sequences WHERE accession = ?", acc)$sequence
  if (length(s_seq) == 0 || length(r_seq) == 0 || !nzchar(s_seq[1]) || !nzchar(r_seq[1])) {
    return(invisible(0L))
  }
  s_seq <- s_seq[1]; r_seq <- r_seq[1]

  s_ann <- app_query_con(con, "SELECT element, start, stop FROM annotations WHERE ID = ?", id)
  r_ann <- app_query_con(con, "SELECT element, start, stop, strand FROM blast_ref_annotations WHERE accession = ?", acc)
  if (nrow(s_ann) == 0 || nrow(r_ann) == 0) return(invisible(0L))

  shared <- intersect(s_ann$element, r_ann$element)
  n <- 0L
  for (el in shared) {
    sa <- s_ann[s_ann$element == el, ][1, ]
    ra <- r_ann[r_ann$element == el, ][1, ]
    s_sub <- substr(s_seq, sa$start, sa$stop)
    r_sub <- substr(r_seq, ra$start, ra$stop)
    if (identical(ra$strand, "-")) {
      r_sub <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r_sub)))
    }
    if (!nzchar(s_sub) || !nzchar(r_sub)) next
    aln <- tryCatch(pairwise_nt(s_sub, r_sub), error = function(e) NULL)
    if (is.null(aln)) next
    s_aln <- as.character(pwalign_pattern(aln))
    r_aln <- as.character(pwalign_subject(aln))
    pid <- pwalign_pid(aln)
    DBI::dbAppendTable(con, "blast_ref_alignment", data.frame(ID = id, accession = acc,
      element = el, pid = round(pid, 1), aln_len = nchar(s_aln),
      sample_aln = s_aln, ref_aln = r_aln, time_stamp = ts, stringsAsFactors = FALSE))
    n <- n + 1L
  }
  invisible(n)
}

# Accessors that work whether pairwiseAlignment came from pwalign or Biostrings.
#' @noRd
pwalign_pattern <- function(aln) {
  fn <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::pattern else getExportedValue("Biostrings", "pattern")
  fn(aln)
}
#' @noRd
pwalign_subject <- function(aln) {
  fn <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::subject else getExportedValue("Biostrings", "subject")
  fn(aln)
}
#' @noRd
pwalign_pid <- function(aln) {
  fn <- if (requireNamespace("pwalign", quietly = TRUE)) pwalign::pid else getExportedValue("Biostrings", "pid")
  fn(aln)
}

#' Query helper bound to an existing connection (parameterized)
#' @noRd
app_query_con <- function(con, sql, ...) {
  DBI::dbGetQuery(con, sql, params = list(...))
}
