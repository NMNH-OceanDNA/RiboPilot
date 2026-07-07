# WF1 reference alignment: per-feature nucleotide alignment of a sample's
# rRNA/ITS elements against its top BLAST reference, for the app's alignment
# view. Pure R (Biostrings/pwalign); no external aligner needed.

#' Live per-feature alignment of a sample's element against its reference,
#' computed from the element's CURRENT boundaries (so the app's alignment view
#' redraws as edits move the boundaries). Aligns the two short, homologous
#' subsequences (reliable global alignment). Returns
#' list(accession, pid, sample_aln, ref_aln) or NULL.
#' @noRd
compute_feature_alignment <- function(con, id, element) {
  acc <- app_query_con(con, "SELECT reference_accession FROM stitched WHERE ID = ?", id)$reference_accession
  if (length(acc) == 0 || is.na(acc[1]) || !nzchar(acc[1])) return(NULL)
  acc <- acc[1]
  s_seq <- app_query_con(con, "SELECT sequence FROM stitched WHERE ID = ?", id)$sequence
  r_seq <- app_query_con(con, "SELECT sequence FROM blast_ref_sequences WHERE accession = ?", acc)$sequence
  sa <- app_query_con(con, "SELECT start, stop FROM annotations WHERE ID = ? AND element = ?", id, element)
  ra <- app_query_con(con, "SELECT start, stop, strand FROM blast_ref_annotations WHERE accession = ? AND element = ?", acc, element)
  if (length(s_seq) == 0 || length(r_seq) == 0 || nrow(sa) == 0 || nrow(ra) == 0) return(NULL)
  s_sub <- substr(s_seq[1], sa$start[1], sa$stop[1])
  r_sub <- substr(r_seq[1], ra$start[1], ra$stop[1])
  if (identical(ra$strand[1], "-")) {
    r_sub <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(r_sub)))
  }
  if (!nzchar(s_sub) || !nzchar(r_sub)) return(NULL)
  aln <- tryCatch(pairwise_nt(s_sub, r_sub), error = function(e) NULL)
  if (is.null(aln)) return(NULL)
  list(accession = acc, pid = pwalign_pid(aln),
    sample_aln = as.character(pwalign_pattern(aln)),
    ref_aln = as.character(pwalign_subject(aln)))
}

#' Multiple alignment of a sample's element against ALL of its BLAST-hit
#' references (the ranked blast_ref_candidates that have a fetched sequence + a
#' parsed feature for the element). Uses the sample's CURRENT boundaries, so the
#' view redraws on edits. Returns an aligned DNAStringSet (sample first) or NULL.
#' @noRd
compute_feature_msa <- function(con, id, element, max_hits = Inf) {
  s_seq <- app_query_con(con, "SELECT sequence FROM stitched WHERE ID = ?", id)$sequence
  sa <- app_query_con(con, "SELECT start, stop FROM annotations WHERE ID = ? AND element = ?", id, element)
  if (length(s_seq) == 0 || nrow(sa) == 0) return(NULL)
  s_sub <- substr(s_seq[1], sa$start[1], sa$stop[1])
  if (!nzchar(s_sub)) return(NULL)

  cands <- app_query_con(con,
    "SELECT accession, species FROM blast_ref_candidates WHERE ID = ? ORDER BY rank", id)
  seqs <- c(s_sub); names(seqs) <- paste0("sample|", element)
  for (i in seq_len(nrow(cands))) {
    if (length(seqs) > max_hits) break
    acc <- cands$accession[i]
    rseq <- app_query_con(con, "SELECT sequence FROM blast_ref_sequences WHERE accession = ?", acc)$sequence
    ra <- app_query_con(con, "SELECT start, stop, strand FROM blast_ref_annotations WHERE accession = ? AND element = ?", acc, element)
    if (length(rseq) == 0 || nrow(ra) == 0) next
    sub <- substr(rseq[1], ra$start[1], ra$stop[1])
    if (identical(ra$strand[1], "-")) sub <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(sub)))
    if (!nzchar(sub)) next
    nm <- if (!is.na(cands$species[i]) && nzchar(cands$species[i])) paste0(acc, " ", cands$species[i]) else acc
    seqs[nm] <- sub
  }
  if (length(seqs) < 2) return(NULL)   # need the sample + at least one reference
  dna <- Biostrings::DNAStringSet(seqs)
  tryCatch(DECIPHER::AlignSeqs(dna, verbose = FALSE), error = function(e) NULL)
}

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
