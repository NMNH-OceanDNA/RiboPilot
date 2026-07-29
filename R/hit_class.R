# Three-state classification of a GenBank defline.
#
# This answers "does this record describe nuclear rDNA". It does NOT answer "is
# this the host": a fungal ITS record and a human 18S record are both correctly
# rdna. Deciding that either is a contaminant is the taxonomy channel's job,
# from rank distance and the identity envelope.
#
# The third state exists because a chromosome-level or whole-genome record is
# not evidence AGAINST rDNA. Counting it as such is what made genuine host
# scaffolds read as "not a nuclear rDNA record".

#' The three hit classes, in no particular order.
#' @noRd
HIT_CLASSES <- c("rdna", "other", "uninformative")

#' Classify GenBank defline(s) into rdna / other / uninformative
#'
#' Precedence is fixed: the organellar veto wins over any rDNA phrase, because a
#' mitochondrial LSU record contains "large subunit ribosomal RNA"; the
#' prokaryotic test comes next so that "16S rRNA gene" cannot match the bare
#' "rRNA gene" construction; then nuclear rDNA phrasing; then assembly-level
#' phrasing; and anything unrecognised stays "other", as it was before.
#'
#' @param title Character vector of defline titles.
#' @return Character vector of the same length, each one of [HIT_CLASSES].
#' @noRd
classify_hit_title <- function(title) {
  t <- as.character(title)
  if (length(t) == 0) return(character(0))
  t[is.na(t)] <- ""
  out <- rep("other", length(t))
  is_rdna <- grepl(NUCLEAR_RRNA_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_prok <- grepl(PROKARYOTIC_RRNA_RE, t, ignore.case = TRUE, perl = TRUE) &
    !grepl(EUKARYOTIC_RRNA_MARKER_RE, t, ignore.case = TRUE, perl = TRUE)
  is_asm <- grepl(ASSEMBLY_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_org <- grepl(ORGANELLAR_TITLE_RE, t, ignore.case = TRUE)
  # assigned low precedence first, so each later line overrides the earlier
  out[is_asm] <- "uninformative"
  out[is_rdna] <- "rdna"
  out[is_prok] <- "other"
  out[is_org] <- "other"
  out[!nzchar(t)] <- "uninformative"
  out
}
