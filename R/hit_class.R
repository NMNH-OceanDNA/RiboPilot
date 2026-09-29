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
#' non-nuclear subunit test comes next so that "16S rRNA gene" cannot match the
#' bare "rRNA gene" construction; then nuclear rDNA phrasing; then assembly-level
#' phrasing; and anything unrecognised stays "other", as it was before.
#'
#' @param title Character vector of defline titles.
#' @return Character vector of the same length, each one of `HIT_CLASSES`.
#' @noRd
classify_hit_title <- function(title) {
  t <- as.character(title)
  if (length(t) == 0) return(character(0))
  t[is.na(t)] <- ""
  out <- rep("other", length(t))
  is_rdna <- grepl(NUCLEAR_RRNA_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_non_nuc <- grepl(NON_NUCLEAR_RRNA_RE, t, ignore.case = TRUE, perl = TRUE) &
    !grepl(NUCLEAR_RRNA_MARKER_RE, t, ignore.case = TRUE, perl = TRUE)
  is_asm <- grepl(ASSEMBLY_TITLE_RE, t, ignore.case = TRUE, perl = TRUE)
  is_org <- grepl(ORGANELLAR_TITLE_RE, t, ignore.case = TRUE)
  # assigned low precedence first, so each later line overrides the earlier
  out[is_asm] <- "uninformative"
  out[is_rdna] <- "rdna"
  out[is_non_nuc] <- "other"
  out[is_org] <- "other"
  out[!nzchar(t)] <- "uninformative"
  out
}

#' Decide a unit's nuclear-rDNA flag from its hit distribution
#'
#' Votes over hits within 10 percent of the top bitscore, taking at least three
#' where that many exist, which is the same defence against one mislabelled
#' GenBank accession that [bitscore_weighted_lca()] uses. Uninformative hits are
#' dropped from the vote rather than counted against rDNA.
#'
#' @param hit_class Character vector of `HIT_CLASSES`.
#' @param bitscore Numeric vector, same length.
#' @return `1L` (nuclear rDNA), `0L` (not), or `NA_integer_` when no hit is
#'   informative.
#' @noRd
consensus_hit_class <- function(hit_class, bitscore) {
  # parallel vectors: a mismatch either crashes on a poisoned max() or silently
  # votes on misaligned rows, so fail loudly instead
  stopifnot(length(hit_class) == length(bitscore))
  hc <- as.character(hit_class)
  bs <- suppressWarnings(as.numeric(bitscore))
  ok <- !is.na(hc) & !is.na(bs)
  hc <- hc[ok]; bs <- bs[ok]
  if (length(bs) == 0) return(NA_integer_)
  sel <- bs >= 0.9 * max(bs)
  if (sum(sel) < 3L) {
    sel <- seq_along(bs) %in% order(bs, decreasing = TRUE)[seq_len(min(3L, length(bs)))]
  }
  hc <- hc[sel]; bs <- bs[sel]
  keep <- hc %in% c("rdna", "other")
  if (!any(keep)) return(NA_integer_)
  hc <- hc[keep]; bs <- bs[keep]
  # tie resolves to rdna: a false 0 raises a spurious triage flag, a false 1
  # only fails to raise one
  if (sum(bs[hc == "rdna"]) >= sum(bs[hc == "other"])) 1L else 0L
}
