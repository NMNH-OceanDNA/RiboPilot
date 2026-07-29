# Intragenomic rDNA variants: assembly-graph isomers versus genuine paralogs.
#
# A sample routinely yields several assembly units and the curator is asked to
# adjudicate between them with no information about why they differ. Two very
# different situations look identical: units that are the same molecule
# assembled two ways, differing only in how many IGS subrepeats the graph
# traversed, and units that are genuinely different rDNA copies from one
# individual. The first is a decision that does not matter; the second
# determines what GenBank receives.
#
# Only units the screen calls on_target are compared, so any difference between
# them is intragenomic by construction and no threshold has to separate a
# paralog from a contaminant (doc 12 section 1). The finding is ORTHOGONAL to
# the verdict: a paralogous variant is emphatically on_target.
#
# Only the ANCHORS are compared. The spacers are what varies in both cases and
# so cannot discriminate between them. RiboPilot's feature model has no IGS
# element, so the isomer test is "anchors identical, unannotated remainder
# differs in length" rather than a measurement of subrepeat copy number.

#' The anchors, in transcription order. Never coordinate order: a circularly
#' permuted assembly would otherwise yield a different string for one molecule.
#' @noRd
SCREEN_ANCHOR_ELEMENTS <- c("18S", "5.8S", "28S")

#' One unit's annotated anchors concatenated into a single comparable string
#'
#' A unit missing all three anchors is not comparable and returns NA, which is
#' distinct from being a `singleton`: NA says nothing was measured.
#'
#' @param seq The unit's `assemblies.sequence`.
#' @param ann Its `annotations` rows, with `element`, `start`, `stop`, `strand`.
#' @return character(1), or `NA_character_` when the unit is not comparable.
#' @noRd
screen_anchor_string <- function(seq, ann) {
  if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) return(NA_character_)
  if (is.null(ann) || nrow(ann) == 0 || !"element" %in% names(ann)) return(NA_character_)
  a <- ann[ann$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (nrow(a) == 0) return(NA_character_)
  ok <- !is.na(a$start) & !is.na(a$stop) & a$stop >= a$start
  a <- a[ok, , drop = FALSE]
  if (nrow(a) == 0) return(NA_character_)

  strand <- if ("strand" %in% names(a)) a$strand else rep(NA_character_, nrow(a))
  strand <- strand[!is.na(strand) & nzchar(strand)]
  # concatenating a sense and an antisense span produces a string that means
  # nothing, and a mixed-strand unit is a chimera signal rather than a paralog
  if (length(strand) && length(unique(strand)) > 1) return(NA_character_)
  minus <- length(strand) > 0 && identical(unique(strand), "-")

  ord <- order(match(a$element, SCREEN_ANCHOR_ELEMENTS))
  a <- a[ord, , drop = FALSE]
  parts <- substr(rep(seq[1], nrow(a)), a$start, a$stop)
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NA_character_)
  # each span reverse-complemented is that anchor in sense; the element sort
  # above already puts them in transcription order, so nothing is reversed twice
  if (minus) {
    parts <- vapply(parts, function(p)
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(p))),
      character(1), USE.NAMES = FALSE)
  }
  paste(parts, collapse = "")
}
