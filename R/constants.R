# Canonical nuclear ribosomal repeat (nrDNA) element names, shared across the
# package (annotation reconciliation, feature-table export, app selectors).
# The repeat is 5'-[ETS]-18S-ITS1-5.8S-ITS2-28S-[IGS]-3', tandemly repeated.

#' Standard nuclear rRNA genes (the conserved anchors)
#' @noRd
RRNA_GENES <- c("18S", "5.8S", "28S")

#' Internal transcribed spacers (the variable, boundary-defined products)
#' @noRd
RRNA_SPACERS <- c("ITS1", "ITS2")

#' Optional external spacers (usually not recovered/curated by default)
#' @noRd
RRNA_EXTRA <- c("ETS", "IGS")

#' Canonical 5'->3' order of elements in the cistron. Used to order features
#' and derive spacers by subtraction between conserved anchors.
#' @noRd
RRNA_ELEMENT_ORDER <- c("ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "IGS")

#' Flat list of element names for use as selectInput choices.
#' @noRd
RRNA_ELEMENT_CHOICES <- c(RRNA_GENES, RRNA_SPACERS, RRNA_EXTRA)

#' GenBank feature type for an nrDNA element name.
#'
#' rRNA genes -> "rRNA"; transcribed spacers -> "misc_RNA".
#'
#' @param element element name (e.g. "18S", "ITS1")
#' @noRd
rrna_element_type <- function(element) {
  if (element %in% RRNA_GENES) {
    "rRNA"
  } else {
    "misc_RNA"
  }
}

#' GenBank product / note string for an nrDNA element name.
#' @noRd
rrna_element_product <- function(element) {
  switch(element,
    "18S"  = "18S ribosomal RNA",
    "5.8S" = "5.8S ribosomal RNA",
    "28S"  = "28S ribosomal RNA",
    "ITS1" = "internal transcribed spacer 1",
    "ITS2" = "internal transcribed spacer 2",
    "ETS"  = "external transcribed spacer",
    "IGS"  = "intergenic spacer",
    element
  )
}
