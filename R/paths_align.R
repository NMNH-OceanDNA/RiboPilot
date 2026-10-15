# Paths alignment panel: homologous contigs across a sample's assembly paths,
# aligned with DECIPHER and shown in the sequence viewer's alignment mode.

#' Homologous sets of contigs: connected components of "shares an anchor"
#' @param feats data.frame(path, scaffold, element), one row per annotation.
#' @return list of list(label, units = data.frame(path, scaffold)); only sets
#'   of two or more contigs, largest first.
#' @noRd
paths_homology_sets <- function(feats) {
  f <- feats[feats$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (!nrow(f)) return(list())
  f$unit <- paste(f$path, f$scaffold)
  units <- unique(f$unit)
  comp <- stats::setNames(seq_along(units), units)
  for (el in unique(f$element)) {
    members <- unique(f$unit[f$element == el])
    comp[comp %in% comp[members]] <- min(comp[members])
  }
  sets <- split(names(comp), comp)
  sets <- Filter(function(u) length(u) >= 2, sets)
  if (!length(sets)) return(list())
  out <- lapply(sets, function(u) {
    pu <- do.call(rbind, strsplit(u, " ", fixed = TRUE))
    units <- data.frame(path = as.integer(pu[, 1]), scaffold = as.integer(pu[, 2]))
    units <- units[order(units$path, units$scaffold), , drop = FALSE]
    rownames(units) <- NULL
    anchors <- SCREEN_ANCHOR_ELEMENTS[SCREEN_ANCHOR_ELEMENTS %in% f$element[f$unit %in% u]]
    list(label = sprintf("%s (%d contigs)", paste(anchors, collapse = "-"), nrow(units)), units = units)
  })
  n <- vapply(out, function(s) nrow(s$units), 1L)
  first <- vapply(out, function(s) s$units$path[1] * 1e6 + s$units$scaffold[1], 1)
  unname(out[order(-n, first)])
}
