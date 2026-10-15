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

#' Align the contigs of one set
#' @param seqs Named character vector (unit keys as names).
#' @return Aligned DNAStringSet or NULL.
#' @noRd
paths_align <- function(seqs) {
  seqs <- seqs[!is.na(seqs) & nzchar(seqs)]
  if (length(seqs) < 2) return(NULL)
  # ponytail: no cap on contig count; DECIPHER on a few 10 kb contigs takes seconds
  dna <- Biostrings::DNAStringSet(toupper(seqs))
  tryCatch(DECIPHER::AlignSeqs(dna, verbose = FALSE), error = function(e) NULL)
}

#' Viewer payload for an alignment: one row per contig, features in columns
#' @noRd
paths_payload <- function(aln, feats, labels) {
  seqs <- as.character(aln)
  rows <- lapply(names(seqs), function(u) {
    s <- seqs[[u]]
    cols <- which(strsplit(s, "", fixed = TRUE)[[1]] != "-")
    f <- feats[feats$unit == u & !is.na(feats$start) & !is.na(feats$stop), , drop = FALSE]
    list(unit = u, label = unname(labels[u]) %|NA|% u, seq = s,
      features = lapply(seq_len(nrow(f)), function(i) list(
        row = f$element[i], gene = f$element[i],
        pos1 = as.integer(cols[f$start[i]]), pos2 = as.integer(cols[f$stop[i]]),
        dir = if (identical(f$strand[i], "-")) "-" else "+",
        partial5 = isTRUE(f$partial_start[i] == 1), partial3 = isTRUE(f$partial_stop[i] == 1),
        color = unname(RRNA_ELEMENT_COLORS[f$element[i]] %|NA|% "#888888"))))
  })
  list(len = as.integer(nchar(seqs[1])), rows = rows)
}

#' Columns where the rows disagree (gaps count)
#' @noRd
paths_diff_count <- function(rows) {
  if (length(rows) < 2) return(0L)
  m <- do.call(rbind, lapply(rows, function(r) strsplit(r$seq, "", fixed = TRUE)[[1]]))
  sum(apply(m, 2, function(col) length(unique(col)) > 1))
}
