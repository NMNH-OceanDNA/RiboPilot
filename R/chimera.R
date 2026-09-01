# Chimera detection: per-gene taxonomy, windowed local BLAST, ribotyper,
# read evidence. Design: dev/design/21_CHIMERA_DETECTION.md.

CHIMERA_WIN <- 500L
CHIMERA_STEP <- 250L
CHIMERA_PID_DROP <- 90
CHIMERA_DEPTH_STEP <- 3
CHIMERA_CLIP_FRAC <- 0.3
CHIMERA_JUNCTION_TOL <- 300L
CHIMERA_TAXONOMY_CAP <- 0.7
CHIMERA_LABELS <- c("clean", "suspect", "likely")

#' @noRd
chimera_tile <- function(seq_len, win = CHIMERA_WIN, step = CHIMERA_STEP) {
  seq_len <- as.integer(seq_len)
  if (seq_len <= win) return(data.frame(start = 1L, stop = seq_len))
  starts <- seq.int(1L, seq_len - win + 1L, by = step)
  if (utils::tail(starts, 1) + win - 1L < seq_len) starts <- c(starts, seq_len - win + 1L)
  data.frame(start = as.integer(starts), stop = as.integer(starts + win - 1L))
}

#' @noRd
chimera_windows_fasta <- function(seq, tiles, file) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  subs <- substring(seq, tiles$start, tiles$stop)
  writeLines(as.vector(rbind(paste0(">", ids), subs)), file)
  invisible(ids)
}

#' Bitscore-weighted LCA per annotated element from stored HSPs
#' @noRd
chimera_gene_lca <- function(hits, feats, min_cov = 0.5) {
  out <- lapply(seq_len(nrow(feats)), function(i) {
    s <- feats$start[i]; e <- feats$stop[i]; len <- e - s + 1
    ov <- pmin(hits$qend, e) - pmax(hits$qstart, s) + 1
    h <- hits[!is.na(ov) & ov >= min_cov * len & !is.na(hits$lineage_taxids), , drop = FALSE]
    h <- h[order(-h$bitscore), , drop = FALSE]
    h <- h[!duplicated(h$accession), , drop = FALSE]
    base <- data.frame(element = feats$element[i], taxon = NA_character_, rank = NA_character_,
      taxid = NA_integer_, species = NA_character_, pident = NA_real_, n_hits = nrow(h),
      lineage_taxids = NA_character_, stringsAsFactors = FALSE)
    if (nrow(h) < 3) return(base)
    lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
    if (is.null(lca) || is.na(lca$taxid %||% NA)) return(base)
    base$taxon <- lca$taxon; base$rank <- lca$rank; base$taxid <- as.integer(lca$taxid)
    base$species <- h$species[1]; base$pident <- mean(h$pident, na.rm = TRUE)
    parts <- strsplit(h$lineage_taxids[1], ";", fixed = TRUE)[[1]]
    k <- match(lca$rank, RRNA_LINEAGE_RANKS)
    if (!is.na(k) && k < length(parts)) parts[(k + 1L):length(parts)] <- ""
    base$lineage_taxids <- paste(parts, collapse = ";")
    base
  })
  do.call(rbind, out)
}

#' @noRd
chimera_unavailable <- function(availability = "unavailable") {
  list(score = 0, junctions = integer(0), reason = NA_character_, availability = availability)
}

#' Channel 1: two resolved elements whose LCAs sit in different classes
#' @noRd
chimera_channel_genes <- function(gene_taxa, feats = NULL) {
  g <- gene_taxa[!is.na(gene_taxa$taxon), , drop = FALSE]
  if (nrow(g) < 2) return(chimera_unavailable("uninformative"))
  worst <- NULL
  for (i in seq_len(nrow(g) - 1)) for (j in (i + 1):nrow(g)) {
    d <- rank_distance(g$lineage_taxids[i], g$lineage_taxids[j])
    if (!is.na(d) && d >= 5 && (is.null(worst) || d > worst$d)) worst <- list(i = i, j = j, d = d)
  }
  if (is.null(worst)) return(list(score = 0, junctions = integer(0),
    reason = sprintf("elements agree: %s", paste(sprintf("%s %s", g$element, g$taxon), collapse = ", ")),
    availability = "available"))
  a <- g[worst$i, ]; b <- g[worst$j, ]
  junction <- NA_integer_
  if (!is.null(feats)) {
    fa <- feats[feats$element == a$element, ]; fb <- feats[feats$element == b$element, ]
    if (nrow(fa) && nrow(fb)) junction <- as.integer(round((min(fa$stop, fb$stop) + max(fa$start, fb$start)) / 2))
  }
  list(score = 0.6, junctions = junction[!is.na(junction)],
    reason = sprintf("%s resolves to %s (%s, %.1f%%) but %s to %s (%s, %.1f%%)",
      a$element, a$taxon, a$species, a$pident, b$element, b$taxon, b$species, b$pident),
    availability = "available")
}
