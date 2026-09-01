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

#' Best hit per tile from a parse_blast_hits_file() frame
#' @noRd
chimera_window_hits <- function(hits, tiles, sample_taxids = NA_character_,
                                taxdump_dir = NULL, db_source = NA_character_) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  best <- if (nrow(hits)) hits[order(hits$qseqid, -hits$bitscore), , drop = FALSE] else hits
  best <- best[!duplicated(best$qseqid), , drop = FALSE]
  i <- match(ids, best$qseqid)
  lin <- tryCatch(tk_ranked_lineage(best$staxid, taxdump_dir), error = function(e) data.frame())
  li <- if (nrow(lin)) match(best$staxid[i], lin$taxid) else rep(NA_integer_, length(i))
  lt <- if (nrow(lin)) lin$lineage_taxids[li] else rep(NA_character_, length(i))
  rd <- if (!is.na(sample_taxids) && nzchar(sample_taxids)) {
    vapply(lt, function(x) if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, sample_taxids)),
      integer(1), USE.NAMES = FALSE)
  } else rep(NA_integer_, length(i))
  data.frame(start = tiles$start, stop = tiles$stop,
    accession = best$accession[i], species = extract_species(best$title[i]),
    staxid = best$staxid[i], lineage = if (nrow(lin)) lin$lineage[li] else NA_character_,
    lineage_taxids = lt, rank_distance = rd,
    hit_class = classify_hit_title(best$title[i]),
    pident = best$pident[i], bitscore = best$bitscore[i], db_source = db_source,
    stringsAsFactors = FALSE)
}

#' Name at a rank from a ';'-joined 7-field lineage
#' @noRd
chimera_rank_name <- function(lineage, rank = "class") {
  k <- match(rank, RRNA_LINEAGE_RANKS)
  vapply(strsplit(lineage %||% NA_character_, ";", fixed = TRUE), function(p) {
    if (length(p) >= k && nzchar(p[k])) p[k] else NA_character_
  }, character(1))
}

#' Collapse windows into runs sharing a name at the deciding rank
#' @noRd
chimera_segments <- function(windows, deciding = "class") {
  if (!nrow(windows)) return(data.frame())
  key <- chimera_rank_name(windows$lineage, deciding)
  key[is.na(key)] <- ifelse(is.na(windows$pident[is.na(key)]), "nohit",
    ifelse(windows$pident[is.na(key)] >= CHIMERA_PID_DROP, "high_id", "low_id"))
  run <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
  do.call(rbind, lapply(split(seq_along(run), run), function(ix) {
    w <- windows[ix, , drop = FALSE]
    top <- ix[which.max(w$bitscore %||% rep(0, nrow(w)))]
    data.frame(start = w$start[1], stop = w$stop[nrow(w)],
      taxon = key[ix[1]], rank = deciding,
      species = windows$species[top] %||% NA_character_,
      pident = mean(w$pident, na.rm = TRUE), n_windows = nrow(w),
      db_source = w$db_source[1], lineage_taxids = windows$lineage_taxids[top] %||% NA_character_,
      stringsAsFactors = FALSE)
  }))
}

#' Channel 2: adjacent segments of comparable length in different classes
#' @noRd
chimera_channel_windows <- function(segments, identity_only = FALSE) {
  if (!nrow(segments)) return(chimera_unavailable())
  s <- segments[segments$n_windows >= 2, , drop = FALSE]
  if (nrow(s) < 2) return(list(score = 0, junctions = integer(0),
    reason = sprintf("one lineage along the scaffold: %s", segments$taxon[1]), availability = "available"))
  best <- NULL
  for (i in seq_len(nrow(s) - 1)) {
    a <- s[i, ]; b <- s[i + 1, ]
    if (identical(a$taxon, b$taxon)) next
    if (a$taxon %in% c("nohit") || b$taxon %in% c("nohit")) next
    la <- a$stop - a$start + 1; lb <- b$stop - b$start + 1
    if (min(la, lb) < 0.25 * max(la, lb)) next
    if (!identity_only && rank_distance(a$lineage_taxids, b$lineage_taxids) < 5) next
    best <- list(a = a, b = b, junction = (a$stop + b$start) %/% 2L)
  }
  if (is.null(best)) return(list(score = 0, junctions = integer(0),
    reason = "no class-level switch between comparable segments", availability = "available"))
  desc <- function(x) sprintf("%d-%d %s (%s, %.1f%%)", x$start, x$stop, x$taxon, x$species %||% "", x$pident)
  list(score = if (identity_only) 0 else 0.5, junctions = best$junction,
    reason = paste(desc(best$a), desc(best$b), sep = "; "),
    availability = if (identity_only) "identity_only" else "available")
}
