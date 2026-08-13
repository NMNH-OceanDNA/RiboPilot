# Scaffold / path scoring and taxonomy-aware auto-select.
#
# GetOrganelle nr emits several candidate paths (each several scaffolds). This
# module scores candidates and auto-selects the one that best represents the
# sample's nuclear rDNA, using: depth (from the assembly graph), BLAST identity,
# anchor completeness (18S/5.8S/28S present), length plausibility, and how
# closely the candidate's top-BLAST-hit lineage matches the sample's resolved
# lineage. See dev/design/07_SCAFFOLD_SELECTION.md.

#' Default composite-score weights (sum need not be 1; relative magnitude matters)
#' @noRd
DEFAULT_SCORE_WEIGHTS <- c(tax = 0.30, identity = 0.20, anchors = 0.20, depth = 0.20, length = 0.10)

# ---- Depth from the assembly graph ------------------------------------------

#' Parse per-segment depth from a GetOrganelle/SPAdes GFA
#'
#' Depth = RC/LN per segment (read count over length). GetOrganelle path fasta
#' headers list the segment ids a scaffold traverses, so a scaffold's depth is
#' the length-weighted mean of its segments' depths.
#'
#' @param gfa_path Path to a `.gfa` file.
#' @return data.frame(segment, length, rc, depth)
#' @noRd
parse_gfa_depth <- function(gfa_path) {
  lines <- readLines(gfa_path, warn = FALSE)
  s <- lines[startsWith(lines, "S\t") | startsWith(lines, "S ")]
  if (length(s) == 0) {
    return(data.frame(segment = character(0), length = integer(0), rc = numeric(0), depth = numeric(0)))
  }
  f <- strsplit(gsub(" ", "\t", s), "\t")
  seg <- vapply(f, `[`, character(1), 2)
  tags <- function(p, key) {
    hit <- grep(paste0("^", key, ":i:"), p, value = TRUE)
    if (length(hit) == 0) NA_real_ else as.numeric(sub(paste0(key, ":i:"), "", hit[1]))
  }
  ln <- vapply(f, tags, numeric(1), key = "LN")
  rc <- vapply(f, tags, numeric(1), key = "RC")
  # fall back to sequence length if LN tag missing
  ln[is.na(ln)] <- vapply(f[is.na(ln)], function(p) nchar(p[3] %||% ""), numeric(1))
  data.frame(segment = seg, length = ln, rc = rc, depth = ifelse(ln > 0, rc / ln, NA_real_), stringsAsFactors = FALSE)
}

#' Length-weighted mean depth of the segments a scaffold traverses
#'
#' @param header A GetOrganelle scaffold header, e.g.
#'   `scaffold_1--24280-,32579-,31679-` (segment ids after `--`, orientation sign
#'   stripped). Merged segments carry underscores (`62361_60645`), so only the
#'   trailing `+`/`-` may be removed; the rest of the id is literal.
#' @param gfa_depth data.frame from [parse_gfa_depth()].
#' @return numeric mean depth (NA if no segments matched)
#' @noRd
scaffold_depth <- function(header, gfa_depth) {
  ids <- sub(".*--", "", header)
  ids <- trimws(unlist(strsplit(ids, ",")))
  ids <- sub("[+-]$", "", ids)
  ids <- ids[nzchar(ids)]
  i <- match(ids, gfa_depth$segment)
  if (anyNA(i)) {
    warning("scaffold header names segments missing from the GFA: ",
      paste(unique(ids[is.na(i)]), collapse = ", "), call. = FALSE)
  }
  d <- gfa_depth[i[!is.na(i)], , drop = FALSE]
  d <- d[!is.na(d$depth), , drop = FALSE]
  if (nrow(d) == 0) {
    return(NA_real_)
  }
  stats::weighted.mean(d$depth, d$length)
}

# ---- Taxonomy match ---------------------------------------------------------

#' Taxonomic match between a candidate lineage and the sample lineage
#'
#' Both are rank-aligned lineages (see [parse_rank_names()]). Scored on the
#' DEEPEST rank the two agree at, not on an unbroken prefix: a real lineage that
#' agrees at genus agrees at every rank above it, so a taxon with an unassigned
#' intermediate rank should still score as the close match it is. The prefix walk
#' this replaces broke at the first gap and scored a congeneric hit near zero,
#' which hit the under-described taxa RiboPilot exists to serve hardest.
#'
#' Taxids are preferred when available: they are immune to the homonyms and
#' spelling variants that names carry.
#'
#' @param candidate_lineage,sample_lineage Rank-aligned ";"-joined lineages.
#' @param use_taxids Treat both as taxid strings rather than names.
#' @return numeric in [0, 1]; 1 = agrees at species, 0 = agrees nowhere.
#' @noRd
taxonomic_match <- function(candidate_lineage, sample_lineage, use_taxids = FALSE) {
  parse <- if (use_taxids) parse_rank_taxids else parse_rank_names
  a <- parse(candidate_lineage)
  b <- parse(sample_lineage)
  if (!use_taxids) {
    a <- tolower(a)
    b <- tolower(b)
  }
  deepest_shared_rank(a, b) / length(RRNA_LINEAGE_RANKS)
}

# ---- Composite scoring ------------------------------------------------------

#' Score and rank candidate paths, auto-selecting the best
#'
#' @param candidates data.frame, one row per candidate path, with columns:
#'   `path`, `depth` (mean graph depth), `blast_pident` (0-100), `blast_lineage`
#'   (rank-aligned lineage of the top hit), `anchors` (comma list of distinct
#'   anchors found among 18S/5.8S/28S), `length` (total bp), and optionally
#'   `blast_lineage_taxids`.
#' @param sample_lineage The sample's resolved rank-aligned lineage.
#' @param sample_taxids The sample's rank-aligned lineage taxids. When both this
#'   and `candidates$blast_lineage_taxids` are present the taxonomy term is
#'   scored on taxids instead of names.
#' @param expected_length Numeric length-2 c(min, max) plausible rDNA repeat unit length
#'   for the clade (default c(4000, 12000)).
#' @param weights Named weights (see [DEFAULT_SCORE_WEIGHTS]).
#' @return `candidates` with added sub-score columns (`s_tax`, `s_identity`,
#'   `s_anchors`, `s_depth`, `s_length`), a `score`, a `rank` (1 = best), and an
#'   `auto_select` logical flagging the top candidate.
#' @export
score_candidates <- function(candidates,
                             sample_lineage = NA_character_,
                             sample_taxids = NA_character_,
                             expected_length = c(4000, 12000),
                             weights = DEFAULT_SCORE_WEIGHTS) {
  n <- nrow(candidates)
  if (n == 0) {
    return(candidates)
  }
  anchor_count <- function(a) {
    if (is.na(a) || !nzchar(a)) {
      return(0L)
    }
    length(intersect(RRNA_GENES, trimws(strsplit(a, ",")[[1]])))
  }

  # taxids where both sides have them, names otherwise
  by_taxid <- !is.na(sample_taxids[1]) && !is.null(candidates$blast_lineage_taxids)
  cand_lin <- if (by_taxid) candidates$blast_lineage_taxids else candidates$blast_lineage %||% rep(NA, n)
  s_tax <- vapply(cand_lin,
    function(l) taxonomic_match(l, if (by_taxid) sample_taxids else sample_lineage,
      use_taxids = by_taxid), numeric(1))
  s_identity <- pmin(pmax((candidates$blast_pident %||% rep(NA, n)) / 100, 0), 1)
  s_identity[is.na(s_identity)] <- 0
  s_anchors <- vapply(candidates$anchors %||% rep(NA, n), anchor_count, integer(1)) / length(RRNA_GENES)
  # depth: normalize to the deepest candidate (relative coverage rank)
  dep <- candidates$depth %||% rep(NA, n)
  s_depth <- if (all(is.na(dep))) rep(0, n) else dep / max(dep, na.rm = TRUE)
  s_depth[is.na(s_depth)] <- 0
  # length plausibility: 1 inside [min,max], linear decay outside over 50% margin
  len <- candidates$length %||% rep(NA, n)
  s_length <- vapply(len, function(L) {
    if (is.na(L)) {
      return(0)
    }
    lo <- expected_length[1]
    hi <- expected_length[2]
    if (L >= lo && L <= hi) {
      return(1)
    }
    if (L < lo) max(0, 1 - (lo - L) / (0.5 * lo)) else max(0, 1 - (L - hi) / (0.5 * hi))
  }, numeric(1))

  w <- weights
  score <- w["tax"] * s_tax + w["identity"] * s_identity + w["anchors"] * s_anchors +
    w["depth"] * s_depth + w["length"] * s_length

  out <- candidates
  out$s_tax <- round(s_tax, 3)
  out$s_identity <- round(s_identity, 3)
  out$s_anchors <- round(s_anchors, 3)
  out$s_depth <- round(s_depth, 3)
  out$s_length <- round(s_length, 3)
  out$score <- round(as.numeric(score), 4)
  out$rank <- rank(-out$score, ties.method = "first")
  out$auto_select <- out$rank == 1
  out[order(out$rank), ]
}
