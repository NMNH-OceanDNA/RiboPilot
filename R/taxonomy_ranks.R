# Rank-aware taxonomy primitives for the contamination screen.
#
# Every taxonomic comparison in RiboPilot used to run over a "; "-joined lineage
# string whose empty rank slots had already been collapsed away, so a taxon
# missing (say) an order had its family sitting in the order slot and positional
# comparison silently under-counted shared ranks. That hits hardest on the
# under-described taxa RiboPilot exists to serve, i.e. exactly backwards.
#
# Here a lineage is always SEVEN rank-aligned fields, carried as both names and
# taxids (`taxonkit reformat -F -t`). Missing ranks are filled with an
# "unclassified <parent>" placeholder name and an EMPTY taxid, so "this rank is
# not real" is unambiguous. Distances are computed on the taxids, which are
# immune to synonyms, homonyms and defline formatting; names are for display.
#
# The core is pure and DB-free (mirroring curate_rrna_core); the taxonkit calls
# live in the thin wrappers at the bottom.

#' The seven ranks carried by a RiboPilot lineage, root -> tip.
#' @noRd
RRNA_LINEAGE_RANKS <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

#' taxonkit reformat format string producing [RRNA_LINEAGE_RANKS].
#'
#' `{K}` is kingdom (Metazoa / Fungi / Viridiplantae) on a current NCBI taxdump,
#' where the old superkingdom rank has been renamed to domain (`{d}`).
#' @noRd
RRNA_RANK_FORMAT <- "{K};{p};{c};{o};{f};{g};{s}"

#' Rank-distance label for a distance returned by [rank_distance()].
#' Index i holds the label for distance i - 1.
#' @noRd
RRNA_RANK_DISTANCE_LABELS <- c(RRNA_LINEAGE_RANKS[7:1], "none")

#' Names that carry no usable taxonomic signal
#'
#' `uncultured`/`environmental` records dominate core_nt 18S and resolve to a
#' kingdom-only lineage, which every rank rule would otherwise read as maximum
#' taxonomic distance. `unclassified` is what `taxonkit reformat -F` writes into
#' a rank it had to fill.
#' @param x Character vector of taxon names.
#' @return logical vector.
#' @noRd
is_placeholder_taxon <- function(x) {
  x <- ifelse(is.na(x), "", trimws(as.character(x)))
  !nzchar(x) | grepl("^(unclassified|uncultured|unidentified|environmental|unassigned|metagenome)\\b",
    x, ignore.case = TRUE) | grepl("\\b(environmental samples|metagenome)\\b", x, ignore.case = TRUE)
}

#' Split a 7-field rank-aligned taxid string into an integer vector
#' @param x A ";"-joined taxid string from `taxonkit reformat -t`.
#' @return integer vector of length 7; NA where the rank is absent or filled.
#' @noRd
parse_rank_taxids <- function(x) {
  out <- rep(NA_integer_, length(RRNA_LINEAGE_RANKS))
  if (length(x) == 0 || is.na(x[1]) || !nzchar(trimws(x[1]))) {
    return(stats::setNames(out, RRNA_LINEAGE_RANKS))
  }
  p <- trimws(strsplit(as.character(x[1]), ";", fixed = TRUE)[[1]])
  n <- min(length(p), length(out))
  if (n > 0) out[seq_len(n)] <- suppressWarnings(as.integer(p[seq_len(n)]))
  stats::setNames(out, RRNA_LINEAGE_RANKS)
}

#' Split a 7-field rank-aligned lineage into a named character vector
#' @param x A ";"-joined lineage string.
#' @return character vector of length 7; NA at placeholder or absent ranks.
#' @noRd
parse_rank_names <- function(x) {
  out <- rep(NA_character_, length(RRNA_LINEAGE_RANKS))
  if (length(x) == 0 || is.na(x[1]) || !nzchar(trimws(x[1]))) {
    return(stats::setNames(out, RRNA_LINEAGE_RANKS))
  }
  p <- trimws(strsplit(as.character(x[1]), ";", fixed = TRUE)[[1]])
  n <- min(length(p), length(out))
  if (n > 0) out[seq_len(n)] <- p[seq_len(n)]
  out[is_placeholder_taxon(out)] <- NA_character_
  stats::setNames(out, RRNA_LINEAGE_RANKS)
}

#' Deepest rank position at which two rank-aligned taxid vectors agree
#'
#' A real lineage that matches at species also matches at every higher rank, so
#' the deepest agreeing position is taken directly rather than as an unbroken
#' prefix. That is what lets a taxon with an unassigned intermediate rank still
#' match at the ranks it does have.
#' @return integer position in 1..7, or 0 when nothing is shared.
#' @noRd
deepest_shared_rank <- function(a, b) {
  n <- length(RRNA_LINEAGE_RANKS)
  hit <- which(!is.na(a[seq_len(n)]) & !is.na(b[seq_len(n)]) & a[seq_len(n)] == b[seq_len(n)])
  if (length(hit) == 0) 0L else max(hit)
}

#' Taxonomic distance between a hit and the sample, in rank steps
#'
#' 0 = same species, 1 = same genus, 2 = family, 3 = order, 4 = class,
#' 5 = phylum, 6 = kingdom, 7 = no shared taxon at any rank.
#'
#' Distance is deliberately NOT a measure of sequence divergence. A large
#' distance on its own is the normal state for an under-represented taxon and
#' must never be read as contamination; it only becomes evidence when paired
#' with an implausibly high percent identity (see the screen envelope).
#'
#' @param hit_taxids,sample_taxids Rank-aligned taxid vectors from
#'   [parse_rank_taxids()], or ";"-joined taxid strings.
#' @return integer in 0..7.
#' @export
rank_distance <- function(hit_taxids, sample_taxids) {
  a <- if (is.character(hit_taxids)) parse_rank_taxids(hit_taxids) else hit_taxids
  b <- if (is.character(sample_taxids)) parse_rank_taxids(sample_taxids) else sample_taxids
  pos <- deepest_shared_rank(a, b)
  if (pos == 0L) length(RRNA_LINEAGE_RANKS) else length(RRNA_LINEAGE_RANKS) - pos
}

#' Human-readable label for a [rank_distance()] value ("genus", "class", "none")
#' @param d Integer distance in 0..7.
#' @return character.
#' @export
rank_distance_label <- function(d) {
  d <- suppressWarnings(as.integer(d))
  ifelse(is.na(d) | d < 0 | d > 7, NA_character_, RRNA_RANK_DISTANCE_LABELS[pmin(pmax(d, 0), 7) + 1L])
}

#' Is a hit at least as close as `min_rank`?
#'
#' @param d Distance from [rank_distance()].
#' @param min_rank One of [RRNA_LINEAGE_RANKS]; "kingdom" is the most permissive.
#' @return logical.
#' @export
rank_at_least <- function(d, min_rank = "phylum") {
  i <- match(min_rank, RRNA_LINEAGE_RANKS)
  if (is.na(i)) stop("min_rank must be one of: ", paste(RRNA_LINEAGE_RANKS, collapse = ", "), call. = FALSE)
  !is.na(d) & d <= (length(RRNA_LINEAGE_RANKS) - i)
}

#' Deepest taxon name shared between a hit and the sample
#'
#' The thing a biologist actually wants to read: "shares: Notropis (genus)" vs
#' "shares: Metazoa (kingdom)" vs "shares: nothing".
#' @param hit_names Rank-aligned names from [parse_rank_names()].
#' @param hit_taxids,sample_taxids Rank-aligned taxid vectors.
#' @return list(taxon, rank, distance).
#' @export
shared_taxon <- function(hit_names, hit_taxids, sample_taxids) {
  a <- if (is.character(hit_taxids)) parse_rank_taxids(hit_taxids) else hit_taxids
  b <- if (is.character(sample_taxids)) parse_rank_taxids(sample_taxids) else sample_taxids
  nm <- if (is.character(hit_names) && length(hit_names) == 1) parse_rank_names(hit_names) else hit_names
  pos <- deepest_shared_rank(a, b)
  if (pos == 0L) {
    return(list(taxon = NA_character_, rank = "none", distance = length(RRNA_LINEAGE_RANKS)))
  }
  list(taxon = unname(nm[pos]), rank = RRNA_LINEAGE_RANKS[pos],
    distance = length(RRNA_LINEAGE_RANKS) - pos)
}

#' Bitscore-weighted lowest common ancestor over the top BLAST hits
#'
#' Defence against a single mislabelled GenBank accession, which a one-hit
#' screen cannot survive. Considers only hits within `frac` of the top bitscore
#' (always at least `min_hits` when that many are available) and returns the
#' deepest rank at which they all agree.
#'
#' @param taxid_paths Character vector of ";"-joined rank-aligned taxid strings,
#'   one per hit, or a list of integer vectors.
#' @param bitscores Numeric vector of the same length.
#' @param names_paths Optional matching lineage-name strings, for the returned label.
#' @param frac Bitscore window below the top hit (default 0.10).
#' @param min_hits Always include at least this many top hits when available.
#' @return list(taxid, taxon, rank, n_hits).
#' @export
bitscore_weighted_lca <- function(taxid_paths, bitscores, names_paths = NULL,
                                  frac = 0.10, min_hits = 3L) {
  keep <- !is.na(bitscores)
  if (!any(keep)) {
    return(list(taxid = NA_integer_, taxon = NA_character_, rank = "none", n_hits = 0L))
  }
  ord <- order(bitscores, decreasing = TRUE, na.last = NA)
  top <- bitscores[ord[1]]
  sel <- ord[bitscores[ord] >= top * (1 - frac)]
  if (length(sel) < min_hits) sel <- utils::head(ord, min_hits)

  mat <- vapply(sel, function(i) {
    x <- taxid_paths[[i]]
    if (is.character(x)) parse_rank_taxids(x) else x
  }, integer(length(RRNA_LINEAGE_RANKS)))
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = length(RRNA_LINEAGE_RANKS))

  agree <- apply(mat, 1, function(r) {
    r <- r[!is.na(r)]
    length(r) > 0 && length(unique(r)) == 1
  })
  # a rank only counts when every selected hit resolved it
  resolved <- apply(mat, 1, function(r) all(!is.na(r)))
  pos <- which(agree & resolved)
  if (length(pos) == 0) {
    return(list(taxid = NA_integer_, taxon = NA_character_, rank = "none", n_hits = length(sel)))
  }
  p <- max(pos)
  nm <- NA_character_
  if (!is.null(names_paths)) {
    nmv <- names_paths[[sel[1]]]
    nm <- unname((if (is.character(nmv) && length(nmv) == 1) parse_rank_names(nmv) else nmv)[p])
  }
  list(taxid = unname(mat[p, 1]), taxon = nm, rank = RRNA_LINEAGE_RANKS[p], n_hits = length(sel))
}

# ---- taxonkit-backed resolution ---------------------------------------------

#' Per-process taxdump extraction cache
#'
#' `write_blast_screen()` used to call [extract_taxdump()] per unit, so a
#' 30-unit sample serially unpacked ~500 MB thirty times. Extractions are cached
#' by source directory for the life of the R process and removed at exit.
#' @noRd
.taxdump_cache <- new.env(parent = emptyenv())

#' Extracted taxdump directory for `dir`, unpacking at most once per process
#' @param dir Directory holding a pinned taxdump (tarball or extracted dmps).
#' @return Path to a directory containing names.dmp / nodes.dmp, or NULL.
#' @noRd
taxdump_session <- function(dir) {
  if (!have_taxdump(dir)) {
    return(NULL)
  }
  if (file.exists(file.path(dir, "names.dmp"))) {
    return(dir)
  }
  key <- normalizePath(dir, mustWork = FALSE)
  hit <- .taxdump_cache[[key]]
  if (!is.null(hit) && dir.exists(hit)) {
    return(hit)
  }
  ex <- extract_taxdump(dir)
  assign(key, ex, envir = .taxdump_cache)
  reg.finalizer(.taxdump_cache, function(e) {
    for (p in ls(e)) unlink(get(p, envir = e), recursive = TRUE)
  }, onexit = TRUE)
  ex
}

#' Back-fill `samples.lineage_taxids` from `samples.taxid` where missing
#'
#' The rank-aligned taxid path is what [rank_distance()] compares, but taxonomy
#' is frozen from three places (the review CSV, the Shiny taxonomy tab, and
#' `resolve_project_taxonomy`). Rather than patch each writer, the path is
#' derived on demand and cached back, which also upgrades projects created
#' before schema v9.
#'
#' @param con Open project DB connection.
#' @param ids Optional character vector of sample IDs (default: all).
#' @param taxdump_dir Directory holding a pinned taxdump.
#' @return (invisibly) the number of samples filled.
#' @noRd
ensure_sample_lineage_taxids <- function(con, ids = NULL, taxdump_dir = NULL) {
  if (!("lineage_taxids" %in% DBI::dbListFields(con, "samples"))) {
    return(invisible(0L))
  }
  q <- "SELECT ID, taxid FROM samples WHERE taxid IS NOT NULL
        AND (lineage_taxids IS NULL OR lineage_taxids = '')"
  need <- tryCatch(DBI::dbGetQuery(con, q), error = function(e) data.frame())
  if (!is.null(ids)) need <- need[need$ID %in% ids, , drop = FALSE]
  if (nrow(need) == 0) {
    return(invisible(0L))
  }
  lin <- tk_ranked_lineage(need$taxid, taxdump_dir)
  if (nrow(lin) == 0) {
    return(invisible(0L))
  }
  i <- match(need$taxid, lin$taxid)
  ok <- !is.na(i)
  if (!any(ok)) {
    return(invisible(0L))
  }
  upd <- DBI::dbSendStatement(con,
    "UPDATE samples SET lineage_taxids = ? WHERE ID = ?")
  on.exit(DBI::dbClearResult(upd), add = TRUE)
  for (k in which(ok)) {
    DBI::dbBind(upd, list(lin$lineage_taxids[i[k]], need$ID[k]))
  }
  invisible(sum(ok))
}

#' Rank-aligned lineages (names + taxids) for a set of taxids
#'
#' One batched taxonkit call. Missing ranks are filled (`-F`) so the seven rank
#' slots always line up positionally, and the parallel taxid string (`-t`) marks
#' a filled rank with an empty taxid.
#'
#' @param taxids Integer vector of NCBI taxids.
#' @param dir Directory holding a pinned taxdump.
#' @return data.frame(taxid, lineage, lineage_taxids).
#' @export
tk_ranked_lineage <- function(taxids, dir) {
  taxids <- unique(stats::na.omit(suppressWarnings(as.integer(taxids))))
  empty <- data.frame(taxid = integer(0), lineage = character(0),
    lineage_taxids = character(0), stringsAsFactors = FALSE)
  if (length(taxids) == 0) {
    return(empty)
  }
  dd <- taxdump_session(dir)
  if (is.null(dd)) {
    return(empty)
  }
  ids <- paste(taxids, collapse = "\n")
  lin <- tk_run("lineage", ids, dd)
  if (length(lin) == 0) {
    return(empty)
  }
  out <- tk_run(c("reformat", "-F", "-t", "-f", RRNA_RANK_FORMAT), paste(lin, collapse = "\n"), dd)
  rows <- lapply(strsplit(out, "\t"), function(p) {
    if (length(p) < 4) {
      return(NULL)
    }
    n <- length(p)
    data.frame(taxid = suppressWarnings(as.integer(p[1])),
      lineage = trimws(p[n - 1]), lineage_taxids = trimws(p[n]),
      stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(empty)
  }
  d <- do.call(rbind, rows)
  d[!is.na(d$taxid), , drop = FALSE]
}
