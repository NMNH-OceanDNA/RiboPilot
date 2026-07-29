# The screen's taxonomy channel: does this unit's hit distribution agree on a
# taxon, and is that taxon the host's.
#
# It uses NO identity threshold. The identity envelope was deferred after
# measurement against the reference project (doc 08 section 3.2.1): that corpus
# has 2 hits at phylum distance and none beyond, so it cannot calibrate a bound
# acting there, and genuine host 18S reaches 99.6 percent at class distance, so
# no ceiling can sit below that without flagging real host data.
#
# What this cannot do: separate two metazoans. A copepod and a fish both sit
# under Metazoa, so concordance against the sample's path says nothing. Doc 08
# section 3.5 puts those cases in the depth and cross-sample channels.

#' Scoring version. Bump to force a recompute of stored rows.
#' @noRd
SCREEN_VERSION <- 1L

#' @noRd
SCREEN_VERDICTS <- c("on_target", "off_target", "no_reference", "uncertain")

#' @noRd
SCREEN_CLASSES <- c("bacterial", "archaeal", "fungal", "algal_plant",
  "other_euk", "metazoan_other", "human")

#' @noRd
SCREEN_NO_REF_CODES <- c("no_hits", "no_taxdump", "unresolved_sample_lineage",
  "unresolved_hit_lineages")

#' Confidence weights.
#'
#' Priors, not measurements, which is why the screen rule ships as Flag rather
#' than Drop. Any change here must bump [SCREEN_VERSION].
#' @noRd
SCREEN_WEIGHTS <- list(support = 0.30, resolution = 0.30, separation = 0.40)

#' Confidence below which a verdict is recorded as uncertain.
#' @noRd
SCREEN_MIN_CONFIDENCE <- 0.6

#' Depth of a resolved rank, 0 (nothing) to 1 (species)
#' @noRd
screen_rank_depth <- function(rank) {
  if (length(rank) == 0 || is.na(rank[1])) return(0)
  i <- match(rank[1], RRNA_LINEAGE_RANKS)
  if (is.na(i)) 0 else i / length(RRNA_LINEAGE_RANKS)
}

#' The three confidence terms and their weighted mean
#'
#' `unanimous` means every surviving hit fell on one side of the concordance
#' test, so there is no opposing hit to measure a margin against. That is
#' MAXIMAL separation, not absent separation: a unit whose hits all agree with
#' the host is the cleanest case there is, and scoring it 0 would put the
#' cleanest possible scaffold below the uncertainty gate.
#' @noRd
screen_confidence_terms <- function(lca_n_hits, lca_rank, margin, top_bitscore,
                                    unanimous = FALSE) {
  support <- min(1, (lca_n_hits %||% 0) / 10)
  resolution <- screen_rank_depth(lca_rank)
  separation <- if (isTRUE(unanimous)) {
    1
  } else if (is.na(margin %|NA|% NA) || is.na(top_bitscore %|NA|% NA) ||
             !isTRUE(top_bitscore > 0)) {
    0
  } else {
    min(1, abs(margin) / top_bitscore)
  }
  w <- SCREEN_WEIGHTS
  list(support = support, resolution = resolution, separation = separation,
    confidence = unname(support * w$support + resolution * w$resolution +
      separation * w$separation))
}

#' Well-known kingdom/domain-rank taxids, stable across taxdump versions.
#'
#' Classification keys off the taxid rather than the kingdom-rank NAME:
#' taxids are immune to synonyms, homonyms and defline formatting, exactly the
#' reason [rank_distance()] compares taxids and treats names as display-only.
#' @noRd
SCREEN_KINGDOM_TAXIDS <- c("2" = "bacterial", "2157" = "archaeal",
  "4751" = "fungal", "33090" = "algal_plant", "33208" = "metazoan_other")

#' Off-target subclass from the LCA's rank-aligned taxid path and its taxon label
#'
#' Matches the WHOLE path rather than a single rank slot: whether taxonkit's
#' `{K}` emits superkingdom or kingdom decides whether Bacteria/Archaea or
#' Fungi/Viridiplantae/Metazoa lands in slot 1, and this codebase disagrees
#' with itself about which (doc 08 line 510 vs taxonomy_ranks.R:24). Matching
#' the whole path is correct under either reading.
#' @param lca_taxids Rank-aligned taxid path of the LCA (";"-joined string or
#'   integer vector from [parse_rank_taxids()]).
#' @param lca_taxon The LCA's own taxon label (used only for the human check).
#' @noRd
screen_class_for <- function(lca_taxids, lca_taxon) {
  if (isTRUE(tolower(lca_taxon %|NA|% "") == "homo sapiens")) return("human")
  p <- if (is.character(lca_taxids)) parse_rank_taxids(lca_taxids) else lca_taxids
  p <- p[!is.na(p)]
  if (!length(p)) return(NA_character_)
  # most specific first: a metazoan path also contains Eukaryota, and slot 1
  # may be superkingdom or kingdom depending on the taxdump
  for (k in c("4751", "33090", "33208", "2", "2157")) {
    if (as.integer(k) %in% p) return(unname(SCREEN_KINGDOM_TAXIDS[k]))
  }
  "other_euk"
}

#' Score one unit's hit distribution against its sample's lineage
#'
#' @param hits Frame with `hit_class`, `bitscore`, `lineage_taxids`, `lineage`.
#' @param sample_taxids The sample's rank-aligned taxid path.
#' @return list(verdict, class, confidence, reason, reason_code, terms,
#'   lca_taxid, lca_taxon, lca_rank, lca_n_hits, host_concordant, margin).
#' @noRd
screen_channel_taxonomy <- function(hits, sample_taxids) {
  out <- list(verdict = "no_reference", class = NA_character_, confidence = 0,
    reason = NA_character_, reason_code = "no_hits",
    terms = list(support = 0, resolution = 0, separation = 0),
    lca_taxid = NA_integer_, lca_taxon = NA_character_, lca_rank = NA_character_,
    lca_n_hits = 0L, host_concordant = NA_integer_, margin = NA_real_)

  samp <- parse_rank_taxids(sample_taxids %||% NA_character_)
  if (all(is.na(samp))) {
    out$reason_code <- "unresolved_sample_lineage"
    out$reason <- "the sample's own lineage could not be resolved"
    return(out)
  }
  if (is.null(hits) || nrow(hits) == 0) {
    out$reason <- "no GenBank hits are stored for this scaffold"
    return(out)
  }

  # uncultured/environmental records carry their junk label at the deepest
  # FILLED rank (taxonkit still writes it into the species slot), which every
  # rank rule would otherwise read as a real species; a merely blank rank
  # (nothing filled in yet) is normal and must not be excluded on that basis
  # alone, so the check looks at the deepest non-empty raw field, not at
  # whichever rank position happens to be blank.
  deepest_raw <- vapply(hits$lineage, function(x) {
    f <- trimws(strsplit(x %|NA|% "", ";", fixed = TRUE)[[1]])
    f <- f[nzchar(f)]
    if (length(f) == 0) NA_character_ else f[length(f)]
  }, character(1), USE.NAMES = FALSE)
  keep <- hits$hit_class %in% "rdna" & !is.na(deepest_raw) & !is_placeholder_taxon(deepest_raw)
  h <- hits[keep, , drop = FALSE]
  if (nrow(h) == 0) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "no hit resolved to a usable lineage"
    return(out)
  }

  lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
  if (is.na(lca$taxid)) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "the hits agree on no shared taxon"
    out$lca_n_hits <- lca$n_hits
    return(out)
  }
  out$lca_taxid <- lca$taxid; out$lca_taxon <- lca$taxon
  out$lca_rank <- lca$rank;  out$lca_n_hits <- lca$n_hits
  out$host_concordant <- as.integer(lca$taxid %in% samp[!is.na(samp)])

  # the margin uses a SECOND, per-hit concordance, distinct from the unit-level
  # one above: it compares the two best-supported opposing hits.
  # a hit is on the sample's lineage when their paths share a resolved taxid at
  # PHYLUM or deeper. Sharing only the shallowest slot means "both are animals",
  # which every metazoan contaminant satisfies: three host hits plus two from a
  # different phylum would otherwise read as unanimous and score on_target.
  samp_ok <- samp[-1]
  samp_ok <- samp_ok[!is.na(samp_ok)]
  hit_conc <- vapply(h$lineage_taxids, function(x) {
    p <- parse_rank_taxids(x)[-1]
    any(p[!is.na(p)] %in% samp_ok)
  }, logical(1), USE.NAMES = FALSE)
  best <- function(v) if (any(v)) max(h$bitscore[v], na.rm = TRUE) else NA_real_
  out$margin <- best(hit_conc) - best(!hit_conc)

  # every surviving hit on one side means nothing opposes the winner, which is
  # maximal separation rather than missing separation. A lone hit is unopposed
  # by definition; without the count floor one accession reaches confidence
  # 0.73 and decides a scaffold on its own, which bitscore_weighted_lca()'s own
  # minimum-of-three exists to prevent.
  unanimous <- (all(hit_conc) || !any(hit_conc)) && lca$n_hits >= 2
  tm <- screen_confidence_terms(lca$n_hits, lca$rank, out$margin,
    suppressWarnings(max(h$bitscore, na.rm = TRUE)), unanimous = unanimous)
  out$terms <- tm[c("support", "resolution", "separation")]
  out$confidence <- tm$confidence
  out$reason_code <- NA_character_

  out$verdict <- if (tm$confidence < SCREEN_MIN_CONFIDENCE) {
    "uncertain"
  } else if (identical(out$host_concordant, 1L)) {
    "on_target"
  } else {
    "off_target"
  }
  if (identical(out$verdict, "off_target")) {
    # the whole path, not one rank slot: see screen_class_for()
    out$class <- screen_class_for(h$lineage_taxids[1], lca$taxon)
  }
  out$reason <- sprintf(
    "hits agree on %s (%s, %d voting); %s the sample's lineage. support %.2f w %.2f, resolution %.2f w %.2f, separation %.2f w %.2f -> confidence %.2f",
    lca$taxon %|NA|% "an unnamed taxon", lca$rank, lca$n_hits,
    if (identical(out$host_concordant, 1L)) "on" else "NOT on",
    tm$support, SCREEN_WEIGHTS$support, tm$resolution, SCREEN_WEIGHTS$resolution,
    tm$separation, SCREEN_WEIGHTS$separation, tm$confidence)
  out
}

#' Score and store one unit's screen verdict
#'
#' Skips a row already stored at the current [SCREEN_VERSION] unless `force`, so
#' re-running the pipeline is cheap and a scoring change forces a recompute
#' without a schema migration.
#'
#' @param con Open project DB connection.
#' @param id,path,scaffold Assembly unit.
#' @param force Recompute even when the stored row is current.
#' @return (invisibly) the verdict written, or the stored one when skipped.
#' @export
screen_one <- function(con, id, path = 1L, scaffold = 1L, force = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  cur <- DBI::dbGetQuery(con,
    "SELECT screen_verdict, screen_version FROM screen
     WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
  if (!isTRUE(force) && nrow(cur) == 1 &&
      isTRUE((cur$screen_version[1] %|NA|% -1L) >= SCREEN_VERSION)) {
    return(invisible(cur$screen_verdict[1]))
  }
  hits <- DBI::dbGetQuery(con,
    "SELECT hit_class, bitscore, lineage_taxids, lineage FROM blast_hits
     WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
  samp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids t FROM samples WHERE ID = ?", params = list(id))$t[1],
    error = function(e) NA_character_)
  r <- screen_channel_taxonomy(hits, samp)
  DBI::dbExecute(con,
    "DELETE FROM screen WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  DBI::dbAppendTable(con, "screen", data.frame(
    ID = id, path = path, scaffold = scaffold,
    screen_verdict = r$verdict, screen_class = r$class,
    screen_confidence = r$confidence, screen_reason = r$reason,
    screen_version = SCREEN_VERSION, lca_taxid = r$lca_taxid,
    lca_taxon = r$lca_taxon, lca_rank = r$lca_rank, lca_n_hits = r$lca_n_hits,
    host_concordant = r$host_concordant, margin = r$margin,
    time_stamp = as.integer(Sys.time()), stringsAsFactors = FALSE))
  invisible(r$verdict)
}

#' Screen every unit of one sample and refresh its reference support
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param force Recompute rows already at the current version.
#' @return (invisibly) the number of units screened.
#' @export
screen_sample <- function(con, id, force = FALSE) {
  u <- DBI::dbGetQuery(con,
    "SELECT path, scaffold FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
    params = list(id))
  for (i in seq_len(nrow(u))) screen_one(con, id, u$path[i], u$scaffold[i], force = force)
  update_reference_support(con, id)
  invisible(nrow(u))
}

#' Refresh a sample's reference support
#'
#' The minimum rank distance across the sample's ON-TARGET units only. Computed
#' from on-target units because a sample whose only resolvable hits belong to a
#' contaminant would otherwise report excellent support for the wrong organism.
#'
#' Recomputed fresh from each hit's `lineage_taxids` against the sample's own
#' lineage, rather than read from the stored `blast_hits.rank_distance` column:
#' that column is only ever populated at BLAST time when the sample lineage
#' was already resolved, and stays NULL for every hit otherwise, which would
#' silently starve this of a value for exactly the samples still awaiting
#' taxonomy resolution.
#' @noRd
update_reference_support <- function(con, id) {
  samp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids t FROM samples WHERE ID = ?", params = list(id))$t[1],
    error = function(e) NA_character_)
  h <- DBI::dbGetQuery(con,
    "SELECT h.lineage_taxids FROM blast_hits h
       JOIN screen s ON s.ID = h.ID AND s.path = h.path AND s.scaffold = h.scaffold
     WHERE h.ID = ? AND s.screen_verdict = 'on_target'", params = list(id))
  v <- NA_integer_
  if (nrow(h) > 0 && !is.na(samp %|NA|% NA) && nzchar(samp %|NA|% "")) {
    d <- vapply(h$lineage_taxids, function(x) {
      if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, samp))
    }, integer(1), USE.NAMES = FALSE)
    if (!all(is.na(d))) v <- min(d, na.rm = TRUE)
  }
  DBI::dbExecute(con, "UPDATE samples SET reference_support = ? WHERE ID = ?",
    params = list(v, id))
  invisible(v)
}
