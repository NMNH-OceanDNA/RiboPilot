# The screen's taxonomy channel: does this unit's hit distribution agree on a
# taxon, and is that taxon the host's.
#
# It uses NO identity threshold. The identity envelope was deferred after
# measurement against the reference project (doc 08 section 3.2.1): that corpus
# has 2 hits at phylum distance and none beyond, so it cannot calibrate a bound
# acting there, and genuine host 18S reaches 99.6 percent at class distance, so
# no ceiling can sit below that without flagging real host data.
#
# What this cannot do: separate two metazoans within the same phylum. Concordance
# is tested at phylum rank, so a copepod on a fish is caught (Arthropoda vs
# Chordata) but a human on a fish is not (both Chordata). Doc 08 section 3.5 puts
# those cases in the depth and cross-sample channels.

#' Scoring version. Bump to force a recompute of stored rows.
#' @noRd
SCREEN_VERSION <- 5L

#' @noRd
SCREEN_VERDICTS <- c("on_target", "off_target", "no_reference", "uncertain")

#' @noRd
SCREEN_CLASSES <- c("bacterial", "archaeal", "fungal", "algal_plant",
  "other_euk", "metazoan_other", "human")

#' @noRd
SCREEN_NO_REF_CODES <- c("no_hits", "no_taxdump", "unresolved_sample_lineage",
  "unresolved_hit_lineages")

#' Confidence weights, and the gate below which a verdict is uncertain.
#'
#' Both are views on [SCREEN_PARAMS_DEFAULT] (R/params_rrna.R) so a clade
#' ruleset and the shipped default cannot drift apart. Any change must bump
#' [SCREEN_VERSION].
#' @noRd
SCREEN_WEIGHTS <- SCREEN_PARAMS_DEFAULT[c("support", "resolution", "separation")]

#' @noRd
SCREEN_MIN_CONFIDENCE <- SCREEN_PARAMS_DEFAULT$min_confidence

#' The screen tunables for a sample, resolved through its curation ruleset
#'
#' Falls back to the default rather than erroring: a screen must not fail
#' because a project predates the ruleset columns.
#' @noRd
screen_params <- function(con, id) {
  rs <- tryCatch(resolve_ruleset(con, id, 1L, 1L)$ruleset, error = function(e) NULL)
  p <- tryCatch(rrna_params(rs)$screen, error = function(e) NULL)
  p %||% SCREEN_PARAMS_DEFAULT
}

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
#'
#' It is also worth 0.40 of the 0.60 gate on its own, which is why the caller
#' only claims it for at least three DISTINCT accessions: two rows agreeing at
#' species rank would otherwise clear the gate at 0.76 on nothing else.
#' @noRd
screen_confidence_terms <- function(lca_n_hits, lca_rank, margin, top_bitscore,
                                    unanimous = FALSE,
                                    params = SCREEN_PARAMS_DEFAULT) {
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
  w <- params %||% SCREEN_PARAMS_DEFAULT
  list(support = support, resolution = resolution, separation = separation,
    confidence = unname(support * w$support + resolution * w$resolution +
      separation * w$separation))
}

#' Well-known kingdom/domain-rank taxids, stable across taxdump versions.
#'
#' Classification keys off the taxid rather than the kingdom-rank NAME:
#' taxids are immune to synonyms, homonyms and defline formatting, exactly the
#' reason [rank_distance()] compares taxids and treats names as display-only.
#'
#' NCBI's 2024 prokaryote reclassification put a kingdom rank under Bacteria and
#' Archaea, so `{K}` resolves to one of these rather than to 2 / 2157 and a
#' bacterial contaminant would otherwise fall through to "other_euk". Both the
#' domains and their kingdoms are listed, so either taxdump generation classifies.
#' @noRd
SCREEN_BACTERIAL_TAXIDS <- c("2", "1783272", "3379134", "3384189", "3384194")
SCREEN_ARCHAEAL_TAXIDS <- c("2157", "1783275", "1783276", "1935183", "3366610")

SCREEN_KINGDOM_TAXIDS <- c(
  stats::setNames(rep("bacterial", length(SCREEN_BACTERIAL_TAXIDS)), SCREEN_BACTERIAL_TAXIDS),
  stats::setNames(rep("archaeal", length(SCREEN_ARCHAEAL_TAXIDS)), SCREEN_ARCHAEAL_TAXIDS),
  "4751" = "fungal", "33090" = "algal_plant", "33208" = "metazoan_other")

#' Off-target subclass from the LCA's rank-aligned taxid path and its taxon label
#'
#' Matches the WHOLE path rather than a single rank slot: whether the taxdump's
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
  for (k in c("4751", "33090", "33208", SCREEN_BACTERIAL_TAXIDS, SCREEN_ARCHAEAL_TAXIDS)) {
    if (as.integer(k) %in% p) return(unname(SCREEN_KINGDOM_TAXIDS[k]))
  }
  "other_euk"
}

#' One vote per accession: the best-bitscore HSP of each subject
#'
#' `blast_hits` stores one row per HSP and blastn runs with `-max_hsps 5`, so a
#' single accession can otherwise cast up to five votes and alone satisfy both
#' [bitscore_weighted_lca()]'s `min_hits` and the unanimity floor.
#' `parse_blast_candidates()` assigns `hit_rank` per subject, so one `hit_rank`
#' is one accession; a frame without the column is already one row per subject.
#' @noRd
screen_one_vote_per_subject <- function(hits) {
  if (is.null(hits) || nrow(hits) == 0) return(hits)
  hr <- if (!is.null(hits$hit_rank)) as.integer(hits$hit_rank) else seq_len(nrow(hits))
  # an unranked row cannot be pooled with any other, so give it its own key
  if (any(is.na(hr))) hr[is.na(hr)] <- -seq_len(sum(is.na(hr)))
  bs <- hits$bitscore %||% rep(NA_real_, nrow(hits))
  bs[is.na(bs)] <- -Inf
  ord <- order(hr, -bs)
  hits[ord, , drop = FALSE][!duplicated(hr[ord]), , drop = FALSE]
}

#' Largest value, or NA for an empty or all-NA vector
#'
#' `max(na.rm = TRUE)` over an all-NA side returns -Inf and warns into the
#' pipeline log, which would make the margin infinite and force separation to 1.
#' @noRd
screen_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) max(x) else NA_real_
}

#' Score one unit's hit distribution against its sample's lineage
#'
#' @param hits Frame with `hit_class`, `bitscore`, `lineage_taxids`, `lineage`
#'   and optionally `hit_rank` (one rank per subject accession).
#' @param sample_taxids The sample's rank-aligned taxid path.
#' @param params The `screen` sublist of a ruleset (see [rrna_params()]).
#' @return list(verdict, class, confidence, reason, reason_code, terms,
#'   lca_taxid, lca_taxon, lca_rank, lca_n_hits, host_concordant, margin,
#'   hit_class).
#' @noRd
screen_channel_taxonomy <- function(hits, sample_taxids,
                                    params = SCREEN_PARAMS_DEFAULT) {
  params <- params %||% SCREEN_PARAMS_DEFAULT
  out <- list(verdict = "no_reference", class = NA_character_, confidence = 0,
    reason = NA_character_, reason_code = "no_hits",
    terms = list(support = 0, resolution = 0, separation = 0),
    lca_taxid = NA_integer_, lca_taxon = NA_character_, lca_rank = NA_character_,
    lca_n_hits = 0L, host_concordant = NA_integer_, margin = NA_real_,
    hit_class = NA_character_)

  samp <- parse_rank_taxids(sample_taxids %||% NA_character_)
  # concordance is tested at phylum rank or deeper, so a path resolved only at
  # the shallowest slot leaves nothing to compare. That is missing evidence, not
  # contrary evidence: read as discordant it makes every hit unanimous and
  # scores clean host data off_target.
  samp_ok <- samp[-1]
  samp_ok <- samp_ok[!is.na(samp_ok)]
  if (!length(samp_ok)) {
    out$reason_code <- "unresolved_sample_lineage"
    out$reason <- "the sample's own lineage is unresolved below the shallowest rank"
    return(out)
  }
  if (is.null(hits) || nrow(hits) == 0) {
    out$reason <- "no GenBank hits are stored for this scaffold"
    return(out)
  }
  hits <- screen_one_vote_per_subject(hits)

  # uncultured/environmental records carry their junk label at the deepest
  # FILLED rank (a placeholder still sits in the species slot), which every
  # rank rule would otherwise read as a real species; a merely blank rank
  # (nothing filled in yet) is normal and must not be excluded on that basis
  # alone, so the check looks at the deepest non-empty raw field, not at
  # whichever rank position happens to be blank.
  deepest_raw <- vapply(hits$lineage, function(x) {
    f <- trimws(strsplit(x %|NA|% "", ";", fixed = TRUE)[[1]])
    f <- f[nzchar(f)]
    if (length(f) == 0) NA_character_ else f[length(f)]
  }, character(1), USE.NAMES = FALSE)
  # `other` votes alongside `rdna`: a bacterial 16S operon or an organellar
  # rRNA contig is POSITIVE evidence about what the scaffold is, and excluding
  # it left the unit with no hits at all, i.e. reported as unjudged rather than
  # as off-target. `uninformative` stays out, which is what that class is for:
  # a chromosome-level record says nothing about the aligned region, so a unit
  # whose only hits are assembly records is still genuinely unscreened.
  keep <- hits$hit_class %in% c("rdna", "other") &
    !is.na(deepest_raw) & !is_placeholder_taxon(deepest_raw)
  h <- hits[keep, , drop = FALSE]
  if (nrow(h) == 0) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "no hit resolved to a usable lineage"
    return(out)
  }
  # which kind of record carried the verdict, so a curator can see whether it
  # rested on nuclear rDNA or on something else entirely
  hbs <- h$bitscore
  hbs[is.na(hbs)] <- 0
  out$hit_class <- if (sum(hbs[h$hit_class %in% "other"]) >
                       sum(hbs[h$hit_class %in% "rdna"])) "other" else "rdna"

  lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
  if (is.na(lca$taxid)) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "the hits agree on no shared taxon"
    out$lca_n_hits <- lca$n_hits
    return(out)
  }
  out$lca_taxid <- lca$taxid; out$lca_taxon <- lca$taxon
  out$lca_rank <- lca$rank;  out$lca_n_hits <- lca$n_hits

  # host_concordant asks about the unit's SINGLE weighted LCA: does that one
  # taxon sit on the sample's path at phylum rank or deeper. The LCA's own
  # ancestry is the prefix of any hit that voted for it, truncated at its rank.
  #
  # It must not test `lca$taxid %in% samp` alone. That is the too-strict rule
  # the margin already discarded: three hits of one congener put the LCA at the
  # congener's species, which is absent from the sample's path, so a clean unit
  # scores off_target. It fires hardest on under-described taxa.
  lca_pos <- match(lca$rank, RRNA_LINEAGE_RANKS)
  hit_paths <- lapply(h$lineage_taxids, parse_rank_taxids)
  donor <- Find(function(p) isTRUE(p[lca_pos] == lca$taxid), hit_paths)
  lca_path <- if (is.null(donor)) integer(0) else donor[seq_len(lca_pos)][-1]
  lca_path <- lca_path[!is.na(lca_path)]
  out$host_concordant <- as.integer(any(lca_path %in% samp_ok))

  # the margin uses a SECOND, per-hit concordance, distinct from the unit-level
  # one above: it compares the two best-supported opposing hits, asking the same
  # question of each hit individually rather than of the single LCA.
  # a hit is on the sample's lineage when their paths share a resolved taxid at
  # PHYLUM or deeper. Sharing only the shallowest slot means "both are animals",
  # which every metazoan contaminant satisfies: three host hits plus two from a
  # different phylum would otherwise read as unanimous and score on_target.
  hit_conc <- vapply(hit_paths, function(p) {
    p <- p[-1]
    any(p[!is.na(p)] %in% samp_ok)
  }, logical(1), USE.NAMES = FALSE)
  out$margin <- screen_max(h$bitscore[hit_conc]) - screen_max(h$bitscore[!hit_conc])

  # every surviving hit on one side means nothing opposes the winner, which is
  # maximal separation rather than missing separation. The floor is three
  # DISTINCT accessions, matching bitscore_weighted_lca()'s own min_hits: one
  # row is unopposed by definition, and two rows agreeing at species rank would
  # clear the 0.60 gate at 0.76 on unanimity alone.
  unanimous <- (all(hit_conc) || !any(hit_conc)) && lca$n_hits >= 3
  top_bitscore <- screen_max(h$bitscore)
  tm <- screen_confidence_terms(lca$n_hits, lca$rank, out$margin,
    top_bitscore, unanimous = unanimous, params = params)
  out$terms <- tm[c("support", "resolution", "separation")]
  out$confidence <- tm$confidence
  out$reason_code <- NA_character_

  out$verdict <- if (tm$confidence < (params$min_confidence %||% SCREEN_MIN_CONFIDENCE)) {
    "uncertain"
  } else if (identical(out$host_concordant, 1L)) {
    "on_target"
  } else {
    "off_target"
  }
  if (identical(out$verdict, "off_target")) {
    # the best-supported surviving hit, not whichever row the query returned
    # first; and the whole path, not one rank slot: see screen_class_for()
    bs <- h$bitscore
    bs[is.na(bs)] <- -Inf
    out$class <- screen_class_for(h$lineage_taxids[which.max(bs)], lca$taxon)
  }
  out$reason <- sprintf(
    "hits agree on %s (%s, %d voting, mostly %s records); %s the sample's lineage. support %.2f w %.2f, resolution %.2f w %.2f, separation %.2f w %.2f -> confidence %.2f",
    lca$taxon %|NA|% "an unnamed taxon", lca$rank, lca$n_hits,
    if (identical(out$hit_class, "other")) "non-rDNA" else "nuclear rDNA",
    if (identical(out$host_concordant, 1L)) "on" else "NOT on",
    tm$support, params$support, tm$resolution, params$resolution,
    tm$separation, params$separation, tm$confidence)
  out
}

#' Score and store one unit's screen verdict
#'
#' Skips a row already stored at the current [SCREEN_VERSION] unless `force`, so
#' re-running the pipeline is cheap and a scoring change forces a recompute
#' without a schema migration.
#'
#' A stored `no_reference` is always recomputed: it records missing evidence,
#' and the evidence routinely arrives later. A unit screened before BLAST ran,
#' or before the sample's lineage was resolved, would otherwise keep that
#' verdict for the life of the project.
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
      isTRUE((cur$screen_version[1] %|NA|% -1L) >= SCREEN_VERSION) &&
      !identical(cur$screen_verdict[1] %|NA|% "", "no_reference")) {
    return(invisible(cur$screen_verdict[1]))
  }
  hits <- DBI::dbGetQuery(con,
    "SELECT hit_rank, hit_class, bitscore, lineage_taxids, lineage FROM blast_hits
     WHERE ID = ? AND path = ? AND scaffold = ?
     ORDER BY hit_rank, hsp", params = list(id, path, scaffold))
  samp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids t FROM samples WHERE ID = ?", params = list(id))$t[1],
    error = function(e) NA_character_)
  r <- screen_channel_taxonomy(hits, samp, params = screen_params(con, id))
  # unresolved_sample_lineage and no_taxdump are excluded on purpose: those are
  # about the sample's own taxonomy or a missing taxdump, and a remote search
  # would not help either one.
  if (r$reason_code %in% c("no_hits", "unresolved_hit_lineages")) {
    # a missing row means remote.
    searched_source <- DBI::dbGetQuery(con,
      "SELECT source FROM blast_searched WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, path, scaffold))$source
    if (length(searched_source) > 0 && identical(searched_source[1], "local")) {
      r$reason <- if (identical(r$reason_code, "no_hits")) {
        paste("the local reference panel returned no usable hit for this",
          "scaffold; a remote BLAST search may identify it")
      } else {
        paste("the local reference panel returned no hit that resolved to a",
          "usable lineage; a remote BLAST search may identify it")
      }
    }
  }
  # one transaction: a failed append after a committed delete would leave the
  # unit with no verdict at all
  DBI::dbWithTransaction(con, {
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
  })
  invisible(r$verdict)
}

#' Screen every unit of one sample
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param force Recompute rows already at the current version.
#' @param variants Also classify on-target units as isomers or paralogs. The pipeline's screen step passes FALSE; its triage step runs that check once annotation is done.
#' @return (invisibly) the number of units screened.
#' @export
screen_sample <- function(con, id, force = FALSE, variants = TRUE) {
  u <- DBI::dbGetQuery(con,
    "SELECT path, scaffold FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
    params = list(id))
  for (i in seq_len(nrow(u))) screen_one(con, id, u$path[i], u$scaffold[i], force = force)
  # after the loop, never before: screen_one() rewrites a row with DELETE plus
  # append, so a variant column written first would be silently wiped
  if (isTRUE(variants)) screen_variants(con, id)
  invisible(nrow(u))
}
