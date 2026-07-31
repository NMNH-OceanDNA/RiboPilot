# Reference-guided curation of the nuclear rDNA feature model. Follows
# MitoPilot's design: a shared core (curate_rrna_core) driven by per-clade
# params (params_rrna.R), plus thin clade wrappers. The core is pure (no DB);
# curate_one() is the DB driver called by the Nextflow curate process,
# right after reconcile writes the annotations table.
#
# The core never base-aligns spacers. It checks the structural model (anchors
# present, ordered, one strand, sane lengths) and, when reference data is
# available, flags divergent anchors and contaminant lineages. Everything else
# is left to the human reviewer.

#' Curate one reconciled nrDNA feature model against a clade ruleset (pure)
#'
#' @param features data.frame(element, type, start, stop, length, strand, ...)
#'   as produced by [reconcile_features()] / stored in the `annotations` table.
#' @param params A ruleset from [rrna_params()].
#' @param ref_pid Optional named numeric vector (element -> percent identity to
#'   the top reference), from the `blast_ref_alignment` table. Missing elements
#'   are skipped (identity check is opportunistic).
#' @param ref_cov Optional named numeric vector (element -> fraction of the
#'   element covered by the reference alignment, `blast_ref_alignment.aln_cov`).
#'   The identity check is skipped when coverage is below `min_ref_coverage`, so
#'   a partial reference cannot raise a false low_identity flag.
#' @param sample_lineage,hit_lineage Optional lineage strings for the contaminant
#'   check: the sample's assigned lineage vs. the lineage of the top hit over
#'   *all* of core_nt (`blast_screen.lineage`), not the nuclear-rDNA reference
#'   panel. An unknown lineage on either side raises `contaminant_unchecked`
#'   rather than passing silently.
#' @param hit_is_nuclear_rrna Optional 0/1 (`blast_screen.is_nuclear_rrna`):
#'   whether that top hit is a nuclear rDNA record at all. 0 raises
#'   `not_nuclear_rdna`.
#' @param chimeric Optional 0/1/logical (`annotate.chimeric`): whether ITSx
#'   flagged the sequence as chimeric. TRUE raises `chimeric`.
#' @return list(features, flags, n_flags, notes) where `features` gains `status`
#'   ("ok"/"flag") and `flags` (comma-separated per-feature codes), and the
#'   sample-level `flags`/`notes` summarize missing anchors, ordering, etc.
#' @export
curate_rrna_core <- function(features, params = rrna_params(),
                             ref_pid = NULL, ref_cov = NULL,
                             sample_lineage = NA, hit_lineage = NA,
                             hit_is_nuclear_rrna = NA, chimeric = NA) {
  sample_flags <- character(0)
  feat_flags <- list()
  add_feat_flag <- function(el, code) feat_flags[[el]] <<- c(feat_flags[[el]], code)

  has_feats <- !is.null(features) && nrow(features) > 0
  present <- if (has_feats) features$element else character(0)

  # 1. Anchors + spacers present ----
  for (g in RRNA_GENES) if (!(g %in% present)) sample_flags <- c(sample_flags, paste0("missing_", g))
  for (s in RRNA_SPACERS) if (!(s %in% present)) sample_flags <- c(sample_flags, paste0("missing_", s))

  if (has_feats) {
    ord <- order(features$start)
    fo <- features[ord, , drop = FALSE]

    # 2. Order: observed elements must follow the canonical 5'->3' sequence ----
    canon <- RRNA_ELEMENT_ORDER[RRNA_ELEMENT_ORDER %in% fo$element]
    if (!identical(fo$element, canon)) sample_flags <- c(sample_flags, "out_of_order")

    # 3. Overlaps between adjacent features ----
    if (nrow(fo) > 1) {
      for (i in 2:nrow(fo)) {
        if (!is.na(fo$start[i]) && !is.na(fo$stop[i - 1]) && fo$start[i] <= fo$stop[i - 1]) {
          add_feat_flag(fo$element[i], "overlap")
          sample_flags <- c(sample_flags, "overlap")
        }
      }
    }

    # 4. Strand: everything should be on '+' after orientation ----
    if ("strand" %in% names(features)) {
      for (i in seq_len(nrow(features))) {
        st <- features$strand[i]
        if (!is.na(st) && !identical(st, "+")) {
          add_feat_flag(features$element[i], "wrong_strand")
          sample_flags <- c(sample_flags, "wrong_strand")
        }
      }
    }

    # 5. Per-feature length tolerance ----
    for (i in seq_len(nrow(features))) {
      el <- features$element[i]
      rng <- params$lengths[[el]]
      len <- features$length[i]
      if (is.null(rng) || is.na(len)) next
      if (len < rng[1] || len > rng[2]) add_feat_flag(el, "length")
    }

    # 6. Reference identity of conserved anchors (opportunistic) ----
    # Per-element identity floor, gated on the reference covering enough of the
    # element: a partial reference (the norm in sparse clades) covers little, so
    # its identity is not judged rather than flagged as low.
    if (!is.null(ref_pid) && length(ref_pid) > 0) {
      min_cov <- params$min_ref_coverage %||% 0.5
      for (g in RRNA_GENES) {
        pid <- if (g %in% names(ref_pid)) ref_pid[[g]] else NA_real_
        cov <- if (!is.null(ref_cov) && g %in% names(ref_cov)) ref_cov[[g]] else NA_real_
        floor_g <- anchor_pid_floor(params$min_anchor_pid, g)
        cov_ok <- is.na(cov) || cov >= min_cov
        if (!is.na(pid) && !is.na(floor_g) && cov_ok && pid < floor_g) {
          add_feat_flag(g, "low_identity")
        }
      }
    }

    # 6b. Anchor truncation (fraction of the rRNA model recovered). A short 18S
    # may be missing its 3' V9 region, which most metazoan eDNA assays target, so
    # a truncated anchor is flagged rather than exported as a complete gene.
    if ("model_coverage" %in% names(features)) {
      min_cov <- params$min_anchor_coverage %||% 90
      for (i in seq_len(nrow(features))) {
        cov <- features$model_coverage[i]
        if (!is.na(cov) && cov < min_cov) add_feat_flag(features$element[i], "truncated")
      }
    }
  }

  # 7. Contaminant checks (from the unrestricted GenBank screen) ----
  # `hit_lineage`/`hit_is_nuclear_rrna` describe the top hit over ALL of core_nt,
  # not the nuclear-rDNA reference panel. That is the point: only an unfiltered
  # hit can say "this is a symbiont" or "this is not rDNA at all".
  if (isTRUE(params$check_contaminant)) {
    if (lineages_conflict(sample_lineage, hit_lineage)) {
      sample_flags <- c(sample_flags, "contaminant")
    } else if (!lineage_known(sample_lineage) || !lineage_known(hit_lineage)) {
      # Not "clean" -- unknown. Absence of evidence must not read as evidence of
      # absence, or an offline/unresolved sample looks verified.
      sample_flags <- c(sample_flags, "contaminant_unchecked")
    }
    if (identical(hit_is_nuclear_rrna, 0L) || identical(hit_is_nuclear_rrna, FALSE)) {
      # Best match in all of GenBank is not a nuclear rDNA record, so the
      # assembly is probably not the target locus (organellar rRNA, bacterial
      # 16S, a genomic repeat, ...).
      sample_flags <- c(sample_flags, "not_nuclear_rdna")
    }
  }

  # 8. ITSx chimera call (detected at annotation; a chimeric sequence must not
  # enter a reference database) ----
  if (isTRUE(chimeric) || identical(chimeric, 1L)) {
    sample_flags <- c(sample_flags, "chimeric")
  }

  # Assemble per-feature status ----
  if (has_feats) {
    features$flags <- vapply(features$element, function(el) {
      fl <- feat_flags[[el]]
      if (is.null(fl)) "" else paste(unique(fl), collapse = ",")
    }, character(1))
    features$status <- ifelse(nzchar(features$flags), "flag", "ok")
  } else {
    features <- data.frame()
  }

  sample_flags <- unique(sample_flags)
  notes <- if (length(sample_flags)) paste(sample_flags, collapse = "; ") else NA_character_
  list(features = features, flags = sample_flags, n_flags = length(sample_flags), notes = notes)
}

#' Thin actinopterygii curation wrapper.
#' @param features,ref_pid,ref_cov,sample_lineage,hit_lineage,hit_is_nuclear_rrna,chimeric see
#'   [curate_rrna_core()].
#' @return see [curate_rrna_core()].
#' @export
curate_actinopterygii_rrna <- function(features, ref_pid = NULL, ref_cov = NULL,
                                       sample_lineage = NA, hit_lineage = NA,
                                       hit_is_nuclear_rrna = NA, chimeric = NA) {
  curate_rrna_core(features, rrna_params("actinopterygii"),
    ref_pid = ref_pid, ref_cov = ref_cov, sample_lineage = sample_lineage,
    hit_lineage = hit_lineage, hit_is_nuclear_rrna = hit_is_nuclear_rrna,
    chimeric = chimeric)
}

#' Is a lineage string usable for the contaminant comparison?
#' @noRd
lineage_known <- function(x) {
  length(x) > 0 && !is.na(x[1]) && nzchar(trimws(x[1]))
}

#' Do two lineage strings share no taxon in common (gross contaminant)?
#' Conservative: only TRUE when both are present and, after dropping ubiquitous
#' high-rank tokens, they share zero taxa. Catches wrong-kingdom hits (e.g. a
#' fungal hit on a fish) without flagging same-clade samples.
#' @noRd
LINEAGE_STOPWORDS <- c("cellular organisms", "eukaryota", "opisthokonta")
lineages_conflict <- function(a, b) {
  toks <- function(x) {
    if (length(x) == 0 || is.na(x) || !nzchar(x)) return(character(0))
    t <- tolower(trimws(unlist(strsplit(x, "[;,]"))))
    unique(setdiff(t[nchar(t) > 3], LINEAGE_STOPWORDS))
  }
  ta <- toks(a); tb <- toks(b)
  if (length(ta) == 0 || length(tb) == 0) return(FALSE)
  length(intersect(ta, tb)) == 0
}

#' Ensure the curation columns/tables exist (idempotent)
#'
#' Thin wrapper over [migrate_db()], kept for call sites that only need the
#' curation schema (e.g. [curate_one()]). New projects get these in [new_db()].
#' @param con Open DBI connection.
#' @export
ensure_curate_schema <- function(con) {
  migrate_db(con)
  invisible()
}

#' Resolve which clade ruleset governs a unit
#'
#' Precedence: explicit argument > the unit's `curate_opts` row. When that is
#' unset or the generic default, auto-select a clade ruleset from the sample's
#' lineage so a clade sample is not silently curated against the generic ranges.
#' If a lineage is present but matches no clade ruleset, `note` is
#' `no_clade_ruleset` so the fallback is visible rather than silent.
#'
#' Shared by [curate_one()] and the screen, which needs the same answer to reach
#' a clade's `screen` tunables.
#' @return list(ruleset, note)
#' @noRd
resolve_ruleset <- function(con, id, path = 1L, scaffold = 1L, ruleset = NULL) {
  sample_lineage <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage FROM samples WHERE ID = ?", params = list(id))$lineage[1],
    error = function(e) NA_character_)
  if (is.null(ruleset)) {
    ruleset <- tryCatch(DBI::dbGetQuery(con,
      "SELECT co.ruleset FROM annotate an JOIN curate_opts co ON an.curate_opts = co.curate_opts
       WHERE an.ID = ? AND an.path = ? AND an.scaffold = ?",
      params = list(id, path, scaffold))$ruleset[1], error = function(e) NA_character_)
  }
  note <- NA_character_
  if (length(ruleset) == 0 || is.na(ruleset) || !nzchar(ruleset) ||
      identical(ruleset, "metazoa_generic")) {
    ruleset <- ruleset_for_lineage(sample_lineage)
    if (identical(ruleset, "metazoa_generic") && length(sample_lineage) &&
        !is.na(sample_lineage) && nzchar(trimws(sample_lineage))) {
      note <- "no_clade_ruleset"
    }
  }
  list(ruleset = ruleset, note = note)
}

#' Curate one sample's annotations and write results back to the DB
#'
#' Reads the reconciled `annotations`, the sample's `curate_opts` ruleset, any
#' per-feature reference identities (`blast_ref_alignment`), and the sample vs.
#' hit lineages, runs [curate_rrna_core()], then updates each feature's
#' `status`/`flags` and the per-sample `annotate` summary (`n_flags`,
#' `annotate_notes`) in one transaction.
#'
#' @param con Open DBI connection.
#' @param id Sample ID.
#' @param path,scaffold Assembly unit within the sample.
#' @param ruleset Optional ruleset override (else read from `curate_opts`).
#' @return (invisibly) the curate_rrna_core() result.
#' @export
curate_one <- function(con, id, path = 1L, scaffold = 1L, ruleset = NULL) {
  ensure_curate_schema(con)
  path <- as.integer(path); scaffold <- as.integer(scaffold)

  feats <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop, length, strand, model_coverage FROM annotations
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))

  sample_lineage <- DBI::dbGetQuery(con, "SELECT lineage FROM samples WHERE ID = ?",
    params = list(id))$lineage[1]

  rs <- resolve_ruleset(con, id, path, scaffold, ruleset)
  ruleset <- rs$ruleset
  ruleset_note <- rs$note
  params <- rrna_params(ruleset)

  # Per-feature reference identity + coverage (opportunistic; empty when
  # BLAST/ref absent).
  ref_pid <- NULL
  ref_cov <- NULL
  rp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT element, pid, aln_cov FROM blast_ref_alignment WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  if (nrow(rp) > 0) {
    ref_pid <- stats::setNames(rp$pid, rp$element)
    if ("aln_cov" %in% names(rp)) ref_cov <- stats::setNames(rp$aln_cov, rp$element)
  }
  # The contaminant check reads the unrestricted GenBank screen, NOT
  # path_scores.blast_lineage. The latter comes from BLAST against the clade
  # seed DB, so its lineage is host-clade by construction and can never conflict
  # with the sample's -- the check could not fire.
  scr <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage, is_nuclear_rrna FROM blast_screen WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  hit_lineage <- if (nrow(scr) > 0) scr$lineage[1] else NA_character_
  hit_is_nuc <- if (nrow(scr) > 0) as.integer(scr$is_nuclear_rrna[1]) else NA_integer_

  chimeric <- tryCatch(DBI::dbGetQuery(con,
    "SELECT chimeric FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$chimeric[1],
    error = function(e) NA_integer_)

  res <- curate_rrna_core(feats, params, ref_pid = ref_pid, ref_cov = ref_cov,
    sample_lineage = sample_lineage, hit_lineage = hit_lineage,
    hit_is_nuclear_rrna = hit_is_nuc, chimeric = chimeric)

  if (!is.na(ruleset_note)) {
    res$flags <- unique(c(res$flags, ruleset_note))
    res$n_flags <- length(res$flags)
    res$notes <- if (length(res$flags)) paste(res$flags, collapse = "; ") else NA_character_
  }

  ts <- as.integer(Sys.time())
  DBI::dbWithTransaction(con, {
    if (nrow(res$features) > 0) {
      for (i in seq_len(nrow(res$features))) {
        DBI::dbExecute(con,
          "UPDATE annotations SET status = ?, flags = ? WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
          params = list(res$features$status[i],
            if (nzchar(res$features$flags[i])) res$features$flags[i] else NA_character_,
            id, path, scaffold, res$features$element[i]))
      }
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET n_flags = ?, annotate_notes = ?, time_stamp = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(res$n_flags, res$notes, ts, id, path, scaffold))
  })
  invisible(res)
}
