# Per-scaffold triage: keep / drop / undecided, with a reason and an owner.
#
# `assemblies.ignore` remains the single downstream gate -- WF2's unit query, the
# Annotate table and export_bundle() all filter on it -- but it is a binary with
# no "not looked at yet" state, no record of why, and no way to tell a machine's
# guess from the curator's decision. This module owns the decision:
#
#   triage = 'keep'  -> ignore = 0     included, reviewed
#   triage = 'drop'  -> ignore = 1     excluded
#   triage = NULL    -> ignore = 0     UNDECIDED: included, but flagged
#
# Undecided deliberately still flows downstream. The cost of wrongly dropping a
# scaffold is a lost record; the cost of wrongly keeping one is an extra row the
# curator looks at. Nothing is ever excluded silently.
#
# triage_source distinguishes 'auto' from 'user' and a user decision is sticky:
# re-running auto-triage never overwrites one.

#' Valid triage verdicts. NA / NULL is the third state (undecided).
#' @noRd
TRIAGE_VERDICTS <- c("keep", "drop")

#' Record a triage decision for one assembly unit
#'
#' Writes the verdict, its justification and its owner, derives
#' `assemblies.ignore` from the verdict, and reconciles the unit's annotate row
#' (held on drop, restored on keep) through [reconcile_ignore_unit()].
#'
#' @param con Open project DB connection.
#' @param id,path,scaffold Assembly unit.
#' @param verdict "keep", "drop", or NA/NULL for undecided.
#' @param reason Short plain-English justification, shown to the curator.
#' @param source "user" or "auto".
#' @return (invisibly) the derived ignore flag.
#' @export
set_unit_triage <- function(con, id, path, scaffold, verdict,
                            reason = NA_character_, source = c("user", "auto")) {
  source <- match.arg(source)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  if (length(verdict) == 0 || is.na(verdict[1]) || !nzchar(verdict[1])) {
    verdict <- NA_character_
  } else if (verdict[1] %in% TRIAGE_VERDICTS) {
    # exact, not match.arg(): partial matching would accept "k"/"d" and then
    # report the abbreviation back to the caller as if it were the verdict
    verdict <- verdict[1]
  } else {
    stop("verdict must be one of: ", paste(TRIAGE_VERDICTS, collapse = ", "),
      ", or NA for undecided", call. = FALSE)
  }
  ig <- as.integer(identical(verdict, "drop"))
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage = ?, triage_reason = ?, triage_source = ?, ignore = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(verdict, reason %||% NA_character_, source, ig, id, path, scaffold))
  reconcile_ignore_unit(con, id, path, scaffold, ig)
  invisible(ig)
}

#' Current triage state of a sample's units, with the evidence auto-triage reads
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return data.frame, one row per unit; empty when the sample has none.
#' @noRd
unit_triage_state <- function(con, id) {
  DBI::dbGetQuery(con,
    "SELECT a.path, a.scaffold, a.length, COALESCE(a.ignore,0) AS ignore,
       a.triage, a.triage_reason, a.triage_source,
       COALESCE(c.median_depth, c.mean_depth) AS depth,
       b.species AS screen_species, b.lineage_taxids AS screen_taxids,
       b.rank_distance, b.is_nuclear_rrna,
       s.lineage_taxids AS sample_taxids, s.sci_name,
       (SELECT COUNT(*) FROM annotations n
          WHERE n.ID = a.ID AND n.path = a.path AND n.scaffold = a.scaffold) AS n_features,
       (SELECT COALESCE(an.annotate_lock, 0) FROM annotate an
          WHERE an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold) AS annotate_lock,
       (SELECT COALESCE(an.reviewed, 'no') FROM annotate an
          WHERE an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold) AS reviewed
     FROM assemblies a
     LEFT JOIN coverage c     ON c.ID = a.ID AND c.path = a.path AND c.scaffold = a.scaffold
     LEFT JOIN blast_screen b ON b.ID = a.ID AND b.path = a.path AND b.scaffold = a.scaffold
     LEFT JOIN samples s      ON s.ID = a.ID
     WHERE a.ID = ? ORDER BY a.path, a.scaffold", params = list(id))
}

#' Default auto-triage settings
#'
#' `min_depth_frac` only ever raises a flag, never drops: within-sample depth is
#' diluted by multi-mapping between near-identical isomer paths, so a low ratio
#' is a reason to look rather than a reason to exclude.
#' @noRd
TRIAGE_DEFAULTS <- list(min_depth_frac = 0.05)

#' Propose a triage verdict for every unit of one sample
#'
#' Pure: takes the state frame, returns verdicts. The rules are deliberately
#' lopsided. A DROP requires POSITIVE evidence that the scaffold is not the
#' target locus; everything suspicious but unproven becomes UNDECIDED, which
#' still flows downstream and merely asks to be looked at.
#'
#' In particular nothing is dropped for being taxonomically distant. Low identity
#' at a large rank distance is the normal state for an under-represented taxon --
#' the reason the tool exists -- and dropping on it would fail exactly the samples
#' it is meant to serve. Distance only becomes evidence when paired with an
#' implausibly high identity, which needs the identity envelope (not yet built),
#' so for now it can only raise a flag.
#'
#' @param d Frame from [unit_triage_state()].
#' @param opts See [TRIAGE_DEFAULTS].
#' @return `d` with `new_triage` and `new_reason` columns.
#' @noRd
triage_rules <- function(d, opts = TRIAGE_DEFAULTS) {
  n <- nrow(d)
  out <- rep(NA_character_, n)
  why <- rep(NA_character_, n)
  if (n == 0) {
    d$new_triage <- out; d$new_reason <- why
    return(d)
  }
  # depth is comparable only within a sample, and only where coverage ran
  dmax <- suppressWarnings(max(d$depth, na.rm = TRUE))
  frac <- if (is.finite(dmax) && dmax > 0) d$depth / dmax else rep(NA_real_, n)

  sib <- unit_concordance(data.frame(
    screen_taxids = d$screen_taxids %||% rep(NA_character_, n),
    screen_lineage = rep(NA_character_, n), stringsAsFactors = FALSE))

  for (i in seq_len(n)) {
    feats <- d$n_features[i] %||% 0
    nuc <- d$is_nuclear_rrna[i]
    has_screen <- !is.na(d$screen_species[i] %||% NA)

    if (isTRUE(feats == 0) && isTRUE(nuc == 0)) {
      # POSITIVE evidence on both channels available today: the annotation tools
      # found no rRNA or ITS feature at all, AND the best match anywhere in
      # core_nt is not a nuclear rDNA record. Neither alone would justify a drop.
      out[i] <- "drop"
      why[i] <- sprintf("no rRNA or ITS features found, and the best GenBank match (%s) is not a nuclear rDNA record",
        d$screen_species[i] %|NA|% "unknown")
    } else if (isTRUE(feats == 0)) {
      why[i] <- if (has_screen) {
        sprintf("no rRNA or ITS features found; best GenBank match is %s", d$screen_species[i] %|NA|% "unknown")
      } else {
        "no rRNA or ITS features found, and the GenBank screen has not run"
      }
    } else if (isTRUE(nuc == 0)) {
      # Has a feature model, but the best match anywhere in core_nt is not a
      # nuclear rDNA record. Not enough to drop on its own -- is_nuclear_rrna is
      # decided by a title regex over one hit -- but it is the single strongest
      # "this may not be the target locus" signal available, so it must not pass
      # silently.
      why[i] <- sprintf("best GenBank match (%s) is not a nuclear rDNA record",
        d$screen_species[i] %|NA|% "unknown")
    } else if (isTRUE(sib$odd_one_out[i])) {
      why[i] <- sprintf("top GenBank hit (%s) differs from this sample's other scaffolds (%s at %s)",
        d$screen_species[i] %|NA|% "unknown", sib$sib_taxon[i] %|NA|% "the majority",
        sib$sib_rank[i] %|NA|% "an agreed rank")
    } else if (!is.na(frac[i]) && frac[i] < (opts$min_depth_frac %||% 0.05)) {
      why[i] <- sprintf("median depth %sx is %d%% of this sample's deepest scaffold",
        round(d$depth[i], 1), as.integer(round(100 * frac[i])))
    } else if (!has_screen) {
      why[i] <- "GenBank screen has not run for this scaffold"
    } else {
      out[i] <- "keep"
    }
  }
  d$new_triage <- out
  d$new_reason <- why
  d
}

#' Apply automatic triage to one sample
#'
#' Skips any unit a human has already decided (`triage_source = 'user'`) and any
#' unit whose annotation is locked, so re-running is safe.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param opts See [TRIAGE_DEFAULTS].
#' @param overwrite_user Also re-decide units a curator set by hand.
#' @return (invisibly) the number of units changed.
#' @export
auto_triage_sample <- function(con, id, opts = TRIAGE_DEFAULTS, overwrite_user = FALSE) {
  d <- unit_triage_state(con, id)
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  d <- triage_rules(d, opts)
  changed <- 0L
  for (i in seq_len(nrow(d))) {
    if (!isTRUE(overwrite_user) && identical(d$triage_source[i], "user")) next
    if (isTRUE(d$annotate_lock[i] == 1)) next
    same <- identical(d$triage[i] %|NA|% "", d$new_triage[i] %|NA|% "") &&
      identical(d$triage_reason[i] %|NA|% "", d$new_reason[i] %|NA|% "")
    if (same) next
    set_unit_triage(con, id, d$path[i], d$scaffold[i], d$new_triage[i],
      reason = d$new_reason[i], source = "auto")
    changed <- changed + 1L
  }
  invisible(changed)
}

#' What would bulk triage do, without doing it
#'
#' A bulk write across a whole project is the highest-blast-radius action in the
#' app, so it is never applied without showing the count first -- including how
#' many samples would end up with NOTHING kept, which is the outcome that
#' silently loses a sample.
#'
#' @param con Open project DB connection.
#' @param ids Sample IDs to consider.
#' @param opts See [TRIAGE_DEFAULTS].
#' @param overwrite_user Include units a curator decided by hand.
#' @return list(n_units, n_drop, n_undecided, n_keep, n_samples, emptied,
#'   n_user_skipped, n_locked_skipped).
#' @export
triage_preview <- function(con, ids, opts = TRIAGE_DEFAULTS, overwrite_user = FALSE) {
  res <- list(n_units = 0L, n_drop = 0L, n_undecided = 0L, n_keep = 0L,
    n_samples = 0L, emptied = character(0), n_user_skipped = 0L, n_locked_skipped = 0L)
  for (id in unique(ids)) {
    d <- unit_triage_state(con, id)
    if (nrow(d) == 0) next
    d <- triage_rules(d, opts)
    skip_user <- !isTRUE(overwrite_user) & identical_user(d$triage_source)
    skip_lock <- !is.na(d$annotate_lock) & d$annotate_lock == 1
    res$n_user_skipped <- res$n_user_skipped + sum(skip_user & !skip_lock)
    res$n_locked_skipped <- res$n_locked_skipped + sum(skip_lock)
    # a skipped unit keeps whatever it already had
    eff <- ifelse(skip_user | skip_lock, d$triage, d$new_triage)
    res$n_units <- res$n_units + nrow(d)
    res$n_drop <- res$n_drop + sum(eff %in% "drop")
    res$n_keep <- res$n_keep + sum(eff %in% "keep")
    res$n_undecided <- res$n_undecided + sum(is.na(eff))
    res$n_samples <- res$n_samples + 1L
    if (all(eff %in% "drop")) res$emptied <- c(res$emptied, id)
  }
  res
}

#' Vectorised "is this row owned by a user"
#' @noRd
identical_user <- function(x) !is.na(x) & x == "user"
