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
#' @param rule Key of the rule that fired, or NA for a user decision.
#' @return (invisibly) the derived ignore flag.
#' @export
set_unit_triage <- function(con, id, path, scaffold, verdict,
                            reason = NA_character_, source = c("user", "auto"),
                            rule = NA_character_) {
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
    "UPDATE assemblies SET triage = ?, triage_reason = ?, triage_source = ?,
       triage_rule = ?, ignore = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(verdict, reason %||% NA_character_, source,
      rule %||% NA_character_, ig, id, path, scaffold))
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
       a.triage, a.triage_reason, a.triage_source, a.triage_rule,
       COALESCE(c.median_depth, c.mean_depth) AS depth,
       b.species AS screen_species, b.lineage_taxids AS screen_taxids,
       b.rank_distance, b.is_nuclear_rrna,
       s.lineage_taxids AS sample_taxids, s.sci_name,
       (SELECT COUNT(*) FROM annotations n
          WHERE n.ID = a.ID AND n.path = a.path AND n.scaffold = a.scaffold) AS n_features,
       (SELECT an.n_features FROM annotate an
          WHERE an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold) AS annotated_n,
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

#' Apply automatic triage to one sample
#'
#' Skips any unit a human has already decided (`triage_source = 'user'`) and any
#' unit whose annotation is locked, so re-running is safe.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param opts Named list of per-rule actions; NULL reads the sample's
#'   triage_opts preset.
#' @param overwrite_user Also re-decide units a curator set by hand.
#' @return (invisibly) the number of units changed.
#' @export
auto_triage_sample <- function(con, id, opts = NULL, overwrite_user = FALSE) {
  d <- unit_triage_state(con, id)
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
  d <- triage_rules(d, act)
  changed <- 0L
  for (i in seq_len(nrow(d))) {
    if (!isTRUE(overwrite_user) && identical(d$triage_source[i], "user")) next
    if (isTRUE(d$annotate_lock[i] == 1)) next
    same <- identical(d$triage[i] %|NA|% "", d$new_triage[i] %|NA|% "") &&
      identical(d$triage_reason[i] %|NA|% "", d$new_reason[i] %|NA|% "") &&
      identical(d$triage_rule[i] %|NA|% "", d$new_rule[i] %|NA|% "")
    if (same) next
    set_unit_triage(con, id, d$path[i], d$scaffold[i], d$new_triage[i],
      reason = d$new_reason[i], source = "auto", rule = d$new_rule[i])
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
#' @param opts Named list of per-rule actions; NULL reads the sample's
#'   triage_opts preset.
#' @param overwrite_user Include units a curator decided by hand.
#' @return list(n_units, n_drop, n_undecided, n_keep, n_samples, emptied,
#'   n_user_skipped, n_locked_skipped, by_rule).
#' @export
triage_preview <- function(con, ids, opts = NULL, overwrite_user = FALSE) {
  res <- list(n_units = 0L, n_drop = 0L, n_undecided = 0L, n_keep = 0L,
    n_samples = 0L, emptied = character(0), n_user_skipped = 0L,
    n_locked_skipped = 0L, by_rule = NULL)
  acc <- NULL
  for (id in unique(ids)) {
    d <- unit_triage_state(con, id)
    if (nrow(d) == 0) next
    act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
    cnt <- triage_rule_counts(d, act)
    acc <- if (is.null(acc)) cnt else {
      acc$n_decided <- acc$n_decided + cnt$n_decided
      acc$n_matching <- acc$n_matching + cnt$n_matching
      acc
    }
    d <- triage_rules(d, act)
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
  res$by_rule <- acc %||% triage_rule_counts(
    data.frame()[0, ], normalize_triage_actions(opts))
  res
}

#' Vectorised "is this row owned by a user"
#' @noRd
identical_user <- function(x) !is.na(x) & x == "user"
