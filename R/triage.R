# Per-scaffold triage: keep / drop / undecided, with a reason and an owner.
#
# `assemblies.ignore` remains the single downstream gate -- the annotate unit query, the
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
       b.species AS screen_species, b.lineage AS screen_lineage,
       b.lineage_taxids AS screen_taxids,
       b.rank_distance, b.is_nuclear_rrna,
       sc.screen_verdict, sc.screen_class, sc.screen_confidence, sc.lca_taxon,
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
     LEFT JOIN screen sc      ON sc.ID = a.ID AND sc.path = a.path AND sc.scaffold = a.scaffold
     LEFT JOIN samples s      ON s.ID = a.ID
     WHERE a.ID = ? ORDER BY a.path, a.scaffold", params = list(id))
}

#' Apply automatic triage to one sample
#'
#' Skips any unit a human has already decided (`triage_source = 'user'`) and any
#' unit whose annotation is locked, so re-running is safe. If every rule is set
#' to "off", this is a complete no-op and writes nothing: see
#' dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3 for why that special
#' case exists.
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param opts Named list of per-rule actions; NULL reads the sample's
#'   triage_opts preset.
#' @param overwrite_user Also re-decide units a curator set by hand.
#' @return (invisibly) the number of units changed.
#' @export
auto_triage_sample <- function(con, id, opts = NULL, overwrite_user = FALSE) {
  # triage is the last per-sample step of WF1, so arriving here is what closes
  # the sample's state -- not what triage goes on to decide, which is why this
  # runs before the early returns below.
  refresh_sample_state(con, id)
  d <- unit_triage_state(con, id)
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
  # every rule off is a deliberate no-op, not "nothing matches, so keep": see
  # dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3
  if (triage_all_off(act)) {
    return(invisible(0L))
  }
  d <- triage_rules(d, act)
  changed <- 0L
  for (i in seq_len(nrow(d))) {
    if (!isTRUE(overwrite_user) && identical(d$triage_source[i], "user")) next
    if (isTRUE(d$annotate_lock[i] == 1)) next
    # A rule that can only flag must not CLEAR a drop. The flag rules assert
    # that a channel is missing, and dropping a scaffold is one way to make a
    # channel missing: parking it means annotation never runs for it, so
    # "annotation has not run" fires and its undecided verdict would put the
    # scaffold straight back. Only a human clears a drop.
    if (identical(d$triage[i], "drop") &&
        (length(d$new_triage[i]) == 0 || is.na(d$new_triage[i]))) next
    if (isTRUE(triage_unit_same(d$triage[i], d$triage_reason[i], d$triage_rule[i],
      d$new_triage[i], d$new_reason[i], d$new_rule[i]))) next
    set_unit_triage(con, id, d$path[i], d$scaffold[i], d$new_triage[i],
      reason = d$new_reason[i], source = "auto", rule = d$new_rule[i])
    changed <- changed + 1L
  }
  invisible(changed)
}

#' Does a proposed verdict match the one already stored
#'
#' Vectorised, with NA and "" the same absent value. The preview's headline is
#' only true while it uses the identical comparison the writer uses, so both
#' call this and neither reimplements it.
#' @noRd
triage_unit_same <- function(stored_triage, stored_reason, stored_rule,
                             new_triage, new_reason, new_rule) {
  navl <- function(x) ifelse(is.na(x), "", as.character(x))
  (navl(stored_triage) == navl(new_triage)) &
    (navl(stored_reason) == navl(new_reason)) &
    (navl(stored_rule) == navl(new_rule))
}

#' What would bulk triage do, without doing it
#'
#' A bulk write across a whole project is the highest-blast-radius action in the
#' app, so it is never applied without showing the count first -- including how
#' many samples would end up with NOTHING kept, which is the outcome that
#' silently loses a sample.
#'
#' If every rule is "off" for a sample, that sample is treated as fully skipped:
#' the preview reports its CURRENT stored state unchanged, with zero counted in
#' every n_new_* counter and in n_changed. See
#' dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3.
#'
#' @param con Open project DB connection.
#' @param ids Sample IDs to consider.
#' @param opts Named list of per-rule actions; NULL reads the sample's
#'   triage_opts preset.
#' @param overwrite_user Include units a curator decided by hand.
#' @return list with, for the END STATE after Apply: n_units, n_drop, n_undecided,
#'   n_keep, n_samples, emptied, n_user_skipped, n_locked_skipped, by_rule; and,
#'   for what Apply would actually CHANGE: n_changed (units that are not skipped
#'   and whose stored triage/reason/rule differ from the proposed ones, by the
#'   same comparison auto_triage_sample() uses to decide what to write), broken
#'   down by proposed verdict into n_new_drop, n_new_keep and n_new_flag, which
#'   therefore sum to n_changed.
#' @export
triage_preview <- function(con, ids, opts = NULL, overwrite_user = FALSE) {
  res <- list(n_units = 0L, n_drop = 0L, n_undecided = 0L, n_keep = 0L,
    n_samples = 0L, emptied = character(0), n_user_skipped = 0L,
    n_locked_skipped = 0L, by_rule = NULL,
    n_new_drop = 0L, n_new_keep = 0L, n_new_flag = 0L, n_changed = 0L)
  acc <- NULL
  for (id in unique(ids)) {
    d <- unit_triage_state(con, id)
    if (nrow(d) == 0) next
    act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
    d <- triage_rules(d, act)
    skip_user <- !isTRUE(overwrite_user) & identical_user(d$triage_source)
    skip_lock <- !is.na(d$annotate_lock) & d$annotate_lock == 1
    # every rule off skips everything: see triage_all_off() and
    # dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3. Deliberately not
    # folded into skip_user/skip_lock below, which must keep counting only
    # hand-decisions and locks
    all_off <- triage_all_off(act)
    skipped <- skip_user | skip_lock | all_off
    # count only the units Apply would touch, so the per-rule table cannot
    # disagree with the totals rendered beside it; the rules still see the whole
    # sample, as Apply does, because sibling concordance depends on it.
    # all_off is deliberately NOT in this mask: Matching is the "what am I
    # giving up" number, and all-off is exactly when it is wanted. Decided is 0
    # under all-off anyway, since no rule can fire.
    cnt <- triage_rule_counts(d, act, mask = !(skip_user | skip_lock))
    if (is.null(acc)) {
      acc <- cnt
    } else {
      acc$n_decided <- acc$n_decided + cnt$n_decided
      acc$n_matching <- acc$n_matching + cnt$n_matching
      # samples may follow different presets; never report one sample's action
      # as if it governed the whole selection
      acc$action[acc$action != cnt$action] <- "mixed"
    }
    res$n_user_skipped <- res$n_user_skipped + sum(skip_user & !skip_lock)
    res$n_locked_skipped <- res$n_locked_skipped + sum(skip_lock)
    eff <- ifelse(skipped, d$triage, d$new_triage)
    res$n_units <- res$n_units + nrow(d)
    res$n_drop <- res$n_drop + sum(eff %in% "drop")
    res$n_keep <- res$n_keep + sum(eff %in% "keep")
    res$n_undecided <- res$n_undecided + sum(is.na(eff))
    res$n_samples <- res$n_samples + 1L
    if (all(eff %in% "drop")) res$emptied <- c(res$emptied, id)
    # what Apply actually writes: units NOT skipped whose proposed triple
    # differs from what is stored. "Newly dropped/kept/flagged" is a breakdown
    # of that set, not of every proposal, or the breakdown outruns the headline
    same <- triage_unit_same(d$triage, d$triage_reason, d$triage_rule,
      d$new_triage, d$new_reason, d$new_rule)
    chg <- !skipped & !same
    res$n_new_drop <- res$n_new_drop + sum(chg & d$new_triage %in% "drop")
    res$n_new_keep <- res$n_new_keep + sum(chg & d$new_triage %in% "keep")
    res$n_new_flag <- res$n_new_flag + sum(chg & is.na(d$new_triage))
    res$n_changed <- res$n_changed + sum(chg)
  }
  res$by_rule <- acc %||% triage_rule_counts(
    data.frame(), normalize_triage_actions(opts %||% triage_opts_row(con)))
  res
}

#' Vectorised "is this row owned by a user"
#' @noRd
identical_user <- function(x) !is.na(x) & x == "user"
