# Named, ordered triage rules over the two independent evidence channels:
# annotation (barrnap + ITSx, local) and the GenBank screen (remote). Their
# independence is the safety argument for dropping on their conjunction, so a
# rule may only drop when BOTH channels are present. Rules 1-2 assert that a
# channel is missing and therefore must precede the combination rules.

#' Actions a triage rule may be set to.
#' @noRd
TRIAGE_ACTIONS <- c("off", "flag", "drop")

#' Ordered triage rule definitions
#'
#' Each entry is list(key, label, meaning, actions, default, test), where `test`
#' is a pure function (d, i, ctx) returning a reason string or NA.
#' @noRd
TRIAGE_RULE_DEFS <- list(
  list(
    key = "not_annotated",
    label = "Annotation has not run",
    meaning = "Neither annotation tool has reported on this scaffold yet.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (isTRUE(ctx$annotated[i])) return(NA_character_)
      "annotation has not run for this scaffold"
    }
  ),
  list(
    key = "no_screen",
    label = "GenBank screen has not run",
    meaning = "No usable core_nt hit is recorded, so the scaffold cannot be judged against GenBank.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (isTRUE(ctx$screened[i])) return(NA_character_)
      "no GenBank hit says whether this scaffold is nuclear rDNA"
    }
  ),
  list(
    key = "no_feat_non_rdna",
    label = "No annotations, non-rDNA hit",
    meaning = "Annotation found no rRNA or ITS feature, and the best GenBank match is not a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "drop",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] == 0) || !isTRUE(ctx$nuc[i] == 0)) return(NA_character_)
      sprintf(paste("annotation found no rRNA or ITS feature, and the best GenBank match",
        "(%s) is not a nuclear rDNA record"), ctx$sp[i])
    }
  ),
  list(
    key = "no_feat_rdna_hit",
    label = "No annotations but rDNA hit",
    meaning = "Annotation found no rRNA or ITS feature, though the best GenBank match is a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] == 0) || !isTRUE(ctx$nuc[i] == 1)) return(NA_character_)
      sprintf(paste("annotation found no rRNA or ITS feature, though the best GenBank match",
        "(%s) is a nuclear rDNA record"), ctx$sp[i])
    }
  ),
  list(
    key = "feat_non_rdna_hit",
    label = "Has annotations but non-rDNA hit",
    meaning = "The scaffold has a feature model, but its best GenBank match is not a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] > 0) || !isTRUE(ctx$nuc[i] == 0)) return(NA_character_)
      sprintf("best GenBank match (%s) is not a nuclear rDNA record", ctx$sp[i])
    }
  ),
  list(
    key = "screen_off_target",
    label = "GenBank taxonomy disagrees with the sample",
    meaning = "The scaffold's BLAST hits agree on a taxon that is not on this sample's lineage.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!identical(ctx$screen[i], "off_target")) return(NA_character_)
      sprintf("the scaffold's GenBank hits agree on %s, which does not sit on this sample's lineage",
        ctx$screen_taxon[i] %|NA|% "another taxon")
    }
  ),
  list(
    key = "odd_one_out",
    label = "Disagrees with sibling scaffolds",
    meaning = "This scaffold's top hit disagrees with the majority of this sample's other scaffolds.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$screened[i]) || !isTRUE(ctx$sib$odd_one_out[i])) return(NA_character_)
      sprintf("top GenBank hit (%s) differs from this sample's other scaffolds (%s at %s)",
        ctx$sp[i], ctx$sib$sib_taxon[i] %|NA|% "the majority",
        ctx$sib$sib_rank[i] %|NA|% "an agreed rank")
    }
  ),
  list(
    key = "chimeric",
    label = "Likely chimeric assembly",
    meaning = "Two or more independent signals place a junction between different lineages or array copies inside this scaffold.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (!identical(ctx$chimera[i], "likely")) return(NA_character_)
      ctx$chimera_reason[i] %|NA|% "chimera detector scored this scaffold as likely"
    }
  )
)

#' Rule keys in evaluation order
#' @noRd
triage_rule_keys <- function() {
  vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "key")
}

#' The shipped action for every rule
#' @noRd
triage_actions_default <- function() {
  stats::setNames(lapply(TRIAGE_RULE_DEFS, `[[`, "default"), triage_rule_keys())
}

#' TRUE when every rule's normalized action is "off"
#' @noRd
triage_all_off <- function(actions) {
  all(vapply(triage_rule_keys(), function(k) identical(actions[[k]], "off"), logical(1)))
}

#' Coerce a partial or hand-edited action set into a complete, legal one
#'
#' An action a rule cannot support is softened to "flag" rather than honoured or
#' rejected: a hand-edited database must not be able to arm a drop the evidence
#' cannot support, and must not be able to make the app error either.
#' @noRd
normalize_triage_actions <- function(x) {
  out <- triage_actions_default()
  if (is.null(x)) return(out)
  for (def in TRIAGE_RULE_DEFS) {
    v <- x[[def$key]]
    if (is.null(v) || length(v) == 0 || is.na(v[1])) next
    v <- tolower(trimws(as.character(v[1])))
    if (!v %in% TRIAGE_ACTIONS) next
    if (!v %in% def$actions) v <- "flag"
    out[[def$key]] <- v
  }
  out
}

#' Shared per-sample context every rule test reads
#' @noRd
triage_rule_ctx <- function(d) {
  n <- nrow(d)
  nuc <- d$is_nuclear_rrna %||% rep(NA_integer_, n)
  sp <- d$screen_species %||% rep(NA_character_, n)
  list(
    # annotate.n_features is written by reconcile_one() on every run, including
    # one that found nothing, and only by reconcile, so it survives a re-queue.
    annotated = !is.na(d$annotated_n %||% rep(NA_integer_, n)),
    # a screen row with no classifiable hit is not a screen for rule purposes
    screened = !is.na(nuc),
    n_feat = d$n_features %||% rep(0L, n),
    nuc = nuc,
    sp = ifelse(is.na(sp), "unknown", sp),
    screen = d$screen_verdict %||% rep(NA_character_, n),
    screen_taxon = d$lca_taxon %||% rep(NA_character_, n),
    sib = unit_concordance(data.frame(
      screen_taxids = d$screen_taxids %||% rep(NA_character_, n),
      screen_lineage = d$screen_lineage %||% rep(NA_character_, n),
      stringsAsFactors = FALSE)),
    chimera = d$chimera_label %||% rep(NA_character_, n),
    chimera_reason = d$chimera_reason %||% rep(NA_character_, n)
  )
}

#' Every rule's reason for every unit, independent of order and of settings
#' @return character matrix nrow(d) x length(TRIAGE_RULE_DEFS), colnames = keys.
#' @noRd
triage_rule_matrix <- function(d) {
  keys <- triage_rule_keys()
  n <- nrow(d)
  m <- matrix(NA_character_, nrow = n, ncol = length(keys),
    dimnames = list(NULL, keys))
  if (n == 0) return(m)
  ctx <- triage_rule_ctx(d)
  for (j in seq_along(TRIAGE_RULE_DEFS)) {
    tst <- TRIAGE_RULE_DEFS[[j]]$test
    for (i in seq_len(n)) m[i, j] <- tst(d, i, ctx)
  }
  m
}

#' Propose a triage verdict for every unit of one sample
#'
#' First enabled rule that matches wins. Nothing matching means keep.
#'
#' @param d Frame from [unit_triage_state()].
#' @param opts Named list key -> action; see [normalize_triage_actions()].
#' @param m Rule matrix for `d`, for callers that already built one.
#' @return `d` with `new_triage`, `new_reason` and `new_rule` columns.
#' @noRd
triage_rules <- function(d, opts = triage_actions_default(), m = triage_rule_matrix(d)) {
  if (is.null(d)) d <- data.frame()
  n <- nrow(d)
  d$new_triage <- rep(NA_character_, n)
  d$new_reason <- rep(NA_character_, n)
  d$new_rule <- rep(NA_character_, n)
  if (n == 0) return(d)
  # catches a wrong row count, stripped dimnames and truncated columns. It
  # cannot catch a row-permuted or foreign same-shape m, which would attribute a
  # verdict to the wrong scaffold: the matrix carries no row identity to check
  # against. Safe only because the sole caller builds m from this same d.
  stopifnot(nrow(m) == n, identical(colnames(m), triage_rule_keys()))
  act <- normalize_triage_actions(opts)
  keys <- colnames(m)
  for (i in seq_len(n)) {
    hit <- NA_integer_
    for (j in seq_along(keys)) {
      if (identical(act[[keys[j]]], "off")) next
      if (is.na(m[i, j])) next
      hit <- j
      break
    }
    if (is.na(hit)) {
      d$new_triage[i] <- "keep"
      next
    }
    d$new_rule[i] <- keys[hit]
    d$new_reason[i] <- m[i, hit]
    if (identical(act[[keys[hit]]], "drop")) d$new_triage[i] <- "drop"
  }
  d
}

#' Per-rule counts for the bulk triage modal
#'
#' `n_decided` is first-match and partitions the units. `n_matching` counts
#' every unit meeting the rule, including ones a rule above acted on first, so
#' a rule set to off still reports what turning it back on would reach.
#'
#' `mask` restricts the COUNTING only. The rules stay evaluated over the whole
#' frame because sibling concordance is a property of the full set of scaffolds,
#' so counting a subset by re-running the rules on it would report a different
#' verdict than the one the caller is about to write.
#' @noRd
triage_rule_counts <- function(d, opts = triage_actions_default(), mask = NULL) {
  act <- normalize_triage_actions(opts)
  keys <- triage_rule_keys()
  out <- data.frame(
    key = keys,
    label = vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "label"),
    meaning = vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "meaning"),
    action = vapply(keys, function(k) act[[k]], character(1)),
    n_decided = 0L, n_matching = 0L, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  if (nrow(d) == 0) return(out)
  m <- triage_rule_matrix(d)
  fired <- triage_rules(d, act, m = m)$new_rule
  if (!is.null(mask)) {
    m <- m[mask, , drop = FALSE]
    fired <- fired[mask]
  }
  out$n_matching <- as.integer(colSums(!is.na(m)))
  out$n_decided <- as.integer(vapply(keys, function(k) sum(fired %in% k), integer(1)))
  out
}
