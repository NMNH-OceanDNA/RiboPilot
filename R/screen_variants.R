# Intragenomic rDNA variants: assembly-graph isomers versus genuine paralogs.
#
# A sample routinely yields several assembly units and the curator is asked to
# adjudicate between them with no information about why they differ. Two very
# different situations look identical: units that are the same molecule
# assembled two ways, differing only in how many IGS subrepeats the graph
# traversed, and units that are genuinely different rDNA copies from one
# individual. The first is a decision that does not matter; the second
# determines what GenBank receives.
#
# Only units the screen calls on_target are compared, so any difference between
# them is intragenomic by construction and no threshold has to separate a
# paralog from a contaminant (doc 12 section 1). The finding is ORTHOGONAL to
# the verdict: a paralogous variant is emphatically on_target.
#
# Only the ANCHORS are compared. The spacers are what varies in both cases and
# so cannot discriminate between them. RiboPilot's feature model has no IGS
# element, so the isomer test is "anchors identical, unannotated remainder
# differs in length" rather than a measurement of subrepeat copy number.

#' The anchors, in transcription order. Never coordinate order: a circularly
#' permuted assembly would otherwise yield a different string for one molecule.
#' @noRd
SCREEN_ANCHOR_ELEMENTS <- c("18S", "5.8S", "28S")

#' Separator introducing the variant clause inside `screen_reason`.
#'
#' The clause is appended to the verdict's own reason rather than given a column
#' of its own, so everything from this token onward is stripped before a rewrite:
#' a rerun must not stack duplicates.
#' @noRd
VARIANT_REASON_SEP <- " | variant: "

#' Percent identity between two anchor strings
#'
#' Identical strings short-circuit to 100 rather than aligning. That is the
#' common isomer case, and an overlap alignment over roughly 4 kb of anchor is
#' O(n*m). A failed alignment returns NA rather than erroring: one unalignable
#' pair must not sink a sample's whole screen.
#' @return numeric percent on the 0-100 scale [pwalign_pid()] uses, or NA.
#' @noRd
screen_anchor_pid <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  if (is.na(a[1]) || is.na(b[1]) || !nzchar(a[1]) || !nzchar(b[1])) return(NA_real_)
  if (identical(a[1], b[1])) return(100)
  aln <- tryCatch(pairwise_nt(a[1], b[1], type = "overlap"), error = function(e) NULL)
  if (is.null(aln)) return(NA_real_)
  tryCatch(as.numeric(pwalign_pid(aln)), error = function(e) NA_real_)
}

#' Group a sample's comparable units by anchor identity
#'
#' Grouping is transitive-by-threshold: units mutually at or above `min_pid`
#' form one group, and the chain A-B-C where A does not match C still forms one
#' group. `variant_identity` then makes the spread visible rather than hiding
#' it. Chained groups are expected to be rare at the identities involved and are
#' not worth a clustering algorithm.
#'
#' `variant_class` is a per-SAMPLE property once two groups exist: telling the
#' curator of one isomer "pick either" would be false when the sample also holds
#' a distinct copy. `variant_group` is what says which units are interchangeable.
#'
#' @param units data.frame(key, anchor, total_length). An NA `anchor` is not
#'   comparable and comes back all-NA, which is distinct from `singleton`.
#' @param min_pid The `min_isomer_anchor_pid` in force for this sample.
#' @return data.frame(key, variant_class, variant_group, variant_identity,
#'   len_delta), one row per input unit, in input order.
#' @noRd
screen_variant_groups <- function(units, min_pid) {
  n <- if (is.null(units)) 0L else nrow(units)
  out <- data.frame(key = if (n) units$key else character(0),
    variant_class = rep(NA_character_, n), variant_group = rep(NA_character_, n),
    variant_identity = rep(NA_real_, n), len_delta = rep(NA_integer_, n),
    stringsAsFactors = FALSE)
  if (n == 0) return(out)

  # %|NA|% is scalar only and errors on a vector, so coalesce explicitly
  anc <- units$anchor
  anc[is.na(anc)] <- ""
  cmp <- which(nzchar(anc))
  if (!length(cmp)) return(out)
  if (length(cmp) == 1) {
    # no second sequence exists, so there is no identity to record
    out$variant_class[cmp] <- "singleton"
    out$variant_group[cmp] <- "g1"
    return(out)
  }

  k <- length(cmp)
  pid <- matrix(NA_real_, k, k)
  for (i in seq_len(k - 1)) {
    for (j in seq(i + 1, k)) {
      pid[i, j] <- pid[j, i] <- screen_anchor_pid(units$anchor[cmp[i]],
        units$anchor[cmp[j]])
    }
  }

  # connected components over the at-or-above-threshold edges
  comp <- rep(NA_integer_, k)
  nxt <- 0L
  for (i in seq_len(k)) {
    if (!is.na(comp[i])) next
    nxt <- nxt + 1L
    queue <- i
    while (length(queue)) {
      v <- queue[1]; queue <- queue[-1]
      if (!is.na(comp[v])) next
      comp[v] <- nxt
      linked <- which(!is.na(pid[v, ]) & pid[v, ] >= min_pid & is.na(comp))
      queue <- c(queue, linked)
    }
  }

  len <- as.integer(units$total_length[cmp])
  rank_len <- len
  rank_len[is.na(rank_len)] <- -1L
  # representative = longest member, ties broken by key so the label is stable
  reps <- vapply(seq_len(nxt), function(g) {
    m <- which(comp == g)
    m[order(-rank_len[m], units$key[cmp[m]])][1]
  }, integer(1))
  # g1 is the group whose representative is longest, again tie-broken by key
  gord <- order(-rank_len[reps], units$key[cmp[reps]])
  label <- character(nxt)
  label[gord] <- paste0("g", seq_len(nxt))
  cls <- if (nxt == 1L) "isomer" else "paralog"

  for (i in seq_len(k)) {
    r <- reps[comp[i]]
    out$variant_class[cmp[i]] <- cls
    out$variant_group[cmp[i]] <- label[comp[i]]
    if (i != r) {
      out$variant_identity[cmp[i]] <- pid[i, r]
      out$len_delta[cmp[i]] <- abs(len[i] - len[r])
    }
  }
  out
}

#' One unit's annotated anchors concatenated into a single comparable string
#'
#' A unit missing all three anchors is not comparable and returns NA, which is
#' distinct from being a `singleton`: NA says nothing was measured.
#'
#' @param seq The unit's `assemblies.sequence`.
#' @param ann Its `annotations` rows, with `element`, `start`, `stop`, `strand`.
#' @return character(1), or `NA_character_` when the unit is not comparable.
#' @noRd
screen_anchor_string <- function(seq, ann) {
  if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) return(NA_character_)
  if (is.null(ann) || nrow(ann) == 0 || !"element" %in% names(ann)) return(NA_character_)
  a <- ann[ann$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (nrow(a) == 0) return(NA_character_)
  ok <- !is.na(a$start) & !is.na(a$stop) & a$stop >= a$start
  a <- a[ok, , drop = FALSE]
  if (nrow(a) == 0) return(NA_character_)

  strand <- if ("strand" %in% names(a)) a$strand else rep(NA_character_, nrow(a))
  strand <- strand[!is.na(strand) & nzchar(strand)]
  # concatenating a sense and an antisense span produces a string that means
  # nothing, and a mixed-strand unit is a chimera signal rather than a paralog
  if (length(strand) && length(unique(strand)) > 1) return(NA_character_)
  minus <- length(strand) > 0 && identical(unique(strand), "-")

  ord <- order(match(a$element, SCREEN_ANCHOR_ELEMENTS))
  a <- a[ord, , drop = FALSE]
  parts <- substr(rep(seq[1], nrow(a)), a$start, a$stop)
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(NA_character_)
  # each span reverse-complemented is that anchor in sense; the element sort
  # above already puts them in transcription order, so nothing is reversed twice
  if (minus) {
    parts <- vapply(parts, function(p)
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(p))),
      character(1), USE.NAMES = FALSE)
  }
  paste(parts, collapse = "")
}

#' The curator-facing sentence for one unit's variant finding
#' @noRd
screen_variant_clause <- function(cls, group, identity, len_delta) {
  if (is.na(cls)) return(NA_character_)
  if (identical(cls, "singleton")) {
    return("singleton, the sample's only comparable on-target unit")
  }
  lead <- sprintf("%s, group %s", cls, group)
  if (is.na(identity)) {
    return(paste0(lead, ", this is the group representative"))
  }
  # the measured identity always prints, so a curator sees the number rather
  # than only the label
  paste0(lead, sprintf(
    ", anchors %.2f percent identical to the group representative", identity),
    if (is.na(len_delta)) "" else sprintf(", total length differs by %d bp", len_delta))
}

#' Strip any previous variant clause from a stored reason
#' @noRd
screen_reason_base <- function(reason) {
  if (is.na(reason)) return(NA_character_)
  # fixed match, never a regex: the separator contains '|', which as a pattern is
  # alternation and matches a bare space anywhere in the verdict's own reason
  i <- regexpr(VARIANT_REASON_SEP, reason, fixed = TRUE)[1]
  if (i > 0) substr(reason, 1, i - 1) else reason
}

#' Classify one sample's on-target units as isomers, paralogs or a singleton
#'
#' Runs AFTER the per-unit loop in [screen_sample()], never before: [screen_one()]
#' rewrites a row with DELETE + append, so a variant column written first would
#' be silently wiped. Writes with UPDATE for the same reason, so the pass cannot
#' destroy the verdict columns screen_one() owns.
#'
#' Only `on_target` units are compared; everything else is a contamination
#' question the screen already owns, and a unit that has changed verdict has any
#' stale finding cleared.
#'
#' No version guard: the only caller is the per-sample `screen` pipeline process,
#' once per run, and identical anchors short-circuit without aligning.
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return (invisibly) the number of units given a variant class.
#' @noRd
screen_variants <- function(con, id) {
  u <- DBI::dbGetQuery(con,
    "SELECT s.path, s.scaffold, s.screen_verdict, s.screen_reason, s.variant_class,
            a.sequence, a.length
       FROM screen s
       LEFT JOIN assemblies a
         ON a.ID = s.ID AND a.path = s.path AND a.scaffold = s.scaffold
      WHERE s.ID = ? ORDER BY s.path, s.scaffold", params = list(id))
  if (nrow(u) == 0) return(invisible(0L))

  keep <- !is.na(u$screen_verdict) & u$screen_verdict == "on_target"
  # a unit that is no longer on_target must not keep a stale finding
  for (i in which(!keep)) {
    DBI::dbExecute(con,
      "UPDATE screen SET variant_class = NULL, variant_group = NULL,
         variant_identity = NULL, screen_reason = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(screen_reason_base(u$screen_reason[i]), id, u$path[i], u$scaffold[i]))
  }
  ui <- which(keep)
  if (!length(ui)) return(invisible(0L))

  anchors <- vapply(ui, function(i) {
    ann <- DBI::dbGetQuery(con,
      "SELECT element, start, stop, strand FROM annotations
        WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, u$path[i], u$scaffold[i]))
    screen_anchor_string(u$sequence[i], ann)
  }, character(1), USE.NAMES = FALSE)

  len <- as.integer(u$length[ui])
  # fall back to the stored sequence when assemblies.length is unset, so a
  # missing column cannot make every unit tie and scramble the representative
  nolen <- is.na(len)
  if (any(nolen)) len[nolen] <- nchar(u$sequence[ui][nolen])

  min_pid <- screen_params(con, id)$min_isomer_anchor_pid %||%
    SCREEN_PARAMS_DEFAULT$min_isomer_anchor_pid
  # as.integer before %d: sprintf errors on a non-whole double, and nothing
  # guarantees the driver hands back integers for these columns
  g <- screen_variant_groups(data.frame(
    key = sprintf("%d.%d", as.integer(u$path[ui]), as.integer(u$scaffold[ui])),
    anchor = anchors, total_length = len, stringsAsFactors = FALSE), min_pid)

  n <- 0L
  for (j in seq_along(ui)) {
    i <- ui[j]
    clause <- screen_variant_clause(g$variant_class[j], g$variant_group[j],
      g$variant_identity[j], g$len_delta[j])
    base <- screen_reason_base(u$screen_reason[i])
    reason <- if (is.na(clause)) base else {
      paste0(base %|NA|% "", VARIANT_REASON_SEP, clause)
    }
    DBI::dbExecute(con,
      "UPDATE screen SET variant_class = ?, variant_group = ?, variant_identity = ?,
         screen_reason = ?
       WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(g$variant_class[j], g$variant_group[j],
        round(g$variant_identity[j], 2), reason, id, u$path[i], u$scaffold[i]))
    if (!is.na(g$variant_class[j])) n <- n + 1L
  }
  invisible(n)
}
