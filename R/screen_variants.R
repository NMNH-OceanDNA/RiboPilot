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

#' Percent identity between two anchor strings, with the coverage each side got
#'
#' Identical strings short-circuit to 100 rather than aligning. That is the
#' common isomer case, and an overlap alignment over roughly 6 kb of anchor is
#' O(n*m). A failed alignment returns NA rather than erroring: one unalignable
#' pair must not sink a sample's whole screen.
#'
#' `cov_a`/`cov_b` are the fraction of EACH input the alignment actually covered.
#' Identity alone is not enough: `pairwise_nt(type = "overlap")` allows free end
#' gaps and pwalign::pid() divides matches by the aligned length, so a fragment
#' matching a substring of a longer string scores 100 while covering half of it.
#' @return list(pid, cov_a, cov_b), any of which may be NA.
#' @noRd
screen_anchor_pid <- function(a, b) {
  out <- list(pid = NA_real_, cov_a = NA_real_, cov_b = NA_real_)
  if (length(a) == 0 || length(b) == 0) return(out)
  if (is.na(a[1]) || is.na(b[1]) || !nzchar(a[1]) || !nzchar(b[1])) return(out)
  if (identical(a[1], b[1])) return(list(pid = 100, cov_a = 1, cov_b = 1))
  aln <- tryCatch(pairwise_nt(a[1], b[1], type = "overlap"), error = function(e) NULL)
  if (is.null(aln)) return(out)
  tryCatch({
    out$pid <- as.numeric(pwalign::pid(aln))
    out$cov_a <- alignment_coverage(as.character(pwalign::pattern(aln)), nchar(a[1]))
    out$cov_b <- alignment_coverage(as.character(pwalign::subject(aln)), nchar(b[1]))
  }, error = function(e) NULL)
  out
}

#' Compare two units on the anchors they BOTH annotate
#'
#' Comparing one unit's 18S+5.8S+28S against another's 28S alone is not a
#' comparison, however high the identity comes out. Only shared elements are
#' concatenated, in transcription order, and the pair is comparable only when the
#' shared span is long enough to mean something and the alignment covers enough of
#' BOTH sides.
#' @param pa,pb Named anchor vectors from [screen_anchor_parts()].
#' @return list(pid, comparable, shared, n_bp, coverage)
#' @noRd
screen_anchor_pair <- function(pa, pb, params = SCREEN_PARAMS_DEFAULT) {
  params <- params %||% SCREEN_PARAMS_DEFAULT
  out <- list(pid = NA_real_, comparable = FALSE, shared = character(0),
    n_bp = NA_integer_, coverage = NA_real_)
  shared <- intersect(names(pa), names(pb))
  shared <- SCREEN_ANCHOR_ELEMENTS[SCREEN_ANCHOR_ELEMENTS %in% shared]
  if (!length(shared)) return(out)
  a <- paste(pa[shared], collapse = "")
  b <- paste(pb[shared], collapse = "")
  out$shared <- shared
  out$n_bp <- as.integer(min(nchar(a), nchar(b)))
  r <- screen_anchor_pid(a, b)
  out$pid <- r$pid
  out$coverage <- suppressWarnings(min(r$cov_a, r$cov_b))
  if (!is.finite(out$coverage)) out$coverage <- NA_real_
  out$comparable <- !is.na(out$pid) && !is.na(out$coverage) &&
    out$n_bp >= (params$min_isomer_anchor_bp %||% 0L) &&
    out$coverage >= (params$min_isomer_anchor_coverage %||% 0)
  out
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
#' A unit with anchors but no COMPARABLE counterpart is `NA`, not `paralog`.
#' Landing it in a group of its own would make the sample read as carrying
#' several rDNA copies on the strength of a measurement that was never made.
#'
#' @param units data.frame(key, total_length), one row per unit.
#' @param parts List of named anchor vectors from [screen_anchor_parts()], same
#'   order and length as `units`. An empty entry is not comparable.
#' @param params The `screen` sublist in force for this sample.
#' @return data.frame(key, variant_class, variant_group, variant_identity,
#'   len_delta, is_rep), one row per input unit, in input order.
#' @noRd
screen_variant_groups <- function(units, parts, params = SCREEN_PARAMS_DEFAULT) {
  params <- params %||% SCREEN_PARAMS_DEFAULT
  min_pid <- params$min_isomer_anchor_pid %||% 99.5
  n <- if (is.null(units)) 0L else nrow(units)
  out <- data.frame(key = if (n) units$key else character(0),
    variant_class = rep(NA_character_, n), variant_group = rep(NA_character_, n),
    variant_identity = rep(NA_real_, n), len_delta = rep(NA_integer_, n),
    is_rep = rep(FALSE, n), stringsAsFactors = FALSE)
  if (n == 0) return(out)

  cand <- which(vapply(parts, function(p) length(p) > 0, logical(1)))
  if (!length(cand)) return(out)
  singleton <- function(i) {
    # no second sequence to compare against, so no identity to record
    out$variant_class[i] <- "singleton"; out$variant_group[i] <- "g1"
    out$is_rep[i] <- TRUE
    out
  }
  if (length(cand) == 1) return(singleton(cand))

  k <- length(cand)
  pid <- matrix(NA_real_, k, k)
  cmpbl <- matrix(FALSE, k, k)
  for (i in seq_len(k - 1)) {
    for (j in seq(i + 1, k)) {
      r <- screen_anchor_pair(parts[[cand[i]]], parts[[cand[j]]], params)
      pid[i, j] <- pid[j, i] <- r$pid
      cmpbl[i, j] <- cmpbl[j, i] <- isTRUE(r$comparable)
    }
  }

  # drop units nothing could be compared against BEFORE counting groups
  sel <- which(apply(cmpbl, 1, any))
  if (!length(sel)) return(out)
  if (length(sel) == 1) return(singleton(cand[sel]))
  idx <- cand[sel]
  pid <- pid[sel, sel, drop = FALSE]
  cmpbl <- cmpbl[sel, sel, drop = FALSE]
  k <- length(sel)

  # connected components over comparable, at-or-above-threshold edges
  edge <- cmpbl & !is.na(pid) & pid >= min_pid
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
      queue <- c(queue, which(edge[v, ] & is.na(comp)))
    }
  }

  len <- as.integer(units$total_length[idx])
  rank_len <- len
  rank_len[is.na(rank_len)] <- -1L
  # representative = longest member, ties broken by key so the label is stable
  reps <- vapply(seq_len(nxt), function(g) {
    m <- which(comp == g)
    m[order(-rank_len[m], units$key[idx[m]])][1]
  }, integer(1))
  # g1 is the group whose representative is longest, again tie-broken by key
  label <- character(nxt)
  label[order(-rank_len[reps], units$key[idx[reps]])] <- paste0("g", seq_len(nxt))
  cls <- if (nxt == 1L) "isomer" else "paralog"

  for (i in seq_len(k)) {
    r <- reps[comp[i]]
    out$variant_class[idx[i]] <- cls
    out$variant_group[idx[i]] <- label[comp[i]]
    if (i == r) {
      out$is_rep[idx[i]] <- TRUE
    } else {
      # NA when this member reached its group through a chain rather than a
      # direct comparison with the representative; that spread is real
      out$variant_identity[idx[i]] <- pid[i, r]
      out$len_delta[idx[i]] <- abs(len[i] - len[r])
    }
  }
  out
}

#' One unit's annotated anchors, kept SEPARATE and named by element
#'
#' Separate rather than pre-concatenated so a pair can be compared on the
#' elements both units actually annotate. A unit with no usable anchor returns
#' `character(0)`, which is not comparable, and is distinct from being a
#' `singleton`: it says nothing was measured.
#'
#' @param seq The unit's `assemblies.sequence`.
#' @param ann Its `annotations` rows, with `element`, `start`, `stop`, `strand`.
#' @return named character vector in transcription order, possibly empty.
#' @noRd
screen_anchor_parts <- function(seq, ann) {
  none <- character(0)
  if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) return(none)
  if (is.null(ann) || nrow(ann) == 0 || !"element" %in% names(ann)) return(none)
  a <- ann[ann$element %in% SCREEN_ANCHOR_ELEMENTS, , drop = FALSE]
  if (nrow(a) == 0) return(none)
  ok <- !is.na(a$start) & !is.na(a$stop) & a$stop >= a$start
  a <- a[ok, , drop = FALSE]
  if (nrow(a) == 0) return(none)

  strand <- if ("strand" %in% names(a)) a$strand else rep(NA_character_, nrow(a))
  strand <- strand[!is.na(strand) & nzchar(strand)]
  # a sense and an antisense span cannot be compared with each other, and a
  # mixed-strand unit is a chimera signal rather than a paralog
  if (length(strand) && length(unique(strand)) > 1) return(none)
  minus <- length(strand) > 0 && identical(unique(strand), "-")

  # transcription order, never coordinate order, and one row per element
  a <- a[order(match(a$element, SCREEN_ANCHOR_ELEMENTS)), , drop = FALSE]
  a <- a[!duplicated(a$element), , drop = FALSE]
  parts <- substr(rep(seq[1], nrow(a)), a$start, a$stop)
  names(parts) <- a$element
  parts <- parts[nzchar(parts)]
  if (!length(parts)) return(none)
  # each span reverse-complemented is that anchor in sense; the element sort
  # above already puts them in transcription order, so nothing is reversed twice
  if (minus) {
    nm <- names(parts)
    parts <- vapply(parts, function(p)
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(p))),
      character(1), USE.NAMES = FALSE)
    names(parts) <- nm
  }
  parts
}

#' All of one unit's anchors as a single string (display and length only)
#'
#' NOT used for comparison: see [screen_anchor_pair()], which compares only the
#' elements two units share.
#' @noRd
screen_anchor_string <- function(seq, ann) {
  p <- screen_anchor_parts(seq, ann)
  if (!length(p)) return(NA_character_)
  paste(p, collapse = "")
}

#' The curator-facing sentence for one unit's variant finding
#' @noRd
screen_variant_clause <- function(cls, group, identity, len_delta, is_rep = FALSE) {
  if (is.na(cls)) return(NA_character_)
  if (identical(cls, "singleton")) {
    return("singleton, the sample's only comparable on-target unit")
  }
  lead <- sprintf("%s, group %s", cls, group)
  if (isTRUE(is_rep)) return(paste0(lead, ", this is the group representative"))
  # a non-representative with no identity reached its group through a chain, not
  # a direct comparison; say so rather than implying a measurement
  if (is.na(identity)) {
    return(paste0(lead, ", grouped through an intermediate member, not compared",
      " with the representative directly"))
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

  parts <- lapply(ui, function(i) {
    ann <- DBI::dbGetQuery(con,
      "SELECT element, start, stop, strand FROM annotations
        WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, u$path[i], u$scaffold[i]))
    screen_anchor_parts(u$sequence[i], ann)
  })

  len <- as.integer(u$length[ui])
  # fall back to the stored sequence when assemblies.length is unset, so a
  # missing column cannot make every unit tie and scramble the representative
  nolen <- is.na(len)
  if (any(nolen)) len[nolen] <- nchar(u$sequence[ui][nolen])

  # as.integer before %d: sprintf errors on a non-whole double, and nothing
  # guarantees the driver hands back integers for these columns
  g <- screen_variant_groups(data.frame(
    key = sprintf("%d.%d", as.integer(u$path[ui]), as.integer(u$scaffold[ui])),
    total_length = len, stringsAsFactors = FALSE), parts, screen_params(con, id))

  n <- 0L
  for (j in seq_along(ui)) {
    i <- ui[j]
    clause <- screen_variant_clause(g$variant_class[j], g$variant_group[j],
      g$variant_identity[j], g$len_delta[j], g$is_rep[j])
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
