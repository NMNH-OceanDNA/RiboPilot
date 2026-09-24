# Annotation core: predict the nuclear rDNA feature model (18S / ITS1 / 5.8S /
# ITS2 / 28S) by combining barrnap (Rfam SSU/LSU anchors) and ITSx (spacer
# boundaries), then reconciling into one gap-free, 5'->3' feature table.
#
# Split so the tools (barrnap, ITSx) can run as Nextflow processes while R
# does the parsing + reconciliation: parse_barrnap()/parse_itsx() consume tool
# output; reconcile_features() combines them; run_barrnap()/run_itsx()
# run the tools from R.

#' Parse barrnap GFF3 text into rRNA gene features
#'
#' Also captures the fraction of the rRNA model that aligned: (py)barrnap marks a
#' short hit `product=... (partial)` and adds `Note=aligned only NN percent of
#' the ...`. That coverage is what tells a complete 18S from one missing its 3'
#' V9 region (the terminus most metazoan eDNA assays target), so it is kept
#' rather than discarded.
#' @return data.frame(seqid, element, start, stop, strand, model_cov) for 18S/5.8S/28S,
#'   where `model_cov` is the percent of the model aligned (100 for a full hit,
#'   NA when the tool reported partial without a percentage).
#' @noRd
parse_barrnap <- function(gff_lines) {
  gff <- gff_lines[!startsWith(gff_lines, "#")]
  if (length(gff) == 0) {
    return(data.frame())
  }
  f <- Filter(function(x) length(x) >= 9, strsplit(gff, "\t"))
  if (length(f) == 0) {
    return(data.frame())
  }
  attrs <- vapply(f, `[`, character(1), 9)
  nm <- sub(".*Name=([^;]+).*", "\\1", attrs)
  map <- c("18S_rRNA" = "18S", "5_8S_rRNA" = "5.8S", "5.8S_rRNA" = "5.8S", "28S_rRNA" = "28S")
  el <- unname(map[nm])
  keep <- !is.na(el)
  if (!any(keep)) {
    return(data.frame())
  }
  attrs <- attrs[keep]
  partial <- grepl("partial", attrs, ignore.case = TRUE)
  has_pct <- grepl("aligned only [0-9.]+ percent", attrs)
  model_cov <- rep(NA_real_, length(attrs))
  model_cov[has_pct] <- suppressWarnings(as.numeric(
    sub(".*aligned only ([0-9.]+) percent.*", "\\1", attrs[has_pct])))
  model_cov[!has_pct & !partial] <- 100 # full-length hit, no partial marker
  data.frame(
    seqid = vapply(f, `[`, character(1), 1)[keep],
    element = el[keep],
    start = as.integer(vapply(f, `[`, character(1), 4))[keep],
    stop = as.integer(vapply(f, `[`, character(1), 5))[keep],
    strand = vapply(f, `[`, character(1), 7)[keep],
    model_cov = model_cov,
    stringsAsFactors = FALSE
  )
}

#' Parse an ITSx `.positions.txt` into region features
#' @return data.frame(seqid, element, start, stop) for SSU/ITS1/5.8S/ITS2/LSU found
#' @noRd
parse_itsx <- function(pos_lines) {
  if (length(pos_lines) == 0) {
    return(data.frame())
  }
  pos_lines <- pos_lines[nzchar(trimws(pos_lines))]
  labmap <- c("SSU" = "18S", "ITS1" = "ITS1", "5.8S" = "5.8S", "ITS2" = "ITS2", "LSU" = "28S")
  out <- list()
  for (ln in pos_lines) {
    fields <- strsplit(ln, "\t")[[1]]
    if (length(fields) == 0) next
    sid <- fields[1]
    for (fd in fields) {
      m <- regmatches(fd, regexec("^\\s*(SSU|ITS1|5\\.8S|ITS2|LSU):\\s*(\\d+)-(\\d+)", fd))[[1]]
      if (length(m) == 4) {
        out[[length(out) + 1]] <- data.frame(
          seqid = sid, element = unname(labmap[m[2]]),
          start = as.integer(m[3]), stop = as.integer(m[4]),
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

#' Did ITSx flag the sequence as chimeric?
#'
#' ITSx appends a literal `Chimeric!` token to the `.positions.txt` line when the
#' detected ITS/anchor domains come from different origins (ITSx source, POS
#' output). A stitched assembly that concatenated a non-target scaffold, or a
#' host/symbiont mix, can be chimeric; such a sequence must never enter a
#' reference database, so the signal is surfaced rather than dropped.
#' @param pos_lines Lines of an ITSx `.positions.txt`.
#' @return logical scalar.
#' @noRd
itsx_is_chimeric <- function(pos_lines) {
  length(pos_lines) > 0 && any(grepl("Chimeric!", pos_lines, fixed = TRUE))
}

#' Reconcile barrnap + ITSx features into a gap-free 5'->3' model
#'
#' barrnap is authoritative for 18S/28S extent; ITSx for ITS1/5.8S/ITS2. Missing
#' spacers are derived by subtraction between adjacent anchors.
#' @param bar,its data.frames from [parse_barrnap()] / [parse_itsx()].
#' @return data.frame(element, type, start, stop, length, strand, tool)
#' @noRd
reconcile_features <- function(bar, its) {
  pick <- function(df, el) if (!is.null(df) && nrow(df) > 0 && any(df$element == el)) df[df$element == el, ][1, ] else NULL
  if (is.null(bar)) bar <- data.frame()
  if (nrow(bar) > 0 && "strand" %in% names(bar)) {
    st <- table(bar$strand[bar$strand %in% c("+", "-")])
    if (length(st) == 2 && st[[1]] != st[[2]]) bar <- bar[bar$strand %in% names(which.max(st)), , drop = FALSE]
  }
  # 18S and 28S are the copies nearest the ITS region, upstream and downstream
  core <- if (!is.null(its) && nrow(its) > 0) its[its$element %in% c("ITS1", "5.8S", "ITS2"), , drop = FALSE] else data.frame()
  if (nrow(core) == 0 && nrow(bar) > 0 && any(bar$element == "5.8S")) {
    core <- bar[bar$element == "5.8S", , drop = FALSE]
    core <- core[which.max(core$stop - core$start), ]
  }
  hub <- if (nrow(core) > 0) (min(core$start) + max(core$stop)) / 2 else NA_real_
  near <- function(el, side) {
    d <- if (nrow(bar) > 0) bar[bar$element == el, , drop = FALSE] else data.frame()
    if (nrow(d) == 0) return(NULL)
    if (is.na(hub)) return(d[which.max(d$stop - d$start), ])
    if (side == "up") d <- d[d$start < hub, , drop = FALSE]
    if (side == "down") d <- d[d$stop > hub, , drop = FALSE]
    if (nrow(d) == 0) return(NULL)
    d[which.min(switch(side, up = hub - d$start, down = d$start - hub,
      abs((d$start + d$stop) / 2 - hub))), ]
  }
  feats <- list()
  add <- function(element, start, stop, tool, strand = NA_character_, model_cov = NA_real_) {
    if (is.null(start) || is.null(stop) || is.na(start) || is.na(stop) || stop < start) {
      return(invisible())
    }
    feats[[length(feats) + 1]] <<- data.frame(
      element = element, type = rrna_element_type(element),
      start = as.integer(start), stop = as.integer(stop), strand = strand,
      model_cov = as.numeric(model_cov), tool = tool, stringsAsFactors = FALSE)
  }

  b18 <- near("18S", "up"); i18 <- pick(its, "18S")
  b58 <- near("5.8S", "at"); i58 <- pick(its, "5.8S")
  b28 <- near("28S", "down"); i28 <- pick(its, "28S")
  its1 <- pick(its, "ITS1"); its2 <- pick(its, "ITS2")

  s18 <- b18 %||% i18
  s58 <- i58 %||% b58
  s28 <- b28 %||% i28

  # barrnap reports strand; ITSx does not. Carry the barrnap strand where it is
  # the source; leave NA otherwise and fill from the anchor consensus below, so
  # a reverse-oriented assembly is recorded as '-' (and caught by curation's
  # wrong_strand check) instead of being silently asserted '+'.
  bstrand <- function(x) if (!is.null(x) && "strand" %in% names(x)) x$strand else NA_character_
  bcov <- function(x) if (!is.null(x) && "model_cov" %in% names(x)) x$model_cov else NA_real_

  if (!is.null(s18)) add("18S", s18$start, s18$stop,
    if (!is.null(b18)) "barrnap" else "ITSx",
    if (!is.null(b18)) bstrand(b18) else NA_character_,
    if (!is.null(b18)) bcov(b18) else NA_real_)
  if (!is.null(its1)) {
    add("ITS1", its1$start, its1$stop, "ITSx")
  } else if (!is.null(s18) && !is.null(s58)) {
    add("ITS1", s18$stop + 1, s58$start - 1, "derived")
  }
  if (!is.null(s58)) add("5.8S", s58$start, s58$stop,
    if (!is.null(i58)) "ITSx" else "barrnap",
    if (!is.null(i58)) NA_character_ else bstrand(b58),
    if (!is.null(i58)) NA_real_ else bcov(b58))
  if (!is.null(its2)) {
    add("ITS2", its2$start, its2$stop, "ITSx")
  } else if (!is.null(s58) && !is.null(s28)) {
    add("ITS2", s58$stop + 1, s28$start - 1, "derived")
  }
  if (!is.null(s28)) add("28S", s28$start, s28$stop,
    if (!is.null(b28)) "barrnap" else "ITSx",
    if (!is.null(b28)) bstrand(b28) else NA_character_,
    if (!is.null(b28)) bcov(b28) else NA_real_)

  if (length(feats) == 0) {
    return(data.frame())
  }
  out <- do.call(rbind, feats)
  # Consensus strand from the anchors barrnap actually called; features with no
  # strand of their own (ITSx anchors, derived/ITSx spacers) inherit it.
  anch <- out$strand[out$strand %in% c("+", "-")]
  consensus <- if (length(anch) == 0) "+" else names(sort(table(anch), decreasing = TRUE))[1]
  out$strand[!(out$strand %in% c("+", "-"))] <- consensus

  # Clip each spacer to sit strictly between its flanking anchors. barrnap (anchor
  # extent) and ITSx (spacer extent) can disagree at a boundary, leaving e.g. an
  # ITS2 that overlaps the 5.8S; the anchor extent is authoritative, so the spacer
  # is trimmed to the gap. A spacer clipped out of existence is dropped.
  clip_spacer <- function(out, spacer, left, right) {
    si <- which(out$element == spacer)
    if (length(si) != 1) return(out)
    li <- which(out$element == left)
    ri <- which(out$element == right)
    if (length(li) == 1) out$start[si] <- max(out$start[si], out$stop[li] + 1L)
    if (length(ri) == 1) out$stop[si] <- min(out$stop[si], out$start[ri] - 1L)
    out
  }
  out <- clip_spacer(out, "ITS1", "18S", "5.8S")
  out <- clip_spacer(out, "ITS2", "5.8S", "28S")
  out <- out[out$stop >= out$start, , drop = FALSE]

  # A gene end up to 5 bp short of the ITSx spacer beside it is moved to meet
  # it: ITSx places these boundaries, as it does at NCBI.
  for (pr in list(c("18S", "ITS1"), c("ITS1", "5.8S"), c("5.8S", "ITS2"), c("ITS2", "28S"))) {
    li <- which(out$element == pr[1])
    ri <- which(out$element == pr[2])
    if (length(li) != 1 || length(ri) != 1) next
    gap <- out$start[ri] - out$stop[li] - 1L
    if (gap < 1L || gap > 5L) next
    if (pr[1] %in% RRNA_GENES) {
      out$stop[li] <- out$start[ri] - 1L
      out$tool[li] <- paste0(out$tool[li], "+ITSx")
    } else {
      out$start[ri] <- out$stop[li] + 1L
      out$tool[ri] <- paste0(out$tool[ri], "+ITSx")
    }
  }

  out$length <- out$stop - out$start + 1
  out <- out[order(out$start),
    c("element", "type", "start", "stop", "length", "strand", "model_cov", "tool")]
  rownames(out) <- NULL
  out
}

#' Translate a kingdom to the vocabulary of the tool that will receive it
#'
#' The two callers do not share a kingdom vocabulary: pybarrnap takes
#' bac/arc/euk/all, while the bioconda Perl barrnap fork takes bac/arc/fun,
#' where "fun" selects its eukaryote Rfam covariance models and is NOT a
#' fungi-only mode (dev/design/03_ANNOTATION_CURATION.md). This is the single
#' place that mapping lives on the R side; the shell branch of
#' inst/nextflow/modules/pipeline_r.nf mirrors the same rule.
#' @param kingdom Stored kingdom (annotate_opts.rrna_kingdom), e.g. "euk".
#' @param tool "barrnap" or "pybarrnap".
#' @return The kingdom string to pass to `tool`.
#' @noRd
rrna_kingdom_for <- function(kingdom, tool) {
  k <- if (length(kingdom) != 1 || is.na(kingdom) || !nzchar(kingdom)) "euk" else tolower(kingdom)
  if (identical(tool, "barrnap")) {
    if (k %in% c("euk", "all")) "fun" else k
  } else {
    if (identical(k, "fun")) "euk" else k
  }
}

#' Locate the rRNA anchor caller: pybarrnap preferred, Perl barrnap fallback
#' @return list(tool, bin) or NULL when neither is installed.
#' @noRd
rrna_caller <- function() {
  py <- unname(Sys.which("pybarrnap"))
  if (nzchar(py)) {
    return(list(tool = "pybarrnap", bin = py))
  }
  bar <- unname(Sys.which("barrnap"))
  if (nzchar(bar)) {
    return(list(tool = "barrnap", bin = bar))
  }
  NULL
}

#' Run an rRNA caller and return its GFF3 lines
#'
#' stderr is captured rather than discarded: a wrong kingdom makes the tool exit
#' non-zero, and silently dropping that turned a hard failure into "no anchors
#' found". Reported as a warning, not an error, because annotation is
#' best-effort per unit.
#' @noRd
run_rrna_tool <- function(bin, tool, fasta, kingdom, cpu = NULL) {
  args <- if (identical(tool, "barrnap")) {
    c("--kingdom", kingdom, "--quiet", shQuote(fasta))
  } else {
    c("-k", kingdom, if (!is.null(cpu)) c("-t", cpu), "-q", shQuote(fasta))
  }
  err <- tempfile("rp_rrna_err_")
  on.exit(unlink(err), add = TRUE)
  gff <- suppressWarnings(system2(bin, args, stdout = TRUE, stderr = err))
  status <- attr(gff, "status")
  failed <- !is.null(status) && !identical(as.integer(status), 0L)
  if (failed || length(gff[!startsWith(gff, "#")]) == 0) {
    msg <- tryCatch(paste(utils::tail(readLines(err, warn = FALSE), 3L), collapse = " | "),
      error = function(e) "")
    warning(sprintf("%s (kingdom '%s') %s on %s%s", tool, kingdom,
      if (failed) "failed" else "found no rRNA records", basename(fasta),
      if (nzchar(msg)) paste0(": ", msg) else ""), call. = FALSE)
  }
  gff
}

#' Run barrnap (Rfam eukaryote models) on a fasta and parse rRNA features
#' @noRd
run_barrnap <- function(fasta, kingdom = "euk") {
  bin <- unname(Sys.which("barrnap"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  parse_barrnap(run_rrna_tool(bin, "barrnap", fasta, rrna_kingdom_for(kingdom, "barrnap")))
}

#' Largest record ITSx will accept.
#'
#' ITSx runs HMMER's protein pipeline, which refuses a target over 100K and
#' aborts the whole run rather than skipping the record, writing no positions
#' file at all. pybarrnap has no such limit.
#' @noRd
ITSX_MAX_BP <- 100000

#' Run ITSx on a fasta and parse region positions
#'
#' @param complement NULL keeps ITSx's own default (`--complement T`, searches
#'   both strands and reports sense-frame coordinates regardless of input
#'   orientation). FALSE adds `--complement F`, for callers that have already
#'   oriented the input and want "no coordinates" rather than mirrored ones on
#'   an unanchored antisense record.
#' @noRd
run_itsx <- function(fasta, cpu = 2, taxon = "all", complement = NULL) {
  bin <- unname(Sys.which("ITSx"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  pre <- tempfile("rp_itsx_")
  err <- tempfile("rp_itsx_err_")
  on.exit(unlink(err), add = TRUE)
  comp_arg <- if (is.null(complement)) NULL else c("--complement", if (complement) "T" else "F")
  st <- suppressWarnings(system2(bin, c("-i", shQuote(fasta), "-o", shQuote(pre),
    "--save_regions", "none", "--graphical", "F", "--cpu", cpu,
    "-t", taxon, "--preserve", "T", comp_arg),
    stdout = FALSE, stderr = err))
  tail_msg <- function() {
    tryCatch(paste(utils::tail(readLines(err, warn = FALSE), 3L), collapse = " | "),
      error = function(e) "")
  }
  if (!identical(as.integer(st), 0L)) {
    warning(sprintf("ITSx failed on %s: %s", basename(fasta), tail_msg()), call. = FALSE)
    return(data.frame())
  }
  pf <- paste0(pre, ".positions.txt")
  if (!file.exists(pf)) {
    warning(sprintf("ITSx wrote no positions file for %s: %s", basename(fasta), tail_msg()),
      call. = FALSE)
    return(data.frame())
  }
  parse_itsx(readLines(pf, warn = FALSE))
}

#' Call rRNA anchors with the pipeline's tool (pybarrnap preferred, barrnap
#' fallback)
#'
#' pybarrnap is the pipeline default (`annotate_opts.rrna_tool`) and, unlike the
#' Perl barrnap on some installs, reliably calls eukaryote 18S/5.8S/28S. Both
#' emit the same GFF3 `Name=` tags, so [parse_barrnap()] handles either.
#' @return data.frame(element, start, stop, strand) or an empty data.frame.
#' @noRd
call_rrna_anchors <- function(fasta, kingdom = "euk") {
  py <- unname(Sys.which("pybarrnap"))
  if (nzchar(py)) {
    gff <- tryCatch(run_rrna_tool(py, "pybarrnap", fasta,
      rrna_kingdom_for(kingdom, "pybarrnap")), error = function(e) character(0))
    out <- parse_barrnap(gff)
    if (!is.null(out) && nrow(out) > 0) return(out)
  }
  if (nzchar(unname(Sys.which("barrnap")))) {
    return(tryCatch(run_barrnap(fasta, kingdom = kingdom), error = function(e) data.frame()))
  }
  data.frame()
}

#' Orient a stitched nrDNA sequence to the sense (5'->3') strand
#'
#' The rDNA repeat is transcribed 5'-18S-ITS1-5.8S-ITS2-28S-3'; an assembly can
#' come out in either orientation. Calls the rRNA anchors and, if the conserved
#' anchors are predominantly on the minus strand, reverse-complements the whole
#' sequence so downstream annotation/curation/export see a consistent sense
#' sequence with plus-strand features. Idempotent: a sense sequence is returned
#' unchanged. A no-op when no anchor caller is available or none are found, in
#' which case a residual reverse orientation is still caught downstream by
#' curation's `wrong_strand` flag.
#'
#' @param seq Character scalar nucleotide sequence.
#' @param kingdom anchor-caller kingdom (default "euk").
#' @return list(seq = possibly reverse-complemented sequence, flipped = logical).
#' @noRd
orient_rrna_sense <- function(seq, kingdom = "euk") {
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) {
    return(list(seq = seq, flipped = FALSE))
  }
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  writeLines(c(">seq", seq), fa)
  bar <- call_rrna_anchors(fa, kingdom = kingdom)
  anch <- if (!is.null(bar) && nrow(bar) > 0 && "strand" %in% names(bar)) {
    bar$strand[bar$element %in% RRNA_GENES & bar$strand %in% c("+", "-")]
  } else {
    character(0)
  }
  # Flip only on a clear minus-strand majority; a tie or no anchors is left as-is.
  if (length(anch) > 0 && sum(anch == "-") > sum(anch == "+")) {
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
    return(list(seq = rc, flipped = TRUE))
  }
  list(seq = seq, flipped = FALSE)
}

#' Rotate a circular rDNA unit so it starts in the IGS
#'
#' Position 1 becomes the middle of the stretch from the end of 28S to the
#' start of the next 18S. Anchors are called on the unit written twice, so a
#' gene cut by the current start point is still seen whole once. Returned
#' unchanged when no plus-strand 28S precedes an 18S.
#' @param seq Sense-strand sequence of one circular unit.
#' @param kingdom Anchor-caller kingdom (default "euk").
#' @return The rotated sequence, same length.
#' @noRd
rotate_to_igs <- function(seq, kingdom = "euk") {
  L <- nchar(seq)
  fa <- tempfile(fileext = ".fasta")
  on.exit(unlink(fa), add = TRUE)
  writeLines(c(">seq", paste0(seq, seq)), fa)
  bar <- call_rrna_anchors(fa, kingdom = kingdom)
  if (is.null(bar) || nrow(bar) == 0) return(seq)
  bar <- bar[bar$strand %in% "+", , drop = FALSE]
  end28 <- bar$stop[bar$element == "28S"]
  s18 <- bar[bar$element == "18S" & vapply(bar$start, function(x) any(end28 < x), logical(1)), , drop = FALSE]
  if (nrow(s18) == 0) return(seq)
  g <- s18[which.max(s18$stop - s18$start), ]
  e <- max(end28[end28 < g$start])
  s <- (e + (g$start - e - 1L) %/% 2L) %% L + 1L
  if (s == 1L) seq else paste0(substr(seq, s, L), substr(seq, 1L, s - 1L))
}

#' Sequence ids short enough for ITSx
#' @param lens Named integer vector of sequence lengths.
#' @noRd
itsx_eligible <- function(lens) {
  names(lens)[!is.na(lens) & lens <= ITSX_MAX_BP]
}

#' Map a seqid/start/stop frame between a sequence and its reverse complement
#'
#' Self-inverse: applying this twice with the same `lens` returns the original
#' coordinates. Used both to move barrnap's input-frame anchors into the
#' revcomp frame ITSx reports on for a minus-strand record, and to move a
#' finished reconciled model back out of that frame.
#' @param d A data.frame with seqid, start, stop columns (any other columns
#'   are left untouched).
#' @param lens Named integer vector of input-frame sequence lengths.
#' @noRd
flip_coords <- function(d, lens) {
  if (is.null(d) || nrow(d) == 0) return(d)
  L <- lens[d$seqid]
  ok <- !is.na(L)
  if (!any(ok)) return(d)
  s <- d$start; e <- d$stop
  d$start[ok] <- as.integer(L[ok] - e[ok] + 1L)
  d$stop[ok] <- as.integer(L[ok] - s[ok] + 1L)
  d
}

#' Reconcile per record
#' @param bar,its parse_barrnap()/parse_itsx() output carrying a seqid column.
#' @return reconcile_features() output with a leading seqid column.
#' @noRd
reconcile_records <- function(bar, its) {
  ids <- unique(c(
    if (!is.null(bar) && nrow(bar) > 0) bar$seqid else character(0),
    if (!is.null(its) && nrow(its) > 0) its$seqid else character(0)))
  if (length(ids) == 0) return(data.frame())
  out <- list()
  drop_seqid <- function(d) {
    if (is.null(d) || nrow(d) == 0) return(data.frame())
    d[, setdiff(names(d), "seqid"), drop = FALSE]
  }
  for (sid in ids) {
    b <- drop_seqid(if (!is.null(bar) && nrow(bar) > 0) bar[bar$seqid == sid, , drop = FALSE] else NULL)
    i <- drop_seqid(if (!is.null(its) && nrow(its) > 0) its[its$seqid == sid, , drop = FALSE] else NULL)
    r <- reconcile_features(b, i)
    if (is.null(r) || nrow(r) == 0) next
    out[[length(out) + 1]] <- cbind(seqid = sid, r, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

#' Annotate every record in a multi-record fasta
#'
#' Runs the rRNA anchor caller once over the whole file and ITSx once per strand
#' group, then reconciles per record. Records longer than [ITSX_MAX_BP] are kept
#' but get anchors only, because ITSx aborts the entire run on an oversized
#' target rather than skipping it.
#'
#' barrnap reports coordinates in the input frame; on a genuinely minus-strand
#' record those anchors run 28S -> 5.8S -> 18S in ascending position, the
#' opposite of the ascending order [reconcile_features()]'s spacer-clipping
#' assumes. ITSx, run on a reverse-complemented copy of those records, reports
#' its regions in the frame where the anchors DO ascend. So minus-strand
#' records are reconciled in that revcomp frame -- barrnap's anchors are
#' flipped into it, ITSx's output is used as reported -- and the finished
#' model is flipped back to the input frame afterward.
#'
#' @param fasta Path to a multi-record fasta.
#' @param kingdom Anchor-caller kingdom.
#' @param cpu Threads for ITSx.
#' @param taxon ITSx `-t` profile set.
#' @return data.frame(seqid, element, type, start, stop, length, strand, model_cov, tool).
#' @export
annotate_records <- function(fasta, kingdom = "euk", cpu = 2, taxon = "all") {
  seqs <- Biostrings::readDNAStringSet(fasta)
  if (length(seqs) == 0) return(data.frame())
  names(seqs) <- sub("\\s.*$", "", names(seqs))
  lens <- stats::setNames(as.integer(Biostrings::width(seqs)), names(seqs))

  caller <- rrna_caller()
  bar <- data.frame()
  if (!is.null(caller)) {
    gff <- run_rrna_tool(caller$bin, caller$tool, fasta,
      rrna_kingdom_for(kingdom, caller$tool), cpu = cpu)
    bar <- parse_barrnap(gff)
  }

  minus <- if (nrow(bar) > 0) unique(bar$seqid[bar$strand == "-"]) else character(0)
  plus <- setdiff(names(seqs), minus)
  eligible <- itsx_eligible(lens)

  run_on <- function(ids, revcomp) {
    if (length(ids) == 0) return(data.frame())
    sub <- seqs[ids]
    if (revcomp) sub <- Biostrings::reverseComplement(sub)
    f <- tempfile(fileext = ".fa")
    on.exit(unlink(f), add = TRUE)
    Biostrings::writeXStringSet(sub, f)
    # both groups are already oriented by construction; sense-frame-only
    # coordinates from ITSx are correct here (see run_itsx() complement doc)
    run_itsx(f, cpu = cpu, taxon = taxon, complement = FALSE)
  }
  bar_by <- function(ids) if (nrow(bar) == 0) data.frame() else bar[bar$seqid %in% ids, , drop = FALSE]

  its_plus <- run_on(intersect(plus, eligible), FALSE)
  out_plus <- reconcile_records(bar_by(plus), its_plus)

  its_minus <- run_on(intersect(minus, eligible), TRUE)     # ITSx's native revcomp frame
  bar_minus_rc <- flip_coords(bar_by(minus), lens)           # barrnap -> revcomp frame
  out_minus <- flip_coords(reconcile_records(bar_minus_rc, its_minus), lens) # -> input frame
  # minus-set membership is barrnap's call; overrides reconcile's "+" default
  if (nrow(out_minus) > 0) out_minus$strand <- "-"

  rbind(out_plus, out_minus)
}
