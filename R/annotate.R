# Annotation core: predict the nuclear rDNA feature model (18S / ITS1 / 5.8S /
# ITS2 / 28S) by combining barrnap (Rfam SSU/LSU anchors) and ITSx (spacer
# boundaries), then reconciling into one gap-free, 5'->3' feature table.
#
# Split so the tools (barrnap, ITSx) can run in the WF2 Nextflow workflow while R
# does the parsing + reconciliation: parse_barrnap()/parse_itsx() consume tool
# output; reconcile_features() combines them; run_barrnap()/run_itsx()/
# annotate_seq() are the all-in-R convenience path.

#' Parse barrnap GFF3 text into rRNA gene features
#'
#' Also captures the fraction of the rRNA model that aligned: (py)barrnap marks a
#' short hit `product=... (partial)` and adds `Note=aligned only NN percent of
#' the ...`. That coverage is what tells a complete 18S from one missing its 3'
#' V9 region (the terminus most metazoan eDNA assays target), so it is kept
#' rather than discarded.
#' @return data.frame(element, start, stop, strand, model_cov) for 18S/5.8S/28S,
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
    element = el[keep],
    start = as.integer(vapply(f, `[`, character(1), 4))[keep],
    stop = as.integer(vapply(f, `[`, character(1), 5))[keep],
    strand = vapply(f, `[`, character(1), 7)[keep],
    model_cov = model_cov,
    stringsAsFactors = FALSE
  )
}

#' Parse an ITSx `.positions.txt` into region features
#' @return data.frame(element, start, stop) for SSU/ITS1/5.8S/ITS2/LSU found
#' @noRd
parse_itsx <- function(pos_lines) {
  if (length(pos_lines) == 0) {
    return(data.frame())
  }
  fields <- strsplit(pos_lines[1], "\t")[[1]]
  labmap <- c("SSU" = "18S", "ITS1" = "ITS1", "5.8S" = "5.8S", "ITS2" = "ITS2", "LSU" = "28S")
  out <- list()
  for (fd in fields) {
    m <- regmatches(fd, regexec("^\\s*(SSU|ITS1|5\\.8S|ITS2|LSU):\\s*(\\d+)-(\\d+)", fd))[[1]]
    if (length(m) == 4) {
      out[[length(out) + 1]] <- data.frame(
        element = unname(labmap[m[2]]), start = as.integer(m[3]), stop = as.integer(m[4]),
        stringsAsFactors = FALSE)
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

  b18 <- pick(bar, "18S"); i18 <- pick(its, "18S")
  b58 <- pick(bar, "5.8S"); i58 <- pick(its, "5.8S")
  b28 <- pick(bar, "28S"); i28 <- pick(its, "28S")
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

  out$length <- out$stop - out$start + 1
  out <- out[order(out$start),
    c("element", "type", "start", "stop", "length", "strand", "model_cov", "tool")]
  rownames(out) <- NULL
  out
}

#' Run barrnap (Rfam eukaryote models) on a fasta and parse rRNA features
#' @noRd
run_barrnap <- function(fasta, kingdom = "fun") {
  bin <- unname(Sys.which("barrnap"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  gff <- suppressWarnings(system2(bin, c("--kingdom", kingdom, "--quiet", shQuote(fasta)),
    stdout = TRUE, stderr = FALSE))
  parse_barrnap(gff)
}

#' Run ITSx on a fasta and parse region positions
#' @noRd
run_itsx <- function(fasta, cpu = 2) {
  bin <- unname(Sys.which("ITSx"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  pre <- tempfile("rp_itsx_")
  suppressWarnings(system2(bin, c("-i", shQuote(fasta), "-o", shQuote(pre),
    "--save_regions", "none", "--graphical", "F", "--cpu", cpu, "-t", "all", "--preserve", "T"),
  stdout = FALSE, stderr = FALSE))
  pf <- paste0(pre, ".positions.txt")
  if (!file.exists(pf)) {
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
    gff <- tryCatch(suppressWarnings(system2(py, c("-k", kingdom, "-q", shQuote(fasta)),
      stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
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

#' Annotate the nuclear rDNA feature model on a single sequence (all in R)
#'
#' @param fasta Path to a single-sequence fasta (e.g. a stitched path).
#' @param cpu CPUs for ITSx.
#' @return data.frame(element, type, start, stop, length, strand, tool).
#' @export
annotate_seq <- function(fasta, cpu = 2) {
  reconcile_features(run_barrnap(fasta), run_itsx(fasta, cpu = cpu))
}
