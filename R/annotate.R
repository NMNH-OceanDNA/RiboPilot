# Annotation core: predict the nuclear rDNA feature model (18S / ITS1 / 5.8S /
# ITS2 / 28S) by combining barrnap (Rfam SSU/LSU anchors) and ITSx (spacer
# boundaries), then reconciling into one gap-free, 5'->3' feature table.
#
# Split so the tools (barrnap, ITSx) can run in the WF2 Nextflow workflow while R
# does the parsing + reconciliation: parse_barrnap()/parse_itsx() consume tool
# output; reconcile_features() combines them; run_barrnap()/run_itsx()/
# annotate_seq() are the all-in-R convenience path.

#' Parse barrnap GFF3 text into rRNA gene features
#' @return data.frame(element, start, stop, strand) for 18S/5.8S/28S
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
  nm <- vapply(f, function(x) sub(".*Name=([^;]+).*", "\\1", x[9]), character(1))
  map <- c("18S_rRNA" = "18S", "5_8S_rRNA" = "5.8S", "5.8S_rRNA" = "5.8S", "28S_rRNA" = "28S")
  el <- unname(map[nm])
  keep <- !is.na(el)
  if (!any(keep)) {
    return(data.frame())
  }
  data.frame(
    element = el[keep],
    start = as.integer(vapply(f, `[`, character(1), 4))[keep],
    stop = as.integer(vapply(f, `[`, character(1), 5))[keep],
    strand = vapply(f, `[`, character(1), 7)[keep],
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
  add <- function(element, start, stop, tool) {
    if (is.null(start) || is.null(stop) || is.na(start) || is.na(stop) || stop < start) {
      return(invisible())
    }
    feats[[length(feats) + 1]] <<- data.frame(
      element = element, type = rrna_element_type(element),
      start = as.integer(start), stop = as.integer(stop), tool = tool, stringsAsFactors = FALSE)
  }

  b18 <- pick(bar, "18S"); i18 <- pick(its, "18S")
  b58 <- pick(bar, "5.8S"); i58 <- pick(its, "5.8S")
  b28 <- pick(bar, "28S"); i28 <- pick(its, "28S")
  its1 <- pick(its, "ITS1"); its2 <- pick(its, "ITS2")

  s18 <- b18 %||% i18
  s58 <- i58 %||% b58
  s28 <- b28 %||% i28

  if (!is.null(s18)) add("18S", s18$start, s18$stop, if (!is.null(b18)) "barrnap" else "ITSx")
  if (!is.null(its1)) {
    add("ITS1", its1$start, its1$stop, "ITSx")
  } else if (!is.null(s18) && !is.null(s58)) {
    add("ITS1", s18$stop + 1, s58$start - 1, "derived")
  }
  if (!is.null(s58)) add("5.8S", s58$start, s58$stop, if (!is.null(i58)) "ITSx" else "barrnap")
  if (!is.null(its2)) {
    add("ITS2", its2$start, its2$stop, "ITSx")
  } else if (!is.null(s58) && !is.null(s28)) {
    add("ITS2", s58$stop + 1, s28$start - 1, "derived")
  }
  if (!is.null(s28)) add("28S", s28$start, s28$stop, if (!is.null(b28)) "barrnap" else "ITSx")

  if (length(feats) == 0) {
    return(data.frame())
  }
  out <- do.call(rbind, feats)
  out$length <- out$stop - out$start + 1
  out$strand <- "+"
  out <- out[order(out$start), c("element", "type", "start", "stop", "length", "strand", "tool")]
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

#' Annotate the nuclear rDNA feature model on a single sequence (all in R)
#'
#' @param fasta Path to a single-sequence fasta (e.g. a stitched path).
#' @param cpu CPUs for ITSx.
#' @return data.frame(element, type, start, stop, length, strand, tool).
#' @export
annotate_seq <- function(fasta, cpu = 2) {
  reconcile_features(run_barrnap(fasta), run_itsx(fasta, cpu = cpu))
}
