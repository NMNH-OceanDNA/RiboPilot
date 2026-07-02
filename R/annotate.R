# Annotation core: predict the nuclear rDNA feature model (18S / ITS1 / 5.8S /
# ITS2 / 28S) on an assembled sequence by combining barrnap (Rfam SSU/LSU
# anchors) and ITSx (spacer boundaries), then reconciling into one gap-free,
# 5'->3' feature table. See docs/03_ANNOTATION_CURATION.md.

#' Run barrnap (Rfam eukaryote models) and parse rRNA gene features
#'
#' The bioconda barrnap fork runs Infernal against Rfam SSU/LSU/5_8S models under
#' `--kingdom fun` (despite the name, these are the eukaryote models).
#' @return data.frame(element, start, stop, strand) for 18S/5.8S/28S found
#' @noRd
run_barrnap <- function(fasta, kingdom = "fun") {
  bin <- unname(Sys.which("barrnap"))
  if (!nzchar(bin)) {
    return(data.frame())
  }
  gff <- suppressWarnings(system2(bin, c("--kingdom", kingdom, "--quiet", shQuote(fasta)),
    stdout = TRUE, stderr = FALSE))
  gff <- gff[!startsWith(gff, "#")]
  if (length(gff) == 0) {
    return(data.frame())
  }
  f <- strsplit(gff, "\t")
  f <- Filter(function(x) length(x) >= 9, f)
  if (length(f) == 0) {
    return(data.frame())
  }
  nm <- vapply(f, function(x) sub(".*Name=([^;]+).*", "\\1", x[9]), character(1))
  map <- c("18S_rRNA" = "18S", "5_8S_rRNA" = "5.8S", "5.8S_rRNA" = "5.8S", "28S_rRNA" = "28S")
  el <- unname(map[nm])
  keep <- !is.na(el)
  data.frame(
    element = el[keep],
    start = as.integer(vapply(f, `[`, character(1), 4))[keep],
    stop = as.integer(vapply(f, `[`, character(1), 5))[keep],
    strand = vapply(f, `[`, character(1), 7)[keep],
    stringsAsFactors = FALSE
  )
}

#' Run ITSx and parse the region positions
#'
#' @return data.frame(element, start, stop) for SSU/ITS1/5.8S/ITS2/LSU found
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
  line <- readLines(pf, warn = FALSE)
  if (length(line) == 0) {
    return(data.frame())
  }
  # e.g. "seqid\t12057 bp.\tSSU: 1-1831\tITS1: 1832-2707\t5.8S: 2708-2865\tITS2: 2866-3414\tLSU: 3415-12057\t..."
  fields <- strsplit(line[1], "\t")[[1]]
  labmap <- c("SSU" = "18S", "ITS1" = "ITS1", "5.8S" = "5.8S", "ITS2" = "ITS2", "LSU" = "28S")
  out <- list()
  for (fd in fields) {
    m <- regmatches(fd, regexec("^\\s*(SSU|ITS1|5\\.8S|ITS2|LSU):\\s*(\\d+)-(\\d+)", fd))[[1]]
    if (length(m) == 4) {
      out[[length(out) + 1]] <- data.frame(
        element = unname(labmap[m[2]]), start = as.integer(m[3]), stop = as.integer(m[4]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) data.frame() else do.call(rbind, out)
}

#' Annotate the nuclear rDNA feature model on a single assembled sequence
#'
#' Combines barrnap (18S/5.8S/28S anchors; bounds the 28S extent) and ITSx
#' (ITS1/5.8S/ITS2 boundaries), preferring each tool where it is authoritative,
#' then derives any missing spacer by subtraction between adjacent anchors so the
#' model is gap-free and ordered 5'->3'.
#'
#' @param fasta Path to a single-sequence nucleotide fasta (e.g. a stitched
#'   path). If multiple sequences are present, only the first is annotated.
#' @param cpu CPUs for ITSx.
#' @return data.frame(element, type, start, stop, length, strand, tool) ordered
#'   by start, or an empty data.frame if nothing was found.
#' @export
annotate_seq <- function(fasta, cpu = 2) {
  bar <- run_barrnap(fasta)
  its <- run_itsx(fasta, cpu = cpu)

  pick <- function(df, el) if (nrow(df) > 0 && any(df$element == el)) df[df$element == el, ][1, ] else NULL
  feats <- list()
  add <- function(element, start, stop, tool) {
    if (is.null(start) || is.null(stop) || is.na(start) || is.na(stop) || stop < start) {
      return(invisible())
    }
    feats[[length(feats) + 1]] <<- data.frame(
      element = element, type = rrna_element_type(element),
      start = as.integer(start), stop = as.integer(stop), tool = tool,
      stringsAsFactors = FALSE
    )
  }

  # anchors: barrnap authoritative for 18S/28S extent; ITSx as fallback
  b18 <- pick(bar, "18S")
  i18 <- pick(its, "18S")
  b58 <- pick(bar, "5.8S")
  i58 <- pick(its, "5.8S")
  b28 <- pick(bar, "28S")
  i28 <- pick(its, "28S")
  its1 <- pick(its, "ITS1")
  its2 <- pick(its, "ITS2")

  s18 <- b18 %||% i18
  s58 <- i58 %||% b58 # ITSx more reliable for the tiny 5.8S
  s28 <- b28 %||% i28

  if (!is.null(s18)) add("18S", s18$start, s18$stop, if (!is.null(b18)) "barrnap" else "ITSx")
  # ITS1: ITSx, else derive between 18S and 5.8S
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
