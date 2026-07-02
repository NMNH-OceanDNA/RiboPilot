# Level-1 scaffold stitching.
#
# GetOrganelle nr often fragments one rDNA repeat unit into several scaffolds (e.g. an 18S
# piece + two 28S pieces). Before scoring/annotation we order and orient a path's
# scaffolds along a reference rDNA sequence into one canonical 5'->3'
# 18S-ITS1-5.8S-ITS2-28S sequence, inserting gap runs where scaffolds don't abut.
# This mirrors MitoPilot's scaffold_join, driven here by BLAST to a reference.

#' Order + orient a path's scaffolds along a reference and stitch them
#'
#' @param path_fasta Fasta of one path's scaffolds.
#' @param reference_fasta Fasta with a single reference rDNA sequence (e.g. the top
#'   BLAST hit's full 18S-ITS-28S sequence).
#' @param max_gap Cap (bp) on the N-gap inserted between non-abutting scaffolds
#'   (default 2000); larger reference gaps are clamped to this.
#' @param min_hit Minimum BLAST alignment length (bp) to place a scaffold
#'   (default 100); unplaced scaffolds are appended at the end and flagged.
#' @param min_new Minimum new reference coverage (bp) a scaffold must add to be
#'   kept in the tiling; scaffolds contained within a larger one are dropped as
#'   redundant (default 200).
#' @return list(sequence = character stitched sequence,
#'   layout = data.frame(scaffold, ref_start, ref_end, strand, placed, length,
#'   redundant, order), n_gaps, n_unplaced, n_redundant)
#' @export
stitch_scaffolds <- function(path_fasta, reference_fasta, max_gap = 2000, min_hit = 100, min_new = 200) {
  if (!nzchar(unname(Sys.which("blastn"))) || !nzchar(unname(Sys.which("makeblastdb")))) {
    stop("blastn/makeblastdb not found (see environment.yml).", call. = FALSE)
  }
  scaffolds <- Biostrings::readDNAStringSet(path_fasta)
  names(scaffolds) <- sub("\\s.*", "", names(scaffolds))

  db <- tempfile("rp_ref_")
  on.exit(unlink(Sys.glob(paste0(db, "*"))), add = TRUE)
  system2("makeblastdb", c("-in", shQuote(reference_fasta), "-dbtype", "nucl", "-out", shQuote(db)),
    stdout = FALSE, stderr = FALSE)

  out <- suppressWarnings(system2("blastn", c(
    "-query", shQuote(path_fasta), "-db", shQuote(db),
    "-outfmt", shQuote("6 qseqid qlen sstart send length bitscore")
  ), stdout = TRUE, stderr = FALSE))

  # all HSPs per scaffold
  hsp <- if (length(out) > 0) {
    f <- strsplit(out, "\t")
    data.frame(
      scaffold = vapply(f, `[`, character(1), 1),
      sstart = as.integer(vapply(f, `[`, character(1), 3)),
      send = as.integer(vapply(f, `[`, character(1), 4)),
      length = as.integer(vapply(f, `[`, character(1), 5)),
      bitscore = as.numeric(vapply(f, `[`, character(1), 6)),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  # build a layout row per scaffold using its FULL reference footprint (all HSPs)
  layout <- data.frame(
    scaffold = names(scaffolds),
    ref_start = NA_integer_, ref_end = NA_integer_, strand = NA_character_,
    placed = FALSE, length = as.integer(Biostrings::width(scaffolds)),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(layout))) {
    h <- hsp[hsp$scaffold == layout$scaffold[i], ]
    if (nrow(h) == 0 || sum(h$length) < min_hit) next
    # footprint = span of all HSPs on the reference
    layout$ref_start[i] <- min(h$sstart, h$send)
    layout$ref_end[i] <- max(h$sstart, h$send)
    # strand from the highest-bitscore HSP
    top <- h[which.max(h$bitscore), ]
    layout$strand[i] <- if (top$send >= top$sstart) "+" else "-"
    layout$placed[i] <- TRUE
  }

  # greedy non-redundant tiling: keep scaffolds (largest footprint first) only
  # if they add >= min_new bp of new reference coverage; drop the rest as
  # redundant (e.g. a short 28S fragment contained in a near-complete scaffold).
  layout$redundant <- FALSE
  placed_all <- layout[layout$placed, ]
  placed_all <- placed_all[order(-(placed_all$ref_end - placed_all$ref_start)), ]
  covered <- list()
  new_cov <- function(a, b) {
    total <- b - a + 1
    for (iv in covered) {
      ov <- max(0, min(b, iv[2]) - max(a, iv[1]) + 1)
      total <- total - ov
    }
    max(total, 0)
  }
  keep <- logical(nrow(placed_all))
  for (i in seq_len(nrow(placed_all))) {
    if (new_cov(placed_all$ref_start[i], placed_all$ref_end[i]) >= min_new) {
      keep[i] <- TRUE
      covered <- c(covered, list(c(placed_all$ref_start[i], placed_all$ref_end[i])))
    }
  }
  redundant_scaffs <- placed_all$scaffold[!keep]
  layout$redundant[layout$scaffold %in% redundant_scaffs] <- TRUE

  placed <- placed_all[keep, ]
  placed <- placed[order(placed$ref_start), ]
  unplaced <- layout[!layout$placed, ]

  # orient + concatenate placed scaffolds with reference-sized (capped) N gaps
  pieces <- character(0)
  n_gaps <- 0L
  prev_end <- NA_integer_
  for (i in seq_len(nrow(placed))) {
    s <- scaffolds[[placed$scaffold[i]]]
    if (placed$strand[i] == "-") s <- Biostrings::reverseComplement(s)
    if (!is.na(prev_end)) {
      gap <- placed$ref_start[i] - prev_end - 1
      if (gap > 0) {
        pieces <- c(pieces, strrep("N", min(gap, max_gap)))
        n_gaps <- n_gaps + 1L
      }
    }
    pieces <- c(pieces, as.character(s))
    prev_end <- placed$ref_end[i]
  }
  # append any unplaced scaffolds (flagged) after a gap
  for (i in seq_len(nrow(unplaced))) {
    pieces <- c(pieces, strrep("N", 100), as.character(scaffolds[[unplaced$scaffold[i]]]))
  }

  # assemble a full layout: kept (in concat order), then redundant, then unplaced
  order_map <- setNames(seq_len(nrow(placed)), placed$scaffold)
  out_layout <- layout
  out_layout$order <- unname(order_map[out_layout$scaffold])
  out_layout <- out_layout[order(is.na(out_layout$order), out_layout$order,
    !out_layout$placed, out_layout$scaffold), ]

  list(
    sequence = paste(pieces, collapse = ""),
    layout = out_layout,
    n_gaps = n_gaps,
    n_unplaced = nrow(unplaced),
    n_redundant = sum(out_layout$redundant)
  )
}
