#!/usr/bin/env Rscript
# Prototype of the Assemble-tab scaffold-stitch visualization (basis for the
# app's plotting function). Reference-anchored gggenes track: reference cistron
# element blocks (18S/ITS/5.8S/28S) on the bottom lane, then one lane per
# candidate path with each scaffold as an orientation arrow at its reference
# position. Overlapping scaffolds are packed onto tight sub-rows; redundant
# (dropped) scaffolds are greyed. The reference GenBank accession is shown.
#
# Run from repo root in the ribopilot env:
#   Rscript dev/prototype_stitch_viz.R <path1.fasta> [<path2.fasta> ...] \
#       --ref <reference.fasta> --out stitch_viz.png
# (defaults target the SRR22396758 dev assembly if present.)

suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(ggplot2); library(gggenes); library(Biostrings)
})

# --- inputs (edit or pass on the command line) ---
AO <- "assembly_out/SRR22396758/go_out"
path_fastas <- Sys.glob(file.path(AO, "*path_sequence.fasta"))
reference_fasta <- "/tmp/rp_ref.fasta" # a single top-hit reference cistron
out_png <- "dev/stitch_viz.png"

stopifnot(length(path_fastas) > 0, file.exists(reference_fasta))

# reference accession/species + element track (barrnap Rfam on the reference)
refhdr <- names(readDNAStringSet(reference_fasta))[1]
acc <- sub("_.*", "", refhdr)
sp <- gsub("_", " ", sub("^[^_]+_", "", refhdr))
reflab <- paste0("reference\n", acc)

gff <- system2("barrnap", c("--kingdom", "fun", "--quiet", shQuote(reference_fasta)),
  stdout = TRUE, stderr = FALSE)
gff <- gff[!startsWith(gff, "#")]
gf <- do.call(rbind, lapply(strsplit(gff, "\t"), function(p) {
  data.frame(start = as.integer(p[4]), end = as.integer(p[5]),
    name = sub(".*Name=([^;]+).*", "\\1", p[9]), stringsAsFactors = FALSE)
}))
gf$gene <- ifelse(grepl("18S", gf$name), "18S", ifelse(grepl("28S", gf$name), "28S", "5.8S"))
if (any(gf$gene == "18S") && any(gf$gene == "28S")) {
  gf <- rbind(gf, data.frame(start = max(gf$end[gf$gene == "18S"]) + 1,
    end = min(gf$start[gf$gene == "28S"]) - 1, name = "ITS", gene = "ITS"))
}
ref <- data.frame(molecule = reflab, start = gf$start, end = gf$end,
  orientation = 1, cat = gf$gene, lab = gf$gene)

mk <- function(fa, label) {
  st <- stitch_scaffolds(fa, reference_fasta)
  L <- st$layout[st$layout$placed, ]
  data.frame(molecule = label, start = L$ref_start, end = L$ref_end,
    orientation = ifelse(L$strand == "+", 1, -1),
    cat = ifelse(L$redundant, "redundant (dropped)", "scaffold"),
    lab = paste0(ifelse(L$redundant, "~ ", "* "), sub("--.*", "", L$scaffold)))
}
paths <- Map(mk, path_fastas, sprintf("path %d", seq_along(path_fastas)))
df <- do.call(rbind, c(list(ref), paths))

# pack overlapping features per molecule onto tight sub-rows (abutting share a row)
pack <- function(d) {
  d <- d[order(d$start), ]
  lane_end <- numeric(0)
  d$lane <- NA_integer_
  for (i in seq_len(nrow(d))) {
    pl <- FALSE
    for (l in seq_along(lane_end)) {
      if (d$start[i] >= lane_end[l]) { d$lane[i] <- l; lane_end[l] <- d$end[i]; pl <- TRUE; break }
    }
    if (!pl) { lane_end <- c(lane_end, d$end[i]); d$lane[i] <- length(lane_end) }
  }
  d
}
mols <- c(reflab, sprintf("path %d", seq_along(path_fastas)))
STEP <- 0.85; GAPB <- 0.75; y0 <- 0; df$y <- NA_real_; labpos <- c(); lab <- c()
for (m in mols) {
  ord <- order(df$start[df$molecule == m])
  sub <- pack(df[df$molecule == m, ])
  nl <- max(sub$lane)
  df$y[df$molecule == m][ord] <- y0 + (sub$lane - 1) * STEP
  labpos <- c(labpos, y0 + (nl - 1) * STEP / 2); lab <- c(lab, m)
  y0 <- y0 + (nl - 1) * STEP + GAPB + STEP
}

pal <- c("18S" = "#0072B2", "ITS" = "#E69F00", "5.8S" = "#F0E442",
  "28S" = "#009E73", "scaffold" = "#56B4E9", "redundant (dropped)" = "#BBBBBB")

p <- ggplot(df, aes(xmin = start, xmax = end, y = y, fill = cat, forward = orientation == 1, label = lab)) +
  geom_gene_arrow(arrowhead_height = grid::unit(4, "mm"), arrow_body_height = grid::unit(3.2, "mm")) +
  geom_gene_label(align = "left", grow = FALSE, height = grid::unit(3.2, "mm")) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(breaks = labpos, labels = lab) +
  labs(title = "Scaffold stitch layout vs reference cistron",
    subtitle = paste0("reference ", acc, "  (", sp, ")   |   arrow = orientation; * kept, ~ redundant"),
    x = "reference position (bp)", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

ggsave(out_png, p, width = 11, height = 4.5, dpi = 140, bg = "white")
cat("wrote", out_png, "\n")
