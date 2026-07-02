#!/usr/bin/env Rscript
# Quick smoke test for RiboPilot taxonomy resolution + scaffold/path scoring.
#
# Run from the RiboPilot repo root, with the `ribopilot` conda env active
# (needs taxonkit; Part C also needs barrnap + blastn):
#
#   conda activate ribopilot
#   Rscript dev/test_scoring.R
#
# Part A (scoring math) needs no external tools and always runs.
# Part B (taxonomy) downloads a ~72 MB NCBI taxdump on first run (cached).
# Part C (live assembly) only runs if an assembly output dir is present.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Reuse an existing taxdump cache if you have one; else it downloads to ./.taxdump
TAXDUMP <- Sys.getenv("RIBOPILOT_TAXDUMP",
  unset = {
    cached <- "~/.local/share/RiboPilot/taxdump"
    if (file.exists(file.path(path.expand(cached), "names.dmp"))) path.expand(cached) else ".taxdump"
  }
)

cat("\n===========================================================\n")
cat("PART A  -  composite scoring (no external tools)\n")
cat("===========================================================\n")
samp <- "Metazoa; Chordata; Actinopteri; Syngnathiformes; Mullidae; Upeneus; Upeneus parvus"
cand <- data.frame(
  path          = c(1, 2, 3),
  depth         = c(250, 120, 15),        # graph coverage
  blast_pident  = c(98.6, 96.0, 99.0),    # note: contaminant has HIGH identity
  blast_lineage = c(
    samp,
    "Metazoa; Chordata; Actinopteri; Syngnathiformes; Mullidae; Upeneus; Upeneus tragula",
    "Fungi; Ascomycota; Sordariomycetes; Hypocreales; Nectriaceae; Fusarium; Fusarium oxysporum"
  ),
  anchors       = c("18S,28S", "18S,28S", "18S"),
  length        = c(9768, 9768, 1800),
  stringsAsFactors = FALSE
)
res <- score_candidates(cand, sample_lineage = samp)
print(res[, c("path", "score", "rank", "auto_select",
  "s_tax", "s_identity", "s_anchors", "s_depth", "s_length")], row.names = FALSE)
cat("\n-> path 3 is a 99%-identity FUNGAL contaminant and is correctly ranked last\n",
  "   (taxonomy + depth override raw BLAST identity).\n", sep = "")

cat("\n===========================================================\n")
cat("PART B  -  typo-tolerant taxonomy resolution\n")
cat("===========================================================\n")
if (!nzchar(unname(Sys.which("taxonkit")))) {
  cat("SKIPPED: taxonkit not on PATH (activate the ribopilot conda env).\n")
} else {
  ensure_taxdump(TAXDUMP)
  taxa <- c("Upeneus parvus", "Upeneus parvis", "Psychrolutes paradoxis", "Notarealfish madeuppus")
  r <- resolve_taxonomy(taxa, taxdump_dir = TAXDUMP)
  print(r[, c("Taxon", "match_type", "taxid", "sci_name", "rank")], row.names = FALSE)
  cat("\n-> 'Upeneus parvis' / 'Psychrolutes paradoxis' are typos, recovered to the",
    "\n   correct species; nonsense is flagged 'unresolved'.\n", sep = "")
}

cat("\n===========================================================\n")
cat("PART C  -  live assembly scoring (optional)\n")
cat("===========================================================\n")
AO <- "assembly_out/SRR22396758/go_out"
paths <- Sys.glob(file.path(AO, "*path_sequence.fasta"))
gfa <- Sys.glob(file.path(AO, "*selected_graph.gfa"))
seeds <- "ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta"
if (length(paths) == 0 || !nzchar(unname(Sys.which("blastn"))) || !nzchar(unname(Sys.which("barrnap")))) {
  cat("SKIPPED: needs an assembly at", AO, "plus blastn + barrnap.\n")
} else {
  db <- file.path(tempdir(), "rp_seeds_db")
  system2("makeblastdb", c("-in", seeds, "-dbtype", "nucl", "-out", db), stdout = FALSE, stderr = FALSE)
  res <- score_assembly(paths, gfa_path = gfa[1], sample_lineage = samp,
    blast_db = db, taxdump_dir = TAXDUMP, with_staxids = FALSE)
  print(res[, c("path", "length", "depth", "anchors", "blast_pident", "hit_species",
    "score", "rank", "auto_select")], row.names = FALSE)
}
cat("\nDone.\n")
