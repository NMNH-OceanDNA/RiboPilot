# Canonical nuclear ribosomal repeat (nrDNA) element names, shared across the
# package (annotation reconciliation, feature-table export, app selectors).
# The repeat is 5'-[ETS]-18S-ITS1-5.8S-ITS2-28S-[IGS]-3', tandemly repeated.

#' Standard nuclear rRNA genes (the conserved anchors)
#' @noRd
RRNA_GENES <- c("18S", "5.8S", "28S")

#' Internal transcribed spacers (the variable, boundary-defined products)
#' @noRd
RRNA_SPACERS <- c("ITS1", "ITS2")

#' Optional external spacers (usually not recovered/curated by default)
#' @noRd
RRNA_EXTRA <- c("ETS", "IGS")

#' Canonical 5'->3' order of elements in the rDNA repeat unit. Used to order features
#' and derive spacers by subtraction between conserved anchors.
#' @noRd
RRNA_ELEMENT_ORDER <- c("ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "IGS")

#' Flat list of element names for use as selectInput choices.
#' @noRd
RRNA_ELEMENT_CHOICES <- c(RRNA_GENES, RRNA_SPACERS, RRNA_EXTRA)

#' Default GetOrganelle nrDNA reference databases (broad Actinopterygii).
#'
#' Single source of truth for the packaged defaults; these are only resolvable
#' once the RiboPilot repository is public. Non-fish projects should build their
#' own with [build_reference_db()] rather than rely on these.
#' @noRd
DEFAULT_SEEDS_DB <- paste0("https://raw.githubusercontent.com/dmacguigan/RiboPilot/",
  "main/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta")

#' @noRd
DEFAULT_LABELS_DB <- paste0("https://raw.githubusercontent.com/dmacguigan/RiboPilot/",
  "main/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta")

#' Executors that run processes inside a container image.
#'
#' Everything except the conda profile, which uses the active environment's
#' tools and sets no `container` directive.
#' @noRd
CONTAINER_EXECUTORS <- c("docker", "local", "awsbatch", "slurm", "pbs", "sge",
                         "lsf", "NMNH_Hydra", "NOAA_SEDNA")

#' Titles naming a nuclear rDNA locus.
#'
#' Broadened from a literal-phrase list because GenBank writes the same locus
#' many ways: "28S ribosomal RNA (ACC), rRNA" carries no "gene", and a
#' chromosome-level host record may name the locus in passing. Safe to include
#' the bare "rRNA gene" construction only because NON_NUCLEAR_RRNA_RE is tested
#' first.
#'
#' The ITS token is case-sensitive (the inline `(?-i:...)`) because the rest of
#' the pattern runs case-insensitively, where it would match the English "its".
#' @noRd
NUCLEAR_RRNA_TITLE_RE <- paste0(
  "internal transcribed spacer|external transcribed spacer|(?-i:\\bITS[12]?\\b)|",
  "\\b(18S|28S|5\\.8S|5S)\\b[^;]{0,40}?(ribosomal RNA|rRNA)|",
  "small subunit ribosomal RNA|large subunit ribosomal RNA|",
  "\\b(SSU|LSU)\\b[^;]{0,20}?rRNA|",
  "\\b45S\\b|",
  "ribosomal RNA gene|\\brRNA gene\\b|ribosomal DNA|\\brDNA\\b")

#' Titles naming a ribosomal subunit that is not nuclear.
#'
#' 16S and 23S are prokaryotic; 12S is the metazoan mitochondrial SSU, whose
#' records often carry no organelle word for [ORGANELLAR_TITLE_RE] to catch.
#' @noRd
NON_NUCLEAR_RRNA_RE <- "\\b(12S|16S|23S)\\b"

#' Nuclear rDNA subunit names, used to spare a record that names both.
#'
#' Deliberately excludes the spacer tokens: bacterial 16S-23S internal
#' transcribed spacer records are common in core_nt and must stay vetoed.
#' @noRd
NUCLEAR_RRNA_MARKER_RE <- "\\b(18S|28S|5\\.8S)\\b"

#' Titles describing an assembly rather than a locus.
#'
#' A chromosome record is not evidence against rDNA, it is no evidence either
#' way. Several markers are matched because real deflines contain typos:
#' "unplaced-saffold_25" occurs in GenBank.
#' @noRd
ASSEMBLY_TITLE_RE <- paste0(
  "chromosome|unplaced|scaffold|contig|whole genome shotgun|\\bWGS\\b|haplotype")

#' Titles indicating an organellar (non-nuclear) record.
#' @noRd
ORGANELLAR_TITLE_RE <- "mitochondri|chloroplast|plastid|apicoplast|kinetoplast"

#' GenBank feature type for an nrDNA element name.
#'
#' rRNA genes -> "rRNA"; transcribed spacers -> "misc_RNA".
#'
#' @param element element name (e.g. "18S", "ITS1")
#' @noRd
rrna_element_type <- function(element) {
  if (element %in% RRNA_GENES) {
    "rRNA"
  } else {
    "misc_RNA"
  }
}

#' GenBank product / note string for an nrDNA element name.
#' @noRd
rrna_element_product <- function(element) {
  switch(element,
    "18S"  = "18S ribosomal RNA",
    "5.8S" = "5.8S ribosomal RNA",
    "28S"  = "28S ribosomal RNA",
    "ITS1" = "internal transcribed spacer 1",
    "ITS2" = "internal transcribed spacer 2",
    "ETS"  = "external transcribed spacer",
    "IGS"  = "intergenic spacer",
    element
  )
}

#' Hover help for each visible column of the assemble modal's scaffold table
#'
#' One place for the wording so the header tooltips and any prose describing the
#' table cannot drift apart.
#' @noRd
SCAFFOLD_COL_HELP <- list(
  unit = paste("Which assembly this row is: p<path>.s<scaffold>. One candidate path",
    "can carry several scaffolds, and each is kept or dropped on its own."),
  triage = paste("Keep, drop, or leave undecided. This is the only decision that counts:",
    "dropped scaffolds are left out of annotation and export."),
  length = "Scaffold length in base pairs.",
  path_score = paste("Rates the whole candidate path this scaffold came from, from its length,",
    "depth, anchor genes and GenBank hit. The same value repeats across every",
    "scaffold of one path."),
  depth = paste("Median read depth per position (a mean is inflated by the conserved",
    "18S/28S cores). The bar compares it with this sample's deepest scaffold:",
    "meaningful within a path, not across paths."),
  depth_track = paste("Read depth along the scaffold, start to end. A deep dip in the middle",
    "can mark a place where the assembler joined two pieces that do not belong together."),
  elements = "The rRNA and ITS features annotation found, in the order they sit on the scaffold.",
  screen_verdict = paste("Contamination check: whether the best GenBank match looks like the",
    "sample's own organism (on target), something else (off target), or the",
    "evidence is too thin to say."),
  variant_class = paste("Whether this scaffold looks like an alternative arrangement of another",
    "scaffold in the sample, or a second copy of the rDNA. Information only:",
    "neither one is a fault."),
  screen_species = "Organism named by the best GenBank match, searched against all of GenBank.",
  shares = paste("The most specific level of classification the best GenBank match shares with",
    "this sample. \"no reference\" is normal for a group with little in GenBank."),
  odd_one_out = paste("Marks a scaffold whose best GenBank match disagrees with the rest of the",
    "sample's scaffolds, which is the usual signature of contamination."),
  screen_accession = "GenBank accession of the best match. Click it to open the record at NCBI.",
  screen_pident = "How much of the scaffold matched the best GenBank record, as a percentage.",
  triage_reason = "Why triage decided the way it did, or what left it undecided."
)
