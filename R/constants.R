# Canonical nuclear ribosomal repeat (nrDNA) element names, shared across the
# package (annotation reconciliation, feature-table export, app selectors).
# The repeat is 5'-[ETS]-18S-ITS1-5.8S-ITS2-28S-[IGS]-3', tandemly repeated.

#' Standard nuclear rRNA genes (the conserved anchors)
#' @noRd
RRNA_GENES <- c("18S", "5.8S", "28S")

#' Internal transcribed spacers (the variable, boundary-defined products)
#' @noRd
RRNA_SPACERS <- c("ITS1", "ITS2")

#' Canonical 5'->3' order of elements in the rDNA repeat unit. Used to order features
#' and derive spacers by subtraction between conserved anchors.
#' @noRd
RRNA_ELEMENT_ORDER <- c("ETS", "18S", "ITS1", "5.8S", "ITS2", "28S", "IGS")

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
  source_contig = paste("What this contig was called in the assembly you supplied, so a",
    "scaffold here can be traced back to your own file."),
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
  screen_reason = paste("Why the screen check reached that verdict, including when the panel",
    "found no usable hit."),
  variant_class = paste("Whether this scaffold looks like an alternative arrangement of another",
    "scaffold in the sample, or a second copy of the rDNA. Information only:",
    "neither one is a fault."),
  screen_species = "Organism named by the best GenBank match, searched against all of GenBank.",
  shares = paste("The most specific level of classification the best GenBank match shares with",
    "this sample. \"no reference\" is normal for a group with little in GenBank."),
  odd_one_out = paste("Marks a scaffold whose best GenBank match disagrees with the rest of the",
    "sample's scaffolds, which is the usual signature of contamination."),
  screen_accession = "GenBank accession of the best match. Click it to open the record at NCBI.",
  screen_pident = paste("How closely the matched stretch resembles the best GenBank record,",
    "as a percentage of identical bases."),
  triage_reason = "Why triage decided the way it did, or what left it undecided.",
  chimera = paste("Chimera check: whether different parts of this scaffold match different",
    "lineages or show read-level breakpoints. clean, suspect, or likely. Hover for the",
    "lineage segments and junction positions. Warns only; never blocks export.")
)

#' Toggle the annotation lock on one unit (used by the bulk lock button and by
#' the "lock anyway" answer to the two-paths warning).
#' @noRd
LOCK_UNIT_SQL <- "UPDATE annotate SET annotate_lock = 1 - COALESCE(annotate_lock, 0) WHERE ID = ?"

#' Column help for the chimera segment and window tables (annotate modal)
#' @noRd
CHIMERA_COL_HELP <- list(
  start = "First base of this stretch of the scaffold (1-based).",
  stop = "Last base of this stretch of the scaffold.",
  taxon = "Best-hit lineage at class rank, shared by every window in the segment. 'nohit' means no window in the stretch matched the local rDNA panel or decoys.",
  species = "Species of the strongest single hit inside the segment.",
  pident = "Mean percent identity of the windows' best hits.",
  n_windows = "Number of 500 bp windows (step 250) collapsed into this segment.",
  win_species = "Species of the best hit for this window.",
  lineage = "Best-hit lineage, last three ranks shown; hover for the full ranked lineage.",
  win_pident = "Percent identity of the window's best hit.",
  bitscore = "BLAST bit score of the best hit.",
  rank_distance = "Ranks between the hit and this sample's own lineage: 0 = same species, 7 = nothing shared.",
  depth_median = "Median read depth across the window.",
  clip_frac = "Largest fraction of reads soft-clipped at one position in the window; high values mark a read-level breakpoint."
)
