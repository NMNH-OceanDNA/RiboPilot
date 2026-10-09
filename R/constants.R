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

#' Features a unit needs before a batch Trim will cut it.
#' @noRd
BATCH_TRIM_FEATURES <- setdiff(RRNA_ELEMENT_ORDER, c("ETS", "IGS"))

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
CONTAINER_EXECUTORS <- c("docker", "local", "slurm", "pbs", "sge", "lsf",
                         "NMNH_Hydra", "NOAA_SEDNA")

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
#' records often carry no organelle word for `ORGANELLAR_TITLE_RE` to catch.
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

#' Hover help for each visible column of the assemble modal's contig table
#'
#' One place for the wording so the header tooltips and any prose describing the
#' table cannot drift apart.
#' @noRd
SCAFFOLD_COL_HELP <- list(
  unit = paste("Which assembly this row is: p<path>.s<scaffold>. A path is one candidate",
    "assembly of the sample; a contig is one sequence in it. One path can carry",
    "several contigs, and each is kept or dropped on its own."),
  triage = paste("Triage is the keep-or-drop call on a contig: keep, drop, or leave undecided.",
    "This is the only decision that counts:",
    "dropped contigs are left out of annotation and export."),
  length = "Contig length in base pairs.",
  source_contig = paste("What this contig was called in the assembly you supplied, so a",
    "contig here can be traced back to your own file."),
  path_score = paste("Rates the whole candidate path this contig came from, from its length,",
    "depth, anchor genes and BLAST hit. The same value repeats across every",
    "contig of one path."),
  depth = paste("Median read depth per position (a mean is inflated by the conserved",
    "18S/28S cores). The bar compares it with this sample's deepest contig:",
    "meaningful within a path, not across paths."),
  depth_track = paste("Read depth along the contig, start to end. Reads that map equally well",
    "to two places are left out, so a repeat inside the contig shows as a dip. A lasting step",
    "in depth, not a short dip, can mark a place where the assembler joined two pieces that do",
    "not belong together."),
  elements = paste("The rRNA and ITS features annotation found, in the order they sit on the contig.",
    "Orange features have curation warnings; hover to list them."),
  itsx_warnings = paste("What stock ITSx (the ITS detection tool; default settings, no organism group, the way NCBI's",
    "Submission Portal runs it) flags on this contig: chimera, gap inside the ITS,",
    "missing or partial 5.8S, a boundary far from RiboPilot's own, or an organism group",
    "outside the sample's lineage, shown with its score next to the sample's own group",
    "('tie' or 'near tie' means that group call is unreliable). Filled in by",
    "the annotation step, and by Update on Annotation Review after a trim or an edit.",
    "Empty means clean; 'not run' means ITSx gave no answer (hover for why);",
    "'not checked' means no ITSx check is stored yet."),
  screen_verdict = paste("BLAST screen: whether the best matches look like the",
    "sample's own organism (on target), something else (off target), or the",
    "evidence is too thin to say."),
  screen_reason = paste("Why the BLAST screen reached that verdict, including when the bundled",
    "rDNA database found no usable hit."),
  variant_class = paste("Whether this contig looks like an alternative arrangement of another",
    "contig in the sample, or a second copy of the rDNA. Information only:",
    "neither one is a fault."),
  screen_species = paste("Organism named by the best match of the BLAST screen: the bundled rDNA",
    "database first, then NCBI core_nt for a contig it cannot match. The evidence panel below says which."),
  shares = paste("The most specific level of classification the best GenBank match shares with",
    "this sample. \"no reference\" is normal for a group with little in GenBank."),
  odd_one_out = paste("Marks a contig whose best GenBank match disagrees with the rest of the",
    "sample's contigs, which is the usual signature of contamination."),
  screen_accession = "GenBank accession of the best match. Click it to open the record at NCBI.",
  screen_pident = paste("How closely the matched stretch resembles the best GenBank record,",
    "as a percentage of identical bases."),
  triage_reason = "Why triage decided the way it did, or what left it undecided.",
  warnings = paste("Number of curation warnings on the contig as a whole, such as a missing gene,",
    "features out of order, or contamination. Hover for the list. Warnings on one",
    "feature, such as length or truncated, are in the contig's Curation warnings panel."),
  chimera = paste("Chimera check: whether different parts of this contig match different",
    "lineages or show read-level breakpoints. severe: parts match different phyla or",
    "kingdoms; the export flags the record as chimeric and asks about it. minor: weaker",
    "evidence, listed in the export report only. clean: nothing found. Hover for the evidence.")
)

#' Column help for the chimera segment and window tables (annotate modal)
#' @noRd
CHIMERA_COL_HELP <- list(
  start = "First base of this stretch of the contig (1-based).",
  stop = "Last base of this stretch of the contig.",
  taxon = "Best-hit lineage at class rank, shared by every window in the segment. 'nohit' means no window in the stretch matched the bundled rDNA database. 'unresolved' means the top hits tie across classes (a conserved stretch), so the windows say nothing about origin.",
  species = "The strongest single hit inside the segment: a species, or the clade tied hits share.",
  pident = "Mean percent identity of the windows' best hits.",
  n_windows = "Number of 500 bp windows (step 250) collapsed into this segment.",
  win_species = "The window's best hit: a species, or, when several hits score about the same, the clade they share and the number of tied hits. The lineage is cut back to that clade.",
  lineage = "Best-hit lineage, last three ranks shown; hover for the full ranked lineage.",
  win_pident = "Percent identity of the window's best hit.",
  bitscore = "BLAST bit score of the best hit.",
  rank_distance = "Ranks between the hit and this sample's own lineage: 0 = same species, 7 = nothing shared.",
  depth_median = "Median read depth across the window.",
  position = "Estimated join position (bp), on the same scale as the plot.",
  checks = "Checks that put a junction near this position: genes (the 18S, 5.8S, and 28S BLAST hits), windows (500 bp BLAST windows), and reads (read depth and clipped reads).",
  why = "What each check saw here.",
  clip_frac = "Largest fraction of reads soft-clipped at one position in the window; high values mark a read-level breakpoint."
)

# ---------------------------------------------------------------------------
# Shared UI vocabularies (state registry, column names, colors, lock wording)
# ---------------------------------------------------------------------------
# constants.R collates last, so read these inside app_ui()/app_server() only;
# never derive a top-level constant from them in another R/ file.

#' Pipeline state registry. Order is lifecycle order, not numeric order.
#' @noRd
MP_STATE_META <- list(
  `5` = list(label = "Waiting for Assemble Update", icon = "fa fa-clock",
             tip = "The main pipeline has not finished this contig. Press Update on the Assemble step first."),
  `0` = list(label = "On hold",      icon = "fa fa-hourglass",
             tip = "Skipped by the next run until you set it to Ready to run."),
  `1` = list(label = "Ready to run", icon = "fa fa-play",
             tip = "Queued: the next run will process it."),
  `4` = list(label = "In progress",  icon = "fa fa-circle-half-stroke",
             tip = "Being processed by the running pipeline, or queued for the next run."),
  `2` = list(label = "Success",      icon = "fa fa-circle-check",
             tip = "Finished without errors."),
  `3` = list(label = "Failed",       icon = "fa fa-triangle-exclamation",
             tip = "Failed. Open the sample's details to see why and what to try next.")
)

#' Codes each module can display. Annotation Review's 5 is display-only.
#' @noRd
MP_STATE_CODES <- list(
  taxonomy = c("0", "1", "4", "2", "3"),
  assemble = c("0", "1", "4", "2", "3"),
  annotate = c("5", "0", "1", "4", "2", "3")
)

#' Codes a user may set by hand.
#' @noRd
MP_STATE_SETTABLE <- c("0", "1", "4", "2", "3")

#' Named code vector for a state picker: names are labels, values are codes.
#' @noRd
mp_state_choices <- function(module) {
  mp_state_field(module, "label", named = TRUE)
}

#' Icon classes for a state column, in module code order.
#' @noRd
mp_state_icons <- function(module) mp_state_field(module, "icon")

#' State labels, in module code order.
#' @noRd
mp_state_labels <- function(module) mp_state_field(module, "label")

#' State tooltips, in module code order.
#' @noRd
mp_state_tips <- function(module) mp_state_field(module, "tip")

#' Pull one field of MP_STATE_META for a module, keyed by code.
#' @noRd
mp_state_field <- function(module, field, named = FALSE) {
  codes <- MP_STATE_CODES[[match.arg(module, names(MP_STATE_CODES))]]
  out <- vapply(codes, function(k) MP_STATE_META[[k]][[field]], character(1))
  if (!named) return(out)
  names(codes) <- out
  codes
}

#' State filter choices shared by the pickers and the visible-row mask.
#' @noRd
STATE_CHOICES <- mp_state_choices("assemble")

#' Hex colors that cannot live in CSS because SweetAlert takes them as values.
#' @noRd
MP_COLORS <- c(primary = "#337ab7", danger = "#b02a37", grey = "#6c757d")

#' Canonical table column headers, keyed by data-frame column name.
#' @noRd
MP_COL_NAMES <- c(
  taxonomy_lock    = "Lock",
  assemble_lock    = "Lock",
  annotate_lock    = "Lock",
  taxonomy_switch  = "State",
  reference_switch = "State",
  assemble_switch  = "State",
  annotate_switch  = "State",
  ID               = "ID",
  path             = "Path #",
  scaffold         = "Contig #",
  Taxon            = "Taxon",
  sci_name         = "Assigned Taxon",
  taxid            = "TaxID",
  rank             = "Rank",
  match_type       = "Confidence",
  reviewed         = "Reviewed",
  problematic      = "Excluded",
  partial          = "Partial",
  checks           = "Checks",
  n_flags          = "Warnings",
  itsx_warnings    = "ITSx Warnings",
  screen_verdict   = "Screen",
  chimera_label    = "Chimera",
  chimera_reason   = "Chimera Reason",
  reference        = "Reference",
  blast_species    = "BLAST Species",
  blast_accession  = "BLAST Accession",
  blast_pident     = "BLAST % Identity",
  blast_qcovs      = "BLAST % Coverage",
  blast_source     = "Hit Source",
  details          = "Details",
  output           = "Output",
  n_features       = "# Features",
  missing          = "Missing",
  length           = "Length (bp)",
  paths            = "# Paths",
  units            = "# Contigs",
  kept             = "# Kept",
  sel_path         = "Selected Path",
  best_score       = "Best Score",
  export_group     = "Export Group",
  exported         = "Exported",
  pre_opts         = "Preprocess Options",
  assemble_opts    = "Assembly Options",
  blast_opts       = "BLAST Options",
  annotate_opts    = "Annotate Options",
  curate_opts      = "Curate Options",
  notes            = "Notes",
  annotate_notes   = "User Notes",
  user_notes       = "User Notes",
  time_stamp       = "Last Updated",
  ref_group        = "Reference Group",
  user_group       = "User Group",
  set_by           = "Set By",
  status           = "Status",
  clade_name       = "Clade",
  clade_rank       = "Rank",
  n_samples        = "# Samples",
  n_records        = "# GenBank Records",
  database         = "Database",
  fail_reason      = "Note",
  action           = "Action",
  tax_alternatives = "Other Matches"
)

#' Header tooltips for the columns whose name does not explain itself.
#' Same keys as MP_COL_NAMES. Sentence case, no terminal period.
#' @noRd
MP_COL_TIPS <- c(
  taxonomy_switch  = "Pipeline state; hover a row icon for its meaning",
  assemble_switch  = "Pipeline state; hover a row icon for its meaning",
  annotate_switch  = "Pipeline state; hover a row icon for its meaning",
  taxonomy_lock    = "Locked rows are finished with this step",
  assemble_lock    = "Locked rows are finished with this step",
  annotate_lock    = "Locked rows are finished with this step",
  reviewed         = "Checked by a curator; Lock and Export do not ask again about its contamination flags",
  problematic      = "Excluded contigs are left out of every export",
  partial          = "A partial contig is described as 'partial sequence'; gene ends are partial at untrimmed contig ends, in truncated genes, and where a trim cut a gene, unless set with 5'/3' in the contig window",
  path             = "Assembly path number for this contig. * = the sample keeps more than one path; type * in the filter to list them",
  scaffold         = "Contig number within the path. * = the sample keeps more than one contig in this path; type * in the filter to list them",
  sci_name         = "Taxon assigned from the NCBI taxdump",
  match_type       = "How the taxon name was matched",
  checks           = paste("Red: Lock and Export ask about it (off target, contaminant, severe chimera,",
    "ITSx chimeric). Amber: worth a look (minor chimera, uncertain screen, ITSx, curation warnings, feature flags).",
    "ok: on target and nothing flagged. On Assemble, covers the sample's kept contigs"),
  n_flags          = SCAFFOLD_COL_HELP$warnings,
  itsx_warnings    = SCAFFOLD_COL_HELP$itsx_warnings,
  screen_verdict   = SCAFFOLD_COL_HELP$screen_verdict,
  chimera_label    = SCAFFOLD_COL_HELP$chimera,
  pre_opts         = "Named set of read preprocessing parameters",
  assemble_opts    = "Named set of assembly parameters",
  blast_opts       = "Named set of BLAST search parameters",
  annotate_opts    = "Named set of annotation parameters",
  curate_opts      = "Named set of curation parameters",
  kept             = "Contigs kept after triage; 0 means nothing to annotate",
  sel_path         = "Assembly path chosen for annotation",
  best_score       = "Highest contig score in this sample",
  export_group     = "Named group this contig is exported with",
  exported         = "When this contig was last exported",
  time_stamp       = "When this row was last changed",
  n_records        = "GenBank records available for this clade",
  tax_alternatives = "Every NCBI taxon the name matched, with its kingdom; click the ID to pick one"
)

# Lock vocabulary shared by the toolbar tooltip, the Lock column header and
# the details-window banner.
MP_UNIT <- c(taxonomy = "sample", reference = "group", assemble = "sample",
             annotate = "contig", export = "contig")
MP_UNITS <- c(taxonomy = "samples", reference = "groups", assemble = "samples",
              annotate = "contigs", export = "contigs")
MP_NEXT <- c(taxonomy = "Assemble", assemble = "Annotation Review", annotate = "Export")
MP_LOCK_DEF <- function(mod = c("taxonomy", "assemble", "annotate")) {
  mod <- match.arg(mod)
  sprintf(
    paste(
      "Locked %s are finished with this step: their options for this step cannot be changed,",
      "the next update will not re-run them, and they move on to %s.",
      "You can unlock at any time."
    ),
    MP_UNITS[[mod]], MP_NEXT[[mod]]
  )
}

#' Why an annotation-locked unit cannot be edited, and how to change that
#' @noRd
MP_UNIT_LOCKED <- paste("This contig is locked, so it is read only here.",
  "To change it, select it in the Annotation Review table and press Unlock.")

#' Batch Trim / Restore tooltips; `what` is what the table's rows are.
#' @noRd
MP_TRIM_TIP <- function() {
  c(trim = paste0("Trim the selected contigs to their gene ends, or gene ends plus a flank. ",
      "Only contigs with all five features (", paste(BATCH_TRIM_FEATURES, collapse = ", "),
      ") are trimmed; locked contigs are skipped. Trimmed contigs show their checks as not run ",
      "until you press Update."),
    restore = paste0("Restore the selected contigs to their originals: trims, boundary edits, and ",
      "5'/3' choices are deleted. Contigs with no edits and locked contigs are skipped."))
}

#' Curator marks: what each toolbar button and window toggle says it will do,
#' what the toast says it did, and what the mark means for export.
#' @noRd
MP_FLAG_WORDS <- list(
  reviewed = c(set = "Mark Reviewed", clear = "Clear Reviewed",
    did_set = "marked reviewed", did_clear = "cleared of reviewed",
    tip = paste("Record that you checked this contig. Lock and Export do not ask again about",
      "contamination flags on reviewed contigs.")),
  problematic = c(set = "Exclude & Lock", clear = "Include in export",
    did_set = "excluded from export and locked", did_clear = "included in export again",
    tip = paste("Exclude & Lock leaves the contig out of every export and locks it as done.",
      "The contig stays in the project; unlock it and include it again to export it.")),
  partial = c(set = "Mark Partial", clear = "Clear Partial",
    did_set = "marked partial", did_clear = "cleared of partial",
    tip = paste("Export describes a partial contig as 'partial sequence'. Gene ends are",
      "marked partial (< and >) at untrimmed contig ends, in truncated genes, and where",
      "a trim cut a gene; the 5' and 3' buttons in the contig window override that.")))
