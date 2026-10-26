# Shared DDL for tables new_db() creates, kept here so the statement lives in one
# place when more than one caller needs it.

#' DDL for blast_hits: the full retained BLAST hit distribution per unit.
#'
#' One row per (subject, HSP). The pipeline already pays for a remote core_nt
#' search; everything a contamination call needs -- the hit distribution, the
#' on/off-target bitscore margin, and the positional structure that separates a
#' chimera from a clean scaffold -- is kept rather than collapsed to one row.
#' @noRd
BLAST_HITS_DDL <- "CREATE TABLE IF NOT EXISTS blast_hits (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  hit_rank INTEGER NOT NULL,
  hsp INTEGER NOT NULL DEFAULT 1,
  accession TEXT,
  title TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  aln_len INTEGER,
  qstart INTEGER,
  qend INTEGER,
  sstart INTEGER,
  send INTEGER,
  slen INTEGER,
  is_nuclear_rrna INTEGER,
  hit_class TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, hit_rank, hsp)
);"

#' DDL for blast_searched: which units have been searched, and by which search.
#'
#' Keyed on a hash of the query sequence and the search settings (blast.nf), so
#' a re-assembled or trimmed unit, or a changed BLAST option, is searched again
#' rather than keeping the previous hits. A unit that found nothing still gets a
#' row, which is what stops it being re-searched on every retry and every re-run.
#' `remote_fallback_failed` marks a unit whose escalation to NCBI failed, so the
#' next run tries it again.
#' `db` is the database that answered: `core_nt` for a remote search, the local
#' panel's name and release (from its manifest) for a local one, "" when a local
#' panel has no manifest.
#' @noRd
BLAST_SEARCHED_DDL <- "CREATE TABLE IF NOT EXISTS blast_searched (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  query_md5 TEXT NOT NULL,
  n_hits INTEGER,
  source TEXT,
  db TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for the triage_opts table (per-rule action presets).
#' @noRd
TRIAGE_OPTS_DDL <- "CREATE TABLE IF NOT EXISTS triage_opts (
  triage_opts            TEXT NOT NULL,
  rule_not_annotated     TEXT,
  rule_no_screen         TEXT,
  rule_no_feat_non_rdna  TEXT,
  rule_no_feat_rdna_hit  TEXT,
  rule_feat_non_rdna_hit TEXT,
  rule_odd_one_out       TEXT,
  rule_screen_off_target TEXT,
  rule_chimeric          TEXT,
  PRIMARY KEY (triage_opts)
);"

#' DDL for the screen table (machine contamination verdict, one row per unit).
#' @noRd
SCREEN_DDL <- "CREATE TABLE IF NOT EXISTS screen (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  screen_verdict    TEXT,
  screen_class      TEXT,
  screen_confidence REAL,
  screen_reason     TEXT,
  screen_version    INTEGER,
  lca_taxid         INTEGER,
  lca_taxon         TEXT,
  lca_rank          TEXT,
  lca_n_hits        INTEGER,
  host_concordant   INTEGER,
  margin            REAL,
  variant_class     TEXT,
  variant_group     TEXT,
  variant_identity  REAL,
  time_stamp        INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for chimera: one row per unit, lineage kept as text at every level.
#' @noRd
CHIMERA_DDL <- "CREATE TABLE IF NOT EXISTS chimera (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  score REAL,
  label TEXT,
  junctions TEXT,
  reason TEXT,
  segments TEXT,
  gene_taxa TEXT,
  channels TEXT,
  ribotyper_codes TEXT,
  db_source TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for chimera_windows: best local hit per sliding window.
#' @noRd
CHIMERA_WINDOWS_DDL <- "CREATE TABLE IF NOT EXISTS chimera_windows (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  start INTEGER NOT NULL,
  stop INTEGER,
  accession TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  hit_class TEXT,
  pident REAL,
  bitscore REAL,
  depth_median REAL,
  clip_frac REAL,
  db_source TEXT,
  PRIMARY KEY (ID, path, scaffold, start)
);"

#' DDL for feature_blast: one local panel BLAST call per core feature of a unit.
#' @noRd
FEATURE_BLAST_DDL <- "CREATE TABLE IF NOT EXISTS feature_blast (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  element TEXT NOT NULL,
  start INTEGER,
  stop INTEGER,
  length INTEGER,
  task TEXT,
  low_conf INTEGER,
  status TEXT,
  note TEXT,
  accession TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  n_tied INTEGER,
  match_rank TEXT,
  rank_distance INTEGER,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  db_source TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, element)
);"

#' DDL for feature_blast_hits: top hits per feature, one row per subject.
#' @noRd
FEATURE_BLAST_HITS_DDL <- "CREATE TABLE IF NOT EXISTS feature_blast_hits (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  element TEXT NOT NULL,
  hit_rank INTEGER NOT NULL,
  accession TEXT,
  title TEXT,
  staxid INTEGER,
  lineage TEXT,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  hit_class TEXT,
  PRIMARY KEY (ID, path, scaffold, element, hit_rank)
);"

#' DDL for annotations: the nrDNA feature model per unit.
#'
#' status/flags are filled by curation (curate_one); tool/edited by annotation.
#' partial_start/partial_stop: 5'/3' end partial in the feature's own orientation.
#' @noRd
ANNOTATIONS_DDL <- "CREATE TABLE IF NOT EXISTS annotations (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  element TEXT NOT NULL,
  type TEXT,
  start INTEGER,
  stop INTEGER,
  length INTEGER,
  strand TEXT,
  model_coverage REAL,
  partial_start INTEGER,
  partial_stop INTEGER,
  tool TEXT,
  status TEXT,
  flags TEXT,
  edited INTEGER,
  notes TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, element)
);"

#' DDL for blast_screen: the unfiltered top remote hit the contaminant check reads.
#' @noRd
BLAST_SCREEN_DDL <- "CREATE TABLE IF NOT EXISTS blast_screen (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  accession TEXT,
  title TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  is_nuclear_rrna INTEGER,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for blast_ref_alignment: per-feature sample-vs-reference alignment.
#' @noRd
BLAST_REF_ALIGNMENT_DDL <- "CREATE TABLE IF NOT EXISTS blast_ref_alignment (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  accession TEXT NOT NULL,
  element TEXT NOT NULL,
  pid REAL,
  aln_cov REAL,
  aln_len INTEGER,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, accession, element)
);"

#' DDL for the reference table (one row per clade reference-database group).
#' @noRd
REFERENCE_DDL <- "CREATE TABLE IF NOT EXISTS reference (
  ref_id           TEXT NOT NULL,
  clade_taxid      INTEGER,
  clade_name       TEXT,
  clade_rank       TEXT,
  query            TEXT,
  n_samples        INTEGER,
  n_records        INTEGER,
  reference_switch INTEGER,
  source           TEXT,
  seeds_db         TEXT,
  labels_db        TEXT,
  log_file         TEXT,
  pid              INTEGER,
  preset_name      TEXT,
  fail_reason      TEXT,
  time_stamp       INTEGER,
  PRIMARY KEY (ref_id)
);"

#' DDL for reference_config: one row of project-level reference settings.
#'
#' `ref_dir` is the folder built databases are written to and reused from; NULL
#' means fall back to RIBOPILOT_REF_DIR and then ribopilot_ref_dir(). The record
#' window bounds the clade suggestion walk.
#' @noRd
REFERENCE_CONFIG_DDL <- "CREATE TABLE IF NOT EXISTS reference_config (
  rowid_key   INTEGER NOT NULL,
  ref_dir     TEXT,
  min_records INTEGER,
  max_records INTEGER,
  blastdb     TEXT,
  blastdb_version TEXT,
  blastdb_info TEXT,
  PRIMARY KEY (rowid_key)
);"

#' DDL for rrna_candidates: the rDNA contig search evidence for user assemblies.
#'
#' One row per contig barrnap hit in a user-supplied assembly, kept or dropped,
#' so the app can say what the search found rather than only that it failed.
#' Contigs with no rRNA hit at all are not recorded: a draft genome would
#' otherwise write a row per scaffold.
#' @noRd
RRNA_CANDIDATES_DDL <- "CREATE TABLE IF NOT EXISTS rrna_candidates (
  ID TEXT NOT NULL,
  contig TEXT NOT NULL,
  length INTEGER,
  rrna_features TEXT,
  rrna_bases INTEGER,
  rank INTEGER,
  selected INTEGER,
  reason TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, contig)
);"

META_RECORDS_DDL <- "CREATE TABLE IF NOT EXISTS meta_records (
  ID TEXT NOT NULL, source TEXT NOT NULL, level TEXT NOT NULL, depth INTEGER NOT NULL,
  ref TEXT, field TEXT NOT NULL, value TEXT, PRIMARY KEY (ID, source, depth, field));"

META_STATUS_DDL <- "CREATE TABLE IF NOT EXISTS meta_status (
  ID TEXT NOT NULL, source TEXT NOT NULL, ref TEXT, status TEXT NOT NULL, message TEXT,
  fetched_at INTEGER, PRIMARY KEY (ID, source));"

META_EXPORT_FIELDS_DDL <- "CREATE TABLE IF NOT EXISTS meta_export_fields (
  key TEXT NOT NULL, PRIMARY KEY (key));"

META_CSV_MAP_DDL <- "CREATE TABLE IF NOT EXISTS meta_csv_map (
  concept TEXT NOT NULL, column TEXT, PRIMARY KEY (concept));"

META_OPTIONS_DDL <- "CREATE TABLE IF NOT EXISTS meta_options (
  key TEXT NOT NULL, value TEXT, PRIMARY KEY (key));"

META_LINKS_DDL <- "CREATE TABLE IF NOT EXISTS meta_links (
  ID TEXT NOT NULL, source TEXT NOT NULL, ref TEXT, via TEXT, note TEXT,
  PRIMARY KEY (ID, source));"

META_VIEW_FIELDS_DDL <- "CREATE TABLE IF NOT EXISTS meta_view_fields (
  key TEXT NOT NULL, shown INTEGER NOT NULL, PRIMARY KEY (key));"

META_SAMPLE_COLS <- c("GEOME_BCID", "GBIF_ID", "BioSample")

#' Create the sample metadata tables and add missing samples columns
#' @noRd
meta_apply_schema <- function(con) {
  for (ddl in c(META_RECORDS_DDL, META_STATUS_DDL, META_EXPORT_FIELDS_DDL, META_CSV_MAP_DDL,
                META_OPTIONS_DDL, META_LINKS_DDL, META_VIEW_FIELDS_DDL)) {
    DBI::dbExecute(con, ddl)
  }
  have <- DBI::dbListFields(con, "samples")
  for (col in META_SAMPLE_COLS) {
    if (col %in% have) next
    hit <- have[tolower(have) == tolower(col)]
    if (length(hit)) {
      DBI::dbExecute(con, paste0("ALTER TABLE samples RENAME COLUMN \"", hit[1], "\" TO ", col))
    } else {
      DBI::dbExecute(con, paste0("ALTER TABLE samples ADD COLUMN ", col, " TEXT"))
    }
  }
  .nmnh_ensure_table(con)
  invisible()
}

#' DDL for review_units: a contig's review copy; a row here means it has one.
#'
#' trim_start/trim_end: the kept window in original coordinates. check_switch:
#' 2 done, 1 not run, 4 running, 3 failed, 0 on hold. rev starts at the time in
#' milliseconds and goes up by one on every invalidation.
#' @noRd
REVIEW_UNITS_DDL <- "CREATE TABLE IF NOT EXISTS review_units (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL,
  scaffold INTEGER NOT NULL,
  sequence TEXT,
  length INTEGER,
  trim_start INTEGER,
  trim_end INTEGER,
  edited INTEGER NOT NULL DEFAULT 0,
  check_switch INTEGER,
  check_note TEXT,
  rev INTEGER NOT NULL,
  checked_at INTEGER,
  n_features INTEGER,
  missing TEXT,
  chimeric INTEGER,
  n_flags INTEGER,
  curate_flags TEXT,
  itsx_warnings TEXT,
  itsx_origin TEXT,
  itsx_scores TEXT,
  blast_accession TEXT,
  blast_species TEXT,
  blast_pident REAL,
  blast_qcovs REAL,
  blast_evalue REAL,
  blast_source TEXT,
  time_stamp INTEGER,
  genbank_check TEXT,
  PRIMARY KEY (ID, path, scaffold)
);"

#' Review mirror of each original per-contig table
#' @noRd
REVIEW_MIRRORS <- c(annotations = "review_annotations", blast_hits = "review_blast_hits",
  blast_screen = "review_blast_screen", screen = "review_screen", chimera = "review_chimera",
  chimera_windows = "review_chimera_windows", blast_ref_alignment = "review_ref_alignment",
  feature_blast = "review_feature_blast", feature_blast_hits = "review_feature_blast_hits")

#' The review tables, review_units first
#' @noRd
REVIEW_TABLES <- c("review_units", unname(REVIEW_MIRRORS))

#' An original's DDL with the table renamed to its mirror
#' @noRd
review_mirror_ddl <- function(ddl, orig) {
  sub(paste0("CREATE TABLE IF NOT EXISTS ", orig, " ("),
    paste0("CREATE TABLE IF NOT EXISTS ", REVIEW_MIRRORS[[orig]], " ("), ddl, fixed = TRUE)
}

REVIEW_ANNOTATIONS_DDL <- sub("  PRIMARY KEY",
  "  partial_start_user INTEGER,\n  partial_stop_user INTEGER,\n  PRIMARY KEY",
  review_mirror_ddl(ANNOTATIONS_DDL, "annotations"), fixed = TRUE)
REVIEW_BLAST_HITS_DDL <- review_mirror_ddl(BLAST_HITS_DDL, "blast_hits")
REVIEW_BLAST_SCREEN_DDL <- review_mirror_ddl(BLAST_SCREEN_DDL, "blast_screen")
REVIEW_SCREEN_DDL <- review_mirror_ddl(SCREEN_DDL, "screen")
REVIEW_CHIMERA_DDL <- review_mirror_ddl(CHIMERA_DDL, "chimera")
REVIEW_CHIMERA_WINDOWS_DDL <- review_mirror_ddl(CHIMERA_WINDOWS_DDL, "chimera_windows")
REVIEW_REF_ALIGNMENT_DDL <- review_mirror_ddl(BLAST_REF_ALIGNMENT_DDL, "blast_ref_alignment")
REVIEW_FEATURE_BLAST_DDL <- review_mirror_ddl(FEATURE_BLAST_DDL, "feature_blast")
REVIEW_FEATURE_BLAST_HITS_DDL <- review_mirror_ddl(FEATURE_BLAST_HITS_DDL, "feature_blast_hits")

#' Column names of a CREATE TABLE constant, in order
#' @noRd
ddl_cols <- function(ddl) {
  x <- sub("[ ,].*", "", trimws(strsplit(ddl, "\n", fixed = TRUE)[[1]][-1]))
  x[grepl("^[A-Za-z_]", x) & x != "PRIMARY"]
}

REVIEW_COPY_SQL <- "EXISTS (SELECT 1 FROM review_units r
  WHERE r.ID = x.ID AND r.path = x.path AND r.scaffold = x.scaffold)"

#' CREATE VIEW: mirror rows of contigs with a copy, original rows of the rest
#' @noRd
review_view_ddl <- function(view, orig, cols, fallback = paste0("x.", cols)) {
  paste0("CREATE VIEW IF NOT EXISTS ", view, " AS\nSELECT ", paste0("x.", cols, collapse = ", "),
    "\nFROM ", REVIEW_MIRRORS[[orig]], " x\nWHERE ", REVIEW_COPY_SQL,
    "\nUNION ALL\nSELECT ", paste(fallback, collapse = ", "),
    "\nFROM ", orig, " x\nWHERE NOT ", REVIEW_COPY_SQL)
}

#' The review layer Annotation Review and export read
#'
#' A contig nobody touched reads its originals. In review_units_v such a contig
#' reads as not run (1) when it is done (2) but has no stored ITSx value.
#' @noRd
REVIEW_VIEWS_DDL <- c(
  review_units_v = paste0("CREATE VIEW IF NOT EXISTS review_units_v AS
SELECT x.ID, x.path, x.scaffold, 1 AS has_copy, a.topology, a.length AS orig_length,
  x.sequence, x.length, x.trim_start, x.trim_end, x.edited, x.check_switch, x.check_note,
  x.rev, x.checked_at, x.n_features, x.missing, x.chimeric, x.n_flags, x.curate_flags,
  x.itsx_warnings, x.itsx_origin, x.itsx_scores, x.blast_accession, x.blast_species,
  x.blast_pident, x.blast_qcovs, x.blast_evalue, x.blast_source, x.time_stamp, x.genbank_check
FROM review_units x
JOIN assemblies a ON a.ID = x.ID AND a.path = x.path AND a.scaffold = x.scaffold
UNION ALL
SELECT x.ID, x.path, x.scaffold, 0, x.topology, x.length,
  x.sequence, x.length, 1, x.length, 0,
  CASE WHEN an.annotate_switch = 2 AND an.itsx_warnings IS NULL THEN 1 ELSE an.annotate_switch END,
  NULL, NULL, NULL, an.n_features, an.missing, an.chimeric, an.n_flags, an.curate_flags,
  an.itsx_warnings, an.itsx_origin, an.itsx_scores, x.blast_accession, x.blast_species,
  x.blast_pident, x.blast_qcovs, x.blast_evalue, x.blast_source, an.time_stamp, an.genbank_check
FROM assemblies x
JOIN annotate an ON an.ID = x.ID AND an.path = x.path AND an.scaffold = x.scaffold
WHERE NOT ", REVIEW_COPY_SQL),
  review_annotations_v = review_view_ddl("review_annotations_v", "annotations",
    c(ddl_cols(ANNOTATIONS_DDL), "partial_start_user", "partial_stop_user"),
    c(paste0("x.", ddl_cols(ANNOTATIONS_DDL)), "NULL AS partial_start_user", "NULL AS partial_stop_user")),
  review_screen_v = review_view_ddl("review_screen_v", "screen", ddl_cols(SCREEN_DDL)),
  review_chimera_v = review_view_ddl("review_chimera_v", "chimera", ddl_cols(CHIMERA_DDL)),
  review_chimera_windows_v = review_view_ddl("review_chimera_windows_v", "chimera_windows",
    ddl_cols(CHIMERA_WINDOWS_DDL)))

#' Create the review tables, then the views over them (idempotent)
#' @noRd
review_apply_schema <- function(con) {
  views <- REVIEW_VIEWS_DDL
  # released step 4 runs this before step 6 adds genbank_check; step 6 rebuilds the view
  if (!"genbank_check" %in% DBI::dbListFields(con, "annotate")) {
    views[["review_units_v"]] <- gsub("(x|an)\\.genbank_check", "NULL AS genbank_check", views[["review_units_v"]])
  }
  for (ddl in c(REVIEW_UNITS_DDL, REVIEW_ANNOTATIONS_DDL, REVIEW_BLAST_HITS_DDL,
                REVIEW_BLAST_SCREEN_DDL, REVIEW_SCREEN_DDL, REVIEW_CHIMERA_DDL,
                REVIEW_CHIMERA_WINDOWS_DDL, REVIEW_REF_ALIGNMENT_DDL, views)) {
    DBI::dbExecute(con, ddl)
  }
  invisible()
}

FEATURE_BLAST_VIEWS_DDL <- c(
  review_feature_blast_v = review_view_ddl("review_feature_blast_v", "feature_blast",
    ddl_cols(FEATURE_BLAST_DDL)),
  review_feature_blast_hits_v = review_view_ddl("review_feature_blast_hits_v", "feature_blast_hits",
    ddl_cols(FEATURE_BLAST_HITS_DDL)))

#' Create feature_blast, feature_blast_hits, their review mirrors and views
#'
#' Kept apart from review_apply_schema(), which released migration step 4 calls.
#' @noRd
feature_blast_apply_schema <- function(con) {
  for (ddl in c(FEATURE_BLAST_DDL, FEATURE_BLAST_HITS_DDL, REVIEW_FEATURE_BLAST_DDL,
                REVIEW_FEATURE_BLAST_HITS_DDL, FEATURE_BLAST_VIEWS_DDL)) {
    DBI::dbExecute(con, ddl)
  }
  invisible()
}
