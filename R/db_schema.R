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
#' Keyed on a hash of the query sequence, so a re-assembled or trimmed unit is
#' searched again rather than inheriting the previous contig's hits. A unit that
#' found nothing still gets a row, which is what stops it being re-searched on
#' every retry and every re-run.
#' @noRd
BLAST_SEARCHED_DDL <- "CREATE TABLE IF NOT EXISTS blast_searched (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  query_md5 TEXT NOT NULL,
  n_hits INTEGER,
  source TEXT,
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

#' DDL for the assembly_backup table (original pre-trim state, for undo/restore).
#' @noRd
ASSEMBLY_BACKUP_DDL <- "CREATE TABLE IF NOT EXISTS assembly_backup (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  sequence TEXT,
  length INTEGER,
  annotations TEXT,
  coverage TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
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
