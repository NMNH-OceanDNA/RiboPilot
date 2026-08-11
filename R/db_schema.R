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
  is_nuclear_rrna INTEGER,
  hit_class TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, hit_rank, hsp)
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
  PRIMARY KEY (triage_opts)
);"

#' Seed row for triage_opts, matching TRIAGE_RULE_DEFS defaults.
#' @noRd
TRIAGE_OPTS_SEED <- "INSERT INTO triage_opts
  (triage_opts, rule_not_annotated, rule_no_screen, rule_no_feat_non_rdna,
   rule_no_feat_rdna_hit, rule_feat_non_rdna_hit, rule_odd_one_out,
   rule_screen_off_target)
  VALUES ('default', 'flag', 'flag', 'drop', 'flag', 'flag', 'flag', 'flag');"

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
  PRIMARY KEY (rowid_key)
);"

#' Back-fill `coverage.median_depth` from the stored per-position tracks
#'
#' Computed lazily rather than at write time: a depth track is one comma-joined
#' integer per base, so a large project holds hundreds of MB of them and parsing
#' every row would make opening the app crawl. Rows with a NULL median fall back
#' to the mean until this runs.
#'
#' @param con Open project DB connection.
#' @param limit Maximum rows to fill in one call (NULL = all).
#' @return (invisibly) the number of rows filled.
#' @export
backfill_coverage_median <- function(con, limit = NULL) {
  if (!("median_depth" %in% DBI::dbListFields(con, "coverage"))) {
    return(invisible(0L))
  }
  q <- "SELECT ID, path, scaffold, depths FROM coverage
        WHERE median_depth IS NULL AND depths IS NOT NULL AND depths != ''"
  if (!is.null(limit)) q <- paste(q, "LIMIT", as.integer(limit))
  d <- tryCatch(DBI::dbGetQuery(con, q), error = function(e) data.frame())
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  upd <- DBI::dbSendStatement(con,
    "UPDATE coverage SET median_depth = ? WHERE ID = ? AND path = ? AND scaffold = ?")
  on.exit(DBI::dbClearResult(upd), add = TRUE)
  n <- 0L
  for (i in seq_len(nrow(d))) {
    v <- suppressWarnings(as.numeric(strsplit(d$depths[i], ",", fixed = TRUE)[[1]]))
    v <- v[!is.na(v)]
    if (!length(v)) next
    DBI::dbBind(upd, list(round(stats::median(v), 2), d$ID[i], d$path[i], d$scaffold[i]))
    n <- n + 1L
  }
  invisible(n)
}
