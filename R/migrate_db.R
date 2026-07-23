# Backwards-compatible DB migration. Brings a project database created by an
# older RiboPilot up to the current schema. Migrations are additive and
# idempotent (guarded ALTER/CREATE), tracked by SQLite's `PRAGMA user_version`.
# new_db() stamps the current version; migrate_db() upgrades older projects.

#' Current project-DB schema version. Bump when adding a migration step below.
#' @noRd
RIBOPILOT_SCHEMA_VERSION <- 6L

#' DDL for the assembly_backup table (original pre-trim state, for undo/restore).
#' Shared by new_db() and the v6 migration so both agree on the schema.
#' @noRd
ASSEMBLY_BACKUP_DDL <- "CREATE TABLE IF NOT EXISTS assembly_backup (
  ID TEXT NOT NULL,
  sequence TEXT,
  length INTEGER,
  annotations TEXT,
  coverage TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID)
);"

#' The pre-v4 restrictive remote-BLAST entrez_query.
#'
#' Retired because it restricted the search space, not just the ranking, so the
#' contamination screen could only ever see nuclear rDNA. Matched literally so
#' the v4 migration only clears it where the user never customized it.
#' @noRd
LEGACY_ENTREZ_QUERY <- paste(
  "(internal transcribed spacer[Title] OR 18S ribosomal RNA[Title]",
  "OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title]",
  "OR small subunit ribosomal RNA[Title]",
  "OR large subunit ribosomal RNA[Title])",
  "NOT (mitochondrion[Title] OR chloroplast[Title] OR plastid[Title])")

#' Bring a project database up to the current schema (idempotent, additive)
#'
#' Adds any tables/columns introduced after a project was created, without
#' touching existing data. Safe to call repeatedly. New columns are only ever
#' added, never dropped or retyped.
#'
#' @param con Open DBI connection to a project database.
#' @return (invisibly) the schema version after migration.
#' @export
migrate_db <- function(con) {
  ver <- tryCatch(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], error = function(e) 0L)
  has_tbl <- function(tbl) tbl %in% DBI::dbListTables(con)
  has_col <- function(tbl, col) has_tbl(tbl) && col %in% DBI::dbListFields(con, tbl)
  add_col <- function(tbl, col, decl) {
    if (has_tbl(tbl) && !has_col(tbl, col)) {
      DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", tbl, col, decl))
    }
  }

  # --- schema v1: reference-guided curation (curate_opts + status/flags) ---
  if (!has_tbl("curate_opts")) {
    DBI::dbExecute(con,
      "CREATE TABLE curate_opts (
        curate_opts TEXT NOT NULL,
        ruleset TEXT,
        min_anchor_pid REAL,
        check_contaminant INTEGER,
        PRIMARY KEY (curate_opts)
      );")
  }
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM curate_opts")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO curate_opts (curate_opts, ruleset, min_anchor_pid, check_contaminant)
       VALUES ('default', 'metazoa_generic', 85, 1)")
  }
  add_col("annotate", "curate_opts", "TEXT DEFAULT 'default'")
  add_col("annotate", "n_flags", "INTEGER")
  add_col("annotations", "status", "TEXT")
  add_col("annotations", "flags", "TEXT")

  # --- schema v2: editable annotation options (annotate_opts) ---
  if (!has_tbl("annotate_opts")) {
    DBI::dbExecute(con,
      "CREATE TABLE annotate_opts (
        annotate_opts TEXT NOT NULL,
        cpus INTEGER,
        memory INTEGER,
        rrna_tool TEXT,
        rrna_kingdom TEXT,
        itsx_taxon TEXT,
        itsx_opts TEXT,
        PRIMARY KEY (annotate_opts)
      );")
  }
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate_opts")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO annotate_opts (annotate_opts, cpus, memory, rrna_tool,
         rrna_kingdom, itsx_taxon, itsx_opts)
       VALUES ('default', 4, 8, 'pybarrnap', 'euk', 'all', '')")
  }
  add_col("annotate", "annotate_opts", "TEXT DEFAULT 'default'")

  # --- schema v3: in-app taxonomy stage (display + lock/state on samples) ---
  had_tax <- has_col("samples", "taxonomy_lock")
  add_col("samples", "sci_name", "TEXT")
  add_col("samples", "rank", "TEXT")
  add_col("samples", "match_type", "TEXT")
  add_col("samples", "taxonomy_switch", "INTEGER DEFAULT 1")
  add_col("samples", "taxonomy_lock", "INTEGER DEFAULT 0")
  if (!had_tax && has_tbl("samples")) {
    # Grandfather samples from projects that predate the taxonomy stage so they
    # are not hidden from Assemble; mark them assigned + locked.
    DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1, taxonomy_switch = 2")
  }

  # --- schema v4: contamination screen over the unrestricted remote BLAST ---
  if (!has_tbl("blast_screen")) {
    DBI::dbExecute(con,
      "CREATE TABLE blast_screen (
        ID TEXT NOT NULL,
        accession TEXT,
        title TEXT,
        species TEXT,
        staxid INTEGER,
        lineage TEXT,
        pident REAL,
        qcovs REAL,
        evalue REAL,
        is_nuclear_rrna INTEGER,
        time_stamp INTEGER,
        PRIMARY KEY (ID)
      );")
  }
  # ITSx chimera flag surfaced through curation.
  add_col("annotate", "chimeric", "INTEGER")
  # Fraction of the rRNA model an anchor covers (truncation signal).
  add_col("annotations", "model_coverage", "REAL")
  # Fraction of a sample element covered by its reference alignment.
  add_col("blast_ref_alignment", "aln_cov", "REAL")
  # Curator flags beyond `reviewed`: problematic (exclude from export) + partial.
  add_col("annotate", "problematic", "TEXT")
  add_col("annotate", "partial", "TEXT")
  # Retire the restrictive entrez_query, but only where it is still verbatim the
  # old default -- a query the user chose is theirs to keep.
  if (has_tbl("blast_opts")) {
    DBI::dbExecute(con,
      "UPDATE blast_opts SET entrez_query = '' WHERE TRIM(COALESCE(entrez_query,'')) = TRIM(?)",
      params = list(LEGACY_ENTREZ_QUERY))
    # 5 was tuned for a pre-filtered search; the client-side filter needs depth.
    DBI::dbExecute(con,
      "UPDATE blast_opts SET max_target_seqs = 50 WHERE COALESCE(max_target_seqs,5) = 5")
  }

  # --- schema v5: taxonomy-stage timestamp on samples ("Last Updated" column) ---
  add_col("samples", "time_stamp", "INTEGER")

  # --- schema v6: assembly_backup (original pre-trim state, for restore) ---
  if (!has_tbl("assembly_backup")) DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)

  if (ver < RIBOPILOT_SCHEMA_VERSION) {
    DBI::dbExecute(con, sprintf("PRAGMA user_version = %d", RIBOPILOT_SCHEMA_VERSION))
  }
  invisible(RIBOPILOT_SCHEMA_VERSION)
}

#' Copy a project's `.sqlite` into `<path>/.old_sqlite_dbs/` before an upgrade
#'
#' Migrations `ALTER` tables in place, which makes the DB unopenable by an older
#' RiboPilot the user may still have on a cluster; keep a copy to fall back to.
#' @noRd
backup_project_db <- function(path) {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) return(invisible(NULL))
  bdir <- file.path(path, ".old_sqlite_dbs")
  dir.create(bdir, showWarnings = FALSE)
  dest <- file.path(bdir, sprintf("sqlite.%d", as.integer(Sys.time())))
  tryCatch(file.copy(db, dest, overwrite = FALSE),
    error = function(e) warning("Could not back up project DB: ", conditionMessage(e)))
  invisible(dest)
}

#' Migrate the project database at a path to the current schema
#'
#' @param path Project directory containing `.sqlite`.
#' @return (invisibly) the schema version after migration.
#' @export
migrate_project <- function(path = ".") {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) {
    stop("No project database found at ", db, call. = FALSE)
  }
  # Back up first when an actual upgrade is pending (older schema version).
  ver <- tryCatch({
    c0 <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(c0), add = TRUE)
    DBI::dbGetQuery(c0, "PRAGMA user_version")[[1]]
  }, error = function(e) 0L)
  if (isTRUE(ver < RIBOPILOT_SCHEMA_VERSION)) backup_project_db(path)

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
}
