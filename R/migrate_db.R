# Backwards-compatible DB migration. Brings a project database created by an
# older RiboPilot up to the current schema. Migrations are additive and
# idempotent (guarded ALTER/CREATE), tracked by SQLite's `PRAGMA user_version`.
# new_db() stamps the current version; migrate_db() upgrades older projects.

#' Current project-DB schema version. Bump when adding a migration step below.
#' @noRd
RIBOPILOT_SCHEMA_VERSION <- 2L

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

  if (ver < RIBOPILOT_SCHEMA_VERSION) {
    DBI::dbExecute(con, sprintf("PRAGMA user_version = %d", RIBOPILOT_SCHEMA_VERSION))
  }
  invisible(RIBOPILOT_SCHEMA_VERSION)
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
  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
}
