# Project database versioning and upgrades.
#
# Any table or column change edits new_db() AND appends a step here. Step k
# takes a version k database to version k + 1. Released steps are never edited.
# Views (REVIEW_VIEWS_DDL) are CREATE VIEW IF NOT EXISTS: a step that changes one
# drops it first. While a view is broken every ALTER TABLE RENAME fails, so a
# step that rebuilds a table a view reads drops the views and reruns
# review_apply_schema() afterwards.

#' Upgrade steps, oldest first: each is list(note = "...", fn = function(con) ...)
#' @noRd
db_migrations <- function() {
  list(
    list(note = "per-tool rRNA caller flags (annotate_opts.barrnap_opts, pybarrnap_opts)",
      fn = function(con) {
        DBI::dbExecute(con, "ALTER TABLE annotate_opts ADD COLUMN barrnap_opts TEXT")
        DBI::dbExecute(con, "ALTER TABLE annotate_opts ADD COLUMN pybarrnap_opts TEXT")
      }),
    list(note = "sample metadata (GEOME/GBIF/NCBI) tables and columns", fn = meta_apply_schema),
    list(note = "when a unit was last exported (annotate.export_time_stamp)",
      fn = function(con) {
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN export_time_stamp INTEGER")
      }),
    list(note = "review copies for edited contigs; ITSx origin and scores",
      fn = function(con) {
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_origin TEXT")
        DBI::dbExecute(con, "ALTER TABLE annotate ADD COLUMN itsx_scores TEXT")
        review_apply_schema(con)
        review_upgrade(con)
        DBI::dbExecute(con, "DROP TABLE assembly_backup")
      }),
    list(note = "per-feature BLAST (feature_blast, feature_blast_hits) and review copies",
      fn = feature_blast_apply_schema),
    list(note = "GenBank BLAST check (annotate.genbank_check, review_units.genbank_check)",
      fn = function(con) {
        DBI::dbExecute(con, "DROP VIEW IF EXISTS review_units_v")
        # step 4 makes review_units from the current DDL, so it may have the column already
        for (t in c("annotate", "review_units")) {
          if (!"genbank_check" %in% DBI::dbListFields(con, t)) {
            DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN genbank_check TEXT", t))
          }
        }
        DBI::dbExecute(con, REVIEW_VIEWS_DDL[["review_units_v"]])
      })
  )
}

#' Schema version new_db() writes and backwards_compatibility() upgrades to
#' @noRd
current_schema_version <- function() 1L + length(db_migrations())

SCHEMA_HISTORY_DDL <- "CREATE TABLE IF NOT EXISTS schema_history (
  event TEXT NOT NULL,
  schema_version INTEGER NOT NULL,
  ribopilot_version TEXT,
  time_stamp INTEGER,
  note TEXT
);"

#' Append one row to schema_history
#' @noRd
record_schema_event <- function(con, event, version, note = NA_character_) {
  DBI::dbExecute(con, SCHEMA_HISTORY_DDL)
  DBI::dbExecute(con, "INSERT INTO schema_history VALUES (?, ?, ?, ?, ?)",
    params = list(event, as.integer(version), as.character(utils::packageVersion("RiboPilot")),
      as.integer(Sys.time()), note))
}

#' Schema version of a project database, or NA when it has no history
#' @noRd
project_schema_version <- function(db) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
  on.exit(DBI::dbDisconnect(con))
  if (!DBI::dbExistsTable(con, "schema_history")) return(NA_integer_)
  v <- DBI::dbGetQuery(con,
    "SELECT schema_version FROM schema_history ORDER BY rowid DESC LIMIT 1")$schema_version
  if (length(v)) as.integer(v) else NA_integer_
}

#' Stop when a database is unversioned or newer than this RiboPilot
#' @noRd
stop_unsupported_schema <- function(v) {
  if (is.na(v)) {
    stop("This project was made before RiboPilot 0.1.0 and cannot be updated. ",
      "Create it again with new_project().", call. = FALSE)
  }
  if (v > current_schema_version()) {
    stop("This project was made by a newer RiboPilot (database version ", v,
      "; this RiboPilot reads version ", current_schema_version(), "). ",
      "Install a newer RiboPilot to open it.", call. = FALSE)
  }
}

#' Update an old project to the installed RiboPilot
#'
#' Upgrades the project database to the current schema one step at a time and
#' points the project's `.config` at the matching container. Each change is
#' logged in the database's `schema_history` table. The database is copied to
#' `.sqlite.bak.<timestamp>` before any step runs, and `.config` to
#' `.config.bak.<timestamp>` before it is edited. Only the container line of
#' `.config` changes, and only for the official `macguigand/ribopilot` image.
#' Projects made before RiboPilot 0.1.0 are not supported.
#'
#' @param path Project directory (default: working directory).
#' @param dry_run Report what would change without changing anything.
#' @return Invisibly, TRUE when something changed (or would change).
#' @export
backwards_compatibility <- function(path = ".", dry_run = FALSE) {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) stop("No RiboPilot project (.sqlite) found in: ", path, call. = FALSE)
  from <- project_schema_version(db)
  stop_unsupported_schema(from)
  to <- current_schema_version()
  steps <- db_migrations()[seq_len(to - from) + from - 1L]
  cgap <- container_version_gap(path)
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")

  if (!length(steps) && is.null(cgap)) {
    message("Project is up to date (database version ", from, ").")
    return(invisible(FALSE))
  }
  if (isTRUE(dry_run)) {
    for (i in seq_along(steps)) message("Would upgrade to version ", from + i, ": ", steps[[i]]$note)
    if (!is.null(cgap)) message("Would change container ", cgap$configured, " to ", cgap$expected)
    return(invisible(TRUE))
  }

  con <- open_project_db(db)
  on.exit(DBI::dbDisconnect(con))
  if (length(steps)) {
    bak <- paste0(db, ".bak.", stamp)
    DBI::dbExecute(con, "VACUUM INTO ?", params = list(bak))
    for (i in seq_along(steps)) {
      tryCatch(DBI::dbWithTransaction(con, {
        steps[[i]]$fn(con)
        record_schema_event(con, "upgraded", from + i, steps[[i]]$note)
      }), error = function(e) {
        stop("Upgrade to database version ", from + i, " failed: ", conditionMessage(e),
          "\nEarlier steps were kept. The original database is at ", bak, call. = FALSE)
      })
      message("Upgraded to database version ", from + i, ": ", steps[[i]]$note)
    }
  }
  if (!is.null(cgap)) {
    cfg <- file.path(path, ".config")
    lines <- readLines(cfg, warn = FALSE)
    file.copy(cfg, paste0(cfg, ".bak.", stamp))
    i <- grep("^\\s*container\\s*=", lines, perl = TRUE)[1]
    lines[i] <- sub(cgap$configured, cgap$expected, lines[i], fixed = TRUE)
    writeLines(lines, cfg)
    note <- paste(cgap$configured, "to", cgap$expected)
    record_schema_event(con, "config", to, note)
    message("Container changed from ", note)
  }
  invisible(TRUE)
}
