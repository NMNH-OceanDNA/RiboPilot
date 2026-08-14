# Cross-project store of fetched reference sequences and their locally derived
# element coordinates. Lives under the user-chosen reference folder so a second
# project in the same clade neither re-fetches nor re-annotates.

#' Path to the store inside a reference folder
#' @noRd
ref_store_path <- function(ref_dir) {
  file.path(path.expand(ref_dir), "refstore", "refstore.sqlite")
}

#' @noRd
REF_STORE_SEQ_DDL <- "CREATE TABLE IF NOT EXISTS ref_sequences (
  accession TEXT NOT NULL,
  species TEXT,
  title TEXT,
  sequence TEXT,
  length INTEGER,
  kingdom TEXT,
  taxon TEXT,
  tool TEXT,
  stored_utc TEXT,
  PRIMARY KEY (accession)
);"

#' @noRd
REF_STORE_ANN_DDL <- "CREATE TABLE IF NOT EXISTS ref_annotations (
  accession TEXT NOT NULL,
  element TEXT NOT NULL,
  type TEXT,
  start INTEGER,
  stop INTEGER,
  strand TEXT,
  PRIMARY KEY (accession, element)
);"

#' Open the reference store under a reference folder
#'
#' @param ref_dir Reference folder, as resolved by [ref_config()].
#' @param write TRUE to create the store if absent and open it writable.
#' @return A DBI connection, or NULL when the store does not exist and `write`
#'   is FALSE.
#' @noRd
open_ref_store <- function(ref_dir, write = FALSE) {
  p <- ref_store_path(ref_dir)
  if (!write) {
    if (!file.exists(p)) return(NULL)
    return(open_sqlite_ro(p))
  }
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "PRAGMA busy_timeout=10000;")
  DBI::dbExecute(con, sprintf("PRAGMA journal_mode=%s;", sqlite_journal_mode(p)))
  DBI::dbExecute(con, REF_STORE_SEQ_DDL)
  DBI::dbExecute(con, REF_STORE_ANN_DDL)
  con
}

#' Fetch stored records for a set of accessions
#' @return list(seqs, ann); zero-row frames when nothing matches.
#' @noRd
ref_store_get <- function(store, accessions) {
  empty <- list(seqs = data.frame(), ann = data.frame())
  if (is.null(store) || length(accessions) == 0) return(empty)
  ph <- paste(rep("?", length(accessions)), collapse = ",")
  seqs <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_sequences WHERE accession IN (%s)", ph),
    params = as.list(accessions))
  if (nrow(seqs) == 0) return(empty)
  ann <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_annotations WHERE accession IN (%s)", ph),
    params = as.list(accessions))
  list(seqs = seqs, ann = ann)
}

#' Write records into the store
#'
#' One transaction per accession, so a failure partway leaves earlier accessions
#' stored and a concurrent writer contends for a short lock rather than a long
#' one.
#' @return (invisibly) the number of accessions written.
#' @noRd
ref_store_put <- function(store, seqs, ann, kingdom, taxon, tool) {
  if (is.null(store) || is.null(seqs) || nrow(seqs) == 0) return(invisible(0L))
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  n <- 0L
  for (i in seq_len(nrow(seqs))) {
    acc <- seqs$accession[i]
    rows <- if (!is.null(ann) && nrow(ann) > 0) ann[ann$accession == acc, , drop = FALSE] else NULL
    DBI::dbWithTransaction(store, {
      DBI::dbExecute(store, "DELETE FROM ref_sequences WHERE accession = ?", params = list(acc))
      DBI::dbAppendTable(store, "ref_sequences", data.frame(
        accession = acc, species = seqs$species[i], title = seqs$title[i],
        sequence = seqs$sequence[i], length = seqs$length[i],
        kingdom = kingdom, taxon = taxon, tool = tool, stored_utc = ts,
        stringsAsFactors = FALSE))
      DBI::dbExecute(store, "DELETE FROM ref_annotations WHERE accession = ?", params = list(acc))
      if (!is.null(rows) && nrow(rows) > 0) {
        DBI::dbAppendTable(store, "ref_annotations", data.frame(
          accession = acc, element = rows$element, type = rows$type,
          start = rows$start, stop = rows$stop, strand = rows$strand,
          stringsAsFactors = FALSE))
      }
    })
    n <- n + 1L
  }
  invisible(n)
}
