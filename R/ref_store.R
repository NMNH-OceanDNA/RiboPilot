# Cross-project store of fetched reference sequences and their locally derived
# element coordinates. Lives under the user-chosen reference folder so a second
# project in the same clade neither re-fetches nor re-annotates.

#' Path to the store inside a reference folder
#' @noRd
ref_store_path <- function(ref_dir) {
  file.path(path.expand(ref_dir), "refstore", "refstore.sqlite")
}

#' Whether ref_dir is usable as a real path
#'
#' NA or blank routes through file.path() into a bogus but literal path
#' ("NA/refstore/..." relative, or "/refstore/..." absolute for "") that can
#' silently create or read a store in the wrong place instead of failing
#' closed.
#' @noRd
ref_dir_ok <- function(ref_dir) {
  length(ref_dir) == 1 && !is.na(ref_dir) && nzchar(trimws(ref_dir))
}

#' Coordinates in ref_sequences are only comparable within one (kingdom, taxon)
#' combination, since a different barrnap kingdom or ITSx taxon can produce
#' different start/stop for the same accession. The primary key includes both
#' so the store can hold more than one combination per accession.
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
  ribopilot_version TEXT,
  stored_utc TEXT,
  PRIMARY KEY (accession, kingdom, taxon)
);"

#' `tool` here is the per-element caller (barrnap, ITSx, or derived), as
#' returned by annotate_records(); distinct from ref_sequences.tool, the
#' per-accession pipeline tag.
#' @noRd
REF_STORE_ANN_DDL <- "CREATE TABLE IF NOT EXISTS ref_annotations (
  accession TEXT NOT NULL,
  element TEXT NOT NULL,
  type TEXT,
  start INTEGER,
  stop INTEGER,
  strand TEXT,
  tool TEXT,
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
  # Always "delete", not sqlite_journal_mode(): this is a low-write cache with
  # short transactions, busy_timeout already covers contention, and "delete" is
  # the only mode that survives being read from a read-only directory (the
  # recommended HPC deployment) -- WAL needs to create a -shm file even to read.
  DBI::dbExecute(con, "PRAGMA journal_mode=delete;")
  DBI::dbExecute(con, REF_STORE_SEQ_DDL)
  DBI::dbExecute(con, REF_STORE_ANN_DDL)
  con
}

#' Fetch stored records for a set of accessions under one kingdom/taxon
#'
#' Coordinates are only comparable within the kingdom/taxon they were derived
#' under (see REF_STORE_SEQ_DDL), so a caller under different settings must not
#' match a row written under different ones.
#' @return list(seqs, ann); zero-row frames when nothing matches.
#' @noRd
ref_store_get <- function(store, accessions, kingdom, taxon) {
  empty <- list(seqs = data.frame(), ann = data.frame())
  if (is.null(store) || length(accessions) == 0) return(empty)
  ph <- paste(rep("?", length(accessions)), collapse = ",")
  seqs <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_sequences WHERE accession IN (%s) AND kingdom = ? AND taxon = ?", ph),
    params = c(as.list(accessions), list(kingdom, taxon)))
  if (nrow(seqs) == 0) return(empty)
  ph2 <- paste(rep("?", nrow(seqs)), collapse = ",")
  ann <- DBI::dbGetQuery(store,
    sprintf("SELECT * FROM ref_annotations WHERE accession IN (%s)", ph2),
    params = as.list(seqs$accession))
  list(seqs = seqs, ann = ann)
}

#' Write records into the store
#'
#' One transaction per accession, so a failure partway leaves earlier accessions
#' stored and a concurrent writer contends for a short lock rather than a long
#' one. A failure on one accession stops the loop there rather than throwing,
#' so the caller still gets the count of what was actually written.
#' @return (invisibly) the number of accessions written.
#' @noRd
ref_store_put <- function(store, seqs, ann, kingdom, taxon, tool) {
  if (is.null(store) || is.null(seqs) || nrow(seqs) == 0) return(invisible(0L))
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  version <- as.character(utils::packageVersion("RiboPilot"))
  n <- 0L
  for (i in seq_len(nrow(seqs))) {
    acc <- seqs$accession[i]
    rows <- if (!is.null(ann) && nrow(ann) > 0) ann[ann$accession == acc, , drop = FALSE] else NULL
    ok <- tryCatch({
      DBI::dbWithTransaction(store, {
        DBI::dbExecute(store, "DELETE FROM ref_sequences WHERE accession = ? AND kingdom = ? AND taxon = ?",
          params = list(acc, kingdom, taxon))
        DBI::dbAppendTable(store, "ref_sequences", data.frame(
          accession = acc, species = seqs$species[i], title = seqs$title[i],
          sequence = seqs$sequence[i], length = seqs$length[i],
          kingdom = kingdom, taxon = taxon, tool = tool, ribopilot_version = version,
          stored_utc = ts, stringsAsFactors = FALSE))
        DBI::dbExecute(store, "DELETE FROM ref_annotations WHERE accession = ?", params = list(acc))
        if (!is.null(rows) && nrow(rows) > 0) {
          DBI::dbAppendTable(store, "ref_annotations", data.frame(
            accession = acc, element = rows$element, type = rows$type,
            start = rows$start, stop = rows$stop, strand = rows$strand,
            tool = if ("tool" %in% names(rows)) rows$tool else NA_character_,
            stringsAsFactors = FALSE))
        }
      })
      TRUE
    }, error = function(e) FALSE)
    if (!ok) break
    n <- n + 1L
  }
  invisible(n)
}

#' Open the store and fetch records for the wanted accessions, as one guarded step
#'
#' `open_sqlite_ro()` only sets `busy_timeout`; it never reads a page, so a
#' file that exists but is not a valid SQLite database "opens" fine and only
#' fails once a query runs. Both the open and the read have to be inside the
#' same guard for a corrupt store to degrade cleanly.
#' @return list(seqs, ann), or NULL when the store does not exist.
#' @noRd
prime_read_store <- function(ref_dir, want, kingdom, taxon) {
  store <- open_ref_store(ref_dir, write = FALSE)
  if (is.null(store)) return(NULL)
  on.exit(DBI::dbDisconnect(store), add = TRUE)
  ref_store_get(store, want, kingdom, taxon)
}

#' Copy one accession's stored record into the project DB
#'
#' The pipeline runs one priming call per sample id and Nextflow parallelises
#' across ids against the same project database, so two processes racing to
#' insert the same accession is expected. `dbWithTransaction()` already rolls
#' back cleanly on a constraint violation; this only needs to keep that
#' failure from aborting the accessions still queued behind it.
#' @return TRUE if copied, FALSE if skipped because another process already
#'   holds this accession.
#' @noRd
prime_copy_accession <- function(con, acc, seq_row, ann_rows, ts) {
  tryCatch({
    DBI::dbWithTransaction(con, {
      DBI::dbAppendTable(con, "blast_ref_sequences", data.frame(
        accession = acc, species = seq_row$species, title = seq_row$title,
        taxid = NA_character_, lineage = NA_character_, sequence = seq_row$sequence,
        length = seq_row$length, time_stamp = ts, stringsAsFactors = FALSE))
      if (!is.null(ann_rows) && nrow(ann_rows) > 0) {
        DBI::dbAppendTable(con, "blast_ref_annotations", data.frame(
          accession = acc, element = ann_rows$element, type = ann_rows$type,
          start = ann_rows$start, stop = ann_rows$stop, strand = ann_rows$strand,
          time_stamp = ts, stringsAsFactors = FALSE))
      }
    })
    TRUE
  }, error = function(e) {
    if (!grepl("constraint failed", conditionMessage(e), ignore.case = TRUE)) stop(e)
    warning("skipping accession primed by a concurrent run: ", acc, call. = FALSE)
    FALSE
  })
}

#' Resolve the kingdom/taxon a sample's coordinates were (or would be) derived
#' under, the same join blast_ref.nf uses for its own annotate_records() call
#'
#' @param con Open project DB connection.
#' @param id Sample ID, or NULL to resolve from any row (used when priming is
#'   not scoped to one sample).
#' @return list(kingdom, taxon); c("euk", "all") when the query errors or
#'   returns nothing, matching the COALESCE defaults blast_ref.nf falls back to.
#' @noRd
prime_kingdom_taxon <- function(con, id = NULL) {
  sql <- "SELECT COALESCE(ao.rrna_kingdom,'euk'), COALESCE(ao.itsx_taxon,'all')
          FROM annotate an JOIN annotate_opts ao ON an.annotate_opts = ao.annotate_opts"
  r <- tryCatch({
    if (is.null(id)) {
      DBI::dbGetQuery(con, paste(sql, "LIMIT 1"))
    } else {
      DBI::dbGetQuery(con, paste(sql, "WHERE an.ID = ? LIMIT 1"), params = list(id))
    }
  }, error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(list(kingdom = "euk", taxon = "all"))
  list(kingdom = r[[1]][1], taxon = r[[2]][1])
}

#' Copy stored references for a project's candidates into the project database
#'
#' The fetch step selects only candidate accessions absent from
#' `blast_ref_sequences`, so anything copied here is not fetched again. An
#' accession the project already holds is left alone. Priming is grouped by
#' sample id because kingdom/taxon -- and so which store row is a valid match
#' -- can differ per id; with no ids given, everything is primed as one group
#' resolved from any row.
#'
#' @param con Open project DB connection.
#' @param ref_dir Reference folder, as resolved by [ref_config()].
#' @param ids Optional sample IDs to limit to; default is every candidate.
#' @return (invisibly) the number of accessions copied.
#' @export
prime_from_ref_store <- function(con, ref_dir, ids = NULL) {
  if (!ref_dir_ok(ref_dir)) return(invisible(0L))
  base_sql <- "SELECT DISTINCT c.accession FROM blast_ref_candidates c
          LEFT JOIN blast_ref_sequences s ON s.accession = c.accession
          WHERE s.accession IS NULL"
  want_for <- function(id) {
    if (is.null(id)) {
      DBI::dbGetQuery(con, base_sql)$accession
    } else {
      DBI::dbGetQuery(con, paste(base_sql, "AND c.ID = ?"), params = list(id))$accession
    }
  }

  groups <- if (is.null(ids)) list(NULL) else as.list(ids)
  n <- 0L
  for (id in groups) {
    want <- want_for(id)
    if (length(want) == 0) next
    kt <- prime_kingdom_taxon(con, id)

    got <- tryCatch(prime_read_store(ref_dir, want, kt$kingdom, kt$taxon), error = function(e) {
      warning("reference store at ", ref_store_path(ref_dir), " could not be read, skipping: ",
        conditionMessage(e), call. = FALSE)
      NULL
    })
    if (is.null(got) || nrow(got$seqs) == 0) next

    ts <- as.integer(Sys.time())
    for (i in seq_len(nrow(got$seqs))) {
      acc <- got$seqs$accession[i]
      rows <- if (nrow(got$ann) > 0) got$ann[got$ann$accession == acc, , drop = FALSE] else NULL
      if (prime_copy_accession(con, acc, got$seqs[i, , drop = FALSE], rows, ts)) n <- n + 1L
    }
  }
  invisible(n)
}

#' Write fetched records and their annotations into the shared store
#'
#' Best effort. A reference folder that cannot be written leaves the run
#' unaffected: the project keeps its own copy either way.
#'
#' @param ref_dir Reference folder.
#' @param fasta_file The multi-record fasta that was fetched.
#' @param ann annotate_records() output, or NULL.
#' @param kingdom,taxon Provenance for the coordinates.
#' @param tool Pipeline-level provenance label for this write, e.g.
#'   `"ribopilot-annotate"`. Not the name of one external binary: anchors come
#'   from whichever rRNA caller is configured (pybarrnap or barrnap) and
#'   spacers from ITSx, so a single per-accession string cannot name both.
#' @return (invisibly) the number of accessions written.
#' @export
push_to_ref_store <- function(ref_dir, fasta_file, ann, kingdom, taxon, tool) {
  if (!ref_dir_ok(ref_dir)) return(invisible(0L))
  store <- tryCatch(open_ref_store(ref_dir, write = TRUE), error = function(e) NULL)
  if (is.null(store)) return(invisible(0L))
  on.exit(DBI::dbDisconnect(store), add = TRUE)

  seqs <- tryCatch(Biostrings::readDNAStringSet(fasta_file), error = function(e) NULL)
  if (is.null(seqs) || length(seqs) == 0) return(invisible(0L))
  ids <- parse_defline_ids(names(seqs))
  seq_chars <- as.character(seqs)
  keep <- nzchar(seq_chars)
  if (!any(keep)) return(invisible(0L))
  sdf <- data.frame(accession = ids$acc[keep], species = extract_species(ids$title[keep]),
    title = ids$title[keep], sequence = seq_chars[keep],
    length = as.integer(Biostrings::width(seqs))[keep], stringsAsFactors = FALSE)

  adf <- NULL
  if (!is.null(ann) && nrow(ann) > 0) {
    adf <- data.frame(accession = ann$seqid, element = ann$element, type = ann$type,
      start = ann$start, stop = ann$stop, strand = ann$strand,
      tool = if ("tool" %in% names(ann)) ann$tool else NA_character_, stringsAsFactors = FALSE)
  }
  tryCatch(ref_store_put(store, sdf, adf, kingdom, taxon, tool),
    error = function(e) invisible(0L))
}
