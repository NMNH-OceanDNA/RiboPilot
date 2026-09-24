# Pure-R NCBI taxonomy backend.
#
# Taxonomy lookups used to shell out to `taxonkit`, an external binary users had
# to install themselves. On a shared cluster the R session almost never has it
# on PATH, so taxonomy assignment failed outright. The project's pinned taxdump
# is now parsed once into a cached SQLite index inside the project
# (`<project>/.taxonomy/tax_index.sqlite`) and every lookup is a query against
# that index, so nothing outside the R package is required.
#
# The index is built from the same pinned tarball as before, so a project still
# resolves against one frozen taxonomy version. Building it also removes the
# repeated ~500 MB extraction of the dmp files that every session used to do.

#' Schema version of the cached taxonomy index; bump to force a rebuild.
#' @noRd
TAX_INDEX_VERSION <- 2L

#' names.dmp classes worth indexing
#'
#' "authority" and "type material" are most of the file and are never what a
#' user (or a BLAST defline) hands us.
#' @noRd
TAX_NAME_CLASSES <- c("scientific name", "synonym", "equivalent name",
  "genbank synonym", "common name", "genbank common name")

#' Path to a taxdump's cached index
#' @noRd
tax_index_path <- function(taxdump_dir) file.path(taxdump_dir, "tax_index.sqlite")

#' Does `idx` exist at the current schema version?
#' @noRd
tax_index_ok <- function(idx) {
  if (!file.exists(idx)) {
    return(FALSE)
  }
  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), idx), error = function(e) NULL)
  if (is.null(con)) {
    return(FALSE)
  }
  on.exit(DBI::dbDisconnect(con))
  v <- tryCatch(DBI::dbGetQuery(con, "SELECT value FROM tax_meta WHERE key = 'version'")$value,
    error = function(e) NULL)
  length(v) == 1 && identical(suppressWarnings(as.integer(v)), TAX_INDEX_VERSION)
}

#' Stream a dmp file in chunks, handing selected columns to `fn`
#'
#' dmp fields are separated by "\\t|\\t" and lines end with "\\t|", so splitting
#' on tabs puts the payload in every other column and the callers pick those
#' positions. Reading in chunks keeps peak memory low: names.dmp is ~4 million
#' rows and a login node is a bad place to hold it all at once.
#' @noRd
read_dmp_chunked <- function(file, cols, col_names, fn, chunk = 500000L) {
  cb <- readr::SideEffectChunkCallback$new(function(x, pos) {
    if (ncol(x) < max(cols)) {
      return(invisible(NULL))
    }
    d <- as.data.frame(x[, cols, drop = FALSE])
    names(d) <- col_names
    fn(d)
    invisible(NULL)
  })
  readr::read_delim_chunked(file, cb, chunk_size = chunk, delim = "\t",
    quote = "", escape_double = FALSE, escape_backslash = FALSE,
    col_names = FALSE, col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE)
  invisible(NULL)
}

#' Build (once) the cached taxonomy index for a taxdump
#'
#' @param taxdump_dir Directory with a pinned taxdump (tarball or extracted dmps).
#' @param rebuild Rebuild even if a current index exists.
#' @return Path to the SQLite index, or `NULL` if there is no taxdump to index.
#' @noRd
build_tax_index <- function(taxdump_dir, rebuild = FALSE) {
  if (is.null(taxdump_dir) || !have_taxdump(taxdump_dir)) {
    return(NULL)
  }
  taxdump_dir <- path.expand(taxdump_dir)
  idx <- tax_index_path(taxdump_dir)
  if (!rebuild && tax_index_ok(idx)) {
    return(idx)
  }
  if (file.exists(file.path(taxdump_dir, "names.dmp"))) {
    data_dir <- taxdump_dir
  } else {
    data_dir <- extract_taxdump(taxdump_dir)
    on.exit(unlink(data_dir, recursive = TRUE), add = TRUE)
  }
  message("Indexing the NCBI taxonomy for this project (one time, a few minutes) ...")
  tmp <- paste0(idx, ".tmp")
  unlink(tmp)
  con <- DBI::dbConnect(RSQLite::SQLite(), tmp)
  ok <- tryCatch({
    DBI::dbExecute(con, "PRAGMA journal_mode = OFF")
    DBI::dbExecute(con, "PRAGMA synchronous = OFF")
    DBI::dbExecute(con, "CREATE TABLE tax_nodes (taxid INTEGER PRIMARY KEY, parent INTEGER, rank TEXT)")
    DBI::dbExecute(con, "CREATE TABLE tax_names (taxid INTEGER, name TEXT, sci INTEGER)")
    DBI::dbExecute(con, "CREATE TABLE tax_merged (old INTEGER PRIMARY KEY, new INTEGER)")

    read_dmp_chunked(file.path(data_dir, "nodes.dmp"), c(1, 3, 5),
      c("taxid", "parent", "rank"), function(d) {
        d$taxid <- suppressWarnings(as.integer(d$taxid))
        d$parent <- suppressWarnings(as.integer(d$parent))
        d <- d[!is.na(d$taxid), , drop = FALSE]
        if (nrow(d) > 0) DBI::dbAppendTable(con, "tax_nodes", d)
      })

    read_dmp_chunked(file.path(data_dir, "names.dmp"), c(1, 3, 7),
      c("taxid", "name", "class"), function(d) {
        d <- d[!is.na(d$class) & d$class %in% TAX_NAME_CLASSES, , drop = FALSE]
        if (nrow(d) == 0) {
          return(invisible(NULL))
        }
        out <- data.frame(
          taxid = suppressWarnings(as.integer(d$taxid)),
          name = trimws(d$name),
          sci = as.integer(d$class == "scientific name"),
          stringsAsFactors = FALSE)
        out <- out[!is.na(out$taxid) & nzchar(out$name), , drop = FALSE]
        if (nrow(out) > 0) DBI::dbAppendTable(con, "tax_names", out)
      })

    merged <- file.path(data_dir, "merged.dmp")
    if (file.exists(merged)) {
      read_dmp_chunked(merged, c(1, 3), c("old", "new"), function(d) {
        d$old <- suppressWarnings(as.integer(d$old))
        d$new <- suppressWarnings(as.integer(d$new))
        d <- d[!is.na(d$old) & !is.na(d$new), , drop = FALSE]
        if (nrow(d) > 0) DBI::dbAppendTable(con, "tax_merged", d)
      })
    }

    n_nodes <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM tax_nodes")$n
    n_names <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM tax_names")$n
    if (n_nodes < getOption("ribopilot.tax_min_nodes", 2e6) || n_names == 0) {
      stop("only ", n_nodes, " taxa were read; the taxdump is incomplete", call. = FALSE)
    }
    DBI::dbExecute(con, "CREATE INDEX i_names_name ON tax_names(name COLLATE NOCASE)")
    DBI::dbExecute(con, "CREATE INDEX i_names_taxid ON tax_names(taxid)")
    DBI::dbExecute(con, "CREATE INDEX i_nodes_parent ON tax_nodes(parent)")
    DBI::dbExecute(con, "CREATE TABLE tax_meta (key TEXT, value TEXT)")
    DBI::dbExecute(con, "INSERT INTO tax_meta VALUES ('version', ?)",
      params = list(as.character(TAX_INDEX_VERSION)))
    TRUE
  }, error = function(e) {
    warning("Could not index the taxonomy: ", conditionMessage(e), call. = FALSE)
    FALSE
  })
  DBI::dbDisconnect(con)
  if (!ok) {
    unlink(tmp)
    return(NULL)
  }
  file.rename(tmp, idx)
  idx
}

#' Open the cached index for a taxdump, building it if needed
#' @return An open connection, or `NULL` if there is no usable taxdump.
#' @noRd
tax_index_con <- function(taxdump_dir) {
  idx <- build_tax_index(taxdump_dir)
  if (is.null(idx)) {
    return(NULL)
  }
  DBI::dbConnect(RSQLite::SQLite(), idx)
}

#' Run `fn` against an open taxonomy index, or return `default`
#'
#' One connection per call site rather than per lookup: name resolution makes
#' several queries in a row and reopening the index for each is wasteful.
#' @noRd
with_tax_index <- function(taxdump_dir, fn, default = NULL) {
  con <- tryCatch(tax_index_con(taxdump_dir), error = function(e) NULL)
  if (is.null(con)) {
    return(default)
  }
  on.exit(DBI::dbDisconnect(con))
  fn(con)
}

#' Run an "IN (...)" query in batches, under SQLite's bound-parameter limit
#' @noRd
tax_in_query <- function(con, sql, ids, chunk = 400L) {
  ids <- unique(ids)
  if (length(ids) == 0) {
    return(NULL)
  }
  parts <- split(ids, ceiling(seq_along(ids) / chunk))
  out <- lapply(parts, function(p) {
    ph <- paste(rep("?", length(p)), collapse = ",")
    DBI::dbGetQuery(con, sprintf(sql, ph), params = as.list(p))
  })
  do.call(rbind, out)
}

#' Map retired taxids onto their replacements (merged.dmp)
#' @noRd
td_resolve_merged <- function(con, taxids) {
  m <- tax_in_query(con, "SELECT old, new FROM tax_merged WHERE old IN (%s)", taxids)
  if (is.null(m) || nrow(m) == 0) {
    return(taxids)
  }
  i <- match(taxids, m$old)
  ifelse(is.na(i), taxids, m$new[i])
}

#' Ancestor path (root -> tip) for each requested taxid
#'
#' @return Named list keyed by the REQUESTED taxid (as character), each a
#'   data.frame(taxid, rank, name) ordered root first. Taxids absent from the
#'   taxdump are simply missing from the list.
#' @noRd
td_paths <- function(con, taxids) {
  taxids <- unique(stats::na.omit(suppressWarnings(as.integer(taxids))))
  if (length(taxids) == 0) {
    return(list())
  }
  live <- td_resolve_merged(con, taxids)
  sql <- "WITH RECURSIVE up(seed, taxid, parent, depth) AS (
            SELECT taxid, taxid, parent, 0 FROM tax_nodes WHERE taxid IN (%s)
            UNION ALL
            SELECT up.seed, n.taxid, n.parent, up.depth + 1
              FROM up JOIN tax_nodes n ON n.taxid = up.parent
             WHERE up.taxid <> up.parent
          )
          SELECT up.seed AS seed, up.taxid AS taxid, up.depth AS depth, n.rank AS rank,
                 (SELECT name FROM tax_names WHERE taxid = up.taxid AND sci = 1 LIMIT 1) AS name
            FROM up JOIN tax_nodes n ON n.taxid = up.taxid
           ORDER BY up.seed, up.depth DESC"
  rows <- tax_in_query(con, sql, live)
  if (is.null(rows) || nrow(rows) == 0) {
    return(list())
  }
  by_seed <- split(rows, rows$seed)
  out <- list()
  for (k in seq_along(taxids)) {
    d <- by_seed[[as.character(live[k])]]
    if (is.null(d) || nrow(d) == 0) next
    d <- d[order(-d$depth), c("taxid", "rank", "name")]
    rownames(d) <- NULL
    out[[as.character(taxids[k])]] <- d
  }
  out
}

#' Lay an ancestor path onto the seven rank-aligned slots
#'
#' Missing ranks are filled the way `taxonkit reformat -F -t` filled them: an
#' "unclassified <nearest known ancestor> <rank>" placeholder name and an EMPTY
#' taxid, so the seven slots always line up positionally and a filled rank stays
#' distinguishable from a real one ([parse_rank_names()] reads the placeholder
#' back as NA).
#' @param path data.frame(taxid, rank, name), root first.
#' @return list(names, taxids), each a ";"-joined 7-field string.
#' @noRd
td_rank_align <- function(path) {
  n <- length(RRNA_LINEAGE_RANKS)
  if (is.null(path) || nrow(path) == 0) {
    return(list(names = NA_character_, taxids = NA_character_))
  }
  i <- match(RRNA_LINEAGE_RANKS, path$rank)
  nm <- ifelse(is.na(i), NA_character_, path$name[i])
  tx <- ifelse(is.na(i), NA_integer_, path$taxid[i])
  known <- path$name[!is.na(path$name)]
  anchor <- if (length(known) > 0) known[length(known)] else "root"
  for (k in seq_len(n)) {
    if (is.na(nm[k]) || !nzchar(nm[k])) {
      nm[k] <- paste("unclassified", anchor, RRNA_LINEAGE_RANKS[k])
    } else {
      anchor <- nm[k]
    }
  }
  list(
    names = paste(nm, collapse = ";"),
    taxids = paste(ifelse(is.na(tx), "", as.character(tx)), collapse = ";")
  )
}

#' Exact (case-insensitive) name -> taxid lookup, synonyms included
#'
#' When a name is some taxon's scientific name, only those taxa are returned.
#' @return Named list: input name -> integer vector of taxids (ascending).
#' @noRd
td_name2taxid <- function(con, names) {
  names <- unique(names)
  out <- stats::setNames(vector("list", length(names)), names)
  if (length(names) == 0) {
    return(out)
  }
  rows <- tax_in_query(con,
    "SELECT n.name AS name, n.taxid AS taxid, MAX(n.sci) AS sci FROM tax_names n
      WHERE n.name COLLATE NOCASE IN (%s) GROUP BY n.name, n.taxid ORDER BY n.taxid", names)
  if (is.null(rows) || nrow(rows) == 0) {
    return(out)
  }
  low <- tolower(rows$name)
  for (nm in names) {
    r <- rows[low == tolower(nm), , drop = FALSE]
    if (nrow(r) == 0) next
    if (any(r$sci == 1)) r <- r[r$sci == 1, , drop = FALSE]
    out[[nm]] <- unique(as.integer(r$taxid))
  }
  out
}

#' Species-rank descendants of a taxid (for fuzzy epithet matching)
#' @return data.frame(taxid, name)
#' @noRd
td_descendant_species <- function(con, taxid) {
  empty <- data.frame(taxid = integer(0), name = character(0), stringsAsFactors = FALSE)
  taxid <- suppressWarnings(as.integer(taxid))
  if (is.na(taxid)) {
    return(empty)
  }
  rows <- DBI::dbGetQuery(con,
    "WITH RECURSIVE down(taxid) AS (
       SELECT ?
       UNION ALL
       SELECT n.taxid FROM tax_nodes n JOIN down d ON n.parent = d.taxid
        WHERE n.taxid <> n.parent
     )
     SELECT t.taxid AS taxid,
            (SELECT name FROM tax_names WHERE taxid = t.taxid AND sci = 1 LIMIT 1) AS name
       FROM tax_nodes t JOIN down ON down.taxid = t.taxid
      WHERE t.rank = 'species'",
    params = list(taxid))
  if (nrow(rows) == 0) {
    return(empty)
  }
  rows <- rows[!is.na(rows$name), , drop = FALSE]
  data.frame(taxid = as.integer(rows$taxid), name = as.character(rows$name),
    stringsAsFactors = FALSE)
}
