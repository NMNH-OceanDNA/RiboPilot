# Named FASTA-header configurations for export, the token data they are resolved
# against, and the checks that run before anything is written.
#
# Mirrors MitoPilot's export_opts: a template is text with {column} tokens
# ("{seqid} [organism={Taxon}] ..."), stored under a name so a project can keep
# several and switch between them. Tokens resolve against the sample's own
# metadata, so any extra column of the user's mapping file is usable with no
# extra plumbing.

EXPORT_OPTS_DDL <- "CREATE TABLE IF NOT EXISTS export_opts (
  export_opts TEXT NOT NULL,
  fasta_header TEXT,
  fasta_header_marker TEXT,
  PRIMARY KEY (export_opts)
);"

#' Package default header templates
#' @noRd
DEFAULT_FASTA_HEADER <- "{seqid} [organism={organism}] [moltype=genomic DNA] {organism} {rrna_definition}"

#' @noRd
DEFAULT_FASTA_HEADER_MARKER <- "{seqid}_{element} [organism={organism}] {organism} {product}"

#' Create the export_opts table and seed its 'default' row (idempotent)
#' @noRd
ensure_export_opts <- function(con) {
  DBI::dbExecute(con, EXPORT_OPTS_DDL)
  DBI::dbExecute(con,
    "INSERT OR IGNORE INTO export_opts (export_opts, fasta_header, fasta_header_marker)
     VALUES ('default', ?, ?)",
    params = list(DEFAULT_FASTA_HEADER, DEFAULT_FASTA_HEADER_MARKER))
  invisible()
}

#' Names of the saved header configurations ('default' always first)
#' @noRd
list_export_templates <- function(con) {
  nm <- tryCatch(DBI::dbGetQuery(con, "SELECT export_opts FROM export_opts")$export_opts,
    error = function(e) character(0))
  unique(c("default", sort(nm[nzchar(nm) & !is.na(nm)])))
}

#' One saved configuration, falling back to the package defaults per field
#'
#' Defensive on purpose: a project made before this table existed, or a row with
#' a blanked-out box, must still export rather than fail.
#' @noRd
get_export_opts <- function(con, name = "default") {
  r <- tryCatch(DBI::dbGetQuery(con,
    "SELECT fasta_header, fasta_header_marker FROM export_opts WHERE export_opts = ?",
    params = list(name %||% "default")), error = function(e) data.frame())
  pick <- function(v, dflt) {
    if (length(v) == 0 || is.na(v[1]) || !nzchar(trimws(v[1]))) dflt else v[1]
  }
  list(
    fasta_header = pick(r$fasta_header, DEFAULT_FASTA_HEADER),
    fasta_header_marker = pick(r$fasta_header_marker, DEFAULT_FASTA_HEADER_MARKER))
}

#' Save (create or update) one named header configuration
#' @noRd
set_export_opts <- function(con, fasta_header, fasta_header_marker, name = "default") {
  tryCatch({
    ensure_export_opts(con)
    DBI::dbExecute(con,
      "INSERT INTO export_opts (export_opts, fasta_header, fasta_header_marker)
       VALUES (?, ?, ?)
       ON CONFLICT(export_opts) DO UPDATE SET
         fasta_header = excluded.fasta_header,
         fasta_header_marker = excluded.fasta_header_marker",
      params = list(name %||% "default", fasta_header, fasta_header_marker))
  }, error = function(e) warning("Could not save header template: ", conditionMessage(e), call. = FALSE))
  invisible()
}

#' GenBank definition phrase for one unit's annotated elements
#'
#' Reproduces the NCBI wording, e.g. "18S ribosomal RNA gene, partial sequence;
#' internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal
#' transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, partial
#' sequence." Consecutive elements sharing a completeness are phrased as one run,
#' which is what makes the sentence read as GenBank's does.
#'
#' @param elements Element names in the order they sit on the scaffold.
#' @param partial Logical of the same length: is that element partial?
#' @noRd
rrna_definition_text <- function(elements, partial) {
  if (length(elements) == 0) return("nuclear ribosomal RNA gene region")
  partial <- as.logical(partial)
  partial[is.na(partial)] <- FALSE
  nm <- vapply(elements, function(el) {
    if (el %in% RRNA_GENES) paste0(el, " ribosomal RNA gene") else rrna_element_product(el)
  }, character(1), USE.NAMES = FALSE)
  # split into consecutive runs of equal completeness
  run <- cumsum(c(TRUE, partial[-1] != partial[-length(partial)]))
  parts <- vapply(unique(run), function(r) {
    ix <- which(run == r)
    n <- nm[ix]
    joined <- if (length(n) == 1) n
      else if (length(n) == 2) paste(n, collapse = " and ")
      else paste0(paste(n[-length(n)], collapse = ", "), ", and ", n[length(n)])
    paste0(joined, if (partial[ix[1]]) ", partial sequence" else ", complete sequence")
  }, character(1))
  if (length(parts) > 1) parts[length(parts)] <- paste0("and ", parts[length(parts)])
  paste0(paste(parts, collapse = "; "), ".")
}

#' Organism name to export: the Taxonomy step's name when it matched the whole
#' Taxon label, or when the curator picked a species there; otherwise the Taxon
#' label itself, so labels like "Upeneus sp. A" are kept
#' @noRd
export_organism <- function(taxon, sci_name = NA, match_type = NA, rank = NA) {
  taxon <- as.character(taxon)
  norm <- function(x) tolower(trimws(gsub("[_[:space:]]+", " ", x)))
  whole <- vapply(taxon, function(t) !is.na(t) && norm(clean_taxon(t)) == norm(t),
    logical(1), USE.NAMES = FALSE)
  use <- !is.na(sci_name) & nzchar(sci_name) &
    ((match_type %in% "exact" & whole) |
      (match_type %in% "manual" & rank %in% c("species", "subspecies", "varietas", "forma")))
  ifelse(use, sci_name, taxon)
}

#' Token data for the header templates: one row per exported unit
#'
#' Everything a template may reference. `samples` is taken whole so that every
#' extra column of the user's mapping file is a usable token; the unit's own
#' facts and four computed fields are added on top.
#'
#' @param con Open project DB connection.
#' @param units data.frame(ID, path, scaffold, seqid).
#' @return data.frame, one row per unit, in the order given.
#' @noRd
export_header_data <- function(con, units) {
  rows <- lapply(seq_len(nrow(units)), function(i) {
    id <- units$ID[i]; p <- as.integer(units$path[i]); s <- as.integer(units$scaffold[i])
    smp <- DBI::dbGetQuery(con, "SELECT * FROM samples WHERE ID = ?", params = list(id))
    if (nrow(smp) == 0) smp <- data.frame(ID = id, stringsAsFactors = FALSE)
    smp <- smp[, !(names(smp) %in% c("R1", "R2")), drop = FALSE]
    asm <- DBI::dbGetQuery(con,
      "SELECT length, topology, sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    ann <- DBI::dbGetQuery(con,
      "SELECT n_features, missing, partial, reviewed, export_group
         FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))
    feat <- DBI::dbGetQuery(con,
      "SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      params = list(id, p, s))
    seq_len_bp <- if (nrow(asm)) asm$length[1] else NA_integer_
    unit_partial <- isTRUE(nrow(ann) > 0 && identical(ann$partial[1], "yes"))
    # A unit marked partial changes only the completeness wording.
    ends <- partial_ends(feat, seq_len_bp)
    part <- ends$p5 | ends$p3
    extra <- data.frame(
      path = p, scaffold = s,
      length = seq_len_bp,
      topology = if (nrow(asm)) asm$topology[1] else NA_character_,
      n_features = if (nrow(ann)) ann$n_features[1] else NA_integer_,
      missing = if (nrow(ann)) ann$missing[1] else NA_character_,
      partial = if (nrow(ann)) ann$partial[1] else NA_character_,
      reviewed = if (nrow(ann)) ann$reviewed[1] else NA_character_,
      export_group = if (nrow(ann)) ann$export_group[1] else NA_character_,
      seqid = units$seqid[i],
      organism = export_organism(smp$Taxon %||% NA, smp$sci_name %||% NA,
        smp$match_type %||% NA, smp$rank %||% NA),
      completeness = if (unit_partial || any(part)) "partial sequence" else "complete sequence",
      elements = paste(feat$element, collapse = "-"),
      rrna_definition = rrna_definition_text(feat$element, part),
      stringsAsFactors = FALSE)
    cbind(smp[, !(names(smp) %in% names(extra)), drop = FALSE], extra)
  })
  cols <- Reduce(intersect, lapply(rows, names))
  out <- do.call(rbind, lapply(rows, function(r) r[, cols, drop = FALSE]))
  # glue renders an NA field as the four characters "NA", so a sample with no
  # Taxon would ship "[organism=NA]", which reads like a real organism name.
  # Blank them instead; export_summary_df flags the unit as missing_organism.
  for (nm in names(out)) {
    if (is.character(out[[nm]])) out[[nm]][is.na(out[[nm]])] <- ""
  }
  out
}

#' First unmatched brace in a template, or NULL
#' @noRd
find_unmatched_brace <- function(x) {
  ch <- strsplit(x, "")[[1]]
  open <- integer(0)
  for (i in seq_along(ch)) {
    if (ch[i] == "{") open <- c(open, i)
    else if (ch[i] == "}") {
      if (length(open) == 0) return(list(brace = "}", pos = i))
      open <- open[-length(open)]
    }
  }
  if (length(open) > 0) list(brace = "{", pos = open[1]) else NULL
}

#' Check a header template before anything is written
#'
#' A bad template is not a cosmetic problem: a header without `{seqid}` produces
#' deflines that disagree with the `>Feature` lines of the feature table, and
#' NCBI rejects the submission. Returns a level rather than stopping, so the app
#' can show it live while the user is still typing.
#'
#' @param template The header template.
#' @param data Token data (from `export_header_data()`); the first row is used
#'   for a dry run. NULL skips the dry run.
#' @param require_seqid Demand `{seqid}` (an error when any sample exports more
#'   than one unit, a warning otherwise).
#' @return list(ok, level = "ok"/"warn"/"error", message)
#' @noRd
validate_fasta_header <- function(template, data = NULL, require_seqid = TRUE) {
  ok <- function(msg = "Template is valid.") list(ok = TRUE, level = "ok", message = msg)
  warn <- function(msg) list(ok = TRUE, level = "warn", message = msg)
  err <- function(msg) list(ok = FALSE, level = "error", message = msg)

  if (is.null(template) || length(template) == 0 || !nzchar(trimws(template))) {
    return(err("Header is empty."))
  }
  br <- find_unmatched_brace(template)
  if (!is.null(br)) {
    near <- substr(template, max(1, br$pos - 12), min(nchar(template), br$pos + 12))
    return(err(sprintf("Unmatched '%s' at position %d (near \"%s\").", br$brace, br$pos, near)))
  }
  if (!is.null(data) && nrow(data) > 0) {
    bad <- tryCatch({
      fill_header(template, data[1, , drop = FALSE])
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(bad)) return(err(bad))
  }
  if (grepl("\\[\\s*isolate\\s*=", template, ignore.case = TRUE) &&
      any(!is.na(data$seqid) & data$seqid != data$ID)) {
    return(err(paste("Take [isolate=...] out of the template for this group. RiboPilot sets",
      "[isolate=<ID>_p<path>_s<scaffold>] on each record of a sample that exports more than",
      "one scaffold, so all of them can go into one GenBank submission.")))
  }
  if (require_seqid && !grepl("{seqid}", template, fixed = TRUE)) {
    multi <- !is.null(data) && any(duplicated(data$ID))
    msg <- paste("No {seqid}. Records would not be named the way the feature table names them,",
      "and NCBI would reject the submission.")
    return(if (multi) err(msg) else warn(msg))
  }
  if (!is.null(data) && nrow(data) > 0) {
    em <- empty_modifiers(vapply(seq_len(nrow(data)), function(i) {
      fill_header(template, data[i, , drop = FALSE])
    }, character(1)))
    hit <- lengths(em) > 0
    if (any(hit)) {
      return(warn(sprintf(paste("Empty %s on %d of %d record(s) (%s). NCBI does not accept an",
        "empty modifier: fill in the column or take the modifier out."),
        paste0("[", unique(unlist(em)), "=]", collapse = ", "), sum(hit), length(hit),
        paste(utils::head((data$seqid %||% data$ID)[hit], 3), collapse = ", "))))
    }
  }
  ok()
}

#' Names of the source modifiers left empty (or "NA") in each header
#' @noRd
empty_modifiers <- function(headers) {
  m <- regmatches(headers, gregexpr("\\[([^]=[]+)=\\s*(NA)?\\s*\\]", headers))
  lapply(m, function(x) unique(sub("^\\[([^=]+)=.*$", "\\1", x)))
}

#' Fill a template's {column} tokens from one row of token data
#'
#' A plain lookup: nothing inside braces is ever run, because a template can
#' come from a project database someone else wrote.
#' @return the filled string; an error names the first token that is not a column.
#' @noRd
fill_header <- function(template, row) {
  p <- regmatches(template, gregexpr("\\{[^{}]*\\}", template), invert = NA)[[1]]
  tok <- seq_along(p) %% 2L == 0L
  nm <- substr(p[tok], 2L, nchar(p[tok]) - 1L)
  miss <- setdiff(nm, names(row))
  if (length(miss)) {
    stop(sprintf("Column \"%s\" is not in the database. Use one of the available columns listed below.",
      miss[1]), call. = FALSE)
  }
  p[tok] <- vapply(nm, function(n) as.character(row[[n]][1]), character(1))
  paste(p, collapse = "")
}

#' Resolve a template against token data
#'
#' A sample that exports more than one record gets `[isolate=<seqid>]` appended,
#' so all of its records can go into one GenBank submission. `validate_fasta_header()`
#' refuses a template's own `[isolate=]` for such samples.
#'
#' @param template Header template. @param data rows from `export_header_data()`.
#' @return character vector of headers, one per row, without the leading ">".
#' @noRd
render_fasta_headers <- function(template, data) {
  h <- vapply(seq_len(nrow(data)), function(i) {
    fill_header(template, data[i, , drop = FALSE])
  }, character(1))
  multi <- !is.na(data$seqid) & data$seqid != data$ID
  h[multi] <- paste0(h[multi], " [isolate=", data$seqid[multi], "]")
  h
}
