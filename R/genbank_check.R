# GenBank BLAST check: the two BLAST flags NCBI's rRNA/ITS curation raises
# before assigning accessions (top-hit query coverage under 75%, more than one
# HSP on the top hit), run the way NCBI runs them (default megablast against
# ITS_eukaryote_sequences and the SSU plus LSU eukaryote sets) on the sequence
# export writes. Stored per contig as annotate.genbank_check or, on a review
# copy, review_units.genbank_check, like the stock ITSx check.

GENBANK_CHECK_DBS <- c(ITS = "ITS_eukaryote_sequences", SSU = "SSU_eukaryote_rRNA", LSU = "LSU_eukaryote_rRNA")
GENBANK_CHECK_MIN_QCOVS <- 75
GENBANK_CHECK_MAX_HSPS <- 1L
GENBANK_CHECK_EVALUE <- 0.05
GENBANK_CHECK_MAX_TARGET_SEQS <- 100L

#' Top hit of one search: accession, species, qcovs, and its HSP spans in query order
#' @noRd
genbank_top_hit <- function(hits) {
  if (!NROW(hits)) return(NULL)
  h <- hits[hits$accession == hits$accession[1], , drop = FALSE]
  h <- h[order(h$qstart), , drop = FALSE]
  list(accession = h$accession[1], species = extract_species(h$title[1]),
    qcovs = round(h$qcovs[1]), spans = paste0(h$qstart, "-", h$qend))
}

#' The stored text for one contig from its two searches
#'
#' `its_hits` and `rrna_hits` are parse_blast_hits_file() frames for one query
#' (possibly empty). The rRNA search's coverage is reported, never flagged.
#' @return character(1): "" when clean, else the flags and a top-hits note.
#' @noRd
genbank_check_calls <- function(its_hits, rrna_hits) {
  its <- genbank_top_hit(its_hits)
  rr <- genbank_top_hit(rrna_hits)
  multi <- function(top, db) if (length(top$spans) > GENBANK_CHECK_MAX_HSPS)
    sprintf("genbank_multi_hsp(%s %d HSPs: %s)", db, length(top$spans), paste(top$spans, collapse = ", "))
  flags <- c(
    if (is.null(its)) "genbank_no_hit(ITS)",
    if (!is.null(its) && its$qcovs < GENBANK_CHECK_MIN_QCOVS) sprintf("genbank_low_coverage(ITS %d%%)", its$qcovs),
    if (!is.null(its)) multi(its, "ITS"),
    if (is.null(rr)) "genbank_no_hit(rRNA)",
    if (!is.null(rr)) multi(rr, "rRNA"))
  if (!length(flags)) return("")
  if (is.null(its) && is.null(rr)) return(paste(flags, collapse = "; "))
  note <- function(top, db) if (is.null(top)) paste(db, "none") else
    sprintf("%s %s %s %d%%", db, top$accession, top$species, top$qcovs)
  paste0(paste(flags, collapse = "; "), " top hits: ", note(its, "ITS"), ", ", note(rr, "rRNA"))
}

#' Short form of a stored check for the export window
#' @noRd
genbank_check_display <- function(text) {
  vapply(text, function(x) {
    if (is.na(x) || !nzchar(x) || startsWith(x, "not_run(")) return(x %|NA|% "")
    x <- sub(" top hits: .*$", "", x)
    tok <- strsplit(x, "; ", fixed = TRUE)[[1]]
    tok <- sub("^genbank_low_coverage\\((.*)\\)$", "low coverage \\1", tok)
    tok <- sub("^genbank_multi_hsp\\((\\S+ \\d+ HSPs):.*\\)$", "\\1", tok)
    tok <- sub("^genbank_no_hit\\((.*)\\)$", "no \\1 hit", tok)
    paste(tok, collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

#' Directory of the NCBI databases from an rp_dbarg.R panel string, or NULL
#' @noRd
genbank_check_db_dir <- function(dbarg) {
  if (is.null(dbarg) || length(dbarg) != 1 || is.na(dbarg) || !nzchar(trimws(dbarg))) return(NULL)
  dirname(strsplit(trimws(dbarg), "[[:space:]]+")[[1]][1])
}

#' @noRd
genbank_blast_run <- function(query, dbarg, out, threads = 1L, blastn = "blastn") {
  status <- system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query),
    "-outfmt", shQuote(RECHECK_BLAST$outfmt), "-evalue", GENBANK_CHECK_EVALUE,
    "-max_target_seqs", GENBANK_CHECK_MAX_TARGET_SEQS, "-num_threads", threads),
    stdout = out, stderr = FALSE)
  if (!identical(as.integer(status), 0L)) stop(sprintf("blastn exited with status %s", status))
  parse_blast_hits_file(out)
}

#' Store a contig's GenBank check text on the annotate row or its review copy
#' @noRd
genbank_store <- function(con, id, path, scaffold, text, review = FALSE) {
  sql <- sprintf("UPDATE %s SET genbank_check = ? WHERE ID = ? AND path = ? AND scaffold = ?",
    if (review) "review_units" else "annotate")
  tryCatch(DBI::dbExecute(con, sql, params = list(text, id, as.integer(path), as.integer(scaffold))),
    error = function(e) if (review) stop(e))
  invisible(text)
}

#' Run the GenBank BLAST check on one contig and store the result
#'
#' @param db_dir Directory holding the three NCBI databases (the panel's
#'   directory). NULL stores `not_run(no local BLAST)`.
#' @param review TRUE checks the contig's review copy (what export writes) and
#'   stores on it.
#' @return (invisibly) the stored text.
#' @noRd
genbank_check_unit <- function(con, id, path = 1L, scaffold = 1L, db_dir = NULL,
                               workdir = tempfile("genbank_check_"), threads = 1L,
                               review = FALSE, blastn = "blastn") {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  seq <- DBI::dbGetQuery(con, sprintf("SELECT sequence FROM %s WHERE ID = ? AND path = ? AND scaffold = ?",
    if (review) "review_units_v" else "assemblies"), params = list(id, path, scaffold))$sequence[1]
  text <- if (is.na(seq %||% NA) || !nzchar(seq)) {
    "not_run(no sequence)"
  } else if (is.null(db_dir) || !nzchar(db_dir %|NA|% "")) {
    "not_run(no local BLAST)"
  } else {
    pre <- stats::setNames(file.path(db_dir, GENBANK_CHECK_DBS), names(GENBANK_CHECK_DBS))
    missing <- GENBANK_CHECK_DBS[!vapply(pre, blast_db_exists, logical(1))]
    if (length(missing)) {
      sprintf("not_run(no %s database)", missing[1])
    } else {
      dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
      q <- file.path(workdir, "query.fa")
      writeLines(c(">query", seq), q)
      tryCatch({
        its <- genbank_blast_run(q, pre[["ITS"]], file.path(workdir, "its.txt"), threads, blastn)
        rr <- genbank_blast_run(q, paste(pre[c("SSU", "LSU")], collapse = " "),
          file.path(workdir, "rrna.txt"), threads, blastn)
        genbank_check_calls(its, rr)
      }, error = function(e) paste0("not_run(BLAST failed: ", conditionMessage(e), ")"))
    }
  }
  genbank_store(con, id, path, scaffold, text, review = review)
}
