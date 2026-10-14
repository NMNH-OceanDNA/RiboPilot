# Export review report: sequence gathering + per-feature alignment helpers, and
# the rmarkdown entry point used by export_bundle(). The heavy lifting lives
# here rather than in inst/export_report.Rmd so it stays unit-testable.

#' Every annotated element name across a set of exported units, in canonical order
#'
#' Elements outside `RRNA_ELEMENT_ORDER` are kept (sorted) at the end so a
#' project with custom element names still gets an alignment section.
#' @noRd
report_elements <- function(con, units) {
  if (is.null(units) || nrow(units) == 0) return(character(0))
  els <- unlist(lapply(seq_len(nrow(units)), function(i) {
    DBI::dbGetQuery(con,
      "SELECT DISTINCT element FROM review_annotations_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(units$ID[i], as.integer(units$path[i]),
        as.integer(units$scaffold[i])))$element
  }), use.names = FALSE)
  els <- unique(els[!is.na(els) & nzchar(els)])
  known <- RRNA_ELEMENT_ORDER[RRNA_ELEMENT_ORDER %in% els]
  c(known, sort(setdiff(els, known)))
}

#' Slice one element's sequence out of every unit that annotates it
#'
#' @param con Open project DB connection.
#' @param units data.frame(ID, path, scaffold, seqid), one row per exported unit.
#' @param element Element name, e.g. "18S".
#' @return A named `Biostrings::DNAStringSet` (names = seqid). Units lacking the
#'   feature, an assembly sequence, or a usable span are skipped.
#' @noRd
gather_feature_seqs <- function(con, units, element) {
  empty <- Biostrings::DNAStringSet(character(0))
  if (is.null(units) || nrow(units) == 0) return(empty)
  seqs <- character(0)
  nms <- character(0)
  for (i in seq_len(nrow(units))) {
    id <- units$ID[i]
    p <- as.integer(units$path[i])
    s <- as.integer(units$scaffold[i])
    ann <- DBI::dbGetQuery(con,
      "SELECT start, stop FROM review_annotations_v
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ? ORDER BY start",
      params = list(id, p, s, element))
    if (nrow(ann) == 0 || is.na(ann$start[1]) || is.na(ann$stop[1])) next
    if (ann$stop[1] < ann$start[1]) next
    seq <- DBI::dbGetQuery(con,
      "SELECT sequence FROM review_units_v WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, p, s))$sequence
    if (length(seq) == 0 || is.na(seq[1]) || !nzchar(seq[1])) next
    sub <- substr(seq[1], ann$start[1], ann$stop[1])
    if (!nzchar(sub)) next
    seqs <- c(seqs, sub)
    nms <- c(nms, if (!is.null(units$seqid)) units$seqid[i] else id)
  }
  if (length(seqs) == 0) return(empty)
  names(seqs) <- nms
  Biostrings::DNAStringSet(seqs)
}

#' Align one element's sequences across units and summarise them
#'
#' Rows come back ordered by complete-linkage clustering on the DECIPHER
#' distance matrix, so near-identical sequences sit next to each other in the
#' widget instead of in whatever order the units were exported.
#'
#' @param seqs A named `Biostrings::DNAStringSet`.
#' @return list(alignment, n, min_len, max_len, mean_identity), or NULL when
#'   fewer than two sequences were supplied or the alignment failed.
#' @noRd
feature_alignment <- function(seqs) {
  if (is.null(seqs) || length(seqs) < 2) return(NULL)
  aln <- tryCatch(DECIPHER::AlignSeqs(seqs, verbose = FALSE), error = function(e) NULL)
  if (is.null(aln) || length(aln) < 2) return(NULL)
  d <- tryCatch(DECIPHER::DistanceMatrix(aln, includeTerminalGaps = FALSE,
    verbose = FALSE), error = function(e) NULL)
  if (!is.null(d) && nrow(d) == length(aln)) {
    ord <- tryCatch(stats::hclust(stats::as.dist(d), method = "complete")$order,
      error = function(e) NULL)
    if (!is.null(ord)) aln <- aln[ord]
    pid <- 100 * (1 - mean(d[upper.tri(d)], na.rm = TRUE))
  } else {
    pid <- NA_real_
  }
  lens <- Biostrings::width(seqs)
  list(alignment = aln, n = length(aln), min_len = min(lens), max_len = max(lens),
    mean_identity = if (is.finite(pid)) pid else NA_real_)
}

#' Render the export review report
#'
#' Renders `inst/export_report.Rmd` to `out_file`, self-contained. Rendering
#' happens in a temporary directory so knitr's intermediates never land in the
#' export bundle. If rmarkdown is missing or the render fails, falls back to the
#' plain hand-built report so an export never ends up without one.
#'
#' @param con Open project DB connection.
#' @param units data.frame(ID, path, scaffold, seqid), one row per exported unit.
#' @param out_file Destination .html path.
#' @param label Export group / bundle label.
#' @param alignments Include the per-feature alignment section.
#' @return (invisibly) `out_file`.
#' @noRd
export_report <- function(con, units, out_file, label, alignments = TRUE) {
  fallback <- function(reason) {
    if (!is.null(reason)) message("Falling back to the plain export report: ", reason)
    export_report_html(con, units, out_file, label)
    invisible(out_file)
  }
  rmd <- system.file("export_report.Rmd", package = "RiboPilot")
  if (!requireNamespace("rmarkdown", quietly = TRUE) || !nzchar(rmd)) {
    return(fallback(if (nzchar(rmd)) "rmarkdown is not installed" else "report template not found"))
  }
  db_path <- tryCatch(DBI::dbGetInfo(con)$dbname, error = function(e) NULL)
  if (is.null(db_path) || !nzchar(db_path) || !file.exists(db_path)) {
    return(fallback("could not locate the project database"))
  }

  work <- file.path(tempdir(), paste0("ribopilot_report_", as.integer(Sys.time()),
    "_", sample.int(1e6, 1)))
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  units_rds <- file.path(work, "units.rds")
  saveRDS(units, units_rds)
  rmd_copy <- file.path(work, "export_report.Rmd")
  file.copy(rmd, rmd_copy, overwrite = TRUE)

  ok <- tryCatch({
    rmarkdown::render(rmd_copy,
      output_file = "report.html", output_dir = work, intermediates_dir = work,
      knit_root_dir = work, quiet = TRUE,
      params = list(db_path = normalizePath(db_path), label = label,
        units_rds = units_rds, alignments = isTRUE(alignments)),
      envir = new.env(parent = globalenv()))
    TRUE
  }, error = function(e) conditionMessage(e))
  built <- file.path(work, "report.html")
  if (!isTRUE(ok) || !file.exists(built)) {
    return(fallback(if (isTRUE(ok)) "the report template produced no output" else ok))
  }
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  file.copy(built, out_file, overwrite = TRUE)
  invisible(out_file)
}
