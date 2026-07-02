# GenBank submission bundle export (mirrors MitoPilot's group-based export).
#
# Samples are organized into export groups (samples.export_group). Each group /
# selection is exported as a combined <label>.fasta (source-modifier deflines) +
# <label>.tbl (concatenated 5-column feature tables), per-sample .fsa/.tbl,
# per-marker fastas (18S/28S/ITS...), a per-sample summary CSV, and a
# self-contained HTML report. table2asn is run if available.

#' Assign selected samples to an export group
#' @param path Project dir. @param ids Sample IDs. @param group Group name.
#' @export
assign_export_group <- function(path, ids, group) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  for (id in ids) DBI::dbExecute(con, "UPDATE samples SET export_group = ? WHERE ID = ?", params = list(group, id))
  invisible()
}

#' Clear the export group of selected samples
#' @export
clear_export_group <- function(path, ids) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  for (id in ids) DBI::dbExecute(con, "UPDATE samples SET export_group = NULL WHERE ID = ?", params = list(id))
  invisible()
}

#' Feature table (.tbl) text for one sample from an open connection
#' @noRd
feature_table_text_con <- function(con, id) {
  d <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
  hdr <- paste0(">Feature ", id)
  if (nrow(d) == 0) {
    return(paste0(hdr, "\n"))
  }
  seqlen <- DBI::dbGetQuery(con, "SELECT length FROM stitched WHERE ID = ?", params = list(id))$length
  seqlen <- if (length(seqlen) == 0 || is.na(seqlen)) max(d$stop) else seqlen
  lines <- hdr
  for (i in seq_len(nrow(d))) {
    s <- if (i == 1 && d$start[i] <= 1) paste0("<", d$start[i]) else as.character(d$start[i])
    e <- if (i == nrow(d) && d$stop[i] >= seqlen) paste0(">", d$stop[i]) else as.character(d$stop[i])
    lines <- c(lines, paste0(s, "\t", e, "\t", d$type[i]))
    lines <- c(lines, paste0("\t\t\tproduct\t", rrna_element_product(d$element[i])))
    if (d$type[i] == "misc_RNA") lines <- c(lines, paste0("\t\t\tnote\t", d$element[i]))
  }
  paste(lines, collapse = "\n")
}

#' Per-sample export summary data.frame
#' @noRd
export_summary_df <- function(con, ids) {
  DBI::dbGetQuery(con, sprintf(
    "SELECT s.ID, s.Taxon AS organism, s.export_group, st.reference_accession AS reference,
       st.length, an.n_features, an.missing, an.reviewed
     FROM samples s LEFT JOIN stitched st ON s.ID = st.ID LEFT JOIN annotate an ON s.ID = an.ID
     WHERE s.ID IN (%s)", paste(sprintf("'%s'", ids), collapse = ",")))
}

#' Inline SVG feature-map for one sample (for the HTML report)
#' @noRd
feature_map_svg <- function(con, id) {
  d <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
  if (nrow(d) == 0) {
    return("<em>no annotations</em>")
  }
  d$molecule <- id
  p <- ggplot2::ggplot(d, ggplot2::aes(xmin = start, xmax = stop, y = molecule, fill = element, label = element)) +
    gggenes::geom_gene_arrow(arrowhead_height = grid::unit(5, "mm"), arrow_body_height = grid::unit(4.5, "mm")) +
    gggenes::geom_gene_label(align = "centre", grow = FALSE) +
    ggplot2::scale_fill_manual(values = RRNA_ELEMENT_COLORS, guide = "none") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank())
  tf <- tempfile(fileext = ".svg")
  ggplot2::ggsave(tf, p, device = svglite::svglite, width = 9, height = 1.1)
  svg <- paste(readLines(tf, warn = FALSE), collapse = "\n")
  unlink(tf)
  sub("(?s).*?(<svg)", "\\1", svg, perl = TRUE) # strip XML/doctype header
}

#' Write a self-contained HTML report for an export group / selection
#' @noRd
export_report_html <- function(con, ids, out_html, label) {
  smry <- export_summary_df(con, ids)
  hdr <- c("ID", "organism", "reference", "length", "n_features", "missing", "reviewed")
  rows <- apply(smry[, intersect(hdr, names(smry)), drop = FALSE], 1, function(r) {
    paste0("<tr>", paste0("<td>", ifelse(is.na(r), "", r), "</td>", collapse = ""), "</tr>")
  })
  maps <- vapply(ids, function(id) {
    paste0("<div class='sample'><h3>", id, " <small>",
      smry$organism[smry$ID == id][1] %||% "", "</small></h3>",
      feature_map_svg(con, id), "</div>")
  }, character(1))
  html <- paste0(
    "<!doctype html><html><head><meta charset='utf-8'><title>RiboPilot export: ", label,
    "</title><style>body{font-family:system-ui,Arial,sans-serif;margin:2em;color:#222;}",
    "table{border-collapse:collapse;margin:1em 0;}th,td{border:1px solid #ddd;padding:4px 8px;font-size:0.9em;}",
    "th{background:#f4f6f8;text-align:left;}.sample{margin:1.2em 0;border-top:1px solid #eee;padding-top:0.5em;}",
    "small{color:#888;font-weight:normal;}svg{max-width:100%;height:auto;}</style></head><body>",
    "<h1>RiboPilot export report</h1><p><strong>Group:</strong> ", label,
    " &middot; <strong>Samples:</strong> ", length(ids),
    " &middot; ", format(Sys.time(), "%Y-%m-%d %H:%M"), "</p>",
    "<h2>Summary</h2><table><thead><tr>",
    paste0("<th>", intersect(hdr, names(smry)), "</th>", collapse = ""), "</tr></thead><tbody>",
    paste(rows, collapse = ""), "</tbody></table>",
    "<h2>Feature models</h2>", paste(maps, collapse = ""),
    "</body></html>")
  writeLines(html, out_html)
}

#' Build a GenBank submission bundle for a set of samples (one export group)
#'
#' @param path Project directory.
#' @param ids Sample IDs to export (default: all with a stitched sequence). If
#'   `group` is given, restricted to that group.
#' @param group Export-group name to filter by and label the bundle.
#' @param out_dir Output directory (default `<path>/export`).
#' @param label Bundle label (default: `group` or "export").
#' @param markers,report Also write per-marker fastas / an HTML report.
#' @param run_table2asn Run table2asn if on PATH.
#' @param zip Also write `<out_dir>/<label>.zip`.
#' @return (invisibly) the group output directory.
#' @export
export_bundle <- function(path = ".", ids = NULL, group = NULL, out_dir = file.path(path, "export"),
                          label = NULL, markers = TRUE, report = TRUE, run_table2asn = FALSE, zip = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  have <- DBI::dbGetQuery(con, "SELECT ID FROM stitched")$ID
  if (!is.null(group)) {
    gids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE export_group = ?", params = list(group))$ID
    ids <- if (is.null(ids)) gids else intersect(ids, gids)
  }
  ids <- if (is.null(ids)) have else intersect(ids, have)
  if (length(ids) == 0) {
    stop("No samples with an assembled/annotated sequence to export.", call. = FALSE)
  }
  label <- label %||% group %||% "export"

  gdir <- file.path(out_dir, label)
  unlink(gdir, recursive = TRUE)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  if (markers) dir.create(file.path(gdir, "markers"), showWarnings = FALSE)

  group_fsa <- character(0)
  group_tbl <- character(0)
  marker_pool <- list()

  for (id in ids) {
    row <- DBI::dbGetQuery(con,
      "SELECT st.sequence AS seq, s.Taxon FROM stitched st JOIN samples s ON st.ID = s.ID WHERE st.ID = ?",
      params = list(id))
    seq <- row$seq[1]
    taxon <- row$Taxon[1]
    defline <- sprintf(">%s [organism=%s] [moltype=genomic DNA]", id, taxon)
    writeLines(c(defline, seq), file.path(gdir, paste0(id, ".fsa")))
    tbl <- feature_table_text_con(con, id)
    writeLines(tbl, file.path(gdir, paste0(id, ".tbl")))
    group_fsa <- c(group_fsa, defline, seq)
    group_tbl <- c(group_tbl, tbl)

    if (markers) {
      ann <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
      wm <- function(el, s, e) {
        if (is.na(s) || is.na(e) || e < s) {
          return(invisible())
        }
        rec <- c(sprintf(">%s_%s [organism=%s]", id, el, taxon), substr(seq, s, e))
        writeLines(rec, file.path(gdir, "markers", paste0(id, "_", gsub("[^A-Za-z0-9]", "", el), ".fasta")))
        marker_pool[[el]] <<- c(marker_pool[[el]], rec)
      }
      for (i in seq_len(nrow(ann))) wm(ann$element[i], ann$start[i], ann$stop[i])
      i1 <- ann[ann$element == "ITS1", ]
      i2 <- ann[ann$element == "ITS2", ]
      if (nrow(i1) && nrow(i2)) wm("ITS", i1$start[1], i2$stop[1])
    }
  }

  # combined group files
  writeLines(group_fsa, file.path(gdir, paste0(label, ".fasta")))
  writeLines(paste(group_tbl, collapse = "\n"), file.path(gdir, paste0(label, ".tbl")))
  if (markers) {
    for (el in names(marker_pool)) {
      writeLines(marker_pool[[el]], file.path(gdir, "markers", paste0("all_", gsub("[^A-Za-z0-9]", "", el), ".fasta")))
    }
  }
  # per-sample summary CSV
  utils::write.csv(export_summary_df(con, ids), file.path(gdir, paste0(label, "_summary.csv")), row.names = FALSE)
  # HTML report
  if (isTRUE(report)) {
    tryCatch(export_report_html(con, ids, file.path(gdir, paste0(label, "_report.html")), label),
      error = function(e) message("Report generation skipped: ", conditionMessage(e)))
  }
  # optional table2asn
  if (isTRUE(run_table2asn) && nzchar(unname(Sys.which("table2asn")))) {
    sbt <- Sys.glob(file.path(gdir, "*.sbt"))
    args <- c("-indir", shQuote(gdir), "-a", "s")
    if (length(sbt) > 0) args <- c(args, "-t", shQuote(sbt[1]))
    system2("table2asn", args, stdout = TRUE, stderr = TRUE)
  }

  message("Exported ", length(ids), " sample(s) to ", gdir)
  if (isTRUE(zip)) {
    zf <- paste0(gdir, ".zip")
    old <- setwd(out_dir)
    on.exit(setwd(old), add = TRUE)
    utils::zip(zf, label, flags = "-rq")
  }
  invisible(gdir)
}

#' Export every assigned export group in a project
#' @param path Project directory. @param out_dir Output dir.
#' @return (invisibly) the export groups exported.
#' @export
export_project <- function(path = ".", out_dir = file.path(path, "export")) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  groups <- DBI::dbGetQuery(con, "SELECT DISTINCT export_group FROM samples WHERE export_group IS NOT NULL")$export_group
  DBI::dbDisconnect(con)
  if (length(groups) == 0) {
    stop("No export groups assigned. Assign samples to a group first.", call. = FALSE)
  }
  for (g in groups) export_bundle(path, group = g, out_dir = out_dir, zip = TRUE)
  invisible(groups)
}
