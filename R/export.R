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

# Valid INSDC source modifiers that may be emitted as `[key=value]` deflines.
# This is an ALLOWLIST, not a blocklist: only map-file columns whose (lowercased)
# name is a real GenBank source modifier are emitted. A blocklist would leak the
# internal bookkeeping columns on the `samples` table (taxonomy_switch,
# taxonomy_lock, rank, ...) into every defline as bogus qualifiers, which
# table2asn rejects. `organism`/`moltype` are emitted separately.
SOURCE_MOD_ALLOW <- c(
  # specimen / material
  "specimen_voucher", "culture_collection", "bio_material", "isolate", "strain",
  "clone", "cultivar", "ecotype", "haplotype", "breed", "variety", "serotype",
  "serovar", "sub_species", "subspecies", "sub_strain", "isolation_source",
  # geography / date
  "country", "geo_loc_name", "lat_lon", "altitude", "collection_date",
  "collected_by", "identified_by",
  # organism biology
  "host", "sex", "tissue_type", "dev_stage", "cell_type", "cell_line",
  "organelle", "environmental_sample", "pop_variant", "biotype",
  # free text
  "note")

#' Build GenBank source-modifier tags from a sample's metadata row
#'
#' Extra map-file columns are emitted as `[column=value]` defline modifiers, but
#' only when the column name is a recognized INSDC source modifier
#' ([SOURCE_MOD_ALLOW]); name the map-file columns after the source modifiers you
#' want (e.g. `specimen_voucher`, `country`/`geo_loc_name`, `collection_date`).
#' Unrecognized columns are ignored so internal columns cannot become qualifiers.
#' @noRd
source_modifiers <- function(row) {
  mods <- character(0)
  for (col in names(row)) {
    key <- tolower(gsub("\\s+", "_", trimws(col)))
    if (!(key %in% SOURCE_MOD_ALLOW)) next
    v <- row[[col]][1]
    if (is.null(v) || length(v) == 0 || is.na(v) || !nzchar(trimws(as.character(v)))) next
    mods <- c(mods, sprintf("[%s=%s]", key, trimws(as.character(v))))
  }
  mods
}

# Minimum run of Ns annotated as an assembly_gap; shorter runs are treated as
# ordinary ambiguous bases, not gaps (per INSDC convention).
MIN_ASSEMBLY_GAP <- 10L

#' 1-based start/stop of every N-run of at least `min_len` in a sequence
#' @noRd
find_n_runs <- function(seq, min_len = MIN_ASSEMBLY_GAP) {
  empty <- data.frame(start = integer(0), stop = integer(0))
  if (is.null(seq) || is.na(seq) || !nzchar(seq)) return(empty)
  m <- gregexpr(sprintf("[Nn]{%d,}", min_len), seq, perl = TRUE)[[1]]
  if (length(m) == 1 && m[1] == -1) return(empty)
  len <- attr(m, "match.length")
  data.frame(start = as.integer(m), stop = as.integer(m + len - 1))
}

#' Feature table (.tbl) text for one sample from an open connection
#'
#' Emits the rRNA/spacer features plus an `assembly_gap` feature for each N-run
#' (>= [MIN_ASSEMBLY_GAP]) so table2asn does not reject fabricated gap runs, and
#' marks the first/last real feature partial (`<`/`>`) when it reaches a
#' sequence end or when `force_partial` is set (a curator-flagged partial repeat,
#' so it is not submitted as a complete cistron).
#' @noRd
feature_table_text_con <- function(con, id, force_partial = FALSE) {
  d <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
  hdr <- paste0(">Feature ", id)
  if (nrow(d) == 0) {
    return(paste0(hdr, "\n"))
  }
  st <- DBI::dbGetQuery(con, "SELECT sequence, length FROM stitched WHERE ID = ?", params = list(id))
  seq <- if (nrow(st)) st$sequence[1] else NA_character_
  seqlen <- if (nrow(st) && !is.na(st$length[1])) st$length[1] else max(d$stop)

  gaps <- find_n_runs(seq)
  feats <- rbind(
    data.frame(start = d$start, stop = d$stop, type = d$type, element = d$element,
      kind = "rna", stringsAsFactors = FALSE),
    if (nrow(gaps) > 0) data.frame(start = gaps$start, stop = gaps$stop,
      type = "assembly_gap", element = NA_character_, kind = "gap", stringsAsFactors = FALSE))
  feats <- feats[order(feats$start), , drop = FALSE]
  rna <- which(feats$kind == "rna")
  first_rna <- if (length(rna)) min(rna) else NA_integer_
  last_rna <- if (length(rna)) max(rna) else NA_integer_

  lines <- hdr
  for (i in seq_len(nrow(feats))) {
    if (identical(feats$kind[i], "gap")) {
      lines <- c(lines,
        paste0(feats$start[i], "\t", feats$stop[i], "\tassembly_gap"),
        paste0("\t\t\testimated_length\t", feats$stop[i] - feats$start[i] + 1),
        "\t\t\tgap_type\twithin scaffold",
        "\t\t\tlinkage_evidence\talign_genus")
      next
    }
    s <- if (i == first_rna && (feats$start[i] <= 1 || isTRUE(force_partial))) {
      paste0("<", feats$start[i])
    } else {
      as.character(feats$start[i])
    }
    e <- if (i == last_rna && (feats$stop[i] >= seqlen || isTRUE(force_partial))) {
      paste0(">", feats$stop[i])
    } else {
      as.character(feats$stop[i])
    }
    lines <- c(lines, paste0(s, "\t", e, "\t", feats$type[i]),
      paste0("\t\t\tproduct\t", rrna_element_product(feats$element[i])))
    if (identical(feats$type[i], "misc_RNA")) lines <- c(lines, paste0("\t\t\tnote\t", feats$element[i]))
  }
  paste(lines, collapse = "\n")
}

#' Portal-parity submission-risk flags for one export sequence
#'
#' NCBI's Submission Portal accepts FASTA only for rRNA-ITS and runs ITSx (stock
#' parameters) on it, silently removing sequences it judges chimeric, gapped
#' between ITS regions, or misassembled
#' (https://submit.ncbi.nlm.nih.gov/help/genbank/#annotation_rRNAITS_faq). This
#' re-checks each sequence the same way so a curator sees the risk at export
#' rather than as a rejection. Warn-only: nothing is dropped from the bundle.
#' @return list(chimeric, gap_in_rdna, misassembly, note); `note` is a
#'   `;`-joined string of the flags that fired ("" when clean).
#' @noRd
submission_risks <- function(con, id) {
  chim <- isTRUE(DBI::dbGetQuery(con,
    "SELECT chimeric FROM annotate WHERE ID = ?", params = list(id))$chimeric[1] == 1L)
  ann <- DBI::dbGetQuery(con,
    "SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
  gap_in_rdna <- FALSE
  misassembly <- FALSE
  if (nrow(ann) > 0) {
    st <- DBI::dbGetQuery(con, "SELECT sequence FROM stitched WHERE ID = ?", params = list(id))
    seq <- if (nrow(st)) st$sequence[1] else NA_character_
    gaps <- find_n_runs(seq)
    if (nrow(gaps) > 0) {
      # an assembly gap inside the transcribed span is NCBI's "gap between ITS regions"
      gap_in_rdna <- any(gaps$stop >= min(ann$start) & gaps$start <= max(ann$stop))
    }
    # features must run strictly 5'->3' in canonical order and not overlap
    idx <- match(ann$element, RRNA_ELEMENT_ORDER)
    misassembly <- any(diff(idx) <= 0, na.rm = TRUE) ||
      (nrow(ann) > 1 && any(ann$start[-1] <= ann$stop[-nrow(ann)]))
  }
  flags <- c(if (chim) "chimeric", if (gap_in_rdna) "gap_in_rDNA", if (misassembly) "misassembly")
  list(chimeric = chim, gap_in_rdna = gap_in_rdna, misassembly = misassembly,
    note = paste(flags, collapse = ";"))
}

#' README text mapping the bundle to NCBI's two submission routes
#' @noRd
submission_readme <- function(label) {
  paste(c(
    paste0("RiboPilot GenBank submission bundle: ", label),
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    "HOW NCBI HANDLES rRNA-ITS SUBMISSIONS",
    "-------------------------------------",
    "The NCBI Submission Portal (https://submit.ncbi.nlm.nih.gov) accepts FASTA",
    "ONLY for rRNA-ITS sequences. It runs ITSx (stock parameters, no organism",
    "tuning) and writes the annotation itself as a SINGLE misc_RNA feature over",
    "the whole region, e.g. (accession PZ692210):",
    "",
    "  misc_RNA  <1..>N",
    "     /note=\"contains small subunit ribosomal RNA, internal transcribed",
    "            spacer 1, 5.8S ribosomal RNA, internal transcribed spacer 2,",
    "            and large subunit ribosomal RNA\"",
    "",
    "Per-element annotation is NOT submitted through the Portal, and the .tbl",
    "feature tables in this bundle are NOT used there.",
    "",
    paste0("  -> Portal route: upload ", label, ".fasta only; NCBI annotates for you."),
    "",
    "The .tbl / .sbt / table2asn output are for the alternate table2asn / BankIt",
    "route only, where a feature table IS accepted.",
    "",
    "SUBMISSION-RISK CHECKS",
    "----------------------",
    "NCBI's ITSx evaluation may silently remove sequences it judges chimeric,",
    "gapped between ITS regions, or misassembled. RiboPilot re-checks each",
    paste0("sequence the same way; see the submission_risk column in ", label, "_summary.csv."),
    "Flagged sequences are still exported (warn-only); review before submitting:",
    "  chimeric     ITSx flagged mixed origins",
    "  gap_in_rDNA  an assembly-gap (N-run) falls inside the 18S..28S span",
    "  misassembly  features overlap or are out of 5'->3' order",
    ""), collapse = "\n")
}

#' Per-sample export summary data.frame
#'
#' Includes a `submission_risk` column ([submission_risks]) so the CSV and HTML
#' report both surface portal-parity risk flags.
#' @noRd
export_summary_df <- function(con, ids) {
  d <- DBI::dbGetQuery(con, sprintf(
    "SELECT s.ID, s.Taxon AS organism, s.export_group, st.reference_accession AS reference,
       st.length, an.n_features, an.missing, an.reviewed, an.partial, an.problematic
     FROM samples s LEFT JOIN stitched st ON s.ID = st.ID LEFT JOIN annotate an ON s.ID = an.ID
     WHERE s.ID IN (%s)", paste(sprintf("'%s'", ids), collapse = ",")))
  d$submission_risk <- vapply(d$ID, function(id) submission_risks(con, id)$note, character(1))
  d
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
  hdr <- c("ID", "organism", "reference", "length", "n_features", "missing", "reviewed", "submission_risk")
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
#' @param template Optional path to a GenBank submission template (`.sbt`);
#'   copied into the bundle as `<label>.sbt` so table2asn picks it up.
#' @param run_table2asn Run table2asn if on PATH.
#' @param zip Also write `<out_dir>/<label>.zip`.
#' @return (invisibly) the group output directory.
#' @export
export_bundle <- function(path = ".", ids = NULL, group = NULL, out_dir = file.path(path, "export"),
                          label = NULL, markers = TRUE, report = TRUE, template = NULL,
                          run_table2asn = FALSE, zip = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # exclude any sample a curator flagged problematic ("do not submit")
  have <- DBI::dbGetQuery(con,
    "SELECT st.ID FROM stitched st LEFT JOIN annotate an ON st.ID = an.ID
     WHERE COALESCE(an.problematic, '') != 'yes'")$ID
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
      "SELECT st.sequence AS seq, s.* FROM stitched st JOIN samples s ON st.ID = s.ID WHERE st.ID = ?",
      params = list(id))
    seq <- row$seq[1]
    taxon <- row$Taxon[1]
    defline <- paste(c(
      sprintf(">%s [organism=%s] [moltype=genomic DNA]", id, taxon),
      source_modifiers(row)), collapse = " ")
    writeLines(c(defline, seq), file.path(gdir, paste0(id, ".fsa")))
    is_partial <- identical(DBI::dbGetQuery(con,
      "SELECT partial FROM annotate WHERE ID = ?", params = list(id))$partial[1], "yes")
    tbl <- feature_table_text_con(con, id, force_partial = is_partial)
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
  # per-sample summary CSV (+ warn-only portal-parity submission-risk gate)
  smry <- export_summary_df(con, ids)
  risky <- which(nzchar(smry$submission_risk))
  if (length(risky) > 0) {
    warning("Submission-risk flags on ", length(risky), " of ", nrow(smry),
      " sample(s) (still exported; see submission_risk column): ",
      paste(sprintf("%s[%s]", smry$ID[risky], smry$submission_risk[risky]), collapse = ", "),
      call. = FALSE)
  }
  utils::write.csv(smry, file.path(gdir, paste0(label, "_summary.csv")), row.names = FALSE)
  writeLines(submission_readme(label), file.path(gdir, "README_SUBMISSION.txt"))
  # HTML report
  if (isTRUE(report)) {
    tryCatch(export_report_html(con, ids, file.path(gdir, paste0(label, "_report.html")), label),
      error = function(e) message("Report generation skipped: ", conditionMessage(e)))
  }
  # submission template (.sbt): copy in so table2asn (-t) picks it up
  if (!is.null(template)) {
    if (!file.exists(template)) {
      warning("Submission template not found: ", template)
    } else {
      file.copy(template, file.path(gdir, paste0(label, ".sbt")), overwrite = TRUE)
    }
  }

  # optional table2asn
  if (isTRUE(run_table2asn) && nzchar(unname(Sys.which("table2asn")))) {
    sbt <- Sys.glob(file.path(gdir, "*.sbt"))
    args <- c("-indir", shQuote(gdir), "-a", "s")
    if (length(sbt) > 0) args <- c(args, "-t", shQuote(sbt[1]))
    system2("table2asn", args, stdout = TRUE, stderr = TRUE)
  }

  message("Exported ", length(ids), " sample(s) to ", gdir,
    "\n  GenBank Portal takes FASTA only: upload ", label,
    ".fasta; NCBI auto-annotates the rRNA-ITS region (see README_SUBMISSION.txt).")
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
