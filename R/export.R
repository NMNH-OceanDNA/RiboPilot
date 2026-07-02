# GenBank submission bundle export.
#
# For each sample: a source-modifier fasta (.fsa) of the assembled nuclear rDNA
# region, a 5-column feature table (.tbl) with the 18S/ITS1/5.8S/ITS2/28S model,
# and per-marker fastas (18S / 28S / full ITS / ITS1 / ITS2) for barcode
# submissions. These are the inputs NCBI's table2asn consumes; if table2asn is
# installed it is run to produce a .sqn.

#' Feature table (.tbl) text for one sample from an open connection
#' @noRd
feature_table_text_con <- function(con, id) {
  d <- DBI::dbGetQuery(con,
    "SELECT element, type, start, stop FROM annotations WHERE ID = ? ORDER BY start",
    params = list(id))
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
    if (d$type[i] == "misc_RNA") {
      lines <- c(lines, paste0("\t\t\tnote\t", d$element[i]))
    }
  }
  paste(lines, collapse = "\n")
}

#' Build a GenBank submission bundle for selected samples
#'
#' @param path Project directory.
#' @param ids Sample IDs to export (default: all with a stitched sequence).
#' @param out_dir Output directory (default `<path>/export`).
#' @param markers Also write per-marker fastas (18S/28S/ITS/ITS1/ITS2).
#' @param run_table2asn If TRUE and `table2asn` is on PATH, run it to produce a
#'   `.sqn` (a template `.sbt` in `out_dir` is used if present).
#' @param zip If TRUE, also write `<out_dir>.zip`.
#' @return (invisibly) the output directory path.
#' @export
export_bundle <- function(path = ".", ids = NULL, out_dir = file.path(path, "export"),
                          markers = TRUE, run_table2asn = FALSE, zip = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  have <- DBI::dbGetQuery(con, "SELECT ID FROM stitched")$ID
  ids <- if (is.null(ids)) have else intersect(ids, have)
  if (length(ids) == 0) {
    stop("No samples with an assembled/annotated sequence to export.", call. = FALSE)
  }

  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (markers) dir.create(file.path(out_dir, "markers"), showWarnings = FALSE)

  marker_pool <- list() # element -> character vector of fasta records (across samples)

  for (id in ids) {
    row <- DBI::dbGetQuery(con,
      "SELECT st.sequence AS seq, s.Taxon FROM stitched st JOIN samples s ON st.ID = s.ID WHERE st.ID = ?",
      params = list(id))
    seq <- row$seq[1]
    taxon <- row$Taxon[1]

    # .fsa with source modifiers in the defline
    defline <- sprintf(">%s [organism=%s] [moltype=genomic DNA]", id, taxon)
    writeLines(c(defline, seq), file.path(out_dir, paste0(id, ".fsa")))

    # .tbl feature table
    writeLines(feature_table_text_con(con, id), file.path(out_dir, paste0(id, ".tbl")))

    # per-marker fastas
    if (markers) {
      ann <- DBI::dbGetQuery(con,
        "SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start", params = list(id))
      write_marker <- function(el, s, e) {
        if (is.na(s) || is.na(e) || e < s) {
          return(invisible())
        }
        sub <- substr(seq, s, e)
        rec <- c(sprintf(">%s_%s [organism=%s]", id, el, taxon), sub)
        writeLines(rec, file.path(out_dir, "markers", paste0(id, "_", gsub("[^A-Za-z0-9]", "", el), ".fasta")))
        marker_pool[[el]] <<- c(marker_pool[[el]], rec)
      }
      for (i in seq_len(nrow(ann))) write_marker(ann$element[i], ann$start[i], ann$stop[i])
      # full ITS = ITS1 start .. ITS2 stop (if both present)
      its1 <- ann[ann$element == "ITS1", ]
      its2 <- ann[ann$element == "ITS2", ]
      if (nrow(its1) && nrow(its2)) write_marker("ITS", its1$start[1], its2$stop[1])
    }
  }

  # combined per-marker multifastas
  if (markers) {
    for (el in names(marker_pool)) {
      writeLines(marker_pool[[el]],
        file.path(out_dir, "markers", paste0("all_", gsub("[^A-Za-z0-9]", "", el), ".fasta")))
    }
  }

  # optional table2asn
  if (isTRUE(run_table2asn) && nzchar(unname(Sys.which("table2asn")))) {
    sbt <- Sys.glob(file.path(out_dir, "*.sbt"))
    args <- c("-indir", shQuote(out_dir), "-a", "s")
    if (length(sbt) > 0) args <- c(args, "-t", shQuote(sbt[1]))
    system2("table2asn", args, stdout = TRUE, stderr = TRUE)
  }

  message("Exported ", length(ids), " sample(s) to ", out_dir)
  if (isTRUE(zip)) {
    zf <- paste0(out_dir, ".zip")
    old <- setwd(dirname(out_dir))
    on.exit(setwd(old), add = TRUE)
    utils::zip(zf, basename(out_dir), flags = "-rq")
    message("Zipped bundle: ", zf)
  }
  invisible(out_dir)
}
