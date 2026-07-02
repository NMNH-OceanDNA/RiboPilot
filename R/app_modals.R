# Shiny modal dialogs (Assemble / Annotate & Curate / Export) + their helpers.

#' Assemble details modal: candidate paths, scores, scaffolds
#' @noRd
assemble_modal <- function(id) {
  shiny::modalDialog(
    title = paste0("Assemble — ", id),
    size = "l", easyClose = TRUE, footer = shiny::modalButton("Close"),
    shiny::h5("Candidate paths (scored; auto-selected in green)"),
    reactable::reactableOutput("assemble_scores"),
    shiny::hr(),
    shiny::h5("Scaffolds"),
    reactable::reactableOutput("assemble_scaffolds")
  )
}

#' Annotate & curate modal: feature map + editable feature table
#' @noRd
annotate_modal <- function(id) {
  shiny::modalDialog(
    title = paste0("Annotate & Curate — ", id),
    size = "l", easyClose = TRUE, footer = shiny::modalButton("Close"),
    shiny::plotOutput("feature_map", height = "160px"),
    shiny::hr(),
    shiny::h5("Feature model (18S / ITS1 / 5.8S / ITS2 / 28S)"),
    reactable::reactableOutput("feature_tbl")
  )
}

#' Export modal: GenBank feature table preview + downloads
#' @noRd
export_modal <- function(id) {
  shiny::modalDialog(
    title = paste0("Export — ", id),
    size = "l", easyClose = TRUE,
    footer = shiny::tagList(
      shiny::downloadButton("download_tbl", "Download .tbl"),
      shiny::downloadButton("download_fsa", "Download .fsa"),
      shiny::modalButton("Close")
    ),
    shiny::h5("GenBank feature table (5-column .tbl)"),
    shiny::tags$pre(
      style = "background:#f6f8fa; padding:0.75em; max-height:340px; overflow:auto;",
      shiny::textOutput("export_tbl_preview", container = shiny::tags$code)
    )
  )
}

#' Element colors shared by the app feature map
#' @noRd
RRNA_ELEMENT_COLORS <- c(
  "18S" = "#0072B2", "ITS1" = "#E69F00", "5.8S" = "#F0E442",
  "ITS2" = "#D55E00", "28S" = "#009E73", "ETS" = "#BBBBBB", "IGS" = "#BBBBBB"
)

#' gggenes feature-map of a sample's annotations
#' @noRd
plot_feature_map <- function(id) {
  if (is.null(id)) {
    return(NULL)
  }
  d <- app_query(
    "SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start",
    params = list(id))
  if (nrow(d) == 0) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", 0, 0, label = "no annotations yet") +
      ggplot2::theme_void())
  }
  d$molecule <- id
  ggplot2::ggplot(d, ggplot2::aes(xmin = start, xmax = stop, y = molecule,
    fill = element, label = element)) +
    gggenes::geom_gene_arrow(arrowhead_height = grid::unit(5, "mm"),
      arrow_body_height = grid::unit(4.5, "mm")) +
    gggenes::geom_gene_label(align = "centre", grow = FALSE) +
    ggplot2::scale_fill_manual(values = RRNA_ELEMENT_COLORS, guide = "none") +
    ggplot2::labs(x = "position (bp)", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Build the GenBank 5-column feature table (.tbl) text for a sample
#'
#' 18S/5.8S/28S -> `rRNA` with /product; ITS1/ITS2 -> `misc_RNA` with /product +
#' /note. Terminal features touching the sequence ends are flagged partial
#' (`<`/`>`), matching NCBI conventions for rRNA-ITS submissions.
#' @noRd
build_feature_table <- function(id) {
  if (is.null(id)) {
    return("")
  }
  d <- app_query(
    "SELECT element, type, start, stop FROM annotations WHERE ID = ? ORDER BY start",
    params = list(id))
  hdr <- paste0(">Feature ", id)
  if (nrow(d) == 0) {
    return(paste0(hdr, "\n(no annotations)"))
  }
  seqlen <- app_query("SELECT length FROM stitched WHERE ID = ?", params = list(id))$length
  seqlen <- if (length(seqlen) == 0 || is.na(seqlen)) max(d$stop) else seqlen

  lines <- hdr
  for (i in seq_len(nrow(d))) {
    s <- d$start[i]
    e <- d$stop[i]
    s_txt <- if (i == 1 && s <= 1) paste0("<", s) else as.character(s)
    e_txt <- if (i == nrow(d) && e >= seqlen) paste0(">", e) else as.character(e)
    lines <- c(lines, paste0(s_txt, "\t", e_txt, "\t", d$type[i]))
    lines <- c(lines, paste0("\t\t\tproduct\t", rrna_element_product(d$element[i])))
    if (d$type[i] == "misc_RNA") {
      lines <- c(lines, paste0("\t\t\tnote\t", d$element[i]))
    }
  }
  paste(lines, collapse = "\n")
}
