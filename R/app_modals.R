# Shiny modal dialogs (Assemble / Annotate) + their plot/table helpers.

#' Assemble details modal: pick the selected path; view scores + scaffolds
#' @noRd
assemble_modal <- function(id) {
  modalDialog(
    title = paste0("Assemble — ", id),
    size = "l", easyClose = TRUE,
    footer = tagList(
      actionButton("set_selected", "Set as selected path", class = "btn-primary btn-sm"),
      modalButton("Close")
    ),
    uiOutput("sel_path_ui"),
    h5("Candidate paths (auto-selected in green)"),
    reactable::reactableOutput("assemble_scores"),
    hr(),
    h5("Scaffolds"),
    reactable::reactableOutput("assemble_scaffolds")
  )
}

#' 5'/3' nucleotide-boundary nudge control (mirrors MitoPilot rRNA nudge)
#' @noRd
nt_nudge <- function(id, sym) {
  actionButton(id, sym, class = "btn-sm",
    style = "min-width:34px; padding:1px 8px; font-weight:bold;")
}

#' Annotate & curate modal: feature map, per-feature nucleotide editor + review
#' @noRd
annotate_modal <- function(id) {
  modalDialog(
    title = paste0("Annotate - ", id),
    size = "l", easyClose = TRUE,
    footer = tagList(
      actionButton("mark_reviewed", "Mark reviewed", class = "btn-sm"),
      modalButton("Close")
    ),
    h5("Feature / reference synteny map"),
    plotOutput("feature_map", height = "160px"),
    h5("Read coverage"),
    plotOutput("coverage_plot", height = "130px"),
    hr(),
    h5("Edit feature boundaries"),
    div(style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
      div(style = "min-width:170px;", uiOutput("nudge_element_ui")),
      shinyWidgets::radioGroupButtons("nudge_step", "Step (nt)",
        choices = c(1, 3, 10, 50), selected = 1, size = "sm"),
      div(
        tags$label("5' end", class = "control-label"),
        div(style = "display:flex; gap:4px;",
          nt_nudge("nt_5_out", "+"), nt_nudge("nt_5_in", "–"))),
      div(
        tags$label("3' end", class = "control-label"),
        div(style = "display:flex; gap:4px;",
          nt_nudge("nt_3_in", "–"), nt_nudge("nt_3_out", "+")))
    ),
    div(style = "color:#888; font-size:0.85em; margin-top:4px;",
      "+ extends the boundary outward, – trims it inward."),
    div(class = "code-block", style = "margin-top:6px; white-space:normal; word-break:break-all;",
      uiOutput("nudge_seq")),
    hr(),
    h5("Feature model"),
    reactable::reactableOutput("feature_tbl")
  )
}

#' Update / run modal: launch the Nextflow workflow with live progress
#' @noRd
run_modal <- function(mode, n) {
  wf <- if (mode == "Annotate") "Annotation (WF2)" else "Assembly (WF1)"
  modalDialog(
    title = div(style = "display:flex; justify-content:space-between; align-items:center; height:38px;",
      span(sprintf("%s — %d sample(s) to update", wf, n)),
      span(id = "rp_gears", class = "gears paused")),
    size = "l",
    shinyWidgets::prettyCheckbox("nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:"),
    div(class = "code-block", textOutput("nf_cmd_text")),
    shinyjs::hidden(div(id = "nf_progress",
      hr(), h5("Progress:"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        textOutput("nf_prog_header"),
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit;",
          textOutput("nf_prog_processes", container = tags$code)),
        textOutput("nf_prog_summary"),
        textOutput("nf_prog_footer")))),
    footer = tagList(
      actionButton("nf_start", "Run from App", class = "btn-success btn-sm"),
      shinyjs::hidden(actionButton("nf_stop", "Stop / Interrupt", class = "btn-sm")),
      actionButton("nf_close", "Close", class = "btn-sm"))
  )
}

#' Assign-export-group modal: name the group for the selected samples
#' @noRd
group_modal <- function() {
  modalDialog(
    title = "Assign export group", size = "s", easyClose = TRUE,
    textInput("group_name", "Group name", placeholder = "e.g. batch_2026_01"),
    footer = tagList(
      actionButton("assign_group_go", "Assign", class = "btn-primary btn-sm"),
      modalButton("Cancel"))
  )
}

#' Element colors shared by the app feature map
#' @noRd
RRNA_ELEMENT_COLORS <- c(
  "18S" = "#0072B2", "ITS1" = "#E69F00", "5.8S" = "#F0E442",
  "ITS2" = "#D55E00", "28S" = "#009E73", "ETS" = "#BBBBBB", "IGS" = "#BBBBBB"
)

#' Read-coverage track for a sample, with feature regions shaded
#' @noRd
plot_coverage <- function(id) {
  if (is.null(id)) {
    return(NULL)
  }
  cov <- app_query("SELECT depths FROM coverage WHERE ID = ?", params = list(id))
  if (nrow(cov) == 0 || is.na(cov$depths[1]) || !nzchar(cov$depths[1])) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", 0, 0, label = "no coverage data (run WF1)") +
      ggplot2::theme_void())
  }
  dep <- as.numeric(strsplit(cov$depths[1], ",")[[1]])
  df <- data.frame(pos = seq_along(dep), depth = dep)
  feats <- app_query("SELECT element, start, stop FROM annotations WHERE ID = ? ORDER BY start",
    params = list(id))
  g <- ggplot2::ggplot(df, ggplot2::aes(pos, depth))
  if (nrow(feats) > 0) {
    g <- g + ggplot2::geom_rect(data = feats, inherit.aes = FALSE,
      ggplot2::aes(xmin = start, xmax = stop, ymin = 0, ymax = Inf, fill = element),
      alpha = 0.15) +
      ggplot2::scale_fill_manual(values = RRNA_ELEMENT_COLORS, guide = "none")
  }
  g + ggplot2::geom_area(fill = "#444444", alpha = 0.85) +
    ggplot2::labs(x = "position (bp)", y = "depth") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}

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
