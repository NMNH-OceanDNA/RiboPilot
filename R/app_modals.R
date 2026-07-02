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

#' Annotate & curate modal: feature map + editable boundaries + review
#' @noRd
annotate_modal <- function(id) {
  modalDialog(
    title = paste0("Annotate — ", id),
    size = "l", easyClose = TRUE,
    footer = tagList(
      actionButton("save_annot", "Save edits", class = "btn-primary btn-sm"),
      actionButton("mark_reviewed", "Mark reviewed", class = "btn-sm"),
      modalButton("Close")
    ),
    plotOutput("feature_map", height = "150px"),
    hr(),
    h5("Edit feature boundaries"),
    uiOutput("annot_edit"),
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
