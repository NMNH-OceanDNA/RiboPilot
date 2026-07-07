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
#' A collapsible <details> section (open by default unless open = FALSE).
#' @noRd
rp_details <- function(summary, ..., open = TRUE) {
  d <- tags$details(tags$summary(tags$b(summary)), div(style = "padding:6px 2px 2px;", ...))
  if (isTRUE(open)) d$attribs$open <- NA
  d
}

annotate_modal <- function(id) {
  modalDialog(
    title = paste0("Annotate - ", id),
    size = "l", easyClose = FALSE,
    footer = tagList(
      actionButton("mark_reviewed", "Mark reviewed", class = "btn-sm"),
      modalButton("Close")
    ),
    # ---- gene list (drives the active feature for the editor + alignment) ----
    h5("Feature model", style = "margin:2px 0 4px;"),
    div(style = "color:#888; font-size:0.82em; margin-bottom:4px;",
      "Click a feature to make it active for editing + reference alignment below."),
    reactable::reactableOutput("feature_tbl"),
    div(style = "margin:6px 0;",
      actionButton("trim_btn", "Trim assembly...", icon = icon("scissors"), class = "btn-sm")),
    tags$hr(style = "margin:10px 0;"),

    # ---- collapsible detail panels (coverage above synteny, per request) ----
    rp_details("Read coverage", plotOutput("coverage_plot", height = "130px")),

    rp_details("Edit active feature boundaries",
      div(style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
        div(style = "font-weight:bold; align-self:center;", textOutput("feat_edit_active", inline = TRUE)),
        div(class = "rp-step-box", style = "width:70px;",
          numericInput("feat_step", "Step (nt)", value = 1, min = 1, max = 50, step = 1, width = "70px")),
        shinyjs::hidden(div(id = "feat_edit_ctrls",
          style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
          div(tags$label("5' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_5_out", "+"), nt_nudge("nt_5_in", "-"))),
          div(tags$label("3' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_3_in", "-"), nt_nudge("nt_3_out", "+"))))),
        div(style = "align-self:center; display:flex; gap:6px;",
          actionButton("feat_edit", "Edit", icon = icon("pen"), class = "btn-sm btn-primary"),
          shinyjs::hidden(actionButton("feat_save", "Save", icon = icon("floppy-disk"), class = "btn-sm btn-success")),
          shinyjs::hidden(actionButton("feat_reset", "Reset", icon = icon("rotate-left"), class = "btn-sm")))
      ),
      div(style = "color:#888; font-size:0.82em; margin-top:4px;",
        "Click Edit to adjust boundaries; + extends outward, - trims inward. Save to keep, Reset to discard."),
      div(class = "code-block", style = "margin-top:6px; white-space:normal; word-break:break-all;",
        uiOutput("nudge_seq"))),

    rp_details("Reference alignment (active feature)",
      uiOutput("msa_header"),
      uiOutput("feat_align_ui")),

    rp_details("Curation flags", open = FALSE,
      uiOutput("curate_notes"),
      reactable::reactableOutput("flag_tbl")),

    rp_details("Notes", open = FALSE,
      textAreaInput("annotate_notes", label = NULL, width = "100%",
        placeholder = "Notes for this sample's annotation"))
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
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("nf_cmd_text")),
    shinyjs::hidden(div(id = "nf_progress",
      hr(style = "margin:12px 0;"), h5("Progress:", style = "margin:0 0 4px;"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit; white-space:pre-wrap;",
          textOutput("nf_prog_all", container = tags$code))))),
    footer = tagList(
      actionButton("nf_start", "Run from App", class = "btn-success btn-sm"),
      shinyjs::hidden(actionButton("nf_stop", "Stop / Interrupt", class = "btn-sm")),
      actionButton("nf_close", "Close", class = "btn-sm"))
  )
}

#' Trim-assembly modal: cut the stitched assembly to the annotated feature span
#' or a custom range (shifts annotations + coverage).
#' @noRd
trim_modal <- function(n, a_from, a_to) {
  modalDialog(
    title = "Trim assembly", size = "m",
    tags$p(class = "text-muted", style = "font-size:0.85em;",
      sprintf("Current assembly length: %d bp. Annotated span: %d..%d bp.", n, a_from, a_to)),
    shinyWidgets::prettyRadioButtons("trim_mode", "Trim to:",
      choices = c("Annotated feature boundaries" = "annotated", "Custom range" = "custom"),
      selected = "annotated", status = "primary"),
    conditionalPanel("input.trim_mode == 'custom'",
      div(style = "display:flex; gap:1em;",
        numericInput("trim_from", "From (bp)", value = a_from, min = 1, max = n),
        numericInput("trim_to", "To (bp)", value = a_to, min = 1, max = n))),
    opts_help("Trimming cuts the flanks off the stitched assembly and shifts every ",
      "feature coordinate + the coverage track to match. Features that fall entirely ",
      "outside the kept range are removed."),
    footer = tagList(
      actionButton("trim_apply", "Trim", icon = icon("scissors"), class = "btn-danger btn-sm"),
      modalButton("Cancel"))
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

#' State-picker modal: choose the state to set on the selected samples.
#' Codes: 0 hold/wait (won't run), 1 ready to run, 2 complete, 3 failed,
#' 4 in progress. Only 1 and 4 are picked up by the next pipeline run.
#' @noRd
state_modal <- function(mode, current = NULL) {
  verb <- if (mode == "annotate") "annotate" else "assemble"
  modalDialog(
    title = "Set state for selected samples", size = "m",
    shinyWidgets::prettyRadioButtons("new_state", label = NULL,
      choices = c("Hold / wait (won't run)" = "0",
        `Ready to run` = "1", "In progress" = "4",
        "Complete" = "2", "Failed / problem" = "3"),
      selected = as.character(current %||% 1L), shape = "square", status = "primary"),
    opts_help("Only 'Ready to run' and 'In progress' samples are picked up the ",
      paste0("next time you run ", verb, "; 'Hold / wait' keeps them out of the run.")),
    footer = tagList(
      actionButton("update_state", "Update", class = "btn-primary btn-sm"),
      modalButton("Cancel"))
  )
}

#' Render a per-feature nucleotide alignment as highlighted monospace HTML
#' (mismatches in red), wrapped at 80 columns.
#' @noRd
render_alignment_html <- function(sample_aln, ref_aln) {
  sa <- strsplit(sample_aln, "")[[1]]
  ra <- strsplit(ref_aln, "")[[1]]
  n <- max(length(sa), length(ra))
  length(sa) <- n; length(ra) <- n
  sa[is.na(sa)] <- "-"; ra[is.na(ra)] <- "-"
  hl <- function(ch, mm) ifelse(mm, sprintf("<span style='color:#d00;font-weight:bold'>%s</span>", ch), ch)
  width <- 80
  blocks <- vapply(seq(1, n, by = width), function(start) {
    idx <- start:min(start + width - 1, n)
    mm <- sa[idx] != ra[idx]
    paste0(
      "sample  ", paste(hl(sa[idx], mm), collapse = ""), "\n",
      "        ", paste(ifelse(mm, " ", "|"), collapse = ""), "\n",
      "ref     ", paste(hl(ra[idx], mm), collapse = ""), "\n")
  }, character(1))
  shiny::HTML(paste0("<pre style='margin:0; font-size:0.8em;'>",
    paste(blocks, collapse = "\n"), "</pre>"))
}

#' Element colors shared by the app feature map
#' @noRd
# Colorblind-friendly (Okabe-Ito) palette, chosen to keep the 5 features
# maximally distinct: 18S blue and 28S vermillion are now clearly different
# (the old sky-blue 28S read too close to the 18S blue).
RRNA_ELEMENT_COLORS <- c(
  "18S" = "#0072B2", "ITS1" = "#E69F00", "5.8S" = "#009E73",
  "ITS2" = "#CC79A7", "28S" = "#D55E00", "ETS" = "#BBBBBB", "IGS" = "#BBBBBB"
)

#' Darken a set of (named) hex colors by a factor (for the coverage bands).
#' @noRd
darken_colors <- function(cols, f = 0.72) {
  m <- grDevices::col2rgb(cols) * f
  stats::setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255), names(cols))
}

#' Pale tint of (named) hex colors (mix f toward white) for table row shading.
#' @noRd
tint_colors <- function(cols, f = 0.80) {
  m <- grDevices::col2rgb(cols)
  m <- m + (255 - m) * f
  stats::setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255), names(cols))
}

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
  ymax <- max(df$depth, 1, na.rm = TRUE)
  dark <- darken_colors(RRNA_ELEMENT_COLORS)
  g <- ggplot2::ggplot(df, ggplot2::aes(pos, depth))
  if (nrow(feats) > 0) {
    feats$mid <- (feats$start + feats$stop) / 2
    g <- g +
      # darker feature bands for visibility
      ggplot2::geom_rect(data = feats, inherit.aes = FALSE,
        ggplot2::aes(xmin = start, xmax = stop, ymin = 0, ymax = Inf, fill = element),
        alpha = 0.28) +
      ggplot2::scale_fill_manual(values = dark, guide = "none")
  }
  g <- g + ggplot2::geom_area(fill = "#444444", alpha = 0.85)
  if (nrow(feats) > 0) {
    # gene names above the corresponding boxes
    g <- g + ggplot2::geom_text(data = feats, inherit.aes = FALSE,
      ggplot2::aes(x = mid, y = ymax * 1.06, label = element, colour = element),
      size = 3.2, fontface = "bold", vjust = 0) +
      ggplot2::scale_colour_manual(values = dark, guide = "none") +
      ggplot2::coord_cartesian(ylim = c(0, ymax * 1.16), clip = "off")
  }
  g + ggplot2::labs(x = "position (bp)", y = "depth") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(14, 6, 6, 6))
}

#' gggenes synteny map: the sample's annotations over its BLAST reference's
#' (when the reference has been fetched + annotated), aligned by element color.
#' @noRd
