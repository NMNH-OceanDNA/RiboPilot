# Shiny modal dialogs (Assemble / Annotate) + their plot/table helpers.

#' Coloured yes/no/na status badge (mirrors MitoPilot's details-modal badges).
#'
#' `invert` flips the colour mapping so "yes" reads as a warning (used for
#' Problem). `neutral_no` renders "no" in grey rather than green (used for
#' Partial, where "no" just means complete, not "good").
#' @noRd
status_badge <- function(label, state, invert = FALSE, neutral_no = FALSE) {
  val <- if (is.null(state) || length(state) == 0 || is.na(state)) "na" else as.character(state)
  color_val <- if (neutral_no && val == "no") "na" else val
  bg <- if (color_val == "yes") { if (invert) "#fde8d0" else "#d4edda" }
    else if (color_val == "no") { if (invert) "#d4edda" else "#fde8d0" } else "#e9ecef"
  fg <- if (color_val == "yes") { if (invert) "#7d4a1e" else "#2d6a4f" }
    else if (color_val == "no") { if (invert) "#2d6a4f" else "#7d4a1e" } else "#6c757d"
  span(style = paste0("background:", bg, "; color:", fg,
      "; border-radius:3px; padding:2px 8px; font-size:0.75em; font-weight:600; white-space:nowrap;"),
    paste0(label, ": ", toupper(val)))
}

#' Sample taxon header line for the detail modals (organism + taxid + rank).
#' @noRd
modal_taxon_header <- function(id) {
  s <- app_query("SELECT Taxon, sci_name, taxid, rank FROM samples WHERE ID = ?", params = list(id))
  if (nrow(s) == 0) {
    return(NULL)
  }
  name <- if (!is.na(s$sci_name[1]) && nzchar(s$sci_name[1])) s$sci_name[1] else s$Taxon[1]
  meta <- c(
    if (!is.na(s$taxid[1]) && nzchar(as.character(s$taxid[1])))
      sprintf("<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s' target='_blank'>taxid %s</a>",
        s$taxid[1], s$taxid[1]),
    if (!is.na(s$rank[1]) && nzchar(s$rank[1])) s$rank[1])
  div(style = "margin:-6px 0 8px; color:#333; font-size:0.92em;",
    tags$b(name %||% id),
    if (length(meta) > 0) HTML(paste0(" &middot; ", paste(meta, collapse = " &middot; "))),
    if (!is.na(name) && !is.na(s$Taxon[1]) && !identical(s$Taxon[1], name))
      tags$span(style = "color:#999;", sprintf(" (input: %s)", s$Taxon[1])))
}

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
    modal_taxon_header(id),
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

annotate_modal <- function(id, path = 1L, scaffold = 1L) {
  n_units <- tryCatch(app_query("SELECT COUNT(*) n FROM annotate WHERE ID = ?", params = list(id))$n[1],
    error = function(e) 1L)
  unit_lbl <- if (isTRUE(n_units > 1)) sprintf(" [path %s / scaffold %s]", path, scaffold) else ""
  # Seed the notes box from the DB. Re-opening the same unit assigns an
  # identical() value to sel_unit, which does not invalidate the reactiveVal, so
  # the server-side reset observer never fires; a value-less textAreaInput then
  # reports "" and the write observer persists it over the stored note.
  notes <- tryCatch(app_query(
    "SELECT annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))$annotate_notes[1],
    error = function(e) NA_character_)
  modalDialog(
    title = paste0("Annotate - ", id, unit_lbl),
    size = "l", easyClose = FALSE,
    footer = tagList(
      actionButton("mark_reviewed", "Mark reviewed", class = "btn-sm"),
      actionButton("mark_problematic", "Problem", class = "btn-sm"),
      actionButton("mark_partial", "Partial", class = "btn-sm"),
      actionButton("lock_and_close", "Lock & Close", class = "btn-primary btn-sm"),
      modalButton("Close")
    ),
    modal_taxon_header(id),
    div(style = "margin:-2px 0 8px; display:flex; gap:6px; flex-wrap:wrap;",
      uiOutput("modal_status_badges", inline = TRUE)),
    # ---- gene list (drives the active feature for the editor + alignment) ----
    h5("Feature model", style = "margin:2px 0 4px;"),
    div(style = "color:#888; font-size:0.82em; margin-bottom:4px;",
      "Click a feature to make it active for editing + reference alignment below."),
    reactable::reactableOutput("feature_tbl"),
    div(style = "margin:6px 0; display:flex; gap:6px;",
      actionButton("trim_btn", "Trim assembly...", icon = icon("scissors"), class = "btn-sm"),
      actionButton("restore_btn", "Restore assembly", icon = icon("rotate-left"), class = "btn-sm")),
    tags$hr(style = "margin:10px 0;"),

    # ---- collapsible detail panels (coverage above synteny, per request) ----
    rp_details("Read coverage", plotOutput("coverage_plot", height = "130px")),

    rp_details("Edit active feature boundaries",
      # normalize control heights so the step box, +/- buttons, and Save/Reset
      # sit on one baseline
      tags$style(HTML(paste0(
        "#feat_edit_row .form-group{margin-bottom:0;} ",
        "#feat_edit_row .control-label{margin-bottom:3px; font-size:0.85em;} ",
        "#feat_edit_row .pretty{margin:0;}"))),
      div(id = "feat_edit_row",
        style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
        div(style = "align-self:flex-end; padding-bottom:6px; font-weight:bold;",
          textOutput("feat_edit_active", inline = TRUE)),
        div(style = "align-self:flex-end; padding-bottom:6px;",
          shinyWidgets::prettyCheckbox("feat_edit", "Edit", value = FALSE,
            status = "primary", inline = TRUE)),
        shinyjs::hidden(div(id = "feat_edit_ctrls",
          style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
          div(class = "rp-step-box", style = "width:70px;",
            numericInput("feat_step", "Step (nt)", value = 1, min = 1, max = 50, step = 1, width = "70px")),
          div(tags$label("5' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_5_out", "+"), nt_nudge("nt_5_in", "-"))),
          div(tags$label("3' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_3_in", "-"), nt_nudge("nt_3_out", "+"))),
          div(style = "display:flex; gap:6px; align-items:flex-end;",
            actionButton("feat_save", "Save", icon = icon("floppy-disk"), class = "btn-sm btn-success"),
            actionButton("feat_reset", "Reset", icon = icon("rotate-left"), class = "btn-sm"))))
      ),
      div(style = "color:#888; font-size:0.82em; margin-top:4px;",
        "Check Edit to adjust boundaries; + extends outward, - trims inward. Save to keep, Reset to discard."),
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
        value = (notes %||% NA_character_) %|NA|% "",
        placeholder = "Notes for this annotation unit"),
      actionButton("annotate_notes_clear", "Clear notes", class = "btn-sm"))
  )
}

#' Update / run modal: launch the Nextflow workflow with live progress
#'
#' Interactive (local browser): "Run from App" spawns Nextflow on this host.
#' Headless (served to a remote browser, e.g. an SSH-tunneled cluster login
#' node): running on the host would burn the login node, so offer an editable
#' cluster job script to Submit as job / Save instead.
#' @noRd
run_modal <- function(mode, n, headless = FALSE, sched = NULL, script = "") {
  wf <- if (mode == "Annotate") "Annotation (WF2)" else "Assembly (WF1)"
  cluster_section <- if (headless) {
    tagList(
      hr(style = "margin:10px 0;"),
      h5("Submit to cluster", style = "margin:6px 0 4px;"),
      div(class = "text-muted", style = "font-size:0.82em; margin-bottom:4px;",
        if (is.null(sched)) {
          "No scheduler (sbatch/qsub/bsub) found on PATH. Save the script and submit it yourself; running from the login node is discouraged."
        } else {
          sprintf("%s detected. Edit the job script if needed, then Submit as job (or Save script).", sched$name)
        }),
      tags$textarea(id = "submit_script_text", rows = 11,
        style = "width:100%; font-family:monospace; font-size:0.78em;", script))
  }
  run_progress <- if (!headless) {
    shinyjs::hidden(div(id = "nf_progress",
      hr(style = "margin:12px 0;"), h5("Progress:", style = "margin:0 0 4px;"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit; white-space:pre-wrap;",
          textOutput("nf_prog_all", container = tags$code)))))
  }
  footer <- if (headless) {
    tagList(
      if (!is.null(sched)) actionButton("submit_job", "Submit as job", class = "btn-success btn-sm"),
      actionButton("save_script", "Save script", class = "btn-sm"),
      actionButton("nf_close", "Close", class = "btn-sm"))
  } else {
    tagList(
      actionButton("nf_start", "Run from App", class = "btn-success btn-sm"),
      shinyjs::hidden(actionButton("nf_stop", "Stop / Interrupt", class = "btn-sm")),
      actionButton("nf_close", "Close", class = "btn-sm"))
  }
  modalDialog(
    title = div(style = "display:flex; justify-content:space-between; align-items:center; height:38px;",
      span(sprintf("%s — %d sample(s) to update", wf, n)),
      span(id = "rp_gears", class = "gears paused")),
    size = "l",
    shinyWidgets::prettyCheckbox("nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("nf_cmd_text")),
    run_progress,
    cluster_section,
    footer = footer
  )
}

#' Trim-assembly modal: cut the stitched assembly to the annotated feature span
#' or a custom range (shifts annotations + coverage).
#' @noRd
trim_modal <- function(n, a_from, a_to) {
  modalDialog(
    title = "Trim assembly", size = "l",
    tags$p(class = "text-muted", style = "font-size:0.85em;",
      sprintf("Current assembly length: %d bp. Annotated span: %d..%d bp.", n, a_from, a_to)),
    shinyWidgets::prettyRadioButtons("trim_mode", "Trim to:",
      choices = c("Annotated feature boundaries" = "annotated", "Custom range" = "custom"),
      selected = "annotated", status = "primary", inline = TRUE),
    tags$p(class = "text-muted", style = "font-size:0.8em; margin:2px 0;",
      "Drag across the plot to set the kept window (grey = trimmed); drag the box ",
      "edges to adjust. The From/To values track your selection."),
    plotOutput("trim_plot", height = "230px",
      brush = brushOpts(id = "trim_brush", direction = "x", resetOnNew = FALSE,
        fill = "#4C9BE8", opacity = 0.15)),
    conditionalPanel("input.trim_mode == 'custom'",
      div(style = "display:flex; gap:1em; margin-top:6px;",
        numericInput("trim_from", "From (bp)", value = a_from, min = 1, max = n),
        numericInput("trim_to", "To (bp)", value = a_to, min = 1, max = n))),
    opts_help("Trimming cuts the flanks off the stitched assembly and shifts every ",
      "feature coordinate + the coverage track to match. Features that fall entirely ",
      "outside the kept range are removed. Use \"Restore assembly\" to undo."),
    footer = tagList(
      actionButton("trim_apply", "Trim", icon = icon("scissors"), class = "btn-danger btn-sm"),
      # actionButton (not modalButton) so Cancel reopens the annotate details
      # modal instead of dropping back to the main Annotate table.
      actionButton("trim_cancel", "Cancel", class = "btn-sm"))
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

#' Export modal: pick which export group to bundle (or all locked samples), then
#' write it to the project's "out/export" folder. Without this, assigned groups
#' were ignored and every bundle collapsed into one undifferentiated "export".
#' @noRd
export_modal <- function(groups, n_locked) {
  choices <- c("All locked units (single bundle)" = "__all__",
    stats::setNames(groups, groups))
  modalDialog(
    title = "Export GenBank bundle", size = "m", easyClose = TRUE,
    if (length(groups) == 0) {
      tags$p(class = "text-muted", style = "font-size:0.9em;",
        sprintf("No export groups assigned. Exporting %d annotation-locked unit(s) as one bundle. Use \"Assign Group\" to batch submissions.", n_locked))
    } else {
      shinyWidgets::pickerInput("export_group_sel", "Export group:",
        choices = choices, selected = "__all__", width = "100%")
    },
    tags$p(class = "text-muted", style = "font-size:0.82em; margin-top:0.6em;",
      "Writes to the project's out/export folder: per-sample .fsa/.tbl, combined FASTA + feature table, per-marker fastas, summary CSV, and an HTML report."),
    footer = tagList(
      actionButton("export_go", "Export", icon = icon("download"), class = "btn-primary btn-sm"),
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

#' Read-coverage track for one assembly unit, with feature regions shaded
#' @noRd
plot_coverage <- function(id, path = 1L, scaffold = 1L) {
  if (is.null(id)) {
    return(NULL)
  }
  cov <- app_query("SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(cov) == 0 || is.na(cov$depths[1]) || !nzchar(cov$depths[1])) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", 0, 0, label = "no coverage data (run WF1)") +
      ggplot2::theme_void())
  }
  dep <- as.numeric(strsplit(cov$depths[1], ",")[[1]])
  df <- data.frame(pos = seq_along(dep), depth = dep)
  feats <- app_query("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
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

#' Coverage + annotation plot for the trim modal, with the kept window marked
#'
#' Reuses [plot_coverage()] and overlays the trim window: the flanks that would
#' be cut are greyed out and the two boundaries drawn as dashed lines.
#' @noRd
plot_trim <- function(id, path = 1L, scaffold = 1L, from = NA, to = NA) {
  g <- plot_coverage(id, path, scaffold)
  if (is.null(g)) {
    return(NULL)
  }
  if (!is.na(from) && !is.na(to) && to >= from) {
    g <- g +
      ggplot2::annotate("rect", xmin = -Inf, xmax = from, ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::annotate("rect", xmin = to, xmax = Inf, ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::geom_vline(xintercept = c(from, to), linetype = "dashed",
        colour = "#111111", linewidth = 0.7)
  }
  g
}

#' gggenes synteny map: the sample's annotations over its BLAST reference's
#' (when the reference has been fetched + annotated), aligned by element color.
#' @noRd
