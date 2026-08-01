# Filter choices, shared by the pickers and the server-side visible-row mask so
# the two cannot drift (that drift is what made "Hold" a no-op).
LOCK_CHOICES <- c("Unlocked" = "0", "Locked" = "1")
STATE_CHOICES <- c("Hold" = "0", "Ready" = "1", "In progress" = "4",
  "Done" = "2", "Failed" = "3")

#' Lock + state status filter pickers (mirrors MitoPilot); values match the
#' assemble/annotate lock (0/1) and switch (0/1/2/3/4) codes. Each mode gets its
#' own picker ids (`<mode>_lock_filter` / `<mode>_state_filter`) so the Assemble
#' and Annotate filters do not share a DOM id and leak across modes.
#' @noRd
filter_pickers <- function(mode) {
  picker_opts <- list(`actions-box` = TRUE, `select-all-text` = "All",
    `deselect-all-text` = "None", `selected-text-format` = "count > 0", width = "140px")
  shiny::tagList(
    shinyWidgets::pickerInput(paste0(mode, "_lock_filter"), label = "Lock:", width = "140px",
      choices = LOCK_CHOICES, selected = unname(LOCK_CHOICES),
      multiple = TRUE, options = picker_opts),
    shinyWidgets::pickerInput(paste0(mode, "_state_filter"), label = "State:", width = "140px",
      choices = STATE_CHOICES, selected = unname(STATE_CHOICES),
      multiple = TRUE, options = picker_opts)
  )
}

#' RiboPilot app UI
#' @import shiny
#' @noRd
app_ui <- function(request) {
  bttn <- function(id, label, icon = NULL) {
    shinyWidgets::actionBttn(id, label = label, icon = if (!is.null(icon)) shiny::icon(icon) else NULL,
      style = "material-flat", size = "sm")
  }
  ctrls <- function(...) shiny::div(style = "display:flex; flex-flow:row nowrap; align-items:center; gap:1em;", ...)

  shiny::tagList(
    add_external_resources(),
    shiny::uiOutput("filter_css"),
    shiny::fluidPage(
      shiny::div(
        style = "display:flex; flex-flow:column; height:96vh;",
        # ---- control bar (mode picker + refresh + per-mode buttons) ----
        shiny::div(
          style = "display:flex; flex-flow:row nowrap; align-items:center; gap:1em; padding:0.5em 0;",
          # label = "" (not NULL) reserves label space above the select, which
          # counterbalances the form-group bottom margin so the picker centers
          # with the buttons (mirrors MitoPilot; no CSS override needed).
          # Named vector: the LABEL is what the user reads, the VALUE is what
          # every mode branch keys on. "Annotate" stays the value so the tab can
          # be renamed without touching the server.
          shinyWidgets::pickerInput("mode", label = "", width = 175,
            choices = c(Taxonomy = "Taxonomy", Reference = "Reference", Assemble = "Assemble",
              "Annotation Review" = "Annotate", Export = "Export")),
          bttn("refresh", NULL, "sync"),
          shiny::conditionalPanel("input.mode == 'Taxonomy'",
            ctrls(bttn("state_btn", "State"), bttn("lock_btn", "Lock"),
              bttn("edit_taxon_btn", "Edit Taxon"),
              bttn("run_btn", "Update", "circle-play"))),
          shiny::conditionalPanel("input.mode == 'Reference'",
            ctrls(bttn("ref_suggest_btn", "Suggest", "wand-magic-sparkles"),
              bttn("ref_edit_btn", "Edit Clade"), bttn("ref_existing_btn", "Use Existing"),
              bttn("ref_skip_btn", "Skip"), bttn("ref_build_btn", "Build", "circle-play"),
              bttn("ref_dir_btn", "Folder", "folder-open"))),
          shiny::conditionalPanel("input.mode == 'Assemble'",
            ctrls(bttn("state_btn", "State"), bttn("lock_btn", "Lock"),
              bttn("triage_btn", "Triage", "filter"),
              bttn("run_btn", "Update", "circle-play"))),
          # No Update button here: one workflow, launched from Assemble. This tab
          # is a review surface -- edit the feature model, trim, mark, lock.
          shiny::conditionalPanel("input.mode == 'Annotate'",
            ctrls(bttn("state_btn", "State"), bttn("lock_btn", "Lock"),
              bttn("reviewed_btn", "Reviewed"), bttn("problematic_btn", "Problem"),
              bttn("partial_btn", "Partial"))),
          shiny::conditionalPanel("input.mode == 'Export'",
            ctrls(bttn("assign_group_btn", "Assign Group"), bttn("clear_group_btn", "Clear Group"),
              bttn("export_btn", "Export Data", "download"),
              shiny::span(shiny::textOutput("export_sel_n", inline = TRUE),
                style = "color:#888; font-size:0.85em;"))),
          shiny::div(style = "margin-left:auto; display:flex; gap:0.9em; align-items:center;",
            shiny::span(shiny::textOutput("n_selected", inline = TRUE),
              style = "color:#888; font-size:0.85em;"),
            shiny::downloadButton("table_csv", "CSV", class = "btn-sm"),
            shinyWidgets::actionBttn("workdirs_btn", "Work Dirs",
              icon = shiny::icon("folder-tree"), style = "material-flat", size = "sm"))
        ),
        # ---- filter row (status pickers, below the control bar like MitoPilot) ----
        shiny::conditionalPanel("input.mode == 'Assemble'",
          shiny::div(style = "display:flex; flex-flow:row wrap; gap:1em; align-items:flex-end; margin:0.6em 0 0.7em;",
            filter_pickers("Assemble"),
            shinyWidgets::pickerInput("col_groups", "Columns:", width = "170px",
              choices = c("Stats", "Options"),
              selected = c("Stats", "Options"), multiple = TRUE,
              options = list(`actions-box` = TRUE, `select-all-text` = "All",
                `deselect-all-text` = "None", `selected-text-format` = "count > 0")))),
        shiny::conditionalPanel("input.mode == 'Annotate'",
          shiny::div(style = "display:flex; flex-flow:row wrap; gap:1em; align-items:flex-end; margin:0.6em 0 0.7em;",
            filter_pickers("Annotate"))),
        shiny::div(
          style = "flex:1 1 auto; min-height:0; overflow:hidden;",
          reactable::reactableOutput("tbl", height = "100%")
        )
      )
    )
  )
}
