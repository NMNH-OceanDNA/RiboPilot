#' Lock + state status filter pickers (mirrors MitoPilot); values match the
#' assemble/annotate lock (0/1) and switch (1/2/3/4) codes.
#' @noRd
filter_pickers <- function(mode) {
  state_choices <- if (mode == "annotate") {
    c("Ready" = "1", "Done" = "2")
  } else {
    c("Ready" = "1", "Done" = "2", "Failed" = "3", "Needs BLAST" = "4")
  }
  picker_opts <- list(`actions-box` = TRUE, `select-all-text` = "All",
    `deselect-all-text` = "None", `selected-text-format` = "count > 0", width = "140px")
  shiny::tagList(
    shinyWidgets::pickerInput("lock_filter", label = "Lock:", width = "140px",
      choices = c("Unlocked" = "0", "Locked" = "1"), selected = c("0", "1"),
      multiple = TRUE, options = picker_opts),
    shinyWidgets::pickerInput("state_filter", label = "State:", width = "140px",
      choices = state_choices, selected = unname(state_choices),
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
          shinyWidgets::pickerInput("mode", label = "", width = 150,
            choices = c("Assemble", "Annotate", "Export")),
          bttn("refresh", NULL, "sync"),
          shiny::conditionalPanel("input.mode == 'Assemble'",
            ctrls(bttn("state_btn", "State"), bttn("lock_btn", "Lock"),
              bttn("run_btn", "Update", "circle-play"))),
          shiny::conditionalPanel("input.mode == 'Annotate'",
            ctrls(bttn("state_btn", "State"), bttn("lock_btn", "Lock"),
              bttn("reviewed_btn", "Mark Reviewed"), bttn("run_btn", "Update", "circle-play"))),
          shiny::conditionalPanel("input.mode == 'Export'",
            ctrls(bttn("assign_group_btn", "Assign Group"), bttn("clear_group_btn", "Clear Group"),
              shiny::downloadButton("export_dl", "Export Data", class = "btn-sm"),
              shiny::span(shiny::textOutput("export_sel_n", inline = TRUE),
                style = "color:#888; font-size:0.85em;"))),
          shiny::div(style = "margin-left:auto;",
            shinyWidgets::actionBttn("workdirs_btn", "Work Dirs",
              icon = shiny::icon("folder-tree"), style = "material-flat", size = "sm"))
        ),
        # ---- filter row (status pickers, below the control bar like MitoPilot) ----
        shiny::conditionalPanel("input.mode == 'Assemble'",
          shiny::div(style = "display:flex; flex-flow:row wrap; gap:1em; align-items:flex-end; margin:0.6em 0 0.7em;",
            filter_pickers("assemble"))),
        shiny::conditionalPanel("input.mode == 'Annotate'",
          shiny::div(style = "display:flex; flex-flow:row wrap; gap:1em; align-items:flex-end; margin:0.6em 0 0.7em;",
            filter_pickers("annotate"))),
        shiny::div(
          style = "flex:1 1 auto; min-height:0; overflow:hidden;",
          reactable::reactableOutput("tbl", height = "100%")
        )
      )
    )
  )
}
