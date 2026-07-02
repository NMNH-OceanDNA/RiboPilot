#' Lock + state status filter pickers (mirrors MitoPilot); values match the
#' assemble/annotate lock (0/1) and switch (1/2/3/4) codes.
#' @noRd
filter_pickers <- function(mode) {
  state_choices <- if (mode == "annotate") {
    c("Ready" = "1", "Done" = "2")
  } else {
    c("Ready" = "1", "Done" = "2", "Failed" = "3", "Needs BLAST" = "4")
  }
  shiny::tagList(
    shinyWidgets::pickerInput("lock_filter", label = NULL, width = "115px",
      choices = c("Unlocked" = "0", "Locked" = "1"), selected = c("0", "1"), multiple = TRUE,
      options = list(`actions-box` = TRUE, `selected-text-format` = "count > 2", title = "lock")),
    shinyWidgets::pickerInput("state_filter", label = NULL, width = "115px",
      choices = state_choices, selected = unname(state_choices), multiple = TRUE,
      options = list(`actions-box` = TRUE, `selected-text-format` = "count > 2", title = "state"))
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
          shinyWidgets::pickerInput("mode", label = NULL, width = "150px",
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
          shiny::span(shiny::textOutput("proj_label"),
            style = "margin-left:auto; color:#888; font-size:0.85em;")
        ),
        # ---- filter row (status pickers, below the control bar like MitoPilot) ----
        shiny::conditionalPanel("input.mode == 'Assemble'",
          shiny::div(style = "display:flex; gap:0.6em; padding:0 0 0.4em;", filter_pickers("assemble"))),
        shiny::conditionalPanel("input.mode == 'Annotate'",
          shiny::div(style = "display:flex; gap:0.6em; padding:0 0 0.4em;", filter_pickers("annotate"))),
        shiny::div(
          style = "flex:1 1 auto; min-height:0; overflow:hidden;",
          reactable::reactableOutput("tbl", height = "100%")
        )
      )
    )
  )
}
