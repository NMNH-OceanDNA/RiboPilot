#' RiboPilot app UI
#' @import shiny
#' @noRd
app_ui <- function(request) {
  bttn <- function(id, label, icon = NULL) {
    shinyWidgets::actionBttn(id, label = label, icon = if (!is.null(icon)) shiny::icon(icon) else NULL,
      style = "material-flat", size = "sm")
  }
  shiny::tagList(
    add_external_resources(),
    shiny::fluidPage(
      shiny::div(
        style = "display:flex; flex-flow:column; height:96vh;",
        # ---- control bar (mirrors MitoPilot) ----
        shiny::div(
          style = "display:flex; flex-flow:row nowrap; align-items:center; gap:0.8em; padding:0.5em 0;",
          shiny::tags$img(src = "www/favicon.ico", height = "26px"),
          shiny::strong("RiboPilot", style = "font-size:1.15em;"),
          shinyWidgets::pickerInput("mode", label = NULL, width = "185px",
            choices = c("Assemble", "Annotate & Curate", "Export")),
          bttn("refresh", NULL, "sync"),
          shiny::conditionalPanel(
            condition = "input.mode == 'Assemble'",
            shiny::div(style = "display:flex; gap:0.6em;",
              bttn("lock_btn", "Lock"), bttn("unlock_btn", "Unlock"),
              bttn("requeue_btn", "Set State: Ready"))
          ),
          shiny::conditionalPanel(
            condition = "input.mode == 'Annotate & Curate'",
            shiny::div(style = "display:flex; gap:0.6em;",
              bttn("lock_btn", "Lock"), bttn("unlock_btn", "Unlock"),
              bttn("reviewed_btn", "Mark Reviewed"))
          ),
          shiny::conditionalPanel(
            condition = "input.mode == 'Export'",
            shiny::div(style = "display:flex; align-items:center; gap:0.8em;",
              shiny::downloadButton("export_dl", "Export Selected", class = "btn-sm",
                icon = shiny::icon("file-export")),
              shiny::span(shiny::textOutput("export_sel_n", inline = TRUE),
                style = "color:#888; font-size:0.85em;"))
          ),
          shiny::span(shiny::textOutput("proj_label"),
            style = "margin-left:auto; color:#888; font-size:0.85em;")
        ),
        shiny::div(
          style = "flex:1 1 auto; min-height:0; overflow:hidden;",
          reactable::reactableOutput("tbl", height = "100%")
        )
      )
    )
  )
}
