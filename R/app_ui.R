#' RiboPilot app UI
#' @import shiny
#' @noRd
app_ui <- function(request) {
  shiny::tagList(
    add_external_resources(),
    shiny::fluidPage(
      shiny::div(
        style = "display:flex; flex-flow:column;",
        # ---- control bar (mirrors MitoPilot) ----
        shiny::div(
          style = "display:flex; flex-flow:row nowrap; align-items:center; gap:1em; padding:0.5em 0;",
          shiny::tags$img(src = "www/favicon.ico", height = "26px"),
          shiny::strong("RiboPilot", style = "font-size:1.15em;"),
          shinyWidgets::pickerInput(
            inputId = "mode", label = NULL, width = "180px",
            choices = c("Assemble", "Annotate & Curate", "Export")
          ),
          shinyWidgets::actionBttn("refresh", label = NULL, icon = shiny::icon("sync"),
            style = "material-flat", size = "sm"),
          # Export-mode controls
          shiny::conditionalPanel(
            condition = "input.mode == 'Export'",
            shiny::div(
              style = "display:flex; flex-flow:row nowrap; align-items:center; gap:1em;",
              shiny::downloadButton("export_dl", "Export Selected",
                class = "btn-sm", icon = shiny::icon("file-export")),
              shiny::span(shiny::textOutput("export_sel_n", inline = TRUE),
                style = "color:#888; font-size:0.85em;")
            )
          ),
          shiny::span(shiny::textOutput("proj_label"), style = "margin-left:auto; color:#888; font-size:0.85em;")
        ),
        shiny::div(
          style = "padding:0.5em 0;",
          reactable::reactableOutput("tbl")
        )
      )
    )
  )
}
