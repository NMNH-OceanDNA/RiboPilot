#' RiboPilot app UI
#' @import shiny
#' @noRd
app_ui <- function(request) {
  shiny::tagList(
    add_external_resources(),
    shiny::fluidPage(
      shiny::div(
        style = "display:flex; flex-flow:row nowrap; align-items:center; gap:1em; margin:0.5em 0 1em;",
        shiny::tags$img(src = "www/favicon.ico", height = "28px"),
        shiny::h3("RiboPilot", style = "margin:0;"),
        shinyWidgets::pickerInput(
          inputId = "mode", label = NULL, width = "200px",
          choices = c("Assemble", "Annotate & Curate", "Export")
        ),
        shinyWidgets::actionBttn("refresh", label = NULL, icon = shiny::icon("sync"),
          style = "material-flat", size = "sm"),
        shiny::span(shiny::textOutput("proj_label"), style = "color:#888; font-size:0.85em;")
      ),
      reactable::reactableOutput("tbl")
    )
  )
}
