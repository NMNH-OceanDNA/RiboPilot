# Shiny app helpers.

#' Attach external resources (CSS/JS/favicon), mirroring MitoPilot styling
#' @noRd
add_external_resources <- function() {
  golem::add_resource_path("www", app_sys("app", "www"))
  shiny::tags$head(
    golem::favicon(),
    golem::bundle_resources(path = app_sys("app", "www"), app_title = "RiboPilot"),
    shinyjs::useShinyjs()
  )
}

#' Active project directory (from golem options, else cwd)
#' @noRd
app_project <- function() {
  golem::get_golem_options("project") %||% getwd()
}

#' Open a connection to the project database
#' @noRd
app_con <- function() {
  DBI::dbConnect(RSQLite::SQLite(), file.path(app_project(), ".sqlite"))
}

#' Query the project database and disconnect
#' @noRd
app_query <- function(sql, params = NULL) {
  con <- app_con()
  on.exit(DBI::dbDisconnect(con))
  if (is.null(params)) DBI::dbGetQuery(con, sql) else DBI::dbGetQuery(con, sql, params = params)
}
