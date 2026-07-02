#' Open the RiboPilot GUI
#'
#' Launches the R-Shiny app for reviewing and curating a RiboPilot project.
#' Run from (or point `project` at) a project directory containing `.sqlite`.
#'
#' @param project Path to the project directory (default: current dir).
#' @param host,port,launch.browser Passed to Shiny (see [run_app()]).
#' @param ... Additional arguments passed to `run_app()`.
#' @export
RiboPilot <- function(project = getwd(), host = NULL, port = NULL, launch.browser = NULL, ...) {
  run_app(project = project, host = host, port = port, launch.browser = launch.browser, ...)
}

#' Run the RiboPilot Shiny application
#'
#' @param project Project directory containing `.sqlite`.
#' @param onStart,options,enableBookmarking,uiPattern Shiny app options.
#' @param host,port,launch.browser Shiny server options.
#' @param ... Passed to golem options.
#' @export
#' @importFrom shiny shinyApp
run_app <- function(project = getwd(), onStart = NULL, options = list(),
                    enableBookmarking = NULL, uiPattern = "/",
                    host = NULL, port = NULL, launch.browser = NULL, ...) {
  project <- normalizePath(project, mustWork = TRUE)
  if (!file.exists(file.path(project, ".sqlite"))) {
    stop("No RiboPilot project (.sqlite) found in: ", project, call. = FALSE)
  }
  if (is.null(launch.browser)) launch.browser <- interactive()
  options[["launch.browser"]] <- launch.browser
  if (!is.null(host)) options[["host"]] <- host
  if (!is.null(port)) options[["port"]] <- port

  golem::with_golem_options(
    app = shiny::shinyApp(
      ui = app_ui, server = app_server, onStart = onStart,
      options = options, enableBookmarking = enableBookmarking, uiPattern = uiPattern
    ),
    golem_opts = list(project = project, ...)
  )
}
