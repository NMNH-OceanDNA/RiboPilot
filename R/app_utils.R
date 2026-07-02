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

#' Enumerate Nextflow task work directories for one sample
#'
#' Parses the project's `.nextflow.log*` files (the app launches Nextflow from
#' the project dir) for completed-task records, which Nextflow writes for every
#' attempt (success or failure), each with its own work dir. Mirrors MitoPilot's
#' work-dir browser.
#' @noRd
find_workdirs <- function(project_dir, sample_id) {
  empty <- data.frame(process = character(0), status = character(0),
    modified = character(0), workdir = character(0), exists = logical(0),
    stringsAsFactors = FALSE)
  if (is.null(project_dir) || is.null(sample_id) || !nzchar(sample_id)) return(empty)
  logs <- list.files(project_dir, pattern = "^\\.nextflow\\.log", full.names = TRUE, all.files = TRUE)
  if (length(logs) == 0) return(empty)
  lines <- unlist(lapply(logs, function(f) tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))))
  th <- grep("Task completed > TaskHandler\\[", lines, value = TRUE)
  if (length(th) == 0) return(empty)

  grab1 <- function(x, pat) vapply(regmatches(x, regexec(pat, x)),
    function(v) if (length(v) >= 2) v[2] else NA_character_, character(1))
  name <- grab1(th, "name:\\s*(.*?);")
  exit <- grab1(th, "exit:\\s*(.*?);")
  workdir <- grab1(th, "workDir:\\s*([^]\\s]+)")
  pt <- regmatches(name, regexec("([^:[:space:]]+)\\s*\\(([^)]+)\\)\\s*$", name))
  process <- vapply(pt, function(v) if (length(v) >= 3) v[2] else NA_character_, character(1))
  sample <- sub("\\.[0-9]+$", "", vapply(pt, function(v) if (length(v) >= 3) v[3] else NA_character_, character(1)))

  keep <- !is.na(workdir) & !is.na(process) & sample == sample_id
  if (!any(keep)) return(empty)
  out <- data.frame(process = process[keep], exit = exit[keep],
    workdir = trimws(workdir[keep]), stringsAsFactors = FALSE)
  out <- out[!duplicated(out$workdir), , drop = FALSE]
  out$status <- ifelse(out$exit == "0", "success", "failed")
  out$exists <- dir.exists(out$workdir)
  mt <- file.info(out$workdir)$mtime
  out$modified <- ifelse(is.na(mt), "-", format(mt, "%Y-%m-%d %H:%M"))
  out[, c("process", "status", "modified", "workdir", "exists")]
}

#' Open a path in the OS file manager (best-effort)
#' @noRd
open_path <- function(path) {
  if (is.null(path) || !nzchar(path)) return(invisible())
  cmd <- if (Sys.info()[["sysname"]] == "Darwin") "open" else
    if (.Platform$OS.type == "windows") "explorer" else "xdg-open"
  try(system2(cmd, shQuote(path), stdout = FALSE, stderr = FALSE, wait = FALSE), silent = TRUE)
  invisible()
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
