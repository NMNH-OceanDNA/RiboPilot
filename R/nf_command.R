# Shared Nextflow command builder used by run_assembly()/run_annotate() (CLI)
# and the app's Update modal (processx). Returns the argument vector that
# follows the `nextflow` binary, so both a `bash -c "cd .. && nextflow <args>"`
# call and `processx::process$new("nextflow", args, wd = path)` can use it.

#' Build the Nextflow argument vector for a project workflow
#'
#' @param path Project directory (contains `.sqlite`, `.config`).
#' @param workflow "assemble" (WF1), "annotate" (WF2), or "fallback".
#' @param resume Add `-resume`.
#' @param work Nextflow work dir.
#' @param ribopilot_src Package source for R-script processes (default: loaded pkg).
#' @param method Fallback method ("bait" or "denovo"), passed as `--method`.
#' @return character vector of args (starting with "run").
#' @noRd
nf_command <- function(path, workflow = c("assemble", "annotate", "fallback"), resume = TRUE,
                       work = file.path(path, "work"), ribopilot_src = NULL, method = "bait") {
  workflow <- match.arg(workflow)
  main_nf <- app_sys("nextflow", switch(workflow,
    assemble = "main.nf", annotate = "main_annotate.nf", fallback = "main_fallback.nf"))
  cfg <- file.path(path, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg) else character(0)
  raw_dir <- config_param(cfg_lines, "rawDir") %||% ""
  min_depth <- config_param(cfg_lines, "minDepth") %||% "0"
  api_key <- config_param(cfg_lines, "ncbi_api_key") %||% ""
  if (length(api_key) == 0 || is.na(api_key) || api_key %in% c("NA", "NULL")) api_key <- ""
  src <- ribopilot_src %||% app_sys()
  taxdump <- file.path(path, ".taxonomy")

  c("run", main_nf, "-c", cfg, "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    "--db", file.path(path, ".sqlite"),
    "--projectDir", path,
    "--rawDir", raw_dir,
    "--minDepth", min_depth,
    "--ribopilot_src", src,
    "--taxdump_dir", taxdump,
    if (workflow == "fallback") c("--method", method),
    if (nzchar(api_key)) c("--ncbi_api_key", api_key))
}

#' Count samples eligible for a workflow (for the Update modal header)
#' @noRd
nf_pending_count <- function(path, workflow = c("assemble", "annotate")) {
  workflow <- match.arg(workflow)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  sql <- if (workflow == "assemble") {
    "SELECT COUNT(*) n FROM assemble WHERE assemble_switch IN (1,4) AND assemble_lock = 0"
  } else {
    "SELECT COUNT(*) n FROM annotate an JOIN assemble a ON a.ID = an.ID
       WHERE an.annotate_switch IN (1,4) AND an.annotate_lock = 0 AND a.assemble_lock = 1
       AND EXISTS (SELECT 1 FROM assemblies x WHERE x.ID = an.ID)"
  }
  DBI::dbGetQuery(con, sql)$n
}
