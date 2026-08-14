# Shared Nextflow command builder used by run_assembly() (CLI) and the app's
# Update modal (processx). Returns the argument vector that follows the
# `nextflow` binary, so both a `bash -c "cd .. && nextflow <args>"` call and
# `processx::process$new("nextflow", args, wd = path)` can use it.
#
# There is ONE workflow. Re-annotation used to be WF2 (main_annotate.nf) and is
# now an arm of WF1: dev/design/15_WF2_DECOMMISSION.md.

#' Build the Nextflow argument vector for the project workflow
#'
#' @param path Project directory (contains `.sqlite`, `.config`).
#' @param resume Add `-resume`.
#' @param work Nextflow work dir.
#' @param ribopilot_src Package source for R-script processes (default: loaded pkg).
#' @return character vector of args (starting with "run").
#' @noRd
nf_command <- function(path, resume = TRUE,
                       work = file.path(path, "work"), ribopilot_src = NULL) {
  main_nf <- app_sys("nextflow", "main.nf")
  cfg <- file.path(path, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg) else character(0)
  # an absent key comes back as NA, and an unfilled config token lands as the
  # literal string "NA"; either would be passed to Nextflow verbatim
  cfg_val <- function(key, default) {
    v <- config_param(cfg_lines, key) %||% NA_character_ %|NA|% default
    if (v %in% c("NA", "NULL")) default else v
  }
  raw_dir <- cfg_val("rawDir", "")
  min_depth <- cfg_val("minDepth", "0")
  api_key <- cfg_val("ncbi_api_key", "")
  if (!nzchar(raw_dir)) {
    stop("No rawDir set in ", cfg, ".\n",
      "Without it the pipeline has nowhere to look for reads and every sample ",
      "fails. Set params.rawDir in the project config.", call. = FALSE)
  }
  src <- ribopilot_src %||% ribopilot_src_default()
  taxdump <- file.path(path, ".taxonomy")
  ref_dir <- nf_ref_dir(path)

  c("run", main_nf, "-c", cfg, "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    "--db", file.path(path, ".sqlite"),
    "--projectDir", path,
    "--rawDir", raw_dir,
    "--minDepth", min_depth,
    "--ribopilot_src", src,
    "--taxdump_dir", taxdump,
    "--ref_dir", ref_dir,
    if (nzchar(api_key)) c("--ncbi_api_key", api_key),
    # the config only mounts what project creation knew about; a database the
    # Reference tab wired later lives wherever the user put it, and the
    # assembler cannot read a staged symlink whose target is unmounted
    container_bind_opts(cfg_lines,
      reference_mount_dirs(path, c(path, raw_dir, config_mounted_dirs(cfg_lines)))))
}

#' The project's reference folder, resolved the same way `ref_config()` does
#'
#' A missing or unreadable project database must not fail command-building
#' the way a missing rawDir does; an empty string reaches Nextflow instead,
#' and priming/pushing degrade to a silent no-op from there.
#' @noRd
nf_ref_dir <- function(path) {
  tryCatch({
    con <- open_project_db(file.path(path, ".sqlite"))
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    ref_config(con)$ref_dir
  }, error = function(e) "")
}

#' Where the pipeline's R processes should load RiboPilot from
#'
#' `app_sys()` returns `<root>/inst` in a source checkout, which has no
#' `DESCRIPTION`, so passing it made `rp_env.R` skip its source branch and fall
#' back to the container's installed copy every time -- exactly the stale-install
#' shadowing that branch exists to prevent. Pass the package root instead when
#' one is there.
#' @noRd
ribopilot_src_default <- function() {
  s <- app_sys()
  root <- dirname(s)
  if (file.exists(file.path(root, "DESCRIPTION")) && dir.exists(file.path(root, "R"))) root else s
}

#' Count the work one run would do (for the Update modal header)
#'
#' Samples queued for assembly, plus units queued for annotation only, matching
#' what `query_samples` and `query_annotate_units` select. One button runs both,
#' so it is one number: a sample queued for assembly is counted once even though
#' its units will be annotated too, because assembly re-seeds them.
#' @noRd
nf_pending_count <- function(path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(path, ".sqlite"))
  on.exit(DBI::dbDisconnect(con))
  n_asm <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM assemble WHERE assemble_switch IN (1,4) AND assemble_lock = 0")$n
  n_ann <- DBI::dbGetQuery(con,
    "SELECT COUNT(*) n FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       JOIN assemble a ON a.ID = an.ID
      WHERE an.annotate_switch IN (1,4) AND an.annotate_lock = 0
        AND COALESCE(asm.ignore,0) = 0
        AND NOT (a.assemble_switch IN (1,4) AND COALESCE(a.assemble_lock,0) = 0)")$n
  as.integer(n_asm + n_ann)
}
