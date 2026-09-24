#' RiboPilot user config directory
#'
#' Persistent, per-user directory where saved cluster profiles
#' (`config.<name>` files) are stored. Created on demand by
#' [generate_config()]. Uses [tools::R_user_dir()] so it follows the platform
#' convention.
#'
#' @return Path to the RiboPilot config directory (not guaranteed to exist).
#' @noRd
ribopilot_config_dir <- function() {
  tools::R_user_dir("RiboPilot", which = "config")
}

#' Build a Nextflow container-engine block
#'
#' The RiboPilot pipeline opens the project SQLite database directly by path, so
#' the container must see the project directory. Every engine block bind-mounts
#' `${params.projectDir}` (resolved by Nextflow at run time from `--projectDir`);
#' new_project() adds the read, assembly, and reference folders through
#' `<<EXTRA_MOUNTS>>`.
#'
#' @param engine One of "docker", "singularity", "apptainer".
#' @param cache Optional cacheDir (singularity/apptainer only).
#' @param run_options Optional extra runOptions appended to the bind mounts.
#' @noRd
container_engine_block <- function(engine, cache = NULL, run_options = NULL) {
  extra <- if (!is.null(run_options) && nzchar(run_options)) run_options else ""
  # <<EXTRA_MOUNTS>> is filled by new_project() with the read, assembly,
  # and local reference folders (empty otherwise).
  if (engine == "docker") {
    # -u / fixOwnership stop SQLite leaving root-owned .sqlite-wal / .sqlite-shm
    # files in the project directory, which the Shiny app then cannot open.
    binds <- paste0(
      "-u \\$(id -u):\\$(id -g) ",
      "-v ${params.projectDir}:${params.projectDir}<<EXTRA_MOUNTS>>")
    ro <- trimws(paste(binds, extra))
    return(paste0("docker {\n  enabled = true\n  fixOwnership = true\n  runOptions = \"", ro, "\"\n}"))
  }
  binds <- "-B ${params.projectDir}<<EXTRA_MOUNTS>>"
  ro <- trimws(paste(binds, extra))
  lines <- c(paste0(engine, " {"), "  enabled = true", "  autoMounts = true", "  pullTimeout = '60 min'")
  if (!is.null(cache) && nzchar(cache)) {
    lines <- c(lines, paste0("  cacheDir = '", cache, "'"))
  }
  lines <- c(lines, paste0("  runOptions = \"", ro, "\""))
  lines <- c(lines, "}")
  paste(lines, collapse = "\n")
}

#' Substitute `<<PLACEHOLDER>>` tokens in config template lines
#'
#' `NULL` values are skipped, leaving the token intact for a later pass.
#'
#' @param lines Character vector of config file lines.
#' @param replacements Named list; name `FOO` replaces token `<<FOO>>`.
#' @noRd
#' Config tokens filled later, by new_project() rather than generate_config()
#' @noRd
PROJECT_CONFIG_TOKENS <- c("CONTAINER_ID", "RAW_DIR", "ASMB_DIR", "MIN_DEPTH",
                           "GENETIC_CODE", "NCBI_API_KEY", "EXTRA_MOUNTS",
                           "HOST_CPUS", "HOST_MEMORY")

#' Placeholder tokens still present in config lines
#' @noRd
config_tokens <- function(lines) {
  unique(unlist(regmatches(lines, gregexpr("<<[A-Z_]+>>", lines))))
}

fill_config <- function(lines, replacements) {
  for (nm in names(replacements)) {
    val <- replacements[[nm]]
    if (is.null(val)) next
    token <- paste0("<<", nm, ">>")
    lines <- stringr::str_replace_all(lines, stringr::fixed(token), as.character(val))
  }
  lines
}

#' Resolve an executor name to a config template path
#'
#' Resolution order: a saved user profile (`<profile_dir>/config.<executor>`)
#' first, then a package built-in (`inst/config.<executor>`).
#'
#' @param executor Executor / profile name.
#' @param profile_dir User profile directory (see [ribopilot_config_dir()]).
#' @noRd
resolve_config <- function(executor, profile_dir = ribopilot_config_dir()) {
  user_path <- file.path(profile_dir, paste0("config.", executor))
  if (file.exists(user_path)) {
    return(user_path)
  }
  pkg_path <- app_sys(paste0("config.", executor))
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }
  available <- list_configs(profile_dir = profile_dir)
  stop(
    "No config found for executor '", executor, "'.\nAvailable: ",
    paste(available$name, collapse = ", "),
    call. = FALSE
  )
}

#' List available executor configs
#'
#' Reports the package built-in templates plus any saved cluster profiles. The
#' `name` column is what you pass to `new_project(executor = ...)`.
#'
#' @param profile_dir User profile directory (see [ribopilot_config_dir()]).
#'
#' @return A data.frame with columns `name`, `type` ("builtin" or "saved"),
#'   and `path`.
#' @noRd
list_configs <- function(profile_dir = ribopilot_config_dir()) {
  builtin_files <- Sys.glob(file.path(app_sys(), "config.*"))
  saved_files <- if (dir.exists(profile_dir)) {
    Sys.glob(file.path(profile_dir, "config.*"))
  } else {
    character(0)
  }
  to_df <- function(files, type) {
    if (length(files) == 0) {
      return(NULL)
    }
    data.frame(
      name = sub("^config\\.", "", basename(files)),
      type = type,
      path = files,
      stringsAsFactors = FALSE
    )
  }
  out <- rbind(to_df(builtin_files, "builtin"), to_df(saved_files, "saved"))
  if (is.null(out)) {
    out <- data.frame(
      name = character(0), type = character(0), path = character(0),
      stringsAsFactors = FALSE
    )
  }
  out
}

#' Generate and save a reusable cluster config profile
#'
#' Configure your HPC cluster once and reuse it for every project. Mirrors
#' MitoPilot's `generate_config()`. Builds a Nextflow config from a generic
#' scheduler template, fills in cluster-level settings, and saves it as a named
#' profile in [ribopilot_config_dir()]. Afterwards
#' `new_project(executor = "<name>")` will find and use it automatically.
#'
#' @param name Profile name. Saved as `config.<name>`; pass as `executor`.
#' @param scheduler Base template: "slurm", "sge", "pbs", "lsf", or "local"
#'   (Docker on one machine, the same as executor = "local").
#' @param container_engine "auto" (default) picks singularity for the schedulers
#'   and docker for "local". "singularity", "apptainer", and "docker" apply to
#'   slurm, sge, pbs, and lsf; "local" always uses Docker.
#' @param container_cache Folder where Singularity or Apptainer keeps the pulled
#'   image, shared by every project that uses this profile; it must be visible
#'   from the compute nodes. Default: `NXF_SINGULARITY_CACHEDIR` (or
#'   `NXF_APPTAINER_CACHEDIR`) when set, otherwise a "containers" folder in the
#'   reference folder (`RIBOPILOT_REF_DIR`, see ribopilot_ref_dir()). Pass "" to
#'   leave it to Nextflow (the project's work/singularity folder).
#' @param container_run_options Optional extra runOptions for the container engine (any engine), appended to the bind mounts.
#' @param queue Partition / queue name. If `NULL`, the directive is omitted.
#' @param walltime Per-job walltime, in a form your scheduler accepts (e.g.
#'   "24h", "48:00:00"). Many SLURM and PBS sites reject a job that requests no
#'   walltime instead of applying a default. If `NULL`, the directive is
#'   omitted and the site default applies.
#' @param account Optional accounting string, folded into clusterOptions.
#' @param cluster_options Optional raw clusterOptions string appended verbatim.
#' @param penv SGE parallel environment (default "mthread").
#' @param profile_dir Directory to save the profile in.
#' @param interactive If `TRUE` and interactive, prompt for unset values.
#' @param overwrite Overwrite an existing profile (default `FALSE`).
#'
#' @return (invisibly) the path to the saved profile.
#' @export
generate_config <- function(
    name,
    scheduler = c("slurm", "sge", "pbs", "lsf", "local"),
    container_engine = c("auto", "singularity", "apptainer", "docker"),
    container_cache = NULL,
    container_run_options = NULL,
    queue = NULL,
    walltime = NULL,
    account = NULL,
    cluster_options = NULL,
    penv = "mthread",
    profile_dir = ribopilot_config_dir(),
    interactive = FALSE,
    overwrite = FALSE) {
  if (missing(name) || !is.character(name) || length(name) != 1 || !nzchar(name)) {
    stop("`name` must be a single non-empty string.", call. = FALSE)
  }
  scheduler <- match.arg(scheduler)
  container_engine <- match.arg(container_engine)
  if (scheduler == "local" && container_engine %nin% c("auto", "docker")) {
    stop("A \"local\" profile runs the Docker image. container_engine = \"", container_engine,
      "\" works with the slurm, sge, pbs, and lsf schedulers.", call. = FALSE)
  }

  if (container_engine == "auto") {
    container_engine <- if (scheduler == "local") "docker" else "singularity"
  }
  if (container_engine != "docker" && is.null(container_cache)) {
    env <- Sys.getenv(c("NXF_SINGULARITY_CACHEDIR", "NXF_APPTAINER_CACHEDIR"))
    container_cache <- c(env[nzchar(env)], file.path(ribopilot_ref_dir(), "containers"))[[1]]
  }

  if (isTRUE(interactive) && base::interactive()) {
    prompt_if_null <- function(val, msg) {
      if (!is.null(val)) {
        return(val)
      }
      ans <- trimws(readline(msg))
      if (nzchar(ans)) ans else NULL
    }
    if (scheduler != "local") {
      queue <- prompt_if_null(queue, "Queue / partition (blank = cluster default): ")
      walltime <- prompt_if_null(walltime, "Walltime per job, e.g. 24h (blank = site default): ")
      account <- prompt_if_null(account, "Account / project (blank = none): ")
      cluster_options <- prompt_if_null(cluster_options, "Extra clusterOptions (blank = none): ")
    }
  }

  acct_flag <- ""
  if (!is.null(account) && nzchar(account)) {
    acct_flag <- switch(scheduler,
      slurm = paste0("--account=", account),
      pbs = paste0("-A ", account),
      lsf = paste0("-P ", account),
      sge = paste0("-P ", account),
      ""
    )
  }
  cluster_options <- trimws(paste(acct_flag, cluster_options %||% ""))

  template <- app_sys(paste0("config.", scheduler))
  if (!nzchar(template) || !file.exists(template)) {
    stop("Base template not found for scheduler '", scheduler, "'.", call. = FALSE)
  }
  lines <- readLines(template)

  if (is.null(queue) || !nzchar(queue)) {
    lines <- lines[!grepl("<<QUEUE>>", lines, fixed = TRUE)]
  }
  if (is.null(walltime) || !nzchar(walltime)) {
    lines <- lines[!grepl("<<WALLTIME>>", lines, fixed = TRUE)]
  }

  lines <- fill_config(lines, list(
    CONTAINER_ENGINE = container_engine_block(container_engine, container_cache, container_run_options),
    QUEUE = queue,
    WALLTIME = walltime,
    CLUSTER_OPTIONS = cluster_options,
    PENV = if (scheduler == "sge") penv else NULL
  ))

  # a token neither round fills (sge with penv = NULL) would stop
  # new_project() later; say so now
  manual <- setdiff(config_tokens(lines), paste0("<<", PROJECT_CONFIG_TOKENS, ">>"))
  if (length(manual)) {
    warning(glue::glue(
      "Profile '{name}' still needs these filled in by hand before new_project(): ",
      "{paste(manual, collapse = ', ')}."), call. = FALSE)
  }

  if (!dir.exists(profile_dir)) {
    dir.create(profile_dir, recursive = TRUE)
  }
  out_path <- file.path(profile_dir, paste0("config.", name))
  if (file.exists(out_path) && !overwrite) {
    stop(
      "Profile '", name, "' already exists at ", out_path,
      ". Use overwrite = TRUE to replace it.",
      call. = FALSE
    )
  }
  writeLines(lines, out_path)

  message("Saved cluster profile '", name, "' to ", out_path)
  if (container_engine != "docker" && nzchar(container_cache)) {
    message("Container images are cached in ", container_cache,
      " (shared by every project that uses this profile).")
  }
  message("Use it with: new_project(..., executor = '", name, "')")
  invisible(out_path)
}
