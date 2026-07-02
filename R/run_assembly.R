# Assembly driver (WF1). Self-contained Nextflow: the pipeline reads samples
# from the project DB (sqlite3), assembles, and writes results back (R-script +
# tool processes). This driver just launches Nextflow with the right params; no
# R samplesheet prep and no R load.

#' Resolve a conda env (name or path) to its `bin/` directory
#' @noRd
resolve_conda_bin <- function(env) {
  if (dir.exists(env) && basename(env) == "bin") {
    return(normalizePath(env))
  }
  if (dir.exists(file.path(env, "bin"))) {
    return(normalizePath(file.path(env, "bin")))
  }
  roots <- c(
    file.path(dirname(dirname(Sys.getenv("CONDA_PREFIX", ""))), "envs"),
    file.path(path.expand("~"), c("miniconda3", "miniforge3", "mambaforge", "anaconda3"), "envs")
  )
  for (r in roots) {
    cand <- file.path(r, env, "bin")
    if (dir.exists(cand)) {
      return(normalizePath(cand))
    }
  }
  stop("Could not locate conda env '", env, "'. Pass a full path to the env or its bin/.", call. = FALSE)
}

#' Read a filled `params.X = '...'` value from the project .config
#' @noRd
config_param <- function(config_lines, key) {
  hit <- grep(paste0("^\\s*", key, "\\s*="), config_lines, value = TRUE)
  if (length(hit) == 0) {
    return(NA_character_)
  }
  val <- sub(paste0("^\\s*", key, "\\s*=\\s*"), "", hit[1])
  val <- sub("//.*$", "", val)
  trimws(gsub("['\"]", "", val))
}

#' Run the assembly workflow (WF1) for a project
#'
#' Launches the self-contained Nextflow pipeline: it reads samples from the
#' project database, runs fastp + GetOrganelle, then scores, stitches, and
#' annotates each sample, writing everything back to the database.
#'
#' @param path Project directory (contains `.sqlite` and `.config`).
#' @param nextflow Path to the `nextflow` binary.
#' @param work Nextflow work directory (default `<path>/work`).
#' @param conda_env Optional conda env providing the pipeline tools; its `bin/`
#'   is prepended to PATH for the run.
#' @param ribopilot_src Path to the RiboPilot package source for the R-script
#'   processes (only needed when the package is not installed; defaults to the
#'   loaded package location).
#' @return (invisibly) the Nextflow exit status.
#' @export
run_assembly <- function(path = ".", nextflow = "nextflow", work = file.path(path, "work"),
                         conda_env = NULL, ribopilot_src = NULL) {
  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  path <- normalizePath(path)
  db <- file.path(path, ".sqlite")
  cfg <- file.path(path, ".config")
  main_nf <- app_sys("nextflow", "main.nf")
  stopifnot(file.exists(db), file.exists(cfg), file.exists(main_nf))

  cfg_lines <- readLines(cfg)
  raw_dir <- config_param(cfg_lines, "rawDir")
  min_depth <- config_param(cfg_lines, "minDepth") %||% "0"
  src <- ribopilot_src %||% app_sys() # package root (source under load_all)
  taxdump <- file.path(path, ".taxonomy")

  message("Running Nextflow assembly (WF1) — self-contained DB I/O...")
  nf_cmd <- paste(c(
    shQuote(nextflow), "run", shQuote(main_nf), "-c", shQuote(cfg),
    "-work-dir", shQuote(work), "-resume",
    "--db", shQuote(db),
    "--projectDir", shQuote(path),
    "--rawDir", shQuote(raw_dir),
    "--minDepth", min_depth,
    "--ribopilot_src", shQuote(src),
    "--taxdump_dir", shQuote(taxdump)
  ), collapse = " ")
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) {
    warning("Nextflow exited with status ", status)
  } else {
    message("WF1 complete (assemble -> score -> stitch -> annotate) written to the database.")
  }
  invisible(status)
}
