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
#' project database, runs fastp + GetOrganelle, then scores and annotates each
#' assembly and computes the assembly-review aids (reference BLAST + coverage),
#' writing everything back to the database. Annotation and curation are handled
#' by re-running this workflow after you review and lock
#' the assemblies in the app.
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
  stopifnot(file.exists(file.path(path, ".sqlite")), file.exists(file.path(path, ".config")))
  check_nextflow_version("run_assembly", nextflow = nextflow)

  message("Running Nextflow assembly (WF1) - self-contained DB I/O...")
  args <- nf_command(path, resume = TRUE, work = work, ribopilot_src = ribopilot_src)
  nf_record_blastdb_launch(path)
  nf_cmd <- paste0(nxf_ver_prefix(nextflow), paste(c(shQuote(nextflow), shQuote(args)), collapse = " "))
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) {
    warning("Nextflow exited with status ", status)
  } else {
    message("WF1 complete (assemble -> score -> BLAST + coverage -> annotate + curate) ",
      "written to the database. Every scaffold is annotated as it is assembled; ",
      "review and lock assemblies in the app, then curate per unit in Annotation Review.")
  }
  invisible(status)
}
