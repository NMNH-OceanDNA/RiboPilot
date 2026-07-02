# Annotation driver (WF2), self-contained: re-annotate samples whose selected
# path changed. Like WF1, the Nextflow pipeline reads/writes the DB itself
# (sqlite3 + R-script processes); this driver just launches it.
#
# Queue a sample for (re-)annotation by setting annotate.annotate_switch = 1
# (done automatically when the selected path is changed in the app).

#' Run the annotation workflow (WF2) for a project
#'
#' @param path Project directory.
#' @param nextflow Path to the nextflow binary.
#' @param work Nextflow work dir (default `<path>/work`).
#' @param conda_env Optional conda env providing the tools (see [run_assembly()]).
#' @param ribopilot_src Package source for the R-script processes.
#' @return (invisibly) the Nextflow exit status.
#' @export
run_annotate <- function(path = ".", nextflow = "nextflow", work = file.path(path, "work"),
                         conda_env = NULL, ribopilot_src = NULL) {
  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }
  path <- normalizePath(path)
  stopifnot(file.exists(file.path(path, ".sqlite")), file.exists(file.path(path, ".config")))

  if (nf_pending_count(path, "annotate") == 0) {
    message("No samples queued for annotation (set annotate_switch = 1 to (re-)annotate).")
    return(invisible(0L))
  }
  message("Running Nextflow annotation (WF2) — self-contained DB I/O...")
  args <- nf_command(path, "annotate", resume = TRUE, work = work, ribopilot_src = ribopilot_src)
  nf_cmd <- paste(c(shQuote(nextflow), shQuote(args)), collapse = " ")
  status <- system2("bash", c("-c", shQuote(paste0("cd ", shQuote(path), " && ", nf_cmd))), wait = TRUE)
  if (status != 0) {
    warning("Nextflow (annotate) exited with status ", status)
  } else {
    message("WF2 complete (stitch -> annotate -> reconcile) written to the database.")
  }
  invisible(status)
}
