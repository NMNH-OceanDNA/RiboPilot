# NMNH Hydra cluster session setup.
#
# An RStudio Server session on Hydra starts with a stripped PATH: no Univa Grid
# Engine, no Java, no ~/bin. Nextflow therefore cannot find `qsub`, and every
# task fails at submission with "Error submitting process ... for execution"
# while `qstat` stays empty and the job log looks like nothing happened.
# Mirrors MitoPilot's hydra_setup().

#' Directories an RStudio Server session on Hydra is missing from PATH
#' @noRd
HYDRA_BIN_DIRS <- c(
  "/share/apps/tools/java/21.0.2/bin",
  "/cm/shared/apps/uge/8.8.1/bin/lx-amd64"
)

#' Are we on the NMNH Hydra cluster?
#'
#' Greps `/etc/hosts` for "hydra", the same check MitoPilot uses.
#' @return `TRUE` if running on Hydra, otherwise `FALSE`.
#' @noRd
is_hydra_cluster <- function() {
  hosts <- tryCatch(readLines("/etc/hosts", warn = FALSE), error = function(e) NULL)
  !is.null(hosts) && any(grepl("hydra", hosts, ignore.case = TRUE))
}

#' Configure the R session environment for the NMNH Hydra cluster
#'
#' RStudio Server sessions on Hydra start with a stripped `PATH` that omits the
#' Univa Grid Engine, Java, and the user's `~/bin` directories, so Nextflow
#' cannot find `qsub` and every job submission fails. Call this once per session,
#' before launching the app, to put those directories back on `PATH`. Has no
#' effect (and warns) when not running on Hydra.
#'
#' @return Invisibly, `TRUE` if the Hydra environment was applied, otherwise
#'   `FALSE`.
#' @export
#' @examples
#' \dontrun{
#' RiboPilot::hydra_setup()
#' RiboPilot::RiboPilot(project = "/scratch/me/my_project")
#' }
hydra_setup <- function() {
  if (!is_hydra_cluster()) {
    warning("Not running on the NMNH Hydra cluster; environment unchanged.", call. = FALSE)
    return(invisible(FALSE))
  }
  dirs <- c(file.path(Sys.getenv("HOME"), "bin"), HYDRA_BIN_DIRS)
  cur <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  # force these to the front, dropping any lower-priority copies, so the UGE
  # qsub and Java 21 always win
  Sys.setenv(PATH = paste(c(dirs, cur[!cur %in% dirs]), collapse = .Platform$path.sep))

  check_nextflow_version("hydra_setup", on_too_old = "warn")
  pin <- nf_pin_version()
  if (!is.na(pin)) {
    Sys.setenv(NXF_VER = pin)
    message("Pinned NXF_VER=", pin, " for this session.")
  }
  sched <- detect_scheduler()
  message("Hydra environment applied. Scheduler: ",
    if (is.null(sched)) "still not found on PATH" else paste0(sched$cmd, " (", sched$name, ")"))
  invisible(TRUE)
}

#' Warn when a Hydra session has not been set up
#'
#' The failure this prevents is silent and expensive: Nextflow submits nothing,
#' retries three times per task, gives up under `errorStrategy 'ignore'`, and the
#' app reports an empty result instead of a missing `qsub`.
#' @noRd
warn_if_hydra_unconfigured <- function() {
  if (!is_hydra_cluster() || !is.null(detect_scheduler())) {
    return(invisible(FALSE))
  }
  message("This looks like NMNH Hydra, but no scheduler (qsub) is on PATH, so ",
    "Nextflow job submission will fail.\n  Run RiboPilot::hydra_setup() before ",
    "launching the app.")
  invisible(TRUE)
}
