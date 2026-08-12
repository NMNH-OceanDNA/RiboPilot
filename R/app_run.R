#' Is the app running headless (served to a remote browser, e.g. a tunneled
#' cluster login node)?
#' @noRd
app_headless <- function() isTRUE(getOption("RiboPilot.headless", FALSE))

#' Should the app treat this session as a remote/cluster one?
#'
#' Headless decides whether the run modal offers a cluster job script or a local
#' "Run from App" button. It used to be `!interactive()`, which is wrong exactly
#' where it matters most: an RStudio Server session on a login node IS
#' interactive, so the app called itself a desktop and hid the submit-to-cluster
#' UI from the only users who need it. RStudio Server and the presence of a batch
#' scheduler are the two signals that this R session is on a cluster.
#' @param launch_browser Whether Shiny will open a local browser.
#' @noRd
detect_headless <- function(launch_browser) {
  if (!isTRUE(launch_browser)) {
    return(TRUE)
  }
  identical(Sys.getenv("RSTUDIO_PROGRAM_MODE"), "server") || !is.null(detect_scheduler())
}

#' Detect a batch-scheduler submit command on PATH (sbatch/qsub/bsub).
#' @return list(cmd, name) or NULL.
#' @noRd
detect_scheduler <- function() {
  for (s in list(c("sbatch", "Slurm"), c("qsub", "SGE/PBS"), c("bsub", "LSF"))) {
    if (nzchar(unname(Sys.which(s[1])))) return(list(cmd = s[1], name = s[2]))
  }
  NULL
}

#' Project subdirectories holding cluster job files
#'
#' Every submission drops three files (a script, the scheduler's job log, and
#' Nextflow's own log). In the project root they bury the handful of files a user
#' actually opens, so each kind gets its own directory. Visible, not dotted: the
#' logs are what you read when a run fails.
#' @noRd
JOB_DIRS <- list(script = "jobs", log = "logs", nextflow = "nextflow_logs")

#' Path inside one of the job directories, creating the directory on demand
#'
#' The scheduler will not create a missing directory for its `-o` file: it fails
#' the job instead, before anything has run and with nowhere to say why.
#' @noRd
job_path <- function(project, which, file) {
  d <- file.path(project, JOB_DIRS[[which]])
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.path(d, file)
}

#' Timestamped base name shared by one submission's script and its logs
#'
#' One base per submission, so a run's script, scheduler log, and Nextflow log
#' sort together and carry the same name. The fixed names this replaces meant
#' every run overwrote the last one's script and Nextflow log, which is exactly
#' the evidence a failed run is diagnosed from. Mirrors MitoPilot.
#' @noRd
job_basename <- function(kind, time = Sys.time()) {
  paste0(kind, "_", format(time, "%Y-%m-%d_%H-%M-%S"))
}

#' Wrap a Nextflow command in a scheduler job script
#'
#' A multi-hour run belongs in the queue, not on the login node. The script is a
#' starting point the user can edit in the run modal before saving/submitting.
#' @noRd
#' Closing lines of every job script
#'
#' With `-ansi-log false` Nextflow's completion summary never reaches the job's
#' output file: the log just stops after the last "Submitted process" line, so a
#' run that succeeded reads exactly like one that was killed. Say how it ended,
#' and hand the status back to the scheduler.
#' @noRd
NF_EXIT_TRAILER <- c(
  "status=$?",
  'echo "[RiboPilot] nextflow exited with status $status"',
  "exit $status")

job_script <- function(project, cmd, sched, jobname) {
  if (is_hydra_cluster()) {
    return(hydra_job_script(project, cmd, jobname))
  }
  # log paths are absolute: the app submits from wherever its R session happens
  # to be, which is not the project, so a relative -o lands somewhere arbitrary
  log <- function(suffix) job_path(project, "log", paste0(jobname, suffix))
  hdr <- if (is.null(sched)) {
    character(0)
  } else {
    switch(sched$cmd,
      sbatch = c(sprintf("#SBATCH -J %s", jobname), "#SBATCH -c 2",
        "#SBATCH --mem=8G", "#SBATCH -t 48:00:00",
        sprintf("#SBATCH -o %s", log("_%%j.log"))),
      qsub = c(sprintf("#$ -N %s", jobname), "#$ -pe mthread 2",
        "#$ -l h_vmem=8G", "#$ -j y", "#$ -S /bin/bash",
        sprintf("#$ -o %s", log(".$JOB_ID.log"))),
      bsub = c(sprintf("#BSUB -J %s", jobname), "#BSUB -n 2",
        "#BSUB -M 8000", sprintf("#BSUB -o %s", log(".%%J.log"))),
      character(0))
  }
  pin <- nf_pin_version()
  pin_line <- if (is.na(pin)) character(0) else paste0("export NXF_VER=", pin)
  paste(c("#!/bin/bash", hdr, "", pin_line, sprintf("cd %s", shQuote(project)), cmd,
    NF_EXIT_TRAILER, ""), collapse = "\n")
}

#' Give the run its own Nextflow log
#'
#' `-log` is a top-level option, so it goes BEFORE `run`. Without it Nextflow
#' writes .nextflow.log into the launch directory and rotates it, so the log of
#' the run you want is gone a couple of runs later.
#' @noRd
nf_log_opt <- function(project, base) {
  c("-log", job_path(project, "nextflow", paste0(base, ".nextflow.log")))
}

#' Cluster job script for the main workflow
#' @param base Timestamped name for this submission's script and logs.
#' @noRd
build_submit_script <- function(project, resume = TRUE, sched = detect_scheduler(),
                                base = job_basename("ribopilot")) {
  # one arg (the container bind options) contains spaces, so quote per argument
  args <- c(nf_log_opt(project, base), nf_command(project, resume = resume),
    "-ansi-log", "false")
  cmd <- paste(c(shQuote("nextflow"), shQuote(args)), collapse = " ")
  job_script(project, cmd, sched, base)
}

#' Cluster job script for a reference BUILD run
#'
#' The manifest must already be on disk: the job reads it when it runs, which is
#' long after the modal that wrote the script has closed.
#' @param base Timestamped name for this submission's script and logs.
#' @noRd
build_reference_submit_script <- function(project, manifest, ref_dir, api_key = "",
                                          refseq = FALSE, resume = TRUE,
                                          sched = detect_scheduler(),
                                          base = job_basename("ribopilot_ref")) {
  args <- nf_reference_command(project, "BUILD", manifest, ref_dir, ref_dir,
    api_key = api_key, refseq = refseq, resume = resume)
  args <- c(nf_log_opt(project, base), args, "-ansi-log", "false")
  job_script(project, paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "),
    sched, base)
}

#' Open the RiboPilot GUI
#'
#' Launches the R-Shiny app for reviewing and curating a RiboPilot project.
#' Run from (or point `project` at) a project directory containing `.sqlite`.
#'
#' @param project Path to the project directory (default: current dir).
#' @param host,port,launch.browser Passed to Shiny (see [run_app()]).
#' @param headless Offer a cluster job script instead of running Nextflow where
#'   the app is served. Detected automatically (RStudio Server or a batch
#'   scheduler on PATH means yes); pass `TRUE`/`FALSE` to override.
#' @param ... Additional arguments passed to `run_app()`.
#' @export
RiboPilot <- function(project = getwd(), host = NULL, port = NULL, launch.browser = NULL,
                      headless = NULL, ...) {
  run_app(project = project, host = host, port = port, launch.browser = launch.browser,
    headless = headless, ...)
}

#' Run the RiboPilot Shiny application
#'
#' @param project Project directory containing `.sqlite`.
#' @param onStart,options,enableBookmarking,uiPattern Shiny app options.
#' @param host,port,launch.browser Shiny server options.
#' @param headless Offer a cluster job script instead of a local run. `NULL`
#'   (default) detects it; see [detect_headless()].
#' @param ... Passed to golem options.
#' @export
#' @importFrom shiny shinyApp
run_app <- function(project = getwd(), onStart = NULL, options = list(),
                    enableBookmarking = NULL, uiPattern = "/",
                    host = NULL, port = NULL, launch.browser = NULL,
                    headless = NULL, ...) {
  project <- normalizePath(project, mustWork = TRUE)
  if (!file.exists(file.path(project, ".sqlite"))) {
    stop("No RiboPilot project (.sqlite) found in: ", project, call. = FALSE)
  }
  check_nextflow_version("RiboPilot app")
  ensure_hydra_setup()
  if (is.null(launch.browser)) launch.browser <- interactive()
  # Headless = this R session is on a remote host (RStudio Server on a login
  # node, or a tunneled session). The run modal then offers a cluster job script
  # instead of running Nextflow where the app is served.
  if (is.null(headless)) headless <- detect_headless(launch.browser)
  headless <- isTRUE(headless)
  options(RiboPilot.headless = headless)
  # only a session with no browser of its own needs tunnel instructions;
  # RStudio Server is headless but serves the app through its own proxy
  if (!isTRUE(launch.browser) && !is.null(port)) {
    message("RiboPilot is headless. Tunnel it from your machine, then open http://localhost:",
      port, " :\n  ssh -N -L ", port, ":localhost:", port, " <user>@<this-host>")
  }
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
