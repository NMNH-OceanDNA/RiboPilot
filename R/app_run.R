#' Is the app running headless (served to a remote browser, e.g. a tunneled
#' cluster login node)?
#' @noRd
app_headless <- function() isTRUE(getOption("RiboPilot.headless", FALSE))

#' Detect a batch-scheduler submit command on PATH (sbatch/qsub/bsub).
#' @return list(cmd, name) or NULL.
#' @noRd
detect_scheduler <- function() {
  for (s in list(c("sbatch", "Slurm"), c("qsub", "SGE/PBS"), c("bsub", "LSF"))) {
    if (nzchar(unname(Sys.which(s[1])))) return(list(cmd = s[1], name = s[2]))
  }
  NULL
}

#' Wrap a Nextflow command in a scheduler job script
#'
#' A multi-hour run belongs in the queue, not on the login node. The script is a
#' starting point the user can edit in the run modal before saving/submitting.
#' @noRd
job_script <- function(project, cmd, sched, jobname) {
  hdr <- if (is.null(sched)) {
    character(0)
  } else {
    switch(sched$cmd,
      sbatch = c(sprintf("#SBATCH -J %s", jobname), "#SBATCH -c 2",
        "#SBATCH --mem=8G", "#SBATCH -t 48:00:00", sprintf("#SBATCH -o %s_%%j.log", jobname)),
      qsub = c(sprintf("#$ -N %s", jobname), "#$ -cwd", "#$ -pe mpi 2",
        "#$ -l h_vmem=8G", "#$ -j y", sprintf("#$ -o %s.$JOB_ID.log", jobname)),
      bsub = c(sprintf("#BSUB -J %s", jobname), "#BSUB -n 2",
        "#BSUB -M 8000", sprintf("#BSUB -o %s.%%J.log", jobname)),
      character(0))
  }
  pin <- nf_pin_version()
  pin_line <- if (is.na(pin)) character(0) else paste0("export NXF_VER=", pin)
  paste(c("#!/bin/bash", hdr, "", pin_line, sprintf("cd %s", shQuote(project)), cmd, ""),
    collapse = "\n")
}

#' Cluster job script for the main workflow
#' @noRd
build_submit_script <- function(project, resume = TRUE, sched = detect_scheduler()) {
  # one arg (the container bind options) contains spaces, so quote per argument
  args <- c(nf_command(project, resume = resume), "-ansi-log", "false")
  cmd <- paste(c(shQuote("nextflow"), shQuote(args)), collapse = " ")
  job_script(project, cmd, sched, "ribopilot")
}

#' Cluster job script for a reference BUILD run
#'
#' The manifest must already be on disk: the job reads it when it runs, which is
#' long after the modal that wrote the script has closed.
#' @noRd
build_reference_submit_script <- function(project, manifest, ref_dir, api_key = "",
                                          refseq = FALSE, resume = TRUE,
                                          sched = detect_scheduler()) {
  args <- nf_reference_command(project, "BUILD", manifest, ref_dir, ref_dir,
    api_key = api_key, refseq = refseq, resume = resume)
  args <- c(args, "-ansi-log", "false")
  job_script(project, paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "),
    sched, "ribopilot_ref")
}

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
  # upgrade older project DBs to the current schema before the app opens
  migrate_project(project)
  check_nextflow_version("RiboPilot app")
  if (is.null(launch.browser)) launch.browser <- interactive()
  # Headless = launched without a local browser, i.e. served on a remote host
  # (typically an SSH-tunneled cluster login node). In that case the run modal
  # submits to the scheduler instead of running Nextflow on the login node.
  headless <- !isTRUE(launch.browser)
  options(RiboPilot.headless = headless)
  if (headless && !is.null(port)) {
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
