#' Is the app in cluster mode (job scripts instead of Run from App)?
#' @noRd
app_headless <- function() isTRUE(getOption("RiboPilot.headless", FALSE))

#' Should the app offer cluster job scripts instead of Run from App?
#'
#' A batch scheduler on PATH means this R session is on a cluster node, where a
#' Nextflow head process belongs in the queue. RStudio Server counts too: on a
#' login node it is interactive, and hydra_setup() only puts qsub on PATH for
#' Hydra.
#' @noRd
detect_headless <- function() {
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

#' Head-job header lines for one scheduler
#'
#' Queue, parallel environment, and cluster options (which carry the account)
#' come from the project config, so the head job goes where the profile sends
#' its tasks. The head job itself is small and long-lived: 2 CPUs, 8 GB, 48 h.
#' @param kind "slurm", "sge", "pbs", "pbspro", or "lsf"; anything else gives no header.
#' @param log Function of a file-name suffix returning the scheduler log path.
#' @noRd
job_header <- function(kind, jobname, log, cfg_lines = character(0)) {
  val <- function(key) {
    v <- config_param(cfg_lines, key)
    if (is.na(v) || !nzchar(v) || grepl("[{}]", v)) NULL else v
  }
  q <- val("queue"); opts <- val("clusterOptions"); penv <- val("penv")
  switch(kind,
    slurm = paste("#SBATCH", c(paste("-J", jobname), "-c 2", "--mem=8G", "-t 48:00:00",
      if (!is.null(q)) paste("-p", q), opts, paste("-o", log("_%j.log")))),
    sge = paste("#$", c(paste("-N", jobname), if (!is.null(penv)) paste("-pe", penv, "2"),
      "-l h_vmem=8G", "-l h_rt=48:00:00", if (!is.null(q)) paste("-q", q), opts,
      "-V", "-j y", "-S /bin/bash", paste("-o", log(".$JOB_ID.log")))),
    pbs = , pbspro = paste("#PBS", c(paste("-N", substr(jobname, 1, 15)),
      if (kind == "pbspro") "-l select=1:ncpus=2:mem=8gb" else "-l nodes=1:ppn=2,mem=8gb",
      "-l walltime=48:00:00", if (!is.null(q)) paste("-q", q), opts,
      "-V", "-j oe", "-S /bin/bash", paste("-o", log(".log")))),
    lsf = paste("#BSUB", c(paste("-J", jobname), "-n 2", "-M 8000", "-W 48:00",
      if (!is.null(q)) paste("-q", q), opts, paste("-o", log(".%J.log")))),
    character(0))
}

#' Wrap a Nextflow command in a scheduler job script
#'
#' A multi-hour run belongs in the queue, not on the login node. The script is a
#' starting point the user can edit in the run modal before saving/submitting.
#' @param sched Scheduler found on PATH (see detect_scheduler()); it picks the
#'   header only when the project config names no grid executor.
#' @noRd
job_script <- function(project, cmd, sched, jobname) {
  if (is_hydra_cluster()) {
    return(hydra_job_script(project, cmd, jobname))
  }
  cfg <- file.path(project, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg, warn = FALSE) else character(0)
  kind <- config_param(cfg_lines, "executor") %|NA|% ""
  if (kind %nin% c("slurm", "sge", "pbs", "pbspro", "lsf")) {
    kind <- switch(sched$cmd %||% "", sbatch = "slurm", qsub = "sge", bsub = "lsf", "")
  }
  # log paths are absolute: the app submits from wherever its R session happens
  # to be, which is not the project, so a relative -o lands somewhere arbitrary
  log <- function(suffix) job_path(project, "log", paste0(jobname, suffix))
  pin <- nf_pin_version()
  pin_line <- if (is.na(pin)) character(0) else paste0("export NXF_VER=", pin)
  paste(c("#!/bin/bash", job_header(kind, jobname, log, cfg_lines), "",
    "export NXF_OPTS='-Xms512m -Xmx4g'", pin_line, sprintf("cd %s", shQuote(project)), cmd,
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
#' @param host,port,launch.browser Passed to Shiny (see [run_app()]). To reach
#'   the app from another computer, pass a `port` and `launch.browser = FALSE`,
#'   then open an SSH tunnel to that port; the console prints the command.
#' @param headless Offer a cluster job script (Save script, Submit as job)
#'   instead of Run from App. By default this is on when a batch scheduler
#'   (sbatch, qsub, or bsub) is on PATH or the session is RStudio Server, and
#'   off otherwise. Pass `TRUE` or `FALSE` to override; `FALSE` runs Nextflow on
#'   the machine that serves the app.
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
#' @param headless Offer a cluster job script instead of Run from App. `NULL`
#'   (default) turns it on when a batch scheduler (sbatch, qsub, or bsub) is on
#'   PATH or the session is RStudio Server; pass `TRUE` or `FALSE` to override.
#' @param ... Passed to golem options.
#' @export
#' @importFrom shiny shinyApp
run_app <- function(project = getwd(), onStart = NULL, options = list(),
                    enableBookmarking = NULL, uiPattern = "/",
                    host = NULL, port = NULL, launch.browser = NULL,
                    headless = NULL, ...) {
  project <- normalizePath(project, mustWork = TRUE)
  if (!file.exists(file.path(project, ".sqlite"))) {
    stop("No RiboPilot project (.sqlite) found in: ", project,
      ". Make one with new_project(), or point `project` at an existing project folder.",
      call. = FALSE)
  }
  check_nextflow_version("RiboPilot app")
  ensure_hydra_setup()
  if (is.null(launch.browser)) launch.browser <- interactive()
  if (is.null(headless)) headless <- detect_headless()
  headless <- isTRUE(headless)
  options(RiboPilot.headless = headless, RiboPilot.local_browser = isTRUE(launch.browser))
  if (!isTRUE(launch.browser) && !is.null(port)) {
    message("Open http://localhost:", port, " in your browser. From another computer, ",
      "open a tunnel first, then use the same address:\n  ssh -N -L ", port,
      ":localhost:", port, " <user>@<this-host>")
  }
  if (headless) {
    message("Cluster mode: the Update window offers a job script to submit instead of ",
      "Run from App. Pass headless = FALSE to run the pipeline on this machine.")
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
