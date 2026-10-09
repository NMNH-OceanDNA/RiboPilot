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

#' Submit a job script; bsub reads #BSUB lines only from stdin
#' @noRd
submit_job <- function(cmd, f) {
  if (identical(cmd, "bsub")) system2(cmd, stdin = f, stdout = TRUE, stderr = TRUE)
  else system2(cmd, shQuote(f), stdout = TRUE, stderr = TRUE)
}

#' Project subdirectories holding run files
#'
#' Every run drops up to three files (a job script, the scheduler's job log, and
#' Nextflow's own log); its report goes to `.runs/reports/`. Older projects
#' keep their old `jobs/`, `logs/`, and `nextflow_logs/` folders.
#' @noRd
JOB_DIRS <- list(script = ".runs/jobs", log = ".runs/logs", nextflow = ".runs/nextflow")

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
  if (!(kind %in% c("slurm", "sge", "pbs", "pbspro", "lsf"))) {
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
                                base = job_basename("pipeline")) {
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
                                          base = job_basename("reference")) {
  args <- nf_reference_command(project, "BUILD", manifest, ref_dir, ref_dir,
    api_key = api_key, refseq = refseq, resume = resume)
  args <- c(nf_log_opt(project, base), args, "-ansi-log", "false")
  cmd <- c(paste("cd", shQuote(file.path(project, ".reference"))),
    paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "))
  job_script(project, cmd, sched, base)
}

#' Cluster job script for a recheck run
#'
#' Launched from its own folder, so the pipeline's -resume history in the
#' project folder is never touched; the job creates that folder itself.
#' @param base Timestamped name for this submission's script and logs.
#' @noRd
build_recheck_submit_script <- function(project, sched = detect_scheduler(),
                                        base = job_basename("recheck")) {
  rc <- file.path(project, ".recheck")
  args <- c(nf_log_opt(project, base), nf_recheck_command(project), "-ansi-log", "false")
  cmd <- c(paste("mkdir -p", shQuote(rc)), paste("cd", shQuote(rc)),
    paste(c(shQuote("nextflow"), shQuote(args)), collapse = " "))
  job_script(project, cmd, sched, base)
}

#' Open the RiboPilot app
#'
#' Launches the R-Shiny app for reviewing and curating a RiboPilot project.
#' Run from (or point `project` at) a project directory containing `.sqlite`.
#' `RiboPilot()` is the same function.
#'
#' @param project Project directory containing `.sqlite` (default: current dir).
#' @param host,port,launch.browser Passed to Shiny. To reach the app from
#'   another computer, pass a `port` and `launch.browser = FALSE`, then open
#'   an SSH tunnel to that port; the console prints the command.
#' @param headless Offer a cluster job script instead of Run from App. `NULL`
#'   (default) turns it on when a batch scheduler (sbatch, qsub, or bsub) is on
#'   PATH or the session is RStudio Server; pass `TRUE` or `FALSE` to override.
#' @export
#' @importFrom shiny shinyApp
run_app <- function(project = getwd(), host = NULL, port = NULL, launch.browser = NULL,
                    headless = NULL) {
  project <- normalizePath(project, mustWork = TRUE)
  if (!file.exists(file.path(project, ".sqlite"))) {
    stop("No RiboPilot project (.sqlite) found in: ", project,
      ". Make one with new_project(), or point `project` at an existing project folder.",
      call. = FALSE)
  }
  v <- project_schema_version(file.path(project, ".sqlite"))
  stop_unsupported_schema(v)
  if (v < current_schema_version()) {
    stop("This project needs an update. Run backwards_compatibility(\"", project, "\") first.",
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
  opts <- list(launch.browser = launch.browser)
  if (!is.null(host)) opts$host <- host
  if (!is.null(port)) opts$port <- port
  app <- shiny::shinyApp(ui = app_ui, server = app_server, options = opts)
  app$appOptions$project <- project
  app
}

#' @rdname run_app
#' @export
RiboPilot <- run_app
#' Wire one Nextflow run window: start, live progress, graceful stop, close
#'
#' Shared by the Assemble step's Update window (`pfx = ""`), the recheck window
#' on Annotation Review (`pfx = "rc_"`), and the Reference build window
#' (`pfx = "ref_"`). Inputs `<pfx>nf_start`, `<pfx>nf_stop`, `<pfx>nf_close`,
#' `<pfx>nf_close_confirm`; outputs `<pfx>nf_prog_all` and the
#' `<pfx>nf_progress` panel. A run holds the project's run lock until Nextflow
#' and every task it started have exited.
#' @param gears Id of the spinner in the window title.
#' @param project function() -> project folder.
#' @param kind Run kind named in the lock.
#' @param workflow Workflow named in the run's files ("pipeline", "recheck",
#'   "reference").
#' @param launch function(nflog) -> a running processx process, or NULL when
#'   there is nothing to run. May stop() with a message for the user.
#' @param on_lines function(lines) called with new output, or NULL.
#' @param on_done function(how) called when the run has ended; `how` is "done",
#'   "stopped", or "killed".
#' @param on_close function() called after the window closes.
#' @noRd
nf_run_window <- function(input, output, session, pfx, gears, project, kind, workflow,
                          launch, on_lines = NULL, on_done = function(how) NULL,
                          on_close = function() NULL) {
  id <- function(x) paste0(pfx, x)
  proc <- shiny::reactiveVal(NULL)
  state <- shiny::reactiveVal(nf_progress_empty())
  run <- new.env()
  output[[id("nf_prog_all")]] <- shiny::renderText(nf_progress_text(state()))

  running <- function(on) {
    shinyjs::toggle(id("nf_start"), condition = !on)
    shinyjs::toggle(id("nf_stop"), condition = on)
    shiny::updateActionButton(session, id("nf_close"), label = if (on) "Stop and Close" else "Close")
    if (on) shinyjs::removeClass(id = gears, class = "paused") else
      shinyjs::addClass(id = gears, class = "paused")
    shinyjs::runjs(if (on) "window.onbeforeunload = function() { return true; };" else
      "window.onbeforeunload = null;")
  }

  shiny::observeEvent(input[[id("nf_start")]], {
    if (!is.null(proc())) return()
    held <- run_lock_holder(project())
    if (!is.null(held)) {
      rp_alert(run_lock_message(held, project()), type = "warning")
      return()
    }
    marker <- NULL
    nflog <- NULL
    p <- tryCatch({
      marker <- run_lock_take(project(), kind)
      nflog <- nf_log_opt(project(), job_basename(workflow))[2]
      launch(nflog)
    }, error = function(e) e)
    if (!is.null(marker)) Sys.unsetenv(marker)
    if (is.null(p) || inherits(p, "error")) {
      if (!is.null(marker)) run_lock_release(project(), marker)
      if (inherits(p, "error")) {
        rp_alert(paste("Could not start the run.", conditionMessage(p)), type = "error")
      }
      return()
    }
    run$marker <- marker; run$log <- nflog
    session$userData$run_live_log <- nflog
    run$stop_at <- NULL; run$exit_at <- NULL; run$killed <- FALSE; run$termed <- FALSE
    state(nf_progress_empty())
    shinyjs::show(id("nf_progress"))
    running(TRUE)
    proc(p)
  })

  waited <- function(t) !is.null(t) && difftime(Sys.time(), t, units = "secs") > NF_STOP_WAIT
  report <- function(p, how) {
    status <- tryCatch(p$get_exit_status() %||% NA_integer_, error = function(e) NA_integer_)
    rep <- tryCatch(nf_run_report(shiny::isolate(project()), run$log, status, how),
      error = function(e) list(type = "error", title = "Run failed", text = "",
        lead = sprintf("The run ended with exit status %s; see %s", status %|NA|% "unknown", run$log)))
    if (!is.null(rep$meta)) try(write_run_report(shiny::isolate(project()), rep), silent = TRUE)
    session$userData$run_live_log <- NULL
    rep
  }

  shiny::observe({
    p <- proc()
    shiny::req(p)
    shiny::invalidateLater(250)
    lines <- tryCatch(c(p$read_output_lines(), p$read_error_lines()), error = function(e) character(0))
    if (length(lines)) {
      state(nf_progress_parse(shiny::isolate(state()), lines))
      if (!is.null(on_lines)) on_lines(lines)
    }
    if (p$is_alive()) {
      if (!is.null(run$stop_at)) nf_escalate(p, run)
      return()
    }
    if (is.null(run$exit_at)) run$exit_at <- Sys.time()
    if (run_busy(run$marker, project()) && !waited(run$exit_at)) return()
    run_clear_containers(project())
    proc(NULL)
    run_lock_release(project(), run$marker)
    running(FALSE)
    how <- if (run$killed) "killed" else if (!is.null(run$stop_at)) "stopped" else "done"
    rep <- report(p, how)
    on_done(how)
    rp_alert(run_report_ui(rep), type = rep$type, title = rep$title)
    if (is.function(session$userData$run_reports_seen)) session$userData$run_reports_seen(rep$meta$base)
  })

  stop_run <- function() {
    p <- proc()
    if (is.null(p) || !is.null(run$stop_at)) return()
    if (p$is_alive()) p$interrupt()
    run$stop_at <- Sys.time()
    shinyjs::hide(id("nf_stop"))
    shiny::updateActionButton(session, id("nf_close"), label = "Close")
    rp_toast("Stopping the run: Nextflow is stopping its tasks. This can take up to a minute.",
      type = "warning", duration = 10)
  }
  shiny::observeEvent(input[[id("nf_stop")]], stop_run())
  shiny::observeEvent(input[[id("nf_close")]], {
    if (!is.null(proc()) && is.null(run$stop_at)) {
      mp_confirm(id("nf_close_confirm"), title = "Stop the run?",
        text = paste("Closing this window stops the run. Nextflow stops its tasks first,",
          "which can take up to a minute."),
        action_label = "Stop and close", danger = TRUE)
      return()
    }
    shiny::removeModal()
    on_close()
  })
  shiny::observeEvent(input[[id("nf_close_confirm")]], {
    shiny::req(isTRUE(input[[id("nf_close_confirm")]]))
    stop_run()
    shiny::removeModal()
    on_close()
  })
  session$onSessionEnded(function() {
    p <- shiny::isolate(proc())
    if (is.null(p)) return()
    clean <- nf_stop_run(p, run$marker, shiny::isolate(project()))
    run_lock_release(shiny::isolate(project()), run$marker)
    how <- if (clean) "stopped" else "killed"
    report(p, how)
    try(on_done(how), silent = TRUE)
  })
  invisible(NULL)
}

#' Alert body for a run report: the lead line, the report, and a copy button
#' @noRd
run_report_ui <- function(rep) {
  if (rep$type == "message" || !nzchar(rep$text)) return(rep$lead)
  shiny::tagList(
    shiny::tags$p(rep$lead),
    shiny::tags$pre(id = "rp_run_report",
      style = "text-align:left; max-height:40vh; overflow:auto; white-space:pre-wrap; word-break:break-all; font-size:11px;",
      rep$text),
    shiny::tags$button(type = "button", class = "btn btn-default btn-sm",
      onclick = "rpCopyEl('rp_run_report'); this.textContent = 'Copied';",
      shiny::icon("clipboard"), " Copy report"))
}
