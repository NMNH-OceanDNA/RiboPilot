# Starting, stopping, and reporting on one Nextflow run.

#' Seconds Nextflow gets to stop its tasks before it is killed
#' @noRd
NF_STOP_WAIT <- 60

#' Seconds after an interrupt before Nextflow is sent SIGTERM
#'
#' A process started with SIGINT ignored (nohup, a non-interactive shell)
#' cannot be interrupted; Nextflow stops its tasks on SIGTERM too.
#' @noRd
NF_TERM_WAIT <- 10

#' Task exit codes nextflow.config retries with more memory
#' @noRd
NF_RETRY_CODES <- c(104L, 134L, 137L, 139L, 140L, 143L, 247L)

#' Environment for a Nextflow run started from R
#'
#' Wide COLUMNS keeps process names whole in the -ansi-log board, so the
#' progress parser's per-process keys stay stable.
#' @noRd
nf_env <- function() {
  env <- c("current", COLUMNS = "500", NXF_ANSI_SUMMARY = "true",
    NXF_DISABLE_CHECK_LATEST = "true")
  pin <- nf_pin_version()
  if (!is.na(pin)) env <- c(env, NXF_VER = pin)
  env
}

#' Where a project's run lock lives
#' @noRd
run_lock_path <- function(project) file.path(project, ".ribopilot_run.lock")

#' Is any process of a run still alive?
#'
#' Every process a run starts inherits the run's marker variable, including the
#' task wrappers that are still stopping containers after Nextflow has exited.
#' @noRd
run_alive <- function(marker) {
  isTRUE(nzchar(marker)) && length(suppressWarnings(ps::ps_find_tree(marker))) > 0
}

#' Output of a docker command, "" when docker is missing or fails
#' @noRd
docker_out <- function(args) {
  if (!nzchar(Sys.which("docker"))) return("")
  tryCatch(processx::run("docker", args, error_on_status = FALSE, timeout = 20)$stdout,
    error = function(e) "")
}

#' Nextflow task containers of a project
#'
#' A project's containers are the nxf-* ones that mount a folder inside it. Only
#' one run holds a project at a time, so these belong to that run.
#' @param all Include stopped containers.
#' @return Container ids.
#' @noRd
run_containers <- function(project, all = FALSE) {
  ids <- strsplit(trimws(docker_out(c("ps", "-q", "--no-trunc", if (all) "-a",
    "--filter", "name=nxf-"))), "\n")[[1]]
  ids <- ids[nzchar(ids)]
  if (!length(ids)) return(character(0))
  x <- strsplit(strsplit(docker_out(c("inspect", "--format",
    "{{.Id}}{{range .Mounts}}\t{{.Source}}{{end}}", ids)), "\n")[[1]], "\t")
  root <- unique(paste0(c(project, normalizePath(project, mustWork = FALSE)), "/"))
  hit <- vapply(x, function(m) any(startsWith(rep(m[-1], each = length(root)), root)), logical(1))
  vapply(x[hit], `[`, "", 1)
}

#' Is a run still busy: any of its processes, or a task container of its project?
#' @noRd
run_busy <- function(marker, project = NULL) {
  run_alive(marker) || (!is.null(project) && length(run_containers(project)) > 0)
}

#' Remove every task container a run left in its project
#' @noRd
run_clear_containers <- function(project) {
  ids <- run_containers(project, all = TRUE)
  if (length(ids)) docker_out(c("rm", "-f", ids))
  invisible(ids)
}

#' Take a project's run lock
#'
#' Sets the run's marker variable in this R session, so the next process it
#' starts carries it; unset it once that process is running.
#' @return The marker.
#' @noRd
run_lock_take <- function(project, kind) {
  marker <- paste0("RIBOPILOT_RUN_", Sys.getpid(), "_", format(Sys.time(), "%Y%m%d%H%M%S"))
  do.call(Sys.setenv, stats::setNames(list("1"), marker))
  write.dcf(data.frame(kind = kind, marker = marker, host = Sys.info()[["nodename"]],
    started = format(Sys.time(), "%Y-%m-%d %H:%M"), stringsAsFactors = FALSE),
    run_lock_path(project))
  marker
}

#' Release a project's run lock, if this run still holds it
#' @noRd
run_lock_release <- function(project, marker) {
  Sys.unsetenv(marker)
  f <- run_lock_path(project)
  held <- if (file.exists(f)) tryCatch(as.list(read.dcf(f)[1, ]), error = function(e) NULL)
  if (identical(held$marker, marker)) unlink(f)
  invisible(NULL)
}

#' The run holding a project's lock, or NULL
#'
#' A lock left by a run on this computer whose processes are all gone is stale
#' and is removed. A lock from another computer cannot be checked from here.
#' @noRd
run_lock_holder <- function(project) {
  f <- run_lock_path(project)
  held <- if (file.exists(f)) tryCatch(as.list(read.dcf(f)[1, ]), error = function(e) NULL)
  if (is.null(held)) return(NULL)
  if (identical(held$host, Sys.info()[["nodename"]]) && !run_alive(held$marker)) {
    unlink(f)
    return(NULL)
  }
  held
}

#' Why a new run cannot start, in words
#' @noRd
run_lock_message <- function(held, project) {
  if (identical(held$host, Sys.info()[["nodename"]])) {
    sprintf(paste("A %s run is already going for this project, or is still stopping its tasks",
      "(started %s). Wait for it to finish, or stop it in the window that started it."),
      held$kind, held$started)
  } else {
    sprintf(paste("A %s run for this project was started on %s at %s, and it cannot be checked",
      "from this computer. If that run is no longer going, delete %s and try again."),
      held$kind, held$host, held$started, run_lock_path(project))
  }
}

#' One step of stopping a run
#'
#' SIGTERM once NF_TERM_WAIT has passed since the interrupt, then the whole
#' process tree once NF_STOP_WAIT has.
#' @param run Environment with `stop_at`; `termed` and `killed` are set here.
#' @noRd
nf_escalate <- function(p, run, term_wait = NF_TERM_WAIT, wait = NF_STOP_WAIT) {
  s <- difftime(Sys.time(), run$stop_at, units = "secs")
  if (s > wait && !isTRUE(run$killed)) {
    tryCatch(p$kill_tree(), error = function(e) p$kill())
    run$killed <- TRUE
  } else if (s > term_wait && !isTRUE(run$termed)) {
    p$signal(tools::SIGTERM)
    run$termed <- TRUE
  }
}

#' Stop a run and wait for it: interrupt, then SIGTERM, kill only if it hangs
#'
#' SIGINT makes Nextflow stop its tasks, remove their containers, and close its
#' cache so -resume works. SIGKILL does none of that. Waits for the run's task
#' processes and containers, then removes containers still left.
#' @return TRUE when Nextflow stopped by itself.
#' @noRd
nf_stop_run <- function(p, marker, project = NULL, wait = NF_STOP_WAIT, term_wait = NF_TERM_WAIT) {
  run <- new.env()
  run$stop_at <- Sys.time()
  if (p$is_alive()) p$interrupt()
  while (p$is_alive()) {
    nf_escalate(p, run, term_wait, wait)
    Sys.sleep(0.25)
  }
  t0 <- Sys.time()
  while (run_busy(marker, project) && difftime(Sys.time(), t0, units = "secs") < wait) Sys.sleep(0.25)
  if (!is.null(project)) run_clear_containers(project)
  !isTRUE(run$killed)
}

#' Why a sample failed and what to try next, from its pipeline note
#' @return list(why, next_step)
#' @noRd
failure_hint <- function(note) {
  note <- note %|NA|% ""
  h <- function(why, next_step) list(why = why, next_step = next_step)
  asm <- paste("Click this sample's Assembly Options cell in the Assemble table, type a new name in",
    "Parameter set name (for example baited), and tick Edit, then ")
  bait <- "choose \"SPAdes, baited to the seed DB\" as the assembler"
  done <- ". Save, then press Update."
  if (grepl("below min_depth", note, fixed = TRUE)) {
    h("Too few reads were left after trimming to assemble this sample.",
      "Add more reads for it, or lower minDepth in the project's .config file. Then set its State to Ready to run and press Update.")
  } else if (grepl("exceed max_paths", note, fixed = TRUE)) {
    h("GetOrganelle found more alternative paths than Max paths allows.",
      paste0(asm, "raise Max paths, or ", bait, done))
  } else if (grepl("exceeds max_scaffolds", note, fixed = TRUE)) {
    h("The assembly came out in more contigs than Max contigs allows.",
      paste0(asm, "raise Max contigs, or ", bait, done))
  } else if (grepl("under min_assembly_length", note, fixed = TRUE)) {
    h("Every assembled contig was shorter than Min assembly length.",
      paste0(asm, bait, ", or lower Min assembly length", done))
  } else if (grepl("rRNA search incomplete", note, fixed = TRUE)) {
    h("Part of the search for rRNA genes in your assembly did not finish.",
      "Set its State to Ready to run and press Update. If it fails again, Work Dirs shows the failed find_rrna_scan step.")
  } else if (grepl("no rRNA genes found in the supplied assembly", note, fixed = TRUE)) {
    h("No rRNA genes were found in the assembly you supplied.",
      "Check that the Assembly file named in the mapping file is the right one. To assemble from reads instead, leave Assembly blank and give R1 and R2.")
  } else if (grepl("no rDNA contig confirmed", note, fixed = TRUE)) {
    h("Contigs in your assembly had rRNA hits, but none had an 18S or 28S gene to anchor it.",
      "See \"Contigs not carried forward\" below. If the rDNA is split across contigs, assemble from reads instead.")
  } else if (grepl("no usable (seed|label) database", note)) {
    h("This sample has no seed or label database the fallback assembler could read, so nothing was assembled.",
      paste0(asm, "set Seed database and Label database to FASTA files you have (the Reference step builds them for a group)", done))
  } else if (grepl("return value: -9|out of memory|bad_alloc", note)) {
    h("The assembler ran out of memory.", paste0(asm, "raise Memory (GB)", done))
  } else {
    h("The assembler returned no rDNA contigs.",
      paste0(asm, bait, done, " Work Dirs shows the assembler's own log."))
  }
}
