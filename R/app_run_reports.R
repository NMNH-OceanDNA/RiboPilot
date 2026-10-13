#' Result cell: icon and label in a pill
#' @noRd
run_result_html <- function(result) {
  tone <- c(finished = "success", failures = "warning", failed = "danger",
    stopped = "neutral", unfinished = "info")
  ic <- c(finished = "circle-check", failures = "triangle-exclamation",
    failed = "circle-xmark", stopped = "circle-stop", unfinished = "hourglass-half")
  vapply(result, function(r) {
    if (is.na(r) || !r %in% names(tone)) return("")
    sprintf("<span class='mp-pill mp-pill-%s'>%s %s</span>",
      tone[[r]], as.character(mp_icon(ic[[r]])), RUN_RESULTS[[r]])
  }, "", USE.NAMES = FALSE)
}

#' Run Reports window and the new-report notice
#'
#' Inputs `rr_open` (toolbar button), `rr_tbl`, `rr_open_loc`, `rr_notice`.
#' Sets `session$userData$run_reports_seen(base)`, which marks a report the
#' app just showed as new without another notice.
#' @param project function() -> project folder.
#' @param refresh Reactive that fires on the table reload button.
#' @noRd
run_reports_handlers <- function(input, output, session, project, refresh) {
  runs <- shiny::reactiveVal(NULL)
  opened <- shiny::reactiveVal(0)
  new_bases <- character(0)

  sync <- function() {
    live <- session$userData$run_live_log
    new <- tryCatch(sync_run_reports(project(), skip = if (!is.null(live)) run_log_base(live)),
      error = function(e) character(0))
    new_bases <<- union(new_bases, new)
    new
  }
  session$userData$run_reports_seen <- function(base) new_bases <<- union(new_bases, base)

  show_reports <- function() {
    sync()
    df <- tryCatch(list_run_reports(project()), error = function(e) NULL)
    runs(df)
    opened(opened() + 1)
    body <- if (is.null(df) || nrow(df) == 0) {
      shiny::p("No run logs found for this project yet.")
    } else {
      shiny::tagList(reactable::reactableOutput("rr_tbl"), shiny::uiOutput("rr_detail"))
    }
    shiny::showModal(shiny::modalDialog(title = mp_modal_title("Run reports"),
      size = "l", easyClose = TRUE, body, footer = mp_footer(dismiss = "Close")))
  }
  shiny::observeEvent(input$rr_open, show_reports())

  output$rr_tbl <- reactable::renderReactable({
    opened()
    df <- runs()
    shiny::req(!is.null(df), nrow(df) > 0)
    disp <- data.frame(
      New = ifelse(df$base %in% new_bases, "New", ""),
      Workflow = ifelse(df$workflow == "reference", "Reference build",
        ifelse(df$workflow == "recheck", "Recheck", "Pipeline")),
      Started = df$started,
      Duration = df$duration,
      Result = run_result_html(df$result),
      Samples = df$n_samples,
      `Failed samples` = df$n_failed_samples,
      `Failed tasks` = df$n_failed_tasks,
      Launched = ifelse(df$launch %in% "job", "Job", "App"),
      check.names = FALSE, stringsAsFactors = FALSE)
    reactable::reactable(disp, selection = "single", onClick = "select",
      height = if (nrow(disp) > 10) 380 else "auto",
      pagination = FALSE, compact = TRUE, highlight = TRUE, wrap = FALSE,
      columns = list(
        New = reactable::colDef(name = "", maxWidth = 60, html = TRUE,
          cell = rt_pill(map = c(New = "info"), empty = "")),
        Started = reactable::colDef(minWidth = 150),
        Result = reactable::colDef(html = TRUE, minWidth = 190)))
  })

  selected_run <- function() {
    i <- reactable::getReactableState("rr_tbl", "selected")
    df <- runs()
    shiny::req(length(i) == 1, !is.null(df), i <= nrow(df))
    df[i, , drop = FALSE]
  }

  output$rr_detail <- shiny::renderUI({
    run <- selected_run()
    rpt <- run$report
    if (is.na(rpt) || !file.exists(rpt)) {
      log <- run$log
      if (is.na(log) || !file.exists(log)) {
        return(shiny::p(style = "margin-top: 1em;", "No report yet: the run has not recorded a finish."))
      }
      shiny::invalidateLater(15000)
      return(shiny::div(style = "margin-top: 1em;",
        shiny::p("No report yet: the run has not recorded a finish. Progress from its Nextflow log, ",
          "updated every 15 seconds. A run with no log update for a long time may have stopped."),
        shiny::tags$pre(style = "max-height: 40vh; overflow: auto; white-space: pre; font-size: 11px;",
          run_progress_text(log))))
    }
    txt <- tryCatch(paste(readLines(rpt, warn = FALSE), collapse = "\n"),
      error = function(e) "The report could not be read.")
    shiny::div(style = "margin-top: 1em;",
      shiny::tags$pre(id = "rr_report_text",
        style = "max-height: 40vh; overflow: auto; white-space: pre-wrap; font-size: 11px;", txt),
      shiny::tags$button(type = "button", class = "btn btn-default btn-sm",
        onclick = "rpCopyEl('rr_report_text'); this.textContent = 'Copied';",
        mp_icon("copy"), " Copy report"),
      shiny::actionButton("rr_open_loc", "Open report location",
        icon = mp_icon("folder-open"), class = "btn-sm"))
  })

  shiny::observeEvent(input$rr_open_loc, {
    rpt <- selected_run()$report
    shiny::req(!is.na(rpt))
    if (!isTRUE(open_path(dirname(rpt)))) {
      rp_alert(paste0("This session has no file browser to open. The report is at: ", rpt),
        type = "message")
    }
  })

  check <- function() {
    new <- sync()
    if (length(new) == 0) return(invisible())
    mp_confirm("rr_notice", title = sprintf("New run reports available (%d)", length(new)),
      text = "Open Run Reports to see how the runs went.",
      action_label = "Open Run Reports", cancel_label = "Dismiss")
  }
  shiny::observeEvent(input$rr_notice, if (isTRUE(input$rr_notice)) show_reports())

  # held until the next reload when a startup alert is up (one alert at a time)
  session$onFlushed(function() {
    if (isTRUE(session$userData$startup_alert)) return(invisible())
    shiny::withReactiveDomain(session, shiny::isolate(check()))
  }, once = TRUE)
  shiny::observeEvent(refresh(), check(), ignoreInit = TRUE)
  invisible(list(runs = runs))
}
