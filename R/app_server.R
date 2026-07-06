#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # Stop the app when the last browser session closes (survives refresh).
  register_app_lifecycle(session)

  refresh <- reactiveVal(0)
  observeEvent(input$refresh, refresh(refresh() + 1))
  sel_id <- reactiveVal(NULL)

  # ---- Work Dirs browser (Nextflow task dirs per sample) ----
  observeEvent(input$workdirs_btn, {
    ids <- app_query("SELECT ID FROM samples ORDER BY ID")$ID
    presel <- selected_ids()
    presel <- if (length(presel)) presel[1] else if (length(ids)) ids[1] else NULL
    showModal(modalDialog(
      title = "Sample work directories", size = "l", easyClose = TRUE,
      tags$p(class = "text-muted", style = "font-size:12px; margin-bottom:0.75em;",
        "Nextflow task work directories for the chosen sample (failed attempts included). ",
        "Use these to inspect intermediate files or troubleshoot a failed process."),
      shinyWidgets::pickerInput("wd_sample", "Sample", choices = ids, selected = presel,
        options = list(`live-search` = TRUE)),
      div(style = "max-height:60vh; overflow:auto;", reactable::reactableOutput("wd_tbl")),
      footer = modalButton("Close")))
  })
  wd_rows <- reactive({
    req(input$wd_sample)
    find_workdirs(app_project(), input$wd_sample)
  })
  output$wd_tbl <- reactable::renderReactable({
    d <- wd_rows()
    if (nrow(d) == 0) {
      return(reactable::reactable(data.frame(` ` = "No work directories found (run the pipeline first).",
        check.names = FALSE)))
    }
    d$open <- "open"
    reactable::reactable(d[, c("process", "status", "modified", "open", "workdir")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 25, wrap = FALSE,
      columns = list(
        status = reactable::colDef(html = TRUE, width = 80, cell = rt_dynamicIcon(
          c(success = "fa fa-circle-check", failed = "fa fa-triangle-exclamation"))),
        open = reactable::colDef(name = "", width = 60, html = TRUE, cell = rt_link("wd_open")),
        workdir = reactable::colDef(minWidth = 340, html = TRUE, cell = rt_longtext())))
  })
  observeEvent(input$wd_open, {
    i <- suppressWarnings(as.integer(input$wd_open))
    d <- wd_rows()
    req(!is.na(i), i >= 1, i <= nrow(d))
    if (!isTRUE(d$exists[i])) {
      showNotification("Work directory not found on this host.", type = "warning")
      return()
    }
    open_path(d$workdir[i])
  })

  # ---- samples table data (per mode; includes state + lock) ----
  samples_data <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    # Column order mirrors MitoPilot: lock, state, ID, Taxon lead the
    # assemble/annotate tables; export leads with ID/Taxon and puts group last.
    if (mode == "Assemble") {
      app_query(
        "SELECT a.assemble_lock, a.assemble_switch, s.ID, s.Taxon, a.paths, a.scaffolds, a.length,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           a.assemble_notes AS notes
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID")
    } else if (mode == "Export") {
      # export only shows samples whose annotation is locked (validated)
      app_query(
        "SELECT s.ID, s.Taxon, st.reference_accession AS reference, st.length, an.n_features,
           an.missing, an.reviewed, an.annotate_lock, s.export_group
         FROM samples s JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID
         WHERE an.annotate_lock = 1")
    } else {
      # annotate only shows samples whose assembly is locked (validated)
      app_query(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, s.Taxon,
           st.reference_accession AS reference, an.n_features, an.missing,
           an.n_flags, an.reviewed
         FROM samples s JOIN assemble a ON s.ID = a.ID
           LEFT JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID
         WHERE a.assemble_lock = 1")
    }
  })

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    d <- samples_data()
    lock_col <- if (mode == "Annotate") "annotate_lock" else "assemble_lock"
    state_col <- if (mode == "Annotate") "annotate_switch" else "assemble_switch"

    cols <- list(
      ID = reactable::colDef(minWidth = 130, sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      Taxon = reactable::colDef(minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      reviewed = reactable::colDef(name = "reviewed", html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      n_flags = reactable::colDef(name = "flags", width = 70, align = "center", html = TRUE,
        filterable = FALSE, cell = rt_flag_badge()),
      reference = reactable::colDef(html = TRUE, cell = rt_ncbi_link()),
      export_group = reactable::colDef(name = "group", minWidth = 110, html = TRUE, cell = rt_longtext()),
      annotate_lock = reactable::colDef(name = "", width = 34, align = "center", filterable = FALSE,
        html = TRUE, cell = rt_dynamicIcon(c(`0` = "fa fa-lock-open", `1` = "fa fa-lock"))),
      notes = reactable::colDef(minWidth = 160, html = TRUE, cell = rt_longtext())
    )
    if (lock_col %in% names(d)) {
      cols[[lock_col]] <- reactable::colDef(name = "", sticky = "left", html = TRUE, width = 34,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(c(`0` = "fa fa-lock-open", `1` = "fa fa-lock")))
    }
    if (state_col %in% names(d)) {
      cols[[state_col]] <- reactable::colDef(name = "", sticky = "left", html = TRUE, width = 32,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(c(`1` = "fa fa-hourglass", `2` = "fa fa-circle-check",
          `3` = "fa fa-triangle-exclamation", `4` = "fa fa-circle-half-stroke")))
    }
    rowclass <- reactable::JS(sprintf(
      "function(rowInfo){ if(!rowInfo||!rowInfo.values) return '';
         return 'rp-lock-'+rowInfo.values['%s']+' rp-state-'+rowInfo.values['%s']; }",
      lock_col, state_col))

    reactable::reactable(
      d,
      resizable = TRUE, filterable = TRUE, striped = TRUE, compact = TRUE,
      searchable = TRUE, showPageSizeOptions = TRUE, wrap = FALSE, height = "100%",
      selection = "multiple", onClick = "select",
      defaultPageSize = 50, pageSizeOptions = c(25, 50, 100, 200),
      rowStyle = rt_highlight_row(), rowClass = rowclass,
      columns = cols[intersect(names(cols), names(d))]
    )
  })

  # ---- open the details modal when an ID cell is clicked ----
  observeEvent(input$open_detail, {
    i <- suppressWarnings(as.integer(input$open_detail))
    if (is.na(i)) {
      return()
    }
    id <- samples_data()$ID[i]
    sel_id(id)
    mode <- input$mode %||% "Assemble"
    if (mode == "Assemble") {
      showModal(assemble_modal(id))
    } else if (mode == "Annotate") {
      locked <- app_query("SELECT assemble_lock FROM assemble WHERE ID = ?", params = list(id))$assemble_lock
      if (!isTRUE(locked[1] == 1)) {
        showNotification("Lock the assembly (Assemble mode) before annotating this sample.",
          type = "warning", duration = 6)
        return()
      }
      showModal(annotate_modal(id))
    }
  })

  selected_ids <- function() {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel)) character(0) else samples_data()$ID[sel]
  }

  # ---- bulk lock / state / reviewed actions ----
  set_flag <- function(sql) {
    ids <- selected_ids()
    if (length(ids) == 0) {
      showNotification("No samples selected.", type = "warning")
      return()
    }
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    for (id in ids) DBI::dbExecute(con, sql, params = list(id))
    refresh(refresh() + 1)
    showNotification(paste("Updated", length(ids), "sample(s)."), type = "message")
  }
  observeEvent(input$lock_btn, {  # toggle lock on selected
    col <- if ((input$mode %||% "") == "Annotate") "annotate" else "assemble"
    set_flag(sprintf("UPDATE %s SET %s_lock = 1 - COALESCE(%s_lock, 0) WHERE ID = ?", col, col, col))
  })
  observeEvent(input$state_btn, {  # set state Ready (re-queue)
    col <- if ((input$mode %||% "") == "Annotate") "annotate" else "assemble"
    set_flag(sprintf("UPDATE %s SET %s_switch = 1 WHERE ID = ?", col, col))
  })
  observeEvent(input$reviewed_btn, set_flag("UPDATE annotate SET reviewed = 'yes' WHERE ID = ?"))

  # ---- Export group actions ----
  observeEvent(input$assign_group_btn, {
    if (length(selected_ids()) == 0) {
      showNotification("No samples selected.", type = "warning"); return()
    }
    showModal(group_modal())
  })
  observeEvent(input$assign_group_go, {
    grp <- trimws(input$group_name %||% "")
    ids <- selected_ids()
    if (!nzchar(grp) || length(ids) == 0) {
      removeModal(); return()
    }
    assign_export_group(app_project(), ids, grp)
    removeModal(); refresh(refresh() + 1)
    showNotification(paste0("Assigned ", length(ids), " sample(s) to group '", grp, "'."), type = "message")
  })
  observeEvent(input$clear_group_btn, {
    ids <- selected_ids()
    if (length(ids) == 0) {
      showNotification("No samples selected.", type = "warning"); return()
    }
    clear_export_group(app_project(), ids)
    refresh(refresh() + 1)
    showNotification(paste("Cleared group for", length(ids), "sample(s)."), type = "message")
  })

  # ---- lock / state row filters (hide via CSS on rp-lock-/rp-state- classes) ----
  output$filter_css <- renderUI({
    all_lock <- c("0", "1")
    all_state <- c("1", "2", "3", "4")
    hide <- character(0)
    if (!is.null(input$lock_filter)) hide <- c(hide, paste0(".rp-lock-", setdiff(all_lock, input$lock_filter)))
    if (!is.null(input$state_filter)) hide <- c(hide, paste0(".rp-state-", setdiff(all_state, input$state_filter)))
    if (length(hide) == 0) {
      return(NULL)
    }
    tags$style(HTML(paste0(paste(hide, collapse = ", "), " { display:none !important; }")))
  })

  # ---- Update modal: run the workflow with live Nextflow progress ----
  wf_for_mode <- function() if ((input$mode %||% "") == "Annotate") "annotate" else "assemble"
  nf_process <- reactiveVal(NULL)
  prog_header <- reactiveVal(character(0))
  prog_proc <- reactiveVal(list())
  prog_summary <- reactiveVal("")
  prog_footer <- reactiveVal(character(0))

  output$nf_cmd_text <- renderText({
    args <- nf_command(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE))
    paste(c("nextflow", args), collapse = " ")
  })

  observeEvent(input$run_btn, {
    prog_header(character(0)); prog_proc(list()); prog_summary(""); prog_footer(character(0))
    showModal(run_modal(input$mode %||% "Assemble", nf_pending_count(app_project(), wf_for_mode())))
  })

  progress_update <- function(lines) {
    lines <- cli::ansi_strip(lines)
    for (ln in lines) {
      if (grepl("^\\[(PIPELINE|WORKDIR)\\]", ln)) {
        prog_header(unique(c(prog_header(), ln)))
      } else if (grepl("^\\[PROCESS ", ln)) {
        key <- sub("^\\[PROCESS \\S+\\]\\s*", "", ln)
        pl <- prog_proc(); pl[[key]] <- ln; prog_proc(pl)
      } else if (grepl("^\\[(SUCCESS|FAILED)\\]", ln)) {
        prog_summary(ln)
      } else if (nzchar(trimws(ln))) {
        prog_footer(c(prog_footer(), ln))
      }
    }
  }

  observeEvent(input$nf_start, {
    prog_header(character(0)); prog_proc(list()); prog_summary(""); prog_footer(character(0))
    shinyjs::hide("nf_start"); shinyjs::show("nf_stop"); shinyjs::show("nf_progress")
    shinyjs::removeClass(id = "rp_gears", class = "paused")
    p <- processx::process$new("nextflow",
      args = nf_command(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE)),
      stdout = "|", stderr = "|", wd = app_project())
    nf_process(p)
  })

  observe({
    req(nf_process())
    invalidateLater(250)
    p <- nf_process()
    grab <- function() {
      out <- tryCatch(p$read_output_lines(), error = function(e) character(0))
      err <- tryCatch(p$read_error_lines(), error = function(e) character(0))
      if (length(c(out, err))) progress_update(c(out, err))
    }
    if (p$is_alive()) {
      grab()
    } else {
      grab()
      nf_process(NULL)
      shinyjs::hide("nf_stop"); shinyjs::show("nf_start")
      shinyjs::addClass(id = "rp_gears", class = "paused")
      refresh(refresh() + 1)
      showNotification("Workflow run finished.", type = "message")
    }
  })

  output$nf_prog_header <- renderText(paste(prog_header(), collapse = "\n"))
  output$nf_prog_processes <- renderText(paste(unlist(prog_proc()), collapse = "\n"))
  output$nf_prog_summary <- renderText(prog_summary())
  output$nf_prog_footer <- renderText(paste(utils::tail(prog_footer(), 40), collapse = "\n"))

  kill_nf <- function() {
    p <- nf_process()
    if (!is.null(p) && p$is_alive()) p$kill()
    nf_process(NULL)
  }
  # don't leave an orphaned Nextflow process if the browser closes mid-run
  session$onSessionEnded(function() {
    p <- isolate(nf_process())
    if (!is.null(p) && p$is_alive()) try(p$kill(), silent = TRUE)
  })
  observeEvent(input$nf_stop, {
    kill_nf()
    shinyjs::hide("nf_stop"); shinyjs::show("nf_start")
    shinyjs::addClass(id = "rp_gears", class = "paused")
    refresh(refresh() + 1)
  })
  observeEvent(input$nf_close, {
    kill_nf()
    removeModal()
  })

  # ---- Export ----
  output$export_sel_n <- renderText({
    refresh()
    paste(length(export_ids()), "locked & ready")
  })
  # gate export to annotation-locked samples (validated)
  export_ids <- function() {
    ids <- selected_ids()
    locked <- app_query("SELECT ID FROM annotate WHERE annotate_lock = 1")$ID
    if (length(ids)) intersect(ids, locked) else locked
  }
  output$export_dl <- downloadHandler(
    filename = function() paste0(basename(app_project()), "_export.zip"),
    content = function(file) {
      ids <- export_ids()
      if (length(ids) == 0) {
        showNotification("No annotation-locked samples to export. Lock annotations first.",
          type = "warning", duration = 6)
        stop("nothing to export")
      }
      td <- file.path(tempdir(), paste0("rp_export_", as.integer(Sys.time())))
      export_bundle(app_project(), ids = ids, out_dir = td, zip = TRUE)
      file.copy(paste0(td, ".zip"), file, overwrite = TRUE)
    }
  )

  # ---- Assemble modal ----
  output$sel_path_ui <- renderUI({
    id <- sel_id()
    req(id)
    d <- app_query("SELECT path, auto_select FROM path_scores WHERE ID = ? ORDER BY path", params = list(id))
    if (nrow(d) == 0) {
      return(em("No scored paths."))
    }
    radioButtons("sel_path", "Selected path:", choices = d$path,
      selected = d$path[d$auto_select == 1][1], inline = TRUE)
  })
  output$assemble_scores <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    d <- app_query(
      "SELECT path, length, round(depth,0) AS depth, anchors, blast_pident AS pident,
         hit_species, round(score,3) AS score, rank, auto_select
       FROM path_scores WHERE ID = ? ORDER BY rank", params = list(id))
    reactable::reactable(d, compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      rowStyle = function(index) if (isTRUE(d$auto_select[index] == 1)) list(background = "#eafaf0") else NULL)
  })
  output$assemble_scaffolds <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    reactable::reactable(
      app_query("SELECT path, scaffold, length, topology FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
        params = list(id)),
      compact = TRUE, resizable = TRUE, defaultPageSize = 10)
  })
  observeEvent(input$set_selected, {
    id <- sel_id()
    req(id, input$sel_path)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "UPDATE path_scores SET auto_select = CASE WHEN path = ? THEN 1 ELSE 0 END WHERE ID = ?",
      params = list(as.integer(input$sel_path), id))
    # queue for re-annotation (WF2 / Annotate mode "Run")
    DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1 WHERE ID = ?", params = list(id))
    refresh(refresh() + 1)
    showNotification(paste0("Set path ", input$sel_path,
      " as selected and queued for re-annotation. Run WF2 (Annotate mode -> Run)."), type = "message")
  })

  # ---- Annotate modal: feature map + nucleotide boundary editor ----
  edit_tick <- reactiveVal(0)  # bumps after each nudge to re-render map/seq/table

  output$feature_map <- renderPlot({
    edit_tick()
    plot_feature_map(sel_id())
  }, height = 155)

  output$coverage_plot <- renderPlot({
    edit_tick()
    plot_coverage(sel_id())
  }, height = 128)

  output$nudge_element_ui <- renderUI({
    id <- sel_id(); req(id)
    els <- app_query("SELECT element FROM annotations WHERE ID = ? ORDER BY start",
      params = list(id))$element
    if (length(els) == 0) {
      return(em("No annotations for this sample."))
    }
    selectInput("nudge_element", "Feature", choices = els, selected = els[1])
  })

  # active feature's stitched-sequence coords + context
  nudge_ctx <- reactive({
    edit_tick()
    id <- sel_id(); el <- input$nudge_element
    req(id, el)
    a <- app_query("SELECT start, stop, length, type FROM annotations WHERE ID = ? AND element = ?",
      params = list(id, el))
    if (nrow(a) == 0) return(NULL)
    seq <- app_query("SELECT sequence FROM stitched WHERE ID = ?", params = list(id))$sequence[1]
    list(el = el, start = a$start[1], stop = a$stop[1], type = a$type[1], seq = seq %||% "")
  })

  output$nudge_seq <- renderUI({
    cx <- nudge_ctx()
    if (is.null(cx) || !nzchar(cx$seq)) {
      return(em("No stitched sequence available for this sample."))
    }
    n <- nchar(cx$seq)
    flank <- 15
    lo <- max(1, cx$start - flank); hi <- min(n, cx$stop + flank)
    pre  <- substr(cx$seq, lo, cx$start - 1)
    feat <- substr(cx$seq, cx$start, cx$stop)
    post <- substr(cx$seq, cx$stop + 1, hi)
    tagList(
      span(sprintf("%s (%s)  %d..%d  len %d bp",
        cx$el, cx$type, cx$start, cx$stop, cx$stop - cx$start + 1),
        style = "display:block; color:#555; margin-bottom:4px;"),
      span(pre, style = "color:#aaa;"),
      span(feat, style = "color:#0072B2; font-weight:bold;"),
      span(post, style = "color:#aaa;")
    )
  })

  output$alignment_view <- renderUI({
    edit_tick()
    id <- sel_id(); el <- input$nudge_element
    if (is.null(id) || is.null(el)) {
      return(em("Select a feature."))
    }
    a <- app_query(
      "SELECT accession, pid, sample_aln, ref_aln FROM blast_ref_alignment WHERE ID = ? AND element = ?",
      params = list(id, el))
    if (nrow(a) == 0 || is.na(a$sample_aln[1])) {
      return(em("No reference alignment for this feature yet (needs WF1 BLAST + reference fetch)."))
    }
    tagList(
      span(sprintf("vs %s  -  %s: %.1f%% identity over %d bp",
        a$accession[1], el, a$pid[1], nchar(a$sample_aln[1])),
        style = "color:#555; display:block; margin-bottom:4px;"),
      render_alignment_html(a$sample_aln[1], a$ref_aln[1])
    )
  })

  output$feature_tbl <- reactable::renderReactable({
    id <- sel_id(); req(id)
    edit_tick()
    d <- app_query("SELECT element, type, start, stop, length, strand, status, flags, edited, tool FROM annotations WHERE ID = ? ORDER BY start",
      params = list(id))
    reactable::reactable(d, compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      columns = list(
        flags = reactable::colDef(html = TRUE, cell = rt_longtext()),
        status = reactable::colDef(width = 70)),
      rowStyle = function(index) if (isTRUE(d$status[index] == "flag")) list(background = "#fdf1e3") else NULL)
  })

  # ---- Curation flags panel (sample notes + per-feature flag list w/ jump) ----
  flag_rows <- reactive({
    edit_tick()
    id <- sel_id(); req(id)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND status = 'flag' ORDER BY start",
      params = list(id))
  })

  output$curate_notes <- renderUI({
    id <- sel_id(); req(id)
    a <- app_query("SELECT n_flags, annotate_notes FROM annotate WHERE ID = ?", params = list(id))
    if (nrow(a) == 0 || is.na(a$n_flags[1]) || a$n_flags[1] == 0) {
      return(span(style = "color:#2d6a4f;", HTML("&#10003; No curation flags for this sample.")))
    }
    tagList(
      span(sprintf("%d sample flag(s): ", a$n_flags[1]), style = "color:#7d4a1e; font-weight:bold;"),
      span(a$annotate_notes[1] %||% "", style = "color:#7d4a1e;")
    )
  })

  output$flag_tbl <- reactable::renderReactable({
    d <- flag_rows()
    if (nrow(d) == 0) {
      return(reactable::reactable(data.frame(` ` = "No flagged features.", check.names = FALSE),
        compact = TRUE))
    }
    d$jump <- "edit"
    reactable::reactable(d[, c("element", "flags", "jump")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      columns = list(
        flags = reactable::colDef(html = TRUE, cell = rt_longtext()),
        jump = reactable::colDef(name = "", width = 60, html = TRUE, cell = rt_link("flag_jump"))))
  })

  # click a flagged feature -> select it in the boundary editor below
  observeEvent(input$flag_jump, {
    i <- suppressWarnings(as.integer(input$flag_jump))
    d <- flag_rows()
    req(!is.na(i), i >= 1, i <= nrow(d))
    updateSelectInput(session, "nudge_element", selected = d$element[i])
  })

  # nudge one boundary by the chosen step, clamp to sequence, persist live
  nudge_boundary <- function(end, dir) {
    cx <- nudge_ctx(); req(cx)
    step <- as.integer(input$nudge_step %||% 1)
    n <- nchar(cx$seq)
    start <- cx$start; stop <- cx$stop
    if (end == "5") {
      start <- start + if (dir == "out") -step else step   # out = extend upstream
    } else {
      stop <- stop + if (dir == "out") step else -step      # out = extend downstream
    }
    start <- max(1, start); stop <- min(if (n > 0) n else stop, stop)
    if (start > stop) {
      showNotification("Boundaries would cross; nudge ignored.", type = "warning")
      return()
    }
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND element = ?",
      params = list(start, stop, stop - start + 1, sel_id(), cx$el))
    edit_tick(edit_tick() + 1)
  }
  observeEvent(input$nt_5_out, nudge_boundary("5", "out"))
  observeEvent(input$nt_5_in,  nudge_boundary("5", "in"))
  observeEvent(input$nt_3_out, nudge_boundary("3", "out"))
  observeEvent(input$nt_3_in,  nudge_boundary("3", "in"))
  observeEvent(input$mark_reviewed, {
    id <- sel_id()
    req(id)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "UPDATE annotate SET reviewed = 'yes' WHERE ID = ?", params = list(id))
    refresh(refresh() + 1)
    showNotification("Marked reviewed.", type = "message")
  })
}
