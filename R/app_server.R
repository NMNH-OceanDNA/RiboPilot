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
  show_workdirs <- function(preselect = NULL) {
    ids <- app_query("SELECT ID FROM samples ORDER BY ID")$ID
    presel <- preselect %||% {
      s <- selected_ids(); if (length(s)) s[1] else if (length(ids)) ids[1] else NULL
    }
    showModal(modalDialog(
      title = "Sample output / work directories", size = "l", easyClose = TRUE,
      tags$p(class = "text-muted", style = "font-size:12px; margin-bottom:0.75em;",
        "Nextflow task work directories for the chosen sample (failed attempts included). ",
        "Use these to inspect intermediate + output files or troubleshoot a failed process."),
      shinyWidgets::pickerInput("wd_sample", "Sample", choices = ids, selected = presel,
        options = list(`live-search` = TRUE)),
      div(style = "max-height:60vh; overflow:auto;", reactable::reactableOutput("wd_tbl")),
      footer = modalButton("Close")))
  }
  observeEvent(input$workdirs_btn, show_workdirs())
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

  # mode -> (table, prefix) for lock/state columns and updates
  mode_meta <- function(mode) {
    switch(mode %||% "Assemble",
      Taxonomy = list(tbl = "samples", pfx = "taxonomy"),
      Annotate = list(tbl = "annotate", pfx = "annotate"),
      list(tbl = "assemble", pfx = "assemble"))
  }

  # ---- samples table data (per mode; includes state + lock) ----
  samples_data <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    # Column order mirrors MitoPilot: lock, state, ID, Taxon lead the tables.
    if (mode == "Taxonomy") {
      app_query(
        "SELECT s.taxonomy_lock, s.taxonomy_switch, s.ID, s.Taxon,
           s.sci_name, s.taxid, s.rank, s.match_type
         FROM samples s")
    } else if (mode == "Assemble") {
      app_query(
        "SELECT a.assemble_lock, a.assemble_switch, s.ID, s.Taxon, s.sci_name, s.taxid,
           pp.pre_opts, a.assemble_opts, a.blast_opts,
           a.paths, a.scaffolds, a.length,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           a.blast_species, a.blast_accession, a.blast_pident, a.blast_qcovs,
           a.assemble_notes AS notes,
           'Details' AS details, 'Output' AS output
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID
           LEFT JOIN preprocess pp ON s.ID = pp.ID
         WHERE s.taxonomy_lock = 1")
    } else if (mode == "Export") {
      # export only shows samples whose annotation is locked (validated)
      app_query(
        "SELECT s.ID, s.Taxon, s.sci_name, s.taxid, st.reference_accession AS reference, st.length,
           an.n_features, an.missing, an.reviewed, an.annotate_lock, s.export_group
         FROM samples s JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID
         WHERE an.annotate_lock = 1")
    } else {
      # annotate only shows samples whose assembly is locked (validated)
      app_query(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, s.Taxon, s.sci_name, s.taxid,
           st.reference_accession AS reference, an.n_features, an.missing,
           an.n_flags, an.annotate_opts, an.curate_opts, an.reviewed
         FROM samples s JOIN assemble a ON s.ID = a.ID
           LEFT JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID
         WHERE a.assemble_lock = 1")
    }
  })

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    d <- samples_data()
    mm <- mode_meta(mode)
    lock_col <- paste0(mm$pfx, "_lock")
    state_col <- paste0(mm$pfx, "_switch")

    cols <- list(
      ID = reactable::colDef(name = "ID", minWidth = 130, sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      Taxon = reactable::colDef(name = "Taxon", minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      sci_name = reactable::colDef(name = "Assigned Taxon", minWidth = 200, html = TRUE,
        cell = rt_taxid_link()),
      taxid = reactable::colDef(show = FALSE),
      rank = reactable::colDef(name = "Rank", minWidth = 90, html = TRUE, cell = rt_longtext()),
      match_type = reactable::colDef(name = "Confidence", minWidth = 110, html = TRUE,
        cell = rt_confidence_badge()),
      reviewed = reactable::colDef(name = "Reviewed", html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      n_flags = reactable::colDef(name = "Flags", width = 70, html = TRUE,
        filterable = FALSE, cell = rt_flag_badge()),
      reference = reactable::colDef(name = "Reference", html = TRUE, cell = rt_ncbi_link()),
      blast_species = reactable::colDef(name = "BLAST Hit", minWidth = 160, html = TRUE, cell = rt_longtext()),
      blast_accession = reactable::colDef(name = "Accession", minWidth = 110, html = TRUE, cell = rt_ncbi_link()),
      blast_pident = reactable::colDef(name = "BLAST %ID", width = 90),
      blast_qcovs = reactable::colDef(name = "BLAST Cov", width = 90),
      details = reactable::colDef(name = "", width = 55, align = "center", filterable = FALSE, html = TRUE,
        cell = rt_icon_link("asm_details", "fa fa-up-right-from-square", "View assembly details")),
      output = reactable::colDef(name = "", width = 55, align = "center", filterable = FALSE, html = TRUE,
        cell = rt_icon_disabled("fa fa-folder-open", "Output browser (disabled for now)")),
      n_features = reactable::colDef(name = "Features", width = 90, html = TRUE),
      missing = reactable::colDef(name = "Missing", minWidth = 120, html = TRUE, cell = rt_longtext()),
      length = reactable::colDef(name = "Length", width = 90),
      paths = reactable::colDef(name = "Paths", width = 80),
      scaffolds = reactable::colDef(name = "Scaffolds", width = 90),
      sel_path = reactable::colDef(name = "Selected Path", width = 110),
      best_score = reactable::colDef(name = "Best Score", width = 100),
      export_group = reactable::colDef(name = "Group", minWidth = 110, html = TRUE, cell = rt_longtext()),
      pre_opts = reactable::colDef(name = "Preprocess", minWidth = 100, html = TRUE,
        filterable = FALSE, cell = rt_link("set_pre_opts")),
      assemble_opts = reactable::colDef(name = "Assemble Opts", minWidth = 100, html = TRUE,
        filterable = FALSE, cell = rt_link("set_assemble_opts")),
      blast_opts = reactable::colDef(name = "BLAST", minWidth = 85, html = TRUE,
        filterable = FALSE, cell = rt_link("set_blast_opts")),
      annotate_opts = reactable::colDef(name = "Annotate Opts", minWidth = 100, html = TRUE,
        filterable = FALSE, cell = rt_link("set_annotate_opts")),
      curate_opts = reactable::colDef(name = "Curate", minWidth = 90, html = TRUE,
        filterable = FALSE, cell = rt_link("set_curate_opts")),
      annotate_lock = reactable::colDef(name = "", width = 34, align = "center", filterable = FALSE,
        html = TRUE, cell = rt_dynamicIcon(c(`0` = "fa fa-lock-open", `1` = "fa fa-lock"))),
      notes = reactable::colDef(name = "Notes", minWidth = 160, html = TRUE, cell = rt_longtext())
    )
    if (lock_col %in% names(d)) {
      cols[[lock_col]] <- reactable::colDef(name = "", sticky = "left", html = TRUE, width = 34,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(c(`0` = "fa fa-lock-open", `1` = "fa fa-lock")))
    }
    if (state_col %in% names(d)) {
      cols[[state_col]] <- reactable::colDef(name = "", sticky = "left", html = TRUE, width = 32,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(c(`0` = "fa fa-hourglass", `1` = "fa fa-person-running",
          `2` = "fa fa-circle-check", `3` = "fa fa-triangle-exclamation",
          `4` = "fa fa-circle-half-stroke")))
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
      # consistent header capitalization (Title Case names above) + left alignment
      defaultColDef = reactable::colDef(align = "left", headerVAlign = "center"),
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
    if (mode == "Taxonomy") {
      open_taxon_review(id)
    } else if (mode == "Assemble") {
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

  # ---- Assemble table "Details" / "Output" buttons ----
  observeEvent(input$asm_details, {
    i <- suppressWarnings(as.integer(input$asm_details))
    req(!is.na(i))
    id <- samples_data()$ID[i]; sel_id(id)
    showModal(assemble_modal(id))
  })
  observeEvent(input$asm_output, {
    i <- suppressWarnings(as.integer(input$asm_output))
    req(!is.na(i))
    show_workdirs(samples_data()$ID[i])
  })

  selected_ids <- function() {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel)) character(0) else samples_data()$ID[sel]
  }

  # ---- named parameter-set editing (preprocess/assemble/blast/annotate/curate) ----
  opts_handlers(input, output, session, list(
    samples_data = samples_data, selected_ids = selected_ids,
    refresh = refresh, con = app_con, mode = function() input$mode %||% "Assemble"
  ))

  # ---- taxonomy review (local search + manual assignment) ----
  app_taxdump <- function() project_taxdump_dir(app_project())
  open_taxon_review <- taxonomy_handlers(input, output, session, list(
    con = app_con, refresh = refresh, taxdump = app_taxdump,
    samples_for = function(ids) {
      if (length(ids) == 0) return(data.frame())
      ph <- paste(rep("?", length(ids)), collapse = ",")
      app_query(sprintf("SELECT ID, Taxon, sci_name, taxid, rank FROM samples WHERE ID IN (%s)", ph),
        params = as.list(ids))
    }
  ))
  observeEvent(input$edit_taxon_btn, open_taxon_review(selected_ids()))

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
    mm <- mode_meta(input$mode)
    set_flag(sprintf("UPDATE %s SET %s_lock = 1 - COALESCE(%s_lock, 0) WHERE ID = ?",
      mm$tbl, mm$pfx, mm$pfx))
  })
  observeEvent(input$state_btn, {  # open a picker to choose the target state
    ids <- selected_ids()
    if (length(ids) == 0) {
      showNotification("No samples selected.", type = "warning"); return()
    }
    mm <- mode_meta(input$mode)
    cur <- app_query(sprintf("SELECT %s_switch AS s FROM %s WHERE ID = ?", mm$pfx, mm$tbl),
      params = list(ids[1]))$s
    showModal(state_modal(mm$pfx, cur))
  })
  observeEvent(input$update_state, {
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state %||% "1")
    set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st))
    removeModal()
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
  prog_executor <- reactiveVal("")
  prog_proc <- reactiveVal(list())
  prog_footer <- reactiveVal(character(0))

  output$nf_cmd_text <- renderText({
    args <- nf_command(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE))
    paste(c("nextflow", args), collapse = " ")
  })

  observeEvent(input$run_btn, {
    # Taxonomy runs locally in the app (lightweight taxdump lookup, no Nextflow).
    if ((input$mode %||% "") == "Taxonomy") { run_taxonomy_update(); return() }
    prog_header(character(0)); prog_executor(""); prog_proc(list()); prog_footer(character(0))
    showModal(run_modal(input$mode %||% "Assemble", nf_pending_count(app_project(), wf_for_mode())))
  })

  # Local (in-app) taxonomy auto-assignment for the ready, unlocked samples.
  # Shows a "hold tight" popup, then defers the (blocking) taxdump resolution so
  # the popup renders first, then swaps in a success/error popup (MitoPilot style).
  run_taxonomy_update <- function() {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ids <- DBI::dbGetQuery(con,
      "SELECT ID FROM samples WHERE taxonomy_lock = 0 AND taxonomy_switch IN (1,4)")$ID
    if (length(ids) == 0) {
      showNotification("No unlocked samples are queued for taxonomy assignment.",
        type = "warning"); return()
    }
    shinyWidgets::sendSweetAlert(
      title = "Assigning taxonomy...",
      text = "Hold tight, resolving organism names against the local NCBI taxonomy. This can take a moment.",
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE
    )
    later::later(function() {
      shiny::withReactiveDomain(session, {
        td <- project_taxdump_dir(app_project())
        res <- tryCatch({
          ensure_taxdump(td)
          c2 <- app_con()
          n <- resolve_samples_taxonomy(c2, ids, td)
          DBI::dbDisconnect(c2)
          list(ok = TRUE, n = n)
        }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
        if (isTRUE(res$ok)) {
          refresh(shiny::isolate(refresh()) + 1)
          shinyWidgets::sendSweetAlert(
            title = "Taxonomy assigned",
            text = sprintf("Auto-assigned taxonomy for %d sample(s). Review and lock.", res$n),
            type = "success")
        } else {
          shinyWidgets::sendSweetAlert(title = "Taxonomy assignment failed",
            text = res$msg, type = "error")
        }
      })
    }, delay = 0.1)
  }

  # Parse Nextflow's ANSI progress board (run with -ansi-log). Each redraw is a
  # full snapshot, so key each process row by its name and overwrite in place;
  # keep the latest executor line; banner lines before the first process are the
  # header, and anything after (warnings, completion) is footer.
  nf_proc_rx <- "^\\[(?:[0-9a-f]{2}/[0-9a-f]+|-\\s*)\\]\\s+(\\S+).*$"
  progress_update <- function(lines) {
    lines <- cli::ansi_strip(lines)
    for (ln in lines) {
      if (grepl("^executor\\s*>", ln)) {
        prog_executor(ln)
      } else if (grepl(nf_proc_rx, ln, perl = TRUE)) {
        key <- sub(nf_proc_rx, "\\1", ln, perl = TRUE)
        pl <- prog_proc(); pl[[key]] <- ln; prog_proc(pl)
      } else if (!nzchar(trimws(ln))) {
        next
      } else if (length(prog_proc()) == 0 && !nzchar(prog_executor())) {
        prog_header(unique(c(prog_header(), ln)))
      } else {
        prog_footer(c(prog_footer(), ln))
      }
    }
  }

  observeEvent(input$nf_start, {
    prog_header(character(0)); prog_executor(""); prog_proc(list()); prog_footer(character(0))
    shinyjs::hide("nf_start"); shinyjs::show("nf_stop"); shinyjs::show("nf_progress")
    shinyjs::removeClass(id = "rp_gears", class = "paused")
    # -ansi-log forces the redrawing progress board; wide COLUMNS keeps Nextflow
    # from truncating process names so the parser's per-process keys stay stable.
    p <- processx::process$new("nextflow",
      args = c(nf_command(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE)),
        "-ansi-log"),
      stdout = "|", stderr = "|", wd = app_project(),
      env = c("current", COLUMNS = "500", NXF_ANSI_SUMMARY = "true"))
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

  # One combined block: skip empty sections so there's no blank header/footer
  # space and no inter-element whitespace inside the pre-wrapped code block.
  output$nf_prog_all <- renderText({
    parts <- c(
      paste(prog_header(), collapse = "\n"),
      prog_executor(),
      paste(unlist(prog_proc()), collapse = "\n"),
      paste(utils::tail(prog_footer(), 40), collapse = "\n")
    )
    paste(parts[nzchar(parts)], collapse = "\n")
  })

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
         hit_species, hit_id AS accession, round(score,3) AS score, rank, auto_select
       FROM path_scores WHERE ID = ? ORDER BY rank", params = list(id))
    reactable::reactable(d, compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      defaultColDef = reactable::colDef(align = "left"),
      columns = list(
        path = reactable::colDef(name = "Path", width = 60),
        length = reactable::colDef(name = "Length", width = 80),
        depth = reactable::colDef(name = "Depth", width = 70),
        anchors = reactable::colDef(name = "Anchors", minWidth = 120),
        pident = reactable::colDef(name = "BLAST %ID", width = 90),
        hit_species = reactable::colDef(name = "BLAST Hit", minWidth = 150, html = TRUE, cell = rt_longtext()),
        accession = reactable::colDef(name = "Accession", minWidth = 110, html = TRUE, cell = rt_ncbi_link()),
        score = reactable::colDef(name = "Score", width = 80),
        rank = reactable::colDef(name = "Rank", width = 60),
        auto_select = reactable::colDef(name = "Auto", width = 60)),
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

  # ---- Annotate modal: gene list + boundary editor (edit/save/reset) + msaR ----
  edit_tick <- reactiveVal(0)         # bumps to re-render map/seq/table
  active_feature <- reactiveVal(NULL) # feature the editor + alignment act on
  edit_active <- reactiveVal(FALSE)   # in edit mode?
  edit_backup <- reactiveVal(NULL)    # annotations snapshot for Reset
  edit_dirty <- reactiveVal(FALSE)    # nudged since entering edit mode?

  # reset per-sample modal state; select the first feature + load notes
  observeEvent(sel_id(), ignoreNULL = FALSE, {
    edit_active(FALSE); edit_backup(NULL); edit_dirty(FALSE)
    shinyjs::hide("feat_edit_ctrls")
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
    id <- sel_id(); active_feature(NULL)
    if (!is.null(id)) {
      els <- app_query("SELECT element FROM annotations WHERE ID = ? ORDER BY start", params = list(id))$element
      if (length(els)) active_feature(els[1])
      a <- app_query("SELECT annotate_notes FROM annotate WHERE ID = ?", params = list(id))
      updateTextAreaInput(session, "annotate_notes", value = a$annotate_notes[1] %||% "")
    }
  })

  output$feat_edit_active <- renderText({
    el <- active_feature(); if (is.null(el)) "No active feature" else paste("Active feature:", el)
  })

  output$coverage_plot <- renderPlot({ edit_tick(); plot_coverage(sel_id()) }, height = 128)

  # active feature's stitched-sequence coords + context
  nudge_ctx <- reactive({
    edit_tick()
    id <- sel_id(); el <- active_feature()
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
      return(em("Select a feature (no stitched sequence available)."))
    }
    n <- nchar(cx$seq); flank <- 15
    lo <- max(1, cx$start - flank); hi <- min(n, cx$stop + flank)
    tagList(
      span(sprintf("%s (%s)  %d..%d  len %d bp",
        cx$el, cx$type, cx$start, cx$stop, cx$stop - cx$start + 1),
        style = "display:block; color:#555; margin-bottom:4px;"),
      span(substr(cx$seq, lo, cx$start - 1), style = "color:#aaa;"),
      span(substr(cx$seq, cx$start, cx$stop), style = "color:#0072B2; font-weight:bold;"),
      span(substr(cx$seq, cx$stop + 1, hi), style = "color:#aaa;")
    )
  })

  # ---- reference alignment: active feature vs ALL BLAST-hit references (msaR) ----
  # Rebuilt from the feature's CURRENT boundaries so it redraws as edits move
  # the boundaries; a DECIPHER multiple alignment of the sample + every candidate
  # reference that has the feature.
  active_msa <- reactive({
    edit_tick()
    id <- sel_id(); el <- active_feature(); req(id, el)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    tryCatch(compute_feature_msa(con, id, el), error = function(e) NULL)
  })
  output$msa_header <- renderUI({
    m <- active_msa()
    if (is.null(m)) return(NULL)
    span(sprintf("%s aligned against %d reference hit(s), %d columns",
      active_feature(), length(m) - 1L, Biostrings::width(m)[1]),
      style = "color:#555; display:block; margin-bottom:4px;")
  })
  output$feat_align_ui <- renderUI({
    m <- active_msa()
    if (is.null(m)) {
      return(em(paste("No reference alignment for this feature. Needs GenBank BLAST hits +",
        "reference fetch (WF1); if those ran, re-run WF1.")))
    }
    msaR::msaROutput("feat_msa", height = "auto")
  })
  output$feat_msa <- msaR::renderMsaR({
    m <- active_msa(); req(!is.null(m))
    fa <- tempfile(fileext = ".fa")
    Biostrings::writeXStringSet(m, fa)
    on.exit(unlink(fa), add = TRUE)
    msaR::msaR(fa, overviewbox = FALSE, seqlogo = FALSE, menu = FALSE, conservation = TRUE,
      colorscheme = "nucleotide", labelNameLength = 220, rowheight = 22,
      alignmentHeight = min(260, 30 + 22 * length(m)))
  })

  # ---- gene list (feature model); click a row -> active feature ----
  feature_rows <- reactive({
    edit_tick(); id <- sel_id(); req(id)
    app_query("SELECT element, type, start, stop, length, strand, status, flags FROM annotations WHERE ID = ? ORDER BY start",
      params = list(id))
  })
  output$feature_tbl <- reactable::renderReactable({
    d <- feature_rows(); act <- active_feature()
    if (nrow(d) == 0) {
      return(reactable::reactable(data.frame(` ` = "No annotations for this sample yet.", check.names = FALSE), compact = TRUE))
    }
    tints <- tint_colors(RRNA_ELEMENT_COLORS)
    reactable::reactable(d, compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      defaultColDef = reactable::colDef(align = "left"),
      columns = list(
        element = reactable::colDef(name = "Feature", html = TRUE, cell = rt_link("feat_select")),
        type = reactable::colDef(name = "Type"),
        start = reactable::colDef(name = "Start", width = 70),
        stop = reactable::colDef(name = "Stop", width = 70),
        length = reactable::colDef(name = "Length", width = 70),
        strand = reactable::colDef(name = "Strand", width = 70),
        status = reactable::colDef(name = "Status", width = 70),
        flags = reactable::colDef(name = "Flags", html = TRUE, cell = rt_longtext())),
      rowStyle = function(index) {
        el <- d$element[index]
        st <- list(background = tints[[el]] %||% "#ffffff")   # row tinted by gene color
        if (isTRUE(el == (act %||% ""))) {                     # active: bold + colored left bar
          st$fontWeight <- "bold"
          st$boxShadow <- paste0("inset 4px 0 0 0 ", RRNA_ELEMENT_COLORS[[el]] %||% "#333")
        }
        st
      })
  })
  observeEvent(input$feat_select, {
    i <- suppressWarnings(as.integer(input$feat_select)); d <- feature_rows()
    req(!is.na(i), i >= 1, i <= nrow(d))
    active_feature(d$element[i])
  })

  # ---- curation flags panel ----
  flag_rows <- reactive({
    edit_tick(); id <- sel_id(); req(id)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND status = 'flag' ORDER BY start", params = list(id))
  })
  output$curate_notes <- renderUI({
    id <- sel_id(); req(id)
    a <- app_query("SELECT n_flags, annotate_notes FROM annotate WHERE ID = ?", params = list(id))
    if (nrow(a) == 0 || is.na(a$n_flags[1]) || a$n_flags[1] == 0) {
      return(span(style = "color:#2d6a4f;", HTML("&#10003; No curation flags for this sample.")))
    }
    tagList(
      span(sprintf("%d sample flag(s): ", a$n_flags[1]), style = "color:#7d4a1e; font-weight:bold;"),
      span(a$annotate_notes[1] %||% "", style = "color:#7d4a1e;"))
  })
  output$flag_tbl <- reactable::renderReactable({
    d <- flag_rows()
    if (nrow(d) == 0) {
      return(reactable::reactable(data.frame(` ` = "No flagged features.", check.names = FALSE), compact = TRUE))
    }
    d$jump <- "select"
    reactable::reactable(d[, c("element", "flags", "jump")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      columns = list(
        element = reactable::colDef(name = "Feature"),
        flags = reactable::colDef(name = "Flags", html = TRUE, cell = rt_longtext()),
        jump = reactable::colDef(name = "", width = 70, html = TRUE, cell = rt_link("flag_jump"))))
  })
  observeEvent(input$flag_jump, {
    i <- suppressWarnings(as.integer(input$flag_jump)); d <- flag_rows()
    req(!is.na(i), i >= 1, i <= nrow(d))
    active_feature(d$element[i])
  })

  # ---- edit / save / reset state machine (Edit is a checkbox) ----
  restore_edit_backup <- function() {
    b <- edit_backup()
    if (!is.null(b) && nrow(b)) {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      for (i in seq_len(nrow(b))) {
        DBI::dbExecute(con,
          "UPDATE annotations SET start = ?, stop = ?, length = ? WHERE ID = ? AND element = ?",
          params = list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, sel_id(), b$element[i]))
      }
    }
  }
  save_edits <- function() {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    tryCatch(curate_one(con, sel_id()), error = function(e) NULL)  # refresh flags after edits
    edit_backup(NULL); edit_dirty(FALSE)
    edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    showNotification("Boundary edits saved.", type = "message")
  }
  observeEvent(input$feat_edit, ignoreInit = TRUE, {
    if (isTRUE(input$feat_edit)) {
      if (is.null(active_feature())) {
        shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
        showNotification("Select a feature first.", type = "warning"); return()
      }
      edit_backup(app_query("SELECT element, start, stop FROM annotations WHERE ID = ?", params = list(sel_id())))
      edit_active(TRUE); edit_dirty(FALSE)
      shinyjs::show("feat_edit_ctrls")
    } else {
      edit_active(FALSE); shinyjs::hide("feat_edit_ctrls")
      if (isTRUE(edit_dirty())) {
        # MitoPilot-style prompt: keep or discard unsaved boundary edits
        shinyWidgets::ask_confirmation("unsaved_confirm",
          title = "Unsaved boundary edits",
          text = "You have unsaved changes to this feature's boundaries. Keep them?",
          type = "warning", btn_labels = c("Discard", "Save"),
          btn_colors = c("#d9534f", "#5cb85c"))
      } else {
        edit_backup(NULL)
      }
    }
  })
  observeEvent(input$unsaved_confirm, {
    if (isTRUE(input$unsaved_confirm)) {
      save_edits()
    } else {
      restore_edit_backup(); edit_backup(NULL); edit_dirty(FALSE); edit_tick(edit_tick() + 1)
      showNotification("Boundary edits discarded.", type = "message")
    }
  }, ignoreNULL = TRUE)
  observeEvent(input$feat_save, {
    save_edits()
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)  # clean exit (not dirty)
  })
  observeEvent(input$feat_reset, {
    restore_edit_backup(); edit_dirty(FALSE); edit_tick(edit_tick() + 1)
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
    showNotification("Boundary edits discarded.", type = "message")
  })

  # ---- trim assembly to annotated / custom boundaries ----
  observeEvent(input$trim_btn, {
    id <- sel_id(); req(id)
    s <- app_query("SELECT length FROM stitched WHERE ID = ?", params = list(id))
    if (nrow(s) == 0 || is.na(s$length[1])) {
      showNotification("No stitched assembly to trim.", type = "warning"); return()
    }
    span <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ?", params = list(id))
    a_from <- if (nrow(span) && !is.na(span$a[1])) span$a[1] else 1L
    a_to <- if (nrow(span) && !is.na(span$b[1])) span$b[1] else s$length[1]
    showModal(trim_modal(s$length[1], a_from, a_to))
  })
  observeEvent(input$trim_apply, {
    id <- sel_id(); req(id)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if ((input$trim_mode %||% "annotated") == "annotated") {
      span <- DBI::dbGetQuery(con, "SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ?", params = list(id))
      from <- span$a[1]; to <- span$b[1]
    } else {
      from <- suppressWarnings(as.integer(input$trim_from)); to <- suppressWarnings(as.integer(input$trim_to))
    }
    if (is.na(from) || is.na(to) || from > to) {
      showNotification("Invalid trim range.", type = "error"); return()
    }
    ok <- tryCatch({
      trim_assembly(con, id, from, to)
      tryCatch(curate_one(con, id), error = function(e) NULL)  # re-curate after the shift
      TRUE
    }, error = function(e) { showNotification(conditionMessage(e), type = "error"); FALSE })
    if (isTRUE(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1); removeModal()
      showNotification("Assembly trimmed; feature coordinates + coverage shifted.", type = "message")
    }
  })

  # ---- notes ----
  observeEvent(input$annotate_notes, ignoreInit = TRUE, {
    id <- sel_id(); req(id)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = ? WHERE ID = ?",
      params = list(input$annotate_notes %||% "", id))
  })

  # nudge one boundary by the chosen step (1-50), clamp to sequence, persist
  nudge_boundary <- function(end, dir) {
    if (!isTRUE(edit_active())) { showNotification("Check Edit before adjusting boundaries.", type = "warning"); return() }
    cx <- nudge_ctx(); req(cx)
    step <- suppressWarnings(as.integer(input$feat_step %||% 1))
    if (is.na(step)) step <- 1
    step <- max(1L, min(50L, step))
    n <- nchar(cx$seq); start <- cx$start; stop <- cx$stop
    if (end == "5") start <- start + if (dir == "out") -step else step
    else stop <- stop + if (dir == "out") step else -step
    start <- max(1, start); stop <- min(if (n > 0) n else stop, stop)
    if (start > stop) { showNotification("Boundaries would cross; nudge ignored.", type = "warning"); return() }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND element = ?",
      params = list(start, stop, stop - start + 1, sel_id(), cx$el))
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
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
