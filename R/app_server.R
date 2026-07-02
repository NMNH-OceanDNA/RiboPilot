#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  refresh <- reactiveVal(0)
  observeEvent(input$refresh, refresh(refresh() + 1))
  sel_id <- reactiveVal(NULL)

  output$proj_label <- renderText(basename(app_project()))

  # ---- samples table data (per mode; includes state + lock) ----
  samples_data <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    if (mode == "Assemble") {
      app_query(
        "SELECT s.ID, s.Taxon, a.assemble_lock, a.assemble_switch, a.paths, a.scaffolds, a.length,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           a.assemble_notes AS notes
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID")
    } else if (mode == "Export") {
      app_query(
        "SELECT s.ID, s.Taxon, an.reviewed, an.n_features, an.missing, st.length,
           st.reference_accession AS reference
         FROM samples s LEFT JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID")
    } else {
      app_query(
        "SELECT s.ID, s.Taxon, an.annotate_lock, an.annotate_switch, an.reviewed, an.n_features,
           an.missing, st.reference_accession AS reference
         FROM samples s LEFT JOIN annotate an ON s.ID = an.ID LEFT JOIN stitched st ON s.ID = st.ID")
    }
  })

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    d <- samples_data()
    lock_col <- if (mode == "Annotate & Curate") "annotate_lock" else "assemble_lock"
    state_col <- if (mode == "Annotate & Curate") "annotate_switch" else "assemble_switch"

    cols <- list(
      ID = reactable::colDef(minWidth = 130, sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      Taxon = reactable::colDef(minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      reviewed = reactable::colDef(name = "reviewed", html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      reference = reactable::colDef(html = TRUE, cell = rt_ncbi_link()),
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
    } else if (mode == "Annotate & Curate") {
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
  observeEvent(input$lock_btn, {
    col <- if ((input$mode %||% "") == "Annotate & Curate") "annotate" else "assemble"
    set_flag(sprintf("UPDATE %s SET %s_lock = 1 WHERE ID = ?", col, col))
  })
  observeEvent(input$unlock_btn, {
    col <- if ((input$mode %||% "") == "Annotate & Curate") "annotate" else "assemble"
    set_flag(sprintf("UPDATE %s SET %s_lock = 0 WHERE ID = ?", col, col))
  })
  observeEvent(input$requeue_btn, set_flag("UPDATE assemble SET assemble_switch = 1 WHERE ID = ?"))
  observeEvent(input$reviewed_btn, set_flag("UPDATE annotate SET reviewed = 'yes' WHERE ID = ?"))

  # ---- Export ----
  output$export_sel_n <- renderText({
    n <- length(selected_ids())
    if (n == 0) "all samples (none checked)" else paste(n, "selected")
  })
  output$export_dl <- downloadHandler(
    filename = function() paste0(basename(app_project()), "_export.zip"),
    content = function(file) {
      ids <- selected_ids()
      td <- file.path(tempdir(), paste0("rp_export_", as.integer(Sys.time())))
      export_bundle(app_project(), ids = if (length(ids)) ids else NULL, out_dir = td, zip = TRUE)
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
    refresh(refresh() + 1)
    showNotification(paste0("Set path ", input$sel_path,
      " as selected. Re-run annotation (run_annotate) to update the feature model."), type = "message")
  })

  # ---- Annotate modal ----
  output$feature_map <- renderPlot(plot_feature_map(sel_id()), height = 150)
  output$annot_edit <- renderUI({
    id <- sel_id()
    req(id)
    d <- app_query("SELECT element, type, start, stop FROM annotations WHERE ID = ? ORDER BY start",
      params = list(id))
    if (nrow(d) == 0) {
      return(em("No annotations for this sample."))
    }
    do.call(tagList, lapply(seq_len(nrow(d)), function(i) {
      el <- d$element[i]
      div(style = "display:flex; gap:0.75em; align-items:center; margin-bottom:2px;",
        strong(el, style = "width:55px;"),
        div(style = "width:120px;", numericInput(paste0("start_", el), "start", d$start[i])),
        div(style = "width:120px;", numericInput(paste0("stop_", el), "stop", d$stop[i])),
        span(paste0("(", d$type[i], ")"), style = "color:#888;"))
    }))
  })
  output$feature_tbl <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    reactable::reactable(
      app_query("SELECT element, type, start, stop, length, strand, tool FROM annotations WHERE ID = ? ORDER BY start",
        params = list(id)),
      compact = TRUE, resizable = TRUE, defaultPageSize = 10)
  })
  observeEvent(input$save_annot, {
    id <- sel_id()
    req(id)
    els <- app_query("SELECT element FROM annotations WHERE ID = ?", params = list(id))$element
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    for (el in els) {
      s <- input[[paste0("start_", el)]]
      e <- input[[paste0("stop_", el)]]
      if (is.null(s) || is.null(e) || is.na(s) || is.na(e)) next
      DBI::dbExecute(con,
        "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND element = ?",
        params = list(as.integer(s), as.integer(e), as.integer(e - s + 1), id, el))
    }
    refresh(refresh() + 1)
    showNotification("Saved annotation edits.", type = "message")
  })
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
