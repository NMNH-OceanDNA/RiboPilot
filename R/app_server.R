#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  refresh <- reactiveVal(0)
  observeEvent(input$refresh, refresh(refresh() + 1))
  sel_id <- reactiveVal(NULL)

  output$proj_label <- renderText(basename(app_project()))

  # ---- samples table (columns per mode) ----
  samples_data <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    if (mode == "Assemble") {
      app_query(
        "SELECT s.ID, s.Taxon, a.paths, a.scaffolds, a.length, a.assemble_notes AS notes,
           (SELECT COUNT(*) FROM path_scores p WHERE p.ID = s.ID) AS scored_paths
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID")
    } else if (mode == "Export") {
      app_query(
        "SELECT s.ID, s.Taxon, an.n_features, an.missing, an.reviewed, st.length AS length,
           st.reference_accession AS reference
         FROM samples s LEFT JOIN annotate an ON s.ID = an.ID
         LEFT JOIN stitched st ON s.ID = st.ID")
    } else {
      app_query(
        "SELECT s.ID, s.Taxon, an.n_features, an.missing, an.reviewed,
           st.reference_accession AS reference
         FROM samples s LEFT JOIN annotate an ON s.ID = an.ID
         LEFT JOIN stitched st ON s.ID = st.ID")
    }
  })

  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    reactable::reactable(
      samples_data(),
      selection = if (mode == "Export") "multiple" else "single",
      onClick = "select", highlight = TRUE,
      defaultPageSize = 15, compact = TRUE, wrap = FALSE,
      theme = reactable::reactableTheme(rowSelectedStyle = list(backgroundColor = "#eef6ff"))
    )
  })

  # ---- Export mode: multi-select + bundle download ----
  output$export_sel_n <- renderText({
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel)) "all samples (none checked)" else paste(length(sel), "selected")
  })
  output$export_dl <- downloadHandler(
    filename = function() paste0(basename(app_project()), "_export.zip"),
    content = function(file) {
      sel <- reactable::getReactableState("tbl", "selected")
      ids <- if (is.null(sel)) NULL else samples_data()$ID[sel]
      td <- file.path(tempdir(), paste0("rp_export_", as.integer(Sys.time())))
      export_bundle(app_project(), ids = ids, out_dir = td, zip = TRUE)
      file.copy(paste0(td, ".zip"), file, overwrite = TRUE)
    }
  )

  # ---- row click opens a modal (single-select modes only) ----
  observeEvent(reactable::getReactableState("tbl", "selected"), {
    mode <- input$mode %||% "Assemble"
    if (mode == "Export") {
      return()
    }
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel) || length(sel) != 1) {
      return()
    }
    id <- samples_data()$ID[sel]
    sel_id(id)
    showModal(if (mode == "Assemble") assemble_modal(id) else annotate_modal(id))
  }, ignoreInit = TRUE)

  # ---- Assemble modal ----
  output$sel_path_ui <- renderUI({
    id <- sel_id()
    req(id)
    d <- app_query("SELECT path, auto_select FROM path_scores WHERE ID = ? ORDER BY path", params = list(id))
    if (nrow(d) == 0) {
      return(em("No scored paths."))
    }
    cur <- d$path[d$auto_select == 1][1]
    radioButtons("sel_path", "Selected path:", choices = d$path, selected = cur, inline = TRUE)
  })
  output$assemble_scores <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    d <- app_query(
      "SELECT path, length, round(depth,0) AS depth, anchors, blast_pident AS pident,
         hit_species, round(score,3) AS score, rank, auto_select
       FROM path_scores WHERE ID = ? ORDER BY rank", params = list(id))
    reactable::reactable(d, compact = TRUE, defaultPageSize = 10,
      rowStyle = function(index) if (isTRUE(d$auto_select[index] == 1)) list(background = "#eafaf0") else NULL)
  })
  output$assemble_scaffolds <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    reactable::reactable(
      app_query("SELECT path, scaffold, length, topology FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
        params = list(id)),
      compact = TRUE, defaultPageSize = 10)
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
      " as selected. Re-run annotation to update the feature model."), type = "message")
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
    rows <- lapply(seq_len(nrow(d)), function(i) {
      el <- d$element[i]
      div(style = "display:flex; gap:0.75em; align-items:center; margin-bottom:2px;",
        strong(el, style = "width:55px;"),
        div(style = "width:120px;", numericInput(paste0("start_", el), "start", d$start[i])),
        div(style = "width:120px;", numericInput(paste0("stop_", el), "stop", d$stop[i])),
        span(paste0("(", d$type[i], ")"), style = "color:#888;"))
    })
    do.call(tagList, rows)
  })
  output$feature_tbl <- reactable::renderReactable({
    id <- sel_id()
    req(id)
    reactable::reactable(
      app_query("SELECT element, type, start, stop, length, strand, tool FROM annotations WHERE ID = ? ORDER BY start",
        params = list(id)),
      compact = TRUE, defaultPageSize = 10)
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
