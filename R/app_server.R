#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  refresh <- shiny::reactiveVal(0)
  shiny::observeEvent(input$refresh, refresh(refresh() + 1))
  sel_id <- shiny::reactiveVal(NULL)

  output$proj_label <- shiny::renderText(basename(app_project()))

  # ---- samples table (columns per mode) ----
  samples_data <- shiny::reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    if (mode == "Assemble") {
      app_query(
        "SELECT s.ID, s.Taxon, a.paths, a.scaffolds, a.length, a.assemble_notes AS notes,
           (SELECT COUNT(*) FROM path_scores p WHERE p.ID = s.ID) AS scored_paths
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID")
    } else if (mode == "Export") {
      app_query(
        "SELECT s.ID, s.Taxon, an.n_features, an.missing, st.length AS stitched_len,
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
    reactable::reactable(
      samples_data(),
      selection = "single", onClick = "select", highlight = TRUE,
      defaultPageSize = 15, compact = TRUE, wrap = FALSE,
      theme = reactable::reactableTheme(
        rowSelectedStyle = list(backgroundColor = "#eef6ff")
      )
    )
  })

  # row click -> open the mode-specific modal
  shiny::observeEvent(reactable::getReactableState("tbl", "selected"), {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel)) {
      return()
    }
    id <- samples_data()$ID[sel]
    sel_id(id)
    mode <- input$mode %||% "Assemble"
    shiny::showModal(switch(mode,
      "Assemble" = assemble_modal(id),
      "Export" = export_modal(id),
      annotate_modal(id)
    ))
  }, ignoreInit = TRUE)

  # ---- Assemble modal content ----
  output$assemble_scores <- reactable::renderReactable({
    id <- sel_id()
    shiny::req(id)
    d <- app_query(
      "SELECT path, length, round(depth,0) AS depth, anchors, blast_pident AS pident,
         hit_species, round(score,3) AS score, rank, auto_select
       FROM path_scores WHERE ID = ? ORDER BY rank", params = list(id))
    reactable::reactable(d, compact = TRUE, defaultPageSize = 10,
      rowStyle = function(index) if (isTRUE(d$auto_select[index] == 1)) list(background = "#eafaf0") else NULL)
  })
  output$assemble_scaffolds <- reactable::renderReactable({
    id <- sel_id()
    shiny::req(id)
    reactable::reactable(
      app_query("SELECT path, scaffold, length, topology FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
        params = list(id)),
      compact = TRUE, defaultPageSize = 10)
  })

  # ---- Annotate modal content ----
  output$feature_map <- shiny::renderPlot(plot_feature_map(sel_id()), height = 160)
  output$feature_tbl <- reactable::renderReactable({
    id <- sel_id()
    shiny::req(id)
    reactable::reactable(
      app_query("SELECT element, type, start, stop, length, strand, tool FROM annotations WHERE ID = ? ORDER BY start",
        params = list(id)),
      compact = TRUE, defaultPageSize = 10)
  })

  # ---- Export modal content ----
  output$export_tbl_preview <- shiny::renderText(build_feature_table(sel_id()))
  output$download_tbl <- shiny::downloadHandler(
    filename = function() paste0(sel_id(), ".tbl"),
    content = function(file) writeLines(build_feature_table(sel_id()), file)
  )
  output$download_fsa <- shiny::downloadHandler(
    filename = function() paste0(sel_id(), ".fsa"),
    content = function(file) {
      id <- sel_id()
      seq <- app_query("SELECT sequence FROM stitched WHERE ID = ?", params = list(id))$sequence
      writeLines(c(paste0(">", id), seq %||% ""), file)
    }
  )
}
