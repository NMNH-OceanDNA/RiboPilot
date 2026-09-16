#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # Stop the app when the last browser session closes (survives refresh).
  register_app_lifecycle(session)

  refresh <- reactiveVal(0)
  observeEvent(input$refresh, refresh(refresh() + 1))
  sel_id <- reactiveVal(NULL)
  sel_unit <- reactiveVal(NULL)  # list(ID, path, scaffold) for the open annotate modal
  # (ID, path, scaffold) of samples_data() row i; path/scaffold default to 1 when
  # the active table has no unit columns (Assemble / Taxonomy modes).
  row_unit <- function(i) {
    d <- samples_data()
    list(ID = d$ID[i],
      path = if ("path" %in% names(d)) as.integer(d$path[i]) else 1L,
      scaffold = if ("scaffold" %in% names(d)) as.integer(d$scaffold[i]) else 1L)
  }
  # positional params (ID, path, scaffold) of the open annotate unit
  uargs <- function() { u <- sel_unit(); list(u$ID, u$path, u$scaffold) }

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
      footer = mp_footer(dismiss = "Close")))
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
    # Never drop unreachable dirs (a failed run's dir may be cleaned up); mark
    # them and show the path so it is still copyable, and only offer Open when the
    # dir exists on this host.
    d$open <- ifelse(d$exists, "open", "")
    d$on_disk <- ifelse(d$exists, "yes", "missing")
    reactable::reactable(d[, c("process", "status", "on_disk", "modified", "open", "workdir")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 25, wrap = FALSE,
      columns = list(
        status = reactable::colDef(html = TRUE, width = 80, cell = rt_dynamicIcon(
          c(success = "fa fa-circle-check mp-fg-success", failed = "fa fa-triangle-exclamation mp-fg-danger"),
          labels = c(success = "Success", failed = "Failed"))),
        on_disk = reactable::colDef(name = "On disk", width = 80),
        open = reactable::colDef(name = "", width = 60, html = TRUE, cell = rt_link("wd_open")),
        workdir = reactable::colDef(name = "Work directory (select to copy)",
          minWidth = 340, html = TRUE, cell = rt_longtext())))
  })
  observeEvent(input$wd_open, {
    i <- suppressWarnings(as.integer(input$wd_open))
    d <- wd_rows()
    req(!is.na(i), i >= 1, i <= nrow(d))
    if (!isTRUE(d$exists[i])) {
      rp_alert("Work directory not found on this host.", type = "warning")
      return()
    }
    if (!open_path(d$workdir[i])) {
      rp_alert(paste0("This session has no file browser to open. The work directory is: ",
        d$workdir[i]), type = "warning")
    }
  })

  output$app_header <- renderUI(div(class = "mp-header",
    span(class = "mp-app-name", "RiboPilot"),
    span(class = "mp-app-ver", paste0("v", utils::packageVersion("RiboPilot"))),
    span(class = "mp-proj", title = app_project(), basename(app_project())),
    a(class = "mp-help", href = "https://github.com/dmacguigan/RiboPilot",
      target = "_blank", rel = "noopener", "Help")))
  # row-count line above the table + CSV export of the table
  output$table_status <- renderUI({
    d <- samples_data()
    mode <- input$mode %||% "Taxonomy"
    noun <- if (mode == "Reference" && ref_view() == "Groups") "group" else MP_UNIT[[tolower(mode)]]
    # the client counts what reactable shows (column filters, search, and the
    # CSS-hidden lock/state rows); the server mask is the fallback before it reports
    shown <- input$tbl_visible_n %||% sum(visible_mask())
    mp_table_status(min(shown, nrow(d)), nrow(d), length(selected_ids()), noun)
  })
  observe(shinyjs::toggleState(selector = ".mp-needs-selection",
    condition = length(selected_ids()) > 0))
  table_csv_handler <- function(suffix, rows) downloadHandler(
    filename = function() paste0(basename(app_project()), "_",
      tolower(input$mode %||% "table"), suffix, ".csv"),
    content = function(file) {
      d <- samples_data()
      d <- d[rows(), setdiff(names(d), c("details", "output", "action")), drop = FALSE]
      utils::write.csv(d, file, row.names = FALSE)
    }
  )
  output$table_csv <- table_csv_handler("", function() TRUE)
  output$table_csv_selected <- table_csv_handler("_selected", selected_rows)

  # mode -> (table, prefix) for lock/state columns and updates
  mode_meta <- function(mode) {
    switch(mode %||% "Assemble",
      Taxonomy = list(tbl = "samples", pfx = "taxonomy"),
      Reference = list(tbl = "reference", pfx = "reference"),
      Annotate = list(tbl = "annotate", pfx = "annotate"),
      list(tbl = "assemble", pfx = "assemble"))
  }

  # Which of the Reference tab's two views is showing.
  ref_view <- reactive(input$ref_view %||% "Samples")

  # ---- samples table data (per mode; includes state + lock) ----
  samples_data <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    # Column order mirrors MitoPilot: lock, state, ID, Taxon lead the tables.
    if (mode == "Taxonomy") {
      app_query(
        "SELECT s.taxonomy_lock, s.taxonomy_switch, s.ID, s.Taxon,
           s.sci_name, s.taxid, s.rank, s.match_type, s.time_stamp
         FROM samples s")
    } else if (mode == "Reference") {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      if (ref_view() == "Samples") reference_sample_table(con) else reference_table(con)
    } else if (mode == "Assemble") {
      app_query(
        "SELECT a.assemble_lock, a.assemble_switch, s.ID, s.Taxon, s.sci_name, s.taxid,
           pp.pre_opts, a.assemble_opts, a.blast_opts,
           -- annotation + curation presets are set here too: one workflow runs
           -- assembly and annotation, so every option a run consumes belongs on
           -- the surface the run is launched from. Taken from the sample's first
           -- unit; opts_apply writes them back across all of the sample's units.
           (SELECT an.annotate_opts FROM annotate an WHERE an.ID = s.ID
              ORDER BY an.path, an.scaffold LIMIT 1) AS annotate_opts,
           (SELECT an.curate_opts FROM annotate an WHERE an.ID = s.ID
              ORDER BY an.path, an.scaffold LIMIT 1) AS curate_opts,
           a.paths, a.length,
           -- unit counts come from `assemblies` so they stay live as scaffolds
           -- are kept/dropped. a.blast_* is not selected: nothing writes it (the
           -- per-unit hits live in `assemblies`), so those columns were blank.
           (SELECT COUNT(*) FROM assemblies x WHERE x.ID = s.ID) AS units,
           (SELECT COUNT(*) FROM assemblies x WHERE x.ID = s.ID
              AND COALESCE(x.ignore,0) = 0) AS kept,
           -- per-sample roll-up of the stock-ITSx check: NULL until any unit
           -- has been checked, else 'n of m units' flagged
           (SELECT CASE WHEN COUNT(x.itsx_warnings) = 0 THEN NULL ELSE
              SUM(CASE WHEN COALESCE(x.itsx_warnings,'') <> ''
                       AND x.itsx_warnings NOT LIKE 'not_run(%' THEN 1 ELSE 0 END)
              || ' of ' || COUNT(x.itsx_warnings) || ' units' END
            FROM annotate x WHERE x.ID = s.ID) AS itsx_warnings,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           -- the rDNA contig search writes its own column because
           -- load_assemblies() overwrites assemble_notes right after it; both
           -- are shown as one Notes cell
           TRIM(COALESCE(a.assemble_notes,'') ||
             CASE WHEN COALESCE(a.find_rrna_notes,'') = '' THEN ''
                  WHEN COALESCE(a.assemble_notes,'') = '' THEN a.find_rrna_notes
                  ELSE '; ' || a.find_rrna_notes END) AS notes, a.time_stamp,
           'details' AS details, 'output' AS output
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID
           LEFT JOIN preprocess pp ON s.ID = pp.ID
         WHERE s.taxonomy_lock = 1")
    } else if (mode == "Export") {
      # export shows annotation-locked UNITS (one row per (ID, path, scaffold))
      app_query(
        "SELECT s.ID, an.path, an.scaffold, s.Taxon, s.sci_name, s.taxid, asm.length,
           an.n_features, an.missing,
           an.n_flags, an.annotate_notes, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason,
           an.reviewed, an.problematic, an.partial,
           an.export_group, an.time_stamp
         FROM annotate an JOIN samples s ON s.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE an.annotate_lock = 1")
    } else {
      # annotate shows one row per non-ignored assembly UNIT of an assembly-locked
      # sample (one (ID, path, scaffold) each, annotated independently)
      app_query(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, an.path, an.scaffold, s.Taxon, s.sci_name, s.taxid,
           asm.length, an.n_features, an.missing,
           -- blast_source says WHICH search named the species: the local seed
           -- BLAST, the local panel BLAST, or the remote core_nt one. Without
           -- it a seed or panel hit reads as a GenBank identification.
           -- blast_qcovs is written by the remote and local-panel paths only
           -- and its column was already defined but never selected.
           asm.blast_species, asm.blast_accession, asm.blast_pident, asm.blast_qcovs,
           asm.blast_source,
           an.n_flags, an.annotate_notes, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason, an.reviewed,
           an.problematic, an.partial, an.time_stamp,
           'details' AS details, 'output' AS output
         FROM annotate an
           JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0")
    }
  })

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  # Re-render only on structural changes (mode / column groups); the data is
  # isolated so a refresh() pushes a data-only update below instead of a full
  # re-render, which would reset the user's column filters, sort, and page.
  # Column defs read the data once; these two re-render the table only when a
  # refresh() changes what the defs depend on (unit columns, group members).
  multi_unit_rv <- reactiveVal(FALSE)
  observe({ d <- samples_data(); multi_unit_rv(any(duplicated(d$ID))) })
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    multi_unit_rv()
    if (mode == "Reference" && ref_view() == "Groups") refresh()
    d <- isolate(samples_data())
    mm <- mode_meta(mode)
    session$sendCustomMessage("rp_table_pfx", mm$pfx)
    lock_col <- paste0(mm$pfx, "_lock")
    state_col <- paste0(mm$pfx, "_switch")

    # Details/Output icon buttons open the mode's detail modal + the work dirs.
    det_input <- if (mode == "Annotate") "ann_details" else "asm_details"
    det_tip <- if (mode == "Annotate") "View annotation details" else "View assembly details"
    out_input <- if (mode == "Annotate") "ann_output" else "asm_output"

    # Path/Scaffold are shown only when the sample has more than one unit; a
    # single-unit project reads exactly like the old per-sample table.
    multi_unit <- all(c("path", "scaffold") %in% names(d)) &&
      nrow(unique(d[, c("ID", "path", "scaffold")])) > length(unique(d$ID))

    cd <- function(f, ...) reactable::colDef(name = mp_col_name(f), header = mp_col_header(f), ...)
    cols <- list(
      ID = cd("ID", width = mp_fit_width(d$ID), sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      path = cd("path", width = 70, align = "center", show = multi_unit),
      scaffold = cd("scaffold", width = 90, align = "center", show = multi_unit),
      Taxon = cd("Taxon", minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      sci_name = cd("sci_name", minWidth = 200, html = TRUE, cell = rt_taxid_link()),
      taxid = reactable::colDef(show = FALSE),
      rank = cd("rank", minWidth = 90, html = TRUE, cell = rt_longtext()),
      match_type = cd("match_type", minWidth = 110, html = TRUE, align = "center",
        cell = rt_confidence_badge()),
      reviewed = cd("reviewed", html = TRUE, align = "center", cell = rt_bool_badge(hide_no = TRUE)),
      problematic = cd("problematic", width = 110, html = TRUE, align = "center",
        cell = rt_bool_badge(hide_no = TRUE, yes_tone = "warning")),
      partial = cd("partial", width = 80, html = TRUE, align = "center",
        cell = rt_bool_badge(hide_no = TRUE, yes_tone = "warning")),
      n_flags = cd("n_flags", width = 100, html = TRUE, filterable = FALSE, align = "center",
        cell = rt_flag_badge()),
      annotate_notes = reactable::colDef(show = FALSE),
      itsx_warnings = cd("itsx_warnings", minWidth = 130, html = TRUE, cell = rt_itsx_badge()),
      screen_verdict = cd("screen_verdict", width = 110, html = TRUE, align = "center",
        cell = rt_screen_badge()),
      chimera_label = cd("chimera_label", width = 100, html = TRUE, align = "center",
        cell = rt_chimera_badge()),
      chimera_reason = reactable::colDef(show = FALSE),
      reference = cd("reference", html = TRUE, cell = rt_ncbi_link()),
      blast_species = cd("blast_species", minWidth = 160, html = TRUE, cell = rt_longtext()),
      blast_accession = cd("blast_accession", minWidth = 130, html = TRUE, cell = rt_ncbi_link()),
      blast_pident = cd("blast_pident", width = 130),
      blast_qcovs = cd("blast_qcovs", width = 140),
      blast_source = cd("blast_source", width = 100, html = TRUE, align = "center",
        cell = rt_blast_source_badge()),
      details = cd("details", width = 80, align = "center", filterable = FALSE, sortable = FALSE,
        sticky = "right", class = "mp-actions-sticky", headerClass = "mp-actions-sticky", html = TRUE,
        cell = rt_icon_bttn_text(det_input, "fas fa-square-arrow-up-right fa-xs", title = det_tip)),
      output = cd("output", width = 80, align = "center", filterable = FALSE, sortable = FALSE,
        sticky = "right", class = "mp-actions-sticky", headerClass = "mp-actions-sticky", html = TRUE,
        cell = rt_icon_bttn_text(out_input, "fas fa-folder-open fa-xs", title = "Work / output directories")),
      n_features = cd("n_features", width = 90, html = TRUE),
      missing = cd("missing", minWidth = 120, html = TRUE, cell = rt_longtext()),
      length = cd("length", width = 100),
      paths = cd("paths", width = 80),
      units = cd("units", width = 70),
      kept = cd("kept", width = 70, html = TRUE, cell = rt_kept_count()),
      sel_path = cd("sel_path", width = 110),
      best_score = cd("best_score", width = 100),
      export_group = cd("export_group", minWidth = 110, html = TRUE, cell = rt_longtext()),
      pre_opts = cd("pre_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_pre_opts", lock_col = "assemble_lock")),
      assemble_opts = cd("assemble_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_assemble_opts", lock_col = "assemble_lock")),
      blast_opts = cd("blast_opts", minWidth = 100, html = TRUE, filterable = FALSE,
        cell = rt_link("set_blast_opts", lock_col = "assemble_lock")),
      annotate_opts = cd("annotate_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_annotate_opts", lock_col = "assemble_lock")),
      curate_opts = cd("curate_opts", minWidth = 100, html = TRUE, filterable = FALSE,
        cell = rt_link("set_curate_opts", lock_col = "assemble_lock")),
      annotate_lock = cd("annotate_lock", width = 34, align = "center", filterable = FALSE,
        html = TRUE, cell = rt_lock_icon()),
      notes = cd("notes", minWidth = 160, html = TRUE, class = "mp-note-cell", cell = rt_longtext()),
      time_stamp = cd("time_stamp", width = 150, filterable = FALSE, cell = rt_ts_date())
    )
    if (lock_col %in% names(d)) {
      cols[[lock_col]] <- cd(lock_col, sticky = "left", html = TRUE, width = 52,
        align = "center", filterable = FALSE, cell = rt_lock_icon())
    }
    if (state_col %in% names(d) && mm$pfx %in% names(MP_STATE_CODES)) {
      cols[[state_col]] <- cd(state_col, sticky = "left", html = TRUE, width = 54,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(icons = mp_state_icon_classes(mm$pfx), labels = mp_state_titles(mm$pfx)))
    }
    if (mode == "Export") {
      cols[["ID"]] <- cd("ID", sticky = "left", width = mp_fit_width(d$ID))
    }
    if (mode == "Reference" && ref_view() == "Samples") {
      cols[["ID"]] <- cd("ID", sticky = "left", width = mp_fit_width(d$ID))
      cols[["Taxon"]] <- cd("Taxon", minWidth = 160)
      cols[["sci_name"]] <- cd("sci_name", minWidth = 160)
      cols[["ref_group"]] <- cd("ref_group", minWidth = 160)
      cols[["clade_rank"]] <- cd("clade_rank", width = 90)
      cols[["user_group"]] <- cd("user_group", width = 110)
      cols[["set_by"]] <- cd("set_by", width = 90)
      cols[[state_col]] <- reactable::colDef(show = FALSE)
    } else if (mode == "Reference") {
      # The group key is the clade name lowercased (ref_id_from_name), so showing
      # both is showing the same thing twice; the key stays in the data for the
      # selection, the CSV and the row classes, and stays visible in the database
      # path, which is named after it.
      cols[["ID"]] <- reactable::colDef(show = FALSE)
      cols[["status"]] <- cd("status", sticky = "left", minWidth = 150)
      cols[["clade_name"]] <- cd("clade_name", minWidth = 160, sticky = "left",
        html = TRUE, cell = rt_link("ref_open_clade"))
      cols[["clade_rank"]] <- cd("clade_rank", width = 90)
      cols[["n_samples"]] <- cd("n_samples", width = 100)
      cols[["n_records"]] <- cd("n_records", width = 160)
      cols[["database"]] <- cd("database", minWidth = 240, html = TRUE, cell = rt_longtext())
      cols[["fail_reason"]] <- cd("fail_reason", minWidth = 200, html = TRUE,
        cell = rt_link("ref_open_log"))
      # the row's own verb, named for the state it is in
      cols[["action"]] <- cd("action", minWidth = 150, sticky = "right", align = "center",
        class = "mp-actions-sticky", headerClass = "mp-actions-sticky",
        html = TRUE, filterable = FALSE, cell = rt_link("ref_row_action"))
      # the state still drives the row class, so it stays in the data, unshown
      cols[[state_col]] <- reactable::colDef(show = FALSE)
    }
    rowclass <- reactable::JS(sprintf(
      "function(rowInfo){ if(!rowInfo||!rowInfo.values) return '';
         return 'rp-lock-'+rowInfo.values['%s']+' rp-state-'+rowInfo.values['%s']; }",
      lock_col, state_col))

    # Explain an empty table: every mode after Taxonomy is gated on the previous
    # stage's lock, and a fresh project has nothing locked, so the default
    # "No rows found" reads like a failed import.
    nodata <- switch(mode,
      Taxonomy = "No samples yet. Create the project with a mapping file (new_project()).",
      Reference = if (ref_view() == "Samples") {
        "No samples yet. Create the project with a mapping file (new_project())."
      } else {
        paste("No reference groups yet. Until a group has its own database, its samples",
          "assemble against RiboPilot's built-in ray-finned fish reference.",
          "Press Propose groups to group your samples by clade from the assigned taxonomy.")
      },
      Assemble = "No taxonomy-locked samples. Assign and lock taxonomy in Taxonomy mode first.",
      Annotate = "No assembly-locked samples. Review and lock assemblies in Assemble mode first.",
      Export   = "No annotation-locked samples. Review and lock annotations in Annotation Review first.",
      "No rows.")

    sorted <- if ("time_stamp" %in% names(d)) list(time_stamp = "desc") else NULL

    # Column groups (Assemble table): collapse wide sections via the "Columns"
    # picker. Columns not in a named group are Core and always shown.
    col_groups <- switch(mode,
      Assemble = list(
        Stats = c("paths", "units", "kept", "sel_path", "best_score"),
        Options = c("pre_opts", "assemble_opts", "blast_opts", "annotate_opts", "curate_opts")),
      Annotate = list(
        BLAST = c("blast_species", "blast_accession", "blast_pident", "blast_qcovs", "blast_source"),
        Checks = c("n_flags", "itsx_warnings", "screen_verdict", "chimera_label")),
      NULL)
    keep_cols <- cols[intersect(names(cols), names(d))]
    keep_cols[[".selection"]] <- reactable::colDef(sticky = "left")
    if (!is.null(col_groups)) {
      sel_g <- input[[if (mode == "Assemble") "col_groups" else "ann_col_groups"]] %||% names(col_groups)
      hidden <- unlist(col_groups[setdiff(names(col_groups), sel_g)], use.names = FALSE)
      for (cn in intersect(hidden, names(d))) keep_cols[[cn]] <- reactable::colDef(show = FALSE)
    }

    # A group row expands to the samples it holds: the count alone never said
    # which sample ended up where, which is the whole reason the tab was opaque.
    ref_details <- NULL
    if (mode == "Reference" && ref_view() == "Groups") {
      con <- app_con(); mem <- reference_members(con); DBI::dbDisconnect(con)
      ref_details <- function(index) {
        reference_member_html(mem[mem$ref_id == d$ID[index], , drop = FALSE])
      }
    }

    reactable::reactable(
      d,
      details = ref_details,
      resizable = TRUE, filterable = TRUE, striped = TRUE, compact = TRUE,
      searchable = TRUE, showPageSizeOptions = TRUE, wrap = FALSE, height = "100%",
      selection = "multiple", onClick = "select",
      defaultPageSize = 50, pageSizeOptions = c(25, 50, 100, 200),
      defaultSorted = sorted,
      language = reactable::reactableLang(noData = nodata),
      # consistent header capitalization (Title Case names above) + left alignment
      defaultColDef = reactable::colDef(align = "left", headerVAlign = "center"),
      rowStyle = rt_highlight_row(), rowClass = rowclass,
      # column order comes from samples_data()
      columns = keep_cols
    )
  })

  # Push data-only updates on refresh() so client-side filters / sort / page
  # survive a lock, state, or group action (a full re-render would clear them).
  # Mode and column-group changes still re-render via output$tbl above.
  observeEvent(refresh(), {
    reactable::updateReactable("tbl", data = samples_data())
  }, ignoreInit = TRUE)

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
        rp_alert("Lock the assembly (Assemble mode) before annotating this sample.",
          type = "warning")
        return()
      }
      u <- row_unit(i); sel_unit(u)
      showModal(annotate_modal(u$ID, u$path, u$scaffold))
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

  # ---- Annotate table "Details" / "Output" buttons ----
  # Rows here are already assembly-locked (table WHERE assemble_lock = 1), so the
  # detail button can open the annotate modal without re-checking the lock.
  observeEvent(input$ann_details, {
    i <- suppressWarnings(as.integer(input$ann_details))
    req(!is.na(i))
    u <- row_unit(i); sel_id(u$ID); sel_unit(u)
    showModal(annotate_modal(u$ID, u$path, u$scaffold))
  })
  observeEvent(input$ann_output, {
    i <- suppressWarnings(as.integer(input$ann_output))
    req(!is.na(i))
    show_workdirs(samples_data()$ID[i])
  })

  # ---- sample navigation in the detail modals (prev / picker / next) ----
  # Assemble rows have no path/scaffold; the modal script reports "1" for both.
  nav_key <- function(d) paste(d$ID, d$path, d$scaffold, sep = "|")
  nav_order <- reactive({
    o <- input$nav_order
    if (length(o$ID) == 0) return(NULL)
    data.frame(ID = as.character(o$ID), path = as.character(o$path),
      scaffold = as.character(o$scaffold), stringsAsFactors = FALSE)
  })
  nav_current <- function() {
    if ((input$mode %||% "Assemble") == "Annotate") {
      u <- sel_unit(); if (is.null(u)) NULL else nav_key(u)
    } else {
      id <- sel_id(); if (is.null(id)) NULL else paste(id, 1, 1, sep = "|")
    }
  }
  output$sample_nav <- renderUI({
    o <- nav_order(); cur <- nav_current(); req(o, cur)
    keys <- nav_key(o); i <- match(cur, keys); req(!is.na(i))
    labels <- if (any(duplicated(o$ID))) paste0(o$ID, " [", o$path, "/", o$scaffold, "]") else o$ID
    arrow <- function(id, dir, off, title) {
      b <- actionButton(id, NULL, icon = icon(paste0("chevron-", dir)), class = "btn-sm",
        title = title, `aria-label` = title)
      if (off) shinyjs::disabled(b) else b
    }
    # A plain select reported only on a user change: a bound input would
    # re-send the old sample while the outgoing and incoming modals overlap.
    picker <- tags$select(class = "form-control input-sm", `aria-label` = "Open another sample",
      onchange = "Shiny.setInputValue('nav_pick', this.value, {priority: 'event'})",
      lapply(seq_len(nrow(o)), function(k) {
        tags$option(value = keys[k], selected = if (k == i) NA, labels[k])
      }))
    div(class = "mp-nav", title = sprintf("%d of %d in the table", i, nrow(o)),
      arrow("nav_prev", "left", i == 1, "Previous sample in the table"),
      picker,
      arrow("nav_next", "right", i == nrow(o), "Next sample in the table"))
  })
  nav_to <- function(key) {
    o <- nav_order(); req(o)
    hit <- o[nav_key(o) == key, , drop = FALSE]
    if (nrow(hit) != 1 || identical(key, nav_current())) return()
    if (edits_block()) return()
    removeModal()
    if ((input$mode %||% "Assemble") == "Annotate") {
      u <- list(ID = hit$ID, path = as.integer(hit$path), scaffold = as.integer(hit$scaffold))
      sel_id(u$ID); sel_unit(u)
      showModal(annotate_modal(u$ID, u$path, u$scaffold))
    } else {
      sel_id(hit$ID)
      showModal(assemble_modal(hit$ID))
    }
  }
  nav_step <- function(by) {
    o <- nav_order(); req(o)
    i <- match(nav_current(), nav_key(o)) + by
    if (!is.na(i) && i >= 1 && i <= nrow(o)) nav_to(nav_key(o[i, ]))
  }
  observeEvent(input$nav_prev, nav_step(-1L))
  observeEvent(input$nav_next, nav_step(1L))
  observeEvent(input$nav_pick, nav_to(input$nav_pick))

  # active mode's lock/state filter values (NULL for modes without pickers, i.e.
  # Taxonomy and Export)
  active_filters <- function() {
    mode <- input$mode %||% "Assemble"
    if (!mode %in% c("Assemble", "Annotate")) {
      return(NULL)
    }
    list(lock = input[[paste0(mode, "_lock_filter")]],
      state = input[[paste0(mode, "_state_filter")]])
  }

  # logical mask of samples_data() rows that pass the active mode's filters
  visible_mask <- reactive({
    d <- samples_data()
    keep <- rep(TRUE, nrow(d))
    f <- active_filters()
    if (is.null(f) || nrow(d) == 0) {
      return(keep)
    }
    mm <- mode_meta(input$mode %||% "Assemble")
    lc <- paste0(mm$pfx, "_lock")
    sc <- paste0(mm$pfx, "_switch")
    if (!is.null(f$lock) && lc %in% names(d)) keep <- keep & as.character(d[[lc]]) %in% f$lock
    if (!is.null(f$state) && sc %in% names(d)) keep <- keep & as.character(d[[sc]]) %in% f$state
    keep
  })

  # Bulk actions must only ever touch VISIBLE samples. Rows filtered out are
  # hidden with CSS but stay in reactable's row model, so intersect the selection
  # with the visible mask here (and prune it in the observer below).
  selected_rows <- function() {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel) || length(sel) == 0) return(integer(0))
    vis <- visible_mask()
    sel[sel >= 1 & sel <= length(vis) & vis[sel]]
  }
  selected_ids <- function() samples_data()$ID[selected_rows()]
  # Visible selected rows as (ID, path, scaffold) units (path/scaffold = 1 when the
  # active table has no unit columns).
  selected_units <- function() {
    sel <- selected_rows()
    if (!length(sel)) return(NULL)
    d <- samples_data()
    data.frame(ID = d$ID[sel],
      path = if ("path" %in% names(d)) as.integer(d$path[sel]) else 1L,
      scaffold = if ("scaffold" %in% names(d)) as.integer(d$scaffold[sel]) else 1L,
      stringsAsFactors = FALSE)
  }

  # Prune hidden rows out of the reactable selection whenever the selection or a
  # filter changes, so a select-all under an active filter cannot silently carry
  # invisible samples into a bulk Lock / State / Reviewed / Group action.
  observeEvent(
    list(reactable::getReactableState("tbl", "selected"),
      input$Assemble_lock_filter, input$Assemble_state_filter,
      input$Annotate_lock_filter, input$Annotate_state_filter, input$mode),
    {
      sel <- reactable::getReactableState("tbl", "selected")
      if (is.null(sel) || length(sel) == 0) {
        return()
      }
      vis <- visible_mask()
      keep <- sel[sel >= 1 & sel <= length(vis) & vis[sel]]
      if (length(keep) != length(sel)) reactable::updateReactable("tbl", selected = keep)
    }, ignoreInit = TRUE)

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

  # ---- Reference databases (suggest / build / reuse / skip) ----
  ref_h <- reference_handlers(input, output, session, list(
    con = app_con, project = app_project, refresh = refresh,
    selected_ids = selected_ids, table_data = samples_data))

  observeEvent(input$ref_move_btn, ref_h$open_move(selected_ids()))
  observeEvent(input$ref_keep_btn, ref_h$do_keep_together(selected_ids()))
  observeEvent(input$ref_release_btn, ref_h$do_dissolve(selected_ids()))

  # A build runs in the background, so the app has to ask the DB for its state.
  # Not gated on the Reference tab being open: a cluster build finishes on its
  # own schedule, and making the user sit on one tab to have it picked up is the
  # kind of thing nobody discovers until their database is silently unwired.
  observe({
    invalidateLater(4000, session)
    if (isTRUE(ref_h$poll() > 0)) refresh(isolate(refresh()) + 1)
  })

  # ---- bulk lock / state / reviewed actions ----
  # In Annotate mode each row is a (path, scaffold) unit, so the bulk update keys
  # on the unit; the annotate `sql` (ending "WHERE ID = ?") is rewritten to add
  # the unit predicate.
  set_flag <- function(sql, verb = "updated") {
    if (identical(input$mode %||% "Assemble", "Annotate")) return(set_unit_flag(sql, verb))
    ids <- selected_ids()
    if (!need_selection(length(ids))) return()
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    for (id in ids) DBI::dbExecute(con, sql, params = list(id))
    refresh(refresh() + 1)
    # refresh() re-renders the table, which drops the client-side selection.
    # Reselect the same samples by ID once the new table has been flushed.
    session$onFlushed(function() {
      idx <- which(isolate(samples_data())$ID %in% ids)
      if (length(idx) > 0) reactable::updateReactable("tbl", selected = idx)
    }, once = TRUE)
    mp_toast(paste(mp_n(length(ids), "sample"), paste0(verb, ".")))
  }
  set_unit_flag <- function(sql, verb = "updated") {
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    usql <- sub("WHERE ID = \\?", "WHERE ID = ? AND path = ? AND scaffold = ?", sql)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    for (i in seq_len(nrow(u))) DBI::dbExecute(con, usql, params = list(u$ID[i], u$path[i], u$scaffold[i]))
    refresh(refresh() + 1)
    session$onFlushed(function() {
      d <- isolate(samples_data())
      pp <- if ("path" %in% names(d)) d$path else rep(1L, nrow(d))
      ss <- if ("scaffold" %in% names(d)) d$scaffold else rep(1L, nrow(d))
      idx <- which(paste(d$ID, pp, ss) %in% paste(u$ID, u$path, u$scaffold))
      if (length(idx) > 0) reactable::updateReactable("tbl", selected = idx)
    }, once = TRUE)
    mp_toast(paste(mp_n(nrow(u), "unit"), paste0(verb, ".")))
  }
  # One direction per click, decided over the whole selection: any locked row
  # means the click unlocks everything selected, else it locks everything.
  # Unlocking drops rows from the next step's table, so it is confirmed first.
  lock_next <- reactiveVal(NULL)
  lock_sql <- function(mm, to) sprintf("UPDATE %s SET %s_lock = %d WHERE ID = ?", mm$tbl, mm$pfx, to)
  next_step <- c(taxonomy = "Assemble", assemble = "Annotation Review", annotate = "Export")
  observeEvent(list(input$tax_lock_btn, input$asmb_lock_btn, input$annot_lock_btn), ignoreInit = TRUE, {
    mm <- mode_meta(input$mode)
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    lc <- paste0(mm$pfx, "_lock")
    d <- samples_data()
    cur <- d[[lc]][selected_rows()]
    to <- if (identical(mp_flag_next(cur, on = 1L, off = 0L), 0L)) 0L else 1L
    n <- nrow(u)
    if (to == 0L) {
      lock_next(list(mm = mm, u = u))
      mp_confirm("unlock_confirm", title = paste("Unlock", mp_n(n, MP_UNIT[[mm$pfx]]), "?"),
        text = sprintf("Unlocked rows leave the %s table until locked again.", next_step[[mm$pfx]]),
        action_label = "Unlock", danger = TRUE)
      return()
    }
    # In Annotate mode locking is per unit, and only ONE assembly path may be
    # locked per sample (paths are alternative resolutions of the same repeat).
    if (identical(input$mode %||% "Assemble", "Annotate")) {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      ids <- unique(u$ID)
      ph <- paste(rep("?", length(ids)), collapse = ",")
      cur <- DBI::dbGetQuery(con, sprintf(
        "SELECT ID, path, scaffold, COALESCE(annotate_lock,0) AS lock FROM annotate WHERE ID IN (%s)", ph),
        params = as.list(ids))
      cur$newlock <- ifelse(paste(cur$ID, cur$path, cur$scaffold) %in% paste(u$ID, u$path, u$scaffold), 1L, cur$lock)
      viol <- tapply(seq_len(nrow(cur)), cur$ID, function(ix)
        length(unique(cur$path[ix][cur$newlock[ix] == 1L])) > 1)
      if (any(unlist(viol))) {
        # asked, not refused: the reviewer answers in the middle of the screen
        showModal(multi_path_lock_modal(names(viol)[unlist(viol)]))
        return()
      }
    }
    set_flag(lock_sql(mm, 1L), verb = "locked")
  })
  observeEvent(input$unlock_confirm, {
    x <- lock_next(); lock_next(NULL)
    req(isTRUE(input$unlock_confirm), x)
    set_flag(lock_sql(x$mm, 0L), verb = "unlocked")
  })
  observeEvent(input$lock_anyway, {
    removeModal()
    set_unit_flag("UPDATE annotate SET annotate_lock = 1 WHERE ID = ?", verb = "locked")
  })
  observeEvent(list(input$tax_state_btn, input$asmb_state_btn, input$annot_state_btn),
    ignoreInit = TRUE, {  # open a picker to choose the target state
    ids <- selected_ids()
    if (!need_selection(length(ids))) return()
    mm <- mode_meta(input$mode)
    # preselect only when every selected row agrees
    cur <- unique(samples_data()[[paste0(mm$pfx, "_switch")]][selected_rows()])
    showModal(mp_state_modal(ids,
      current = if (length(cur) == 1 && !is.na(cur)) as.character(cur) else character(0),
      mode = mm$pfx))
  })
  observeEvent(input$update_state, {
    if (is.null(input$new_state)) {
      mp_toast("Choose a state first.", type = "warning"); return()
    }
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state)
    set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st))
    removeModal()
  })
  # Toolbar flags: the whole selection goes one way (mp_flag_next), same as the
  # Mark / Clear toggles in the annotation window.
  toolbar_flag <- function(col) {
    d <- samples_data(); sel <- selected_rows()
    if (!need_selection(length(sel))) return()
    to <- mp_flag_next(d[[col]][sel])
    set_flag(sprintf("UPDATE annotate SET %s = '%s' WHERE ID = ?", col, to),
      verb = paste(if (to == "yes") "marked" else "cleared", col))
  }
  observeEvent(input$reviewed_btn, toolbar_flag("reviewed"))
  observeEvent(input$problematic_btn, toolbar_flag("problematic"))
  observeEvent(input$partial_btn, toolbar_flag("partial"))

  # ---- Export group actions ----
  observeEvent(input$assign_group_btn, {
    if (!need_selection(length(selected_ids()))) return()
    showModal(group_modal())
  })
  observeEvent(input$assign_group_go, {
    grp <- trimws(input$group_name %||% "")
    u <- selected_units()
    if (is.null(u) || nrow(u) == 0) { removeModal(); return() }
    if (!nzchar(grp)) { mp_toast("Type a group name first.", type = "warning"); return() }
    assign_export_group(app_project(), u, grp)
    removeModal(); refresh(refresh() + 1)
    mp_toast(paste0(mp_n(nrow(u), "unit"), " assigned to group '", grp, "'."))
  })
  observeEvent(input$clear_group_btn, {
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    clear_export_group(app_project(), u)
    refresh(refresh() + 1)
    mp_toast(paste0("Group cleared for ", mp_n(nrow(u), "unit"), "."))
  })

  # ---- lock / state row filters (hide via CSS on rp-lock-/rp-state- classes) ----
  # Only the active mode's filters apply, and the rules are scoped to the main
  # table, so a filter set in Assemble does not blank the Taxonomy/Export tables
  # (which share the #tbl output but have no filter controls). all_state includes
  # "0" so deselecting "Hold" actually hides held rows.
  output$filter_css <- renderUI({
    f <- active_filters()
    if (is.null(f)) {
      return(NULL)
    }
    hide <- character(0)
    if (!is.null(f$lock)) hide <- c(hide, paste0(".rp-lock-", setdiff(unname(LOCK_CHOICES), f$lock)))
    if (!is.null(f$state)) hide <- c(hide, paste0(".rp-state-", setdiff(unname(STATE_CHOICES), f$state)))
    if (length(hide) == 0) {
      return(NULL)
    }
    sel <- paste0("#tbl ", hide)
    tags$style(HTML(paste0(paste(sel, collapse = ", "), " { display:none !important; }")))
  })

  # ---- Update modal: run the workflow with live Nextflow progress ----
  nf_process <- reactiveVal(NULL)
  nf_state <- reactiveVal(nf_progress_empty())
  submit_base <- reactiveVal(NULL)   # timestamp shared by one submission's files

  output$nf_cmd_text <- renderText({
    args <- nf_command(app_project(), resume = isTRUE(input$nf_resume %||% TRUE))
    paste(c("nextflow", args), collapse = " ")
  })

  observeEvent(list(input$tax_run_btn, input$asmb_run_btn), ignoreInit = TRUE, {
    # Taxonomy runs locally in the app (lightweight taxdump lookup, no Nextflow).
    if ((input$mode %||% "") == "Taxonomy") { run_taxonomy_update(); return() }
    # block: a preset naming a database that is not there assembles nothing, and
    # the pipeline reports success while doing it (dangling symlink in, exit 0
    # out), so the run is pure wasted cluster time
    dead <- tryCatch({
      c0 <- app_con(); on.exit(DBI::dbDisconnect(c0), add = TRUE); missing_assembly_dbs(c0)
    }, error = function(e) character(0))
    if (length(dead) > 0) {
      rp_alert(sprintf(paste0(
        "Cannot run: %d reference database file(s) named by the samples queued for assembly ",
        "are missing or empty (%s). Rebuild them in the Reference tab, or point those samples ",
        "at a database that exists."), length(dead), paste(dead, collapse = ", ")), type = "error")
      return()
    }
    # warn, never block: a pending group only means those samples keep the
    # database their options set already names
    pending <- tryCatch({
      c1 <- app_con(); on.exit(DBI::dbDisconnect(c1), add = TRUE); reference_pending(c1)
    }, error = function(e) character(0))
    if (length(pending) > 0) {
      rp_alert(sprintf(
        "%d reference group(s) are not built or skipped (%s). Those samples will assemble with the database their options set already has.",
        length(pending), paste(pending, collapse = ", ")), type = "warning")
    }
    nf_state(nf_progress_empty())
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    # one timestamp per opening of the modal, so the script the user is looking
    # at is saved under the same name as the log paths written inside it
    if (hl) submit_base(job_basename("ribopilot"))
    script <- if (hl) build_submit_script(app_project(), resume = TRUE, sched = sched,
      base = submit_base()) else ""
    showModal(run_modal(nf_pending_count(app_project()),
      headless = hl, sched = sched, script = script))
  })

  # Headless run modal: save / submit the cluster job script (instead of running
  # Nextflow on the login node).
  write_submit_script <- function() {
    base <- submit_base() %||% job_basename("ribopilot")
    typed <- input$submit_script_text %||% ""
    script <- if (nzchar(trimws(typed))) typed else
      build_submit_script(app_project(), resume = isTRUE(input$nf_resume %||% TRUE),
        base = base)
    f <- job_path(app_project(), "script", paste0(base, ".sh"))
    writeLines(script, f)
    try(Sys.chmod(f, "0755"), silent = TRUE)
    f
  }
  observeEvent(input$save_script, {
    f <- write_submit_script()
    rp_alert(paste0("Saved job script to ", f), type = "message")
  })
  observeEvent(input$submit_job, {
    sched <- detect_scheduler()
    if (is.null(sched)) {
      rp_alert("No scheduler (sbatch/qsub/bsub) found on PATH.", type = "warning"); return()
    }
    # write the script HERE: it reads input$submit_script_text, and the later()
    # callback below has a session but no reactive context
    f <- write_submit_script()
    # a busy scheduler can sit on the submit for several seconds, and a dead
    # button reads as a missed click. The popup has no dismiss and covers the
    # modal, so the second click that would submit the job twice cannot land.
    shinyWidgets::sendSweetAlert(title = "Submitting...",
      text = paste0("Handing the job script to ", sched$cmd,
        ". Hold tight, a busy cluster can take a few seconds to answer."),
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        out <- tryCatch(system2(sched$cmd, shQuote(f), stdout = TRUE, stderr = TRUE),
          error = function(e) conditionMessage(e))
        rp_alert(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " ")),
          type = "message")
        removeModal()
      })
    }, delay = 0.1)
  })

  # Local (in-app) taxonomy auto-assignment for the ready, unlocked samples.
  # Shows a "hold tight" popup, then defers the (blocking) taxdump resolution so
  # the popup renders first, then swaps in a success/error popup (MitoPilot style).
  run_taxonomy_update <- function() {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ids <- DBI::dbGetQuery(con,
      "SELECT ID FROM samples WHERE taxonomy_lock = 0 AND taxonomy_switch IN (1,4)")$ID
    if (length(ids) == 0) {
      rp_alert(paste("No unlocked samples are queued for taxonomy assignment.",
        "Select samples and set their State to Ready, and make sure they are not locked."),
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

  progress_update <- function(lines) {
    nf_state(nf_progress_parse(nf_state(), lines))
  }

  observeEvent(input$nf_start, {
    nf_state(nf_progress_empty())
    shinyjs::hide("nf_start"); shinyjs::show("nf_stop"); shinyjs::show("nf_progress")
    shinyjs::removeClass(id = "rp_gears", class = "paused")
    # -ansi-log forces the redrawing progress board; wide COLUMNS keeps Nextflow
    # from truncating process names so the parser's per-process keys stay stable.
    nf_env <- c("current", COLUMNS = "500", NXF_ANSI_SUMMARY = "true")
    nxf_pin <- nf_pin_version()
    if (!is.na(nxf_pin)) nf_env <- c(nf_env, NXF_VER = nxf_pin)
    nf_record_blastdb_launch(app_project())
    p <- processx::process$new("nextflow",
      args = c(nf_command(app_project(), resume = isTRUE(input$nf_resume %||% TRUE)),
        "-ansi-log"),
      stdout = "|", stderr = "|", wd = app_project(),
      env = nf_env)
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
      rp_alert("Workflow run finished.", type = "message")
    }
  })

  output$nf_prog_all <- renderText(nf_progress_text(nf_state()))

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
    p <- nf_process()
    if (!is.null(p) && p$is_alive()) {
      mp_confirm("nf_close_confirm", title = "Stop the run?",
        text = "Closing this window stops the pipeline. Leave it open to keep running.",
        action_label = "Stop and close", danger = TRUE)
      return()
    }
    removeModal()
  })
  observeEvent(input$nf_close_confirm, {
    req(isTRUE(input$nf_close_confirm))
    kill_nf(); removeModal()
  })

  # ---- Export ----
  output$export_sel_n <- renderText({
    refresh()
    paste0(mp_n(nrow(export_units_df()), "unit"), " locked & ready")
  })
  # Exportable UNITS: annotation-locked, kept, not flagged problematic;
  # narrowed to one export group when given.
  export_units_df <- function(group = NULL) {
    locked <- app_query(
      "SELECT an.ID, an.path, an.scaffold, an.export_group
       FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       WHERE an.annotate_lock = 1 AND COALESCE(asm.ignore,0) = 0
         AND COALESCE(an.problematic, '') != 'yes'")
    if (is.null(group) || !nzchar(group)) return(locked)
    locked[!is.na(locked$export_group) & locked$export_group == group, , drop = FALSE]
  }
  export_group_choice <- function() {
    g <- input$export_group_sel %||% ""
    if (nzchar(g)) g else NULL
  }
  # SeqIDs are counted over ALL exportable units of a sample, not this bundle's
  # subset, so the header preview matches what export_bundle() will write.
  app_export_seqids <- function(con, u) {
    all_ids <- DBI::dbGetQuery(con,
      "SELECT an.ID FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       WHERE an.annotate_lock = 1 AND COALESCE(asm.ignore,0) = 0
         AND COALESCE(an.problematic,'') != 'yes'")$ID
    k <- as.integer(table(all_ids)[u$ID])
    ifelse(!is.na(k) & k > 1,
      sprintf("%s_p%s_s%s", u$ID, u$path, u$scaffold), u$ID)
  }
  # Token data the header templates are checked and resolved against.
  export_hdr_data <- reactive({
    refresh()
    u <- export_units_df(export_group_choice())
    if (nrow(u) == 0) return(NULL)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u$seqid <- app_export_seqids(con, u)
    tryCatch(export_header_data(con, u), error = function(e) NULL)
  })

  # Marker templates also resolve {element} / {product}, which belong to a
  # feature rather than to the record; a representative pair is enough to check
  # the template.
  export_hdr_data_marker <- reactive({
    d <- export_hdr_data()
    if (is.null(d)) return(NULL)
    el <- app_query("SELECT element FROM annotations LIMIT 1")$element
    el <- if (length(el)) el[1] else "18S"
    d$element <- el
    d$product <- rrna_element_product(el)
    d
  })

  # ---- header templates ----
  export_template <- reactiveVal("default")

  # "Export Data" opens a modal to pick a group before export. Built by a
  # function so the overwrite dialog's Cancel can bring it back.
  export_modal_current <- function() {
    groups <- app_query(
      "SELECT DISTINCT export_group FROM annotate
       WHERE export_group IS NOT NULL AND export_group != '' ORDER BY export_group")$export_group
    if (length(groups) == 0) {
      mp_alert(title = "No export group is assigned",
        text = paste("Export writes one group at a time. Select rows in the Export",
          "table, press Assign Group, then press Export Data."),
        type = "info")
      return(NULL)
    }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ensure_export_opts(con)
    tmpls <- list_export_templates(con)
    nm <- if (export_template() %in% tmpls) export_template() else "default"
    export_template(nm)
    o <- get_export_opts(con, nm)
    cols <- names(export_hdr_data_marker()) %||% character(0)
    sample_cols <- tryCatch(DBI::dbListFields(con, "samples"), error = function(e) character(0))
    export_modal(groups, tmpls, nm, o$fasta_header, o$fasta_header_marker, cols, sample_cols)
  }
  observeEvent(input$export_btn, {
    m <- export_modal_current()
    if (!is.null(m)) showModal(m)
  })

  # Picking a saved name refills the boxes; nothing is written until Save.
  observeEvent(input$export_template, ignoreInit = TRUE, {
    nm <- trimws(input$export_template %||% "")
    req(nzchar(nm))
    export_template(nm)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (nm %in% list_export_templates(con)) {
      o <- get_export_opts(con, nm)
      updateTextAreaInput(session, "fasta_header", value = o$fasta_header)
      updateTextAreaInput(session, "fasta_header_marker", value = o$fasta_header_marker)
    }
  })

  hdr_rec <- debounce(reactive(input$fasta_header), 400)
  hdr_mrk <- debounce(reactive(input$fasta_header_marker), 400)
  # Both boxes demand {seqid}: without it the per-marker files of different
  # records would collide, and the record deflines would disagree with the
  # feature table.
  hdr_rec_chk <- reactive(validate_fasta_header(hdr_rec(), export_hdr_data(), require_seqid = TRUE))
  hdr_mrk_chk <- reactive(validate_fasta_header(hdr_mrk(), export_hdr_data_marker(), require_seqid = TRUE))
  output$export_hdr_status <- renderUI(header_status_ui(hdr_rec_chk()))
  output$export_hdr_marker_status <- renderUI(header_status_ui(hdr_mrk_chk()))

  # The marker box is hidden unless markers are exported, so it must not gate
  # anything the user cannot see.
  markers_ok <- function() !isTRUE(input$export_markers) || isTRUE(hdr_mrk_chk()$ok)
  observe(shinyjs::toggleState("save_template",
    condition = isTRUE(hdr_rec_chk()$ok) && markers_ok()))

  observeEvent(input$save_template, {
    nm <- trimws(input$export_template %||% "")
    if (!nzchar(nm)) {
      mp_alert(title = "Name the template first",
        text = "Pick a template name, or type a new one, then press Save template.",
        type = "info")
      return()
    }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    set_export_opts(con, input$fasta_header, input$fasta_header_marker, nm)
    updateSelectizeInput(session, "export_template", choices = list_export_templates(con),
      selected = nm, options = list(create = TRUE, maxItems = 1))
    export_template(nm)
    mp_toast(sprintf("Saved header template \"%s\".", nm))
  })

  # What Export will write, for the group currently chosen.
  output$export_summary <- renderUI({
    grp <- export_group_choice(); req(grp)
    refresh()
    n <- nrow(export_units_df(grp))
    path <- file.path(app_project(), "out", "export", grp)
    files <- c(paste0(grp, ".fasta"), paste0(grp, ".tbl"), "samples/",
      paste0(grp, "_summary.csv"), "README.txt")
    if (isTRUE(input$export_markers)) files <- c(files, "markers/")
    if (isTRUE(input$export_report)) files <- c(files, paste0(grp, "_report.html"))
    div(style = paste("font-size:var(--mp-fs-meta); padding:8px 12px; margin-top:12px;",
        "background:var(--mp-surface-alt); border-left:3px solid var(--mp-primary);"),
      div(sprintf("%s in group \"%s\".", mp_n(n, "unit"), grp)),
      div("Written to ", tags$code(class = "mp-path", path), " as: ", paste(files, collapse = ", ")),
      if (dir.exists(path)) div(class = "mp-fg-warning", icon("triangle-exclamation"), " ",
        "This folder already exists. Export deletes it and writes it again."))
  })
  # Write the bundle into the project's out/export folder (no zip download,
  # mirroring MitoPilot), then show where it landed with an "Open folder" button.
  export_done_dir <- reactiveVal(NULL)
  observeEvent(input$export_go, {
    grp <- export_group_choice()
    u <- if (is.null(grp)) NULL else export_units_df(grp)
    if (is.null(u) || nrow(u) == 0) {
      rp_alert("No exportable units in this group. Lock annotations first.",
        type = "warning")
      return()
    }
    rc <- validate_fasta_header(input$fasta_header, export_hdr_data(), require_seqid = TRUE)
    mc <- if (isTRUE(input$export_markers)) {
      validate_fasta_header(input$fasta_header_marker, export_hdr_data_marker(), require_seqid = TRUE)
    } else list(level = "ok")
    bad <- if (identical(rc$level, "error")) list("Record FASTA header", rc$message)
      else if (identical(mc$level, "error")) list("Per-marker FASTA header", mc$message)
    if (!is.null(bad)) {
      shinyWidgets::sendSweetAlert(session, title = "Header cannot be used",
        text = paste0(bad[[1]], ": ", bad[[2]]), type = "error")
      return()
    }
    # export_bundle() rebuilds the group folder from scratch, so anything the
    # curator put there (or an earlier export they have not submitted yet) goes
    # with it. Say so before it happens, not after.
    target <- file.path(app_project(), "out", "export", grp)
    if (dir.exists(target)) {
      showModal(export_overwrite_modal(target,
        length(list.files(target, recursive = TRUE))))
      return()
    }
    run_export(u, grp)
  })
  observeEvent(input$export_overwrite_back, {
    m <- export_modal_current()
    if (!is.null(m)) showModal(m)
  })
  observeEvent(input$export_overwrite_go, {
    grp <- export_group_choice()
    run_export(export_units_df(grp), grp)
  })
  run_export <- function(u, grp) {
    detail <- if (isTRUE(input$export_itsx_check)) {
      "Writing files and running the final ITSx check. Hold tight, this can take a few minutes."
    } else {
      "Writing files. Hold tight."
    }
    gdir <- withProgress(message = "Exporting...", detail = detail, value = 0.5, tryCatch(
      export_bundle(app_project(), units = u, group = grp,
        out_dir = file.path(app_project(), "out"), zip = FALSE,
        fasta_header = input$fasta_header,
        fasta_header_marker = input$fasta_header_marker,
        markers = isTRUE(input$export_markers),
        report = isTRUE(input$export_report),
        itsx_check = isTRUE(input$export_itsx_check)),
      error = function(e) {
        rp_alert(paste("Export failed:", conditionMessage(e)), type = "error")
        NULL
      }))
    if (is.null(gdir)) return(invisible())
    export_done_dir(gdir)
    removeModal()
    showModal(modalDialog(
      title = "Export complete", size = "m", easyClose = TRUE,
      tags$p("Data exported to:"),
      div(class = "code-block", style = "white-space:normal; word-break:break-all;", gdir),
      export_check_summary_ui(attr(gdir, "summary")),
      footer = mp_footer(dismiss = "Close",
        extra = actionButton("open_export_dir", "Open folder", icon = icon("folder-open"), class = "btn-sm"))
    ))
    invisible()
  }
  observeEvent(input$open_export_dir, {
    d <- export_done_dir()
    # a remote session with no RStudio Files pane can show the user nothing;
    # saying where the export is beats a button that appears to do nothing
    if (!open_path(d)) {
      rp_alert(paste0("This session has no file browser to open. The export is at: ", d),
        type = "warning")
    }
  })

  # ---- Assemble modal ----
  # Single source of truth for the scaffold table AND its click handlers. The
  # handlers used to re-query and index by row position, so any change to the
  # display order silently acted on the wrong scaffold; clicks are now keyed on
  # the row's own uid and nothing depends on ordering.
  scaffold_rows <- reactive({
    id <- sel_id(); req(id); refresh()
    unit_evidence(id)
  })
  unit_by_uid <- function(uid) {
    d <- scaffold_rows()
    if (is.null(uid) || !nrow(d)) return(NULL)
    p <- strsplit(as.character(uid), "|", fixed = TRUE)[[1]]
    if (length(p) < 3) return(NULL)
    r <- d[d$path == as.integer(p[2]) & d$scaffold == as.integer(p[3]), , drop = FALSE]
    if (!nrow(r)) NULL else r[1, ]
  }
  output$assemble_scaffolds <- reactable::renderReactable({
    d <- scaffold_rows()
    if (nrow(d) == 0) {
      return(reactable::reactable(
        data.frame(` ` = "No assembled scaffolds for this sample yet.", check.names = FALSE),
        compact = TRUE))
    }
    # reactable renders in data-frame column order, so the row identifier has to
    # lead here; colDef(sticky=) alone does not move it.
    lead <- c("unit", "source_contig", "triage", "length", "path_score", "depth", "depth_track", "elements",
      "screen_verdict", "screen_reason", "variant_class")
    d <- d[, c(intersect(lead, names(d)), setdiff(names(d), lead)), drop = FALSE]
    reactable::reactable(d, compact = TRUE, resizable = TRUE, filterable = TRUE,
      searchable = TRUE, sortable = TRUE, wrap = FALSE, defaultPageSize = 20,
      # tick boxes drive Copy FASTA and, with more than one ticked, a bulk K/D
      # verdict; clicking the row anywhere else opens that scaffold below
      selection = "multiple",
      onClick = rt_examine_row("scaf_examine"),
      meta = list(examined = isolate(examined_uid())),
      defaultColDef = reactable::colDef(align = "left", show = FALSE),
      columns = list(
        # the tick-box column is a built-in column, so defaultColDef's show =
        # FALSE hides it too; it has to opt back in by name. The copy icon says
        # what the ticks are for, since they do NOT drive the panels below.
        .selection = reactable::colDef(show = TRUE, width = 36, sticky = "left",
          html = TRUE, header = as.character(icon("clipboard",
            title = paste("Tick rows to copy them as FASTA, or to keep/drop them",
              "together: with more than one ticked, K or D on any ticked row",
              "applies to all of them. K or D on an unticked row affects only",
              "that row.")))),
        # Every visible column explains itself on hover ("?" in the header), so a
        # reviewer never has to remember which number is comparable with what.
        unit = reactable::colDef(name = "Unit", width = 72, show = TRUE, sticky = "left",
          html = TRUE, header = rt_help_header("Unit", SCAFFOLD_COL_HELP$unit)),
        length = reactable::colDef(name = "Length", width = 80, show = TRUE,
          html = TRUE, header = rt_help_header("Length", SCAFFOLD_COL_HELP$length)),
        # only meaningful for a supplied assembly: for an assembled sample this
        # holds the assembler's own scaffold header, which says nothing a
        # reviewer can act on
        source_contig = reactable::colDef(name = "Your contig", width = 150,
          show = sample_has_user_assembly(sel_id()), html = TRUE,
          header = rt_help_header("Your contig", SCAFFOLD_COL_HELP$source_contig),
          cell = rt_longtext()),
        path_score = reactable::colDef(name = "Path score", width = 115, show = TRUE,
          html = TRUE, header = rt_help_header("Path score", SCAFFOLD_COL_HELP$path_score)),
        # The number is this scaffold's absolute median depth; the bar beside it
        # is that value as a fraction of the sample's deepest scaffold. Naming
        # the column "relative" described the bar and misdescribed the number.
        depth = reactable::colDef(name = "Median depth", width = 125, show = TRUE, html = TRUE,
          header = rt_help_header("Median depth", SCAFFOLD_COL_HELP$depth),
          cell = rt_depth_bar("depth_max")),
        depth_track = reactable::colDef(name = "Coverage", width = 100, show = TRUE,
          filterable = FALSE, sortable = FALSE, html = TRUE,
          header = rt_help_header("Coverage", SCAFFOLD_COL_HELP$depth_track),
          cell = rt_depth_sparkline()),
        elements = reactable::colDef(name = "rRNA features", width = 240, show = TRUE, html = TRUE,
          header = rt_help_header("rRNA features", SCAFFOLD_COL_HELP$elements),
          cell = rt_element_strip()),
        screen_verdict = reactable::colDef(name = "Screen", width = 100, show = TRUE, html = TRUE,
          header = rt_help_header("Screen", SCAFFOLD_COL_HELP$screen_verdict),
          cell = rt_screen_badge()),
        # advisory only, not a primary column: tells a curator why a scaffold has
        # no reference to judge against, including when a local panel search
        # came back empty and a remote search might still resolve it
        screen_reason = reactable::colDef(name = "Screen why", width = 140, show = TRUE, html = TRUE,
          header = rt_help_header("Screen why", SCAFFOLD_COL_HELP$screen_reason),
          cell = rt_longtext()),
        # Information, not a flag: neither an isomer nor a paralog is a defect.
        variant_class = reactable::colDef(name = "Variant", width = 115, show = TRUE, html = TRUE,
          header = rt_help_header("Variant", SCAFFOLD_COL_HELP$variant_class),
          cell = rt_variant_badge()),
        chimera_label = reactable::colDef(name = "Chimera", width = 90, show = TRUE, html = TRUE,
          header = rt_help_header("Chimera", SCAFFOLD_COL_HELP$chimera), cell = rt_chimera_badge()),
        chimera_reason = reactable::colDef(show = FALSE),
        itsx_warnings = reactable::colDef(name = "ITSx Warnings", minWidth = 130, show = TRUE, html = TRUE,
          header = rt_help_header("ITSx Warnings", SCAFFOLD_COL_HELP$itsx_warnings), cell = rt_itsx_badge()),
        screen_species = reactable::colDef(name = "Top GenBank hit", minWidth = 160, show = TRUE,
          html = TRUE, header = rt_help_header("Top GenBank hit", SCAFFOLD_COL_HELP$screen_species),
          cell = rt_longtext()),
        shares = reactable::colDef(name = "Shares", width = 105, show = TRUE, html = TRUE,
          header = rt_help_header("Shares", SCAFFOLD_COL_HELP$shares),
          cell = rt_shares_badge()),
        odd_one_out = reactable::colDef(name = "", width = 55, show = TRUE, html = TRUE,
          filterable = FALSE, header = rt_help_header("", SCAFFOLD_COL_HELP$odd_one_out),
          cell = rt_odd_badge()),
        screen_accession = reactable::colDef(name = "Accession", width = 110, show = TRUE,
          html = TRUE, header = rt_help_header("Accession", SCAFFOLD_COL_HELP$screen_accession),
          cell = rt_ncbi_link()),
        screen_pident = reactable::colDef(name = "%ID", width = 65, show = TRUE,
          html = TRUE, header = rt_help_header("%ID", SCAFFOLD_COL_HELP$screen_pident)),
        triage = reactable::colDef(name = "Keep", width = 68, show = TRUE, html = TRUE,
          filterable = FALSE, header = rt_help_header("Keep", SCAFFOLD_COL_HELP$triage),
          cell = rt_triage_bttn("scaf_triage")),
        triage_reason = reactable::colDef(name = "Why", minWidth = 200, show = TRUE,
          html = TRUE, header = rt_help_header("Why", SCAFFOLD_COL_HELP$triage_reason),
          cell = rt_longtext()),
        uid = reactable::colDef(show = FALSE), depth_max = reactable::colDef(show = FALSE),
        sib_taxon = reactable::colDef(show = FALSE), sib_rank = reactable::colDef(show = FALSE)),
      # JS rather than an R function of the row index: the examined row changes
      # on every click, and re-rendering the table to restyle it would throw away
      # the user's sort, filters and page.
      rowStyle = rt_scaffold_row_style())
  })

  # The ticked scaffolds as (path, scaffold) rows. Two controls read this: Copy
  # FASTA, and a K/D click while more than one row is ticked.
  ticked_units <- function() {
    d <- scaffold_rows()
    sel <- reactable::getReactableState("assemble_scaffolds", "selected")
    sel <- sel[!is.na(sel) & sel >= 1 & sel <= nrow(d)]
    d[sel, c("path", "scaffold"), drop = FALSE]
  }

  # ---- Copy FASTA for the ticked scaffolds ----
  observeEvent(input$copy_fasta, {
    u <- ticked_units()
    if (nrow(u) == 0) {
      rp_alert("Tick one or more scaffolds first.", type = "warning"); return()
    }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    txt <- unit_fasta_text(con, sel_id(), u)
    if (length(txt) == 0) {
      rp_alert("Those scaffolds have no stored sequence.", type = "warning"); return()
    }
    copy_to_clipboard(paste(txt, collapse = "\n"))
    rp_toast(sprintf("Copied %d scaffold(s) as FASTA.", nrow(u)), type = "message")
  })

  # ---- per-unit BLAST evidence: the unfiltered screen vs the rDNA panel ----
  # The examined unit drives both panels below the table. Clicking a row picks
  # it; the tick boxes are for Copy FASTA and bulk K/D, and do not move it.
  sel_hit_unit <- reactiveVal(NULL)
  examined_uid <- reactive({
    h <- sel_hit_unit()
    if (is.null(h)) NA_character_ else paste(sel_id(), h$path, h$scaffold, sep = "|")
  })
  # push the highlight without re-rendering, which would drop sort and filters
  examine_unit <- function(path, scaffold) {
    sel_hit_unit(list(path = path, scaffold = scaffold))
    shinyjs::runjs("document.querySelectorAll('#asm_unit_panels details').forEach(function(d){ d.open = false; });")
    reactable::updateReactable("assemble_scaffolds",
      meta = list(examined = paste(sel_id(), path, scaffold, sep = "|")))
  }
  observeEvent(input$scaf_examine, {
    u <- unit_by_uid(input$scaf_examine); req(!is.null(u))
    examine_unit(u$path, u$scaffold)
  })
  # default to the first kept unit so the panel is never empty on open
  observeEvent(scaffold_rows(), {
    d <- scaffold_rows()
    if (nrow(d) == 0) { sel_hit_unit(NULL); return() }
    cur <- sel_hit_unit()
    if (!is.null(cur) && any(d$path == cur$path & d$scaffold == cur$scaffold)) return()
    k <- which(d$ignore == 0)[1]
    if (is.na(k)) k <- 1L
    examine_unit(d$path[k], d$scaffold[k])
  })
  output$unit_screen <- renderUI({
    h <- sel_hit_unit(); d <- scaffold_rows()
    if (is.null(h) || nrow(d) == 0) return(NULL)
    r <- d[d$path == h$path & d$scaffold == h$scaffold, , drop = FALSE]
    if (nrow(r) == 0) return(NULL)
    r <- r[1, ]
    lbl <- sprintf("Scaffold p%s.s%s", r$path, r$scaffold)
    if (is.na(r$screen_accession)) {
      # Two different states used to share one message: BLAST switched off, and
      # BLAST switched on but the remote search never came back. The second told
      # the user to enable something already enabled.
      return(tagList(
        tags$b(lbl),
        div(style = "color:#888; margin-top:4px;",
          if (blast_enabled(sel_id())) {
            paste("No GenBank screen for this scaffold. The remote BLAST is still",
              "running, or it did not return (NCBI queue, network, or retries used up).")
          } else {
            paste("GenBank BLAST is switched off for this sample (Assemble tab,",
              "BLAST options). Run the pipeline (Assemble > Update) with BLAST enabled to screen this scaffold.")
          }),
        seed_hit_block(r)))
    }
    tagList(
      tags$b(lbl),
      div(style = "margin-top:4px;",
        tags$span(style = "color:#555;", "Best match anywhere in core_nt, NOT restricted to nuclear rDNA:")),
      div(class = "code-block", style = "white-space:normal; margin-top:3px;",
        HTML(sprintf("<b>%s</b> &nbsp; %s &nbsp; %s%% id &nbsp; %s%% cov<br><span style='color:#666;'>%s</span>",
          r$screen_species %|NA|% "unknown",
          sprintf("<a href='https://www.ncbi.nlm.nih.gov/nuccore/%s' target='_blank'>%s</a>",
            r$screen_accession, r$screen_accession),
          r$screen_pident %|NA|% "?", r$screen_qcovs %|NA|% "?",
          r$screen_title %|NA|% ""))),
      div(style = "margin-top:5px;",
        # yes is the GOOD case here (the hit is a nuclear rDNA record), so no invert.
        # NA means no hit said either way and must not render as "no".
        status_badge("Nuclear rDNA",
          if (is.na(r$is_nuclear_rrna)) NA else if (r$is_nuclear_rrna == 1) "yes" else "no"),
        tags$span(style = "margin-left:8px; color:#555;",
          sprintf("shares %s with %s",
            if (is.na(r$shares)) "no resolved rank" else r$shares,
            r$sci_name %|NA|% "this sample"))),
      if (isTRUE(r$odd_one_out)) div(class = "mp-fg-warning", style = "margin-top:5px;",
        sprintf("This scaffold's top hit differs from the rest of this sample (%s at %s).",
          r$sib_taxon %|NA|% "the majority", r$sib_rank %|NA|% "")))
  })
  output$unit_candidates <- reactable::renderReactable({
    refresh()
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    d <- app_query(
      "SELECT rank, accession, species, title, pident, qcovs, bitscore, evalue
       FROM blast_ref_candidates WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY rank",
      params = list(id, h$path, h$scaffold))
    if (nrow(d) == 0) {
      # An empty panel is itself a finding: the nuclear-rDNA title filter matched
      # nothing, which is common in sparse invertebrate clades. But when the
      # screen above is empty too, pointing at it sent the reader in a circle.
      sr <- scaffold_rows()
      sr <- sr[sr$path == h$path & sr$scaffold == h$scaffold, , drop = FALSE]
      screened <- nrow(sr) > 0 && !is.na(sr$screen_accession[1])
      return(reactable::reactable(data.frame(
        ` ` = if (screened) {
          "No nuclear rDNA reference matched this scaffold. See the unrestricted screen above."
        } else {
          "No GenBank result for this scaffold, so there is nothing to rank here yet."
        },
        check.names = FALSE), compact = TRUE))
    }
    reactable::reactable(d, compact = TRUE, resizable = TRUE, filterable = TRUE,
      wrap = FALSE, defaultPageSize = 10, defaultColDef = reactable::colDef(align = "left"),
      columns = list(
        rank = reactable::colDef(name = "#", width = 45),
        accession = reactable::colDef(name = "Accession", width = 115, html = TRUE, cell = rt_ncbi_link()),
        species = reactable::colDef(name = "Species", minWidth = 140, html = TRUE, cell = rt_longtext()),
        title = reactable::colDef(name = "Title", minWidth = 240, html = TRUE, cell = rt_longtext()),
        pident = reactable::colDef(name = "%ID", width = 65),
        qcovs = reactable::colDef(name = "Cov", width = 60),
        bitscore = reactable::colDef(name = "Bits", width = 70),
        evalue = reactable::colDef(name = "E", width = 80)))
  })
  # Panel titles name the highlighted scaffold; one output per title since an
  # output id can render in only one place.
  for (o in c("hl_scaf_ev", "hl_scaf_cov", "hl_scaf_chim")) local({
    oid <- o
    output[[oid]] <- renderText({
      h <- sel_hit_unit(); req(!is.null(h)); sprintf("p%s.s%s", h$path, h$scaffold)
    })
  })
  output$unit_coverage_plot <- renderPlot({
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    grid::grid.draw(fixed_gutter(plot_coverage(id, h$path, h$scaffold)))
  }, height = 130)
  # Triage one scaffold: keep / drop / undecided. `assemblies.ignore` is derived
  # from the verdict by set_unit_triage(), so nothing downstream changes.
  # Dropping is no longer destructive (the annotate row is held, not deleted),
  # so a reviewed unit only warns; an annotation-locked one is still refused,
  # because the Annotate-mode lock guard cannot see this path.
  pending_triage <- reactiveVal(NULL)
  unit_label <- function(u) paste0(u$path, ".", u$scaffold)
  do_scaf_triage <- function(id, units, verdict) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    for (i in seq_len(nrow(units))) {
      set_unit_triage(con, id, units$path[i], units$scaffold[i], verdict,
        reason = NA_character_, source = "user")
    }
    refresh(refresh() + 1)
    word <- if (is.na(verdict)) "undecided" else verdict
    rp_toast(if (nrow(units) == 1) {
      sprintf("Scaffold %s marked %s.", unit_label(units[1, ]), word)
    } else {
      sprintf("%d scaffolds marked %s.", nrow(units), word)
    }, type = "message")
  }
  # Split the units a drop would touch into the ones it must refuse, the ones it
  # must ask about, and the rest. Annotation-locked is a refusal because the
  # Annotate-mode lock guard cannot see this path; reviewed or exported only
  # warns, since the annotate row is held rather than deleted.
  classify_drop <- function(id, units) {
    locked <- character(0); ask <- character(0); keep_idx <- integer(0)
    for (i in seq_len(nrow(units))) {
      p <- units$path[i]; s <- units$scaffold[i]
      st <- app_query(
        "SELECT annotate_lock, reviewed, export_group FROM annotate
         WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, p, s))
      if (nrow(st) && isTRUE(as.integer(st$annotate_lock[1]) == 1L)) {
        locked <- c(locked, paste0(p, ".", s)); next
      }
      if (nrow(st) && (identical(st$reviewed[1], "yes") ||
                       (!is.na(st$export_group[1]) && nzchar(st$export_group[1])))) {
        ask <- c(ask, paste0(p, ".", s))
      }
      keep_idx <- c(keep_idx, i)
    }
    list(locked = locked, ask = ask, units = units[keep_idx, , drop = FALSE])
  }
  observeEvent(input$scaf_triage, {
    id <- sel_id(); req(id)
    parts <- strsplit(as.character(input$scaf_triage), "|", fixed = TRUE)[[1]]
    req(length(parts) >= 4)
    # exact match only: match.arg() would silently expand an abbreviation, so a
    # malformed payload would act rather than fail
    want <- parts[4]
    req(want %in% TRIAGE_VERDICTS)
    u <- unit_by_uid(paste(parts[1:3], collapse = "|")); req(!is.null(u))

    # A click INSIDE a multi-row selection is a bulk verdict over that selection.
    # A click OUTSIDE it stays a single-row action on the row clicked: the ticked
    # rows are not what the user is pointing at, and retargeting the click onto
    # them would change rows they never aimed at while leaving the row they did
    # aim at alone.
    # Bulk SETS rather than toggles: a mixed selection has no single state to
    # toggle out of, so "clear it again" would mean different things per row.
    ticked <- ticked_units()
    on_ticked <- nrow(ticked) > 0 &&
      any(ticked$path == u$path & ticked$scaffold == u$scaffold)
    if (nrow(ticked) > 1 && on_ticked) {
      units <- ticked
      verdict <- want
    } else {
      units <- u[, c("path", "scaffold"), drop = FALSE]
      # clicking the state a row is already in clears it back to undecided
      verdict <- if (identical(u$triage %|NA|% "", want)) NA_character_ else want
    }

    if (identical(verdict, "drop")) {
      cls <- classify_drop(id, units)
      if (length(cls$locked)) {
        rp_alert(sprintf(
          "Scaffold%s %s %s annotation-locked. Unlock %s in Annotation Review first.",
          if (length(cls$locked) == 1) "" else "s", paste(cls$locked, collapse = ", "),
          if (length(cls$locked) == 1) "is" else "are",
          if (length(cls$locked) == 1) "that unit" else "those units"),
          type = if (nrow(cls$units)) "warning" else "error")
      }
      if (nrow(cls$units) == 0) return()
      if (length(cls$ask)) {
        pending_triage(list(id = id, units = cls$units))
        showModal(modalDialog(
          title = if (length(cls$ask) == 1) "Drop a reviewed scaffold?" else "Drop reviewed scaffolds?",
          size = "s",
          tags$p(sprintf(paste("Scaffold%s %s already been reviewed. Dropping",
            "%s from annotation and export. Annotation and curation are kept, so this can be undone."),
            if (length(cls$ask) == 1) paste0(" ", cls$ask) else
              paste0("s ", paste(cls$ask, collapse = ", ")),
            if (length(cls$ask) == 1) "has" else "have",
            if (nrow(cls$units) == 1) "excludes it" else
              sprintf("excludes %d scaffolds", nrow(cls$units)))),
          footer = mp_footer(dismiss = NULL,
            extra = actionButton("scaf_triage_back", "Cancel"),
            danger = actionButton("scaf_triage_ok", "Drop it"))))
        return()
      }
      units <- cls$units
    }
    do_scaf_triage(id, units, verdict)
  })
  observeEvent(input$scaf_triage_back, {
    pt <- pending_triage(); pending_triage(NULL); req(pt)
    removeModal(); showModal(assemble_modal(pt$id))
  })
  observeEvent(input$scaf_triage_ok, {
    pt <- pending_triage(); req(pt); pending_triage(NULL)
    removeModal()
    do_scaf_triage(pt$id, pt$units, "drop")
    showModal(assemble_modal(pt$id))
  })
  # ---- bulk triage over the selected samples ----
  triage_target_ids <- function() {
    ids <- unique(selected_ids())
    if (length(ids)) ids else unique(samples_data()$ID)  # no selection: the whole table
  }
  # the preset the modal was opened with; seeding from the shipped defaults
  # instead would report -- and on Apply write back -- settings the curator
  # never chose, until the browser echoed the controls back
  triage_preset <- reactiveVal(NULL)
  # settings come from the DB so the modal and the WF1 triage process agree; the
  # inputs override them live so the preview reflects what Apply would do.
  triage_actions <- reactive({
    a <- triage_preset() %||% triage_actions_default()
    for (k in triage_rule_keys()) {
      v <- input[[paste0("triage_rule_", k)]]
      if (!is.null(v) && nzchar(v)) a[[k]] <- v
    }
    normalize_triage_actions(a)
  })
  observeEvent(input$triage_btn, {
    ids <- triage_target_ids()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    cur <- if (length(ids)) triage_opts_for(con, ids[1]) else triage_actions_default()
    triage_preset(cur)
    showModal(triage_modal(length(ids), cur))
  })
  triage_pv <- reactive({
    ids <- triage_target_ids()
    if (!length(ids)) return(NULL)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    tryCatch(triage_preview(con, ids, opts = triage_actions(),
      overwrite_user = isTRUE(input$triage_overwrite_user)), error = function(e) NULL)
  })
  # one count output per rule, rendered into the static table so changing a
  # control never rebuilds (and resets) the controls themselves
  for (k in triage_rule_keys()) {
    local({
      key <- k
      output[[paste0("triage_n_", key)]] <- renderUI({
        p <- triage_pv()
        if (is.null(p) || is.null(p$by_rule)) return(span("-"))
        r <- p$by_rule[p$by_rule$key == key, ]
        if (nrow(r) == 0) return(span("-"))
        span(style = "white-space:nowrap;",
          tags$strong(r$n_decided[1]), " / ",
          span(style = "color:#777;", r$n_matching[1]))
      })
    })
  }
  output$triage_preview <- renderUI({
    p <- triage_pv()
    if (is.null(p)) return(em("No samples to triage."))
    no_drop <- !any(triage_actions() %in% "drop")
    all_off <- triage_all_off(triage_actions())
    tagList(
      tags$p(
        if (all_off)
          "Every rule is off, so Apply would change nothing."
        else if (p$n_changed == 0 && p$n_units > 0 &&
            p$n_user_skipped + p$n_locked_skipped == p$n_units)
          # nothing is compared against these rules here, so saying the stored
          # verdicts match them would be a claim about units never examined
          "Every scaffold in this selection is hand-decided or annotation-locked, so Apply would skip them all."
        else if (p$n_changed == 0)
          "Apply would change nothing: every scaffold already has the verdict these rules propose."
        else
          sprintf(paste("Apply would change %d of %d scaffold(s):",
            "%d newly dropped, %d newly kept, %d newly flagged."),
            p$n_changed, p$n_units, p$n_new_drop, p$n_new_keep, p$n_new_flag)),
      tags$p(sprintf("Resulting state across %d sample(s): %d dropped, %d kept, %d undecided.",
        p$n_samples, p$n_drop, p$n_keep, p$n_undecided)),
      # the all-off message above already says Apply is a no-op; do not also
      # print the no-Drop note, or the modal shows two overlapping sentences
      if (all_off) NULL
      else if (no_drop) tags$p(style = "color:#555;",
        "No rule is set to Drop, so Apply will only add and clear reasons.")
      else if (p$n_new_drop == 0) tags$p(style = "color:#555;",
        "A rule is set to Drop, but nothing in this selection will be dropped."),
      if (p$n_user_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d scaffold(s) you decided by hand will be left alone.", p$n_user_skipped)),
      if (p$n_locked_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d annotation-locked scaffold(s) will be left alone.", p$n_locked_skipped)),
      # under all-off nothing is caused by Apply, so this reports the standing
      # state rather than an outcome
      if (length(p$emptied) > 0) div(class = "mp-fg-danger", style = "font-weight:600; margin-top:6px;",
        sprintf(if (all_off) "%d sample(s) currently have NO scaffolds kept: %s"
          else "%d sample(s) would end up with NO scaffolds kept: %s",
          length(p$emptied), paste(utils::head(p$emptied, 12), collapse = ", ")))
      else div(class = "mp-fg-success", style = "margin-top:6px;",
        if (all_off) "No sample currently has nothing kept."
        else "No sample would be left with nothing kept."))
  })
  observeEvent(input$triage_apply, {
    ids <- triage_target_ids(); req(length(ids))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    act <- triage_actions()
    # persist first: a pipeline re-run must reproduce what the curator just saw
    write_triage_opts(con, act)
    n <- 0L
    for (id in ids) {
      n <- n + tryCatch(auto_triage_sample(con, id, opts = act,
        overwrite_user = isTRUE(input$triage_overwrite_user)), error = function(e) 0L)
    }
    removeModal()
    refresh(refresh() + 1)
    rp_alert(sprintf("Triage updated %d scaffold(s) across %d sample(s).", n, length(ids)),
      type = "message")
  })

  # ---- Annotate modal: gene list + boundary editor (edit/save/reset) + msaR ----
  edit_tick <- reactiveVal(0)         # bumps to re-render map/seq/table
  active_feature <- reactiveVal(NULL) # feature the editor + alignment act on
  edit_active <- reactiveVal(FALSE)   # in edit mode?
  edit_backup <- reactiveVal(NULL)    # annotations snapshot for Reset
  edit_dirty <- reactiveVal(FALSE)    # nudged since entering edit mode?

  # reset per-unit modal state; select the first feature + load notes
  observeEvent(sel_unit(), ignoreNULL = FALSE, {
    edit_active(FALSE); edit_backup(NULL); edit_dirty(FALSE)
    shinyjs::hide("feat_edit_ctrls")
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
    u <- sel_unit(); active_feature(NULL)
    if (!is.null(u)) {
      els <- app_query("SELECT element FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
        params = uargs())$element
      if (length(els)) active_feature(els[1])
      a <- app_query("SELECT annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
      updateTextAreaInput(session, "annotate_notes", value = a$annotate_notes[1] %||% "")
    }
  })

  output$feat_edit_active <- renderText({
    el <- active_feature(); if (is.null(el)) "No active feature" else paste("Active feature:", el)
  })

  output$coverage_plot <- renderPlot({
    edit_tick(); u <- sel_unit(); req(u); plot_coverage(u$ID, u$path, u$scaffold)
  }, height = 128)

  # Chimera panels follow the annotate modal's unit, or the assemble modal's
  # highlighted scaffold.
  hit_unit <- reactive({
    if ((input$mode %||% "Assemble") == "Annotate") {
      u <- sel_unit(); req(u); return(u)
    }
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    list(ID = id, path = h$path, scaffold = h$scaffold)
  })
  output$chimera_plot <- renderPlot({
    u <- hit_unit(); grid::grid.draw(fixed_gutter(plot_chimera(u$ID, u$path, u$scaffold)))
  }, height = 296)
  output$chimera_db_source <- renderText({
    refresh(); u <- hit_unit()
    ch <- app_query("SELECT db_source FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$db_source[1]))
    paste("windows searched against:", paste(basename(strsplit(ch$db_source[1], " ", fixed = TRUE)[[1]]), collapse = ", "))
  })
  output$chimera_segments <- reactable::renderReactable({
    refresh()
    u <- hit_unit()
    ch <- app_query("SELECT segments FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$segments[1]))
    s <- jsonlite::fromJSON(ch$segments[1])
    samp <- app_query("SELECT lineage_taxids FROM samples WHERE ID = ?", list(u$ID))
    cols <- chimera_seg_colors(s, if (nrow(samp)) samp$lineage_taxids[1] else NA_character_)
    s <- s[, intersect(c("start", "stop", "taxon", "species", "pident", "n_windows"), names(s))]
    reactable::reactable(s, compact = TRUE, pagination = FALSE,
      defaultColDef = reactable::colDef(minWidth = 70),
      rowStyle = function(index) list(background = grDevices::adjustcolor(cols[index], alpha.f = CHIMERA_SEG_ALPHA)),
      columns = list(
        start = reactable::colDef(width = 80, html = TRUE, header = rt_help_header("Start", CHIMERA_COL_HELP$start)),
        stop = reactable::colDef(width = 80, html = TRUE, header = rt_help_header("Stop", CHIMERA_COL_HELP$stop)),
        taxon = reactable::colDef(minWidth = 140, html = TRUE, header = rt_help_header("Taxon", CHIMERA_COL_HELP$taxon)),
        species = reactable::colDef(minWidth = 220, na = "", html = TRUE, header = rt_help_header("Species", CHIMERA_COL_HELP$species)),
        pident = reactable::colDef(width = 100, format = reactable::colFormat(digits = 1), html = TRUE,
          header = rt_help_header("Identity %", CHIMERA_COL_HELP$pident)),
        n_windows = reactable::colDef(width = 100, html = TRUE, header = rt_help_header("Windows", CHIMERA_COL_HELP$n_windows))))
  })
  output$chimera_windows <- reactable::renderReactable({
    refresh()
    u <- hit_unit()
    w <- app_query("SELECT start, stop, species, lineage, pident, bitscore, rank_distance, depth_median, clip_frac FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      list(u$ID, u$path, u$scaffold))
    req(nrow(w) > 0)
    w$lineage_short <- vapply(w$lineage, function(l) {
      if (is.na(l)) return(NA_character_)
      p <- strsplit(l, ";", fixed = TRUE)[[1]]; p <- p[nzchar(p)]
      paste(utils::tail(p, 3), collapse = " > ")
    }, character(1), USE.NAMES = FALSE)
    w <- w[, c("start", "stop", "species", "lineage_short", "pident", "bitscore", "rank_distance", "depth_median", "clip_frac", "lineage")]
    # each window belongs to the segment whose start is the last one at or before its own
    ch <- app_query("SELECT segments FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    samp <- app_query("SELECT lineage_taxids FROM samples WHERE ID = ?", list(u$ID))
    cols <- if (nrow(ch) == 1 && !is.na(ch$segments[1])) {
      s <- jsonlite::fromJSON(ch$segments[1])
      chimera_seg_colors(s, if (nrow(samp)) samp$lineage_taxids[1] else NA_character_)[findInterval(w$start, s$start)]
    } else rep("#FFFFFF", nrow(w))
    reactable::reactable(w, compact = TRUE, pagination = TRUE, defaultPageSize = 15, searchable = TRUE,
      defaultColDef = reactable::colDef(na = "", minWidth = 70, headerStyle = list(whiteSpace = "nowrap")),
      rowStyle = function(index) list(background = grDevices::adjustcolor(cols[index], alpha.f = CHIMERA_SEG_ALPHA)),
      columns = list(
        start = reactable::colDef(width = 75, html = TRUE, header = rt_help_header("Start", CHIMERA_COL_HELP$start)),
        stop = reactable::colDef(width = 75, html = TRUE, header = rt_help_header("Stop", CHIMERA_COL_HELP$stop)),
        species = reactable::colDef(minWidth = 150, html = TRUE, header = rt_help_header("Species", CHIMERA_COL_HELP$win_species)),
        lineage_short = reactable::colDef(minWidth = 180, html = TRUE, header = rt_help_header("Lineage", CHIMERA_COL_HELP$lineage),
          cell = htmlwidgets::JS("function(c){ if (c.value == null) return ''; var full = c.row.lineage || ''; return '<span title=\"' + String(full).replace(/\"/g, '&quot;') + '\">' + c.value + '</span>'; }")),
        lineage = reactable::colDef(show = FALSE),
        pident = reactable::colDef(width = 100, format = reactable::colFormat(digits = 1), html = TRUE,
          header = rt_help_header("Identity %", CHIMERA_COL_HELP$win_pident)),
        bitscore = reactable::colDef(width = 90, format = reactable::colFormat(digits = 0), html = TRUE,
          header = rt_help_header("Bitscore", CHIMERA_COL_HELP$bitscore)),
        rank_distance = reactable::colDef(width = 95, html = TRUE, header = rt_help_header("Rank dist", CHIMERA_COL_HELP$rank_distance)),
        depth_median = reactable::colDef(width = 80, format = reactable::colFormat(digits = 0), html = TRUE,
          header = rt_help_header("Depth", CHIMERA_COL_HELP$depth_median)),
        clip_frac = reactable::colDef(width = 90, format = reactable::colFormat(digits = 2), html = TRUE,
          header = rt_help_header("Clip frac", CHIMERA_COL_HELP$clip_frac))))
  })

  # Annotate modal header badges: reviewed / problem / partial (live via refresh()).
  # Header flags for the open annotate unit: passive pills mirrored by the
  # toggle buttons beside them (MitoPilot's status_badges / status_toggles).
  modal_flags <- reactive({
    refresh()
    u <- sel_unit(); req(u)
    a <- app_query("SELECT reviewed, problematic, partial, annotate_lock FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    req(nrow(a) == 1)
    a
  })
  output$modal_status_badges <- renderUI({
    a <- modal_flags()
    pill <- function(label, state, yes_tone = "success") {
      val <- if (is.na(state)) "not set" else tolower(as.character(state))
      tone <- if (identical(val, "yes")) yes_tone else "neutral"
      span(class = paste0("mp-pill mp-pill-", tone), paste0(label, ": ", val))
    }
    tagList(
      pill("Reviewed", a$reviewed[1]),
      pill("Problematic", a$problematic[1], yes_tone = "warning"),
      pill("Partial", a$partial[1], yes_tone = "warning"))
  })
  output$modal_status_toggles <- renderUI({
    a <- modal_flags()
    # mp_flag_next() picks the value the click writes and the label that names
    # it, so the two never disagree. Coloured while the flag is set.
    toggle_btn <- function(id, label, state, title, tone = "success") {
      set <- identical(mp_flag_next(state), "yes")
      btn <- actionButton(id, paste(if (set) "Mark" else "Clear", label),
        icon = icon(if (set) "check" else "minus"),
        class = paste("btn-sm btn-default", if (!set) paste0("mp-btn-", tone)), title = title)
      if (isTRUE(a$annotate_lock[1] == 1)) shinyjs::disabled(btn) else btn
    }
    tagList(
      toggle_btn("mark_reviewed", "Reviewed", a$reviewed[1],
        "Record whether this annotation has been reviewed"),
      toggle_btn("mark_problematic", "Problematic", a$problematic[1],
        "Flag this annotation for another look", tone = "warning"),
      toggle_btn("mark_partial", "Partial", a$partial[1],
        "Flag this annotation as incomplete", tone = "warning"))
  })

  # active feature's unit-sequence coords + context
  nudge_ctx <- reactive({
    edit_tick()
    u <- sel_unit(); el <- active_feature()
    req(u, el)
    a <- app_query("SELECT start, stop, length, type FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(uargs(), list(el)))
    if (nrow(a) == 0) return(NULL)
    seq <- app_query("SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$sequence[1]
    list(el = el, start = a$start[1], stop = a$stop[1], type = a$type[1], seq = seq %||% "")
  })

  output$nudge_seq <- renderUI({
    cx <- nudge_ctx()
    if (is.null(cx) || !nzchar(cx$seq)) {
      return(em("Select a feature (no unit sequence available)."))
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
    u <- sel_unit(); el <- active_feature(); req(u, el)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    tryCatch(compute_feature_msa(con, u$ID, u$path, u$scaffold, el), error = function(e) NULL)
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
        "reference fetch from the pipeline run; if those ran, run Update again.")))
    }
    div(class = "rp-msa-box", msaR::msaROutput("feat_msa", height = "auto"))
  })
  output$feat_msa <- msaR::renderMsaR({
    input$msa_relayout  # redraw once the modal is laid out (custom.js)
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
    edit_tick(); u <- sel_unit(); req(u)
    app_query("SELECT element, type, start, stop, length, strand, status, flags FROM annotations
       WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", params = uargs())
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
        flags = reactable::colDef(name = "Warnings", html = TRUE, cell = rt_longtext())),
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
    edit_tick(); u <- sel_unit(); req(u)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND status = 'flag' ORDER BY start",
      params = uargs())
  })
  output$curate_notes <- renderUI({
    u <- sel_unit(); req(u)
    a <- app_query("SELECT n_flags, annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(a) == 0 || is.na(a$n_flags[1]) || a$n_flags[1] == 0) {
      return(span(class = "mp-fg-success", HTML("&#10003; No curation warnings for this sample.")))
    }
    tagList(
      span(sprintf("%d sample warning(s): ", a$n_flags[1]), class = "mp-fg-warning", style = "font-weight:bold;"),
      span(a$annotate_notes[1] %||% "", class = "mp-fg-warning"))
  })
  output$flag_tbl <- reactable::renderReactable({
    d <- flag_rows()
    if (nrow(d) == 0) {
      return(reactable::reactable(data.frame(` ` = "No features with warnings.", check.names = FALSE), compact = TRUE))
    }
    d$jump <- "select"
    reactable::reactable(d[, c("element", "flags", "jump")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 10,
      columns = list(
        element = reactable::colDef(name = "Feature"),
        flags = reactable::colDef(name = "Warnings", html = TRUE, cell = rt_longtext()),
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
      u <- sel_unit()
      for (i in seq_len(nrow(b))) {
        DBI::dbExecute(con,
          "UPDATE annotations SET start = ?, stop = ?, length = ? WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
          params = list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, u$ID, u$path, u$scaffold, b$element[i]))
      }
    }
  }
  save_edits <- function() {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # refresh flags after edits
    edit_backup(NULL); edit_dirty(FALSE)
    edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    rp_toast("Boundary edits saved.", type = "message")
  }
  observeEvent(input$feat_edit, ignoreInit = TRUE, {
    if (isTRUE(input$feat_edit)) {
      if (is.null(active_feature())) {
        shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
        rp_toast("Select a feature first.", type = "warning"); return()
      }
      edit_backup(app_query("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs()))
      edit_active(TRUE); edit_dirty(FALSE)
      shinyjs::show("feat_edit_ctrls")
    } else {
      edit_active(FALSE); shinyjs::hide("feat_edit_ctrls")
      if (isTRUE(edit_dirty())) {
        # MitoPilot-style prompt: keep or discard unsaved boundary edits
        mp_choice("unsaved_confirm", title = "Unsaved boundary edits",
          text = "You have unsaved changes to this feature's boundaries.",
          labels = c("Discard", "Save"))
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
      rp_toast("Boundary edits discarded.", type = "message")
    }
  }, ignoreNULL = TRUE)
  observeEvent(input$feat_save, {
    save_edits()
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)  # clean exit (not dirty)
  })
  observeEvent(input$feat_reset, {
    restore_edit_backup(); edit_dirty(FALSE); edit_tick(edit_tick() + 1)
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
    rp_toast("Boundary edits discarded.", type = "message")
  })

  # helper: (re)open the annotate details modal on the current unit
  show_annotate_unit <- function() {
    u <- sel_unit(); req(u); showModal(annotate_modal(u$ID, u$path, u$scaffold))
  }

  # ---- trim the active unit to annotated / custom boundaries ----
  observeEvent(input$trim_btn, {
    u <- sel_unit(); req(u)
    if (edits_block()) return()
    s <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(s) == 0 || is.na(s$length[1])) {
      rp_alert("No unit sequence to trim.", type = "warning"); return()
    }
    span <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    a_from <- if (nrow(span) && !is.na(span$a[1])) span$a[1] else 1L
    a_to <- if (nrow(span) && !is.na(span$b[1])) span$b[1] else s$length[1]
    showModal(trim_modal(s$length[1], a_from, a_to))
  })

  # From/To are the single source of truth for the kept window. The annotated
  # and flank modes are presets that fill them in; they are editable only in
  # custom mode (brushing the plot switches to custom).
  trim_span <- reactive({
    u <- sel_unit(); req(u)
    sp <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    n <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
    n <- if (length(n)) n[1] else NA_integer_
    c(a = if (nrow(sp) && !is.na(sp$a[1])) sp$a[1] else 1L,
      b = if (nrow(sp) && !is.na(sp$b[1])) sp$b[1] else n, n = n)
  })
  trim_preset <- reactive({
    mode <- input$trim_mode %||% "annotated"
    if (mode == "custom") return(NULL)
    sp <- trim_span()
    x <- if (mode == "flank") suppressWarnings(as.integer(input$trim_flank)) else 0L
    if (is.na(x) || x < 0) x <- 0L
    c(max(1L, sp[["a"]] - x), min(sp[["n"]], sp[["b"]] + x))
  })
  observeEvent(trim_preset(), {
    p <- trim_preset(); req(p)
    updateNumericInput(session, "trim_from", value = p[1])
    updateNumericInput(session, "trim_to", value = p[2])
  })
  observeEvent(input$trim_mode, {
    shinyjs::toggleState("trim_from", condition = input$trim_mode == "custom")
    shinyjs::toggleState("trim_to", condition = input$trim_mode == "custom")
  })
  output$trim_flank_hint <- renderText({
    req((input$trim_mode %||% "annotated") == "flank")
    sp <- trim_span(); p <- trim_preset(); req(p)
    x <- suppressWarnings(as.integer(input$trim_flank))
    if (is.na(x) || x < 1) return("Enter a flank of at least 1 bp.")
    clipped <- c(if (p[1] > sp[["a"]] - x) sprintf("upstream clipped to %d bp (scaffold start)", sp[["a"]] - p[1]),
      if (p[2] < sp[["b"]] + x) sprintf("downstream clipped to %d bp (scaffold end)", p[2] - sp[["b"]]))
    if (length(clipped)) paste0("Flank ", paste(clipped, collapse = "; "), ".") else ""
  })
  trim_window <- reactive({
    c(suppressWarnings(as.integer(input$trim_from)), suppressWarnings(as.integer(input$trim_to)))
  })
  output$trim_plot <- renderPlot({
    edit_tick(); u <- sel_unit(); req(u); w <- trim_window(); plot_trim(u$ID, u$path, u$scaffold, w[1], w[2])
  }, height = 230)
  output$trim_chimera_hint <- renderText({
    u <- sel_unit(); req(u)
    ch <- app_query("SELECT junctions FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    req(nrow(ch) == 1, nzchar(ch$junctions[1] %|NA|% ""))
    sprintf("suggested junctions: %s", ch$junctions[1])
  })
  # Dragging a brush over the plot sets the custom window from the x-range.
  observeEvent(input$trim_brush, {
    b <- input$trim_brush; u <- sel_unit(); req(b, u)
    n <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
    n <- if (length(n)) n[1] else NA
    lo <- max(1L, round(b$xmin)); hi <- round(b$xmax)
    if (!is.na(n)) hi <- min(n, hi)
    if (is.na(hi) || hi <= lo) return()
    shinyWidgets::updatePrettyRadioButtons(session, "trim_mode", selected = "custom")
    updateNumericInput(session, "trim_from", value = lo)
    updateNumericInput(session, "trim_to", value = hi)
  })

  observeEvent(input$trim_apply, {
    u <- sel_unit(); req(u)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    w <- trim_window(); from <- w[1]; to <- w[2]; n <- trim_span()[["n"]]
    if (is.na(from) || is.na(to) || from < 1 || from > to || (!is.na(n) && to > n)) {
      rp_alert("Invalid trim range.", type = "error"); return()
    }
    ok <- tryCatch({
      trim_assembly(con, u$ID, u$path, u$scaffold, from, to)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # re-curate after the shift
      TRUE
    }, error = function(e) { rp_alert(conditionMessage(e), type = "error"); FALSE })
    if (isTRUE(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()  # back to the details modal, now showing trimmed coords
      rp_alert("Unit trimmed; feature coordinates + coverage shifted.", type = "message")
    }
  })

  # Cancel out of the trim modal: return to the annotate details modal.
  observeEvent(input$trim_cancel, { req(sel_unit()); show_annotate_unit() })

  # ---- restore the active unit to its original (pre-trim) length ----
  observeEvent(input$restore_btn, {
    u <- sel_unit(); req(u)
    if (edits_block()) return()
    len <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length[1]
    mp_confirm("restore_confirm", title = "Restore the full assembly?",
      text = sprintf("The unit is %s bp now. Restoring undoes the trim and re-runs curation.", len %||% "?"),
      action_label = "Restore", danger = TRUE)
  })
  observeEvent(input$restore_confirm, {
    req(isTRUE(input$restore_confirm))
    u <- sel_unit(); req(u)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ok <- tryCatch({
      len <- restore_assembly(con, u$ID, u$path, u$scaffold)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # recompute flags for the full unit
      len
    }, error = function(e) { rp_alert(conditionMessage(e), type = "warning"); NULL })
    if (!is.null(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()
      rp_alert(paste0("Unit restored to original length (", ok, " bp)."), type = "message")
    }
  })

  # ---- notes ----
  observeEvent(input$annotate_notes, ignoreInit = TRUE, {
    u <- sel_unit(); req(u)
    new <- input$annotate_notes %||% ""
    cur <- app_query("SELECT annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$annotate_notes[1]
    cur <- (cur %||% NA_character_) %|NA|% ""
    if (identical(new, cur)) return()
    # A freshly rendered modal reports "" before the stored value lands; never
    # let that blank a stored note. Clearing text is done with the Clear button.
    if (!nzchar(new) && nzchar(cur)) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(new), uargs()))
  })
  observeEvent(input$annotate_notes_clear, {
    u <- sel_unit(); req(u)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    updateTextAreaInput(session, "annotate_notes", value = "")
    rp_toast("Notes cleared.", type = "message")
  })

  # nudge one boundary by the chosen step (1-50), clamp to sequence, persist
  nudge_boundary <- function(end, dir) {
    if (!isTRUE(edit_active())) { rp_toast("Check Edit before adjusting boundaries.", type = "warning"); return() }
    cx <- nudge_ctx(); req(cx)
    step <- suppressWarnings(as.integer(input$feat_step %||% 1))
    if (is.na(step)) step <- 1
    step <- max(1L, min(50L, step))
    n <- nchar(cx$seq); start <- cx$start; stop <- cx$stop
    if (end == "5") start <- start + if (dir == "out") -step else step
    else stop <- stop + if (dir == "out") step else -step
    start <- max(1, start); stop <- min(if (n > 0) n else stop, stop)
    if (start > stop) { rp_toast("Boundaries would cross; nudge ignored.", type = "warning"); return() }
    u <- sel_unit()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = list(start, stop, stop - start + 1, u$ID, u$path, u$scaffold, cx$el))
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
  }
  observeEvent(input$nt_5_out, nudge_boundary("5", "out"))
  observeEvent(input$nt_5_in,  nudge_boundary("5", "in"))
  observeEvent(input$nt_3_out, nudge_boundary("3", "out"))
  observeEvent(input$nt_3_in,  nudge_boundary("3", "in"))
  # per-unit curator flag toggles + Lock & Close in the annotate modal
  modal_flag_toggle <- function(col, label) {
    u <- sel_unit(); req(u)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, sprintf(
      "UPDATE annotate SET %s = CASE WHEN %s = 'yes' THEN 'no' ELSE 'yes' END WHERE ID = ? AND path = ? AND scaffold = ?",
      col, col), params = uargs())
    refresh(refresh() + 1)
    rp_toast(paste0("Toggled ", label, "."), type = "message")
  }
  observeEvent(input$mark_reviewed, modal_flag_toggle("reviewed", "reviewed"))
  observeEvent(input$mark_problematic, modal_flag_toggle("problematic", "problematic"))
  observeEvent(input$mark_partial, modal_flag_toggle("partial", "partial"))
  # Unsaved boundary edits stop any exit from the annotation window.
  edits_block <- function() {
    if (!isTRUE(edit_dirty())) return(FALSE)
    rp_alert("Save or reset your boundary edits before leaving this window.", type = "warning")
    TRUE
  }
  observeEvent(input$annotate_close, {
    if (edits_block()) return()
    removeModal()
  })
  observeEvent(input$lock_and_close, {
    u <- sel_unit(); req(u)
    if (edits_block()) return()
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    # One assembly path per sample may be locked (paths are alternative
    # resolutions of the same repeat); scaffolds of that path can all be locked.
    other <- DBI::dbGetQuery(con,
      "SELECT DISTINCT path FROM annotate WHERE ID = ? AND annotate_lock = 1 AND path != ?",
      params = list(u$ID, u$path))
    if (nrow(other) > 0) {
      rp_alert(sprintf(
        "Path %s is already locked for %s. Only one assembly path can be locked per sample; unlock it first.",
        other$path[1], u$ID), type = "error")
      return()
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET annotate_lock = 1, reviewed = 'yes' WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    refresh(refresh() + 1)
    removeModal()
    rp_toast("Locked annotation unit.", type = "message")
  })
}
