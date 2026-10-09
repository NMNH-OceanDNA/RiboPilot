#' RiboPilot app server
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

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

  # Container version check. Warn only: the project still runs, with the
  # RiboPilot code and tools in its image.
  cgap <- container_version_gap(app_project())
  if (!is.null(cgap)) {
    mp_alert(
      title = "Container version does not match RiboPilot",
      text = tags$div(
        tags$p("This project runs the pipeline from:"),
        tags$pre(cgap$configured),
        tags$p("but the installed RiboPilot is ",
          tags$b(as.character(utils::packageVersion("RiboPilot"))), ", which expects:"),
        tags$pre(cgap$expected),
        tags$p("The app works either way, but pipeline runs will use the RiboPilot code ",
          "and tools in the older image."),
        tags$p(tags$b("To switch, run backwards_compatibility() on this project."),
          " It changes only the container line and leaves the rest of your settings alone.")),
      html = TRUE, type = "info")
    session$userData$startup_alert <- TRUE
  }

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
  run_reports_handlers(input, output, session, app_project, refresh = reactive(input$refresh))
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
    mode <- input$mode %||% "Taxonomy"
    noun <- if (mode == "Reference" && ref_view() == "Groups") "group" else MP_UNIT[[tolower(mode)]]
    mp_table_status(min(length(visible_rows()), nrow(samples_data())), nrow(all_rows()),
      length(selected_ids()), noun)
  })
  observe(shinyjs::toggleState(selector = ".mp-needs-selection",
    condition = length(selected_ids()) > 0))
  table_csv_handler <- function(suffix, rows, data = samples_data) downloadHandler(
    filename = function() paste0(basename(app_project()), "_",
      tolower(input$mode %||% "table"), suffix, ".csv"),
    content = function(file) {
      d <- data()
      d <- d[rows(), setdiff(names(d), c("details", "output", "action", "checks_tone", "specimen_message", "specimen_icons")), drop = FALSE]
      utils::write.csv(d, file, row.names = FALSE)
    }
  )
  output$table_csv <- table_csv_handler("", function() TRUE, all_rows)
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
  meta_view_ver <- reactiveVal(0L)
  meta_view_open <- meta_view_setup(input, output, session, con = app_con, ver = meta_view_ver)
  all_rows <- reactive({
    meta_view_ver()
    d <- all_rows_base()
    if ((input$mode %||% "Assemble") %in% c("Taxonomy", "Assemble", "Annotate", "Export") &&
        "Taxon" %in% names(d)) {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      d <- .specimen_status_join(d, con) |>
        dplyr::relocate(specimen, specimen_message, specimen_icons, .after = Taxon)
      d <- meta_view_join(d, con)
    }
    d
  })
  all_rows_base <- reactive({
    refresh()
    mode <- input$mode %||% "Assemble"
    # Column order mirrors MitoPilot: lock, state, ID, Taxon lead the tables.
    if (mode == "Taxonomy") {
      app_query(
        "SELECT s.taxonomy_lock, s.taxonomy_switch, s.ID, s.Taxon,
           s.sci_name, s.taxid, s.rank, s.match_type, s.tax_alternatives, s.time_stamp
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
           -- has been checked, else 'n of m contigs' flagged
           (SELECT CASE WHEN COUNT(x.itsx_warnings) = 0 THEN NULL ELSE
              SUM(CASE WHEN COALESCE(x.itsx_warnings,'') <> ''
                       AND x.itsx_warnings NOT LIKE 'not_run(%' THEN 1 ELSE 0 END)
              || ' of ' || COUNT(x.itsx_warnings) || ' contigs' END
            FROM annotate x WHERE x.ID = s.ID) AS itsx_warnings,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           -- the rDNA contig search writes its own column because
           -- load_assemblies() overwrites assemble_notes right after it; both
           -- are shown as one Notes cell
           TRIM(COALESCE(a.assemble_notes,'') ||
             CASE WHEN COALESCE(a.find_rrna_notes,'') = '' THEN ''
                  WHEN COALESCE(a.assemble_notes,'') = '' THEN a.find_rrna_notes
                  ELSE '; ' || a.find_rrna_notes END) AS notes, a.user_notes, a.time_stamp,
           'details' AS details, 'output' AS output
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID
           LEFT JOIN preprocess pp ON s.ID = pp.ID
         WHERE s.taxonomy_lock = 1") |> with_checks("sample")
    } else if (mode == "Export") {
      # the units an export writes (EXPORTABLE_WHERE), one row per (ID, path, scaffold)
      app_query(paste(
        "SELECT s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ", s.Taxon, s.sci_name, s.taxid, asm.length,
           an.n_features, an.missing,
           an.n_flags, an.curate_flags, an.annotate_notes, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason,
           an.reviewed, an.partial,
           an.export_group,", UNIT_EXPORTED_SQL, ", an.time_stamp
         FROM annotate an JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE", EXPORTABLE_WHERE))
    } else {
      # annotate shows one row per non-ignored assembly UNIT of an assembly-locked
      # sample (one (ID, path, scaffold) each, annotated independently)
      app_query(paste(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, an.path, an.scaffold,", UNIT_MULTI_SQL, ",
           s.Taxon, s.sci_name, s.taxid,
           asm.length, an.n_features, an.missing,
           -- blast_source says WHICH search named the species: the local seed
           -- BLAST, the local panel BLAST, or the remote core_nt one. Without
           -- it a seed or panel hit reads as a GenBank identification.
           -- blast_qcovs is written by the remote and local-panel paths only
           -- and its column was already defined but never selected.
           asm.blast_species, asm.blast_accession, asm.blast_pident, asm.blast_qcovs,
           asm.blast_source,
           an.n_flags, an.curate_flags, an.itsx_warnings, sc.screen_verdict,
           ch.label AS chimera_label, ch.reason AS chimera_reason, an.reviewed,
           an.problematic, an.partial, an.annotate_notes,
           an.export_group,", UNIT_EXPORTED_SQL, ", an.time_stamp,
           'details' AS details, 'output' AS output
         FROM annotate an
           JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
           LEFT JOIN screen sc ON sc.ID = an.ID AND sc.path = an.path AND sc.scaffold = an.scaffold
           LEFT JOIN chimera ch ON ch.ID = an.ID AND ch.path = an.path AND ch.scaffold = an.scaffold
         WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0")) |> with_checks("unit")
    }
  })

  # The rows the table shows: the step's Lock and State filters are applied here,
  # so the table pages and counts only the rows that pass them.
  samples_data <- reactive(all_rows()[visible_mask(), , drop = FALSE])

  # The Checks summary goes right after the row's ID (and path and scaffold).
  with_checks <- function(d, by) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ck <- unit_checks(con, by)
    k <- function(x) if (by == "unit") unit_key(x$ID, x$path, x$scaffold) else x$ID
    i <- match(k(d), k(ck))
    at <- max(match(c("ID", "path", "scaffold"), names(d)), na.rm = TRUE)
    cbind(d[seq_len(at)], checks = ck$checks[i], checks_tone = ck$checks_tone[i], d[-seq_len(at)])
  }

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  # Re-render only on structural changes (mode / column groups); the data is
  # isolated so a refresh() pushes a data-only update below instead of a full
  # re-render, which would reset the user's column filters, sort, and page.
  # Column defs read the data once; these two re-render the table only when a
  # refresh() changes what the defs depend on (unit columns, group members).
  multi_unit_rv <- reactiveVal(FALSE)
  # row keys of the data the table was last given
  tbl_keys <- character(0)
  # TRUE from a re-render until the new table reports its own page and selection
  tbl_stale <- FALSE
  row_keys <- function(d) unit_key(d$ID, d[["path"]] %||% 1L, d[["scaffold"]] %||% 1L)
  observe({ d <- all_rows(); multi_unit_rv(any(duplicated(d$ID))) })
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    multi_unit_rv()
    meta_view_ver()
    if (mode == "Reference" && ref_view() == "Groups") refresh()
    d <- isolate(samples_data())
    tbl_keys <<- row_keys(d)
    st <- isolate(reactable::getReactableState("tbl"))
    tbl_stale <<- length(st$selected) > 0 || !(st$page %||% 1) %in% 1
    mm <- mode_meta(mode)
    lock_col <- paste0(mm$pfx, "_lock")
    state_col <- paste0(mm$pfx, "_switch")

    # Details/Output icon buttons open the mode's detail modal + the work dirs.
    det_input <- if (mode == "Annotate") "ann_details" else "asm_details"
    det_tip <- if (mode == "Annotate") "View annotation details" else "View assembly details"
    out_input <- if (mode == "Annotate") "ann_output" else "asm_output"

    # Path/Scaffold are shown only when the sample has more than one unit; a
    # single-unit project reads exactly like the old per-sample table.
    a <- isolate(all_rows())
    multi_unit <- all(c("path", "scaffold") %in% names(a)) &&
      nrow(unique(a[, c("ID", "path", "scaffold")])) > length(unique(a$ID))

    cd <- function(f, ...) reactable::colDef(name = mp_col_name(f), header = mp_col_header(f), ...)
    cols <- list(
      ID = cd("ID", width = mp_fit_width(a$ID), sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      path = cd("path", width = 70, align = "center", show = multi_unit, html = TRUE,
        cell = rt_multi_star("path_n", "kept paths"), filterMethod = rt_multi_star_filter("path_n")),
      scaffold = cd("scaffold", width = 90, align = "center", show = multi_unit, html = TRUE,
        cell = rt_multi_star("scaffold_n", "kept contigs"),
        filterMethod = rt_multi_star_filter("scaffold_n")),
      path_n = reactable::colDef(show = FALSE),
      scaffold_n = reactable::colDef(show = FALSE),
      Taxon = cd("Taxon", minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      sci_name = cd("sci_name", minWidth = 200, html = TRUE, cell = rt_taxid_link()),
      specimen = specimen_col_def("specimen_open", sticky = "left"),
      specimen_message = reactable::colDef(show = FALSE),
      specimen_icons = reactable::colDef(show = FALSE),
      taxid = reactable::colDef(show = FALSE),
      rank = cd("rank", minWidth = 90, html = TRUE, cell = rt_longtext()),
      match_type = cd("match_type", minWidth = 110, html = TRUE, align = "center",
        cell = rt_confidence_badge()),
      tax_alternatives = cd("tax_alternatives", minWidth = 260, html = TRUE, cell = rt_longtext()),
      checks = cd("checks", minWidth = 180, html = TRUE, cell = rt_checks_badge()),
      checks_tone = reactable::colDef(show = FALSE),
      reviewed = cd("reviewed", html = TRUE, align = "center", cell = rt_bool_badge(hide_no = TRUE)),
      problematic = cd("problematic", width = 110, html = TRUE, align = "center",
        cell = rt_bool_badge(hide_no = TRUE, yes_tone = "warning")),
      partial = cd("partial", width = 80, html = TRUE, align = "center",
        cell = rt_bool_badge(yes_tone = "warning", no_tone = "success")),
      n_flags = cd("n_flags", width = 100, html = TRUE, filterable = FALSE, align = "center",
        cell = rt_flag_badge()),
      annotate_notes = cd("annotate_notes", minWidth = 160, html = TRUE, class = "mp-note-cell", cell = rt_longtext()),
      curate_flags = reactable::colDef(show = FALSE),
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
      units = cd("units", width = 95),
      kept = cd("kept", width = 70, html = TRUE, cell = rt_kept_count()),
      sel_path = cd("sel_path", width = 110),
      best_score = cd("best_score", width = 100),
      export_group = cd("export_group", minWidth = 110, html = TRUE, cell = rt_longtext()),
      exported = cd("exported", width = 170, html = TRUE, cell = rt_text_pill("success")),
      pre_opts = cd("pre_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_pre_opts", lock_col = "assemble_lock")),
      assemble_opts = cd("assemble_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_assemble_opts", lock_col = "assemble_lock")),
      blast_opts = cd("blast_opts", minWidth = 100, html = TRUE, filterable = FALSE,
        cell = rt_link("set_blast_opts", lock_col = "assemble_lock")),
      # guarded by the annotation lock on click (app_opts.R), not the assembly lock
      annotate_opts = cd("annotate_opts", minWidth = 120, html = TRUE, filterable = FALSE,
        cell = rt_link("set_annotate_opts")),
      curate_opts = cd("curate_opts", minWidth = 100, html = TRUE, filterable = FALSE,
        cell = rt_link("set_curate_opts")),
      annotate_lock = cd("annotate_lock", width = 34, align = "center", filterable = FALSE,
        html = TRUE, cell = rt_lock_icon()),
      notes = cd("notes", minWidth = 160, html = TRUE, class = "mp-note-cell", cell = rt_longtext()),
      user_notes = cd("user_notes", minWidth = 160, html = TRUE, class = "mp-note-cell", cell = rt_longtext()),
      time_stamp = cd("time_stamp", width = 150, filterable = FALSE, cell = rt_ts_date())
    )
    if (lock_col %in% names(d)) {
      cols[[lock_col]] <- cd(lock_col, sticky = "left", html = TRUE, width = 52,
        align = "center", filterable = FALSE, cell = rt_lock_icon("lock_row"))
    }
    if (state_col %in% names(d) && mm$pfx %in% names(MP_STATE_CODES)) {
      cols[[state_col]] <- cd(state_col, sticky = "left", html = TRUE, width = 54,
        align = "center", filterable = FALSE,
        cell = rt_dynamicIcon(icons = mp_state_icon_classes(mm$pfx), labels = mp_state_titles(mm$pfx),
          inputId = "state_row"))
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
      cols[[state_col]] <- reactable::colDef(show = FALSE)
    }
    # Explain an empty table: every mode after Taxonomy is gated on the previous
    # stage's lock, and a fresh project has nothing locked, so the default
    # "No rows found" reads like a failed import.
    nodata <- switch(mode,
      Taxonomy = "No samples yet. Create the project with a mapping file (new_project()).",
      Reference = if (ref_view() == "Samples") {
        "No samples yet. Create the project with a mapping file (new_project())."
      } else {
        paste("No reference groups yet. Press Propose groups to group your samples by clade",
          "from the assigned taxonomy, then set up a database for each group.")
      },
      Assemble = paste("No samples to show. A sample is listed once its taxonomy is locked",
        "in the Taxonomy step. The Lock and State filters and the search box also hide rows."),
      Annotate = paste("No contigs to show. A contig is listed once its sample's assembly is locked",
        "in the Assemble step. The Lock and State filters and the search box also hide rows."),
      Export   = paste("No contigs ready to export. A contig is listed once its assembly (Assemble)",
        "and its annotation (Annotation Review) are both locked; contigs excluded from export",
        "are not listed."),
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
        `Check details` = c("n_flags", "itsx_warnings", "screen_verdict", "chimera_label")),
      NULL)
    keep_cols <- cols[intersect(names(cols), names(d))]
    if (mode != "Reference") {
      mcon <- app_con()
      on.exit(DBI::dbDisconnect(mcon), add = TRUE)
      keep_cols <- c(keep_cols, meta_view_table_defs(mcon, d))
    }
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
      rowStyle = rt_highlight_row(),
      # column order comes from samples_data()
      columns = keep_cols
    )
  })

  # Data-only updates keep the search, column filters, and sort; the page and
  # the selection (by row key) are passed back. Mode and column-group changes
  # still re-render via output$tbl above.
  push_tbl <- function(page = NULL) {
    d <- samples_data(); k <- row_keys(d)
    sel <- if (!tbl_stale) reactable::getReactableState("tbl", "selected")
    reactable::updateReactable("tbl", data = d, selected = which(k %in% tbl_keys[sel]),
      page = if (!tbl_stale) page)
    tbl_keys <<- k
    session$sendCustomMessage("rp_tbl_report", TRUE)
  }
  observeEvent(list(reactable::getReactableState("tbl", "selected"),
    reactable::getReactableState("tbl", "page")), tbl_stale <<- FALSE, ignoreInit = TRUE)
  observeEvent(refresh(), push_tbl(reactable::getReactableState("tbl", "page")), ignoreInit = TRUE)
  observeEvent(list(input$Assemble_lock_filter, input$Assemble_state_filter,
    input$Annotate_lock_filter, input$Annotate_state_filter), push_tbl(), ignoreInit = TRUE)

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
        rp_alert("Lock the assembly in the Assemble step before annotating this sample.",
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

  # logical mask of all_rows() that pass the active mode's filters
  visible_mask <- reactive({
    d <- all_rows()
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

  # Rows left by the table's search and column filters, as the browser reports
  # them (custom.js); every row until it has.
  visible_rows <- function() {
    v <- input$tbl_visible
    if (is.null(v)) seq_len(nrow(samples_data())) else as.integer(unlist(v$rows))
  }
  # Bulk actions only ever touch rows the user can see.
  selected_rows <- function() {
    intersect(reactable::getReactableState("tbl", "selected") %||% integer(0), visible_rows())
  }
  selected_ids <- function() samples_data()$ID[selected_rows()]
  # Visible selected rows as (ID, path, scaffold) units (path/scaffold = 1 when the
  # active table has no unit columns).
  selected_units <- function(sel = selected_rows()) {
    if (!length(sel)) return(NULL)
    d <- samples_data()
    data.frame(ID = d$ID[sel],
      path = if ("path" %in% names(d)) as.integer(d$path[sel]) else 1L,
      scaffold = if ("scaffold" %in% names(d)) as.integer(d$scaffold[sel]) else 1L,
      stringsAsFactors = FALSE)
  }

  # Rows the search or a column filter hides leave the selection.
  observeEvent(input$tbl_visible, {
    keep <- selected_rows()
    if (length(keep) != length(reactable::getReactableState("tbl", "selected"))) {
      reactable::updateReactable("tbl", selected = keep)
    }
  })

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
      app_query(sprintf(
        "SELECT ID, Taxon, sci_name, taxid, rank, tax_alternatives FROM samples WHERE ID IN (%s)", ph),
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
  set_flag <- function(sql, verb = "updated", u = selected_units()) {
    if (identical(input$mode %||% "Assemble", "Annotate")) return(set_unit_flag(sql, verb, u))
    ids <- u$ID
    if (!need_selection(length(ids))) return()
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbWithTransaction(con,
      for (id in ids) for (q in sql) DBI::dbExecute(con, q, params = list(id))))) return()
    refresh(refresh() + 1)
    mp_toast(paste(mp_n(length(ids), "sample"), paste0(verb, ".")))
  }
  set_unit_flag <- function(sql, verb = "updated", u = selected_units()) {
    if (!need_selection(NROW(u))) return()
    usql <- sub("WHERE ID = \\?", "WHERE ID = ? AND path = ? AND scaffold = ?", sql)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbWithTransaction(con, for (i in seq_len(nrow(u)))
      DBI::dbExecute(con, usql, params = list(u$ID[i], u$path[i], u$scaffold[i]))))) return()
    refresh(refresh() + 1)
    mp_toast(paste(mp_n(nrow(u), "contig"), paste0(verb, ".")))
  }
  # One direction per click (mp_flag_next): a selection that is all locked
  # unlocks, anything else locks. Unlocking drops rows from the next step's
  # table, so it is confirmed first. lock_next holds what a confirmed prompt does.
  lock_next <- reactiveVal(NULL)
  observeEvent(input$lock_confirm, {
    f <- lock_next(); lock_next(NULL)
    req(isTRUE(input$lock_confirm), f)
    f()
  })
  # Units flagged as possible contamination and not marked Reviewed, asked
  # about before a lock, as the Export window does.
  confirm_flagged_lock <- function(f, then) {
    lock_next(then)
    one <- nrow(f) == 1
    mp_confirm("lock_confirm", title = "Possible contamination or chimera",
      text = paste0(mp_n(nrow(f), "contig"), if (one) " is" else " are",
        " flagged and not marked Reviewed: ",
        mp_id_list(sprintf("%s p%s.s%s (%s)", f$ID, f$path, f$scaffold,
          risk_reasons(f$submission_risk))),
        ". Lock anyway? Export asks about ", if (one) "it" else "them",
        " again unless ", if (one) "it is" else "they are", " marked Reviewed."),
      action_label = "Lock anyway", danger = TRUE)
  }
  # Excluded units are never exported, so Lock does not ask about them.
  lock_flagged <- function(u) {
    f <- screen_flagged(u)
    if (nrow(f)) f[!(f$problematic %in% "yes"), , drop = FALSE] else f
  }
  lock_sql <- function(mm, to) sprintf("UPDATE %s SET %s_lock = %d WHERE ID = ?", mm$tbl, mm$pfx, to)
  next_step <- c(taxonomy = "Assemble", assemble = "Annotation Review", annotate = "Export")
  lock_rows <- function(s) {
    mm <- mode_meta(input$mode)
    u <- selected_units(s)
    if (!need_selection(NROW(u))) return()
    lc <- paste0(mm$pfx, "_lock")
    d <- samples_data()
    cur <- d[[lc]][s]
    to <- if (identical(mp_flag_next(cur, on = 1L, off = 0L), 0L)) 0L else 1L
    n <- nrow(u)
    if (to == 0L) {
      # An unlocked sample leaves Annotation Review, where its units could no
      # longer be unlocked, so their annotation locks are released with it.
      asm <- identical(mm$pfx, "assemble")
      unlock <- function() set_flag(c(lock_sql(mm, 0L),
        if (asm) "UPDATE annotate SET annotate_lock = 0 WHERE ID = ?"), verb = "unlocked", u = u)
      skip <- mp_confirm("lock_confirm", title = paste0("Unlock ", mp_n(n, MP_UNIT[[mm$pfx]]), "?"),
        text = paste0(sprintf("Unlocked rows leave the %s table until locked again.", next_step[[mm$pfx]]),
          if (asm) " Their contigs are unlocked in Annotation Review too."),
        action_label = "Unlock", danger = TRUE, skippable = TRUE)
      if (skip) unlock() else lock_next(unlock)
      return()
    }
    do_lock <- function() set_flag(lock_sql(mm, 1L), verb = "locked", u = u)
    if (identical(mm$pfx, "taxonomy")) {
      weak <- d$ID[s][is.na(d$taxid[s])]
      if (length(weak) > 0) {
        lock_next(do_lock)
        mp_confirm("lock_confirm", title = "Lock samples with no taxon?",
          text = sprintf(paste("%s: %s. Without a taxon a sample gets no reference group,",
            "and the BLAST screen cannot judge it. Cancel and click the ID to pick",
            "a taxon, or lock anyway."),
            mp_n(length(weak), "selected sample has no taxon", "selected samples have no taxon"),
            paste(c(utils::head(weak, 10), if (length(weak) > 10) "..."), collapse = ", ")),
          action_label = "Lock anyway", danger = TRUE)
        return()
      }
    }
    if (identical(input$mode %||% "Assemble", "Annotate") && !need_one_path(u)) return()
    open <- !(cur %in% 1)
    ids <- unique(u$ID[open])
    units <- if (identical(mm$pfx, "annotate")) u[open, , drop = FALSE]
      else if (identical(mm$pfx, "assemble") && length(ids)) app_query(sprintf(
        "SELECT an.ID, an.path, an.scaffold FROM annotate an
         JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
         WHERE COALESCE(asm.ignore, 0) = 0 AND an.ID IN (%s)",
        paste(rep("?", length(ids)), collapse = ",")), params = as.list(ids))
    f <- lock_flagged(units)
    if (nrow(f)) return(confirm_flagged_lock(f, do_lock))
    do_lock()
  }
  observeEvent(list(input$tax_lock_btn, input$asmb_lock_btn, input$annot_lock_btn), ignoreInit = TRUE,
    lock_rows(selected_rows()))
  # Row icons send their index into samples_data(); search, filters, sort and
  # paging leave it unchanged.
  clicked_row <- function(i) {
    i <- suppressWarnings(as.integer(i))
    req(isTRUE(i >= 1 && i <= nrow(samples_data())))
    i
  }
  observeEvent(input$lock_row, lock_rows(clicked_row(input$lock_row)))
  state_units <- reactiveVal(NULL)
  open_state <- function(s) {
    if (!need_selection(length(s))) return()
    mm <- mode_meta(input$mode)
    state_units(selected_units(s))
    # preselect only when every selected row agrees
    cur <- unique(samples_data()[[paste0(mm$pfx, "_switch")]][s])
    showModal(mp_state_modal(samples_data()$ID[s],
      current = if (length(cur) == 1 && !is.na(cur)) as.character(cur) else character(0),
      mode = mm$pfx))
  }
  observeEvent(list(input$tax_state_btn, input$asmb_state_btn, input$annot_state_btn),
    ignoreInit = TRUE, open_state(selected_rows()))
  observeEvent(input$state_row, open_state(clicked_row(input$state_row)))
  observeEvent(input$update_state, {
    if (is.null(input$new_state)) {
      mp_toast("Choose a state first.", type = "warning"); return()
    }
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state)
    set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st), u = state_units())
    removeModal()
  })
  # Toolbar marks: the unlocked part of the selection goes one way
  # (mp_flag_next), as the toggles in the annotation window do; locked units
  # are left alone.
  toolbar_flag <- function(col) {
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    d <- samples_data(); sel <- selected_rows()
    lk <- d$annotate_lock[sel] %in% 1
    if (all(lk)) {
      need_unlocked(sprintf("%s p%s.s%s", u$ID, u$path, u$scaffold), "unit")
      return()
    }
    to <- mp_flag_next(d[[col]][sel][!lk])
    w <- MP_FLAG_WORDS[[col]]
    set_unit_flag(sprintf("UPDATE annotate SET %s = '%s' WHERE ID = ?", col, to),
      verb = w[[if (to == "yes") "did_set" else "did_clear"]], u = u[!lk, , drop = FALSE])
    if (any(lk)) mp_toast(sprintf("Left %s unchanged: unlock before changing marks.",
      mp_n(sum(lk), "locked contig")), type = "warning")
  }
  # Lock and mark buttons say what the next click does (MitoPilot relabels its
  # mark buttons the same way).
  observe({
    mode <- input$mode %||% "Taxonomy"
    lock_btn <- c(Taxonomy = "tax_lock_btn", Assemble = "asmb_lock_btn", Annotate = "annot_lock_btn")
    if (!mode %in% names(lock_btn)) return()
    d <- samples_data(); sel <- selected_rows()
    lk <- d[[paste0(mode_meta(mode)$pfx, "_lock")]][sel]
    updateActionButton(session, lock_btn[[mode]],
      label = if (identical(mp_flag_next(lk, on = 1L, off = 0L), 0L)) "Unlock" else "Lock")
    if (mode != "Annotate") return()
    open <- sel[!(d$annotate_lock[sel] %in% 1)]
    for (col in names(MP_FLAG_WORDS)) {
      updateActionButton(session, paste0(col, "_btn"), label = MP_FLAG_WORDS[[col]][[
        if (identical(mp_flag_next(d[[col]][open]), "yes")) "set" else "clear"]])
    }
  })
  observeEvent(input$reviewed_btn, toolbar_flag("reviewed"))
  # Exclude & Lock: excluded units are done, so the click also locks them
  observeEvent(input$problematic_btn, {
    d <- samples_data(); sel <- selected_rows()
    open <- sel[!(d$annotate_lock[sel] %in% 1)]
    if (!length(open) || !identical(mp_flag_next(d$problematic[open]), "yes")) return(toolbar_flag("problematic"))
    u <- selected_units()
    lk <- d$annotate_lock[sel] %in% 1
    set_unit_flag("UPDATE annotate SET problematic = 'yes', annotate_lock = 1 WHERE ID = ?",
      verb = MP_FLAG_WORDS$problematic[["did_set"]], u = u[!lk, , drop = FALSE])
    if (any(lk)) mp_toast(sprintf("Left %s unchanged: unlock before changing marks.",
      mp_n(sum(lk), "locked contig")), type = "warning")
  })
  observeEvent(input$partial_btn, toolbar_flag("partial"))

  # ---- batch Trim / Restore on ticked units ----
  batch_trim <- reactiveVal(NULL)
  open_batch_trim <- function(restore) {
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    batch_trim(list(u = u, restore = restore))
    showModal(batch_trim_modal(restore))
  }
  observeEvent(input$annot_trim_btn, open_batch_trim(FALSE))
  observeEvent(input$annot_restore_btn, open_batch_trim(TRUE))
  batch_flank <- function() {
    if (!identical(input$batch_trim_mode, "flank")) return(0L)
    suppressWarnings(as.integer(input$batch_trim_flank))[1]
  }
  batch_flank_bad <- function() identical(input$batch_trim_mode, "flank") &&
    (is.na(batch_flank()) || batch_flank() < 1)
  flank_msg <- "Enter a flank of at least 1 bp."
  batch_trim_plan_now <- function() {
    b <- batch_trim(); req(b)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    batch_trim_plan(con, b$u, batch_flank(), b$restore,
      running = !is.null(run_lock_holder(app_project())))
  }
  observe(shinyjs::toggleState("batch_trim_go", condition = !batch_flank_bad()))
  output$batch_trim_summary <- renderUI({
    if (batch_flank_bad()) return(tags$p(style = "color:#c0392b; font-weight:bold;", flank_msg))
    batch_trim_summary(batch_trim_plan_now(), batch_trim()$restore)
  })
  observeEvent(input$batch_trim_go, {
    b <- batch_trim(); req(b)
    if (batch_flank_bad()) return(mp_toast(flank_msg, type = "warning"))
    p <- batch_trim_plan_now()
    ok <- p[is.na(p$reason), , drop = FALSE]
    if (!nrow(ok)) return(mp_toast("Nothing to change: every contig is skipped.", type = "warning"))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    err <- batch_trim_apply(con, ok, b$restore)
    removeModal(); batch_trim(NULL); edit_tick(edit_tick() + 1); refresh(refresh() + 1)
    skipped <- nrow(p) - nrow(ok)
    mp_toast(paste0(mp_n(sum(is.na(err)), "contig"), if (b$restore) " restored" else " trimmed",
      if (skipped) paste0(", ", skipped, " skipped"), "."))
    if (any(!is.na(err))) rp_alert(paste0(mp_n(sum(!is.na(err)), "contig"), " failed: ",
      paste(sprintf("%s p%s s%s (%s)", ok$ID, ok$path, ok$scaffold, err)[!is.na(err)], collapse = "; ")),
      type = "error")
  })

  # ---- Export group actions ----
  observeEvent(input$assign_group_btn, {
    if (!need_selection(length(selected_ids()))) return()
    showModal(group_modal())
  })
  pending_group <- reactiveVal(NULL)
  do_assign_group <- function(u, groups) {
    if (!app_write(assign_export_group(app_project(), u, groups))) return()
    removeModal(); refresh(refresh() + 1)
    mp_toast(if (length(unique(groups)) == 1) {
      paste0(mp_n(nrow(u), "contig"), " assigned to group '", groups[1], "'.")
    } else {
      paste0(mp_n(nrow(u), "contig"), " split into groups ",
        paste0("'", unique(groups), "'", collapse = ", "), ".")
    })
  }
  observeEvent(input$assign_group_go, {
    grp <- trimws(input$group_name %||% "")
    u <- selected_units()
    if (is.null(u) || nrow(u) == 0) { removeModal(); return() }
    if (!nzchar(grp)) { mp_toast("Type a group name first.", type = "warning"); return() }
    bad <- export_group_problem(grp)
    if (!is.null(bad)) { mp_toast(bad, type = "warning"); return() }
    # A GenBank submission records one assembly method, so a group should not
    # mix assemblers. Offer to split "<name>-<assembler>" like MitoPilot's
    # complete/partial split.
    asm <- sample_assembler(app_project(), u$ID)
    if (length(unique(asm)) > 1) {
      pending_group(list(name = grp, units = u, assembler = asm))
      counts <- table(asm)
      showModal(modalDialog(
        title = mp_modal_title("Mixed assemblers"), size = "m", easyClose = FALSE,
        HTML(sprintf(paste0(
          "A GenBank submission records one assembly method, so a group should ",
          "hold samples from one assembler. This selection has %s. Split into %s, ",
          "or keep them as one group?"),
          paste(sprintf("%d from %s", as.integer(counts), names(counts)), collapse = " and "),
          paste0("'", grp, "-", names(counts), "'", collapse = " and "))),
        footer = mp_footer(
          primary = actionButton("group_split", "Split by assembler"), dismiss = NULL,
          extra = tagList(actionButton("group_back", "Cancel"),
            actionButton("group_keep_one", "Keep as one mixed group")))))
      return()
    }
    do_assign_group(u, grp)
  })
  observeEvent(input$group_split, {
    p <- req(pending_group()); pending_group(NULL)
    do_assign_group(p$units, paste0(p$name, "-", p$assembler))
  })
  observeEvent(input$group_keep_one, {
    p <- req(pending_group()); pending_group(NULL)
    do_assign_group(p$units, p$name)
  })
  observeEvent(input$group_back, {
    pending_group(NULL)
    showModal(group_modal())
  })
  observeEvent(input$clear_group_btn, {
    u <- selected_units()
    if (!need_selection(NROW(u))) return()
    if (!app_write(clear_export_group(app_project(), u))) return()
    refresh(refresh() + 1)
    mp_toast(paste0("Group cleared for ", mp_n(nrow(u), "contig"), "."))
  })

  # ---- Update modal: run the workflow with live Nextflow progress ----
  submit_base <- reactiveVal(NULL)   # timestamp shared by one submission's files

  output$nf_cmd_text <- renderText({
    args <- nf_command(app_project(), resume = isTRUE(input$nf_resume %||% TRUE))
    paste(c("nextflow", args), collapse = " ")
  })

  observeEvent(list(input$tax_run_btn, input$asmb_run_btn), ignoreInit = TRUE, {
    # Taxonomy runs locally in the app (lightweight taxdump lookup, no Nextflow).
    if ((input$mode %||% "") == "Taxonomy") { run_taxonomy_update(); return() }
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    n <- nf_pending_count(app_project())
    if (n == 0) {
      rp_alert(paste("No samples or contigs are Ready to run. Samples need a locked",
        "taxonomy, and locked or finished samples are skipped. Set the State of the rows you",
        "want to Ready to run, then press Update."), type = "warning", title = "Nothing to update")
      return()
    }
    # block: a preset naming a database that is not there assembles nothing, and
    # the pipeline reports success while doing it (dangling symlink in, exit 0
    # out), so the run is pure wasted cluster time
    dead <- tryCatch({
      c0 <- app_con(); on.exit(DBI::dbDisconnect(c0), add = TRUE); missing_assembly_dbs(c0)
    }, error = function(e) character(0))
    if (length(dead) > 0) {
      rp_alert(sprintf(paste0(
        "Cannot run: %d reference database file(s) named by the samples queued for assembly ",
        "are missing or empty (%s). Rebuild them in the Reference step, or point those samples ",
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
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    # one timestamp per opening of the modal, so the script the user is looking
    # at is saved under the same name as the log paths written inside it
    if (hl) submit_base(job_basename("pipeline"))
    script <- if (hl) build_submit_script(app_project(), resume = TRUE, sched = sched,
      base = submit_base()) else ""
    showModal(run_modal(n, headless = hl, sched = sched, script = script))
  })

  # Headless run modal: save / submit the cluster job script (instead of running
  # Nextflow on the login node).
  write_submit_script <- function() {
    base <- submit_base() %||% job_basename("pipeline")
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
  job_sent <- reactiveVal(FALSE)
  resubmit_ok <- reactiveVal(FALSE)
  observeEvent(input$submit_job, {
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    if (isTRUE(job_sent()) && !isTRUE(resubmit_ok())) {
      resubmit_ok(TRUE)
      rp_alert(paste("A job for this project was already submitted from this window and may",
        "still be queued or running. Press Submit as job again to submit another."), type = "warning")
      return()
    }
    resubmit_ok(FALSE)
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
        out <- tryCatch(submit_job(sched$cmd, f),
          error = function(e) conditionMessage(e))
        job_sent(TRUE)
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
    td <- project_taxdump_dir(app_project())
    shinyWidgets::sendSweetAlert(
      title = "Assigning taxonomy...",
      text = if (file.exists(taxdump_tarball(td))) {
        "Hold tight, resolving organism names against the local NCBI taxonomy. This can take a moment."
      } else {
        "Hold tight. The first time, this downloads the NCBI taxonomy (about 80 MB) and indexes it, which takes a few minutes."
      },
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE
    )
    later::later(function() {
      shiny::withReactiveDomain(session, {
        res <- tryCatch({
          ensure_taxdump(td)
          c2 <- app_con()
          mt <- resolve_samples_taxonomy(c2, ids, td)
          DBI::dbDisconnect(c2)
          list(ok = TRUE, mt = mt)
        }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
        if (isTRUE(res$ok)) {
          refresh(shiny::isolate(refresh()) + 1)
          open <- any(res$mt %in% c("ambiguous", "unresolved"))
          shinyWidgets::sendSweetAlert(
            title = if (open) "Some samples need a taxon" else "Taxonomy assigned",
            text = tax_count_text(res$mt),
            type = if (open) "warning" else "success")
        } else {
          shinyWidgets::sendSweetAlert(title = "Taxonomy assignment failed",
            text = res$msg, type = "error")
        }
      })
    }, delay = 0.1)
  }

  nf_run_window(input, output, session, pfx = "", gears = "rp_gears", project = app_project,
    kind = "pipeline", workflow = "pipeline",
    launch = function(nflog) {
      p <- nf_spawn("nextflow", c("-log", nflog,
        nf_command(app_project(), resume = isTRUE(input$nf_resume %||% TRUE)), "-ansi-log"),
        app_project(), nf_env())
      nf_record_blastdb_launch(app_project())
      p
    },
    on_done = function(how) refresh(isolate(refresh()) + 1))

  # ---- Recheck: Update on Annotation Review reruns the checks on review copies ----
  rc_submit_base <- reactiveVal(NULL)   # timestamp shared by one submission's files

  output$rc_nf_cmd_text <- renderText({
    input$annot_run_btn
    paste(c("nextflow", nf_recheck_command(app_project())), collapse = " ")
  })

  observeEvent(input$annot_run_btn, ignoreInit = TRUE, {
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    n <- nf_recheck_pending(app_project())
    if (n == 0) {
      w <- nf_recheck_waiting(app_project())
      rp_alert(paste(c("No contigs are Ready to run. A trim or a moved gene boundary queues a",
        "contig, and so does setting its State to Ready to run. Locked contigs are skipped.",
        if (w > 0) paste0("Waiting for the Assemble step's Update: ", mp_n(w, "contig"), ".")),
        collapse = " "), type = "warning", title = "Nothing to update")
      return()
    }
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    if (hl) rc_submit_base(job_basename("recheck"))
    script <- if (hl) build_recheck_submit_script(app_project(), sched = sched,
      base = rc_submit_base()) else ""
    showModal(run_modal(n, headless = hl, sched = sched, script = script, pfx = "rc_",
      title = sprintf("RiboPilot recheck - %d contig(s) to update", n), resume = FALSE))
  })

  write_recheck_script <- function() {
    base <- rc_submit_base() %||% job_basename("recheck")
    typed <- input$rc_submit_script_text %||% ""
    script <- if (nzchar(trimws(typed))) typed else
      build_recheck_submit_script(app_project(), base = base)
    f <- job_path(app_project(), "script", paste0(base, ".sh"))
    writeLines(script, f)
    try(Sys.chmod(f, "0755"), silent = TRUE)
    f
  }
  observeEvent(input$rc_save_script, {
    f <- write_recheck_script()
    rp_alert(paste0("Saved job script to ", f), type = "message")
  })
  rc_job_sent <- reactiveVal(FALSE)
  rc_resubmit_ok <- reactiveVal(FALSE)
  observeEvent(input$rc_submit_job, {
    held <- run_lock_holder(app_project())
    if (!is.null(held)) { rp_alert(run_lock_message(held, app_project()), type = "warning"); return() }
    if (isTRUE(rc_job_sent()) && !isTRUE(rc_resubmit_ok())) {
      rc_resubmit_ok(TRUE)
      rp_alert(paste("A recheck job for this project was already submitted from this window and",
        "may still be queued or running. Press Submit as job again to submit another."),
        type = "warning")
      return()
    }
    rc_resubmit_ok(FALSE)
    sched <- detect_scheduler()
    if (is.null(sched)) {
      rp_alert("No scheduler (sbatch/qsub/bsub) found on PATH.", type = "warning"); return()
    }
    f <- write_recheck_script()
    shinyWidgets::sendSweetAlert(title = "Submitting...",
      text = paste0("Handing the job script to ", sched$cmd,
        ". Hold tight, a busy cluster can take a few seconds to answer."),
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        out <- tryCatch(submit_job(sched$cmd, f), error = function(e) conditionMessage(e))
        rc_job_sent(TRUE)
        rp_alert(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " ")),
          type = "message")
        removeModal()
      })
    }, delay = 0.1)
  })

  # Launched from .recheck so the pipeline's -resume history is never touched.
  # nf_run_window writes the report before on_done, so contigs the run left
  # running are reported as failures before they go back to the queue.
  nf_run_window(input, output, session, pfx = "rc_", gears = "rc_rp_gears",
    project = app_project, kind = "recheck", workflow = "recheck",
    launch = function(nflog) {
      wd <- file.path(app_project(), ".recheck")
      dir.create(wd, recursive = TRUE, showWarnings = FALSE)
      nf_spawn("nextflow", c("-log", nflog, nf_recheck_command(app_project()), "-ansi-log"),
        wd, nf_env())
    },
    on_done = function(how) {
      try(local({
        con <- app_con()
        on.exit(DBI::dbDisconnect(con))
        recheck_release_running(con)
      }), silent = TRUE)
      refresh(isolate(refresh()) + 1)
    })

  # Stop the app when the last browser session closes (survives refresh). After
  # the run windows, so the grace period starts once their runs have stopped.
  register_app_lifecycle(session)

  # ---- Sample metadata: viewer and export fields ----
  specimen_viewer_server("specimen", open = reactive(input$specimen_open),
    on_change = function() refresh(isolate(refresh()) + 1), con = app_con)
  specimen_viewer_server("nmnh_spec", open = reactive(input$nmnh_meta_open),
    on_change = function() refresh(isolate(refresh()) + 1), con = app_con,
    on_back = function() nmnh_show_report(), back_label = "Back to voucher report")

  export_snap <- reactiveVal(NULL)
  # Choose fields inside Export data: keep what is on screen so closing the
  # fields window reopens Export data as it was.
  take_export_snap <- function() {
    export_snap(lapply(stats::setNames(nm = c("export_group_sel", "export_template",
      "fasta_header", "fasta_header_marker", "export_markers", "export_report",
      "export_itsx_check")), function(id) input[[id]]))
  }
  snap_export <- function(source) {
    take_export_snap()
    meta_view_open(source, closed_input = "meta_view_closed")
  }
  observeEvent(input$token_fields_geome, snap_export("GEOME"))
  observeEvent(input$token_fields_gbif, snap_export("GBIF"))
  observeEvent(input$token_fields_ncbi, snap_export("NCBI"))
  observeEvent(input$meta_view_closed, {
    req(export_snap())
    back_to_export()
  })

  # ---- Export ----
  output$export_sel_n <- renderText({
    refresh()
    paste0(mp_n(nrow(export_units_df()), "contig"), " locked & ready")
  })
  # Exportable units (exportable_units()), narrowed to one export group when given.
  export_units_df <- function(group = NULL) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    locked <- exportable_units(con)
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
    all_ids <- exportable_units(con)$ID
    k <- as.integer(table(all_ids)[u$ID])
    ifelse(!is.na(k) & k > 1,
      sprintf("%s_p%s_s%s", u$ID, u$path, u$scaffold), u$ID)
  }
  # Token data the header templates are checked and resolved against.
  export_hdr_data <- reactive({
    refresh()
    meta_view_ver()
    u <- export_units_df(export_group_choice())
    if (nrow(u) == 0) return(NULL)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u$seqid <- app_export_seqids(con, u)
    with_nmnh(tryCatch(export_header_data(con, u), error = function(e) NULL))
  })

  # Marker templates also resolve {element} / {product}, which belong to a
  # feature rather than to the record; a representative pair is enough to check
  # the template.
  with_marker_fields <- function(d) {
    if (is.null(d)) return(NULL)
    el <- app_query("SELECT element FROM annotations LIMIT 1")$element
    el <- if (length(el)) el[1] else "18S"
    d$element <- el
    d$product <- rrna_element_product(el)
    d
  }
  export_hdr_data_marker <- reactive(with_marker_fields(export_hdr_data()))

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
    tmpls <- list_export_templates(con)
    nm <- if (export_last_template(con) %in% tmpls) export_last_template(con) else "default"
    export_template(nm)
    snap <- export_snap()
    export_snap(NULL)
    if (!is.null(snap)) {
      tmpls <- union(tmpls, snap$export_template)
      nm <- snap$export_template %||% nm
      export_template(nm)
    }
    o <- get_export_opts(con, nm)
    sample_cols <- tryCatch(DBI::dbListFields(con, "samples"), error = function(e) character(0))
    ticked <- tryCatch(DBI::dbGetQuery(con, "SELECT key FROM meta_export_fields")$key,
      error = function(e) character(0))
    u <- export_units_df()
    u <- u[!is.na(u$export_group) & nzchar(u$export_group), , drop = FALSE]
    u$seqid <- app_export_seqids(con, u)
    tok <- tryCatch(with_marker_fields(export_header_data(con, u)), error = function(e) NULL)
    if (is.null(tok)) tok <- data.frame()
    export_modal(groups, tmpls, nm, snap$fasta_header %||% nmnh_template_pref(con, o$fasta_header),
      snap$fasta_header_marker %||% o$fasta_header_marker,
      token_groups = export_token_groups(tok, sample_cols, ticked),
      totals = jsonlite::toJSON(as.list(table(tok$export_group)), auto_unbox = TRUE),
      snap = snap, nmnh_cols = nmnh_columns(con),
      nmnh_choices = c("(none: look up from GBIF, NCBI, or GEOME)" = "",
        setdiff(sample_cols, c("ID", "Taxon", "R1", "R2"))))
  }
  observeEvent(input$export_btn, {
    m <- export_modal_current()
    if (!is.null(m)) showModal(m)
  })

  # Picking a saved name refills the boxes; nothing is written until Save.
  observeEvent(input$export_template, ignoreInit = TRUE, {
    nm <- trimws(input$export_template %||% "")
    req(nzchar(nm))
    # the dropdown reports its value when the window opens; only a change reloads
    if (identical(nm, export_template())) return()
    export_template(nm)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    # The template decides NMNH user; the project setting follows it
    if (nm %in% list_export_templates(con)) {
      o <- get_export_opts(con, nm)
      export_set_last_template(con, nm)
      nmnh_set_user_pref(con, nmnh_template_on(o$fasta_header))
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
  output$fasta_header_preview <- renderUI(hdr_preview_ui(hdr_rec(), export_hdr_data()))
  output$fasta_header_marker_preview <- renderUI(hdr_preview_ui(hdr_mrk(), export_hdr_data_marker()))

  # ---- NMNH user ----
  # The switch mirrors the record template: on iff it carries both NMNH tokens.
  nmnh_res <- reactiveVal(NULL)
  with_nmnh <- function(d) {
    if (is.null(d)) return(NULL)
    r <- nmnh_res()
    for (k in c("nmnh_specimen_voucher", "nmnh_voucherURI")) {
      v <- if (is.null(r)) rep(NA_character_, nrow(d)) else r[[k]][match(d$ID, r$ID)]
      d[[k]] <- ifelse(is.na(v), "", v)
    }
    d
  }
  nmnh_on <- reactive(nmnh_template_on(hdr_rec()))
  nmnh_key <- reactiveVal(NULL)
  nmnh_refresh <- function(force = FALSE) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ids <- unique(export_units_df(export_group_choice())$ID)
    key <- paste(c(ids, unlist(nmnh_columns(con))), collapse = "|")
    if (!force && identical(key, nmnh_key())) return(invisible())
    nmnh_res(NULL)
    res <- withProgress(message = "Checking NMNH vouchers", nmnh_resolve(con, ids))
    nmnh_res(res)
    nmnh_key(key)
  }
  nmnh_refresh_ids <- function(ids) {
    r <- nmnh_res()
    ids <- intersect(ids, r$ID)
    if (!length(ids)) return(invisible())
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    r[match(ids, r$ID), ] <- nmnh_resolve(con, ids)
    nmnh_res(r)
    refresh(isolate(refresh()) + 1)
  }
  observeEvent(list(nmnh_on(), input$export_group_sel), {
    if (isTRUE(nmnh_on())) nmnh_refresh()
  })
  observeEvent(hdr_rec(), {
    on <- nmnh_template_on(hdr_rec())
    if (!identical(on, isTRUE(input$nmnh_user))) {
      shinyWidgets::updatePrettyCheckbox(session, "nmnh_user", value = on)
    }
  })
  nmnh_apply <- function() {
    updateTextAreaInput(session, "fasta_header", value = nmnh_template_add(input$fasta_header))
    if (isTRUE(input$export_markers) && !nmnh_template_on(input$fasta_header_marker)) {
      updateTextAreaInput(session, "fasta_header_marker",
        value = nmnh_template_add(input$fasta_header_marker))
    }
  }
  # marker export switched on after NMNH: give the marker header the NMNH fields too
  observeEvent(input$export_markers, ignoreInit = TRUE, {
    if (isTRUE(input$export_markers) && isTRUE(input$nmnh_user) &&
        !nmnh_template_on(input$fasta_header_marker %||% "")) {
      updateTextAreaInput(session, "fasta_header_marker",
        value = nmnh_template_add(input$fasta_header_marker %||% ""))
    }
  })
  observeEvent(input$nmnh_user, ignoreInit = TRUE, {
    on <- isTRUE(input$nmnh_user)
    local({
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      nmnh_set_user_pref(con, on)
    })
    if (on == nmnh_template_on(input$fasta_header)) return()
    if (!on) {
      main <- nmnh_template_remove(input$fasta_header)
      if (identical(main, input$fasta_header)) {
        shinyWidgets::updatePrettyCheckbox(session, "nmnh_user", value = TRUE)
        mp_toast("The NMNH modifiers in the header were edited, so they were left in. Take them out by hand.",
          "warning")
        return()
      }
      updateTextAreaInput(session, "fasta_header", value = main)
      updateTextAreaInput(session, "fasta_header_marker",
        value = nmnh_template_remove(input$fasta_header_marker %||% ""))
      return()
    }
    found <- unique(unlist(lapply(
      c(input$fasta_header, if (isTRUE(input$export_markers)) input$fasta_header_marker),
      nmnh_existing_mods)))
    if (!length(found)) return(nmnh_apply())
    mp_confirm("nmnh_replace", title = "Replace voucher modifiers?",
      text = paste0("The header already has ", paste(found, collapse = ", "),
        ". Replace with NMNH values?"),
      action_label = "Replace")
  })
  observeEvent(input$nmnh_replace, {
    if (isTRUE(input$nmnh_replace)) nmnh_apply()
    else shinyWidgets::updatePrettyCheckbox(session, "nmnh_user", value = FALSE)
  })
  observeEvent(list(input$nmnh_voucher_col, input$nmnh_uri_col), ignoreInit = TRUE, {
    req(!is.null(input$nmnh_voucher_col), !is.null(input$nmnh_uri_col))
    v <- input$nmnh_voucher_col
    u <- input$nmnh_uri_col
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    cur <- nmnh_columns(con)
    if (identical(v, cur$voucher %|NA|% "") && identical(u, cur$uri %|NA|% "")) return()
    nmnh_set_columns(con, v, u)
    if (isTRUE(nmnh_on())) nmnh_refresh()
  })
  observeEvent(input$nmnh_recheck, {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- nmnh_res()$nmnh_voucherURI
    if (DBI::dbExistsTable(con, "nmnh_ark_cache")) {
      for (x in u[!is.na(u)]) {
        DBI::dbExecute(con, "DELETE FROM nmnh_ark_cache WHERE ark = ? OR parent_ark = ?",
          params = list(x, x))
      }
    }
    nmnh_refresh(force = TRUE)
  })

  output$nmnh_status <- renderUI({
    req(nmnh_on())
    r <- nmnh_res()
    if (is.null(r)) return(span(class = "text-muted", "Checking NMNH values..."))
    bad <- sum(!r$ok)
    skipped <- if (any(!is.na(r$nmnh_voucherURI) & !r$checked & !r$not_nmnh)) {
      " GBIF could not be reached, so the specimen check was skipped."
    }
    if (bad == 0) {
      n_not <- sum(r$not_nmnh)
      return(div(class = "mp-nmnh-status mp-fg-success",
        icon("circle-check"), sprintf(" All %s ready%s. ", mp_n(nrow(r), "sample"),
          if (n_not) sprintf(" (%d marked Not NMNH)", n_not) else ""),
        actionLink("nmnh_report", "Review NMNH values"), skipped))
    }
    div(class = "mp-nmnh-status mp-fg-warning",
      icon("triangle-exclamation"),
      sprintf(" %d of %s %s attention. ", bad, mp_n(nrow(r), "sample"), if (bad == 1) "needs" else "need"),
      actionButton("nmnh_report", "Review and fix", class = "btn-sm"), skipped)
  })
  output$nmnh_sources <- renderUI({
    req(nmnh_on())
    r <- nmnh_res()
    if (is.null(r)) return(NULL)
    div(class = "mp-nmnh-sources",
      div(tags$b("Catalog number: "), nmnh_source_text(r$voucher_source)),
      div(tags$b("Specimen link: "), nmnh_source_text(r$uri_source)))
  })

  nmnh_note <- function(src, note) {
    x <- c(if (!is.na(src)) paste("source:", nmnh_source_label(src)), if (nzchar(note)) note)
    paste(x, collapse = "; ")
  }
  nmnh_show_report <- function() {
    take_export_snap()
    r <- nmnh_res()
    nmnh_report_ids(r$ID[!r$ok | r$not_nmnh])
    nmnh_ready_ids(r$ID[r$ok & !r$not_nmnh])
    nmnh_typed(list())
    showModal(modalDialog(
      title = mp_modal_title("NMNH voucher report", close = FALSE), size = "xl",
      p(class = "mp-fg-warning", icon("triangle-exclamation"),
        " NMNH records should not be submitted to GenBank without both a ",
        tags$code("specimen_voucher"), " and a ", tags$code("voucherURI"), "."),
      opts_help("Type a value into a cell to fix it. It is checked when you leave the ",
        "cell and saved with this project's NMNH values; your mapping file is not ",
        "changed. Clear a cell to go back to the automatic value. Pick a field under a ",
        "cell to take that value from another metadata field. Tick Not NMNH for samples ",
        "that are not NMNH specimens: their values are used as found, without checks. ",
        "A missing value leaves its modifier out of the header."),
      div(class = "mp-nmnh-fields", "Fields to pick from: mapping-file columns and the metadata ",
        "fields ticked for export. ", actionLink("nmnh_choose_fields", "Choose fields"),
        tags$br(), "Not NMNH: ", actionLink("nmnh_flag_all", "Select all"), " | ",
        actionLink("nmnh_flag_none", "Clear all")),
      div(class = "mp-nmnh-table-title",
        sprintf("Need attention or marked Not NMNH (%d)", length(nmnh_report_ids()))),
      div(class = "mp-nmnh-scroll", reactable::reactableOutput("nmnh_table")),
      div(class = "mp-nmnh-table-title",
        sprintf("Pass the NMNH check (%d)", length(nmnh_ready_ids()))),
      div(class = "mp-nmnh-scroll", reactable::reactableOutput("nmnh_ready_table")),
      div(class = "mp-nmnh-actions",
        uiOutput("nmnh_copy"),
        actionButton("nmnh_recheck", "Re-check with GBIF", class = "btn-sm",
          title = "Ask GBIF again about this group's specimen links"),
        downloadButton("nmnh_csv", "Download CSV for bulk edit", class = "btn-sm btn-default"),
        fileInput("nmnh_upload", NULL, accept = ".csv", buttonLabel = "Upload edited CSV",
          placeholder = "CSV file")),
      footer = mp_footer(primary = actionButton("nmnh_back", "Back to export"), dismiss = NULL)
    ))
  }
  observeEvent(input$nmnh_report, nmnh_show_report())
  observeEvent(input$nmnh_choose_fields, meta_view_open(closed_input = "nmnh_fields_closed"))
  observeEvent(input$nmnh_fields_closed, {
    nmnh_ver(nmnh_ver() + 1L)
    nmnh_show_report()
  })
  nmnh_flag_many <- function(value) {
    ids <- nmnh_report_ids()
    req(length(ids))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    for (id in ids) nmnh_set_not_nmnh(con, id, value)
    nmnh_typed(list())
    nmnh_refresh_ids(ids)
  }
  observeEvent(input$nmnh_flag_all, nmnh_flag_many(TRUE))
  observeEvent(input$nmnh_flag_none, nmnh_flag_many(FALSE))
  observeEvent(input$nmnh_back, back_to_export())

  nmnh_problems <- reactive({
    r <- req(nmnh_res())
    r[!r$ok, , drop = FALSE]
  })
  # Rows shown in the report: the problems found when it opened, kept after a
  # fix so the entry and its check stay visible. Invalid entries are not saved.
  nmnh_report_ids <- reactiveVal(character(0))
  nmnh_ready_ids <- reactiveVal(character(0))
  nmnh_typed <- reactiveVal(list())
  nmnh_ver <- reactiveVal(0L)
  # One editable table: problem rows, or rows that pass (with an Auto-fixed column)
  nmnh_report_table <- function(ids, ready = FALSE) {
    r <- req(nmnh_res())
    r <- r[match(ids, r$ID), , drop = FALSE]
    r <- r[!is.na(r$ID), , drop = FALSE]
    typed <- nmnh_typed()
    nmnh_ver()
    choices <- local({
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      nmnh_field_choices(con)
    })
    send <- function(input_id, id, rest) {
      sprintf("Shiny.setInputValue('%s', {id: %s, %s}, {priority: 'event'})",
        input_id, jsonlite::toJSON(id, auto_unbox = TRUE), rest)
    }
    pick <- function(field, index) {
      cur <- r[[paste0(field, "_field")]][index] %|NA|% ""
      opts <- c(None = "none", Automatic = "", choices)
      if (nzchar(cur) && !cur %in% opts) opts <- c(opts, stats::setNames(cur, .nmnh_key_label(cur)))
      tags$select(
        class = "form-control input-sm mp-nmnh-pick", `aria-label` = paste(field, "field", r$ID[index]),
        title = "Field the value comes from",
        onchange = send("nmnh_field", r$ID[index], sprintf("field: '%s', key: this.value", field)),
        lapply(seq_along(opts), function(i) {
          tags$option(value = opts[[i]], selected = if (identical(opts[[i]], cur)) NA, names(opts)[i])
        })
      )
    }
    # html = TRUE: reactable drops inline handlers on R tags
    cell <- function(field) function(value, index) {
      bad <- typed[[paste(r$ID[index], field)]]
      src <- r[[paste0(field, "_source")]][index]
      raw <- r[[paste0(field, "_note")]][index] %|NA|% ""
      # short line: verdict and source; the full notes sit behind a help icon
      short <- if (is.na(src)) "" else paste0(" (", nmnh_source_short(src), ")")
      head1 <- strsplit(raw, "; ", fixed = TRUE)[[1]][1] %|NA|% "problem"
      used <- r[[paste0(field, "_from")]][index]
      full <- if (is.null(bad)) nmnh_note(src, raw) else ""
      if (is.null(bad) && !is.na(used)) full <- paste0("Field: ", used, if (nzchar(raw)) paste0("\n", gsub("; ", "\n", raw, fixed = TRUE)))
      if (is.null(bad) && (identical(raw, "missing") || !nzchar(full))) {
        full <- "No value found in the mapping file, GBIF, NCBI, or GEOME."
      }
      if (identical(r[[paste0(field, "_field")]][index], "none")) {
        full <- "Field set to None: no value, so this modifier is left out of the header."
      }
      st <- if (isTRUE(r$not_nmnh[index]) && is.null(bad)) {
        if (is.na(value)) list("muted", "circle-minus", "Missing: modifier left out")
        else list("muted", "circle-minus", paste0("Not checked", short))
      } else if (!is.null(bad)) {
        list("danger", "circle-xmark", bad$note)
      } else if (isTRUE(r[[paste0(field, "_ok")]][index])) {
        list("success", "circle-check", paste0("Valid", short))
      } else {
        list("danger", "circle-xmark", paste0(toupper(substr(head1, 1, 1)), substring(head1, 2), short))
      }
      as.character(tagList(
        tags$input(
          type = "text", class = "form-control input-sm",
          value = if (!is.null(bad)) bad$value else value %|NA|% "",
          `aria-label` = paste(field, r$ID[index]),
          onchange = send("nmnh_edit", r$ID[index], sprintf("field: '%s', value: this.value", field))
        ),
        pick(field, index),
        div(class = if (st[[1]] == "muted") "text-muted" else paste0("mp-fg-", st[[1]]),
            style = "font-size: var(--mp-fs-meta); white-space: normal;",
            icon(st[[2]]), " ", st[[3]],
              if (nzchar(full)) mp_help_tip(gsub("\n", "<br>", htmltools::htmlEscape(full), fixed = TRUE),
                                             label = paste(field, r$ID[index])))
      ))
    }
    d <- data.frame(ID = r$ID, not_nmnh = r$not_nmnh, specimen_voucher = r$nmnh_specimen_voucher,
                    voucherURI = r$nmnh_voucherURI)
    if (ready) d$auto_fixed <- ifelse(r$fixed, "Yes", "")
    reactable::reactable(
      d,
      compact = TRUE, pagination = FALSE,
      language = reactable::reactableLang(
        noData = if (ready) "No samples pass the NMNH check yet." else "No NMNH problems."),
      columns = c(list(
        ID = reactable::colDef(width = 120, html = TRUE, cell = function(value) {
          as.character(tags$a(
            href = "#", title = "View this sample's metadata",
            onclick = sprintf(
              "event.preventDefault(); Shiny.setInputValue('%s', %s, {priority: 'event'})",
              "nmnh_meta_open", jsonlite::toJSON(value, auto_unbox = TRUE)),
            value
          ))
        }),
        not_nmnh = reactable::colDef(
          name = "Not NMNH", width = 110, align = "center", html = TRUE,
          header = rt_header("Not NMNH", "Tick when the sample is not an NMNH specimen. Its values are used as found, without checks."),
          cell = function(value, index) as.character(tags$input(
            type = "checkbox", checked = if (isTRUE(value)) NA, `aria-label` = paste("Not NMNH", r$ID[index]),
            onchange = send("nmnh_flag", r$ID[index], "value: this.checked")))
        ),
        specimen_voucher = reactable::colDef(minWidth = 280, html = TRUE, cell = cell("voucher")),
        voucherURI = reactable::colDef(minWidth = 280, html = TRUE, cell = cell("uri"))
      ), if (ready) list(auto_fixed = reactable::colDef(
        name = "Auto-fixed", width = 100, align = "center",
        header = rt_header("Auto-fixed", "RiboPilot corrected the format of a value, or replaced a tissue link with its specimen. The ? under a cell says what changed.")
      )))
    )
  }
  output$nmnh_table <- reactable::renderReactable(nmnh_report_table(nmnh_report_ids()))
  output$nmnh_ready_table <- reactable::renderReactable(nmnh_report_table(nmnh_ready_ids(), ready = TRUE))
  output$nmnh_copy <- renderUI({
    r <- nmnh_problems()
    d <- data.frame(ID = r$ID, specimen_voucher = r$nmnh_specimen_voucher,
      voucher_note = r$voucher_note, voucherURI = r$nmnh_voucherURI, uri_note = r$uri_note)
    d[is.na(d)] <- ""
    tsv <- paste(c(paste(names(d), collapse = "\t"), do.call(paste, c(d, sep = "\t"))),
      collapse = "\n")
    tags$button(type = "button", class = "btn btn-default btn-sm",
      onclick = paste0(
        "navigator.clipboard.writeText(", jsonlite::toJSON(tsv, auto_unbox = TRUE), ");",
        "var t=this.querySelector('span');",
        "if(t){var o=t.innerText;t.innerText='Copied!';",
        "setTimeout(function(){t.innerText=o;},1200);}"),
      icon("copy"), tags$span(style = "margin-left: 0.3em;", "Copy table"))
  })
  output$nmnh_csv <- downloadHandler(
    filename = function() paste0("nmnh_vouchers_", Sys.Date(), ".csv"),
    content = function(file) {
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      utils::write.csv(nmnh_download_df(con, nmnh_problems()), file, row.names = FALSE, na = "")
    })
  observeEvent(input$nmnh_edit, {
    e <- input$nmnh_edit
    t <- nmnh_typed()
    key <- paste(e$id, e$field)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    note <- nmnh_edit_value(con, e$id, e$field, e$value)
    t[[key]] <- if (!is.null(note)) list(value = e$value, note = note)
    nmnh_typed(t)
    if (is.null(note)) nmnh_refresh_ids(e$id)
  })
  observeEvent(input$nmnh_flag, {
    e <- input$nmnh_flag
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    nmnh_set_not_nmnh(con, e$id, isTRUE(e$value))
    t <- nmnh_typed()
    t[paste(e$id, c("voucher", "uri"))] <- list(NULL)
    nmnh_typed(t)
    nmnh_refresh_ids(e$id)
  })
  observeEvent(input$nmnh_field, {
    e <- input$nmnh_field
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    nmnh_set_field(con, e$id, e$field, e$key)
    t <- nmnh_typed()
    t[[paste(e$id, e$field)]] <- NULL
    nmnh_typed(t)
    nmnh_refresh_ids(e$id)
  })
  observeEvent(input$nmnh_upload, {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    res <- tryCatch(nmnh_upload_csv(con, input$nmnh_upload$datapath), error = function(e) {
      mp_alert(title = "Upload failed", text = conditionMessage(e), type = "error")
      NULL
    })
    if (is.null(res)) return()
    nmnh_refresh_ids(res$ids)
    mp_toast(sprintf("Updated %s from the CSV.", mp_n(length(res$ids), "sample")))
    if (length(res$bad)) {
      mp_alert(title = "Some values were not used", type = "warning", html = TRUE,
        text = tags$ul(style = "text-align: left;", lapply(res$bad, tags$li)))
    }
  })

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
    rc <- validate_fasta_header(input$fasta_header, export_hdr_data(), require_seqid = TRUE)
    mc <- if (isTRUE(input$export_markers)) {
      validate_fasta_header(input$fasta_header_marker, export_hdr_data_marker(), require_seqid = TRUE)
    } else list(ok = TRUE)
    if (!isTRUE(rc$ok) || !isTRUE(mc$ok)) {
      mp_alert(title = "Invalid FASTA header template",
        text = if (!isTRUE(rc$ok)) rc$message else mc$message, type = "error")
      return()
    }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    set_export_opts(con, input$fasta_header, input$fasta_header_marker, nm)
    updateSelectizeInput(session, "export_template", choices = list_export_templates(con),
      selected = nm, options = list(create = TRUE, maxItems = 1))
    export_template(nm)
    export_set_last_template(con, nm)
    mp_toast(sprintf("Saved header template \"%s\".", nm))
  })

  # What Export will write, for the group currently chosen.
  output$export_summary <- renderUI({
    grp <- export_group_choice(); req(grp)
    refresh()
    n <- nrow(export_units_df(grp))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    out <- export_left_out(con, grp)
    path <- file.path(app_project(), "out", "export", grp)
    files <- c(paste0(grp, ".fasta"), paste0(grp, ".tbl"), "samples/",
      paste0(grp, "_summary.csv"), "README.txt")
    if (isTRUE(input$export_markers)) files <- c(files, "markers/")
    if (isTRUE(input$export_report)) files <- c(files, paste0(grp, "_report.html"))
    div(class = "mp-card mp-card-summary",
      div(class = "mp-card-title", sprintf("%s in group \"%s\"", mp_n(n, "contig"), grp)),
      if (nrow(out)) div(class = "mp-fg-warning", icon("triangle-exclamation"), " ",
        sprintf("%s left out: %s.", mp_n(nrow(out), "contig"), mp_id_list(
          sprintf("%s p%s.s%s (%s)", out$ID, out$path, out$scaffold, out$why), max = 10))),
      div("Written to ", tags$code(class = "mp-path", path), " as: ", paste(files, collapse = ", ")),
      if (dir.exists(path)) div(class = "mp-fg-warning", icon("triangle-exclamation"), " ",
        "This folder already exists. Export deletes it and writes it again."))
  })
  # nothing to write: no group chosen, or none of its units can be exported
  observe({
    grp <- export_group_choice(); refresh()
    shinyjs::toggleState("export_go", condition = !is.null(grp) && nrow(export_units_df(grp)) > 0)
  })
  # Write the bundle into the project's out/export folder (no zip download,
  # mirroring MitoPilot), then show where it landed with an "Open folder" button.
  export_done_dir <- reactiveVal(NULL)
  observeEvent(input$export_go, {
    grp <- export_group_choice()
    u <- if (is.null(grp)) NULL else export_units_df(grp)
    if (is.null(u) || nrow(u) == 0) {
      rp_alert("No contig in this group can be exported. The export window lists each contig left out and why.",
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
    r <- if (nmnh_template_on(input$fasta_header)) nmnh_res()
    if (!is.null(r) && any(!r$ok)) {
      export_pending(u)
      take_export_snap()
      removeModal()
      mp_confirm("nmnh_confirm", title = "NMNH values missing or wrong",
        text = sprintf(paste(
          "%d of %s in this group have NMNH problems. NMNH records should not be",
          "submitted to GenBank without both a specimen_voucher and a voucherURI.",
          "Exporting anyway leaves the empty modifiers out of the headers."),
          sum(!r$ok), mp_n(nrow(r), "sample")),
        action_label = "Export anyway", cancel_label = "View report", danger = TRUE)
      return()
    }
    specimen_then_export(u, grp)
  })
  observeEvent(input$nmnh_confirm, ignoreInit = TRUE, {
    if (isTRUE(input$nmnh_confirm)) {
      export_snap(NULL)
      specimen_then_export(export_pending(), export_group_choice())
    } else nmnh_show_report()
  })
  specimen_then_export <- function(u, grp) {
    sp <- tryCatch(specimen_conflict_rows(u, c(input$fasta_header,
        if (isTRUE(input$export_markers)) input$fasta_header_marker)),
      error = function(e) {
        message("Specimen conflict check failed: ", conditionMessage(e))
        NULL
      })
    if (!is.null(sp) && nrow(sp) > 0) {
      export_pending(u)
      take_export_snap()
      removeModal()
      mp_confirm("specimen_confirm", title = "Sample metadata disagrees",
        text = specimen_warning_html(sp), action_label = "Export anyway",
        danger = TRUE, html = TRUE)
      return()
    }
    screen_then_proceed(u, grp)
  }
  screen_then_proceed <- function(u, grp) {
    flagged <- screen_flagged(u)
    export_pending(u)
    if (nrow(flagged)) {
      showModal(export_screen_modal(flagged))
      return()
    }
    proceed_export(u, grp)
  }
  # Conflicts on the specimen items the header templates use; "not checked"
  # items never count and templates without specimen tokens never warn.
  specimen_conflict_rows <- function(u, tmpl) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ids <- unique(u$ID)
    concepts <- specimen_template_concepts(tmpl, specimen_csv_columns(con))
    if (!length(concepts)) return(NULL)
    specimen_export_warnings(specimen_conflicts(con, ids), concepts, ids)
  }
  observeEvent(input$specimen_confirm, ignoreInit = TRUE, {
    if (isTRUE(input$specimen_confirm)) {
      export_snap(NULL)
      screen_then_proceed(export_pending(), export_group_choice())
    } else back_to_export()
  })
  # Units the contamination checks flagged that no curator has signed off on.
  screen_flagged <- function(u) {
    if (NROW(u) == 0) return(data.frame())
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    s <- export_summary_df(con, cbind(u, seqid = u$ID))
    s[has_screen_risk(s$submission_risk) & !(s$reviewed %in% "yes"), , drop = FALSE]
  }
  # export_bundle() rebuilds the group folder from scratch, so anything the
  # curator put there (or an earlier export they have not submitted yet) goes
  # with it. Say so before it happens, not after.
  proceed_export <- function(u, grp) {
    export_pending(u)
    target <- file.path(app_project(), "out", "export", grp)
    if (dir.exists(target)) {
      showModal(export_overwrite_modal(target,
        length(list.files(target, recursive = TRUE))))
      return(invisible())
    }
    run_export(u, grp)
  }
  export_pending <- reactiveVal(NULL)
  back_to_export <- function() {
    m <- export_modal_current()
    if (!is.null(m)) showModal(m)
  }
  observeEvent(input$export_overwrite_back, back_to_export())
  observeEvent(input$export_screen_back, back_to_export())
  observeEvent(input$export_overwrite_go, run_export(export_pending(), export_group_choice()))
  observeEvent(input$export_screen_keep, proceed_export(export_pending(), export_group_choice()))
  observeEvent(input$export_screen_skip, {
    u <- export_pending()
    f <- screen_flagged(u)
    u <- u[!unit_key(u$ID, u$path, u$scaffold) %in% unit_key(f$ID, f$path, f$scaffold), , drop = FALSE]
    if (nrow(u) == 0) {
      back_to_export()
      rp_alert("Every contig in this group was left out, so nothing was exported.", type = "warning")
      return()
    }
    proceed_export(u, export_group_choice())
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
        data.frame(` ` = "No assembled contigs for this sample yet.", check.names = FALSE),
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
        unit = reactable::colDef(name = "Contig", width = 72, show = TRUE, sticky = "left",
          html = TRUE, header = rt_help_header("Contig", SCAFFOLD_COL_HELP$unit)),
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
        screen_species = reactable::colDef(name = "Top BLAST hit", minWidth = 160, show = TRUE,
          html = TRUE, header = rt_help_header("Top BLAST hit", SCAFFOLD_COL_HELP$screen_species),
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
      rp_alert("Tick one or more contigs first.", type = "warning"); return()
    }
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    txt <- unit_fasta_text(con, sel_id(), u)
    if (length(txt) == 0) {
      rp_alert("Those contigs have no stored sequence.", type = "warning"); return()
    }
    copy_to_clipboard(paste(txt, collapse = "\n"))
    rp_toast(sprintf("Copied %d contig(s) as FASTA.", nrow(u)), type = "message")
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
    lbl <- sprintf("Contig p%s.s%s", r$path, r$scaffold)
    if (is.na(r$screen_accession)) {
      # Two different states used to share one message: BLAST switched off, and
      # BLAST switched on but the remote search never came back. The second told
      # the user to enable something already enabled.
      return(tagList(
        tags$b(lbl),
        div(style = "color:#888; margin-top:4px;",
          no_screen_note(
            failed = identical(as.integer(app_query("SELECT assemble_switch FROM assemble WHERE ID = ?",
              list(sel_id()))$assemble_switch[1]), 3L),
            blast_on = blast_enabled(sel_id()),
            running = !is.null(run_lock_holder(app_project())))),
        seed_hit_block(r)))
    }
    tagList(
      tags$b(lbl),
      div(style = "margin-top:4px;",
        tags$span(style = "color:#555;",
          paste0("Best match in ", searched_in(r), ", not restricted to nuclear rDNA:"))),
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
        sprintf("This contig's top hit differs from the rest of this sample (%s at %s).",
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
          "No nuclear rDNA reference matched this contig. See the unrestricted screen above."
        } else {
          "No BLAST result for this contig, so there is nothing to rank here yet."
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
  observe({
    edit_tick(); refresh()
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    try(seqview_send(session, "seqv_asm-canvas", NULL, id, h$path, h$scaffold), silent = TRUE)
  })
  # Triage one scaffold: keep / drop / undecided. `assemblies.ignore` is derived
  # from the verdict by set_unit_triage(), so nothing downstream changes.
  # Dropping is no longer destructive (the annotate row is held, not deleted),
  # so a reviewed unit only warns; an annotation-locked one is still refused,
  # because the Annotate-mode lock guard cannot see this path.
  pending_triage <- reactiveVal(NULL)
  unit_label <- function(u) paste0(u$path, ".", u$scaffold)
  do_scaf_triage <- function(id, units, verdict) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(for (i in seq_len(nrow(units))) {
      set_unit_triage(con, id, units$path[i], units$scaffold[i], verdict,
        reason = NA_character_, source = "user")
    })) return()
    refresh(refresh() + 1)
    word <- if (is.na(verdict)) "undecided" else verdict
    rp_toast(if (nrow(units) == 1) {
      sprintf("Contig %s marked %s.", unit_label(units[1, ]), word)
    } else {
      sprintf("%d contigs marked %s.", nrow(units), word)
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
    if (isTRUE(app_query("SELECT assemble_lock FROM assemble WHERE ID = ?",
      params = list(id))$assemble_lock == 1)) {
      return(need_unlocked(id, how = "Unlock it in the Assemble table first"))
    }

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
          "Contig%s %s %s annotation-locked. Unlock %s in Annotation Review first.",
          if (length(cls$locked) == 1) "" else "s", paste(cls$locked, collapse = ", "),
          if (length(cls$locked) == 1) "is" else "are",
          if (length(cls$locked) == 1) "that contig" else "those contigs"),
          type = "warning")
      }
      if (nrow(cls$units) == 0) return()
      if (length(cls$ask)) {
        pending_triage(list(id = id, units = cls$units))
        showModal(modalDialog(
          title = if (length(cls$ask) == 1) "Drop a reviewed contig?" else "Drop reviewed contigs?",
          size = "s",
          tags$p(sprintf(paste("Contig%s %s already been reviewed. Dropping",
            "%s from annotation and export. Annotation and curation are kept, so this can be undone."),
            if (length(cls$ask) == 1) paste0(" ", cls$ask) else
              paste0("s ", paste(cls$ask, collapse = ", ")),
            if (length(cls$ask) == 1) "has" else "have",
            if (nrow(cls$units) == 1) "excludes it" else
              sprintf("excludes %d contigs", nrow(cls$units)))),
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
  # the unlocked samples selected when the window opened
  triage_ids <- reactiveVal(character(0))
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
    d <- samples_data(); sel <- selected_rows()
    if (!need_selection(length(sel))) return()
    lk <- d$assemble_lock[sel] %in% 1
    if (all(lk)) return(need_unlocked(d$ID[sel]))
    ids <- unique(d$ID[sel][!lk])
    triage_ids(ids)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    cur <- triage_opts_for(con, ids[1])
    triage_preset(cur)
    showModal(triage_modal(length(ids), cur))
    if (any(lk)) mp_toast(sprintf("Left %s unchanged: unlock before triage.",
      mp_n(sum(lk), "locked sample")), type = "warning")
  })
  triage_pv <- reactive({
    ids <- triage_ids()
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
          "Every contig in this selection is hand-decided or annotation-locked, so Apply would skip them all."
        else if (p$n_changed == 0)
          "Apply would change nothing: every contig already has the verdict these rules propose."
        else
          sprintf(paste("Apply would change %d of %d contig(s):",
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
        sprintf("%d contig(s) you decided by hand will be left alone.", p$n_user_skipped)),
      if (p$n_locked_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d annotation-locked contig(s) will be left alone.", p$n_locked_skipped)),
      # under all-off nothing is caused by Apply, so this reports the standing
      # state rather than an outcome
      if (length(p$emptied) > 0) div(class = "mp-fg-danger", style = "font-weight:600; margin-top:6px;",
        sprintf(if (all_off) "%d sample(s) currently have NO contigs kept: %s"
          else "%d sample(s) would end up with NO contigs kept: %s",
          length(p$emptied), paste(utils::head(p$emptied, 12), collapse = ", ")))
      else div(class = "mp-fg-success", style = "margin-top:6px;",
        if (all_off) "No sample currently has nothing kept."
        else "No sample would be left with nothing kept."))
  })
  observeEvent(input$triage_apply, {
    ids <- triage_ids(); req(length(ids))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    act <- triage_actions()
    # persist first: a pipeline re-run must reproduce what the curator just saw
    if (!app_write(write_triage_opts(con, act))) return()
    n <- 0L
    for (id in ids) {
      n <- n + tryCatch(auto_triage_sample(con, id, opts = act,
        overwrite_user = isTRUE(input$triage_overwrite_user)), error = function(e) 0L)
    }
    removeModal()
    refresh(refresh() + 1)
    rp_alert(sprintf("Triage updated %d contig(s) across %d sample(s).", n, length(ids)),
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

  # effective 5'/3' partial state of the unit's features (curator's flag, else automatic)
  unit_partial_rows <- function(d) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    n <- DBI::dbGetQuery(con, "SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length
    e <- partial_ends(d, if (length(n)) n[1] else NA, unit_trimmed(con, u$ID, u$path, u$scaffold))
    d$partial_start <- as.integer(e$p5); d$partial_stop <- as.integer(e$p3)
    d
  }
  partial_state <- function(el) {
    unit_partial_rows(app_query("SELECT * FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = c(uargs(), list(el))))
  }
  active_partial <- reactive({
    edit_tick(); u <- sel_unit(); el <- active_feature(); req(u, el)
    partial_state(el)
  })
  output$feat_edit_active <- renderText({
    el <- active_feature(); if (is.null(el)) return("No active feature")
    p <- active_partial()
    ends <- c(if (isTRUE(p$partial_start[1] == 1)) "5'", if (isTRUE(p$partial_stop[1] == 1)) "3'")
    paste0("Active feature: ", el, if (length(ends)) paste0(" (partial ", paste(ends, collapse = ", "), ")"))
  })
  output$partial_ctrls <- renderUI({
    req(edit_active())
    p <- active_partial()
    btn <- function(id, label, on) actionButton(id, label,
      class = paste("btn-sm", if (on) "btn-warning" else "btn-default"))
    div(tags$label("Partial", class = "control-label"),
      div(style = "display:flex; gap:4px;",
        btn("toggle_partial_start", "5'", isTRUE(p$partial_start[1] == 1)),
        btn("toggle_partial_stop", "3'", isTRUE(p$partial_stop[1] == 1))))
  })

  observe({
    edit_tick(); refresh(); u <- sel_unit(); req(u)
    try(seqview_send(session, "seqv_ann-canvas", "seqv_ann_pick", u$ID, u$path, u$scaffold, active_feature()), silent = TRUE)
  })
  seqview_fit <- function(el) session$sendCustomMessage("rpseq_select", list(id = "seqv_ann-canvas", row = el))
  observeEvent(input$seqv_ann_pick, active_feature(input$seqv_ann_pick$row))

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
    edit_tick(); refresh(); u <- hit_unit(); grid::grid.draw(fixed_gutter(plot_chimera(u$ID, u$path, u$scaffold)))
  }, height = 296)
  output$chimera_result <- renderUI({
    edit_tick(); refresh(); u <- hit_unit()
    ch <- app_query("SELECT label, reason FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$label[1]))
    chimera_result_text(ch$label[1], ch$reason[1])
  })
  output$chimera_db_source <- renderText({
    edit_tick(); refresh(); u <- hit_unit()
    ch <- app_query("SELECT db_source FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$db_source[1]))
    paste("windows searched against:", paste(basename(strsplit(ch$db_source[1], " ", fixed = TRUE)[[1]]), collapse = ", "))
  })
  output$chimera_segments <- reactable::renderReactable({
    edit_tick(); refresh()
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
        species = reactable::colDef(minWidth = 220, na = "", html = TRUE, header = rt_help_header("Best Match", CHIMERA_COL_HELP$species)),
        pident = reactable::colDef(width = 100, format = reactable::colFormat(digits = 1), html = TRUE,
          header = rt_help_header("Identity %", CHIMERA_COL_HELP$pident)),
        n_windows = reactable::colDef(width = 100, html = TRUE, header = rt_help_header("Windows", CHIMERA_COL_HELP$n_windows))))
  })
  output$chimera_junctions <- renderUI({
    edit_tick(); refresh(); u <- hit_unit()
    o <- chimera_original(app_query, u$ID, u$path, u$scaffold)
    req(NROW(o$ch) == 1)
    j <- chimera_junction_table(o$ch$junctions[1], o$ch$channels[1])
    if (!nrow(j)) return(tags$p(class = "mp-help-text", "No possible junction was found."))
    if (!is.null(o$kept)) j$why <- ifelse(j$position <= o$kept[1] | j$position > o$kept[2],
      paste(j$why, "This part is now trimmed away."), j$why)
    reactable::reactable(j, compact = TRUE, pagination = FALSE, wrap = TRUE,
      columns = list(
        position = reactable::colDef(width = 90, html = TRUE, header = rt_help_header("Position", CHIMERA_COL_HELP$position)),
        checks = reactable::colDef(width = 130, html = TRUE, header = rt_help_header("Checks", CHIMERA_COL_HELP$checks)),
        why = reactable::colDef(minWidth = 300, html = TRUE, header = rt_help_header("Why", CHIMERA_COL_HELP$why))))
  })
  output$chimera_windows <- reactable::renderReactable({
    edit_tick(); refresh()
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
    reactable::reactable(w, compact = TRUE, pagination = TRUE, defaultPageSize = 15,
      defaultColDef = reactable::colDef(na = "", minWidth = 70, headerStyle = list(whiteSpace = "nowrap")),
      rowStyle = function(index) list(background = grDevices::adjustcolor(cols[index], alpha.f = CHIMERA_SEG_ALPHA)),
      columns = list(
        start = reactable::colDef(width = 75, html = TRUE, header = rt_help_header("Start", CHIMERA_COL_HELP$start)),
        stop = reactable::colDef(width = 75, html = TRUE, header = rt_help_header("Stop", CHIMERA_COL_HELP$stop)),
        species = reactable::colDef(minWidth = 150, html = TRUE, header = rt_help_header("Best Match", CHIMERA_COL_HELP$win_species)),
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
      pill("Excluded", a$problematic[1], yes_tone = "warning"),
      pill("Partial", a$partial[1], yes_tone = "warning"))
  })
  output$modal_status_toggles <- renderUI({
    a <- modal_flags()
    # mp_flag_next() picks the value the click writes and the label that names
    # it, so the two never disagree. Coloured while the flag is set.
    toggle_btn <- function(col, tone = "success") {
      w <- MP_FLAG_WORDS[[col]]
      set <- identical(mp_flag_next(a[[col]][1]), "yes")
      btn <- actionButton(paste0("mark_", col), w[[if (set) "set" else "clear"]],
        icon = icon(if (set) "check" else "minus"),
        class = paste("btn-sm btn-default", if (!set) paste0("mp-btn-", tone)), title = w[["tip"]])
      if (isTRUE(a$annotate_lock[1] == 1)) tags$span(title = MP_UNIT_LOCKED, shinyjs::disabled(btn)) else btn
    }
    tagList(
      toggle_btn("reviewed"),
      toggle_btn("problematic", tone = "warning"),
      toggle_btn("partial", tone = "warning"))
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
      return(em("Select a feature (no contig sequence available)."))
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
      return(em(paste("No reference alignment for this feature. Needs BLAST hits +",
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
    unit_partial_rows(app_query("SELECT element, type, start, stop, length, strand, partial_start, partial_stop, status, flags
       FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", params = uargs()))
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
        partial_start = reactable::colDef(name = "Partial", width = 80, html = TRUE, align = "center",
          cell = htmlwidgets::JS("
            function(cellInfo) {
              var row = cellInfo.row || {};
              function pill(t) {
                return '<span style=\"background:#FAA34A30;color:#111111;border:1px solid #FAA34A;' +
                  'border-radius:3px;padding:1px 4px;font-size:11px;white-space:nowrap;\">' + t + '</span>';
              }
              var out = [];
              if (row.partial_start == 1) out.push(pill(\"5'\"));
              if (row.partial_stop == 1) out.push(pill(\"3'\"));
              return out.join(' ');
            }")),
        partial_stop = reactable::colDef(show = FALSE),
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
    active_feature(d$element[i]); seqview_fit(d$element[i])
  })

  # ---- curation flags panel ----
  flag_rows <- reactive({
    edit_tick(); u <- sel_unit(); req(u)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND status = 'flag' ORDER BY start",
      params = uargs())
  })
  output$curate_notes <- renderUI({
    u <- sel_unit(); req(u)
    a <- app_query("SELECT n_flags, curate_flags FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(a) == 0 || is.na(a$n_flags[1]) || a$n_flags[1] == 0) {
      return(span(class = "mp-fg-success", HTML("&#10003; No warnings on the contig as a whole.")))
    }
    tagList(
      span(sprintf("%d contig warning(s): ", a$n_flags[1]), class = "mp-fg-warning", style = "font-weight:bold;"),
      span(a$curate_flags[1] %||% "", class = "mp-fg-warning"))
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
    active_feature(d$element[i]); seqview_fit(d$element[i])
  })

  # ---- edit / save / reset state machine (Edit is a checkbox) ----
  restore_edit_backup <- function() {
    b <- edit_backup()
    if (is.null(b) || !nrow(b)) return(TRUE)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    u <- sel_unit()
    app_write(DBI::dbWithTransaction(con, for (i in seq_len(nrow(b))) {
      DBI::dbExecute(con,
        "UPDATE annotations SET start = ?, stop = ?, length = ?, partial_start = ?, partial_stop = ?
         WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
        params = list(b$start[i], b$stop[i], b$stop[i] - b$start[i] + 1, b$partial_start[i], b$partial_stop[i],
          u$ID, u$path, u$scaffold, b$element[i]))
    }))
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
      if (lock_blocked()) {
        shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE); return()
      }
      if (is.null(active_feature())) {
        shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
        rp_toast("Select a feature first.", type = "warning"); return()
      }
      edit_backup(app_query("SELECT element, start, stop, partial_start, partial_stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?",
        params = uargs()))
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
      if (!restore_edit_backup()) return()
      edit_backup(NULL); edit_dirty(FALSE); edit_tick(edit_tick() + 1)
      rp_toast("Boundary edits discarded.", type = "message")
    }
  }, ignoreNULL = TRUE)
  observeEvent(input$feat_save, {
    if (lock_blocked()) return()
    save_edits()
    shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)  # clean exit (not dirty)
  })
  observeEvent(input$feat_reset, {
    if (!restore_edit_backup()) return()
    edit_dirty(FALSE); edit_tick(edit_tick() + 1)
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
    if (edits_block() || lock_blocked()) return()
    s <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(s) == 0 || is.na(s$length[1])) {
      rp_alert("No contig sequence to trim.", type = "warning"); return()
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
    edit_tick(); u <- sel_unit(); req(u)
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
    trim_bounds(sp[["a"]], sp[["b"]], sp[["n"]], if (mode == "flank") input$trim_flank else 0L)
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
    clipped <- c(if (p[1] > sp[["a"]] - x) sprintf("upstream clipped to %d bp (contig start)", sp[["a"]] - p[1]),
      if (p[2] < sp[["b"]] + x) sprintf("downstream clipped to %d bp (contig end)", p[2] - sp[["b"]]))
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
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    w <- trim_window(); from <- w[1]; to <- w[2]
    n <- DBI::dbGetQuery(con, "SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$length[1]
    if (is.na(from) || is.na(to) || from < 1 || from > to || (!is.na(n) && to > n)) {
      rp_alert(sprintf("Invalid trim range. From must be at least 1 and no more than To, and To no more than the contig's current length (%s bp).",
        n %|NA|% "?"), type = "error"); return()
    }
    ok <- tryCatch({
      trim_assembly(con, u$ID, u$path, u$scaffold, from, to)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # re-curate after the shift
      TRUE
    }, error = function(e) { rp_alert(conditionMessage(e), type = "error"); FALSE })
    if (isTRUE(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()  # back to the details modal, now showing trimmed coords
      rp_alert("Contig trimmed; feature coordinates + coverage shifted.", type = "message")
    }
  })

  # Cancel out of the trim modal: return to the annotate details modal.
  observeEvent(input$trim_cancel, { req(sel_unit()); show_annotate_unit() })

  # ---- restore the active unit to its original (pre-trim) length ----
  observeEvent(input$restore_btn, {
    u <- sel_unit(); req(u)
    if (edits_block() || lock_blocked()) return()
    len <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length[1]
    mp_confirm("restore_confirm", title = "Restore the full assembly?",
      text = sprintf("The contig is %s bp now. Restoring undoes the trim and re-runs curation.", len %||% "?"),
      action_label = "Restore", danger = TRUE)
  })
  observeEvent(input$restore_confirm, {
    req(isTRUE(input$restore_confirm))
    u <- sel_unit(); req(u)
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ok <- tryCatch({
      len <- restore_assembly(con, u$ID, u$path, u$scaffold)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # recompute flags for the full unit
      len
    }, error = function(e) { rp_alert(conditionMessage(e), type = "warning"); NULL })
    if (!is.null(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()
      rp_alert(paste0("Contig restored to original length (", ok, " bp)."), type = "message")
    }
  })

  # ---- notes ----
  # the table reloads once the window closes, not on every keystroke
  notes_dirty <- reactiveVal(FALSE)
  observeEvent(input$modal_hidden, {
    if (!notes_dirty()) return()
    notes_dirty(FALSE); refresh(refresh() + 1)
  })
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
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (app_write(DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = ? WHERE ID = ? AND path = ? AND scaffold = ?",
      params = c(list(new), uargs())))) notes_dirty(TRUE)
  })
  observeEvent(input$annotate_notes_clear, {
    u <- sel_unit(); req(u)
    if (lock_blocked()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con, "UPDATE annotate SET annotate_notes = NULL WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs()))) return()
    notes_dirty(TRUE)
    updateTextAreaInput(session, "annotate_notes", value = "")
    rp_toast("Notes cleared.", type = "message")
  })
  observeEvent(input$assemble_user_notes, ignoreInit = TRUE, {
    id <- sel_id(); req(id)
    new <- input$assemble_user_notes %||% ""
    cur <- app_query("SELECT user_notes FROM assemble WHERE ID = ?", params = list(id))$user_notes[1] %|NA|% ""
    if (identical(new, cur) || (!nzchar(new) && nzchar(cur))) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (app_write(DBI::dbExecute(con, "UPDATE assemble SET user_notes = ? WHERE ID = ?", params = list(new, id)))) notes_dirty(TRUE)
  })
  observeEvent(input$assemble_user_notes_clear, {
    id <- sel_id(); req(id)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con, "UPDATE assemble SET user_notes = NULL WHERE ID = ?", params = list(id)))) return()
    notes_dirty(TRUE)
    updateTextAreaInput(session, "assemble_user_notes", value = "")
    rp_toast("Notes cleared.", type = "message")
  })

  # nudge one boundary by the chosen step (1-50), clamp to sequence, persist
  nudge_boundary <- function(end, dir) {
    if (!isTRUE(edit_active())) { rp_toast("Check Edit before adjusting boundaries.", type = "warning"); return() }
    if (lock_blocked()) return()
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
    if (!app_write(DBI::dbExecute(con,
      "UPDATE annotations SET start = ?, stop = ?, length = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      params = list(start, stop, stop - start + 1, u$ID, u$path, u$scaffold, cx$el)))) return()
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
  }
  observeEvent(input$nt_5_out, nudge_boundary("5", "out"))
  observeEvent(input$nt_5_in,  nudge_boundary("5", "in"))
  observeEvent(input$nt_3_out, nudge_boundary("3", "out"))
  observeEvent(input$nt_3_in,  nudge_boundary("3", "in"))
  toggle_partial <- function(col) {
    if (!isTRUE(edit_active()) || lock_blocked()) return()
    el <- active_feature(); req(el)
    to <- if (isTRUE(partial_state(el)[[col]][1] == 1)) 0L else 1L
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con, sprintf(
      "UPDATE annotations SET %s = ?, edited = 1 WHERE ID = ? AND path = ? AND scaffold = ? AND element = ?",
      col), params = c(list(to), uargs(), list(el))))) return()
    edit_dirty(TRUE); edit_tick(edit_tick() + 1)
  }
  observeEvent(input$toggle_partial_start, toggle_partial("partial_start"))
  observeEvent(input$toggle_partial_stop, toggle_partial("partial_stop"))
  # per-unit curator flag toggles + Lock & Close in the annotate modal
  modal_flag_toggle <- function(col) {
    u <- sel_unit(); req(u)
    if (lock_blocked()) return()
    to <- mp_flag_next(modal_flags()[[col]][1])
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con, sprintf("UPDATE annotate SET %s = ? WHERE ID = ? AND path = ? AND scaffold = ?", col),
      params = c(list(to), uargs())))) return()
    refresh(refresh() + 1)
    rp_toast(paste0("1 contig ", MP_FLAG_WORDS[[col]][[if (to == "yes") "did_set" else "did_clear"]], "."),
      type = "message")
  }
  observeEvent(input$mark_reviewed, modal_flag_toggle("reviewed"))
  observeEvent(input$mark_problematic, {
    if (!identical(mp_flag_next(modal_flags()$problematic[1]), "yes")) return(modal_flag_toggle("problematic"))
    u <- sel_unit(); req(u)
    if (lock_blocked() || edits_block()) return()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con,
      "UPDATE annotate SET problematic = 'yes', annotate_lock = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs()))) return()
    refresh(refresh() + 1)
    rp_toast(paste0("1 contig ", MP_FLAG_WORDS$problematic[["did_set"]], "."), type = "message")
    showModal(annotate_modal(u$ID, u$path, u$scaffold))
  })
  observeEvent(input$mark_partial, modal_flag_toggle("partial"))
  # A locked unit is read only in this window (MitoPilot's lock_blocked()).
  lock_blocked <- function() {
    lk <- app_query("SELECT annotate_lock FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())$annotate_lock
    if (!isTRUE(lk[1] == 1)) return(FALSE)
    rp_toast(MP_UNIT_LOCKED, type = "warning")
    TRUE
  }
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
  # `close`: Lock & Close; otherwise Lock, which redraws the window read only.
  lock_unit <- function(close) {
    u <- sel_unit(); req(u)
    if (edits_block()) return()
    if (!need_one_path(as.data.frame(u))) return()
    f <- lock_flagged(as.data.frame(u))
    if (nrow(f)) return(confirm_flagged_lock(f, function() lock_unit_now(close)))
    lock_unit_now(close)
  }
  lock_unit_now <- function(close) {
    u <- sel_unit()
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    if (!app_write(DBI::dbExecute(con,
      "UPDATE annotate SET annotate_lock = 1 WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs()))) return()
    refresh(refresh() + 1)
    rp_toast("1 contig locked - ready to export.", type = "message")
    if (close) removeModal() else showModal(annotate_modal(u$ID, u$path, u$scaffold))
  }
  observeEvent(input$lock_and_close, lock_unit(TRUE))
  observeEvent(input$modal_lock, lock_unit(FALSE))
}
