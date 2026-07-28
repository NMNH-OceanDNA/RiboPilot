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
    # Never drop unreachable dirs (a failed run's dir may be cleaned up); mark
    # them and show the path so it is still copyable, and only offer Open when the
    # dir exists on this host.
    d$open <- ifelse(d$exists, "open", "")
    d$on_disk <- ifelse(d$exists, "yes", "missing")
    reactable::reactable(d[, c("process", "status", "on_disk", "modified", "open", "workdir")],
      compact = TRUE, resizable = TRUE, defaultPageSize = 25, wrap = FALSE,
      columns = list(
        status = reactable::colDef(html = TRUE, width = 80, cell = rt_dynamicIcon(
          c(success = "fa fa-circle-check", failed = "fa fa-triangle-exclamation"))),
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
      showNotification("Work directory not found on this host.", type = "warning")
      return()
    }
    open_path(d$workdir[i])
  })

  # selection count (feedback before a bulk action) + CSV export of the table
  output$n_selected <- renderText({
    n <- length(selected_ids())
    if (n == 0) "" else paste(n, "selected")
  })
  output$table_csv <- downloadHandler(
    filename = function() paste0(basename(app_project()), "_",
      tolower(input$mode %||% "table"), ".csv"),
    content = function(file) {
      d <- samples_data()
      d <- d[, setdiff(names(d), c("details", "output")), drop = FALSE]
      utils::write.csv(d, file, row.names = FALSE)
    }
  )

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
           s.sci_name, s.taxid, s.rank, s.match_type, s.time_stamp
         FROM samples s")
    } else if (mode == "Assemble") {
      app_query(
        "SELECT a.assemble_lock, a.assemble_switch, s.ID, s.Taxon, s.sci_name, s.taxid,
           pp.pre_opts, a.assemble_opts, a.blast_opts,
           a.paths, a.length,
           -- unit counts come from `assemblies` so they stay live as scaffolds
           -- are kept/dropped. a.blast_* is not selected: nothing writes it (the
           -- per-unit hits live in `assemblies`), so those columns were blank.
           (SELECT COUNT(*) FROM assemblies x WHERE x.ID = s.ID) AS units,
           (SELECT COUNT(*) FROM assemblies x WHERE x.ID = s.ID
              AND COALESCE(x.ignore,0) = 0) AS kept,
           (SELECT p.path FROM path_scores p WHERE p.ID = s.ID AND p.auto_select = 1) AS sel_path,
           (SELECT round(MAX(p.score),3) FROM path_scores p WHERE p.ID = s.ID) AS best_score,
           a.assemble_notes AS notes, a.time_stamp,
           'Details' AS details, 'Output' AS output
         FROM samples s LEFT JOIN assemble a ON s.ID = a.ID
           LEFT JOIN preprocess pp ON s.ID = pp.ID
         WHERE s.taxonomy_lock = 1")
    } else if (mode == "Export") {
      # export shows annotation-locked UNITS (one row per (ID, path, scaffold))
      app_query(
        "SELECT s.ID, an.path, an.scaffold, s.Taxon, s.sci_name, s.taxid, asm.length,
           an.n_features, an.missing, an.reviewed, an.problematic, an.partial,
           an.export_group, an.time_stamp
         FROM annotate an JOIN samples s ON s.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
         WHERE an.annotate_lock = 1")
    } else {
      # annotate shows one row per non-ignored assembly UNIT of an assembly-locked
      # sample (one (ID, path, scaffold) each, annotated independently)
      app_query(
        "SELECT an.annotate_lock, an.annotate_switch, s.ID, an.path, an.scaffold, s.Taxon, s.sci_name, s.taxid,
           asm.length, an.n_features, an.missing,
           asm.blast_species, asm.blast_accession, asm.blast_pident,
           an.n_flags, an.annotate_opts, an.curate_opts, an.reviewed,
           an.problematic, an.partial, an.time_stamp,
           'Details' AS details, 'Output' AS output
         FROM annotate an
           JOIN samples s ON s.ID = an.ID
           JOIN assemble a ON a.ID = an.ID
           JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
         WHERE a.assemble_lock = 1 AND COALESCE(asm.ignore, 0) = 0")
    }
  })

  # ---- reactable (MitoPilot-style: resizable, filterable, sticky, icons) ----
  # Re-render only on structural changes (mode / column groups); the data is
  # isolated so a refresh() pushes a data-only update below instead of a full
  # re-render, which would reset the user's column filters, sort, and page.
  output$tbl <- reactable::renderReactable({
    mode <- input$mode %||% "Assemble"
    d <- isolate(samples_data())
    mm <- mode_meta(mode)
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

    cols <- list(
      ID = reactable::colDef(name = "ID", minWidth = 130, sticky = "left", html = TRUE, cell = rt_link("open_detail")),
      path = reactable::colDef(name = "Path", width = 60, show = multi_unit),
      scaffold = reactable::colDef(name = "Scaf", width = 60, show = multi_unit),
      Taxon = reactable::colDef(name = "Taxon", minWidth = 140, sticky = "left", html = TRUE, cell = rt_longtext()),
      sci_name = reactable::colDef(name = "Assigned Taxon", minWidth = 200, html = TRUE,
        cell = rt_taxid_link()),
      taxid = reactable::colDef(show = FALSE),
      rank = reactable::colDef(name = "Rank", minWidth = 90, html = TRUE, cell = rt_longtext()),
      match_type = reactable::colDef(name = "Confidence", minWidth = 110, html = TRUE,
        cell = rt_confidence_badge()),
      reviewed = reactable::colDef(name = "Reviewed", html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      problematic = reactable::colDef(name = "Problem", width = 90, html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      partial = reactable::colDef(name = "Partial", width = 80, html = TRUE, cell = rt_bool_badge(hide_no = TRUE)),
      n_flags = reactable::colDef(name = "Flags", width = 70, html = TRUE,
        filterable = FALSE, cell = rt_flag_badge()),
      reference = reactable::colDef(name = "Reference", html = TRUE, cell = rt_ncbi_link()),
      blast_species = reactable::colDef(name = "BLAST Hit", minWidth = 160, html = TRUE, cell = rt_longtext()),
      blast_accession = reactable::colDef(name = "Accession", minWidth = 110, html = TRUE, cell = rt_ncbi_link()),
      blast_pident = reactable::colDef(name = "BLAST %ID", width = 90),
      blast_qcovs = reactable::colDef(name = "BLAST Cov", width = 90),
      details = reactable::colDef(name = "", width = 55, align = "center", filterable = FALSE, html = TRUE,
        cell = rt_icon_link(det_input, "fa fa-up-right-from-square", det_tip)),
      output = reactable::colDef(name = "", width = 55, align = "center", filterable = FALSE, html = TRUE,
        cell = rt_icon_link(out_input, "fa fa-folder-open", "Work / output directories")),
      n_features = reactable::colDef(name = "Features", width = 90, html = TRUE),
      missing = reactable::colDef(name = "Missing", minWidth = 120, html = TRUE, cell = rt_longtext()),
      length = reactable::colDef(name = "Length", width = 90),
      paths = reactable::colDef(name = "Paths", width = 80),
      units = reactable::colDef(name = "Units", width = 70),
      kept = reactable::colDef(name = "Kept", width = 70, html = TRUE,
        cell = rt_kept_count()),
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
      notes = reactable::colDef(name = "Notes", minWidth = 160, html = TRUE, cell = rt_longtext()),
      time_stamp = reactable::colDef(name = "Last Updated", width = 150, filterable = FALSE,
        cell = rt_ts_date())
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

    # Explain an empty table: every mode after Taxonomy is gated on the previous
    # stage's lock, and a fresh project has nothing locked, so the default
    # "No rows found" reads like a failed import.
    nodata <- switch(mode,
      Taxonomy = "No samples yet. Create the project with a mapping file (new_project()).",
      Assemble = "No taxonomy-locked samples. Assign and lock taxonomy in Taxonomy mode first.",
      Annotate = "No assembly-locked samples. Review and lock assemblies in Assemble mode first.",
      Export   = "No annotation-locked samples. Review and lock annotations in Annotate mode first.",
      "No rows.")

    sorted <- if ("time_stamp" %in% names(d)) list(time_stamp = "desc") else NULL

    # Column groups (Assemble table): collapse wide sections via the "Columns"
    # picker. Columns not in a named group are Core and always shown.
    col_groups <- list(
      Stats = c("paths", "units", "kept", "sel_path", "best_score"),
      Options = c("pre_opts", "assemble_opts", "blast_opts", "annotate_opts", "curate_opts"))
    keep_cols <- cols[intersect(names(cols), names(d))]
    if (mode == "Assemble") {
      sel_g <- input$col_groups %||% names(col_groups)
      hidden <- unlist(col_groups[setdiff(names(col_groups), sel_g)], use.names = FALSE)
      for (cn in intersect(hidden, names(d))) keep_cols[[cn]] <- reactable::colDef(show = FALSE)
    }

    reactable::reactable(
      d,
      resizable = TRUE, filterable = TRUE, striped = TRUE, compact = TRUE,
      searchable = TRUE, showPageSizeOptions = TRUE, wrap = FALSE, height = "100%",
      selection = "multiple", onClick = "select",
      defaultPageSize = 50, pageSizeOptions = c(25, 50, 100, 200),
      defaultSorted = sorted,
      language = reactable::reactableLang(noData = nodata),
      # consistent header capitalization (Title Case names above) + left alignment
      defaultColDef = reactable::colDef(align = "left", headerVAlign = "center"),
      rowStyle = rt_highlight_row(), rowClass = rowclass,
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
        showNotification("Lock the assembly (Assemble mode) before annotating this sample.",
          type = "warning", duration = 6)
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
  selected_ids <- function() {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel) || length(sel) == 0) {
      return(character(0))
    }
    vis <- visible_mask()
    sel <- sel[sel >= 1 & sel <= length(vis) & vis[sel]]
    samples_data()$ID[sel]
  }
  # Visible selected rows as (ID, path, scaffold) units (path/scaffold = 1 when the
  # active table has no unit columns).
  selected_units <- function() {
    sel <- reactable::getReactableState("tbl", "selected")
    if (is.null(sel) || length(sel) == 0) return(NULL)
    vis <- visible_mask()
    sel <- sel[sel >= 1 & sel <= length(vis) & vis[sel]]
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

  # ---- bulk lock / state / reviewed actions ----
  # In Annotate mode each row is a (path, scaffold) unit, so the bulk update keys
  # on the unit; the annotate `sql` (ending "WHERE ID = ?") is rewritten to add
  # the unit predicate.
  set_flag <- function(sql) {
    if (identical(input$mode %||% "Assemble", "Annotate")) return(set_unit_flag(sql))
    ids <- selected_ids()
    if (length(ids) == 0) {
      showNotification("No samples selected.", type = "warning")
      return()
    }
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
    showNotification(paste("Updated", length(ids), "sample(s)."), type = "message")
  }
  set_unit_flag <- function(sql) {
    u <- selected_units()
    if (is.null(u) || nrow(u) == 0) {
      showNotification("No units selected.", type = "warning"); return()
    }
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
    showNotification(paste("Updated", nrow(u), "unit(s)."), type = "message")
  }
  observeEvent(input$lock_btn, {  # toggle lock on selected
    mm <- mode_meta(input$mode)
    # In Annotate mode locking is per unit, and only ONE assembly path may be
    # locked per sample (paths are alternative resolutions of the same repeat).
    # Guard the bulk toggle against locking a second path.
    if (identical(input$mode %||% "Assemble", "Annotate")) {
      u <- selected_units()
      if (is.null(u) || nrow(u) == 0) { showNotification("No units selected.", type = "warning"); return() }
      con <- app_con(); on.exit(DBI::dbDisconnect(con))
      ids <- unique(u$ID)
      ph <- paste(rep("?", length(ids)), collapse = ",")
      cur <- DBI::dbGetQuery(con, sprintf(
        "SELECT ID, path, scaffold, COALESCE(annotate_lock,0) AS lock FROM annotate WHERE ID IN (%s)", ph),
        params = as.list(ids))
      selk <- paste(u$ID, u$path, u$scaffold)
      curk <- paste(cur$ID, cur$path, cur$scaffold)
      cur$newlock <- ifelse(curk %in% selk, 1L - cur$lock, cur$lock)  # toggle selected units
      viol <- tapply(seq_len(nrow(cur)), cur$ID, function(ix)
        length(unique(cur$path[ix][cur$newlock[ix] == 1L])) > 1)
      if (any(unlist(viol))) {
        showNotification(paste("That would lock more than one assembly path for a sample.",
          "Only one path per sample can be locked; unlock the other path first."),
          type = "error", duration = 8)
        return()
      }
      return(set_unit_flag("UPDATE annotate SET annotate_lock = 1 - COALESCE(annotate_lock, 0) WHERE ID = ?"))
    }
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
      params = list(ids[1]))$s[1]  # annotate has many rows per ID; take one for the default
    showModal(state_modal(mm$pfx, cur))
  })
  observeEvent(input$update_state, {
    mm <- mode_meta(input$mode)
    st <- as.integer(input$new_state %||% "1")
    set_flag(sprintf("UPDATE %s SET %s_switch = %d WHERE ID = ?", mm$tbl, mm$pfx, st))
    removeModal()
  })
  observeEvent(input$reviewed_btn, set_flag("UPDATE annotate SET reviewed = 'yes' WHERE ID = ?"))
  observeEvent(input$problematic_btn, set_flag(
    "UPDATE annotate SET problematic = CASE WHEN problematic = 'yes' THEN 'no' ELSE 'yes' END WHERE ID = ?"))
  observeEvent(input$partial_btn, set_flag(
    "UPDATE annotate SET partial = CASE WHEN partial = 'yes' THEN 'no' ELSE 'yes' END WHERE ID = ?"))

  # ---- Export group actions ----
  observeEvent(input$assign_group_btn, {
    if (length(selected_ids()) == 0) {
      showNotification("No samples selected.", type = "warning"); return()
    }
    showModal(group_modal())
  })
  observeEvent(input$assign_group_go, {
    grp <- trimws(input$group_name %||% "")
    u <- selected_units()
    if (!nzchar(grp) || is.null(u) || nrow(u) == 0) {
      removeModal(); return()
    }
    assign_export_group(app_project(), u, grp)
    removeModal(); refresh(refresh() + 1)
    showNotification(paste0("Assigned ", nrow(u), " unit(s) to group '", grp, "'."), type = "message")
  })
  observeEvent(input$clear_group_btn, {
    u <- selected_units()
    if (is.null(u) || nrow(u) == 0) {
      showNotification("No units selected.", type = "warning"); return()
    }
    clear_export_group(app_project(), u)
    refresh(refresh() + 1)
    showNotification(paste("Cleared group for", nrow(u), "unit(s)."), type = "message")
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
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    script <- if (hl) build_submit_script(app_project(), wf_for_mode(), resume = TRUE, sched = sched) else ""
    showModal(run_modal(input$mode %||% "Assemble", nf_pending_count(app_project(), wf_for_mode()),
      headless = hl, sched = sched, script = script))
  })

  # Headless run modal: save / submit the cluster job script (instead of running
  # Nextflow on the login node).
  write_submit_script <- function() {
    script <- input$submit_script_text %||%
      build_submit_script(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE))
    f <- file.path(app_project(), paste0("submit_", wf_for_mode(), ".sh"))
    writeLines(script, f)
    try(Sys.chmod(f, "0755"), silent = TRUE)
    f
  }
  observeEvent(input$save_script, {
    f <- write_submit_script()
    showNotification(paste0("Saved job script to ", f), type = "message", duration = 10)
  })
  observeEvent(input$submit_job, {
    sched <- detect_scheduler()
    if (is.null(sched)) {
      showNotification("No scheduler (sbatch/qsub/bsub) found on PATH.", type = "warning"); return()
    }
    f <- write_submit_script()
    out <- tryCatch(system2(sched$cmd, shQuote(f), stdout = TRUE, stderr = TRUE),
      error = function(e) conditionMessage(e))
    showNotification(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " ")),
      type = "message", duration = 12)
    removeModal()
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
    nf_env <- c("current", COLUMNS = "500", NXF_ANSI_SUMMARY = "true")
    nxf_pin <- nf_pin_version()
    if (!is.na(nxf_pin)) nf_env <- c(nf_env, NXF_VER = nxf_pin)
    p <- processx::process$new("nextflow",
      args = c(nf_command(app_project(), wf_for_mode(), resume = isTRUE(input$nf_resume %||% TRUE)),
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
    n <- nrow(export_units_df())
    paste(n, if (n == 1) "unit" else "units", "locked & ready")
  })
  # Exportable UNITS: annotation-locked, kept, not flagged problematic;
  # optionally narrowed to one group, else to the current selection. Unit grain
  # throughout, so selecting one scaffold of a sample exports that scaffold.
  export_units_df <- function(group = NULL) {
    locked <- app_query(
      "SELECT an.ID, an.path, an.scaffold, an.export_group
       FROM annotate an
       JOIN assemblies asm ON asm.ID = an.ID AND asm.path = an.path AND asm.scaffold = an.scaffold
       WHERE an.annotate_lock = 1 AND COALESCE(asm.ignore,0) = 0
         AND COALESCE(an.problematic, '') != 'yes'")
    if (nrow(locked) == 0) return(locked)
    key <- function(d) paste(d$ID, as.integer(d$path), as.integer(d$scaffold))
    if (!is.null(group) && nzchar(group)) {
      return(locked[!is.na(locked$export_group) & locked$export_group == group, , drop = FALSE])
    }
    u <- selected_units()
    if (is.null(u) || nrow(u) == 0) return(locked)
    locked[key(locked) %in% key(u), , drop = FALSE]
  }
  # "Export Data" opens a modal to pick a group (or all locked) before export.
  observeEvent(input$export_btn, {
    groups <- app_query(
      "SELECT DISTINCT export_group FROM annotate
       WHERE export_group IS NOT NULL AND export_group != '' ORDER BY export_group")$export_group
    showModal(export_modal(groups, nrow(export_units_df())))
  })
  export_group_choice <- function() {
    g <- input$export_group_sel %||% "__all__"
    if (identical(g, "__all__") || !nzchar(g)) NULL else g
  }
  # Write the bundle into the project's out/export folder (no zip download,
  # mirroring MitoPilot), then show where it landed with an "Open folder" button.
  export_done_dir <- reactiveVal(NULL)
  observeEvent(input$export_go, {
    grp <- export_group_choice()
    u <- export_units_df(grp)
    if (nrow(u) == 0) {
      showNotification("No annotation-locked units to export. Lock annotations first.",
        type = "warning", duration = 6)
      return()
    }
    gdir <- tryCatch(
      export_bundle(app_project(), units = u, group = grp,
        out_dir = file.path(app_project(), "out"), zip = FALSE),
      error = function(e) {
        showNotification(paste("Export failed:", conditionMessage(e)), type = "error", duration = 8)
        NULL
      })
    if (is.null(gdir)) return()
    export_done_dir(gdir)
    removeModal()
    showModal(modalDialog(
      title = "Export complete", size = "m", easyClose = TRUE,
      tags$p("Data exported to:"),
      div(class = "code-block", style = "white-space:normal; word-break:break-all;", gdir),
      footer = tagList(
        actionButton("open_export_dir", "Open folder", icon = icon("folder-open"), class = "btn-sm"),
        modalButton("Close"))
    ))
  })
  observeEvent(input$open_export_dir, open_path(export_done_dir()))

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
    d$hits <- "hits"
    # reactable renders in data-frame column order, so the row identifier has to
    # lead here; colDef(sticky=) alone does not move it.
    lead <- c("unit", "triage", "length", "topology", "depth", "elements")
    d <- d[, c(intersect(lead, names(d)), setdiff(names(d), lead)), drop = FALSE]
    reactable::reactable(d, compact = TRUE, resizable = TRUE, filterable = TRUE,
      searchable = TRUE, sortable = TRUE, wrap = FALSE, defaultPageSize = 20,
      defaultColDef = reactable::colDef(align = "left", show = FALSE),
      columns = list(
        unit = reactable::colDef(name = "Unit", width = 72, show = TRUE, sticky = "left"),
        length = reactable::colDef(name = "Length", width = 80, show = TRUE),
        topology = reactable::colDef(name = "Topology", width = 85, show = TRUE),
        # The number is this scaffold's absolute median depth; the bar beside it
        # is that value as a fraction of the sample's deepest scaffold. Naming
        # the column "relative" described the bar and misdescribed the number.
        depth = reactable::colDef(name = "Median depth", width = 125, show = TRUE, html = TRUE,
          cell = rt_depth_bar("depth_max")),
        elements = reactable::colDef(name = "rRNA features", width = 240, show = TRUE, html = TRUE,
          cell = rt_element_strip()),
        screen_species = reactable::colDef(name = "Top GenBank hit", minWidth = 160, show = TRUE,
          html = TRUE, cell = rt_longtext()),
        shares = reactable::colDef(name = "Shares", width = 105, show = TRUE, html = TRUE,
          cell = rt_shares_badge()),
        odd_one_out = reactable::colDef(name = "", width = 55, show = TRUE, html = TRUE,
          filterable = FALSE, cell = rt_odd_badge()),
        screen_accession = reactable::colDef(name = "Accession", width = 110, show = TRUE,
          html = TRUE, cell = rt_ncbi_link()),
        screen_pident = reactable::colDef(name = "%ID", width = 65, show = TRUE),
        hits = reactable::colDef(name = "", width = 55, show = TRUE, html = TRUE, filterable = FALSE,
          cell = rt_unit_link("scaf_hits")),
        triage = reactable::colDef(name = "Keep", width = 68, show = TRUE, html = TRUE,
          filterable = FALSE, cell = rt_triage_bttn("scaf_triage")),
        triage_reason = reactable::colDef(name = "Why", minWidth = 200, show = TRUE,
          html = TRUE, cell = rt_longtext()),
        uid = reactable::colDef(show = FALSE), depth_max = reactable::colDef(show = FALSE),
        sib_taxon = reactable::colDef(show = FALSE), sib_rank = reactable::colDef(show = FALSE)),
      rowStyle = function(index) if (isTRUE(d$ignore[index] == 1)) list(background = "#f3f3f3", color = "#999") else NULL)
  })

  # ---- per-unit BLAST evidence: the unfiltered screen vs the rDNA panel ----
  sel_hit_unit <- reactiveVal(NULL)
  observeEvent(input$scaf_hits, {
    u <- unit_by_uid(input$scaf_hits); req(!is.null(u))
    sel_hit_unit(list(path = u$path, scaffold = u$scaffold))
  })
  # default to the first kept unit so the panel is never empty on open
  observeEvent(scaffold_rows(), {
    d <- scaffold_rows()
    if (nrow(d) == 0) { sel_hit_unit(NULL); return() }
    cur <- sel_hit_unit()
    if (!is.null(cur) && any(d$path == cur$path & d$scaffold == cur$scaffold)) return()
    k <- which(d$ignore == 0)[1]
    if (is.na(k)) k <- 1L
    sel_hit_unit(list(path = d$path[k], scaffold = d$scaffold[k]))
  })
  output$unit_screen <- renderUI({
    h <- sel_hit_unit(); d <- scaffold_rows()
    if (is.null(h) || nrow(d) == 0) return(NULL)
    r <- d[d$path == h$path & d$scaffold == h$scaffold, , drop = FALSE]
    if (nrow(r) == 0) return(NULL)
    r <- r[1, ]
    lbl <- sprintf("Scaffold p%s.s%s", r$path, r$scaffold)
    if (is.na(r$screen_accession)) {
      return(tagList(
        tags$b(lbl),
        div(style = "color:#888; margin-top:4px;",
          "No GenBank screen for this scaffold yet (run WF1 with BLAST enabled).")))
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
        # yes is the GOOD case here (the hit is a nuclear rDNA record), so no invert
        status_badge("Nuclear rDNA", if (isTRUE(r$is_nuclear_rrna == 1)) "yes" else "no"),
        tags$span(style = "margin-left:8px; color:#555;",
          sprintf("shares %s with %s",
            if (is.na(r$shares)) "no resolved rank" else r$shares,
            r$sci_name %|NA|% "this sample"))),
      if (isTRUE(r$odd_one_out)) div(style = "margin-top:5px; color:#7d4a1e;",
        sprintf("This scaffold's top hit differs from the rest of this sample (%s at %s).",
          r$sib_taxon %|NA|% "the majority", r$sib_rank %|NA|% "")))
  })
  output$unit_candidates <- reactable::renderReactable({
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    d <- app_query(
      "SELECT rank, accession, species, title, pident, qcovs, bitscore, evalue
       FROM blast_ref_candidates WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY rank",
      params = list(id, h$path, h$scaffold))
    if (nrow(d) == 0) {
      # An empty panel is itself a finding: the nuclear-rDNA title filter matched
      # nothing, which is common in sparse invertebrate clades.
      return(reactable::reactable(data.frame(
        ` ` = "No nuclear rDNA reference matched this scaffold. See the unrestricted screen above.",
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
  output$unit_coverage_plot <- renderPlot({
    h <- sel_hit_unit(); id <- sel_id(); req(!is.null(h), id)
    plot_coverage(id, h$path, h$scaffold)
  }, height = 130)
  # Triage one scaffold: keep / drop / undecided. `assemblies.ignore` is derived
  # from the verdict by set_unit_triage(), so nothing downstream changes.
  # Dropping is no longer destructive (the annotate row is held, not deleted),
  # so a reviewed unit only warns; an annotation-locked one is still refused,
  # because the Annotate-mode lock guard cannot see this path.
  pending_triage <- reactiveVal(NULL)
  do_scaf_triage <- function(id, p, s, verdict) {
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    set_unit_triage(con, id, p, s, verdict, reason = NA_character_, source = "user")
    refresh(refresh() + 1)
    showNotification(sprintf("Scaffold %s.%s marked %s.", p, s,
      if (is.na(verdict)) "undecided" else switch(verdict, keep = "keep", drop = "drop", verdict)),
      type = "message")
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
    p <- u$path; s <- u$scaffold
    # clicking the state a row is already in clears it back to undecided
    verdict <- if (identical(u$triage %|NA|% "", want)) NA_character_ else want

    if (identical(verdict, "drop")) {
      st <- app_query(
        "SELECT annotate_lock, reviewed, export_group FROM annotate
         WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, p, s))
      if (nrow(st) && isTRUE(as.integer(st$annotate_lock[1]) == 1L)) {
        showNotification(sprintf(
          "Scaffold %s.%s is annotation-locked. Unlock that unit in Annotate mode first.", p, s),
          type = "error", duration = 8)
        return()
      }
      if (nrow(st) && (identical(st$reviewed[1], "yes") ||
                       (!is.na(st$export_group[1]) && nzchar(st$export_group[1])))) {
        pending_triage(list(id = id, p = p, s = s))
        showModal(modalDialog(
          title = "Drop a reviewed scaffold?", size = "s",
          tags$p(sprintf(paste("Scaffold %s.%s has already been reviewed. Dropping it excludes it",
            "from annotation and export. Its annotation and curation are kept, so this can be undone."), p, s)),
          footer = tagList(
            actionButton("scaf_triage_ok", "Drop it", class = "btn-danger btn-sm"),
            modalButton("Cancel"))))
        return()
      }
    }
    do_scaf_triage(id, p, s, verdict)
  })
  observeEvent(input$scaf_triage_ok, {
    pt <- pending_triage(); req(pt); pending_triage(NULL)
    removeModal()
    do_scaf_triage(pt$id, pt$p, pt$s, "drop")
    showModal(assemble_modal(pt$id))
  })
  # ---- bulk triage over the selected samples ----
  triage_target_ids <- function() {
    ids <- unique(selected_ids())
    if (length(ids)) ids else unique(samples_data()$ID)
  }
  # settings come from the DB so the modal and the WF1 triage process agree; the
  # inputs override them live so the preview reflects what Apply would do.
  triage_actions <- reactive({
    a <- triage_actions_default()
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
    tagList(
      tags$p(sprintf("Would drop %d of %d scaffold(s) across %d sample(s); %d kept, %d left undecided.",
        p$n_drop, p$n_units, p$n_samples, p$n_keep, p$n_undecided)),
      if (no_drop) tags$p(style = "color:#555;",
        "No rule is set to Drop, so Apply will only add and clear reasons."),
      if (p$n_user_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d scaffold(s) you decided by hand will be left alone.", p$n_user_skipped)),
      if (p$n_locked_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d annotation-locked scaffold(s) will be left alone.", p$n_locked_skipped)),
      if (length(p$emptied) > 0) div(style = "color:#842029; font-weight:600; margin-top:6px;",
        sprintf("%d sample(s) would end up with NO scaffolds kept: %s",
          length(p$emptied), paste(utils::head(p$emptied, 12), collapse = ", ")))
      else div(style = "color:#2d6a4f; margin-top:6px;", "No sample would be left with nothing kept."))
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
    showNotification(sprintf("Triage updated %d scaffold(s) across %d sample(s).", n, length(ids)),
      type = "message")
  })

  observeEvent(input$set_selected, {
    id <- sel_id()
    req(id, input$sel_path)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    # Records the reviewer's choice for display. It does NOT gate WF2, which
    # selects units by assemblies.ignore, so nothing is re-queued here: the old
    # sample-wide "UPDATE annotate SET annotate_switch = 1 WHERE ID = ?" also
    # re-queued annotate-locked units the reviewer never touched.
    DBI::dbExecute(con,
      "UPDATE path_scores SET auto_select = CASE WHEN path = ? THEN 1 ELSE 0 END WHERE ID = ?",
      params = list(as.integer(input$sel_path), id))
    refresh(refresh() + 1)
    showNotification(paste0("Marked path ", input$sel_path,
      " as the selected path. To annotate only this path, ignore the other paths' scaffolds below."),
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

  # Annotate modal header badges: reviewed / problem / partial (live via refresh()).
  output$modal_status_badges <- renderUI({
    refresh()
    u <- sel_unit(); req(u)
    a <- app_query("SELECT reviewed, problematic, partial FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    if (nrow(a) == 0) {
      return(NULL)
    }
    tagList(
      status_badge("Reviewed", a$reviewed[1]),
      status_badge("Problem", a$problematic[1], invert = TRUE),
      status_badge("Partial", a$partial[1], invert = TRUE, neutral_no = TRUE))
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
    edit_tick(); u <- sel_unit(); req(u)
    app_query("SELECT element, flags FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? AND status = 'flag' ORDER BY start",
      params = uargs())
  })
  output$curate_notes <- renderUI({
    u <- sel_unit(); req(u)
    a <- app_query("SELECT n_flags, annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
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
    showNotification("Boundary edits saved.", type = "message")
  }
  observeEvent(input$feat_edit, ignoreInit = TRUE, {
    if (isTRUE(input$feat_edit)) {
      if (is.null(active_feature())) {
        shinyWidgets::updatePrettyCheckbox(session, "feat_edit", value = FALSE)
        showNotification("Select a feature first.", type = "warning"); return()
      }
      edit_backup(app_query("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs()))
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

  # helper: (re)open the annotate details modal on the current unit
  show_annotate_unit <- function() {
    u <- sel_unit(); req(u); showModal(annotate_modal(u$ID, u$path, u$scaffold))
  }

  # ---- trim the active unit to annotated / custom boundaries ----
  observeEvent(input$trim_btn, {
    u <- sel_unit(); req(u)
    s <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    if (nrow(s) == 0 || is.na(s$length[1])) {
      showNotification("No unit sequence to trim.", type = "warning"); return()
    }
    span <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    a_from <- if (nrow(span) && !is.na(span$a[1])) span$a[1] else 1L
    a_to <- if (nrow(span) && !is.na(span$b[1])) span$b[1] else s$length[1]
    showModal(trim_modal(s$length[1], a_from, a_to))
  })

  # Interactive trim plot: coverage + features with the kept window shaded. The
  # window is the annotated span (annotated mode) or the From/To inputs (custom).
  trim_window <- reactive({
    u <- sel_unit(); req(u)
    if ((input$trim_mode %||% "annotated") == "annotated") {
      sp <- app_query("SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
      n <- app_query("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())$length
      lo <- if (nrow(sp) && !is.na(sp$a[1])) sp$a[1] else 1L
      hi <- if (nrow(sp) && !is.na(sp$b[1])) sp$b[1] else if (length(n)) n[1] else NA
      c(lo, hi)
    } else {
      c(suppressWarnings(as.integer(input$trim_from)), suppressWarnings(as.integer(input$trim_to)))
    }
  })
  output$trim_plot <- renderPlot({
    edit_tick(); u <- sel_unit(); req(u); w <- trim_window(); plot_trim(u$ID, u$path, u$scaffold, w[1], w[2])
  }, height = 230)
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
    if ((input$trim_mode %||% "annotated") == "annotated") {
      span <- DBI::dbGetQuery(con, "SELECT MIN(start) AS a, MAX(stop) AS b FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
      from <- span$a[1]; to <- span$b[1]
    } else {
      from <- suppressWarnings(as.integer(input$trim_from)); to <- suppressWarnings(as.integer(input$trim_to))
    }
    if (is.na(from) || is.na(to) || from > to) {
      showNotification("Invalid trim range.", type = "error"); return()
    }
    ok <- tryCatch({
      trim_assembly(con, u$ID, u$path, u$scaffold, from, to)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # re-curate after the shift
      TRUE
    }, error = function(e) { showNotification(conditionMessage(e), type = "error"); FALSE })
    if (isTRUE(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()  # back to the details modal, now showing trimmed coords
      showNotification("Unit trimmed; feature coordinates + coverage shifted.", type = "message")
    }
  })

  # Cancel out of the trim modal: return to the annotate details modal.
  observeEvent(input$trim_cancel, { req(sel_unit()); show_annotate_unit() })

  # ---- restore the active unit to its original (pre-trim) length ----
  observeEvent(input$restore_btn, {
    u <- sel_unit(); req(u)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    ok <- tryCatch({
      len <- restore_assembly(con, u$ID, u$path, u$scaffold)
      tryCatch(curate_one(con, u$ID, u$path, u$scaffold), error = function(e) NULL)  # recompute flags for the full unit
      len
    }, error = function(e) { showNotification(conditionMessage(e), type = "warning"); NULL })
    if (!is.null(ok)) {
      edit_tick(edit_tick() + 1); refresh(refresh() + 1)
      show_annotate_unit()
      showNotification(paste0("Unit restored to original length (", ok, " bp)."), type = "message")
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
    showNotification("Notes cleared.", type = "message")
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
    showNotification(paste0("Toggled ", label, "."), type = "message")
  }
  observeEvent(input$mark_reviewed, {
    u <- sel_unit(); req(u)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    DBI::dbExecute(con, "UPDATE annotate SET reviewed = 'yes' WHERE ID = ? AND path = ? AND scaffold = ?", params = uargs())
    refresh(refresh() + 1)
    showNotification("Marked reviewed.", type = "message")
  })
  observeEvent(input$mark_problematic, modal_flag_toggle("problematic", "problematic"))
  observeEvent(input$mark_partial, modal_flag_toggle("partial", "partial"))
  observeEvent(input$lock_and_close, {
    u <- sel_unit(); req(u)
    con <- app_con()
    on.exit(DBI::dbDisconnect(con))
    # One assembly path per sample may be locked (paths are alternative
    # resolutions of the same repeat); scaffolds of that path can all be locked.
    other <- DBI::dbGetQuery(con,
      "SELECT DISTINCT path FROM annotate WHERE ID = ? AND annotate_lock = 1 AND path != ?",
      params = list(u$ID, u$path))
    if (nrow(other) > 0) {
      showNotification(sprintf(
        "Path %s is already locked for %s. Only one assembly path can be locked per sample; unlock it first.",
        other$path[1], u$ID), type = "error", duration = 8)
      return()
    }
    DBI::dbExecute(con,
      "UPDATE annotate SET annotate_lock = 1, reviewed = 'yes' WHERE ID = ? AND path = ? AND scaffold = ?",
      params = uargs())
    refresh(refresh() + 1)
    removeModal()
    showNotification("Locked annotation unit.", type = "message")
  })
}
