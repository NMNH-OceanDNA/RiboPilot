# Named parameter-set editing for each pipeline step (preprocess / assemble /
# blast / annotate / curate). Mirrors MitoPilot's opts-modal pattern (a
# selectizeInput with create = TRUE to pick or name a set, an Edit checkbox, tool
# fields with inline + CLI help), adapted to RiboPilot's flat (non-module) app.
# Clicking an opts cell in the samples table opens the relevant modal for the
# selected samples; Update upserts the set and assigns it, re-queuing the step.

# ---- shared modal building blocks ----

# selectize (pick/create a set name) + Edit checkbox
opts_name_row <- function(name_id, edit_id, choices, selected) {
  shiny::div(style = "display:flex; flex-flow:row nowrap; align-items:center; gap:2em;",
    shiny::selectizeInput(name_id, "Parameter set name:", width = "100%",
      choices = choices, selected = selected,
      options = list(create = TRUE, maxItems = 1)),
    shiny::div(class = "form-group shiny-input-container", style = "margin-top:39px;",
      shinyWidgets::prettyCheckbox(edit_id, "Edit", value = FALSE, status = "primary")))
}

`%|NA|%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

opts_modal_shell <- function(title, n, name_row, ..., update_id) {
  shiny::modalDialog(
    title = sprintf("%s for %d sample(s)", title, n), size = "m",
    name_row,
    opts_help("Reusable named set applied to the selected samples; check Edit to ",
      "change values, or type a new name to create a set."),
    ...,
    footer = shiny::tagList(
      shiny::actionButton(update_id, "Update", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

# ---- the five modals ----

pre_opts_modal <- function(cur, choices, n) {
  opts_modal_shell("Preprocessing options", n,
    opts_name_row("po_name", "po_edit", choices, cur$pre_opts),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("po_cpus", "CPUs:",
        value = cur$cpus %|NA|% 4, min = 1, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("po_mem", "Memory (GB):",
        value = cur$memory %|NA|% 8, min = 1, width = "100%")))),
    shinyjs::disabled(shiny::textInput("po_fastp",
      shiny::tagList("fastp options", tool_help_icon("fastp")),
      value = cur$fastp %|NA|% "", width = "100%")),
    opts_help("Command-line flags for fastp (adapter + quality trimming).",
      href = "https://github.com/OpenGene/fastp"),
    update_id = "update_po")
}

assemble_opts_modal <- function(cur, choices, n) {
  opts_modal_shell("Assembly options", n,
    opts_name_row("ao_name", "ao_edit", choices, cur$assemble_opts),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_cpus", "CPUs:",
        value = cur$cpus %|NA|% 4, min = 1, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_mem", "Memory (GB):",
        value = cur$memory %|NA|% 16, min = 1, width = "100%")))),
    shinyjs::disabled(shiny::textInput("ao_getorg",
      shiny::tagList("GetOrganelle options", tool_help_icon("getOrganelle")),
      value = cur$getOrganelle %|NA|% "", width = "100%")),
    opts_help("Command-line flags for GetOrganelle (nrDNA assembly). Seed/label ",
      "databases are set at project creation.",
      href = "https://github.com/Kinggerm/GetOrganelle"),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_minlen",
        "Min assembly length (bp):", value = cur$min_assembly_length %|NA|% 0, min = 0, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_maxpaths",
        "Max paths:", value = cur$max_paths %|NA|% 0, min = 0, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_maxscaf",
        "Max scaffolds:", value = cur$max_scaffolds %|NA|% 0, min = 0, width = "100%")))),
    update_id = "update_ao")
}

blast_opts_modal <- function(cur, choices, n) {
  opts_modal_shell("GenBank BLAST options", n,
    opts_name_row("bo_name", "bo_edit", choices, cur$blast_opts),
    shinyjs::disabled(shinyWidgets::prettyCheckbox("bo_run",
      "Run GenBank BLAST (find the closest reference for review + curation)",
      value = (cur$run_blast %|NA|% 1) == 1, status = "primary")),
    shinyjs::disabled(shiny::textAreaInput("bo_entrez",
      shiny::tagList("Entrez query", tool_help_icon("blastn")),
      value = cur$entrez_query %|NA|% "", width = "100%", height = "90px")),
    opts_help("NCBI Entrez filter restricting the remote BLAST to nuclear ",
      "ribosomal / ITS records."),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:2", shinyjs::disabled(shiny::textInput("bo_extra",
        "Extra blastn options:", value = cur$extra_opts %|NA|% "", width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("bo_mts",
        "Max target seqs:", value = cur$max_target_seqs %|NA|% 5, min = 1, width = "100%")))),
    update_id = "update_bo")
}

annotate_opts_modal <- function(cur, choices, n) {
  tool <- cur$rrna_tool %|NA|% "pybarrnap"
  kingdoms <- if (tool == "barrnap") c("bac", "arc", "fun") else c("bac", "arc", "euk", "all")
  opts_modal_shell("Annotation options", n,
    opts_name_row("no_name", "no_edit", choices, cur$annotate_opts),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("no_cpus", "CPUs:",
        value = cur$cpus %|NA|% 4, min = 1, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("no_mem", "Memory (GB):",
        value = cur$memory %|NA|% 8, min = 1, width = "100%")))),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::selectInput("no_tool",
        shiny::tagList("rRNA anchor tool", tool_help_icon("pybarrnap"), tool_help_icon("barrnap")),
        choices = c("pybarrnap", "barrnap"), selected = tool, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::selectInput("no_kingdom", "Kingdom:",
        choices = kingdoms, selected = cur$rrna_kingdom %|NA|% "euk", width = "100%")))),
    opts_help("pybarrnap (default) is a euk-capable barrnap reimplementation. ",
      "Kingdom values differ by tool (pybarrnap: bac/arc/euk/all; barrnap: bac/arc/fun)."),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::textInput("no_itsxtaxon",
        shiny::tagList("ITSx profile set (-t)", tool_help_icon("ITSx")),
        value = cur$itsx_taxon %|NA|% "all", width = "100%"))),
      shiny::div(style = "flex:2", shinyjs::disabled(shiny::textInput("no_itsxopts",
        "Extra ITSx options:", value = cur$itsx_opts %|NA|% "", width = "100%")))),
    opts_help("ITSx locates the ITS1 / ITS2 spacers.",
      href = "https://microbiology.se/software/itsx/"),
    update_id = "update_no")
}

curate_opts_modal <- function(cur, choices, n) {
  opts_modal_shell("Curation options", n,
    opts_name_row("co_name", "co_edit", choices, cur$curate_opts),
    shinyjs::disabled(shiny::selectInput("co_ruleset", "Clade ruleset:",
      choices = RRNA_RULESETS, selected = cur$ruleset %|NA|% "metazoa_generic", width = "100%")),
    opts_help("Per-clade expected feature lengths + anchor-identity thresholds ",
      "used to flag oddities during curation."),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("co_pid",
        "Min anchor identity (%):", value = cur$min_anchor_pid %|NA|% 85, min = 0, max = 100, width = "100%"))),
      shiny::div(style = "flex:1; margin-top:32px;", shinyjs::disabled(shinyWidgets::prettyCheckbox("co_contam",
        "Flag contaminants", value = (cur$check_contaminant %|NA|% 1) == 1, status = "primary")))),
    update_id = "update_co")
}

# ---- server wiring ----

#' Register the opts-modal observers on the app server.
#'
#' @param input,output,session the moduleless app's shiny objects.
#' @param deps list of shared accessors from app_server: samples_data (reactive
#'   data.frame with an ID column), selected_ids (fn -> character), refresh
#'   (reactiveVal), con (fn -> DBI connection), mode (fn -> "Assemble"/...).
#' @noRd
opts_handlers <- function(input, output, session, deps) {
  # each opts type: fields it edits + where it lives + how it re-queues
  specs <- list(
    po = list(table = "pre_opts",      key = "pre_opts",      link_tbl = "preprocess",
              queue_tbl = "assemble", queue = "assemble_switch", lock = "assemble_lock",
              modal = pre_opts_modal),
    ao = list(table = "assemble_opts", key = "assemble_opts", link_tbl = "assemble",
              queue_tbl = "assemble", queue = "assemble_switch", lock = "assemble_lock",
              modal = assemble_opts_modal),
    bo = list(table = "blast_opts",    key = "blast_opts",    link_tbl = "assemble",
              queue_tbl = "assemble", queue = "assemble_switch", lock = "assemble_lock",
              modal = blast_opts_modal),
    no = list(table = "annotate_opts", key = "annotate_opts", link_tbl = "annotate",
              queue_tbl = "annotate", queue = "annotate_switch", lock = "annotate_lock",
              modal = annotate_opts_modal),
    co = list(table = "curate_opts",   key = "curate_opts",   link_tbl = "annotate",
              queue_tbl = "annotate", queue = "annotate_switch", lock = "annotate_lock",
              modal = curate_opts_modal)
  )
  targets <- shiny::reactiveVal(character(0))   # sample IDs the open modal applies to
  cur_row <- shiny::reactiveVal(NULL)           # base opts row (full, for new-set copy)

  # CLI-help popovers
  for (tl in c("fastp", "getOrganelle", "blastn", "pybarrnap", "barrnap", "ITSx")) {
    register_tool_help(tl, input, session = session)
  }

  # sample IDs among `ids` that are locked for this step
  locked_ids <- function(ids, lock_col, tbl) {
    if (length(ids) == 0) return(character(0))
    con <- deps$con(); on.exit(DBI::dbDisconnect(con))
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbGetQuery(con, sprintf("SELECT ID FROM %s WHERE %s = 1 AND ID IN (%s)",
      tbl, lock_col, ph), params = as.list(ids))$ID
  }

  # click -> which sample IDs (clicked row, or the whole selection if it's in it)
  click_ids <- function(idx) {
    i <- suppressWarnings(as.integer(idx))
    if (is.na(i)) return(character(0))
    d <- deps$samples_data()
    if (i < 1 || i > nrow(d)) return(character(0))
    clicked <- d$ID[i]
    sel <- deps$selected_ids()
    if (length(sel) && clicked %in% sel) sel else clicked
  }

  open_modal <- function(pfx) {
    sp <- specs[[pfx]]
    ids <- targets()
    con <- deps$con(); on.exit(DBI::dbDisconnect(con))
    # current set for the first target, joined to its opts row
    q <- sprintf("SELECT o.* FROM %s x JOIN %s o ON x.%s = o.%s WHERE x.ID = ?",
      sp$link_tbl, sp$table, sp$key, sp$key)
    cur <- tryCatch(DBI::dbGetQuery(con, q, params = list(ids[1])), error = function(e) NULL)
    if (is.null(cur) || nrow(cur) == 0) {
      cur <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE %s = 'default'", sp$table, sp$key))
    }
    cur_row(cur)
    choices <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", sp$key, sp$table, sp$key))[[1]]
    shiny::showModal(sp$modal(as.list(cur), choices, length(ids)))
  }

  # register open + edit-toggle + update observers for each opts type
  field_ids <- list(
    po = c("po_cpus", "po_mem", "po_fastp"),
    ao = c("ao_cpus", "ao_mem", "ao_getorg", "ao_minlen", "ao_maxpaths", "ao_maxscaf"),
    bo = c("bo_run", "bo_entrez", "bo_extra", "bo_mts"),
    no = c("no_cpus", "no_mem", "no_tool", "no_kingdom", "no_itsxtaxon", "no_itsxopts"),
    co = c("co_ruleset", "co_pid", "co_contam")
  )
  for (pfx in names(specs)) local({
    p <- pfx; sp <- specs[[p]]; fids <- field_ids[[p]]
    shiny::observeEvent(input[[paste0("set_", sp$table)]], {
      ids <- click_ids(input[[paste0("set_", sp$table)]])
      if (length(ids) == 0) {
        shiny::showNotification("No sample selected.", type = "warning"); return()
      }
      locked <- locked_ids(ids, sp$lock, sp$queue_tbl)
      if (length(locked) == length(ids)) {
        shiny::showNotification("Selected sample(s) are locked; unlock to change options.",
          type = "warning"); return()
      }
      targets(setdiff(ids, locked))
      open_modal(p)
    })
    # enable fields on Edit; auto-check Edit when a new (unknown) name is typed
    shiny::observeEvent(input[[paste0(p, "_edit")]], ignoreInit = TRUE, {
      for (f in fids) shinyjs::toggleState(f, condition = isTRUE(input[[paste0(p, "_edit")]]))
    })
    shiny::observeEvent(input[[paste0(p, "_name")]], ignoreInit = TRUE, {
      known <- tryCatch({
        con <- deps$con(); on.exit(DBI::dbDisconnect(con))
        DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s", sp$key, sp$table))[[1]]
      }, error = function(e) character(0))
      if (!(input[[paste0(p, "_name")]] %in% known)) {
        shinyWidgets::updatePrettyCheckbox(session, paste0(p, "_edit"), value = TRUE)
        for (f in fids) shinyjs::enable(f)
      }
    })
    shiny::observeEvent(input[[paste0("update_", p)]], {
      nm <- trimws(input[[paste0(p, "_name")]] %||% "")
      chk <- opts_edit_check(sp, nm, targets(), isTRUE(input[[paste0(p, "_edit")]]), deps)
      if (identical(chk$action, "refuse")) {
        shiny::showNotification(chk$msg, type = "error", duration = 10)
        return()
      }
      if (identical(chk$action, "confirm")) {
        shinyWidgets::confirmSweetAlert(session, inputId = paste0(p, "_confirm_edit"),
          title = "Edit beyond selection?", text = chk$msg, type = "warning",
          btn_labels = c("Cancel", "Apply to all"), danger_mode = TRUE)
        return()
      }
      opts_apply(p, sp, input, deps, targets(), cur_row())
    })
    # proceed with a confirmed "edit beyond selection"
    shiny::observeEvent(input[[paste0(p, "_confirm_edit")]], ignoreInit = TRUE, {
      if (isTRUE(input[[paste0(p, "_confirm_edit")]])) {
        opts_apply(p, sp, input, deps, targets(), cur_row())
      }
    })
  })

  # annotate kingdom choices follow the tool selection
  shiny::observeEvent(input$no_tool, ignoreInit = TRUE, {
    ks <- if (input$no_tool == "barrnap") c("bac", "arc", "fun") else c("bac", "arc", "euk", "all")
    sel <- if (input$no_kingdom %in% ks) input$no_kingdom else ks[[if ("euk" %in% ks) which(ks == "euk") else 1]]
    shiny::updateSelectInput(session, "no_kingdom", choices = ks, selected = sel)
  })
}

# Guard editing a shared, named opts set. Editing an EXISTING set rewrites the
# recorded parameters of every sample on that name, not only the selection, so:
#   refuse  - some of those other samples are locked (their stored results were
#             produced under the old parameters; rewriting corrupts provenance);
#   confirm - other unlocked samples share the set (still a surprise);
#   ok      - a new name (new set), or the set is used only by the selection.
opts_edit_check <- function(sp, nm, ids, editing, deps) {
  if (!isTRUE(editing) || !nzchar(nm) || length(ids) == 0) return(list(action = "ok"))
  con <- deps$con(); on.exit(DBI::dbDisconnect(con))
  exists_set <- nrow(DBI::dbGetQuery(con,
    sprintf("SELECT 1 FROM %s WHERE %s = ? LIMIT 1", sp$table, sp$key),
    params = list(nm))) > 0
  if (!exists_set) return(list(action = "ok"))  # new name -> new set, no shared rewrite
  others <- DBI::dbGetQuery(con,
    sprintf("SELECT ID FROM %s WHERE %s = ?", sp$link_tbl, sp$key),
    params = list(nm))$ID
  others <- setdiff(others, ids)
  if (length(others) == 0) return(list(action = "ok"))
  ph <- paste(rep("?", length(others)), collapse = ",")
  locked <- DBI::dbGetQuery(con,
    sprintf("SELECT ID FROM %s WHERE %s = 1 AND ID IN (%s)", sp$queue_tbl, sp$lock, ph),
    params = as.list(others))$ID
  if (length(locked) > 0) {
    return(list(action = "refuse", msg = sprintf(
      "Parameter set '%s' is used by %d locked sample(s), so it cannot be edited (that would rewrite their recorded parameters). Type a new set name to change parameters for the selected samples only.",
      nm, length(locked))))
  }
  list(action = "confirm", msg = sprintf(
    "Parameter set '%s' is also used by %d sample(s) outside your selection. Editing it changes their parameters too and re-queues only the selected ones. Apply to all?",
    nm, length(others)))
}

# build the full opts row (edited fields overriding the base) and upsert it,
# then assign the set to the target samples and re-queue the step.
opts_apply <- function(pfx, sp, input, deps, ids, cur) {
  if (length(ids) == 0) { shiny::removeModal(); return() }
  nm <- trimws(input[[paste0(pfx, "_name")]] %||% "")
  if (!nzchar(nm)) { shiny::showNotification("Give the parameter set a name.", type = "warning"); return() }
  editing <- isTRUE(input[[paste0(pfx, "_edit")]])
  base <- if (!is.null(cur) && nrow(cur)) as.list(cur) else list()

  con <- deps$con(); on.exit(DBI::dbDisconnect(con))
  if (editing) {
    row <- switch(pfx,
      po = data.frame(pre_opts = nm, cpus = as.integer(input$po_cpus),
        memory = as.integer(input$po_mem), fastp = input$po_fastp %||% "",
        stringsAsFactors = FALSE),
      ao = data.frame(assemble_opts = nm, cpus = as.integer(input$ao_cpus),
        memory = as.integer(input$ao_mem), getOrganelle = input$ao_getorg %||% "",
        seeds_db = base$seeds_db %||% NA_character_, labels_db = base$labels_db %||% NA_character_,
        target = base$target %||% NA_character_, max_paths = as.integer(input$ao_maxpaths),
        max_scaffolds = as.integer(input$ao_maxscaf),
        min_assembly_length = as.integer(input$ao_minlen), stringsAsFactors = FALSE),
      bo = data.frame(blast_opts = nm, run_blast = as.integer(isTRUE(input$bo_run)),
        entrez_query = input$bo_entrez %||% "", extra_opts = input$bo_extra %||% "",
        max_target_seqs = as.integer(input$bo_mts), stringsAsFactors = FALSE),
      no = data.frame(annotate_opts = nm, cpus = as.integer(input$no_cpus),
        memory = as.integer(input$no_mem), rrna_tool = input$no_tool,
        rrna_kingdom = input$no_kingdom, itsx_taxon = input$no_itsxtaxon %||% "all",
        itsx_opts = input$no_itsxopts %||% "", stringsAsFactors = FALSE),
      co = data.frame(curate_opts = nm, ruleset = input$co_ruleset,
        min_anchor_pid = as.numeric(input$co_pid),
        check_contaminant = as.integer(isTRUE(input$co_contam)), stringsAsFactors = FALSE)
    )
    dplyr::rows_upsert(dplyr::tbl(con, sp$table), row, by = sp$key, in_place = TRUE, copy = TRUE)
  }

  for (id in ids) {
    DBI::dbExecute(con, sprintf("UPDATE %s SET %s = ? WHERE ID = ?", sp$link_tbl, sp$key),
      params = list(nm, id))
    DBI::dbExecute(con, sprintf("UPDATE %s SET %s = 1 WHERE ID = ?", sp$queue_tbl, sp$queue),
      params = list(id))
  }
  deps$refresh(deps$refresh() + 1)
  shiny::removeModal()
  shiny::showNotification(sprintf("Applied '%s' to %d sample(s); queued to (re)run.", nm, length(ids)),
    type = "message")
}
