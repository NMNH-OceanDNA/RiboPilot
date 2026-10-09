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
    shiny::div(class = "form-group shiny-input-container mp-opts-checkbox",
      mp_checkbox(edit_id, "Edit")))
}

# assembler values written to assemble_opts.assembler; main.nf routes anything
# that is not "getorganelle" to the fallback assembler
ASSEMBLER_CHOICES <- c(
  "GetOrganelle (default)" = "getorganelle",
  "SPAdes, baited to the seed DB" = "bait",
  "SPAdes, de-novo (all reads)" = "denovo")

# which tool-option inputs belong to each assembler. De-novo assembles all
# reads, so it has no baiting step and no bowtie2 flags.
assembler_fields <- function(assembler) {
  all <- c("ao_getorg", "ao_bt2", "ao_spades")
  show <- switch(assembler,
    bait = c("ao_bt2", "ao_spades"),
    denovo = "ao_spades",
    "ao_getorg")
  list(show = show, hide = setdiff(all, show))
}

opts_modal_shell <- function(title, n, name_row, ..., update_id) {
  shiny::modalDialog(
    title = mp_modal_title(paste(title, "for", mp_n(n, "sample"))), size = "m",
    name_row,
    opts_help("Reusable named set applied to the selected samples; check Edit to ",
      "change values, or type a new name to create a set.",
      " Saving queues the selected samples for the next Update."),
    ...,
    footer = mp_footer(primary = shiny::actionButton(update_id, "Save")))
}

#' Where a preset's reference database came from, for the modal's help line
#' @param cur One assemble_opts row (may carry a `clade_name` column).
#' @noRd
ref_db_origin <- function(cur) {
  clade <- cur$clade_name %||% NA_character_
  if (length(clade) != 1 || is.na(clade) || !nzchar(trimws(clade))) {
    return("Set at project creation, or by hand.")
  }
  sprintf("Set by the Reference step for clade %s.", clade)
}

#' "" -> NA, so an emptied path is stored as absent rather than as an empty path
#' @noRd
blank_to_na <- function(x) {
  if (length(x) != 1 || is.na(x) || !nzchar(trimws(x))) NA_character_ else trimws(x)
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
    # assembler choice + only that assembler's tool options, boxed off. Each
    # help line is appended INSIDE its input's container, so hiding the input
    # takes its help with it; a standalone help paragraph would be left behind.
    shiny::div(class = "mp-opts-section",
      shiny::div(style = "font-weight:bold; margin-bottom:8px;", "Assembler"),
      shinyjs::disabled(shiny::selectInput("ao_assembler", label = NULL,
        choices = ASSEMBLER_CHOICES,
        selected = (cur$assembler %||% "getorganelle") %|NA|% "getorganelle",
        width = "100%")),
      opts_help("GetOrganelle is the primary assembler. The SPAdes fallbacks are ",
        "for samples it cannot assemble: 'baited' maps reads to the seed database ",
        "and assembles only what maps, 'de-novo' assembles every read. Both keep ",
        "the contigs that match the label database. The options below belong to ",
        "the selected assembler."),
      shiny::tagAppendChild(
        shinyjs::disabled(shiny::textInput("ao_getorg",
          shiny::tagList("GetOrganelle options", tool_help_icon("getOrganelle")),
          value = cur$getOrganelle %|NA|% "", width = "100%")),
        opts_help("Command-line flags for GetOrganelle (nrDNA assembly). The seed ",
          "and label databases are set below.",
          href = "https://github.com/Kinggerm/GetOrganelle",
          id = "help_ao_getorg", nested = TRUE)),
      shiny::tagAppendChild(
        shinyjs::disabled(shiny::textInput("ao_bt2",
          shiny::tagList("bowtie2 options", tool_help_icon("bowtie2")),
          value = (cur$bowtie2 %||% "--local") %|NA|% "--local", width = "100%")),
        opts_help("Flags for the baiting step that recruits reads against the ",
          "seed database. Output and pairing flags are set by the pipeline.",
          href = "https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml",
          id = "help_ao_bt2", nested = TRUE)),
      shiny::tagAppendChild(
        shinyjs::disabled(shiny::textInput("ao_spades",
          shiny::tagList("SPAdes options", tool_help_icon("spades")),
          value = (cur$spades %||% "--isolate") %|NA|% "--isolate", width = "100%")),
        opts_help("Flags for the SPAdes assembly. If SPAdes rejects the library ",
          "the run is retried once with no extra flags.",
          href = "https://ablab.github.io/spades/",
          id = "help_ao_spades", nested = TRUE))),
    # only reaches samples that supplied their own assembly; the assemblers
    # above are what every other sample uses
    shiny::div(class = "mp-opts-section",
      shiny::div(style = "font-weight:bold; margin-bottom:8px;", "Supplied assemblies"),
      shinyjs::disabled(shiny::checkboxInput("ao_findrrna",
        "Search the supplied assembly for nuclear rDNA contigs",
        value = ((cur$find_rrna %||% 1) %|NA|% 1) == 1)),
      opts_help("Applies only to samples with an assembly in the mapping file. ",
        "The assembly is scanned for rRNA genes and the contigs carrying an 18S ",
        "or 28S anchor go on through the pipeline; the rest are recorded and ",
        "left behind. Turn this off to use the file exactly as supplied, which ",
        "is what a single ITS-only fragment needs, since it carries no anchor ",
        "to find.", id = "help_ao_findrrna")),
    # the Reference tab writes these two paths when it wires a clade database;
    # showing them here is the only place a preset's database is visible
    shiny::div(class = "mp-opts-section",
      shiny::div(style = "font-weight:bold; margin-bottom:8px;", "Reference database"),
      shinyjs::disabled(shiny::textInput("ao_seeds", "Seed database:",
        value = cur$seeds_db %|NA|% "", width = "100%")),
      shinyjs::disabled(shiny::textInput("ao_labels", "Label database:",
        value = cur$labels_db %|NA|% "", width = "100%")),
      opts_help(ref_db_origin(cur), " Editing these points this parameter set at a ",
        "different database; the Reference step overwrites them the next time it ",
        "builds or reuses one for this group.", id = "help_ao_refdb")),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_minlen",
        "Min assembly length (bp):", value = cur$min_assembly_length %|NA|% 0, min = 0, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_maxpaths",
        "Max paths:", value = cur$max_paths %|NA|% 0, min = 0, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("ao_maxscaf",
        "Max contigs:", value = cur$max_scaffolds %|NA|% 0, min = 0, width = "100%")))),
    opts_help("Contigs under the minimum length are dropped from the sample but ",
      "kept in the database, so lowering the limit brings them back without ",
      "re-assembling. Max paths and max contigs (counted within one path) are ",
      "different: exceeding either fails the whole sample. 0 means no limit.",
      id = "help_ao_limits"),
    update_id = "update_ao")
}

# Where a BLAST search goes. Stored as two integers because that is what
# blast.nf reads, but offered as one choice: as two checkboxes the user could
# ask for "NCBI only" AND "fall back to NCBI", which blast.nf silently ignores.
BLAST_SOURCE_CHOICES <- c(
  "Bundled database only" = "bundled",
  "Bundled database, then NCBI for anything it cannot name" = "escalate",
  "NCBI only (core_nt)" = "remote")

#' The stored columns as one choice
#'
#' remote wins over fallback for a row carrying both, because that is what the
#' pipeline does with it.
#' @noRd
blast_source <- function(remote_blast, remote_fallback) {
  r <- remote_blast %|NA|% 0
  f <- remote_fallback %|NA|% 1
  if (r == 1) "remote" else if (f == 1) "escalate" else "bundled"
}

#' The choice as the stored columns
#' @noRd
blast_source_cols <- function(source) {
  switch(source %||% "escalate",
    bundled  = list(remote_blast = 0L, remote_fallback = 0L),
    remote   = list(remote_blast = 1L, remote_fallback = 0L),
    list(remote_blast = 0L, remote_fallback = 1L))
}

#' The bundled-database note, rendered as a callout rather than a field
#'
#' A full border plus a radius is exactly what Bootstrap gives a form control,
#' so in a column of disabled grey inputs the old box read as one more input the
#' user could not type in. A tinted fill, a left accent bar and a lead-in INSIDE
#' the block (every input in this app labels from outside and above) each say
#' "this is text", and none of them relies on the user noticing a colour.
#' @param note list(known, text) from panel_note().
#' @noRd
panel_note_box <- function(note) {
  known <- isTRUE(note$known)
  accent <- if (known) "#5b7fa6" else "#9a7b2e"
  shiny::tags$div(role = "note",
    style = paste0("background:#f5f7fa; border-left:3px solid ", accent, ";",
      "border-radius:0 4px 4px 0; padding:9px 13px; margin:6px 0 14px;",
      "font-size:0.9em; color:#333;"),
    shiny::tags$strong(if (known) "Bundled database" else
      "Bundled database: Not yet recorded"),
    shiny::tags$br(),
    note$text %||% "")
}

#' Attach the reference-database note to a blast_opts row
#'
#' Takes the project explicitly: the app is launched from wherever the user's R
#' session happens to sit, so resolving the config relative to the working
#' directory silently found nothing and left the modal showing its placeholder.
#' @noRd
with_panel_note <- function(cur, con, project = app_project()) {
  if (NROW(cur) < 1) return(cur)
  note <- tryCatch(project_panel_note(con, nf_blastdb(project), project = project),
    error = function(e) panel_note(NULL, ""))
  cur$panel_note <- note$text
  cur$panel_known <- isTRUE(note$known)
  cur
}

blast_opts_modal <- function(cur, choices, n) {
  source <- blast_source(cur$remote_blast, cur$remote_fallback)
  # the Entrez filter is passed to every search that leaves the machine, and an
  # escalated scaffold is one of those
  reaches_ncbi <- source != "bundled"
  opts_modal_shell("BLAST options", n,
    opts_name_row("bo_name", "bo_edit", choices, cur$blast_opts),
    shinyjs::disabled(shinyWidgets::prettyCheckbox("bo_run",
      "Run BLAST (find the closest reference for review + curation)",
      value = (cur$run_blast %|NA|% 1) == 1, status = "primary")),
    panel_note_box(list(known = cur$panel_known %|NA|% FALSE,
      text = cur$panel_note %|NA|% "")),
    shinyjs::disabled(shiny::radioButtons("bo_source", "Search:",
      choices = BLAST_SOURCE_CHOICES, selected = source, width = "100%")),
    opts_help("The bundled database answers in milliseconds; NCBI is far slower, ",
      "rate limited, and needs an internet connection from wherever the pipeline ",
      "runs. Escalating catches contaminants and off-target contigs, the things ",
      "a metazoan rDNA database has no record of."),
    shiny::div(id = "bo_entrez_group",
      style = if (reaches_ncbi) NULL else "display:none;",
      shinyjs::disabled(shiny::textAreaInput("bo_entrez",
        shiny::tagList("Entrez query", tool_help_icon("blastn")),
        value = cur$entrez_query %|NA|% "", width = "100%", height = "90px")),
      opts_help("NCBI filter, applied to any search that goes to NCBI, including ",
        "an escalated contig. The bundled database cannot use it.")),
    shiny::div(style = "display:flex; gap:2em;",
      shiny::div(style = "flex:2", shinyjs::disabled(shiny::textInput("bo_extra",
        "Extra blastn options:", value = cur$extra_opts %|NA|% "", width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::numericInput("bo_mts",
        "Max target seqs:", value = cur$max_target_seqs %|NA|% 5, min = 1, width = "100%")))),
    opts_help("Both apply to whichever search runs."),
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
        shiny::tagList("rRNA anchor tool",
          shiny::span(id = "no_help_pybarrnap", tool_help_icon("pybarrnap")),
          shiny::span(id = "no_help_barrnap", tool_help_icon("barrnap"))),
        choices = c("pybarrnap", "barrnap"), selected = tool, width = "100%"))),
      shiny::div(style = "flex:1", shinyjs::disabled(shiny::selectInput("no_kingdom", "Kingdom:",
        choices = kingdoms, selected = cur$rrna_kingdom %|NA|% "euk", width = "100%")))),
    opts_help("pybarrnap (default) is a euk-capable barrnap reimplementation. ",
      "Kingdom values differ by tool (pybarrnap: bac/arc/euk/all; barrnap: bac/arc/fun)."),
    shinyjs::disabled(shiny::textInput("no_pybarrnap", "Extra pybarrnap options:",
      value = cur$pybarrnap_opts %|NA|% "", width = "100%")),
    shinyjs::disabled(shiny::textInput("no_barrnap", "Extra barrnap options:",
      value = cur$barrnap_opts %|NA|% "", width = "100%")),
    opts_help("Flags for the selected tool only; kingdom and threads are set by the pipeline. ",
      "The two tools take different flags (see the ? icon next to the tool name)."),
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
    opts_help("Per-clade expected feature lengths and anchor-identity thresholds used to flag oddities during curation."),
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
  for (tl in c("fastp", "getOrganelle", "bowtie2", "spades", "blastn", "pybarrnap",
               "barrnap", "ITSx")) {
    register_tool_help(tl, input, session = session)
  }

  # show only the selected assembler's tool options
  toggle_assembler_fields <- function(assembler) {
    f <- assembler_fields(assembler)
    for (id in f$show) shinyjs::show(id = id)
    for (id in f$hide) shinyjs::hide(id = id)
  }
  shiny::observeEvent(input$ao_assembler, ignoreInit = TRUE, {
    toggle_assembler_fields(input$ao_assembler)
  })

  # show only the selected rRNA caller's flag field
  toggle_rrna_fields <- function(tool) {
    for (id in c("no_pybarrnap", "no_help_pybarrnap")) shinyjs::toggle(id = id, condition = !identical(tool, "barrnap"))
    for (id in c("no_barrnap", "no_help_barrnap")) shinyjs::toggle(id = id, condition = identical(tool, "barrnap"))
  }
  shiny::observeEvent(input$no_tool, ignoreInit = TRUE, {
    toggle_rrna_fields(input$no_tool)
  })

  # sample IDs among `ids` that are locked for this step
  locked_ids <- function(ids, lock_col, tbl) {
    if (length(ids) == 0) return(character(0))
    con <- deps$con(); on.exit(DBI::dbDisconnect(con))
    ph <- paste(rep("?", length(ids)), collapse = ",")
    DBI::dbGetQuery(con, sprintf("SELECT DISTINCT ID FROM %s WHERE %s = 1 AND ID IN (%s)",
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

  open_modal <- function(pfx, set = NULL) {
    sp <- specs[[pfx]]
    ids <- targets()
    con <- deps$con(); on.exit(DBI::dbDisconnect(con))
    # the named set, else the first target's current set (for annotate, its
    # first unit, as the table shows)
    q <- if (is.null(set)) {
      sprintf("SELECT o.* FROM %s x JOIN %s o ON x.%s = o.%s WHERE x.ID = ?%s LIMIT 1",
        sp$link_tbl, sp$table, sp$key, sp$key,
        if (sp$link_tbl == "annotate") " ORDER BY x.path, x.scaffold" else "")
    } else {
      sprintf("SELECT * FROM %s WHERE %s = ?", sp$table, sp$key)
    }
    cur <- tryCatch(DBI::dbGetQuery(con, q, params = list(set %||% ids[1])),
      error = function(e) NULL)
    if (is.null(cur) || nrow(cur) == 0) {
      cur <- DBI::dbGetQuery(con, sprintf("SELECT * FROM %s WHERE %s = 'default'", sp$table, sp$key))
    }
    # the Reference tab names the preset it wires, so the modal can say which
    # clade a database came from
    if (pfx == "ao" && nrow(cur) > 0) {
      cur$clade_name <- tryCatch(DBI::dbGetQuery(con,
        "SELECT clade_name FROM reference WHERE preset_name = ? LIMIT 1",
        params = list(cur[[sp$key]][1]))$clade_name[1], error = function(e) NA_character_)
    }
    if (pfx == "bo") cur <- with_panel_note(cur, con)
    cur_row(cur)
    choices <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s ORDER BY %s", sp$key, sp$table, sp$key))[[1]]
    shiny::showModal(sp$modal(as.list(cur), choices, length(ids)))
    if (pfx == "ao") {
      toggle_assembler_fields((cur$assembler[1] %||% "getorganelle") %|NA|% "getorganelle")
    }
    if (pfx == "no") {
      toggle_rrna_fields((cur$rrna_tool[1] %||% "pybarrnap") %|NA|% "pybarrnap")
    }
  }

  # register open + edit-toggle + update observers for each opts type
  field_ids <- list(
    po = c("po_cpus", "po_mem", "po_fastp"),
    ao = c("ao_cpus", "ao_mem", "ao_getorg", "ao_minlen", "ao_maxpaths", "ao_maxscaf",
           "ao_assembler", "ao_bt2", "ao_spades", "ao_seeds", "ao_labels",
           "ao_findrrna"),
    bo = c("bo_run", "bo_source", "bo_entrez", "bo_extra", "bo_mts"),
    no = c("no_cpus", "no_mem", "no_tool", "no_kingdom", "no_itsxtaxon", "no_itsxopts",
           "no_barrnap", "no_pybarrnap"),
    co = c("co_ruleset")
  )
  for (pfx in names(specs)) local({
    p <- pfx; sp <- specs[[p]]; fids <- field_ids[[p]]
    shiny::observeEvent(input[[paste0("set_", sp$table)]], {
      ids <- click_ids(input[[paste0("set_", sp$table)]])
      if (!need_selection(length(ids))) return()
      locked <- locked_ids(ids, sp$lock, sp$queue_tbl)
      what <- if (identical(sp$lock, "annotate_lock")) "annotation" else "assembly"
      if (length(locked) == length(ids)) {
        need_unlocked(locked, "sample", how = if (what == "annotation") {
          paste("Their annotation is locked. Unlock their contigs in Annotation Review first;",
            "a sample is listed there only while its assembly is locked")
        } else "Their assembly is locked. Select them on this tab and press Unlock first")
        return()
      }
      targets(setdiff(ids, locked))
      open_modal(p)
      if (length(locked)) mp_toast(sprintf("Skipping %s whose %s is locked: %s.",
        mp_n(length(locked), "sample"), what, mp_id_list(locked)), type = "warning")
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
      nm <- input[[paste0(p, "_name")]]
      if (!(nm %in% known)) {
        shinyWidgets::updatePrettyCheckbox(session, paste0(p, "_edit"), value = TRUE)
        for (f in fids) shinyjs::enable(f)
      } else if (!identical(nm, cur_row()[[sp$key]][1]) ||
                 isTRUE(input[[paste0(p, "_edit")]])) {
        # a known set shows its own values, read-only
        open_modal(p, set = nm)
      }
    })
    shiny::observeEvent(input[[paste0("update_", p)]], {
      nm <- trimws(input[[paste0(p, "_name")]] %||% "")
      chk <- opts_edit_check(sp, nm, targets(), isTRUE(input[[paste0(p, "_edit")]]), deps)
      if (identical(chk$action, "refuse")) {
        rp_alert(chk$msg, "error", session = session)
        return()
      }
      if (identical(chk$action, "confirm")) {
        if (!mp_confirm(paste0(p, "_confirm_edit"), title = "Edit beyond the selection",
          text = chk$msg, action_label = "Continue", danger = TRUE, skippable = TRUE)) return()
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

  # the Entrez filter only matters for a search that can reach NCBI
  shiny::observeEvent(input$bo_source, ignoreInit = TRUE, {
    shinyjs::toggle(id = "bo_entrez_group",
      condition = !identical(input$bo_source, "bundled"))
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
    sprintf("SELECT DISTINCT ID FROM %s WHERE %s = 1 AND ID IN (%s)", sp$queue_tbl, sp$lock, ph),
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

# one opts-table row from the modal inputs, with `base` supplying the fields the
# modal does not expose. `input` may be shiny's input object or a plain list.
opts_row <- function(pfx, nm, input, base) {
  switch(pfx,
    po = data.frame(pre_opts = nm, cpus = as.integer(input$po_cpus),
      memory = as.integer(input$po_mem), fastp = input$po_fastp %||% "",
      stringsAsFactors = FALSE),
    ao = data.frame(assemble_opts = nm, cpus = as.integer(input$ao_cpus),
      memory = as.integer(input$ao_mem), getOrganelle = input$ao_getorg %||% "",
      seeds_db = blank_to_na(input$ao_seeds %||% base$seeds_db),
      labels_db = blank_to_na(input$ao_labels %||% base$labels_db),
      max_paths = as.integer(input$ao_maxpaths),
      max_scaffolds = as.integer(input$ao_maxscaf),
      min_assembly_length = as.integer(input$ao_minlen),
      assembler = input$ao_assembler %||% "getorganelle",
      find_rrna = as.integer(isTRUE(input$ao_findrrna)),
      bowtie2 = input$ao_bt2 %||% "--local", spades = input$ao_spades %||% "--isolate",
      stringsAsFactors = FALSE),
    bo = data.frame(blast_opts = nm, run_blast = as.integer(isTRUE(input$bo_run)),
      entrez_query = input$bo_entrez %||% "",
      remote_blast = blast_source_cols(input$bo_source)$remote_blast,
      remote_fallback = blast_source_cols(input$bo_source)$remote_fallback,
      extra_opts = input$bo_extra %||% "",
      max_target_seqs = as.integer(input$bo_mts), stringsAsFactors = FALSE),
    no = data.frame(annotate_opts = nm, cpus = as.integer(input$no_cpus),
      memory = as.integer(input$no_mem), rrna_tool = input$no_tool,
      rrna_kingdom = input$no_kingdom, itsx_taxon = input$no_itsxtaxon %||% "all",
      itsx_opts = input$no_itsxopts %||% "", barrnap_opts = input$no_barrnap %||% "",
      pybarrnap_opts = input$no_pybarrnap %||% "", stringsAsFactors = FALSE),
    co = data.frame(curate_opts = nm, ruleset = input$co_ruleset, stringsAsFactors = FALSE)
  )
}

#' Warnings for reference-database paths a saved preset cannot use
#' @param row One assemble_opts row (other opts tables have no database columns).
#' @return character vector of messages, empty when everything is usable.
#' @noRd
unusable_db_warnings <- function(row) {
  out <- character(0)
  for (f in c(seeds_db = "Seed", labels_db = "Label")) {
    col <- names(which(c(seeds_db = "Seed", labels_db = "Label") == f))
    if (!(col %in% names(row))) next
    p <- row[[col]][1]
    if (is.na(p) || !nzchar(trimws(p))) next
    if (!usable_fasta(p)) {
      out <- c(out, sprintf(
        "%s database not found or empty: %s. Assembly will fail for samples on this set until it exists.",
        f, p))
    }
  }
  out
}

# build the full opts row (edited fields overriding the base) and upsert it,
# then assign the set to the target samples and re-queue the step.
opts_apply <- function(pfx, sp, input, deps, ids, cur) {
  if (length(ids) == 0) { shiny::removeModal(); return() }
  nm <- trimws(input[[paste0(pfx, "_name")]] %||% "")
  if (!nzchar(nm)) { rp_alert("Give the parameter set a name.", "warning"); return() }
  editing <- isTRUE(input[[paste0(pfx, "_edit")]])
  base <- if (!is.null(cur) && nrow(cur)) as.list(cur) else list()

  con <- deps$con(); on.exit(DBI::dbDisconnect(con))
  known <- DBI::dbGetQuery(con, sprintf("SELECT %s FROM %s", sp$key, sp$table))[[1]]
  warns <- character(0)
  saved <- app_write({
    if (editing || !(nm %in% known)) {
      row <- opts_row(pfx, nm, input, base)
      dplyr::rows_upsert(dplyr::tbl(con, sp$table), row, by = sp$key, in_place = TRUE, copy = TRUE)
      # a database that is not there yet is saved anyway (it may live on a share
      # that is not mounted, or be built later), but it is said out loud now
      # rather than surfacing as a failed assembly
      warns <- unusable_db_warnings(row)
    }
    for (id in ids) {
      DBI::dbExecute(con, sprintf("UPDATE %s SET %s = ? WHERE ID = ?", sp$link_tbl, sp$key),
        params = list(nm, id))
      DBI::dbExecute(con, sprintf("UPDATE %s SET %s = 1 WHERE ID = ?", sp$queue_tbl, sp$queue),
        params = list(id))
      if (sp$queue_tbl == "annotate") {
        cp <- DBI::dbGetQuery(con, "SELECT path, scaffold FROM review_units WHERE ID = ?",
          params = list(id))
        for (i in seq_len(nrow(cp))) review_follow_original(con, id, cp$path[i], cp$scaffold[i])
      }
    }
  })
  if (!saved) return()
  deps$refresh(deps$refresh() + 1)
  shiny::removeModal()
  # one alert, not one per warning: a second sweet alert replaces the first
  if (length(warns)) {
    msg <- paste0(mp_n(length(ids), "sample"), " switched to ", nm, " and set to Ready to run.")
    rp_alert(shiny::tagList(shiny::tags$p(msg), lapply(warns, shiny::tags$p)), "warning")
  } else {
    mp_opts_saved_toast(length(ids), nm)
  }
}
