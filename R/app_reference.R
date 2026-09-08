# The Reference tab: one row per clade group, a ladder to change the clade, and
# the build / reuse / skip actions.
#
# Every verb lives in the row, not in the control bar: each group is at its own
# point of the lifecycle and admits one sensible action, so a bar of verbs meant
# the user had to guess which rows a click would touch and how many were allowed.
#
# The table's key column is named ID so the app's shared selection, CSV export
# and click handlers work on it unchanged; the value is the group key, not a
# sample.

#' Reference-tab table data
#' @param con Open project DB connection.
#' @return data.frame with the columns the tab renders.
#' @noRd
reference_table <- function(con) {
  d <- tryCatch(DBI::dbGetQuery(con,
    "SELECT reference_switch,
       CASE COALESCE(reference_switch, 1)
         WHEN 2 THEN 'Ready' WHEN 3 THEN 'Build failed' WHEN 4 THEN 'Building'
         WHEN 5 THEN 'Left on current database' ELSE 'Needs a database' END AS status,
       ref_id AS ID, clade_name, clade_rank, n_samples, n_records,
       COALESCE(seeds_db, '') AS database, COALESCE(labels_db, '') AS labels_db,
       COALESCE(fail_reason, '') AS fail_reason, time_stamp,
       CASE COALESCE(reference_switch, 1)
         WHEN 3 THEN 'Try again' WHEN 2 THEN 'Change' WHEN 4 THEN 'Change'
         WHEN 5 THEN 'Change' ELSE 'Set up database' END AS action
     FROM reference ORDER BY ref_id"), error = function(e) NULL)
  if (is.null(d)) {
    d <- data.frame(reference_switch = integer(0), status = character(0), ID = character(0),
      clade_name = character(0), clade_rank = character(0), n_samples = integer(0),
      n_records = integer(0), database = character(0), fail_reason = character(0),
      time_stamp = integer(0), action = character(0), stringsAsFactors = FALSE)
  }
  # "Ready" has to mean the database is still on disk. The reference folder is
  # shared and can be deleted or moved long after a group was wired to it, and
  # nothing else ever revisits that decision, so a row wired to files that are
  # gone would keep reading Ready and offering only "Change".
  if (nrow(d) > 0) {
    gone <- !is.na(d$reference_switch) & d$reference_switch == 2 &
      !(vapply(d$database, reference_db_ok, logical(1)) &
        vapply(d$labels_db, reference_db_ok, logical(1)))
    d$status[gone] <- "Database missing"
    d$action[gone] <- "Set up database"
    d$fail_reason[gone] <- paste("The seed/label FASTAs this group was wired to are no longer",
      "there. Build the database again, or point the group at one that exists.")
  }
  d$labels_db <- NULL
  d
}

#' Reference-tab table data, sample view: one row per sample
#'
#' The group view answers "what will be built"; this one answers "where did my
#' sample end up", which the group view can only show by expanding a row. It is
#' also the surface bulk moves act on, because the app's selection machinery
#' keys on the table's ID column and here that is a sample.
#' @noRd
reference_sample_table <- function(con) {
  d <- tryCatch(DBI::dbGetQuery(con,
    "SELECT COALESCE(r.reference_switch, 1) AS reference_switch, s.ID, s.Taxon, s.sci_name,
       COALESCE(r.clade_name, 'Unassigned') AS ref_group,
       COALESCE(r.clade_rank, '') AS clade_rank,
       COALESCE(s.user_group, '') AS user_group,
       CASE WHEN s.ref_id IS NULL THEN ''
            WHEN s.ref_source = 'manual' THEN 'you' ELSE 'auto' END AS set_by,
       s.time_stamp
     FROM samples s LEFT JOIN reference r ON r.ref_id = s.ref_id
     ORDER BY s.ID"), error = function(e) NULL)
  if (is.null(d)) {
    d <- data.frame(reference_switch = integer(0), ID = character(0), Taxon = character(0),
      sci_name = character(0), ref_group = character(0), clade_rank = character(0),
      user_group = character(0), set_by = character(0), time_stamp = integer(0),
      stringsAsFactors = FALSE)
  }
  d
}

#' Every group's members, for the group view's expandable rows
#' @noRd
reference_members <- function(con) {
  tryCatch(DBI::dbGetQuery(con,
    "SELECT ref_id, ID, COALESCE(Taxon, '') AS Taxon, COALESCE(user_group, '') AS user_group,
       COALESCE(ref_source, '') AS ref_source
     FROM samples WHERE ref_id IS NOT NULL ORDER BY ID"),
    error = function(e) data.frame(ref_id = character(0), ID = character(0), Taxon = character(0),
      user_group = character(0), ref_source = character(0), stringsAsFactors = FALSE))
}

#' The member list shown when a group row is expanded
#'
#' Read-only on purpose: everything that changes an assignment lives in the
#' sample view, so there is one place to look for it rather than two.
#' @noRd
reference_member_html <- function(m) {
  if (nrow(m) == 0) {
    return(shiny::div(style = "padding:0.5em 1em; color:#888;", "No samples in this group."))
  }
  tag <- function(txt, colour) shiny::span(txt,
    style = sprintf("margin-left:0.6em; font-size:0.8em; color:%s;", colour))
  shiny::div(style = "padding:0.4em 1em 0.8em 2.5em;",
    lapply(seq_len(nrow(m)), function(i) shiny::div(
      style = "padding:1px 0;",
      shiny::strong(m$ID[i]),
      shiny::span(m$Taxon[i], style = "margin-left:0.8em; color:#555;"),
      if (nzchar(m$user_group[i])) tag(paste("User group", m$user_group[i]), "#1a7f5a"),
      if (identical(m$ref_source[i], "manual")) tag("assigned by you", "#8a6d1a"))))
}

#' Preview modal: what Propose groups would change, before anything is written
#'
#' @param diff From [reference_group_diff()].
#' @param n_old,n_new Group counts before and after.
#' @param warn Character vector of warning lines (clades over the record ceiling).
#' @noRd
reference_diff_modal <- function(diff, n_old, n_new, warn = character(0)) {
  hand <- diff$ID[diff$hand_set]
  rows <- if (nrow(diff) == 0) {
    shiny::helpText("No sample changes group.")
  } else {
    shiny::div(style = "max-height:45vh; overflow:auto; font-size:0.9em;",
      shiny::tags$table(class = "table table-sm table-striped",
        shiny::tags$thead(shiny::tags$tr(shiny::tags$th("Sample"), shiny::tags$th("Now"),
          shiny::tags$th("Proposed"), shiny::tags$th(""))),
        shiny::tags$tbody(lapply(seq_len(nrow(diff)), function(i) shiny::tags$tr(
          shiny::tags$td(diff$ID[i]),
          shiny::tags$td(diff$from_group[i]),
          shiny::tags$td(shiny::strong(diff$to_group[i])),
          shiny::tags$td(style = "color:#8a6d1a;",
            if (diff$hand_set[i]) "you assigned this" else ""))))))
  }
  shiny::modalDialog(
    title = "Proposed reference groups", size = "l",
    rows,
    shiny::div(style = "margin-top:0.6em;",
      sprintf("%d group(s) now, %d proposed.", n_old, n_new)),
    if (length(warn) > 0) {
      shiny::div(style = "margin-top:0.6em; color:#a33;",
        lapply(warn, function(w) shiny::div(w)))
    },
    if (length(hand) > 0) {
      shiny::tagList(
        shiny::hr(),
        shiny::checkboxGroupInput("ref_diff_keep", "Keep any of these together instead:",
          choices = hand, inline = TRUE),
        opts_help("Making a user group out of them stops any future proposal from ",
          "splitting them up. Nothing is applied; press Propose groups again after."),
        shiny::actionButton("ref_diff_group", "Make a user group", class = "btn-sm"))
    },
    footer = shiny::tagList(
      shiny::actionButton("ref_diff_apply", "Apply", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Pick the group a bulk move sends samples to
#' @param groups data.frame(ref_id, clade_name) of existing groups.
#' @noRd
reference_move_modal <- function(groups, n) {
  choices <- c(stats::setNames(groups$ref_id, groups$clade_name), "Unassigned" = "")
  shiny::modalDialog(
    title = sprintf("Move %d sample(s) to a reference group", n),
    shinyWidgets::prettyRadioButtons("ref_move_pick", NULL, choices = choices,
      shape = "round", status = "primary"),
    opts_help("The move is yours, not the proposer's, so the next Propose groups ",
      "will flag it before changing it. To make it permanent, select the samples ",
      "and press Keep together."),
    footer = shiny::tagList(
      shiny::actionButton("ref_move_apply", "Move", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' One rung of the clade ladder, as a radio choice label
#' @noRd
rung_label <- function(name, rank, n_samples, n_records, chosen) {
  recs <- if (is.na(n_records)) "records unknown" else paste0(format(n_records, big.mark = ","), " records")
  sprintf("%s [%s] - %d sample(s), %s%s", name, rank, n_samples, recs,
    if (isTRUE(chosen)) "  (current)" else "")
}

#' Clade ladder modal: pick a rung, or search the taxdump for anything else
#'
#' @param ref_id Group key. @param clade_name Current clade.
#' @param ladder data.frame from [group_ladder()].
#' @param chosen_taxid The group's current clade taxid.
#' @param n_group Samples in the group.
#' @param index_ready Is the taxdump name-search index already built? When it is
#'   not, the first search pays for building it, which is worth saying up front.
#' @noRd
reference_ladder_modal <- function(ref_id, clade_name, ladder, chosen_taxid, n_group,
                                   index_ready = TRUE) {
  ladder <- ladder[order(-ladder$level), , drop = FALSE]
  choices <- stats::setNames(as.character(ladder$taxid),
    vapply(seq_len(nrow(ladder)), function(i) rung_label(ladder$name[i], ladder$rank[i],
      ladder$n_samples[i], ladder$n_records[i],
      !is.na(chosen_taxid) && ladder$taxid[i] == chosen_taxid), character(1)))
  shiny::modalDialog(
    title = sprintf("Reference clade for '%s' (%d sample(s))", ref_id, n_group), size = "l",
    shinyWidgets::prettyRadioButtons("ref_rung", "Shared lineage (deepest first):",
      choices = choices,
      selected = if (!is.na(chosen_taxid)) as.character(chosen_taxid) else NULL,
      shape = "round", status = "primary"),
    opts_help("A clade with too few records makes a weak seed database; one with ",
      "too many stops being specific. Picking a rung that covers only some of the ",
      "group's samples offers to split the rest into their own group."),
    shiny::textInput("ref_search", "Or search NCBI for another clade:", width = "100%",
      placeholder = "e.g. Echinodermata"),
    if (!isTRUE(index_ready)) {
      opts_help("The first search builds a searchable copy of the local NCBI ",
        "taxonomy, which takes a minute or two. Every search after that is instant.")
    },
    shiny::uiOutput("ref_search_results"),
    footer = shiny::tagList(
      shiny::actionButton("ref_clade_apply", "Use this clade", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Use-an-existing-database modal: manifest entries plus typed paths
#' @param manifest From [read_ref_manifest()].
#' @param label What to call the target in the title: the clade name for one
#'   group, "N selected groups" for several. Never the raw group key, which the
#'   table no longer shows.
#' @noRd
reference_existing_modal <- function(manifest, label) {
  listing <- if (nrow(manifest) == 0) {
    shiny::helpText("No databases recorded in this reference folder yet.")
  } else {
    # %|NA|% is not vectorised (R/utils.R), so a manifest with more than one row
    # (or an NA record count) needs ifelse() here, not the scalar %|NA|% helper.
    rank_txt <- ifelse(is.na(manifest$rank), "?", manifest$rank)
    recs_txt <- ifelse(is.na(manifest$n_records), "records unknown",
      paste0(format(manifest$n_records, big.mark = ","), " records"))
    shinyWidgets::prettyRadioButtons("ref_existing_pick", "Databases in the reference folder:",
      choices = stats::setNames(as.character(seq_len(nrow(manifest))),
        sprintf("%s [%s] - %s, built %s", manifest$clade, rank_txt, recs_txt, manifest$built)),
      shape = "round", status = "primary")
  }
  shiny::modalDialog(
    title = sprintf("Use an existing database for %s", label), size = "l",
    listing,
    shiny::textInput("ref_existing_seeds", "Seed FASTA path:", width = "100%"),
    shiny::textInput("ref_existing_labels", "Label FASTA path:", width = "100%"),
    opts_help("Choosing a listed database fills these paths; you can also type ",
      "paths to a database built somewhere else."),
    footer = shiny::tagList(
      shiny::actionButton("ref_existing_apply", "Use this database", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

# Section headings in the settings modal. Bootstrap's h5 is lighter than a
# .control-label, so a bare h5 above the Minimum/Maximum labels reads as the
# smaller of the two; 600 puts the heading back on top.
SETTINGS_HEAD <- "font-weight:600; margin:2px 0 6px;"

#' Settings modal: the reference folder, and the record window Propose groups uses
#'
#' @param cfg From [ref_config()].
#' @param blastdb_version From [blastdb_provenance()]; "" when nothing has run yet.
#' @noRd
reference_settings_modal <- function(cfg, blastdb_version = "") {
  shiny::modalDialog(
    title = "Reference settings", size = "l",
    shiny::h5("Reference database folder", style = SETTINGS_HEAD),
    shiny::textInput("ref_dir_path", NULL, width = "100%", value = cfg$ref_dir),
    opts_help("Databases built here are reused by every project pointed at this ",
      "folder, including a shared drive. Changing it does not move databases ",
      "already built, and groups already wired keep working."),
    shiny::h5("GenBank records per suggested clade", style = paste0(SETTINGS_HEAD, " margin-top:18px;")),
    # these two are read as a pair, so they sit side by side; naming the measure
    # once in the heading is what lets the labels be one word each
    shiny::div(style = "display:flex; gap:1.25em;",
      shiny::numericInput("ref_min", "Minimum", value = cfg$min_records,
        min = 1, step = 100, width = "150px"),
      shiny::numericInput("ref_max", "Maximum", value = cfg$max_records,
        min = 1, step = 1000, width = "150px")),
    opts_help("Propose groups climbs each sample's lineage and stops at the first ",
      "clade holding at least the minimum. Below that there is too little published ",
      "sequence to seed an assembly; past the maximum a clade is so broad that its ",
      "seeds stop being specific to your organism."),
    shiny::h5("Local BLAST panel", style = SETTINGS_HEAD),
    shiny::div(
      style = "font-size:0.85em; color:#666; margin-bottom:0.4em;",
      "Path prefix of a downloaded panel, without a file extension. ",
      "Leave empty and RiboPilot uses a panel shipped with the image or ",
      "installed locally if it finds one, otherwise it searches NCBI. ",
      "To search NCBI instead, use Remote BLAST in the BLAST options; the ",
      "environment variable RIBOPILOT_NO_LOCAL_BLAST=1 does the same for every ",
      "project when set before launching the app."),
    shiny::textInput("ref_blastdb_path", NULL, width = "100%",
      value = cfg$blastdb %||% "",
      placeholder = "/data/blastdb/2026.08.18/ribopilot_rrna"),
    if (nzchar(blastdb_version %||% "")) {
      shiny::div(style = "font-size:0.85em; color:#666; margin-top:0.4em;",
        paste("Last run searched against:", blastdb_version))
    },
    footer = shiny::tagList(
      shiny::actionButton("ref_dir_apply", "Save", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Row action modal: everything that can be done to a reference group
#'
#' The table's action cell is the only door onto these, so each choice is a full
#' sentence naming its consequence. "Leave these samples..." is the only place in
#' the app that says what skipping actually costs, which is why it is spelled out
#' rather than labelled "Skip".
#'
#' @param ids Group keys the action will apply to.
#' @param label Display name for the title (clade name, or "N selected groups").
#' @noRd
reference_action_modal <- function(ids, label) {
  choice <- function(id, text, help, primary = FALSE) {
    shiny::div(style = "margin-bottom:0.9em;",
      shiny::actionButton(id, text, width = "100%",
        class = if (primary) "btn-primary btn-sm" else "btn-sm"),
      # nested: opts_help's default -6px top margin is tuned to pull the hint up
      # under a form field's own bottom margin. A bare button has none, so the
      # default drags the text into the button.
      opts_help(help, nested = TRUE))
  }
  one <- length(ids) == 1
  m <- shiny::modalDialog(
    title = sprintf("Reference database for %s", label),
    choice("ref_act_build", "Download sequences from GenBank",
      paste(if (one) "Builds a seed database for this clade from GenBank."
        else "Builds a seed database for each of these clades from GenBank.",
        "This is the normal choice."), primary = TRUE),
    choice("ref_act_existing", "Use a database I already have",
      paste("Point these samples at seed and label FASTA files already on disk,",
        "for example one a colleague built.")),
    # the ladder walks one group's shared lineage, so it has no meaning across a
    # multi-group selection
    if (one) {
      choice("ref_act_clade", "Change which clade this group uses",
        "Widen or narrow the group, or pick a different clade entirely.")
    },
    choice("ref_act_skip", "Leave these samples on the built-in fish database",
      paste("Nothing is built. These samples assemble against RiboPilot's",
        "ray-finned fish reference, which is usually wrong outside fish.")),
    footer = shiny::modalButton("Cancel"))
  # This dialog is a short stack of buttons, so the app's 600px modal floor is
  # dead space that stretches every label across it. Narrow it and let the hints
  # wrap. modalDialog() takes no class argument, hence the tag surgery.
  htmltools::tagQuery(m)$find(".modal-dialog")$addClass("rp-modal-narrow")$allTags()
}

#' Wire the Reference tab
#'
#' @param deps list(con = fn -> DBI con, project = fn -> chr, refresh = reactiveVal,
#'   selected_ids = fn -> chr, table_data = reactive -> data.frame).
#' @return list(open_ladder, apply_clade, poll, do_build, do_skip, open_existing,
#'   open_move, do_keep_together, do_dissolve),
#'   each taking the group key(s) it acts on.
#' @noRd
reference_handlers <- function(input, output, session, deps) {
  ladder_target <- shiny::reactiveVal(NULL)
  existing_target <- shiny::reactiveVal(NULL)
  build_targets <- shiny::reactiveVal(character(0))
  action_target <- shiny::reactiveVal(character(0))  # groups the row action applies to
  action_label <- shiny::reactiveVal(NULL)           # what to call them in a title
  picked <- shiny::reactiveVal(NULL)   # taxid chosen from the search box
  ref_proc <- shiny::reactiveVal(NULL)
  ref_log <- shiny::reactiveVal("")
  ref_state <- shiny::reactiveVal(nf_progress_empty())
  headless_open <- shiny::reactiveVal(FALSE)
  gen_scripts <- shiny::reactiveVal(character(0))   # scripts this session generated
  submitted <- shiny::reactiveVal(character(0))     # groups with a cluster job out
  resubmit_ok <- shiny::reactiveVal(FALSE)
  ref_base <- shiny::reactiveVal(NULL)              # timestamp for one submission
  proposal <- shiny::reactiveVal(NULL)              # suggestion awaiting Apply

  with_con <- function(fn) {
    con <- deps$con(); on.exit(DBI::dbDisconnect(con)); fn(con)
  }
  bump <- function() deps$refresh(shiny::isolate(deps$refresh()) + 1)
  # read fresh each time so tests can substitute a runner that spawns nothing
  count_runner <- function() getOption("ribopilot.count_runner", nf_run_blocking)
  notify <- function(msg, type = "message") rp_alert(msg, type, session = session)

  poll_now <- function() {
    tryCatch(with_con(poll_reference_builds), error = function(e) {
      # poll_reference_builds() stop()s (e.g. wire_reference_db() finding no
      # base assemble_opts preset) before it can mark the row it was working
      # on; left at reference_switch = 4 every later tick hits the same error
      # and the user watches "running" forever. Fail every currently-running
      # row with this message instead, so it surfaces and stops spinning.
      msg <- conditionMessage(e)
      tryCatch(with_con(function(con) DBI::dbExecute(con,
        "UPDATE reference SET reference_switch = 3, pid = NULL, fail_reason = ?, time_stamp = ?
         WHERE reference_switch = 4",
        params = list(msg, as.integer(Sys.time())))), error = function(e2) NULL)
      notify(msg, "error")
      invisible(0L)
    })
  }

  # taxonomy writers store the lineage names only; the rank-aligned taxid path
  # this tab walks is derived on demand, so back-fill before every read
  ladders_now <- function(con) {
    tryCatch(ensure_sample_lineage_taxids(con,
      taxdump_dir = project_taxdump_dir(deps$project())), error = function(e) NULL)
    sample_ladders(con)
  }

  # ---- Suggest ----
  shiny::observeEvent(input$ref_suggest_btn, {
    shinyWidgets::sendSweetAlert(title = "Suggesting clades...",
      text = "Reading sample lineages and checking GenBank record counts. This can take a moment.",
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        # bump() and the success alert are INSIDE the tryCatch too: the "Suggesting..."
        # alert has no dismiss button, so an error raised after the DB write (e.g. from
        # bump() or the alert call itself) must still reach the "failed" branch, or the
        # user is stuck looking at a popup they cannot close.
        outcome <- tryCatch({
          res <- with_con(function(con) {
            cfg <- ref_config(con)
            l <- ladders_now(con)
            if (nrow(l) == 0) stop("No samples have a resolved lineage yet. Assign taxonomy first.",
              call. = FALSE)
            # one start rank for both: the walk must only ask about rungs the
            # count run was told to count, or every answer comes back NA
            start_rank <- "family"
            counts <- run_clade_counts(deps$project(), cfg$ref_dir, candidate_rungs(l, start_rank),
              api_key = project_ncbi_key(deps$project()), runner = count_runner())
            s <- suggest_reference_groups(l, cache_count_fn(cfg$ref_dir),
              cfg$min_records, cfg$max_records, start_rank = start_rank,
              user_groups = sample_user_groups(con))
            if (nrow(s$groups) == 0) {
              # the counting run is the only thing that could have failed here,
              # so its log is the message, not a pointer to a file. A task that
              # failed and was ignored says why only in its own stderr, and that
              # is nearly always the real answer, so it leads.
              log <- attr(counts, "log")
              errs <- attr(counts, "errors") %||% character(0)
              tail <- log_tail(log, 5)
              stop("Could not count GenBank records for any clade (the counting run returned no results).",
                if (length(errs) > 0) paste0(" Task error: ", errs[1]) else "",
                if (!is.null(log)) paste0(" Log: ", log, ".") else "",
                if (nzchar(tail)) paste0(" Last lines: ", tail) else "",
                call. = FALSE)
            }
            # nothing is written yet: the user sees what would move first, and a
            # sample they placed by hand gets flagged rather than quietly undone
            over <- s$groups[!is.na(s$groups$n_records) &
                s$groups$n_records > cfg$max_records, , drop = FALSE]
            list(prop = s, diff = reference_group_diff(con, s$groups, s$members),
              n_old = nrow(DBI::dbGetQuery(con, "SELECT ref_id FROM reference")),
              warn = if (nrow(over) == 0) character(0) else sprintf(
                "%s holds %s records, above the maximum of %s.", over$clade_name,
                format(over$n_records, big.mark = ","), format(cfg$max_records, big.mark = ",")))
          })
          proposal(res$prop)
          shinyWidgets::closeSweetAlert(session = session)
          shiny::showModal(reference_diff_modal(res$diff, res$n_old, nrow(res$prop$groups),
            res$warn))
          NULL
        }, error = function(e) e)
        if (inherits(outcome, "error")) {
          shinyWidgets::sendSweetAlert(title = "Suggestion failed",
            text = conditionMessage(outcome), type = "error")
        }
      })
    }, delay = 0.1)
  })

  # ---- apply / refuse a proposal ----
  shiny::observeEvent(input$ref_diff_apply, {
    s <- proposal()
    if (is.null(s)) { notify("Nothing to apply.", "warning"); return() }
    res <- tryCatch(with_con(function(con) write_reference_groups(con, s$groups, s$members)),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    proposal(NULL)
    shiny::removeModal(); bump()
    # land on the Groups view, then go straight to giving the new groups a
    # database: one dialog for all of them, since one Build runs them together
    shinyWidgets::updateRadioGroupButtons(session, "ref_view", selected = "Groups")
    todo <- tryCatch(with_con(function(con) DBI::dbGetQuery(con,
      "SELECT ref_id, clade_name FROM reference WHERE COALESCE(reference_switch, 1) = 1 ORDER BY ref_id")),
      error = function(e) data.frame())
    todo <- todo[todo$ref_id %in% s$groups$ref_id, , drop = FALSE]
    if (nrow(todo) == 0) {
      notify(sprintf("Applied %d reference group(s); every group already has a database.", nrow(s$groups)))
      return()
    }
    lbl <- if (nrow(todo) == 1) todo$clade_name[1] else sprintf("%d new groups", nrow(todo))
    action_label(lbl); action_target(todo$ref_id)
    shiny::showModal(reference_action_modal(todo$ref_id, lbl))
  })

  shiny::observeEvent(input$ref_diff_group, {
    ids <- input$ref_diff_keep %||% character(0)
    if (length(ids) < 2) {
      notify("Pick at least two samples: a user group is what keeps them together.", "warning")
      return()
    }
    res <- tryCatch(with_con(function(con) {
      lab <- next_user_group_label(con)
      set_user_group(con, ids, lab)
      lab
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    proposal(NULL)
    shiny::removeModal(); bump()
    notify(sprintf("User group %s created from %d sample(s). Press Propose groups again.",
      res, length(ids)))
  })

  # ---- sample view: move, keep together, dissolve ----
  move_target <- shiny::reactiveVal(character(0))
  open_move <- function(ids) {
    if (length(ids) == 0) { notify("Select the sample(s) to move.", "warning"); return() }
    move_target(ids)
    res <- tryCatch(with_con(function(con) DBI::dbGetQuery(con,
      "SELECT ref_id, clade_name FROM reference ORDER BY clade_name")), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    if (nrow(res) == 0) {
      notify("No reference groups yet. Press Propose groups first.", "warning"); return()
    }
    shiny::showModal(reference_move_modal(res, length(ids)))
  }
  shiny::observeEvent(input$ref_move_apply, {
    ids <- move_target()
    key <- input$ref_move_pick %||% ""
    res <- tryCatch(with_con(function(con) {
      assign_samples_ref(con, ids, if (nzchar(key)) key else NA_character_)
      if (!nzchar(key)) return(0L)
      # a clade that does not contain a sample makes a seed database that cannot
      # match it, so say so; the move still stands, because the user asked for it
      tid <- DBI::dbGetQuery(con, "SELECT clade_taxid FROM reference WHERE ref_id = ?",
        params = list(key))$clade_taxid
      if (length(tid) == 0 || is.na(tid[1])) return(0L)
      length(setdiff(ids, ladder_coverage(ladders_now(con), tid[1])))
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    if (res > 0) {
      notify(sprintf("Moved, but %d sample(s) are not inside that clade.", res), "warning")
    } else {
      notify(sprintf("Moved %d sample(s).", length(ids)))
    }
  })

  do_keep_together <- function(ids) {
    if (length(ids) < 2) {
      notify("Select at least two samples to keep together.", "warning"); return()
    }
    res <- tryCatch(with_con(function(con) {
      lab <- next_user_group_label(con)
      set_user_group(con, ids, lab)
      lab
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    bump()
    notify(sprintf("User group %s: these %d sample(s) will not be split up.", res, length(ids)))
  }

  do_dissolve <- function(ids) {
    if (length(ids) == 0) { notify("Select the sample(s) to release.", "warning"); return() }
    res <- tryCatch(with_con(function(con) set_user_group(con, ids, NA_character_)),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    bump()
    notify(sprintf("Released %d sample(s); the next proposal may regroup them.", length(ids)))
  }

  # ---- ladder modal ----
  # `picked` tracks whichever of the ladder radio or the search-result radio the
  # user touched LAST; Shiny retains an input's value after its widget is gone
  # (the search results renderUI clears, the modal closes), so reading either
  # input directly risks applying a stale search pick over a fresh rung click.
  shiny::observeEvent(input$ref_rung, ignoreNULL = TRUE, { picked(input$ref_rung) })
  shiny::observeEvent(input$ref_search_pick, ignoreNULL = TRUE, { picked(input$ref_search_pick) })

  open_ladder <- function(ref_id) {
    if (length(ref_id) == 0) { notify("Select a group first.", "warning"); return() }
    ref_id <- ref_id[1]
    picked(NULL)
    ladder_target(ref_id)
    res <- tryCatch(with_con(function(con) {
      row <- DBI::dbGetQuery(con, "SELECT * FROM reference WHERE ref_id = ?", params = list(ref_id))
      if (nrow(row) == 0) stop("Group not found.", call. = FALSE)
      ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE ref_id = ?", params = list(ref_id))$ID
      cfg <- ref_config(con)
      # cache-only: opening the modal must not shell out. Suggest pre-fills the
      # cache for these rungs; anything else renders as "records unknown"
      l <- group_ladder(ladders_now(con), ids, count_fn = cache_count_fn(cfg$ref_dir))
      reference_ladder_modal(ref_id, row$clade_name[1], l, row$clade_taxid[1], length(ids),
        index_ready = file.exists(tax_index_path(project_taxdump_dir(deps$project()))))
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  }
  shiny::observeEvent(input$ref_open_clade, {
    i <- suppressWarnings(as.integer(input$ref_open_clade))
    d <- deps$table_data()
    if (!is.na(i) && i >= 1 && i <= nrow(d)) open_ladder(d$ID[i])
  })

  # taxdump search inside the ladder modal
  q <- shiny::debounce(shiny::reactive(input$ref_search), 300)
  output$ref_search_results <- shiny::renderUI({
    query <- trimws(q() %||% "")
    if (nchar(query) < 2) return(NULL)
    r <- tax_search(project_taxdump_dir(deps$project()), query, 20)
    if (nrow(r) == 0) return(shiny::helpText("No matches in the local taxdump."))
    shinyWidgets::prettyRadioButtons("ref_search_pick", NULL,
      choices = stats::setNames(as.character(r$taxid),
        sprintf("%s   [%s]   (taxid %s)", r$name, ifelse(is.na(r$rank), "?", r$rank), r$taxid)),
      shape = "round", status = "primary")
  })

  shiny::observeEvent(input$ref_clade_apply, {
    tid <- suppressWarnings(as.integer(picked() %||% input$ref_rung %||% NA))
    if (is.na(tid)) { notify("Pick a clade first.", "warning"); return() }
    # the clade check walks the taxdump and asks GenBank for a record count,
    # which takes a few seconds; say so rather than leaving the button dead
    shinyWidgets::sendSweetAlert(title = "Applying clade...",
      text = "Checking which samples it covers and counting GenBank records.",
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    # read the target HERE: withReactiveDomain sets a session, not a reactive
    # context, so a reactiveVal read inside the later() callback would error
    target <- ladder_target()
    later::later(function() {
      shiny::withReactiveDomain(session, apply_clade(tid, target))
    }, delay = 0.1)
  })

  apply_clade <- function(tid, ref_id) {
    res <- tryCatch(with_con(function(con) {
      ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE ref_id = ?",
        params = list(ref_id))$ID
      ladders <- ladders_now(con)
      covered <- intersect(ladder_coverage(ladders, tid), ids)
      # ranks outside the seven-rank ladder (superclass, subphylum, ...) are
      # still valid clades, so fall back to the full ancestor path
      if (length(covered) == 0) {
        covered <- intersect(
          lineage_coverage(con, tid, project_taxdump_dir(deps$project())), ids)
      }
      # a clade that covers none of the group's samples must not be applied to
      # all of them anyway -- that would silently mis-clade every sample
      if (length(covered) == 0) {
        stop("That clade does not cover any sample in this group; pick a shared rung.", call. = FALSE)
      }
      info <- ladders[ladders$taxid == tid, , drop = FALSE]
      if (nrow(info) == 0) {
        i <- tk_taxid_info(tid, project_taxdump_dir(deps$project()))
        info <- data.frame(name = i$sci_name[1], rank = i$rank[1], stringsAsFactors = FALSE)
      }
      cfg <- ref_config(con)
      run_clade_counts(deps$project(), cfg$ref_dir,
        data.frame(ref_id = ref_id_from_name(info$name[1]), taxid = as.integer(tid),
          clade_name = info$name[1], query = clade_query(info$name[1]),
          stringsAsFactors = FALSE),
        api_key = project_ncbi_key(deps$project()), runner = count_runner())
      n <- cache_count_fn(cfg$ref_dir)(tid, info$name[1])
      # a clade change re-keys the group, so the follow-up work uses the key
      # set_reference_clade() actually wrote, not the one the modal opened on
      new_key <- set_reference_clade(con, ref_id = ref_id, clade_taxid = tid,
        clade_name = info$name[1], clade_rank = info$rank[1], n_records = n,
        member_ids = covered)

      # samples the chosen clade does not cover are not left inside a clade
      # that does not contain them: split them into a group of their own,
      # keyed on their own shared lineage (deepest rung all of them pass through)
      leftover <- setdiff(ids, covered)
      left_name <- NULL
      if (length(leftover) > 0) {
        gl <- group_ladder(ladders, leftover)
        if (nrow(gl) > 0) {
          lg <- gl[order(-gl$level), , drop = FALSE][1, ]
          left_key <- unique_split_ref_id(con, lg$name, lg$taxid, avoid = new_key)
          set_reference_clade(con, ref_id = left_key, clade_taxid = lg$taxid, clade_name = lg$name,
            clade_rank = lg$rank, n_records = NA_integer_, member_ids = leftover)
          left_name <- lg$name
        }
      }
      list(name = info$name[1], left_n = length(leftover), left_name = left_name)
    }), error = function(e) e)
    # the "Applying..." alert has no dismiss button, so every exit closes it
    if (inherits(res, "error")) {
      shinyWidgets::sendSweetAlert(title = "Could not apply that clade",
        text = conditionMessage(res), type = "error")
      return(invisible(NULL))
    }
    shinyWidgets::closeSweetAlert(session = session)
    shiny::removeModal()
    bump()
    if (!is.null(res$left_name)) {
      notify(sprintf("%d sample(s) outside %s were split into a new group (%s); open it to pick a clade.",
        res$left_n, res$name, res$left_name), "warning")
    } else {
      notify(sprintf("Group now uses %s.", res$name))
    }
    invisible(NULL)
  }

  # ---- skip ----
  do_skip <- function(ids) {
    if (length(ids) == 0) { notify("Select a group first.", "warning"); return() }
    res <- tryCatch(with_con(function(con) { for (r in ids) skip_reference_group(con, r); length(ids) }),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    bump()
    notify(sprintf("Left %d group(s) on the built-in fish database.", res))
  }

  # ---- use existing ----
  open_existing <- function(ids, label = NULL) {
    if (length(ids) == 0) { notify("Select a group first.", "warning"); return() }
    existing_target(ids)
    res <- tryCatch(with_con(function(con) reference_existing_modal(
      read_ref_manifest(ref_config(con)$ref_dir),
      label %||% sprintf("%d selected groups", length(ids)))), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  }
  shiny::observeEvent(input$ref_existing_pick, {
    i <- suppressWarnings(as.integer(input$ref_existing_pick))
    with_con(function(con) {
      m <- read_ref_manifest(ref_config(con)$ref_dir)
      if (!is.na(i) && i >= 1 && i <= nrow(m)) {
        shiny::updateTextInput(session, "ref_existing_seeds", value = m$seeds[i])
        shiny::updateTextInput(session, "ref_existing_labels", value = m$labels[i])
      }
    })
  })
  shiny::observeEvent(input$ref_existing_apply, {
    res <- tryCatch(with_con(function(con) {
      n <- 0L
      for (r in existing_target()) {
        n <- n + use_existing_reference(con, r,
          input$ref_existing_seeds %||% "", input$ref_existing_labels %||% "")
      }
      n
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    notify(sprintf("Pointed %d sample(s) at the existing database.", res))
  })

  # ---- build: one Nextflow run for every selected group ----
  # the job script is written now but runs much later, so the manifest it reads
  # has to be on disk before the modal closes
  ref_submit_script <- function() {
    # the settings observer below regenerates this script several times while the
    # modal is open; the base must NOT move with it, or the file finally saved is
    # named for a different moment than the log paths written inside it
    base <- shiny::isolate(ref_base()) %||% job_basename("ribopilot_ref")
    s <- with_con(function(con) {
      cfg <- ref_config(con)
      ids <- reference_missing_dbs(cfg$ref_dir, build_targets())
      if (length(ids) == 0) return("")
      man <- write_build_manifest(deps$project(), reference_rows(con, ids))
      build_reference_submit_script(deps$project(), man, cfg$ref_dir,
        api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
        resume = isTRUE(input$ref_nf_resume %||% TRUE), sched = detect_scheduler(),
        base = base)
    })
    if (nzchar(s)) gen_scripts(unique(c(shiny::isolate(gen_scripts()), s)))
    s
  }

  do_build <- function(ids) {
    if (length(ids) == 0) { notify("Select the group(s) to build.", "warning"); return() }
    build_targets(ids)
    ref_state(nf_progress_empty())
    res <- tryCatch(with_con(function(con) {
      cfg <- ref_config(con)
      # count what will actually run: a group whose database is already on disk
      # is wired from it, never built
      todo <- reference_missing_dbs(cfg$ref_dir, ids)
      list(todo = todo, ref_dir = cfg$ref_dir)
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    # a group whose database is now on disk has no job outstanding any more
    submitted(setdiff(shiny::isolate(submitted()), setdiff(ids, res$todo)))
    if (length(res$todo) == 0) {
      wired <- tryCatch(with_con(function(con) start_reference_build(con, deps$project(), ids)),
        error = function(e) e)
      if (inherits(wired, "error")) { notify(conditionMessage(wired), "error"); return() }
      bump()
      notify(sprintf("%d group(s) already had a database in the reference folder; wired those.",
        length(ids)))
      return()
    }
    hl <- app_headless()
    sched <- if (hl) detect_scheduler() else NULL
    # one timestamp per opening of the modal; see ref_submit_script()
    if (hl) ref_base(job_basename("ribopilot_ref"))
    script <- if (hl) tryCatch(ref_submit_script(), error = function(e) "") else ""
    headless_open(hl)
    shiny::showModal(reference_run_modal(length(res$todo), res$ref_dir,
      headless = hl, sched = sched, script = script,
      api_key = project_ncbi_key(deps$project())))
  }

  # ---- row action: one door onto build / reuse / re-clade / skip ----
  # The clicked row folds into the current selection when it is part of it, so a
  # bulk action is still one click and no action is ever refused after the click.
  shiny::observeEvent(input$ref_row_action, {
    i <- suppressWarnings(as.integer(input$ref_row_action))
    d <- deps$table_data()
    if (is.na(i) || i < 1 || i > nrow(d)) return()
    sel <- deps$selected_ids()
    ids <- if (d$ID[i] %in% sel) sel else d$ID[i]
    lbl <- if (length(ids) == 1) d$clade_name[i] else sprintf("%d selected groups", length(ids))
    action_label(lbl)
    action_target(ids)
    shiny::showModal(reference_action_modal(ids, lbl))
  })
  shiny::observeEvent(input$ref_act_build, { shiny::removeModal(); do_build(action_target()) })
  shiny::observeEvent(input$ref_act_existing, {
    shiny::removeModal(); open_existing(action_target(), action_label())
  })
  shiny::observeEvent(input$ref_act_clade, { shiny::removeModal(); open_ladder(action_target()) })
  shiny::observeEvent(input$ref_act_skip, { shiny::removeModal(); do_skip(action_target()) })

  # the script is built before the modal renders, so its settings inputs are
  # still NULL; regenerate it whenever one of them changes, or the saved job
  # would silently ignore the API key and the RefSeq tick
  shiny::observeEvent(list(input$ref_api_key, input$ref_refseq, input$ref_nf_resume), {
    if (!isTRUE(headless_open())) return()
    s <- tryCatch(ref_submit_script(), error = function(e) NULL)
    if (!is.null(s) && nzchar(s)) {
      shiny::updateTextAreaInput(session, "ref_submit_script_text", value = s)
    }
  }, ignoreInit = TRUE)

  build_args <- function(con) {
    nf_reference_command(deps$project(), "BUILD", reference_manifest_path(deps$project()),
      ref_config(con)$ref_dir, ref_config(con)$ref_dir,
      api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
      resume = isTRUE(input$ref_nf_resume %||% TRUE))
  }
  output$ref_nf_cmd_text <- shiny::renderText({
    args <- tryCatch(with_con(build_args), error = function(e) NULL)
    if (is.null(args)) return("")
    paste(c("nextflow", args), collapse = " ")
  })

  # ---- headless: save / submit the job script instead of running locally ----
  write_ref_script <- function() {
    typed <- input$ref_submit_script_text
    typed <- if (is.null(typed)) "" else paste(typed, collapse = "\n")
    # an untouched textarea can still hold a script generated from settings the
    # user has since changed (the regenerated one is echoed back a beat later),
    # and an emptied one holds nothing worth writing: only a hand-edited script
    # beats a freshly generated one
    edited <- nzchar(trimws(typed)) && !(typed %in% shiny::isolate(gen_scripts()))
    fresh <- ref_submit_script()
    script <- if (edited) typed else fresh
    if (!nzchar(trimws(script))) stop("Nothing to build: every selected group already has a database.", call. = FALSE)
    f <- job_path(deps$project(), "script",
      paste0(shiny::isolate(ref_base()) %||% job_basename("ribopilot_ref"), ".sh"))
    writeLines(script, f)
    try(Sys.chmod(f, "0755"), silent = TRUE)
    f
  }
  shiny::observeEvent(input$ref_save_script, {
    res <- tryCatch(write_ref_script(), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    notify(paste0("Saved job script to ", res, ". Run Build again once the job finishes to wire the databases."))
  })
  shiny::observeEvent(input$ref_submit_job, {
    sched <- detect_scheduler()
    if (is.null(sched)) { notify("No scheduler (sbatch/qsub/bsub) found on PATH.", "warning"); return() }
    # two BUILD runs copy into the same seed FASTA, so a second job for a group
    # whose first job is still out can leave a half-written database behind
    dup <- intersect(build_targets(), shiny::isolate(submitted()))
    if (length(dup) > 0 && !isTRUE(shiny::isolate(resubmit_ok()))) {
      resubmit_ok(TRUE)
      notify(sprintf(paste0("A build job was already submitted for %d of these group(s) this session ",
        "and may still be queued. Click Submit again to submit anyway."), length(dup)), "warning")
      return()
    }
    resubmit_ok(FALSE)
    res <- tryCatch(write_ref_script(), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    # read the targets HERE, same reason as the clade apply above: the later()
    # callback gets a session, not a reactive context
    targets <- build_targets()
    # a busy scheduler can sit on the submit for several seconds, and a dead
    # button reads as a missed click. The popup has no dismiss and covers the
    # modal, so the second click that would submit the build twice cannot land.
    shinyWidgets::sendSweetAlert(title = "Submitting...",
      text = paste0("Handing the job script to ", sched$cmd,
        ". Hold tight, a busy cluster can take a few seconds to answer."),
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        out <- tryCatch(system2(sched$cmd, shQuote(res), stdout = TRUE, stderr = TRUE),
          error = function(e) conditionMessage(e))
        submitted(unique(c(shiny::isolate(submitted()), targets)))
        # record the build so the poller wires it when the job's FASTAs appear
        try(with_con(function(con) {
          mark_reference_submitted(con,
            reference_missing_dbs(ref_config(con)$ref_dir, targets))
        }), silent = TRUE)
        headless_open(FALSE)
        shiny::removeModal()
        bump()
        notify(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " "),
          ". The database(s) will be wired automatically once the job finishes."))
      })
    }, delay = 0.1)
  })

  idle_buttons <- function() {
    shinyjs::hide("ref_nf_stop"); shinyjs::show("ref_nf_start")
    shinyjs::addClass(id = "ref_gears", class = "paused")
  }

  shiny::observeEvent(input$ref_nf_start, {
    ids <- build_targets()
    if (length(ids) == 0) { notify("Select the group(s) to build.", "warning"); return() }
    ref_state(nf_progress_empty())
    shinyjs::hide("ref_nf_start"); shinyjs::show("ref_nf_stop"); shinyjs::show("ref_nf_progress")
    shinyjs::removeClass(id = "ref_gears", class = "paused")
    # the process object is kept so its output can be streamed; the pid handed
    # back is what the reference rows store
    launcher <- function(cmd, args, log) {
      nf_env <- c("current", COLUMNS = "500", NXF_ANSI_SUMMARY = "true")
      pin <- nf_pin_version()
      if (!is.na(pin)) nf_env <- c(nf_env, NXF_VER = pin)
      p <- nf_spawn(cmd, c(args, "-ansi-log"), deps$project(), nf_env)
      ref_proc(p); ref_log(log)
      p$get_pid()
    }
    res <- tryCatch(with_con(function(con) start_reference_build(con, deps$project(), ids,
      api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
      resume = isTRUE(input$ref_nf_resume %||% TRUE), launcher = launcher)),
      error = function(e) e)
    if (inherits(res, "error")) {
      # the launcher may already have spawned Nextflow; nothing stamped the rows,
      # so leaving it running would produce results the poll never wires
      kill_ref_nf(); idle_buttons(); notify(conditionMessage(res), "error"); return()
    }
    bump()
    if (is.null(res)) {
      idle_buttons()
      notify("Every selected group already had a database; nothing to run.")
    }
  })

  shiny::observe({
    shiny::req(ref_proc())
    shiny::invalidateLater(250)
    p <- ref_proc()
    grab <- function() {
      out <- tryCatch(p$read_output_lines(), error = function(e) character(0))
      err <- tryCatch(p$read_error_lines(), error = function(e) character(0))
      if (length(c(out, err))) {
        ref_state(nf_progress_parse(ref_state(), c(out, err)))
        # the rows record a log path, so the piped stream is teed into it
        append_run_log(ref_log(), c(out, err))
      }
    }
    if (p$is_alive()) {
      grab()
    } else {
      grab()
      ref_proc(NULL)
      idle_buttons()
      poll_now()
      bump()
      notify("Reference build run finished.")
    }
  })

  output$ref_nf_prog_all <- shiny::renderText(nf_progress_text(ref_state()))

  kill_ref_nf <- function() {
    p <- ref_proc()
    if (!is.null(p) && p$is_alive()) p$kill()
    ref_proc(NULL)
  }
  # don't leave an orphaned Nextflow process if the browser closes mid-run
  session$onSessionEnded(function() {
    p <- shiny::isolate(ref_proc())
    if (!is.null(p) && p$is_alive()) try(p$kill(), silent = TRUE)
  })
  shiny::observeEvent(input$ref_nf_stop, {
    kill_ref_nf(); idle_buttons(); poll_now(); bump()
  })
  shiny::observeEvent(input$ref_nf_close, {
    kill_ref_nf(); headless_open(FALSE); shiny::removeModal(); poll_now(); bump()
  })

  # ---- log tail ----
  shiny::observeEvent(input$ref_open_log, {
    i <- suppressWarnings(as.integer(input$ref_open_log))
    d <- deps$table_data()
    if (is.na(i) || i < 1 || i > nrow(d)) return()
    ref_id <- d$ID[i]
    row <- with_con(function(con) DBI::dbGetQuery(con,
      "SELECT log_file, fail_reason FROM reference WHERE ref_id = ?", params = list(ref_id)))
    log <- row$log_file[1]
    txt <- if (is.na(log) || !file.exists(log)) "No log for this group yet." else
      paste(utils::tail(readLines(log, warn = FALSE), 200), collapse = "\n")
    shiny::showModal(shiny::modalDialog(title = paste("Build log:", ref_id), size = "l",
      easyClose = TRUE,
      shiny::tags$pre(style = "max-height:60vh; overflow:auto; font-size:12px;", txt),
      footer = shiny::modalButton("Close")))
  })

  # ---- settings ----
  shiny::observeEvent(input$ref_dir_btn, {
    res <- tryCatch(with_con(function(con)
        reference_settings_modal(ref_config(con), blastdb_provenance(con))),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  })
  shiny::observeEvent(input$ref_dir_apply, {
    v <- validate_ref_dir(input$ref_dir_path %||% "")
    if (!isTRUE(v$ok)) { notify(v$msg, "error"); return() }
    b <- validate_blastdb_prefix(input$ref_blastdb_path %||% "")
    if (!isTRUE(b$ok)) { notify(b$msg, "error"); return() }
    with_con(function(con) set_ref_config(con, ref_dir = v$msg,
      min_records = input$ref_min, max_records = input$ref_max, blastdb = b$msg))
    shiny::removeModal(); bump()
    notify(paste("Reference folder set to", v$msg))
  })

  # apply_clade is exposed so a test can run it the way later() does: session
  # domain, no reactive context; the row actions are exposed so tests can drive
  # them without a button id, which is what the UI does now too
  list(open_ladder = open_ladder, apply_clade = apply_clade, poll = poll_now,
    do_build = do_build, do_skip = do_skip, open_existing = open_existing,
    open_move = open_move, do_keep_together = do_keep_together, do_dissolve = do_dissolve)
}
