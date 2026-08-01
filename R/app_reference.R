# The Reference tab: one row per clade group, a ladder to change the clade, and
# the build / reuse / skip actions.
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
    "SELECT reference_switch, ref_id AS ID, clade_name, clade_rank, n_samples, n_records,
       COALESCE(source, '') AS source, COALESCE(seeds_db, '') AS database,
       COALESCE(fail_reason, '') AS fail_reason, time_stamp
     FROM reference ORDER BY ref_id"), error = function(e) NULL)
  if (is.null(d)) {
    d <- data.frame(reference_switch = integer(0), ID = character(0), clade_name = character(0),
      clade_rank = character(0), n_samples = integer(0), n_records = integer(0),
      source = character(0), database = character(0), fail_reason = character(0),
      time_stamp = integer(0), stringsAsFactors = FALSE)
  }
  d
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
#' @noRd
reference_ladder_modal <- function(ref_id, clade_name, ladder, chosen_taxid, n_group) {
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
    shiny::uiOutput("ref_search_results"),
    footer = shiny::tagList(
      shiny::actionButton("ref_clade_apply", "Use this clade", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Build modal: the options a build run takes
#' @param ref_ids Groups that will be built. @param ref_dir Destination folder.
#' @noRd
reference_build_modal <- function(ref_ids, ref_dir) {
  shiny::modalDialog(
    title = sprintf("Build %d reference database(s)", length(ref_ids)), size = "m",
    shiny::tags$p(class = "text-muted", style = "font-size:12px;",
      "Downloads GenBank records and writes seed + label FASTAs to:"),
    shiny::tags$pre(style = "font-size:12px;", ref_dir),
    shiny::tags$pre(style = "max-height:120px; overflow:auto; font-size:12px;",
      paste(ref_ids, collapse = "\n")),
    shiny::textInput("ref_api_key", "NCBI API key (optional):", width = "100%",
      value = Sys.getenv("NCBI_API_KEY", "")),
    shiny::textInput("ref_conda", "Conda environment providing Entrez Direct (optional):",
      width = "100%", value = ""),
    shinyWidgets::prettyCheckbox("ref_refseq", "RefSeq records only", value = FALSE,
      status = "primary", shape = "curve"),
    shinyWidgets::prettyCheckbox("ref_overwrite", "Rebuild even if the database exists",
      value = FALSE, status = "primary", shape = "curve"),
    opts_help("Builds run in the background; the tab shows their state and log."),
    footer = shiny::tagList(
      shiny::actionButton("ref_build_go", "Start build", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Use-an-existing-database modal: manifest entries plus typed paths
#' @param manifest From [read_ref_manifest()]. @param ref_id Group key.
#' @noRd
reference_existing_modal <- function(manifest, ref_id) {
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
    title = sprintf("Use an existing database for '%s'", ref_id), size = "l",
    listing,
    shiny::textInput("ref_existing_seeds", "Seed FASTA path:", width = "100%"),
    shiny::textInput("ref_existing_labels", "Label FASTA path:", width = "100%"),
    opts_help("Choosing a listed database fills these paths; you can also type ",
      "paths to a database built somewhere else."),
    footer = shiny::tagList(
      shiny::actionButton("ref_existing_apply", "Use this database", class = "btn-primary btn-sm"),
      shiny::modalButton("Cancel")))
}

#' Wire the Reference tab
#'
#' @param deps list(con = fn -> DBI con, project = fn -> chr, refresh = reactiveVal,
#'   selected_ids = fn -> chr, table_data = reactive -> data.frame).
#' @return list(open_ladder = function(ref_id), poll = function()).
#' @noRd
reference_handlers <- function(input, output, session, deps) {
  ladder_target <- shiny::reactiveVal(NULL)
  existing_target <- shiny::reactiveVal(NULL)
  build_targets <- shiny::reactiveVal(character(0))
  picked <- shiny::reactiveVal(NULL)   # taxid chosen from the search box

  with_con <- function(fn) {
    con <- deps$con(); on.exit(DBI::dbDisconnect(con)); fn(con)
  }
  bump <- function() deps$refresh(shiny::isolate(deps$refresh()) + 1)
  notify <- function(msg, type = "message") shiny::showNotification(msg, type = type, duration = 8)

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
          n <- with_con(function(con) {
            cfg <- ref_config(con)
            l <- sample_ladders(con)
            if (nrow(l) == 0) stop("No samples have a resolved lineage yet. Assign taxonomy first.",
              call. = FALSE)
            cf <- cached_count_fn(cfg$ref_dir,
              api_key = Sys.getenv("NCBI_API_KEY", ""), conda_env = NULL)
            s <- suggest_reference_groups(l, cf, cfg$min_records, cfg$max_records)
            if (nrow(s$groups) == 0) stop("Could not count GenBank records for any clade (is Entrez Direct on PATH?).",
              call. = FALSE)
            write_reference_groups(con, s$groups, s$members)
            nrow(s$groups)
          })
          bump()
          shinyWidgets::sendSweetAlert(title = "Clades suggested",
            text = sprintf("Proposed %d reference group(s). Review, then build or skip.", n),
            type = "success")
          NULL
        }, error = function(e) e)
        if (inherits(outcome, "error")) {
          shinyWidgets::sendSweetAlert(title = "Suggestion failed",
            text = conditionMessage(outcome), type = "error")
        }
      })
    }, delay = 0.1)
  })

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
      l <- group_ladder(sample_ladders(con), ids,
        count_fn = cached_count_fn(cfg$ref_dir, api_key = Sys.getenv("NCBI_API_KEY", "")))
      reference_ladder_modal(ref_id, row$clade_name[1], l, row$clade_taxid[1], length(ids))
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  }
  shiny::observeEvent(input$ref_edit_btn, open_ladder(deps$selected_ids()))
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

  # Slug a leftover group's own clade name into a ref_id that cannot collide
  # with the group being edited (or any other existing group): the derived
  # slug for a clade re-used by name, e.g. the parent clade the split group
  # came from, can otherwise coincide with the literal ref_id string of the
  # row this same apply just rewrote, silently overwriting it back.
  unique_split_ref_id <- function(con, name, taxid, avoid) {
    slug <- ref_id_from_name(name)
    existing <- tryCatch(DBI::dbGetQuery(con, "SELECT ref_id FROM reference")$ref_id,
      error = function(e) character(0))
    if (!(slug %in% existing) && !(slug %in% avoid)) return(slug)
    paste(slug, taxid, sep = "_")
  }

  shiny::observeEvent(input$ref_clade_apply, {
    ref_id <- ladder_target()
    tid <- suppressWarnings(as.integer(picked() %||% input$ref_rung %||% NA))
    if (is.na(tid)) { notify("Pick a clade first.", "warning"); return() }
    res <- tryCatch(with_con(function(con) {
      ids <- DBI::dbGetQuery(con, "SELECT ID FROM samples WHERE ref_id = ?",
        params = list(ref_id))$ID
      ladders <- sample_ladders(con)
      covered <- intersect(ladder_coverage(ladders, tid), ids)
      # a clade that covers none of the group's samples must not be applied to
      # all of them anyway -- that would silently mis-clade every sample
      if (length(covered) == 0) {
        stop("That clade does not cover any sample in this group; pick a shared rung.", call. = FALSE)
      }
      info <- ladders[ladders$taxid == tid, , drop = FALSE]
      if (nrow(info) == 0) {
        info <- with_taxdump_data(project_taxdump_dir(deps$project()), function(dd) {
          i <- tk_taxid_info(tid, dd)
          data.frame(name = i$sci_name[1], rank = i$rank[1], stringsAsFactors = FALSE)
        })
      }
      cfg <- ref_config(con)
      n <- cached_count_fn(cfg$ref_dir, api_key = Sys.getenv("NCBI_API_KEY", ""))(tid, info$name[1])
      set_reference_clade(con, ref_id = ref_id, clade_taxid = tid, clade_name = info$name[1],
        clade_rank = info$rank[1], n_records = n, member_ids = covered)

      # samples the chosen clade does not cover are not left inside a clade
      # that does not contain them: split them into a group of their own,
      # keyed on their own shared lineage (deepest rung all of them pass through)
      leftover <- setdiff(ids, covered)
      left_name <- NULL
      if (length(leftover) > 0) {
        gl <- group_ladder(ladders, leftover)
        if (nrow(gl) > 0) {
          lg <- gl[order(-gl$level), , drop = FALSE][1, ]
          left_key <- unique_split_ref_id(con, lg$name, lg$taxid, avoid = ref_id)
          set_reference_clade(con, ref_id = left_key, clade_taxid = lg$taxid, clade_name = lg$name,
            clade_rank = lg$rank, n_records = NA_integer_, member_ids = leftover)
          left_name <- lg$name
        }
      }
      list(name = info$name[1], left_n = length(leftover), left_name = left_name)
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal()
    bump()
    if (!is.null(res$left_name)) {
      notify(sprintf("%d sample(s) outside %s were split into a new group (%s); open it to pick a clade.",
        res$left_n, res$name, res$left_name), "warning")
    } else {
      notify(sprintf("Group now uses %s.", res$name))
    }
  })

  # ---- skip ----
  shiny::observeEvent(input$ref_skip_btn, {
    ids <- deps$selected_ids()
    if (length(ids) == 0) { notify("Select a group first.", "warning"); return() }
    res <- tryCatch(with_con(function(con) { for (r in ids) skip_reference_group(con, r); length(ids) }),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    bump()
    notify(sprintf("Skipped %d group(s); their samples keep their current database.", res))
  })

  # ---- use existing ----
  shiny::observeEvent(input$ref_existing_btn, {
    ids <- deps$selected_ids()
    if (length(ids) != 1) { notify("Select exactly one group.", "warning"); return() }
    existing_target(ids[1])
    res <- tryCatch(with_con(function(con) reference_existing_modal(
      read_ref_manifest(ref_config(con)$ref_dir), ids[1])), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  })
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
    res <- tryCatch(with_con(function(con) use_existing_reference(con, existing_target(),
      input$ref_existing_seeds %||% "", input$ref_existing_labels %||% "")),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    notify(sprintf("Pointed %d sample(s) at the existing database.", res))
  })

  # ---- build ----
  shiny::observeEvent(input$ref_build_btn, {
    ids <- deps$selected_ids()
    if (length(ids) == 0) { notify("Select the group(s) to build.", "warning"); return() }
    build_targets(ids)
    res <- tryCatch(with_con(function(con) reference_build_modal(ids, ref_config(con)$ref_dir)),
      error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  })
  shiny::observeEvent(input$ref_build_go, {
    ids <- build_targets()
    res <- tryCatch(with_con(function(con) {
      for (r in ids) {
        start_reference_build(con, deps$project(), r,
          api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
          conda_env = input$ref_conda %||% "", overwrite = isTRUE(input$ref_overwrite))
      }
      length(ids)
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::removeModal(); bump()
    notify(sprintf("Started %d build(s); watch the State column.", res))
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

  # ---- reference folder ----
  shiny::observeEvent(input$ref_dir_btn, {
    res <- tryCatch(with_con(function(con) {
      cfg <- ref_config(con)
      shiny::modalDialog(title = "Reference database folder", size = "l",
        shiny::textInput("ref_dir_path", NULL, width = "100%", value = cfg$ref_dir),
        opts_help("Databases built here are reused by every project pointed at this ",
          "folder, including a shared drive. Changing it does not move databases ",
          "already built, and groups already wired keep working."),
        shiny::numericInput("ref_min", "Minimum GenBank records for a suggested clade:",
          value = cfg$min_records, min = 1, step = 100),
        shiny::numericInput("ref_max", "Maximum GenBank records for a suggested clade:",
          value = cfg$max_records, min = 1, step = 1000),
        footer = shiny::tagList(
          shiny::actionButton("ref_dir_apply", "Save", class = "btn-primary btn-sm"),
          shiny::modalButton("Cancel")))
    }), error = function(e) e)
    if (inherits(res, "error")) { notify(conditionMessage(res), "error"); return() }
    shiny::showModal(res)
  })
  shiny::observeEvent(input$ref_dir_apply, {
    v <- validate_ref_dir(input$ref_dir_path %||% "")
    if (!isTRUE(v$ok)) { notify(v$msg, "error"); return() }
    with_con(function(con) set_ref_config(con, ref_dir = v$msg,
      min_records = input$ref_min, max_records = input$ref_max))
    shiny::removeModal(); bump()
    notify(paste("Reference folder set to", v$msg))
  })

  list(open_ladder = open_ladder,
    poll = function() {
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
    })
}
