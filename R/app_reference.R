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
  ref_proc <- shiny::reactiveVal(NULL)
  ref_log <- shiny::reactiveVal("")
  ref_state <- shiny::reactiveVal(nf_progress_empty())
  headless_open <- shiny::reactiveVal(FALSE)
  gen_scripts <- shiny::reactiveVal(character(0))   # scripts this session generated
  submitted <- shiny::reactiveVal(character(0))     # groups with a cluster job out
  resubmit_ok <- shiny::reactiveVal(FALSE)

  with_con <- function(fn) {
    con <- deps$con(); on.exit(DBI::dbDisconnect(con)); fn(con)
  }
  bump <- function() deps$refresh(shiny::isolate(deps$refresh()) + 1)
  # read fresh each time so tests can substitute a runner that spawns nothing
  count_runner <- function() getOption("ribopilot.count_runner", nf_run_blocking)
  notify <- function(msg, type = "message") shiny::showNotification(msg, type = type, duration = 8)

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
          n <- with_con(function(con) {
            cfg <- ref_config(con)
            l <- ladders_now(con)
            if (nrow(l) == 0) stop("No samples have a resolved lineage yet. Assign taxonomy first.",
              call. = FALSE)
            # one start rank for both: the walk must only ask about rungs the
            # count run was told to count, or every answer comes back NA
            start_rank <- "family"
            counts <- run_clade_counts(deps$project(), cfg$ref_dir, candidate_rungs(l, start_rank),
              api_key = Sys.getenv("NCBI_API_KEY", ""), runner = count_runner())
            s <- suggest_reference_groups(l, cache_count_fn(cfg$ref_dir),
              cfg$min_records, cfg$max_records, start_rank = start_rank)
            if (nrow(s$groups) == 0) {
              # the counting run is the only thing that could have failed here,
              # so its log is the message, not a pointer to a file
              log <- attr(counts, "log")
              tail <- log_tail(log, 5)
              stop("Could not count GenBank records for any clade (the counting run returned no results).",
                if (!is.null(log)) paste0(" Log: ", log, ".") else "",
                if (nzchar(tail)) paste0(" Last lines: ", tail) else "",
                call. = FALSE)
            }
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
      # cache-only: opening the modal must not shell out. Suggest pre-fills the
      # cache for these rungs; anything else renders as "records unknown"
      l <- group_ladder(ladders_now(con), ids, count_fn = cache_count_fn(cfg$ref_dir))
      reference_ladder_modal(ref_id, row$clade_name[1], l, row$clade_taxid[1], length(ids),
        index_ready = file.exists(tax_index_path(project_taxdump_dir(deps$project()))))
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
        info <- with_taxdump_data(project_taxdump_dir(deps$project()), function(dd) {
          i <- tk_taxid_info(tid, dd)
          data.frame(name = i$sci_name[1], rank = i$rank[1], stringsAsFactors = FALSE)
        })
      }
      cfg <- ref_config(con)
      run_clade_counts(deps$project(), cfg$ref_dir,
        data.frame(ref_id = ref_id_from_name(info$name[1]), taxid = as.integer(tid),
          clade_name = info$name[1], query = clade_query(info$name[1]),
          stringsAsFactors = FALSE),
        api_key = Sys.getenv("NCBI_API_KEY", ""), runner = count_runner())
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

  # ---- build: one Nextflow run for every selected group ----
  # the job script is written now but runs much later, so the manifest it reads
  # has to be on disk before the modal closes
  ref_submit_script <- function() {
    s <- with_con(function(con) {
      cfg <- ref_config(con)
      ids <- reference_missing_dbs(cfg$ref_dir, build_targets())
      if (length(ids) == 0) return("")
      man <- write_build_manifest(deps$project(), reference_rows(con, ids))
      build_reference_submit_script(deps$project(), man, cfg$ref_dir,
        api_key = input$ref_api_key %||% "", refseq = isTRUE(input$ref_refseq),
        resume = isTRUE(input$ref_nf_resume %||% TRUE), sched = detect_scheduler())
    })
    if (nzchar(s)) gen_scripts(unique(c(shiny::isolate(gen_scripts()), s)))
    s
  }

  shiny::observeEvent(input$ref_build_btn, {
    ids <- deps$selected_ids()
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
    script <- if (hl) tryCatch(ref_submit_script(), error = function(e) "") else ""
    headless_open(hl)
    shiny::showModal(reference_run_modal(length(res$todo), res$ref_dir,
      headless = hl, sched = sched, script = script))
  })

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
    f <- file.path(deps$project(), "submit_ribopilot_reference.sh")
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
    out <- tryCatch(system2(sched$cmd, shQuote(res), stdout = TRUE, stderr = TRUE),
      error = function(e) conditionMessage(e))
    submitted(unique(c(shiny::isolate(submitted()), build_targets())))
    headless_open(FALSE)
    shiny::removeModal()
    notify(paste0("Submitted via ", sched$cmd, ": ", paste(out, collapse = " "),
      ". Run Build again once the job finishes to wire the databases."))
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

  # apply_clade is exposed so a test can run it the way later() does: session
  # domain, no reactive context
  list(open_ladder = open_ladder, apply_clade = apply_clade, poll = poll_now)
}
