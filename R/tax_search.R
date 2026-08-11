# Offline organism-name search over a project's pinned NCBI taxdump, powering the
# in-app taxonomy review. Searches are indexed LIKE lookups against the project's
# cached taxonomy index (see R/tax_index.R), so the app's searchable box stays
# responsive without any network calls.

#' Search organism scientific names in a taxdump's cached index
#'
#' Prefix matches rank first (so "Danio" surfaces the genus and its species),
#' then substring matches; results are ordered shortest-name first.
#'
#' @param taxdump_dir Directory with the taxdump (index built on demand).
#' @param query Search string (>= 2 chars).
#' @param limit Max results (default 25).
#' @return data.frame(taxid, name, rank), possibly empty.
#' @noRd
tax_search <- function(taxdump_dir, query, limit = 25) {
  empty <- data.frame(taxid = integer(0), name = character(0), rank = character(0),
    stringsAsFactors = FALSE)
  query <- trimws(query %||% "")
  if (nchar(query) < 2) return(empty)
  idx <- build_tax_index(taxdump_dir)
  if (is.null(idx)) return(empty)
  con <- DBI::dbConnect(RSQLite::SQLite(), idx)
  on.exit(DBI::dbDisconnect(con))
  esc <- function(s) gsub("([%_])", "\\\\\\1", s)  # escape LIKE wildcards
  q <- esc(query)
  res <- DBI::dbGetQuery(con,
    "SELECT n.taxid AS taxid, n.name AS name, t.rank AS rank,
       CASE WHEN n.name LIKE ? ESCAPE '\\' THEN 0 ELSE 1 END AS pref
     FROM tax_names n JOIN tax_nodes t ON t.taxid = n.taxid
     WHERE n.sci = 1 AND n.name LIKE ? ESCAPE '\\'
     ORDER BY pref, length(n.name), n.name
     LIMIT ?",
    params = list(paste0(q, "%"), paste0("%", q, "%"), as.integer(limit)))
  res$pref <- NULL
  res$taxid <- suppressWarnings(as.integer(res$taxid))
  res
}

#' Auto-assign taxonomy to samples from their `Taxon` string (local, no Nextflow)
#'
#' Resolves each sample's `Taxon` against the local taxdump and writes
#' `taxid`/`sci_name`/`rank`/`lineage`/`match_type` into `samples`. This is what
#' the Taxonomy tab's Update button runs.
#'
#' @param con Open DBI connection to the project database.
#' @param ids Sample IDs to (re)resolve.
#' @param taxdump_dir Directory with the project's pinned taxdump.
#' @return (invisibly) number of samples updated.
#' @noRd
resolve_samples_taxonomy <- function(con, ids, taxdump_dir) {
  if (length(ids) == 0) return(invisible(0L))
  taxdump_dir <- path.expand(taxdump_dir)
  ph <- paste(rep("?", length(ids)), collapse = ",")
  rows <- DBI::dbGetQuery(con, sprintf("SELECT ID, Taxon FROM samples WHERE ID IN (%s)", ph),
    params = as.list(ids))
  res <- resolve_taxonomy(rows$Taxon, taxdump_dir)
  m <- res[match(rows$Taxon, res$Taxon), ]
  for (i in seq_len(nrow(rows))) {
    DBI::dbExecute(con,
      "UPDATE samples SET taxid = ?, sci_name = ?, rank = ?, lineage = ?, match_type = ?,
         taxonomy_switch = 2, time_stamp = ? WHERE ID = ?",
      params = list(
        if (is.na(m$taxid[i])) NA_integer_ else as.integer(m$taxid[i]),
        m$sci_name[i], m$rank[i], m$lineage[i], m$match_type[i],
        as.integer(Sys.time()), rows$ID[i]))
  }
  invisible(nrow(rows))
}

#' Assign a chosen taxid to samples (manual override from the review modal)
#' @noRd
assign_samples_taxid <- function(con, ids, taxid, taxdump_dir) {
  if (length(ids) == 0 || is.na(taxid)) return(invisible(0L))
  info <- tk_taxid_info(as.integer(taxid), taxdump_dir)
  if (nrow(info) == 0) stop("Could not resolve taxid ", taxid, call. = FALSE)
  m <- info[1, ]
  for (id in ids) {
    DBI::dbExecute(con,
      "UPDATE samples SET taxid = ?, sci_name = ?, rank = ?, lineage = ?,
         match_type = 'manual', taxonomy_switch = 2, time_stamp = ? WHERE ID = ?",
      params = list(as.integer(taxid), m$sci_name, m$rank, m$lineage,
        as.integer(Sys.time()), id))
  }
  invisible(length(ids))
}

#' Taxonomy review modal: search the local taxdump + assign to selected samples
#' @noRd
taxon_review_modal <- function(samples_df, n) {
  cur <- utils::head(samples_df, 12)
  rowsum <- paste(apply(cur, 1, function(r) {
    sprintf("%s: \"%s\" -> %s", r[["ID"]], r[["Taxon"]],
      if (is.na(r[["sci_name"]]) || !nzchar(r[["sci_name"]] %||% "")) "unassigned"
      else sprintf("%s [%s, %s]", r[["sci_name"]], r[["taxid"]], r[["rank"]] %||% "?"))
  }), collapse = "\n")
  modalDialog(
    title = sprintf("Review taxonomy for %d sample(s)", n), size = "l",
    tags$pre(style = "max-height:150px; overflow:auto; font-size:12px; background:#f6f6f6; padding:8px; border-radius:4px;",
      rowsum),
    textInput("tax_query", "Search NCBI organism name (species or genus):",
      width = "100%", placeholder = "e.g. Danio rerio, or a genus like Upeneus"),
    opts_help("Type at least 2 letters; results come from the project's local NCBI ",
      "taxdump. Pick a species or genus row, then apply."),
    uiOutput("tax_results"),
    footer = tagList(
      actionButton("tax_apply", "Assign to selected sample(s)", class = "btn-primary btn-sm"),
      modalButton("Cancel"))
  )
}

#' Wire the taxonomy review modal (search results + apply). Returns an opener
#' function `open_review(ids)` for the app server to call.
#'
#' @param deps list with: con (fn -> DBI con), refresh (reactiveVal),
#'   taxdump (fn -> taxdump dir), samples_for (fn(ids) -> data.frame).
#' @noRd
taxonomy_handlers <- function(input, output, session, deps) {
  targets <- shiny::reactiveVal(character(0))
  q <- shiny::debounce(shiny::reactive(input$tax_query), 300)

  output$tax_results <- shiny::renderUI({
    query <- trimws(q() %||% "")
    if (nchar(query) < 2) return(shiny::helpText("Type at least 2 letters to search."))
    r <- tax_search(deps$taxdump(), query, 30)
    if (nrow(r) == 0) return(shiny::helpText("No matches in the local taxdump."))
    choices <- stats::setNames(as.character(r$taxid),
      sprintf("%s   [%s]   (taxid %s)", r$name,
        ifelse(is.na(r$rank), "?", r$rank), r$taxid))
    shinyWidgets::prettyRadioButtons("tax_pick", NULL, choices = choices,
      shape = "round", status = "primary")
  })

  shiny::observeEvent(input$tax_apply, {
    tid <- suppressWarnings(as.integer(input$tax_pick %||% NA))
    if (is.na(tid)) {
      shiny::showNotification("Pick an organism from the results first.", type = "warning")
      return()
    }
    ids <- targets()
    td <- deps$taxdump()
    # "hold tight" popup, then defer the (blocking) taxdump lookup so it renders
    shinyWidgets::sendSweetAlert(
      title = "Assigning taxonomy...",
      text = "Hold tight, applying the assignment to the selected sample(s).",
      type = "info", btn_labels = NA, closeOnClickOutside = FALSE)
    later::later(function() {
      shiny::withReactiveDomain(session, {
        res <- tryCatch({
          con <- deps$con()
          n <- assign_samples_taxid(con, ids, tid, td)
          DBI::dbDisconnect(con)
          list(ok = TRUE, n = n)
        }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
        if (isTRUE(res$ok)) {
          deps$refresh(shiny::isolate(deps$refresh()) + 1)
          shiny::removeModal()
          shinyWidgets::sendSweetAlert(title = "Taxonomy assigned",
            text = sprintf("Assigned taxid %s to %d sample(s).", tid, res$n), type = "success")
        } else {
          shinyWidgets::sendSweetAlert(title = "Assignment failed", text = res$msg, type = "error")
        }
      })
    }, delay = 0.1)
  })

  function(ids) {
    if (length(ids) == 0) {
      shiny::showNotification("No samples selected.", type = "warning"); return()
    }
    targets(ids)
    shiny::showModal(taxon_review_modal(deps$samples_for(ids), length(ids)))
  }
}
