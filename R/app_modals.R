# Shiny modal dialogs (Assemble / Annotate) + their plot/table helpers.

#' Coloured yes/no/na status badge (mirrors MitoPilot's details-modal badges).
#'
#' `invert` flips the colour mapping so "yes" reads as a warning (used for
#' Problem). `neutral_no` renders "no" in grey rather than green (used for
#' Partial, where "no" just means complete, not "good").
#' @noRd
status_badge <- function(label, state, invert = FALSE, neutral_no = FALSE) {
  val <- if (is.null(state) || length(state) == 0 || is.na(state)) "na" else as.character(state)
  color_val <- if (neutral_no && val == "no") "na" else val
  tone <- if (color_val == "yes") { if (invert) "warning" else "success" }
    else if (color_val == "no") { if (invert) "success" else "warning" } else "neutral"
  span(class = paste0("mp-pill mp-pill-", tone), paste0(label, ": ", toupper(val)))
}

#' The stored per-contig BLAST hit, labelled with the search that produced it
#'
#' Three searches can write `assemblies.blast_*`: the remote core_nt one, the
#' local one against a configured rDNA reference panel, and the local one
#' against the GetOrganelle seed file. The seed file is the same set of
#' sequences that baited the assembly, so its hit answers "which seed does this
#' contig resemble", not "what is this"; the other two are real reference
#' searches. Shown where no BLAST screen exists, and never without saying
#' which search it came from.
#' @param r One row of `unit_evidence()`.
#' @noRd
seed_hit_block <- function(r) {
  sp <- r$blast_species %||% NA_character_
  if (length(sp) == 0 || is.na(sp) || !nzchar(sp)) {
    return(NULL)
  }
  acc <- r$blast_accession %|NA|% ""
  src <- r$blast_source %|NA|% ""
  # a legacy row has no source recorded; it is a seed hit, so it keeps the caveat
  from_seeds <- identical(src, "seeds") || !nzchar(src)
  from_local <- identical(src, "local_rrna")
  div(style = "margin-top:8px;",
    tags$span(style = "color:#555;", if (from_seeds) {
      "Closest sequence in the local seed file (a hint, NOT a GenBank identification):"
    } else if (from_local) {
      "Best nuclear rDNA reference from the bundled rDNA database:"
    } else {
      "Best nuclear rDNA reference from the NCBI core_nt search:"
    }),
    div(class = "code-block", style = "white-space:normal; margin-top:3px;",
      HTML(sprintf("<b>%s</b> &nbsp; %s &nbsp; %s%% id", sp,
        if (nzchar(acc)) {
          sprintf("<a href='https://www.ncbi.nlm.nih.gov/nuccore/%s' target='_blank'>%s</a>", acc, acc)
        } else "no accession",
        r$blast_pident %|NA|% "?"))),
    if (from_seeds) tags$p(class = "mp-help-text mp-fg-warning", style = "margin-top:3px;",
      paste("These seeds are what GetOrganelle used to bait the assembly, so a high",
        "identity here says the contig looks like its own bait. A missing or",
        "distant relative in the seed file can put this name on the wrong genus.")))
}

#' Sample taxon header line for the detail modals (organism + taxid + rank).
#' @noRd
modal_taxon_header <- function(id) {
  s <- app_query("SELECT Taxon, sci_name, taxid, rank FROM samples WHERE ID = ?", params = list(id))
  if (nrow(s) == 0) {
    return(NULL)
  }
  name <- if (!is.na(s$sci_name[1]) && nzchar(s$sci_name[1])) s$sci_name[1] else s$Taxon[1]
  meta <- c(
    if (!is.na(s$taxid[1]) && nzchar(as.character(s$taxid[1])))
      sprintf("<a href='https://www.ncbi.nlm.nih.gov/datasets/taxonomy/%s/' target='_blank' title='NCBI taxonomy ID'>taxid %s</a>",
        s$taxid[1], s$taxid[1]),
    if (!is.na(s$rank[1]) && nzchar(s$rank[1])) s$rank[1])
  div(style = "margin:-6px 0 8px; color:#333; font-size:0.92em;",
    tags$b(name %||% id),
    if (length(meta) > 0) HTML(paste0(" &middot; ", paste(meta, collapse = " &middot; "))),
    if (!is.na(name) && !is.na(s$Taxon[1]) && !identical(s$Taxon[1], name))
      tags$span(style = "color:#999;", sprintf(" (input: %s)", s$Taxon[1])))
}

#' Why a failed sample failed and what to try next, or nothing
#' @param id Sample ID.
#' @noRd
failure_banner <- function(id) {
  r <- app_query("SELECT assemble_switch, TRIM(COALESCE(assemble_notes, '') || ' ' ||
      COALESCE(find_rrna_notes, '')) AS note FROM assemble WHERE ID = ?", params = list(id))
  if (nrow(r) == 0 || !identical(as.integer(r$assemble_switch[1]), 3L)) {
    return(NULL)
  }
  h <- failure_hint(r$note[1])
  div(class = "mp-fg-danger",
    style = "border:1px solid currentColor; border-radius:4px; padding:8px 12px; margin:0 0 8px;",
    tags$b("This sample failed. "), h$why, tags$br(), tags$b("Next step: "), h$next_step,
    if (nzchar(r$note[1] %|NA|% "")) div(style = "margin-top:4px; color:#555;",
      "Pipeline note: ", r$note[1]))
}

#' Why a contig has no BLAST screen
#' @noRd
no_screen_note <- function(failed, blast_on, running) {
  if (failed) {
    "No BLAST screen for this contig. BLAST did not run because the assembly failed."
  } else if (!blast_on) {
    paste("BLAST is switched off for this sample (Assemble step, BLAST Options).",
      "Switch Run BLAST on there, then press Update to screen this contig.")
  } else if (running) {
    paste("No BLAST screen for this contig yet. A run is going for this project,",
      "and its BLAST search may still be working on it.")
  } else {
    paste("No BLAST screen for this contig: the search found no match, or it could",
      "not finish (for example, NCBI could not be reached).")
  }
}

#' Which database answered a contig's BLAST search, and when
#' @param r One row of `unit_evidence()`.
#' @noRd
searched_in <- function(r) {
  at <- r$searched_at %||% NA
  when <- if (is.na(at)) "" else format(as.POSIXct(at, origin = "1970-01-01"), ", searched %Y-%m-%d")
  db <- r$searched_db %|NA|% ""
  if (identical(r$searched_source, "local")) {
    paste0("the bundled rDNA database ", if (nzchar(db)) db else "(release not recorded)", when)
  } else {
    paste0("NCBI core_nt", when)
  }
}

#' Prev / picker / next sample navigation for the detail modals
#'
#' The order is whatever the main table shows when the modal opens: its Lock
#' and State filters, column filters, search, and sort. The script snapshots it
#' into input$nav_order; the server renders the controls.
#' @noRd
rp_sample_nav <- function() {
  tagList(
    uiOutput("sample_nav", inline = TRUE),
    tags$script(HTML(
      "setTimeout(function(){
         var rows = window.rpVisibleRows ? rpVisibleRows() : [];
         var col = function(k){ return rows.map(function(r){ return r[k] == null ? '1' : String(r[k]); }); };
         Shiny.setInputValue('nav_order', {ID: col('ID'), path: col('path'), scaffold: col('scaffold')}, {priority: 'event'});
       }, 0);"))
  )
}

#' Chimera check panel body, shared by the assemble and annotate modals
#'
#' Only one modal is open at a time, so both can render the same output ids;
#' the server picks the unit from whichever modal is open.
#' @noRd
chimera_panel <- function() {
  tagList(
    tags$p(style = "margin:2px 0 4px; font-weight:600;", uiOutput("chimera_result", inline = TRUE)),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      "Looks for pieces of different origin joined in one contig, using the genes' BLAST hits, ",
      "500 bp BLAST windows, ribotyper, and read depth and clipping. Severe (different phyla or ",
      "kingdoms) is flagged at export; minor is only reported. Unresolved windows match several ",
      "classes equally and do not count."),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;", textOutput("chimera_db_source", inline = TRUE)),
    plotOutput("chimera_plot", height = "300px"),
    reactable::reactableOutput("chimera_segments"),
    div(class = "rp-sub-list",
      htmltools::tagAppendAttributes(class = "rp-subsection", rp_details("Possible junctions",
        tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
          "Where the checks think two pieces may be joined, which checks point there, and why. ",
          "Positions match the plot above."),
        uiOutput("chimera_junctions"))),
      htmltools::tagAppendAttributes(class = "rp-subsection", rp_details("All windows",
        tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
          "One row per window: its best BLAST hit in the bundled rDNA database, the hit identity, and the read evidence over that stretch."),
        reactable::reactableOutput("chimera_windows"), open = FALSE))))
}

#' Assemble details modal: contigs, their evidence, and keep/drop
#' @noRd
assemble_modal <- function(id) {
  modalDialog(
    title = div(class = "mp-sticky-head", paste0("Assemble - ", id)),
    size = "l", easyClose = TRUE,
    footer = div(class = "mp-sticky-foot", rp_sample_nav(), mp_footer(dismiss = "Close")),
    modal_taxon_header(id),
    failure_banner(id),
    h5("Contigs", style = "margin:2px 0 4px;"),
    tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
      "Click a contig row to see its BLAST evidence, read coverage, and chimera check below."),
    reactable::reactableOutput("assemble_scaffolds"),
    div(style = "margin-top:6px;",
      actionButton("copy_fasta", "Copy FASTA", icon = icon("clipboard"), class = "btn-sm"),
      tags$span(class = "text-muted", style = "margin-left:8px; font-size:0.82em;",
        "Copies the ticked contigs as FASTA, headed by sample, taxon and contig.")),
    tags$hr(style = "margin:10px 0;"),

    # The evidence that settles the keep/drop call, at the point it is made.
    # blast_screen and blast_ref_candidates are both written per unit by WF1 and
    # neither was rendered anywhere: when they name different organisms, that is
    # the contamination signal.
    # All start collapsed, and collapse again when another scaffold is picked
    # (examine_unit), so the modal opens short for every scaffold.
    div(id = "asm_unit_panels",
      rp_details(tagList("BLAST evidence for ", textOutput("hl_scaf_ev", inline = TRUE)), open = FALSE,
        uiOutput("unit_screen"),
        div(style = "margin-top:8px; color:#555;",
          "Best matches among records whose title says nuclear rDNA:"),
        reactable::reactableOutput("unit_candidates")),

      seqview_ui("seqv_asm", tagList("Sequence and read coverage for ", textOutput("hl_scaf_cov", inline = TRUE)),
        open = FALSE, pick = FALSE),

      rp_details(tagList("Chimera check for ", textOutput("hl_scaf_chim", inline = TRUE)), open = FALSE,
        chimera_panel())),

    # Only a supplied assembly has these: the contigs the rDNA search found but
    # did not carry forward. They never become assembly units, so this is the
    # one place they are visible at all.
    rejected_contigs_section(id),

    rp_details("Notes", open = FALSE,
      textAreaInput("assemble_user_notes", label = NULL, width = "100%",
        value = app_query("SELECT user_notes FROM assemble WHERE ID = ?", params = list(id))$user_notes[1] %|NA|% "",
        placeholder = "Notes for this sample"),
      tags$p(class = "mp-help-text", "Saved as you type. Deleting the text does not clear the stored note; use Clear notes."),
      actionButton("assemble_user_notes_clear", "Clear notes", class = "btn-sm"))
  )
}

#' The contigs the rDNA search left behind, or nothing at all
#'
#' Rendered inline rather than as a server output: the rows never change while
#' the modal is open, and a sample that dropped nothing should not show an empty
#' section.
#' @param id Sample ID.
#' @noRd
rejected_contigs_section <- function(id) {
  d <- rejected_rrna_contigs(id)
  if (nrow(d) == 0) {
    return(NULL)
  }
  names(d) <- c("Contig", "Length", "rRNA found", "Why it was left out")
  rp_details(sprintf("Contigs not carried forward (%d)", nrow(d)), open = FALSE,
    tags$p(class = "mp-help-text", style = "margin-bottom:6px;",
      paste("These contigs in the assembly you supplied had rRNA hits but did",
        "not go on through the pipeline. A contig needs an 18S or 28S anchor to",
        "be carried forward; 5.8S on its own is not enough. Contigs with no rRNA",
        "hit at all are not listed.")),
    reactable::reactable(d, compact = TRUE, sortable = TRUE, wrap = FALSE,
      defaultPageSize = 10,
      columns = list(
        Contig = reactable::colDef(minWidth = 160),
        Length = reactable::colDef(width = 90),
        `rRNA found` = reactable::colDef(width = 110),
        `Why it was left out` = reactable::colDef(minWidth = 200))))
}

#' Warn that exporting will replace an export folder that already exists
#'
#' The group folder is rebuilt from scratch on every export, so a previous
#' bundle (and anything filed alongside it) is deleted. Cheap to redo, painful
#' to discover afterwards.
#' @param dir The folder that would be replaced. @param n Files it holds.
#' @noRd
export_overwrite_modal <- function(dir, n) {
  modalDialog(
    title = mp_modal_title("This export already exists"),
    size = "m", easyClose = FALSE,
    div(style = "font-size:1.02em; line-height:1.5;",
      tags$p(sprintf("Exporting replaces this folder and the %d file(s) in it:", n)),
      div(class = "code-block", style = "white-space:normal; word-break:break-all;", dir),
      tags$p(style = "color:#666;",
        "Rename or move it first if you need to keep the earlier bundle.")),
    footer = mp_footer(dismiss = NULL,
      extra = actionButton("export_overwrite_back", "Cancel"),
      danger = actionButton("export_overwrite_go", "Replace it"))
  )
}

#' Plain words for the tags of submission_risk strings, one string per record
#' @noRd
risk_reasons <- function(risk) {
  vapply(risk, function(r) {
    r <- r %|NA|% ""
    arg <- function(tag) regmatches(r, regexec(paste0(tag, "\\(([^)]*)\\)"), r))[[1]][2]
    tx <- arg("off_target")
    nm <- arg("taxon_differs")
    paste(c(if (!is.na(tx)) sprintf("BLAST hits agree on %s", tx),
      if (grepl("(^|;)contaminant(;|$)", r)) "the best BLAST hit is from another phylum",
      if (grepl("(^|;)chimeric\\(ITSx\\)", r)) "ITSx, the tool NCBI runs on submissions, calls it chimeric"
      else if (grepl("(^|;)chimeric", r)) "parts of the contig match organisms from different phyla or kingdoms",
      if (grepl("(^|;)gap_in_rDNA", r)) "an assembly gap inside the rDNA",
      if (grepl("(^|;)misassembly", r)) "features overlap or are out of order",
      if (!is.na(nm)) sprintf("your Taxon differs from %s", nm),
      if (grepl("(^|;)missing_organism", r)) "no organism name"),
      collapse = "; ")
  }, character(1), USE.NAMES = FALSE)
}

#' Ask before exporting units the contamination checks flagged
#'
#' `flagged` is rows of `export_summary_df()` whose submission_risk carries an
#' off_target or contaminant tag and that no curator has marked Reviewed.
#' @noRd
export_screen_modal <- function(flagged) {
  why <- risk_reasons(flagged$submission_risk)
  items <- lapply(seq_len(nrow(flagged)), function(i) tags$li(
    tags$b(sprintf("%s p%s.s%s", flagged$ID[i], flagged$path[i], flagged$scaffold[i])),
    sprintf(" (%s): %s.", flagged$organism[i] %|NA|% "no Taxon", why[i])))
  modalDialog(
    title = mp_modal_title("Possible contamination or chimera"),
    size = "m", easyClose = FALSE,
    div(style = "font-size:1.02em; line-height:1.5;",
      tags$p(sprintf("%s in this group may contain another organism's rDNA, exported under the sample's name:",
        mp_n(nrow(flagged), "contig"))),
      tags$ul(items),
      tags$p(style = "color:#666;",
        "Leave them out unless you have checked them. Marking a contig Reviewed in",
        "Annotation Review means it will not be asked about again.")),
    footer = mp_footer(dismiss = NULL,
      extra = actionButton("export_screen_back", "Cancel"),
      danger = actionButton("export_screen_keep", "Export anyway"),
      primary = actionButton("export_screen_skip", "Leave them out"))
  )
}

#' One rule's three-way action control
#'
#' Off is the disable, so there is no separate tickbox. Drop is offered only on
#' rules that combine both evidence channels; on the others it is absent and the
#' label carries the reason, because a control that looks available and refuses
#' to act is worse than one that is not there.
#' @noRd
triage_rule_control <- function(def, current) {
  choices <- c("Off" = "off", "Flag" = "flag")
  if ("drop" %in% def$actions) choices <- c(choices, "Drop" = "drop")
  shinyWidgets::radioGroupButtons(
    inputId = paste0("triage_rule_", def$key),
    label = NULL, choices = choices, selected = current,
    size = "xs", individual = FALSE)
}

#' Bulk triage modal: re-run automatic triage over the selected samples
#'
#' A bulk write across a project is the highest-blast-radius action in the app,
#' so it always shows what it would do first, per rule, including how many
#' samples would end up with NOTHING kept. Never applied without the counts.
#' @param n_samples Number of samples the action would cover.
#' @param actions Named list key -> action, from `triage_opts_for()`.
#' @noRd
triage_modal <- function(n_samples, actions = triage_actions_default()) {
  rows <- lapply(TRIAGE_RULE_DEFS, function(def) {
    tags$tr(
      tags$td(style = "padding:3px 8px 3px 0; vertical-align:top;",
        tags$div(style = "font-weight:600;", def$label),
        tags$p(class = "mp-help-text", def$meaning)),
      tags$td(style = "padding:3px 8px; vertical-align:top; white-space:nowrap;",
        triage_rule_control(def, actions[[def$key]])),
      tags$td(style = "padding:3px 8px; vertical-align:top; text-align:right;",
        uiOutput(paste0("triage_n_", def$key), inline = TRUE)))
  })
  modalDialog(
    title = paste("Triage contigs for", mp_n(n_samples, "sample")),
    size = "l", easyClose = TRUE,
    tags$p(class = "mp-help-text",
      paste("Proposes keep / drop / undecided per contig from the evidence already",
        "in the database. Depth is shown in the contig table but never decides:",
        "high depth is often k-mer collapse across the conserved core or IGS",
        "subrepeats, and low depth is often a repeat inside the contig, where reads",
        "map to two places and are left out.")),
    tags$table(style = "width:100%; font-size:0.9em;",
      tags$thead(tags$tr(
        tags$th(style = "text-align:left;", "Rule"),
        tags$th(style = "text-align:left;", "Action"),
        tags$th(style = "text-align:right;", "Decided / Matching"))),
      tags$tbody(rows)),
    tags$p(class = "mp-help-text", style = "margin-top:6px;",
      paste("Rules are applied top to bottom and the first one that matches decides",
        "whether to drop or flag the contig. \"Decided\" counts contigs this rule",
        "will drop or flag. \"Matching\" counts every contig that meets the rule,",
        "even if a higher-priority rule acts on it first.")),
    div(style = "margin:0.4em 0;",
      shinyWidgets::prettyCheckbox("triage_overwrite_user",
        "Also re-decide contigs I set by hand", value = FALSE, status = "danger")),
    tags$hr(style = "margin:8px 0;"),
    uiOutput("triage_preview"),
    footer = mp_footer(danger = actionButton("triage_apply", "Apply"))
  )
}

#' 5'/3' nucleotide-boundary nudge control (mirrors MitoPilot rRNA nudge)
#' @noRd
nt_nudge <- function(id, sym, label) {
  actionButton(id, sym, class = "btn-sm", title = label, `aria-label` = label,
    style = "min-width:34px; padding:1px 8px; font-weight:bold;")
}

#' A collapsible <details> section (open by default unless open = FALSE).
#' @noRd
rp_details <- function(summary, ..., open = TRUE) {
  d <- tags$details(tags$summary(tags$b(summary)), div(style = "padding:6px 2px 2px;", ...))
  if (isTRUE(open)) d$attribs$open <- NA
  d
}

annotate_modal <- function(id, path = 1L, scaffold = 1L) {
  # Always name the unit, even when the sample has only one: every edit in here
  # (features, trim, marks, lock) applies to one (ID, path, scaffold), and the
  # title is the only thing on screen that says which.
  unit_lbl <- sprintf(" [path %s / contig %s]", path, scaffold)
  # Seed the notes box from the DB. Re-opening the same unit assigns an
  # identical() value to sel_unit, which does not invalidate the reactiveVal, so
  # the server-side reset observer never fires; a value-less textAreaInput then
  # reports "" and the write observer persists it over the stored note.
  row <- tryCatch(app_query(
    "SELECT annotate_notes, annotate_lock FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold))),
    error = function(e) data.frame())
  notes <- if (nrow(row)) row$annotate_notes[1] else NA_character_
  # A locked unit is read only here, as in MitoPilot's details window; the
  # server refuses the same edits (lock_blocked() in app_server.R).
  locked <- isTRUE(row$annotate_lock[1] == 1)
  off <- function(x) if (locked) shinyjs::disabled(x) else x
  modalDialog(
    # Flags sit in the pinned header beside the badges they mirror, as in
    # MitoPilot, so they stay in reach while the body scrolls.
    title = div(class = "mp-sticky-head",
      paste0("Annotation Review - ", id, unit_lbl),
      div(class = "mp-sticky-row",
        uiOutput("modal_status_badges", inline = TRUE, style = "display: contents;"),
        uiOutput("modal_status_toggles", inline = TRUE, style = "display: contents;"))),
    size = "l", easyClose = FALSE,
    # Close is a server button so unsaved boundary edits can stop it.
    footer = div(class = "mp-sticky-foot", rp_sample_nav(), if (locked) {
      mp_footer(primary = actionButton("annotate_close", "Close"), dismiss = NULL)
    } else {
      mp_footer(primary = actionButton("lock_and_close", "Lock & Close"), dismiss = NULL,
        extra = tagList(actionButton("annotate_close", "Close"), actionButton("modal_lock", "Lock")))
    }),
    modal_taxon_header(id),
    if (locked) div(class = "alert alert-warning", role = "status", icon("lock"), " ", MP_UNIT_LOCKED),
    # ---- gene list (drives the active feature for the editor + alignment) ----
    h5("Feature model", style = "margin:2px 0 4px;"),
    tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
      "Click a feature to make it active for editing + reference alignment below."),
    reactable::reactableOutput("feature_tbl"),
    div(style = "margin:6px 0; display:flex; gap:6px;",
      off(actionButton("trim_btn", "Trim assembly...", icon = icon("scissors"), class = "btn-sm")),
      off(actionButton("restore_btn", "Restore assembly", icon = icon("rotate-left"), class = "btn-sm"))),
    tags$hr(style = "margin:10px 0;"),

    # ---- collapsible detail panels ----
    # Coverage, boundary editing, and the alignment are the review itself, so
    # they open by default; the Assemble window opens short because its panels
    # are per-scaffold evidence the reviewer picks through.
    seqview_ui("seqv_ann", "Sequence and read coverage"),

    rp_details("Edit active feature boundaries", open = FALSE,
      # normalize control heights so the step box, +/- buttons, and Save/Reset
      # sit on one baseline
      tags$style(HTML(paste0(
        "#feat_edit_row .form-group{margin-bottom:0;} ",
        "#feat_edit_row .control-label{margin-bottom:3px; font-size:0.85em;} ",
        "#feat_edit_row .pretty{margin:0;}"))),
      div(id = "feat_edit_row",
        style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
        div(style = "align-self:flex-end; padding-bottom:6px; font-weight:bold;",
          textOutput("feat_edit_active", inline = TRUE)),
        div(style = "align-self:flex-end; padding-bottom:6px;",
          off(shinyWidgets::prettyCheckbox("feat_edit", "Edit", value = FALSE,
            status = "primary", inline = TRUE))),
        shinyjs::hidden(div(id = "feat_edit_ctrls",
          style = "display:flex; flex-flow:row wrap; gap:1.2em; align-items:flex-end;",
          div(class = "rp-step-box", style = "width:70px;",
            numericInput("feat_step", "Step (nt)", value = 1, min = 1, max = 50, step = 1, width = "70px")),
          div(tags$label("5' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_5_out", "+", "Extend 5' end"),
              nt_nudge("nt_5_in", "-", "Trim 5' end"))),
          div(tags$label("3' end", class = "control-label"),
            div(style = "display:flex; gap:4px;", nt_nudge("nt_3_in", "-", "Trim 3' end"),
              nt_nudge("nt_3_out", "+", "Extend 3' end"))),
          uiOutput("partial_ctrls"),
          div(style = "display:flex; gap:6px; align-items:flex-end;",
            actionButton("feat_save", "Save", icon = icon("floppy-disk"), class = "btn-sm btn-success"),
            actionButton("feat_reset", "Reset", icon = icon("rotate-left"), class = "btn-sm"))))
      ),
      tags$p(class = "mp-help-text", style = "margin-top:4px;",
        paste("Check Edit to adjust boundaries; + extends outward, - trims inward.",
          "Partial 5' and 3' show whether an end is partial and set it by hand. Save to keep, Reset to discard.")),
      div(class = "code-block", style = "margin-top:6px; white-space:normal; word-break:break-all;",
        uiOutput("nudge_seq"))),

    rp_details("Reference alignment (active feature)", open = FALSE,
      uiOutput("msa_header"),
      uiOutput("feat_align_ui")),

    rp_details("Curation warnings", open = FALSE,
      uiOutput("curate_notes"),
      reactable::reactableOutput("flag_tbl")),

    rp_details("Chimera check", open = FALSE, chimera_panel()),

    rp_details("Notes", open = FALSE,
      off(textAreaInput("annotate_notes", label = NULL, width = "100%",
        value = (notes %||% NA_character_) %|NA|% "",
        placeholder = "Notes for this contig")),
      tags$p(class = "mp-help-text", "Saved as you type. Deleting the text does not clear the stored note; use Clear notes."),
      off(actionButton("annotate_notes_clear", "Clear notes", class = "btn-sm")))
  )
}

#' Run-modal title: text on the left, an animated gears spinner on the right
#' @noRd
run_modal_title <- function(text, gears_id) {
  div(class = "mp-sticky-head",
    style = "display:flex; justify-content:space-between; align-items:center; height:38px;",
    span(text), span(id = gears_id, class = "gears paused"))
}

#' Update / run modal: launch the Nextflow workflow with live progress
#'
#' Interactive (local browser): "Run from App" spawns Nextflow on this host.
#' Headless (served to a remote browser, e.g. an SSH-tunneled cluster login
#' node): running on the host would burn the login node, so offer an editable
#' cluster job script to Submit as job / Save instead.
#' @noRd
run_modal <- function(n, headless = FALSE, sched = NULL, script = "") {
  wf <- "RiboPilot pipeline"
  cluster_section <- if (headless) {
    tagList(
      hr(style = "margin:10px 0;"),
      h5("Submit to cluster", style = "margin:6px 0 4px;"),
      tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
        if (is.null(sched)) {
          "No scheduler (sbatch, qsub, or bsub) found on PATH. Save the script and submit it from a machine that has one, or reopen the app with headless = FALSE to run here."
        } else {
          sprintf("%s found. Check the queue or partition, account, and walltime in the script, then Submit as job (or Save script). The job log goes to .runs/logs/ and the Nextflow log to .runs/nextflow/ in the project folder. This window does not follow the job: use the reload button above the table to see new results, and Run Reports to see how the job went.", sched$name)
        }),
      tags$textarea(id = "submit_script_text", rows = 11,
        style = "width:100%; font-family:monospace; font-size:0.78em;", script))
  }
  run_progress <- if (!headless) {
    shinyjs::hidden(div(id = "nf_progress",
      hr(style = "margin:12px 0;"), h5("Progress:", style = "margin:0 0 4px;"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit; white-space:pre-wrap;",
          textOutput("nf_prog_all", container = tags$code)))))
  }
  footer <- if (headless) {
    mp_footer(
      primary = if (!is.null(sched)) actionButton("submit_job", "Submit as job"),
      dismiss = NULL,
      extra = tagList(
        actionButton("save_script", "Save script"),
        actionButton("nf_close", "Close")))
  } else {
    mp_footer(
      primary = actionButton("nf_start", "Run from App"), dismiss = NULL,
      extra = tagList(
        shinyjs::hidden(actionButton("nf_stop", "Stop")),
        actionButton("nf_close", "Close")))
  }
  modalDialog(
    title = run_modal_title(sprintf("%s - %d sample(s) or contig(s) to update", wf, n), "rp_gears"),
    size = "l",
    shinyWidgets::prettyCheckbox("nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("nf_cmd_text")),
    run_progress,
    cluster_section,
    footer = footer
  )
}

#' Reference build modal: run the reference workflow with live progress
#'
#' The same shape as `run_modal()` so a database build looks like any other
#' RiboPilot run: command preview, progress board, Run / Stop / Close, and the
#' cluster-submit section when the app is served headless.
#' @param api_key Key to prefill; the project's own key, not just the env var.
#' @noRd
reference_run_modal <- function(n, ref_dir, headless = FALSE, sched = NULL, script = "",
                                api_key = Sys.getenv("NCBI_API_KEY", "")) {
  cluster_section <- if (headless) {
    tagList(
      hr(style = "margin:10px 0;"),
      h5("Submit to cluster", style = "margin:6px 0 4px;"),
      tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
        if (is.null(sched)) {
          "No scheduler (sbatch, qsub, or bsub) found on PATH. Save the script and submit it from a machine that has one, or reopen the app with headless = FALSE to run here."
        } else {
          sprintf("%s found. Check the queue or partition, account, and walltime in the script, then Submit as job (or Save script). The job log goes to .runs/logs/ and the Nextflow log to .runs/nextflow/ in the project folder. This window does not follow the job: use the reload button above the table to see new results, and Run Reports to see how the job went.", sched$name)
        }),
      tags$textarea(id = "ref_submit_script_text", rows = 11,
        style = "width:100%; font-family:monospace; font-size:0.78em;", script))
  }
  run_progress <- if (!headless) {
    shinyjs::hidden(div(id = "ref_nf_progress",
      hr(style = "margin:12px 0;"), h5("Progress:", style = "margin:0 0 4px;"),
      div(class = "code-block", style = "max-height:320px; overflow:auto;",
        tags$pre(style = "margin:0; background:transparent; border:0; color:inherit; white-space:pre-wrap;",
          textOutput("ref_nf_prog_all", container = tags$code)))))
  }
  footer <- if (headless) {
    mp_footer(
      primary = if (!is.null(sched)) actionButton("ref_submit_job", "Submit as job"),
      dismiss = NULL,
      extra = tagList(
        actionButton("ref_save_script", "Save script"),
        actionButton("ref_nf_close", "Close")))
  } else {
    mp_footer(
      primary = actionButton("ref_nf_start", "Run from App"), dismiss = NULL,
      extra = tagList(
        shinyjs::hidden(actionButton("ref_nf_stop", "Stop")),
        actionButton("ref_nf_close", "Close")))
  }
  modalDialog(
    title = run_modal_title(sprintf("Reference databases - %d clade(s) to build", n), "ref_gears"),
    size = "l",
    tags$p(class = "mp-help-text", style = "margin-bottom:6px;",
      sprintf("Databases are written to %s", ref_dir)),
    textInput("ref_api_key", "NCBI API key (optional):", width = "100%",
      value = api_key),
    shinyWidgets::prettyCheckbox("ref_refseq", "RefSeq records only", value = FALSE,
      status = "primary", inline = TRUE),
    shinyWidgets::prettyCheckbox("ref_nf_resume", "Resume previous run?", value = TRUE,
      status = "primary", inline = TRUE),
    h5("Nextflow command:", style = "margin:6px 0 4px;"),
    div(class = "code-block", textOutput("ref_nf_cmd_text")),
    run_progress,
    cluster_section,
    footer = footer
  )
}

#' Trim-assembly modal: cut the contig to the annotated feature span
#' or a custom range (shifts annotations + coverage).
#' @noRd
trim_modal <- function(n, a_from, a_to) {
  modalDialog(
    title = "Trim assembly", size = "l",
    tags$p(class = "mp-help-text",
      sprintf("Current assembly length: %d bp. Annotated span: %d..%d bp.", n, a_from, a_to)),
    shinyWidgets::prettyRadioButtons("trim_mode", "Trim to:",
      choices = c("Annotated feature boundaries" = "annotated", "Annotated + flank" = "flank", "Custom range" = "custom"),
      selected = "annotated", status = "primary", inline = TRUE),
    tags$p(class = "mp-help-text", style = "margin:2px 0;",
      "Drag across the plot to set the kept window (grey = trimmed); drag the box ",
      "edges to adjust. From/To always show the window that will be kept; choose Custom range to edit them."),
    plotOutput("trim_plot", height = "230px",
      brush = brushOpts(id = "trim_brush", direction = "x", resetOnNew = FALSE,
        fill = "#4C9BE8", opacity = 0.15)),
    div(style = "display:flex; gap:1em; margin-top:6px; align-items:flex-end;",
      conditionalPanel("input.trim_mode == 'flank'",
        numericInput("trim_flank", "Flank (bp)", value = 100, min = 1, step = 50)),
      shinyjs::disabled(numericInput("trim_from", "From (bp)", value = a_from, min = 1, max = n)),
      shinyjs::disabled(numericInput("trim_to", "To (bp)", value = a_to, min = 1, max = n))),
    div(style = "color:#c0392b; font-weight:bold; margin:6px 0 14px;", textOutput("trim_flank_hint")),
    textOutput("trim_chimera_hint"),
    opts_help("Trimming cuts the flanks off the contig and shifts every ",
      "feature coordinate + the coverage track to match. Features that fall entirely ",
      "outside the kept range are removed. Use \"Restore assembly\" to undo."),
    footer = mp_footer(
      danger = actionButton("trim_apply", "Trim", icon = icon("scissors")), dismiss = NULL,
      # actionButton (not modalButton) so Cancel reopens the annotate details
      # modal instead of dropping back to the main Annotate table.
      extra = actionButton("trim_cancel", "Cancel"))
  )
}

#' Batch Trim / Restore window; the count and skipped list render into
#' `batch_trim_summary` (see [batch_trim_summary()]).
#' @param restore TRUE for the Restore window.
#' @noRd
batch_trim_modal <- function(restore = FALSE) {
  modalDialog(
    title = if (restore) "Restore trimmed contigs" else "Trim contigs", easyClose = TRUE,
    if (!restore) tagList(
      shinyWidgets::prettyRadioButtons("batch_trim_mode", "Trim to:",
        choices = c("Match gene ends" = "annotated", "Gene ends + flank" = "flank"),
        selected = "annotated", status = "primary", inline = TRUE),
      conditionalPanel("input.batch_trim_mode == 'flank'",
        numericInput("batch_trim_flank", "Flank (bp)", value = 100, min = 1, step = 50))),
    uiOutput("batch_trim_summary"),
    footer = mp_footer(danger = actionButton("batch_trim_go", if (restore) "Restore" else "Trim",
      icon = icon(if (restore) "rotate-left" else "scissors")))
  )
}

#' How many units a batch Trim / Restore acts on, and a scrolling list of the rest
#' @param p [batch_trim_plan()] result. @param restore TRUE for Restore.
#' @noRd
batch_trim_summary <- function(p, restore = FALSE) {
  ok <- sum(is.na(p$reason)); bad <- p[!is.na(p$reason), , drop = FALSE]
  tagList(
    tags$p(paste0(mp_n(ok, "contig"), if (ok == 1) " is" else " are",
      if (restore) " restored to full length and re-curated." else " trimmed and re-curated.")),
    if (nrow(bad)) tagList(
      tags$p(style = "margin:8px 0 4px;", paste0("Skipped (", nrow(bad), "):")),
      div(style ="max-height:200px; overflow-y:auto; border:1px solid #ddd; padding:4px 8px; font-size:0.9em;",
        lapply(sprintf("%s p%s s%s: %s", bad$ID, bad$path, bad$scaffold, bad$reason), div))))
}

#' Assign-export-group modal: name the group for the selected samples
#' @noRd
group_modal <- function() {
  modalDialog(
    title = "Assign export group", size = "s", easyClose = TRUE,
    textInput("group_name", "Group name",
      placeholder = "letters, digits, dot, dash, underscore; e.g. batch_2026_01"),
    footer = mp_footer(primary = actionButton("assign_group_go", "Assign"))
  )
}

#' Export modal: pick which export group to bundle (or all locked samples), then
#' write it to the project's "out/export" folder. Without this, assigned groups
#' were ignored and every bundle collapsed into one undifferentiated "export".
#' @noRd
export_modal <- function(groups, templates = "default", template = "default",
                         fasta_header = "", fasta_header_marker = "",
                         token_groups = NULL, totals = "{}", snap = NULL,
                         nmnh_cols = list(voucher = NA, uri = NA), nmnh_choices = c("(none)" = "")) {
  flag <- function(id, dflt) snap[[id]] %||% dflt
  # status line between the label and the box it describes
  hdr_box <- function(id, label, value, status) {
    htmltools::tagQuery(textAreaInput(id, tagList(label, insert_field_link(id)), value, width = "100%", height = "60px"))$
      find("textarea")$addClass("mp-hl")$addAttrs(`aria-describedby` = status)$
      before(uiOutput(status))$
      after(tags$details(class = "mp-hdr-preview", open = NA,
        tags$summary("Preview"), uiOutput(paste0(id, "_preview"))))$allTags()
  }
  modalDialog(
    title = mp_modal_title("Export data"), size = "l", class = "mp-modal-form",
    div(style = "display:flex; flex-flow:row nowrap; gap:1em;",
      div(style = "flex:1; min-width:0;",
        shinyWidgets::pickerInput("export_group_sel", "Export group:",
          choices = groups, selected = if (isTRUE(snap$export_group_sel %in% groups)) snap$export_group_sel,
          width = "100%")),
      div(style = "flex:1; min-width:0;",
        selectizeInput("export_template", "Header template:", width = "100%",
          choices = templates, selected = template,
          options = list(create = TRUE, maxItems = 1,
            placeholder = "select or type a new template name"))),
      div(class = "mp-opts-checkbox",
        actionButton("save_template", "Save template",
          title = "Store the header text below under this template name"))),
    opts_help("Export group bundles contigs into one output set. Header template is ",
      "a reusable, named set of the FASTA header patterns below: export uses the ",
      "text on screen, and Save template keeps it for next time."),
    uiOutput("export_summary"),
    mp_checkbox("nmnh_user", tagList(tags$img(src = "www/si_sun.svg", alt = "", class = "nmnh-logo"), "NMNH user"), value = nmnh_template_on(fasta_header)),
    opts_help("Adds the ", tags$code("[specimen_voucher=]"), " and ",
      tags$code("[voucherURI=]"), " modifiers required by the NMNH Minimum Genomics ",
      "Metadata Requirements, so each GenBank record links back to its museum specimen. ",
      "Values come from your mapping file, GBIF, NCBI, or GEOME, and are checked against GBIF."),
    conditionalPanel("input.nmnh_user == true",
      div(class = "mp-nmnh-card",
        uiOutput("nmnh_sources"),
        tags$details(class = "mp-nmnh-cols",
          open = if (is.na(nmnh_cols$voucher) && is.na(nmnh_cols$uri)) NA,
          tags$summary("Change which mapping file columns are used"),
          div(style = "display:flex; flex-flow:row nowrap; gap:1em;",
            div(style = "flex:1; min-width:0;",
              selectInput("nmnh_voucher_col", "Catalog number column (e.g. USNM:FISH:12345)",
                nmnh_choices, selected = nmnh_cols$voucher %|NA|% "", width = "100%")),
            div(style = "flex:1; min-width:0;",
              selectInput("nmnh_uri_col", "Specimen link column (ARK, e.g. http://n2t.net/ark:/65665/3...)",
                nmnh_choices, selected = nmnh_cols$uri %|NA|% "", width = "100%")))),
        uiOutput("nmnh_status"))),
    hdr_box("fasta_header", "Record FASTA header:", fasta_header, "export_hdr_status"),
    mp_checkbox("export_markers", "Per-marker FASTA files", value = flag("export_markers", TRUE)),
    opts_help("Writes one FASTA per rRNA element (18S, ITS1, 5.8S, ITS2, 28S) into a ",
      "markers folder beside the group files."),
    conditionalPanel("input.export_markers == true",
      hdr_box("fasta_header_marker", "Per-marker FASTA header:", fasta_header_marker,
        "export_hdr_marker_status"),
      opts_help("Same {column} syntax, plus {element} and {product} for the marker itself.",
        nested = TRUE)),
    tags$details(
      class = "mp-meta-panel",
      meta_panel_summary(token_groups),
      opts_help("Click a column to insert it at the cursor of the header box you last clicked. ",
        "Write a column name in braces to use its value, for example ", tags$code("{Taxon}"),
        ". ", tags$code("{seqid}"), " is the record name RiboPilot gives this contig: the sample ",
        "ID, or ", tags$code("ID_p<path>_s<contig>"), " when one sample exports more than one ",
        "record. Columns from your mapping file work here even when the table does not show ",
        "them. Orange columns are empty for some records in the chosen export group; hover ",
        "one for the count.", nested = TRUE),
      if (!is.null(token_groups)) {
        export_token_ui(token_groups, target_id = "fasta_header", ns = shiny::NS(NULL), totals = totals)
      },
      opts_help(tags$code("{organism}"), " is the name set in the Taxonomy step when ",
        "your Taxon matched it exactly or you picked a species there; otherwise it is ",
        "your Taxon, so a label like \"Upeneus sp. A\" is kept. ", tags$code("{Taxon}"),
        " and ", tags$code("{sci_name}"), " give either name on its own.", nested = TRUE),
      opts_help(tags$b("isolate"), ": each record of a sample that exports more than one ",
        "contig gets ", tags$code("[isolate=ID_p<path>_s<contig>]"), " added, so all of ",
        "them can go into one GenBank submission. Your own ", tags$code("[isolate=...]"),
        " works only in a group where every sample exports one record.", nested = TRUE)),
    mp_checkbox("export_report", "HTML report with alignments", value = flag("export_report", TRUE)),
    opts_help("One HTML page with per-element alignments across the group."),
    mp_checkbox("export_itsx_check", "Final ITSx check", value = flag("export_itsx_check", TRUE)),
    opts_help("Runs stock ITSx over the exported records the way the NCBI Submission ",
      "Portal does and lists what it would flag in the summary CSV. Can take a few minutes."),
    footer = mp_footer(primary = actionButton("export_go", "Export", icon = icon("download")),
      dismiss = "Cancel")
  )
}

#' One-paragraph verdict for the export-complete dialog
#'
#' Reads the summary attached to [export_bundle()]'s return value; names the
#' records the submission-risk gate flagged, and why, and says how many the
#' ITSx check flagged, so a curator sees it before opening the CSV.
#' @param smry data.frame from `export_summary_df()` (possibly with `itsx_check`).
#' @noRd
export_check_summary_ui <- function(smry) {
  if (is.null(smry) || nrow(smry) == 0) return(NULL)
  n <- nrow(smry)
  risky <- nzchar(smry$submission_risk %||% "")
  lines <- list()
  if (any(risky)) {
    lines[[1]] <- tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
      sprintf("Submission-risk flags on %d of %d record(s): %s. See the submission_risk column in the summary CSV.",
        sum(risky), n, mp_id_list(sprintf("%s (%s)", smry$seqid[risky],
          risk_reasons(smry$submission_risk[risky])), max = 10)))
  }
  if (!is.null(smry$itsx_check)) {
    not_run <- startsWith(smry$itsx_check, "not_run(")
    if (all(not_run)) {
      lines[[length(lines) + 1]] <- tags$p(class = "mp-fg-warning",
        "ITSx check not run: ", sub("^not_run\\((.*)\\)$", "\\1", smry$itsx_check[1]), ".")
    } else {
      flagged <- sum(nzchar(smry$itsx_check) & !not_run)
      lines[[length(lines) + 1]] <- if (flagged == 0) {
        tags$p(class = if (!any(risky)) "mp-fg-success", if (!any(risky)) HTML("&#10003; "),
          sprintf("Final ITSx check: nothing flagged on %d record(s).", n))
      } else {
        tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
          sprintf("Final ITSx check flagged %d of %d record(s): %s. See the itsx_check column in the summary CSV.",
            flagged, n, paste(smry$seqid[nzchar(smry$itsx_check) & !not_run], collapse = ", ")))
      }
    }
  }
  if (length(lines) == 0) return(NULL)
  do.call(tagList, lines)
}

#' Coloured one-line verdict for a header template
#' @param chk result of `validate_fasta_header()`.
#' @noRd
header_status_ui <- function(chk) {
  if (is.null(chk)) return(NULL)
  style <- switch(chk$level,
    ok = list("success", "circle-check"),
    warn = list("warning", "triangle-exclamation"),
    list("danger", "circle-xmark"))
  tags$p(class = paste0("mp-help-text mp-fg-", style[[1]]), style = "margin:0 0 2px;",
    icon(style[[2]]), " ", chk$message)
}

#' Element colors shared by the app feature map
#' @noRd
# Colorblind-friendly (Okabe-Ito) palette, chosen to keep the 5 features
# maximally distinct: 18S blue and 28S vermillion are now clearly different
# (the old sky-blue 28S read too close to the 18S blue).
RRNA_ELEMENT_COLORS <- c(
  "18S" = "#0072B2", "ITS1" = "#E69F00", "5.8S" = "#009E73",
  "ITS2" = "#CC79A7", "28S" = "#D55E00", "ETS" = "#BBBBBB", "IGS" = "#BBBBBB"
)

#' Darken a set of (named) hex colors by a factor (for the coverage bands).
#' @noRd
darken_colors <- function(cols, f = 0.72) {
  m <- grDevices::col2rgb(cols) * f
  stats::setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255), names(cols))
}

#' Pale tint of (named) hex colors (mix f toward white) for table row shading.
#' @noRd
tint_colors <- function(cols, f = 0.80) {
  m <- grDevices::col2rgb(cols)
  m <- m + (255 - m) * f
  stats::setNames(grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255), names(cols))
}

#' Pad a plot's left gutter (y title, axis text, strips) to a fixed width so
#' plots stacked in a modal share one x axis. Returns a gtable; draw it with
#' grid::grid.draw().
#' @noRd
fixed_gutter <- function(g, left_cm = 3.2) {
  gt <- ggplot2::ggplotGrob(g)
  l <- min(gt$layout$l[grepl("panel", gt$layout$name)])
  have <- sum(grid::convertWidth(gt$widths[seq_len(l - 1)], "cm", valueOnly = TRUE))
  gt$widths[1] <- gt$widths[1] + grid::unit(max(0, left_cm - have), "cm")
  gt
}

#' Read-coverage track for one assembly unit, with feature regions shaded
#' @noRd
plot_coverage <- function(id, path = 1L, scaffold = 1L) {
  if (is.null(id)) {
    return(NULL)
  }
  cov <- app_query("SELECT depths FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  if (nrow(cov) == 0 || is.na(cov$depths[1]) || !nzchar(cov$depths[1])) {
    # a supplied-assembly sample with no reads will never have coverage
    r1 <- app_query("SELECT R1 FROM preprocess WHERE ID = ?", params = list(id))$R1
    msg <- if (length(r1) && (is.na(r1[1]) || !nzchar(r1[1]))) {
      "coverage not measured (no reads for this sample)"
    } else {
      "no coverage data (run the pipeline from Assemble > Update)"
    }
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", 0, 0, label = msg) +
      ggplot2::theme_void())
  }
  dep <- parse_depths(cov$depths)
  df <- data.frame(pos = seq_along(dep), depth = dep)
  feats <- app_query("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  ymax <- max(df$depth, 1, na.rm = TRUE)
  dark <- darken_colors(RRNA_ELEMENT_COLORS)
  g <- ggplot2::ggplot(df, ggplot2::aes(pos, depth))
  if (nrow(feats) > 0) {
    feats$mid <- (feats$start + feats$stop) / 2
    g <- g +
      # darker feature bands for visibility
      ggplot2::geom_rect(data = feats, inherit.aes = FALSE,
        ggplot2::aes(xmin = start, xmax = stop, ymin = 0, ymax = Inf, fill = element),
        alpha = 0.28) +
      ggplot2::scale_fill_manual(values = dark, guide = "none")
  }
  g <- g + ggplot2::geom_area(fill = "#444444", alpha = 0.85)
  if (nrow(feats) > 0) {
    # gene names above the corresponding boxes
    g <- g + ggplot2::geom_text(data = feats, inherit.aes = FALSE,
      ggplot2::aes(x = mid, y = ymax * 1.06, label = element, colour = element),
      size = 3.2, fontface = "bold", vjust = 0) +
      ggplot2::scale_colour_manual(values = dark, guide = "none")
  }
  # same x extent and padding as plot_chimera so the two stack with aligned axes
  n <- nrow(df)
  g + ggplot2::scale_x_continuous(breaks = seq(0, n, 1000), minor_breaks = seq(0, n, 500), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, n), ylim = c(0, ymax * 1.16), clip = "off") +
    ggplot2::labs(x = "position (bp)", y = "depth") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey78", linewidth = 0.45),
      panel.grid.minor.x = ggplot2::element_line(colour = "grey86", linewidth = 0.3),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(14, 6, 6, 6))
}

#' Coverage + annotation plot for the trim modal, with the kept window marked
#'
#' Reuses `plot_coverage()` and overlays the trim window: the flanks that would
#' be cut are greyed out and the two boundaries drawn as dashed lines.
#' @noRd
plot_trim <- function(id, path = 1L, scaffold = 1L, from = NA, to = NA) {
  g <- plot_coverage(id, path, scaffold)
  if (is.null(g)) {
    return(NULL)
  }
  if (!is.na(from) && !is.na(to) && to >= from) {
    g <- g +
      ggplot2::annotate("rect", xmin = -Inf, xmax = from, ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::annotate("rect", xmin = to, xmax = Inf, ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::geom_vline(xintercept = c(from, to), linetype = "dashed",
        colour = "#111111", linewidth = 0.7)
  }
  g
}

#' Colorblind-safe fills for lineage segments, distinct from the feature colors
#' (Tol light palette minus blue, orange, and pink). Shared by the chimera track
#' and the segments table, which doubles as the legend.
#' @noRd
CHIMERA_LINEAGE_COLORS <- c("#EEDD88", "#99DDFF", "#BBCC33", "#AAAA00", "#44BB99")
CHIMERA_SEG_ALPHA <- 0.6  # shared by the lineage bar and the table row fills

#' One fill per segment: white for no hit, grey for the sample's own lineage,
#' palette color per foreign taxon in order of first appearance.
#' @noRd
chimera_seg_colors <- function(seg, sample_lt) {
  seg_lt <- if (!is.null(seg$lineage_taxids)) seg$lineage_taxids else rep(NA_character_, nrow(seg))
  nohit <- seg$taxon %in% c("nohit", "unresolved", "high_id", "low_id")
  host <- !is.na(sample_lt) & vapply(seg_lt, function(x)
    !is.na(x) && rank_distance(x, sample_lt) < 5, logical(1))
  other <- unique(seg$taxon[!host & !nohit])
  pal <- stats::setNames(rep_len(CHIMERA_LINEAGE_COLORS, length(other)), other)
  ifelse(nohit, "#FFFFFF", ifelse(host, "#B3B3B3", pal[seg$taxon]))
}

#' The chimera call for one unit, in words
#' @noRd
chimera_result_text <- function(label, reason) {
  has <- !is.na(reason) && nzchar(reason)
  tone <- switch(label, clean = "success", severe = "danger", "warning")
  tags$span(title = if (has) reason, "Chimera call: ",
    tags$span(class = paste0("mp-pill mp-pill-", tone), label), if (has) tags$span(
    class = "mp-help-text", style = "font-weight:normal;", " (hover for the full reason)"))
}

#' A unit's chimera check as first run: before its first trim, in untrimmed
#' coordinates, with the current genes moved to match and the kept span
#' @param q Query function (sql, params).
#' @noRd
chimera_original <- function(q, id, path, scaffold) {
  ch <- q("SELECT label, junctions, segments, channels FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  win <- q("SELECT start, stop, pident, lineage FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", list(id, path, scaffold))
  feats <- q("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  n <- q("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))$length
  # a trimmed unit shows the check from before its first trim, in untrimmed coordinates
  kept <- NULL
  bak <- tryCatch(q("SELECT b.sequence bak, b.length, b.chimera, a.sequence cur FROM assembly_backup b
    JOIN assemblies a ON a.ID = b.ID AND a.path = b.path AND a.scaffold = b.scaffold
    WHERE b.ID = ? AND b.path = ? AND b.scaffold = ?", list(id, path, scaffold)), error = function(e) data.frame())
  if (nrow(bak) && !is.na(bak$bak[1]) && !is.na(bak$cur[1])) {
    off <- as.integer(regexpr(bak$cur[1], bak$bak[1], fixed = TRUE)) - 1L
    if (off >= 0L && nchar(bak$cur[1]) < bak$length[1]) {
      kept <- c(off, off + nchar(bak$cur[1]))
      feats$start <- feats$start + off; feats$stop <- feats$stop + off
      n <- bak$length[1]
      if (!is.na(bak$chimera[1])) {
        s <- jsonlite::fromJSON(bak$chimera[1])
        ch <- s$chimera
        win <- if (is.data.frame(s$windows)) s$windows[order(s$windows$start), ] else win[0, ]
      }
    }
  }
  list(ch = ch, win = win, feats = feats, n = n, kept = kept)
}

plot_chimera <- function(id, path = 1L, scaffold = 1L, db = NULL) {
  q <- if (is.null(db)) app_query else {
    con <- open_project_db(db); on.exit(DBI::dbDisconnect(con))
    function(sql, params = NULL) DBI::dbGetQuery(con, sql, params = params)
  }
  empty <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = msg) + ggplot2::theme_void()
  o <- chimera_original(q, id, path, scaffold)
  ch <- o$ch; win <- o$win; feats <- o$feats; n <- o$n; kept <- o$kept
  if (NROW(ch) == 0) return(empty("no chimera check yet (run the pipeline from Assemble > Update)"))
  seg <- if (!is.na(ch$segments[1])) jsonlite::fromJSON(ch$segments[1]) else data.frame()
  junc <- if (nzchar(ch$junctions[1] %|NA|% "")) as.integer(strsplit(ch$junctions[1], ",")[[1]]) else integer(0)
  xmax <- max(c(n, feats$stop, win$stop, seg$stop, 1), na.rm = TRUE)
  # one facet per row; the two bar rows span 0..12 so facet heights come out
  # roughly 1:1:3 against the 68..102 identity panel
  rows <- c("genes", "lineage", "ident")
  fr <- function(df, r) { df$row <- factor(r, levels = rows); df }
  g <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = fr(data.frame(y = c(0, 12, 0, 12, 68, 102)), rep(rows, each = 2)), ggplot2::aes(y = y))
  if (nrow(feats)) {
    dark <- darken_colors(RRNA_ELEMENT_COLORS)
    g <- g + ggplot2::geom_rect(data = fr(feats, "genes"), ggplot2::aes(xmin = start, xmax = stop, ymin = 1, ymax = 9, fill = element),
      alpha = 0.28) +
      ggplot2::scale_fill_manual(values = dark, guide = "none", na.value = "grey70") +
      ggplot2::geom_text(data = fr(feats, "genes"), ggplot2::aes(x = (start + stop) / 2, y = 5, label = element, colour = element),
        size = 3.2, fontface = "bold") +
      ggplot2::scale_colour_manual(values = dark, guide = "none", na.value = "grey40")
  }
  if (nrow(seg)) {
    samp <- q("SELECT lineage_taxids FROM samples WHERE ID = ?", list(id))
    seg$col <- chimera_seg_colors(seg, if (nrow(samp)) samp$lineage_taxids[1] else NA_character_)
    g <- g + ggplot2::geom_rect(data = fr(seg, "lineage"), ggplot2::aes(xmin = start, xmax = stop, ymin = 0.5, ymax = 11.5),
      fill = seg$col, alpha = CHIMERA_SEG_ALPHA, colour = "grey40", linewidth = 0.3, linetype = ifelse(seg$col == "#FFFFFF", "dotted", "solid"))
  }
  if (nrow(win)) {
    # one range line per window at its identity, dot at the midpoint filled by
    # lineage segment; no-hit windows sit hollow just under the 70% floor
    win$nohit <- is.na(win$pident)
    win$y <- ifelse(win$nohit, 68, pmax(win$pident, 70))
    win$col <- if (nrow(seg)) seg$col[findInterval(win$start, seg$start)] else "#B3B3B3"
    win$col[win$nohit] <- "#FFFFFF"
    g <- g + ggplot2::geom_segment(data = fr(win, "ident"), ggplot2::aes(x = start, xend = stop, y = y, yend = y),
        colour = ifelse(win$nohit, "grey70", "grey30"), linewidth = 0.5,
        linetype = ifelse(win$nohit, "dotted", "solid")) +
      ggplot2::geom_point(data = fr(win, "ident"), ggplot2::aes(x = (start + stop) / 2, y = y),
        shape = 21, size = 3, fill = win$col, colour = ifelse(win$nohit, "grey60", "grey30"))
  }
  sev <- identical(ch$label[1], "severe")
  jcol <- if (sev) "#c0392b" else "#b9770e"
  if (length(junc)) g <- g +
    ggplot2::geom_vline(xintercept = junc, colour = jcol, linetype = "dashed", linewidth = 0.7) +
    ggplot2::geom_label(data = fr(data.frame(x = junc, y = 12), "genes"), ggplot2::aes(x = x, y = y),
      label = sprintf("%s\njunction at %d bp", if (sev) "chimera" else "possible", junc), vjust = 0,
      size = 3.8, lineheight = 0.9, colour = jcol, fill = "white", linewidth = 0)
  if (!is.null(kept)) {
    cuts <- c(if (kept[1] > 0) kept[1], if (kept[2] < n) kept[2])
    g <- g +
      ggplot2::annotate("rect", xmin = c(0, kept[2]), xmax = c(kept[1], n), ymin = -Inf, ymax = Inf,
        fill = "#5b5b5b", alpha = 0.32) +
      ggplot2::geom_vline(xintercept = cuts, linetype = "dashed", colour = "#111111", linewidth = 0.7) +
      ggplot2::geom_text(data = fr(data.frame(x = cuts, y = 12), "genes"), ggplot2::aes(x = x, y = y),
        label = "trim", vjust = -0.3, size = 3.2, colour = "#111111")
  }
  g + ggplot2::facet_grid(row ~ ., scales = "free_y", space = "free_y", switch = "y",
      labeller = ggplot2::as_labeller(c(genes = "", lineage = "lineage", ident = "% identity"))) +
    ggplot2::scale_x_continuous(breaks = seq(0, xmax, 1000), minor_breaks = seq(0, xmax, 500), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = c(70, 80, 90, 100), expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), clip = "off") +
    ggplot2::labs(x = "position (bp)", y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "grey78", linewidth = 0.45),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(colour = "grey78", linewidth = 0.45),
      panel.grid.minor.x = ggplot2::element_line(colour = "grey86", linewidth = 0.3),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold", colour = "grey30", size = 9.5),
      axis.text = ggplot2::element_text(size = 10), panel.spacing.y = grid::unit(3, "pt"),
      plot.margin = ggplot2::margin(36, 6, 6, 6))
}
