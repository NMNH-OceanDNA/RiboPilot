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

#' The stored per-scaffold BLAST hit, labelled with the search that produced it
#'
#' Three searches can write `assemblies.blast_*`: the remote core_nt one, the
#' local one against a configured rDNA reference panel, and the local one
#' against the GetOrganelle seed file. The seed file is the same set of
#' sequences that baited the assembly, so its hit answers "which seed does this
#' scaffold resemble", not "what is this"; the other two are real reference
#' searches. Shown where no GenBank screen exists, and never without saying
#' which search it came from.
#' @param r One row of [unit_evidence()].
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
      "Best nuclear rDNA reference from the local reference panel search:"
    } else {
      "Best nuclear rDNA reference from the GenBank search:"
    }),
    div(class = "code-block", style = "white-space:normal; margin-top:3px;",
      HTML(sprintf("<b>%s</b> &nbsp; %s &nbsp; %s%% id", sp,
        if (nzchar(acc)) {
          sprintf("<a href='https://www.ncbi.nlm.nih.gov/nuccore/%s' target='_blank'>%s</a>", acc, acc)
        } else "no accession",
        r$blast_pident %|NA|% "?"))),
    if (from_seeds) tags$p(class = "mp-help-text mp-fg-warning", style = "margin-top:3px;",
      paste("These seeds are what GetOrganelle used to bait the assembly, so a high",
        "identity here says the scaffold looks like its own bait. A missing or",
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

#' Prev / picker / next sample navigation for the detail modals
#'
#' The order is whatever the main table shows when the modal opens: sort,
#' column filters, search, and the lock/state rows the filter CSS hides. The
#' script snapshots it into input$nav_order; the server renders the controls.
#' @param pfx "assemble" or "annotate": names the lock/state row classes.
#' @noRd
rp_sample_nav <- function(pfx) {
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
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      "A chimera is one scaffold built from pieces of different origin: two organisms, ",
      "or rDNA copies joined in the wrong place. Four checks look for one: the 18S, 5.8S, and 28S ",
      "GenBank hits must agree on a lineage; 500 bp windows BLASTed against the local rDNA panel ",
      "are grouped into segments by class; ribotyper checks the 18S and 28S structure; and reads ",
      "clipped at one position, or a lasting step in read depth, mark a possible join."),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
      "Severe: two parts match different phyla or kingdoms. The export calls the record chimeric ",
      "and asks before including it. Minor: weaker evidence, such as a class-level switch or read ",
      "evidence alone. It is listed in the export report but not flagged to GenBank. Windows in ",
      "conserved 18S or 28S that match several classes equally well are shown as unresolved and ",
      "do not count."),
    tags$p(class = "mp-help-text", style = "margin:2px 0 4px;", textOutput("chimera_db_source", inline = TRUE)),
    plotOutput("chimera_plot", height = "300px"),
    reactable::reactableOutput("chimera_segments"),
    rp_details("All windows",
      tags$p(class = "mp-help-text", style = "margin:2px 0 4px;",
        "One row per window: its best BLAST hit in the local panel, the hit identity, and the read evidence over that stretch."),
      reactable::reactableOutput("chimera_windows"), open = FALSE))
}

#' Assemble details modal: scaffolds, their evidence, and keep/drop
#' @noRd
assemble_modal <- function(id) {
  modalDialog(
    title = div(class = "mp-sticky-head", paste0("Assemble - ", id)),
    size = "l", easyClose = TRUE,
    footer = div(class = "mp-sticky-foot", rp_sample_nav("assemble"), mp_footer(dismiss = "Close")),
    modal_taxon_header(id),
    h5("Scaffolds", style = "margin:2px 0 4px;"),
    tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
      paste("Click on a scaffold row to examine GenBank evidence and coverage below")),
    reactable::reactableOutput("assemble_scaffolds"),
    div(style = "margin-top:6px;",
      actionButton("copy_fasta", "Copy FASTA", icon = icon("clipboard"), class = "btn-sm"),
      tags$span(class = "text-muted", style = "margin-left:8px; font-size:0.82em;",
        "Copies the ticked scaffolds as FASTA, headed by sample, taxon and unit.")),
    tags$hr(style = "margin:10px 0;"),

    # The evidence that settles the keep/drop call, at the point it is made.
    # blast_screen and blast_ref_candidates are both written per unit by WF1 and
    # neither was rendered anywhere: when they name different organisms, that is
    # the contamination signal.
    # All start collapsed, and collapse again when another scaffold is picked
    # (examine_unit), so the modal opens short for every scaffold.
    div(id = "asm_unit_panels",
      rp_details(tagList("GenBank evidence for ", textOutput("hl_scaf_ev", inline = TRUE)), open = FALSE,
        uiOutput("unit_screen"),
        div(style = "margin-top:8px; color:#555;",
          "Best matches among records whose title says nuclear rDNA:"),
        reactable::reactableOutput("unit_candidates")),

      rp_details(tagList("Read coverage for ", textOutput("hl_scaf_cov", inline = TRUE)), open = FALSE,
        plotOutput("unit_coverage_plot", height = "140px")),

      rp_details(tagList("Chimera check for ", textOutput("hl_scaf_chim", inline = TRUE)), open = FALSE,
        chimera_panel())),

    # Only a supplied assembly has these: the contigs the rDNA search found but
    # did not carry forward. They never become assembly units, so this is the
    # one place they are visible at all.
    rejected_contigs_section(id)
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

#' Warn that locking would leave a sample with two locked assembly paths
#'
#' A corner-of-the-screen notification was too easy to miss for something that
#' silently does nothing, so this blocks until the reviewer answers. Locking two
#' paths is unusual rather than impossible -- alternative resolutions of the same
#' repeat are normally competing answers, not both right -- so "Lock anyway"
#' stays available.
#' @param ids Sample IDs that would end up with more than one locked path.
#' @noRd
multi_path_lock_modal <- function(ids) {
  modalDialog(
    title = mp_modal_title("More than one assembly path would be locked"),
    size = "m", easyClose = FALSE,
    div(style = "font-size:1.02em; line-height:1.5;",
      tags$p(paste("Paths are competing reconstructions of the same rDNA repeat, so a",
        "sample normally carries one locked path. This would leave more than one",
        "locked for:")),
      tags$p(style = "font-weight:600; margin:6px 0;", paste(ids, collapse = ", ")),
      tags$p(style = "color:#666;",
        paste("Unlock the other path first, or proceed if you mean to keep both",
          "(both will be annotated and exported)."))),
    footer = mp_footer(danger = actionButton("lock_anyway", "Lock anyway"))
  )
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

#' Ask before exporting units the contamination checks flagged
#'
#' `flagged` is rows of [export_summary_df()] whose submission_risk carries an
#' off_target or contaminant tag and that no curator has marked Reviewed.
#' @noRd
export_screen_modal <- function(flagged) {
  why <- vapply(flagged$submission_risk, function(r) {
    tx <- regmatches(r, regexec("off_target\\(([^)]*)\\)", r))[[1]][2]
    paste(c(if (!is.na(tx)) sprintf("GenBank hits agree on %s", tx),
      if (grepl("(^|;)contaminant(;|$)", r)) "the best GenBank hit is from another phylum",
      if (grepl("(^|;)chimeric\\(ITSx\\)", r)) "ITSx, the tool NCBI runs on submissions, calls it chimeric"
      else if (grepl("(^|;)chimeric", r)) "parts of the scaffold match organisms from different phyla or kingdoms"),
      collapse = "; ")
  }, character(1))
  items <- lapply(seq_len(nrow(flagged)), function(i) tags$li(
    tags$b(sprintf("%s p%s.s%s", flagged$ID[i], flagged$path[i], flagged$scaffold[i])),
    sprintf(" (%s): %s.", flagged$organism[i] %|NA|% "no Taxon", why[i])))
  modalDialog(
    title = mp_modal_title("Possible contamination or chimera"),
    size = "m", easyClose = FALSE,
    div(style = "font-size:1.02em; line-height:1.5;",
      tags$p(sprintf("%s in this group may contain another organism's rDNA, exported under the sample's name:",
        mp_n(nrow(flagged), "unit"))),
      tags$ul(items),
      tags$p(style = "color:#666;",
        "Leave them out unless you have checked them. Marking a unit Reviewed in",
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
#' @param actions Named list key -> action, from [triage_opts_for()].
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
    title = paste("Triage scaffolds for", mp_n(n_samples, "sample")),
    size = "l", easyClose = TRUE,
    tags$p(class = "mp-help-text",
      paste("Proposes keep / drop / undecided per scaffold from the evidence already",
        "in the database. Depth is shown in the scaffold table but never decides:",
        "high depth is often k-mer collapse across the conserved core or IGS",
        "subrepeats, and low depth is often a repeat inside the scaffold, where reads",
        "map to two places and are left out.")),
    tags$table(style = "width:100%; font-size:0.9em;",
      tags$thead(tags$tr(
        tags$th(style = "text-align:left;", "Rule"),
        tags$th(style = "text-align:left;", "Action"),
        tags$th(style = "text-align:right;", "Decided / Matching"))),
      tags$tbody(rows)),
    tags$p(class = "mp-help-text", style = "margin-top:6px;",
      paste("Rules are applied top to bottom and the first one that matches decides",
        "whether to drop or flag the scaffold. \"Decided\" counts scaffolds this rule",
        "will drop or flag. \"Matching\" counts every scaffold that meets the rule,",
        "even if a higher-priority rule acts on it first.")),
    div(style = "margin:0.4em 0;",
      shinyWidgets::prettyCheckbox("triage_overwrite_user",
        "Also re-decide scaffolds I set by hand", value = FALSE, status = "danger")),
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
  unit_lbl <- sprintf(" [path %s / scaffold %s]", path, scaffold)
  # Seed the notes box from the DB. Re-opening the same unit assigns an
  # identical() value to sel_unit, which does not invalidate the reactiveVal, so
  # the server-side reset observer never fires; a value-less textAreaInput then
  # reports "" and the write observer persists it over the stored note.
  notes <- tryCatch(app_query(
    "SELECT annotate_notes FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, as.integer(path), as.integer(scaffold)))$annotate_notes[1],
    error = function(e) NA_character_)
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
    footer = div(class = "mp-sticky-foot", rp_sample_nav("annotate"), mp_footer(
      primary = actionButton("lock_and_close", "Lock & Close"), dismiss = NULL,
      extra = actionButton("annotate_close", "Close"))),
    modal_taxon_header(id),
    # ---- gene list (drives the active feature for the editor + alignment) ----
    h5("Feature model", style = "margin:2px 0 4px;"),
    tags$p(class = "mp-help-text", style = "margin-bottom:4px;",
      "Click a feature to make it active for editing + reference alignment below."),
    reactable::reactableOutput("feature_tbl"),
    div(style = "margin:6px 0; display:flex; gap:6px;",
      actionButton("trim_btn", "Trim assembly...", icon = icon("scissors"), class = "btn-sm"),
      actionButton("restore_btn", "Restore assembly", icon = icon("rotate-left"), class = "btn-sm")),
    tags$hr(style = "margin:10px 0;"),

    # ---- collapsible detail panels ----
    # Coverage, boundary editing, and the alignment are the review itself, so
    # they open by default; the Assemble window opens short because its panels
    # are per-scaffold evidence the reviewer picks through.
    rp_details("Read coverage", plotOutput("coverage_plot", height = "130px")),

    rp_details("Edit active feature boundaries",
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
          shinyWidgets::prettyCheckbox("feat_edit", "Edit", value = FALSE,
            status = "primary", inline = TRUE)),
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
          div(style = "display:flex; gap:6px; align-items:flex-end;",
            actionButton("feat_save", "Save", icon = icon("floppy-disk"), class = "btn-sm btn-success"),
            actionButton("feat_reset", "Reset", icon = icon("rotate-left"), class = "btn-sm"))))
      ),
      tags$p(class = "mp-help-text", style = "margin-top:4px;",
        "Check Edit to adjust boundaries; + extends outward, - trims inward. Save to keep, Reset to discard."),
      div(class = "code-block", style = "margin-top:6px; white-space:normal; word-break:break-all;",
        uiOutput("nudge_seq"))),

    rp_details("Reference alignment (active feature)",
      uiOutput("msa_header"),
      uiOutput("feat_align_ui")),

    rp_details("Curation warnings", open = FALSE,
      uiOutput("curate_notes"),
      reactable::reactableOutput("flag_tbl")),

    rp_details("Chimera check", open = FALSE, chimera_panel()),

    rp_details("Notes", open = FALSE,
      textAreaInput("annotate_notes", label = NULL, width = "100%",
        value = (notes %||% NA_character_) %|NA|% "",
        placeholder = "Notes for this annotation unit"),
      tags$p(class = "mp-help-text", "Saved as you type. Deleting the text does not clear the stored note; use Clear notes."),
      actionButton("annotate_notes_clear", "Clear notes", class = "btn-sm"))
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
          "No scheduler (sbatch/qsub/bsub) found on PATH. Save the script and submit it yourself; running from the login node is discouraged."
        } else {
          sprintf("%s detected. Edit the job script if needed, then Submit as job (or Save script).", sched$name)
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
        shinyjs::hidden(actionButton("nf_stop", "Stop / Interrupt")),
        actionButton("nf_close", "Close")))
  }
  modalDialog(
    title = run_modal_title(sprintf("%s - %d sample(s) to update", wf, n), "rp_gears"),
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
#' The same shape as [run_modal()] so a database build looks like any other
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
          "No scheduler (sbatch/qsub/bsub) found on PATH. Save the script and submit it yourself; running from the login node is discouraged."
        } else {
          sprintf("%s detected. Edit the job script if needed, then Submit as job (or Save script).", sched$name)
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
        shinyjs::hidden(actionButton("ref_nf_stop", "Stop / Interrupt")),
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

#' Trim-assembly modal: cut the scaffold to the annotated feature span
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
    opts_help("Trimming cuts the flanks off the scaffold and shifts every ",
      "feature coordinate + the coverage track to match. Features that fall entirely ",
      "outside the kept range are removed. Use \"Restore assembly\" to undo."),
    footer = mp_footer(
      danger = actionButton("trim_apply", "Trim", icon = icon("scissors")), dismiss = NULL,
      # actionButton (not modalButton) so Cancel reopens the annotate details
      # modal instead of dropping back to the main Annotate table.
      extra = actionButton("trim_cancel", "Cancel"))
  )
}

#' Assign-export-group modal: name the group for the selected samples
#' @noRd
group_modal <- function() {
  modalDialog(
    title = "Assign export group", size = "s", easyClose = TRUE,
    textInput("group_name", "Group name", placeholder = "e.g. batch_2026_01"),
    footer = mp_footer(primary = actionButton("assign_group_go", "Assign"))
  )
}

#' Export modal: pick which export group to bundle (or all locked samples), then
#' write it to the project's "out/export" folder. Without this, assigned groups
#' were ignored and every bundle collapsed into one undifferentiated "export".
#' @noRd
export_modal <- function(groups, templates = "default", template = "default",
                         fasta_header = "", fasta_header_marker = "",
                         columns = character(0), sample_cols = character(0)) {
  yours <- sort(intersect(columns, sample_cols))
  ours <- sort(setdiff(columns, yours))
  hdr_box <- function(id, label, value) {
    htmltools::tagQuery(textAreaInput(id, label, value, width = "100%", height = "60px"))$
      find("textarea")$addAttrs(`aria-describedby` = paste0(id, "_status"))$allTags()
  }
  modalDialog(
    title = mp_modal_title("Export data"), size = "l", class = "mp-modal-form",
    div(style = "display:flex; flex-flow:row nowrap; gap:1em;",
      div(style = "flex:1; min-width:0;",
        shinyWidgets::pickerInput("export_group_sel", "Export group:",
          choices = groups, width = "100%")),
      div(style = "flex:1; min-width:0;",
        selectizeInput("export_template", "Header template:", width = "100%",
          choices = templates, selected = template,
          options = list(create = TRUE, maxItems = 1,
            placeholder = "select or type a new template name"))),
      div(class = "mp-opts-checkbox",
        actionButton("save_template", "Save template",
          title = "Store the header text below under this template name"))),
    opts_help("Export group bundles units into one output set. Header template is ",
      "a reusable, named set of the FASTA header patterns below: export uses the ",
      "text on screen, and Save template keeps it for next time."),
    hdr_box("fasta_header", "Record FASTA header:", fasta_header),
    uiOutput("export_hdr_status"),
    tags$details(
      tags$summary("Available columns"),
      opts_help("Write a column name in braces to use its value, for example ",
        tags$code("{Taxon}"), ". ", tags$code("{seqid}"), " is the record name ",
        "RiboPilot gives this unit: the sample ID, or ",
        tags$code("ID_p<path>_s<scaffold>"), " when one sample exports more than ",
        "one record. Columns from your mapping file work here even when the ",
        "table does not show them.", nested = TRUE),
      opts_help(tags$b("isolate"), " is filled by RiboPilot: records of a sample ",
        "that exports more than one unit get ", tags$code("[isolate=ID_p<path>_s<scaffold>]"),
        " appended, because GenBank rejects several sequences with the same ID and ",
        "voucher. Leave it out of the template.", nested = TRUE),
      tags$p(tags$b("Your columns: "), if (length(yours)) paste(yours, collapse = ", ") else "none"),
      tags$p(tags$b("RiboPilot columns: "), if (length(ours)) paste(ours, collapse = ", ") else "none")),
    mp_checkbox("export_markers", "Per-marker FASTA files", value = TRUE),
    opts_help("Writes one FASTA per rRNA element (18S, ITS1, 5.8S, ITS2, 28S) into a ",
      "markers folder beside the group files."),
    conditionalPanel("input.export_markers == true",
      hdr_box("fasta_header_marker", "Per-marker FASTA header:", fasta_header_marker),
      uiOutput("export_hdr_marker_status"),
      opts_help("Same {column} syntax, plus {element} and {product} for the marker itself.",
        nested = TRUE)),
    mp_checkbox("export_report", "HTML report with alignments", value = TRUE),
    opts_help("One HTML page with per-element alignments across the group."),
    mp_checkbox("export_itsx_check", "Final ITSx check", value = TRUE),
    opts_help("Runs stock ITSx over the exported records the way the NCBI Submission ",
      "Portal does and lists what it would flag in the summary CSV. Can take a few minutes."),
    uiOutput("export_summary"),
    footer = mp_footer(primary = actionButton("export_go", "Export", icon = icon("download")),
      dismiss = "Cancel")
  )
}

#' One-paragraph verdict for the export-complete dialog
#'
#' Reads the summary attached to [export_bundle()]'s return value and says how
#' many records the ITSx check and the submission-risk gate flagged, so a
#' curator sees it before opening the CSV.
#' @param smry data.frame from [export_summary_df()] (possibly with `itsx_check`).
#' @noRd
export_check_summary_ui <- function(smry) {
  if (is.null(smry) || nrow(smry) == 0) return(NULL)
  n <- nrow(smry)
  risk <- sum(nzchar(smry$submission_risk %||% ""))
  lines <- list()
  if (!is.null(smry$itsx_check)) {
    not_run <- startsWith(smry$itsx_check, "not_run(")
    if (all(not_run)) {
      lines[[length(lines) + 1]] <- tags$p(class = "mp-fg-warning",
        "ITSx check not run: ", sub("^not_run\\((.*)\\)$", "\\1", smry$itsx_check[1]), ".")
    } else {
      flagged <- sum(nzchar(smry$itsx_check) & !not_run)
      lines[[length(lines) + 1]] <- if (flagged == 0) {
        tags$p(class = "mp-fg-success", HTML("&#10003; "),
          sprintf("Final ITSx check: nothing flagged on %d record(s).", n))
      } else {
        tags$p(class = "mp-fg-warning", style = "font-weight:bold;",
          sprintf("Final ITSx check flagged %d of %d record(s): %s. See the itsx_check column in the summary CSV.",
            flagged, n, paste(smry$seqid[nzchar(smry$itsx_check) & !not_run], collapse = ", ")))
      }
    }
  }
  if (risk > 0) {
    lines[[length(lines) + 1]] <- tags$p(class = "mp-fg-warning",
      sprintf("Submission-risk flags on %d record(s) (submission_risk column).", risk))
  }
  if (length(lines) == 0) return(NULL)
  do.call(tagList, lines)
}

#' Coloured one-line verdict for a header template
#' @param chk result of [validate_fasta_header()].
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

#' Read-coverage track for one assembly unit, with feature regions shaded
#' @noRd
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
  dep <- as.numeric(strsplit(cov$depths[1], ",")[[1]])
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
#' Reuses [plot_coverage()] and overlays the trim window: the flanks that would
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

#' Lineage segments, window identity, and junctions under the depth track
#' @noRd
plot_chimera <- function(id, path = 1L, scaffold = 1L, db = NULL) {
  q <- if (is.null(db)) app_query else {
    con <- open_project_db(db); on.exit(DBI::dbDisconnect(con))
    function(sql, params = NULL) DBI::dbGetQuery(con, sql, params = params)
  }
  empty <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = msg) + ggplot2::theme_void()
  ch <- q("SELECT label, junctions, segments, channels FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  if (nrow(ch) == 0) return(empty("no chimera check yet (run the pipeline from Assemble > Update)"))
  seg <- if (!is.na(ch$segments[1])) jsonlite::fromJSON(ch$segments[1]) else data.frame()
  win <- q("SELECT start, stop, pident, lineage FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", list(id, path, scaffold))
  feats <- q("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  junc <- if (nzchar(ch$junctions[1] %|NA|% "")) as.integer(strsplit(ch$junctions[1], ",")[[1]]) else integer(0)
  n <- q("SELECT length FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))$length
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

#' Synteny map: the sample's annotations over its BLAST reference's
#' (when the reference has been fetched + annotated), aligned by element color.
#' @noRd
