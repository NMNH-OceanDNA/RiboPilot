# Filter choices, shared by the pickers and the server-side visible-row mask so
# the two cannot drift (that drift is what made "Hold" a no-op).
LOCK_CHOICES <- c("Unlocked" = "0", "Locked" = "1")

#' Lock + state status filter pickers (mirrors MitoPilot); values match the
#' assemble/annotate lock (0/1) and switch (0/1/2/3/4) codes. Each mode gets its
#' own picker ids (`<mode>_lock_filter` / `<mode>_state_filter`) so the Assemble
#' and Annotate filters do not share a DOM id and leak across modes.
#' @noRd
filter_pickers <- function(mode) {
  shiny::tagList(
    mp_filter_picker(paste0(mode, "_lock_filter"), "Lock:", LOCK_CHOICES, width = "140px"),
    mp_filter_picker(paste0(mode, "_state_filter"), "State:", STATE_CHOICES, width = "140px")
  )
}

#' RiboPilot app UI
#' @import shiny
#' @noRd
app_ui <- function(request) {
  btn <- mp_toolbar_button
  # curator marks; app_server.R relabels them to what the next click does
  mark <- function(col) btn(paste0(col, "_btn"), MP_FLAG_WORDS[[col]][["set"]],
    title = paste(MP_FLAG_WORDS[[col]][["tip"]], "Acts on the selected units; locked units are skipped."),
    needs_selection = TRUE)
  step <- function(n, lab) shiny::HTML(sprintf("<span class='mp-step-n'>%d</span> %s", n, lab))
  h1 <- function(mode, lab) shiny::conditionalPanel(sprintf("input.mode == '%s'", mode),
    shiny::tags$h1(class = "mp-module-h1 sr-only", lab))

  structure(lang = "en", shiny::tagList(
    add_external_resources(),
    shiny::fluidPage(
      shiny::tags$a(href = "#mp-content", class = "sr-only sr-only-focusable", "Skip to table"),
      shiny::div(
        style = "display:flex; flex-direction:column;",
        shiny::tags$header(
          role = "banner",
          shiny::uiOutput("app_header"),
          shiny::div(
            class = "mp-toolbar",
            # "Annotate" stays the value so the tab can be renamed without
            # touching the server.
            shinyWidgets::radioGroupButtons("mode", label = "Step:",
              choiceNames = list(step(1, "Taxonomy"), step(2, "Reference"), step(3, "Assemble"),
                step(4, "Annotation Review"), step(5, "Export")),
              choiceValues = c("Taxonomy", "Reference", "Assemble", "Annotate", "Export"),
              selected = "Taxonomy", size = "sm"),
            btn("refresh", label = NULL, icon = mp_icon("arrows-rotate"),
              title = "Reload the table from the database"),
            shiny::conditionalPanel("input.mode == 'Taxonomy'", id = "tax_ctrls", class = "mp-toolbar",
              # Toolbar ids are per mode: every conditionalPanel stays in the DOM, so
              # a shared id would bind several buttons to one input.
              btn("tax_state_btn", "State", title = "Set the pipeline state of the selected samples",
                needs_selection = TRUE),
              btn("tax_lock_btn", "Lock", title = MP_LOCK_DEF("taxonomy"), needs_selection = TRUE),
              btn("tax_run_btn", "Assign Taxonomy", icon = mp_icon("circle-play"), emphasis = "primary",
                title = "Assign taxonomy to every unlocked sample that is Ready to run"),
              btn("edit_taxon_btn", "Edit Taxonomy", icon = mp_icon("pen-to-square"),
                title = "Search the NCBI taxonomy and assign a taxon to the selected samples",
                needs_selection = TRUE)),
            # Only the two controls that act on the whole tab live here. Everything
            # that acts on ONE group is in that group's own row (see app_reference.R).
            shiny::conditionalPanel("input.mode == 'Reference'", id = "ref_ctrls", class = "mp-toolbar",
              htmltools::tagAppendAttributes(
                shinyWidgets::radioGroupButtons("ref_view", label = NULL, size = "sm",
                  choices = c("Samples", "Groups"), selected = "Samples"),
                style = "margin-bottom:0;"),
              btn("ref_suggest_btn", "Propose groups", icon = mp_icon("wand-magic-sparkles"),
                emphasis = "primary", title = "Group samples by clade from the assigned taxonomy"),
              btn("ref_dir_btn", "Settings", icon = mp_icon("gear"),
                title = "Reference folder and record limits")),
            shiny::conditionalPanel("input.mode == 'Reference' && input.ref_view == 'Samples'",
              id = "ref_sample_ctrls", class = "mp-toolbar",
              btn("ref_move_btn", "Move", title = "Move the selected samples to a reference group",
                needs_selection = TRUE),
              btn("ref_keep_btn", "Keep together",
                title = "Make a user group from the selected samples", needs_selection = TRUE),
              btn("ref_release_btn", "Release",
                title = "Return the selected samples to automatic grouping", needs_selection = TRUE)),
            shiny::conditionalPanel("input.mode == 'Assemble'", id = "asmb_ctrls", class = "mp-toolbar",
              btn("asmb_state_btn", "State", title = "Set the pipeline state of the selected samples",
                needs_selection = TRUE),
              btn("asmb_lock_btn", "Lock", title = MP_LOCK_DEF("assemble"), needs_selection = TRUE),
              btn("triage_btn", "Triage", icon = mp_icon("traffic-light"),
                title = "Re-run automatic scaffold triage on the selected samples. Locked samples are skipped.",
                needs_selection = TRUE),
              btn("asmb_trim_btn", "Trim", icon = mp_icon("scissors"),
                title = MP_TRIM_TIP("samples")[["trim"]], needs_selection = TRUE),
              btn("asmb_restore_btn", "Restore", icon = mp_icon("rotate-left"),
                title = MP_TRIM_TIP("samples")[["restore"]], needs_selection = TRUE),
              btn("asmb_run_btn", "Update", icon = mp_icon("circle-play"), emphasis = "primary",
                title = "Review and launch the RiboPilot pipeline")),
            # No Update button here: one workflow, launched from Assemble.
            shiny::conditionalPanel("input.mode == 'Annotate'", id = "annot_ctrls", class = "mp-toolbar",
              btn("annot_state_btn", "State", title = "Set the pipeline state of the selected units",
                needs_selection = TRUE),
              btn("annot_lock_btn", "Lock", title = MP_LOCK_DEF("annotate"), needs_selection = TRUE),
              mark("reviewed"), mark("problematic"), mark("partial"),
              btn("annot_trim_btn", "Trim", icon = mp_icon("scissors"),
                title = MP_TRIM_TIP("units")[["trim"]], needs_selection = TRUE),
              btn("annot_restore_btn", "Restore", icon = mp_icon("rotate-left"),
                title = MP_TRIM_TIP("units")[["restore"]], needs_selection = TRUE)),
            shiny::conditionalPanel("input.mode == 'Export'", id = "export_ctrls", class = "mp-toolbar",
              btn("assign_group_btn", "Assign Group",
                title = "Assign the selected units to an export group", needs_selection = TRUE),
              btn("clear_group_btn", "Clear Group",
                title = "Remove the selected units from their export group", needs_selection = TRUE),
              btn("export_btn", "Export Data", icon = mp_icon("download"), emphasis = "primary",
                title = "Open the export window"),
              shiny::span(shiny::textOutput("export_sel_n", inline = TRUE),
                style = "color:#888; font-size:0.85em;")),
            shiny::div(style = "margin-left:auto; display:flex; gap:0.9em; align-items:center;",
              # per-SAMPLE Nextflow task dirs: only meaningful on Assemble
              shiny::conditionalPanel("input.mode == 'Assemble'",
                btn("workdirs_btn", "Work Dirs", icon = mp_icon("folder-tree"),
                  title = "Nextflow work directories of the selected sample",
                  needs_selection = TRUE)))
          )
        ),
        shiny::tags$main(
          id = "mp-content", style = "padding: 1em;",
          h1("Taxonomy", "Taxonomy"), h1("Reference", "Reference"), h1("Assemble", "Assemble"),
          h1("Annotate", "Annotation Review"), h1("Export", "Export"),
          shiny::conditionalPanel("input.mode == 'Assemble'",
            shiny::div(class = "mp-filter-row",
              filter_pickers("Assemble"),
              shiny::div(class = "mp-filter-cols",
                mp_filter_picker("col_groups", "Columns:", c("Stats", "Options"), width = "170px")))),
          shiny::conditionalPanel("input.mode == 'Annotate'",
            shiny::div(class = "mp-filter-row",
              filter_pickers("Annotate"),
              shiny::div(class = "mp-filter-cols",
                mp_filter_picker("ann_col_groups", "Columns:", c("BLAST", "Check details"), width = "170px")))),
          shiny::uiOutput("table_status"),
          shiny::div(class = "mp-table-resize", reactable::reactableOutput("tbl", height = "100%")),
          # CSV downloads sit under the table, as in MitoPilot.
          shiny::div(class = "mp-csv-row",
            shiny::span(class = "mp-needs-selection-wrap",
              title = "Select one or more rows in the table first",
              shiny::downloadButton("table_csv_selected", "Download selected rows",
                class = "btn-sm btn-default mp-needs-selection")),
            shiny::downloadButton("table_csv", "Download all rows", class = "btn-sm btn-default"))
        )
      )
    )
  ))
}
