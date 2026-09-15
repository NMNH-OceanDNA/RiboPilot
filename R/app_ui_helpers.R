# Shared UI building blocks. One toolbar button, one filter picker, one modal
# footer, one modal title, one CSV row, one pluraliser. See dev/ui_review
# themes T05, T10, T13, T16, T23.

#' Toolbar button with the shared class set and a tooltip.
#'
#' A disabled `<button>` fires no pointer events, so a `title` on it never
#' shows. When `needs_selection` is TRUE the button is wrapped in a span that
#' carries the requirement tooltip and the button keeps its own action tooltip.
#'
#' @noRd
mp_toolbar_button <- function(id, label, icon = NULL,
                              emphasis = c("default", "primary", "danger"),
                              title = NULL, needs_selection = FALSE) {
  emphasis <- match.arg(emphasis)
  btn <- actionButton(id, label, icon = icon, title = title)
  btn$attribs$class <- paste0("btn btn-sm action-button mp-toolbar-btn btn-", emphasis)
  if (!is.null(title) && (is.null(label) || (is.character(label) && !nzchar(label)))) {
    btn <- tagAppendAttributes(btn, `aria-label` = title)
  }
  if (!needs_selection) {
    return(btn)
  }
  tags$span(
    class = "mp-needs-selection-wrap",
    title = "Select one or more rows in the table first",
    tagAppendAttributes(btn, class = "mp-needs-selection")
  )
}

#' Multi-select filter picker for a table filter row.
#' @noRd
mp_filter_picker <- function(id, label, choices, selected = choices, width = "150px") {
  shinyWidgets::pickerInput(
    inputId = id, label = label, choices = choices, selected = selected,
    multiple = TRUE, width = width,
    options = list(
      `actions-box` = TRUE,
      `select-all-text` = "All",
      `deselect-all-text` = "None",
      `selected-text-format` = "count > 1",
      `count-selected-text` = "{0} of {1}"
    )
  )
}

#' Modal footer: extra, danger, dismiss, primary (primary right-most).
#'
#' @noRd
mp_footer <- function(primary = NULL, dismiss = "Cancel", danger = NULL, extra = NULL) {
  emph <- function(x, cls) if (is.null(x)) NULL else tagAppendAttributes(x, class = cls)
  parts <- list(
    extra,
    emph(danger, "btn-danger"),
    if (is.null(dismiss)) NULL else modalButton(dismiss),
    emph(primary, "btn-primary")
  )
  do.call(tagList, parts[!vapply(parts, is.null, logical(1))])
}

#' Modal title with a Bootstrap close X and an optional subtitle.
#'
#' `modalDialog()` already wraps whatever it is handed in its own
#' `h4.modal-title`, so this returns a plain div: one heading per dialog.
#'
#' @noRd
mp_modal_title <- function(text, subtitle = NULL, close = TRUE) {
  div(
    class = "mp-modal-title",
    if (close) {
      tags$button(
        type = "button", class = "close", `data-dismiss` = "modal",
        `aria-label` = "Close", tags$span(`aria-hidden` = "true", HTML("&times;"))
      )
    },
    text,
    if (is.null(subtitle)) NULL else tags$p(class = "text-muted mp-modal-subtitle", subtitle)
  )
}

#' Count plus noun, pluralised: "1 sample", "3 samples", "2 assemblies".
#'
#' Errors on `NA` rather than putting "NA" on screen.
#'
#' @noRd
mp_n <- function(n, noun, plural = NULL) {
  if (length(n) != 1 || is.na(n)) stop("mp_n() needs one non-missing count")
  if (n == 1) {
    return(paste(n, noun))
  }
  if (is.null(plural)) {
    plural <- if (grepl("[^aeiou]y$", noun)) sub("y$", "ies", noun) else paste0(noun, "s")
  }
  paste(n, plural)
}

#' Checkbox with a visible tick and an accessible name.
#'
#' `prettyCheckbox()` draws a filled square with no glyph and never binds its
#' label to the input; the check icon and the aria-label fix both (theme T18).
#'
#' @noRd
mp_checkbox <- function(inputId, label, value = FALSE, ...) {
  x <- shinyWidgets::prettyCheckbox(
    inputId, label = label, value = value, status = "primary", icon = icon("check"), ...
  )
  lab <- if (is.character(label)) label else paste(htmltools::doRenderTags(label), collapse = " ")
  htmltools::tagQuery(x)$find("input")$addAttrs(`aria-label` = lab)$allTags()
}


#' The row-count and selection line above a sample table
#' @noRd
mp_table_status <- function(n_visible, n_total, n_selected, noun = "sample") {
  div(
    class = "mp-table-status", role = "status", `aria-live` = "polite",
    span(
      "Showing ", tags$b(n_visible), " of ", tags$b(n_total), " ",
      sub("^\\S+ ", "", mp_n(n_total, noun))
    ),
    span(class = "mp-sep", `aria-hidden` = "true", HTML("&middot;")),
    span(tags$b(n_selected), " selected"),
    # Filtered-out rows are hidden by CSS, so reactable's own noData never fires.
    if (n_visible == 0) {
      span(class = "mp-empty-state", "No rows match the current filters.")
    }
  )
}

#' Canonical header text for a table column
#' @noRd
mp_col_name <- function(field) unname(MP_COL_NAMES[[field]])

#' Canonical header, with its tooltip when the name needs one
#' @noRd
mp_col_header <- function(field, tip = NULL) {
  if (is.null(tip)) {
    tip <- unname(MP_COL_TIPS[field])
    if (is.na(tip)) tip <- NULL
  }
  rt_header(mp_col_name(field), tip)
}

#' Column width that fits the longest value in `x`
#' @noRd
mp_fit_width <- function(x, floor = 120, cap = 260) {
  n <- suppressWarnings(max(nchar(as.character(x)), na.rm = TRUE))
  if (!is.finite(n)) {
    return(floor)
  }
  max(floor, min(cap, ceiling(8.5 * n) + 26))
}

#' Toast after an options set is applied to rows
#' @noRd
mp_opts_saved_toast <- function(n, set, noun = "sample",
                                session = getDefaultReactiveDomain()) {
  mp_toast(
    paste0(mp_n(n, noun), " switched to ", set, " and set to Ready to run."),
    session = session
  )
}

#' State glyphs for a table icon column, tone class included
#' @noRd
mp_state_icon_classes <- function(module) {
  codes <- MP_STATE_CODES[[module]]
  tone <- c(`0` = "mp-fg-neutral", `1` = "mp-fg-neutral", `4` = "mp-fg-neutral",
            `2` = "mp-fg-success", `3` = "mp-fg-danger")
  stats::setNames(paste(mp_state_icons(module), tone[codes]), codes)
}

#' Accessible name for each state glyph: the label plus its meaning.
#' @noRd
mp_state_titles <- function(module) {
  codes <- MP_STATE_CODES[[module]]
  stats::setNames(
    paste0(mp_state_labels(module), " - ", mp_state_tips(module)),
    codes
  )
}

#' A decorative glyph: the button text or aria-label is the name.
#' @noRd
mp_icon <- function(name) {
  i <- icon(name)
  i$attribs$`aria-label` <- NULL
  i$attribs$role <- NULL
  tagAppendAttributes(i, `aria-hidden` = "true")
}

#' Modal that sets the pipeline state of the selected rows
#' @noRd
mp_state_modal <- function(ids, current = character(0), mode = "assemble") {
  codes <- MP_STATE_SETTABLE
  effect <- c(
    `0` = "Skipped by the next Update.",
    `1` = "Processed by the next Update; results this module already stored are replaced.",
    `4` = "Picked up by the next run; normally set by the pipeline itself.",
    `2` = "Treated as finished; the next Update skips it.",
    `3` = "Treated as finished with an error; the next Update skips it."
  )
  choice_names <- lapply(codes, function(k) {
    tagList(
      tags$strong(MP_STATE_META[[k]]$label),
      tags$span(class = "mp-help-text", paste0(" - ", effect[[k]]))
    )
  })
  modalDialog(
    title = mp_modal_title(paste("Set state for", mp_n(length(ids), MP_UNIT[[mode]]))),
    tags$p(class = "text-muted", mp_id_list(ids)),
    if (length(current) == 1 && current == "4") {
      tags$p(class = "text-muted",
             paste("These are In progress. The pipeline sets that state;",
                   "choose another to take them out of it."))
    } else if (length(current) == 0) {
      tags$p(class = "text-muted", "The selected rows are not all in the same state.")
    },
    shinyWidgets::prettyRadioButtons(
      "new_state",
      label = "New state",
      choiceValues = codes,
      choiceNames = choice_names,
      selected = if (length(current) == 1 && current %in% codes) current else character(0),
      shape = "round",
      status = "primary"
    ),
    size = "m",
    easyClose = TRUE,
    footer = mp_footer(primary = actionButton("update_state", "Set state"))
  )
}
