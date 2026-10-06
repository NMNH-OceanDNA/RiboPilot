# CLI help + inline option help for the opts modals. Mirrors MitoPilot's
# help_utils.R: a muted one-line hint under a field (opts_help), and a "?" icon
# that toggles an in-place popover with a tool's bundled `--help` text
# (tool_help_icon + register_tool_help), read from inst/tool_help/<tool>.txt.

#' Brief inline help text shown directly under an options-modal field
#' @noRd
opts_help <- function(..., href = NULL, link_text = "learn more", id = NULL,
                      nested = FALSE) {
  inner <- list(...)
  if (!is.null(href)) {
    inner <- c(inner, list(
      " (",
      shiny::tags$a(href = href, target = "_blank", rel = "noopener", link_text),
      ")"
    ))
  }
  margin_top <- if (nested) "2px" else "-6px"
  shiny::tags$p(
    id = id,
    class = "text-muted mp-help-text",
    style = paste0("margin-top: ", margin_top, "; margin-bottom: 14px;"),
    inner
  )
}

#' Read a bundled tool help text file (inst/tool_help/<tool>.txt)
#' @noRd
read_tool_help <- function(tool) {
  f <- app_sys("tool_help", paste0(tool, ".txt"))
  if (!nzchar(f) || !file.exists(f)) {
    return(paste0("No bundled help for '", tool, "'. Regenerate with ",
      "tools/capture_tool_help.sh against the RiboPilot image."))
  }
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

#' Small inline "?" icon that toggles an in-place popover with the tool's
#' bundled `--help` text. The popover is an absolutely-positioned overlay
#' (styles inline, so it floats regardless of external CSS) shown/hidden
#' client-side, so opening help never tears down the options modal and never
#' discards unsaved edits. Same approach as MitoPilot.
#' @noRd
tool_help_icon <- function(tool,
                           label = paste(tool, "options"),
                           session = shiny::getDefaultReactiveDomain()) {
  ns <- session$ns
  input_id <- paste0("help_", tool)
  panel_id <- paste0("help_panel_", tool)
  shiny::tags$span(
    style = "position: relative; display: inline-block;",
    shiny::actionLink(
      ns(input_id), label = NULL, icon = shiny::icon("circle-question"),
      title = paste("Show", label, "documentation"),
      style = "color: #888; margin-left: 4px;"
    ),
    shiny::tags$div(
      id = ns(panel_id), class = "rp-tool-help-popover",
      style = paste(
        "display: none; position: absolute; z-index: 1080;",
        "left: 0; top: 1.6em; width: 720px; max-width: 90vw;",
        "max-height: 50vh; overflow-y: auto; text-align: left;",
        "font-weight: normal; background: #fff; color: #333; border: 1px solid #ccc;",
        "border-radius: 4px; box-shadow: 0 4px 16px rgba(0,0,0,0.2);",
        "padding: 10px;"
      ),
      shiny::tags$pre(
        style = "white-space: pre-wrap; font-size: 12px; margin: 0; background: transparent; border: none; padding: 0;",
        read_tool_help(tool)
      )
    )
  )
}

#' Register an observer that toggles a tool-help popover when its icon is
#' clicked. Call once per tool_help_icon() in the server.
#' @noRd
register_tool_help <- function(tool, input,
                               session = shiny::getDefaultReactiveDomain()) {
  input_id <- paste0("help_", tool)
  panel_id <- paste0("help_panel_", tool)
  shiny::observeEvent(input[[input_id]], ignoreInit = TRUE, {
    shinyjs::toggle(id = session$ns(panel_id))
  })
}

#' Small "?" icon that opens a help popover on click or keyboard focus
#'
#' The one help affordance: a grey `circle-question` always means "click for
#' a popover". The popover itself is wired up once in `custom.js` for every
#' `[data-toggle="mp-popover"]` on the page, so this works inside modals and
#' inside dynamically rendered UI with no server code.
#'
#' @param text help text, one or two sentences (HTML allowed).
#' @param label optional name of the thing being explained, used in the
#'   accessible name.
#' @return a focusable button carrying the popover content.
#' @noRd
mp_help_tip <- function(text, label = NULL) {
  name <- if (is.null(label)) "Show help" else paste("Show help for", label)
  shiny::tags$button(
    type = "button",
    class = "mp-help-icon",
    `data-toggle` = "mp-popover",
    `data-content` = as.character(text),
    `aria-label` = name,
    shiny::tags$i(class = "fa-solid fa-circle-question", `aria-hidden` = "true")
  )
}

#' A field or column label with a help popover beside it
#'
#' @param label the visible label text.
#' @param tip help text for the popover.
#' @noRd
mp_help_label <- function(label, tip) {
  shiny::tagList(label, mp_help_tip(tip, label = label))
}
