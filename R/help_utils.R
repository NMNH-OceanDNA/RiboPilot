# CLI help + inline option help for the opts modals. Mirrors MitoPilot's
# help_utils.R: a muted one-line hint under a field (opts_help), and a "?" icon
# that toggles an in-place popover with a tool's bundled `--help` text
# (tool_help_icon), read from inst/tool_help/<tool>.txt. Both "?" kinds share
# the .mp-help-icon look and are opened and closed in custom.js.

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
    # a span, not a button: shinyjs::disabled() on a field would switch a button off
    shiny::tags$span(
      role = "button", tabindex = "0", id = ns(input_id), class = "mp-help-icon",
      `data-toggle` = "mp-tool-help", `aria-expanded` = "false",
      `aria-label` = paste("Show", label, "documentation"),
      shiny::tags$i(class = "fa-solid fa-circle-question", `aria-hidden` = "true")
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
#' @return a focusable element carrying the popover content.
#' @noRd
mp_help_tip <- function(text, label = NULL) {
  name <- if (is.null(label)) "Show help" else paste("Show help for", label)
  # a span, not a button: shinyjs::disabled() on a field would switch a button off
  shiny::tags$span(
    role = "button", tabindex = "0",
    class = "mp-help-icon",
    `data-toggle` = "mp-popover",
    `data-content` = as.character(text),
    `aria-label` = name,
    shiny::tags$i(class = "fa-solid fa-circle-question", `aria-hidden` = "true")
  )
}

#' Help popover content as a short bullet list
#' @param ... one string (or tag) per bullet.
#' @noRd
mp_help_bullets <- function(...) {
  as.character(shiny::tags$ul(class = "mp-help-list", lapply(list(...), shiny::tags$li)))
}

#' A field or column label with a help popover beside it
#'
#' @param label the visible label text.
#' @param tip help text for the popover.
#' @noRd
mp_help_label <- function(label, tip) {
  shiny::tagList(label, mp_help_tip(tip, label = label))
}
