# Where the app puts messages.
#
# Center sweet alert: the outcome of something the user deliberately did, so a
# result, a refusal, or a failure. The default.
# Bottom-right toast: incidental chatter the user can miss at no cost, e.g. a
# state toggle they can already see change in the table.
#
# The annotation editor is the deliberate exception: its messages fire every few
# seconds during boundary work, so they stay toasts.

#' Message type -> (sweetalert type, default title)
#'
#' Keyed by the `shiny::showNotification()` type names the app already used, so
#' call sites keep reading the same way.
#' @noRd
ALERT_TYPES <- list(
  message = c("success", "Done"),
  success = c("success", "Done"),
  warning = c("warning", "Heads up"),
  error = c("error", "Error"),
  info = c("info", "Note")
)

#' Center-screen alert reporting the outcome of a user action
#'
#' Sweet alerts do not queue: a second one replaces the first rather than
#' stacking. Collapse multiple messages into one call rather than looping.
#' @noRd
rp_alert <- function(msg, type = "message", title = NULL,
                     session = shiny::getDefaultReactiveDomain()) {
  spec <- ALERT_TYPES[[type]] %||% ALERT_TYPES$message
  # pass a shiny tag as `msg` when one alert has to carry several paragraphs
  shinyWidgets::sendSweetAlert(session = session, title = title %||% spec[2],
    text = msg, type = spec[1], html = !is.character(msg))
}

#' Bottom-right toast for incidental messages
#' @noRd
rp_toast <- function(msg, type = "message", duration = 5) {
  shiny::showNotification(msg, type = type, duration = duration)
}
