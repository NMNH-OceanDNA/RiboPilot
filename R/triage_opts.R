# Per-rule triage actions, read from the triage_opts preset table so the Shiny
# app's Apply and the WF1 triage process cannot disagree. All resolution is
# defensive: an old project, a missing preset or a hand-edited value falls back
# to the shipped defaults rather than erroring or arming an unsupported drop.

#' Column name holding a rule's action
#' @noRd
triage_opts_col <- function(key) paste0("rule_", key)

#' Rule actions from one preset
#' @param con Open project DB connection.
#' @param preset Preset name.
#' @return Named list key -> action.
#' @noRd
triage_opts_row <- function(con, preset = "default") {
  r <- tryCatch(
    DBI::dbGetQuery(con, "SELECT * FROM triage_opts WHERE triage_opts = ?",
      params = list(preset %||% "default")),
    error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(triage_actions_default())
  x <- list()
  for (k in triage_rule_keys()) {
    col <- triage_opts_col(k)
    if (col %in% names(r)) x[[k]] <- r[[col]][1]
  }
  normalize_triage_actions(x)
}

#' Rule actions for one sample, following assemble.triage_opts
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return Named list key -> action.
#' @noRd
triage_opts_for <- function(con, id) {
  p <- tryCatch(
    DBI::dbGetQuery(con, "SELECT triage_opts AS t FROM assemble WHERE ID = ?",
      params = list(id))$t,
    error = function(e) character(0))
  triage_opts_row(con, if (length(p) && !is.na(p[1]) && nzchar(p[1])) p[1] else "default")
}

#' Write rule actions into a preset, creating it if needed
#' @param con Open project DB connection.
#' @param actions Named list key -> action; missing keys take the default.
#' @param preset Preset name.
#' @return (invisibly) the preset name.
#' @noRd
write_triage_opts <- function(con, actions, preset = "default") {
  a <- normalize_triage_actions(actions)
  keys <- triage_rule_keys()
  cols <- vapply(keys, triage_opts_col, character(1))
  # one transaction: a failed insert after a committed delete would silently
  # discard a customised preset and leave the project on the defaults
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM triage_opts WHERE triage_opts = ?", params = list(preset))
    DBI::dbExecute(con,
      sprintf("INSERT INTO triage_opts (triage_opts, %s) VALUES (%s)",
        paste(cols, collapse = ", "),
        paste(rep("?", length(cols) + 1), collapse = ", ")),
      params = c(list(preset), unname(a[keys])))
  })
  invisible(preset)
}
