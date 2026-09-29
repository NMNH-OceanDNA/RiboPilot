# Click-time guards. Each returns TRUE to continue, or explains why the click
# cannot proceed and returns FALSE. Replaces bare req() calls that used to abort
# silently (theme T01).

#' Shorten an ID list for a message: five IDs, then "and N more".
#' @noRd
mp_id_list <- function(ids, max = 5) {
  ids <- unique(ids)
  if (length(ids) <= max) return(paste(ids, collapse = ", "))
  paste0(paste(ids[seq_len(max)], collapse = ", "), ", and ", length(ids) - max, " more")
}

#' At least one row must be selected.
#' @noRd
need_selection <- function(n, session = getDefaultReactiveDomain()) {
  if (length(n) > 0 && !is.na(n[1]) && n[1] > 0) return(TRUE)
  mp_toast("Select one or more rows in the table first.", type = "warning", session = session)
  FALSE
}

#' None of the selected units may be locked. `locked_ids` are the IDs that are.
#' @noRd
need_unlocked <- function(locked_ids, noun = "sample",
                          how = "Select them and press Unlock first",
                          session = getDefaultReactiveDomain()) {
  if (length(locked_ids) == 0) return(TRUE)
  plural <- sub("^\\S+ ", "", mp_n(2, noun))
  mp_alert(
    title = paste("Locked", plural, "cannot be edited"),
    text = paste0("Locked: ", mp_id_list(locked_ids), ". ", how, "."),
    type = "warning", session = session
  )
  FALSE
}

#' Only one assembly path per sample may be annotation-locked (MitoPilot rule).
#' `u` is data.frame(ID, path, scaffold) of the units about to be locked.
#' @noRd
need_one_path <- function(u, session = getDefaultReactiveDomain()) {
  ids <- unique(u$ID)
  cur <- app_query(sprintf(
    "SELECT ID, path, scaffold, COALESCE(annotate_lock,0) AS lock, problematic FROM annotate WHERE ID IN (%s)",
    paste(rep("?", length(ids)), collapse = ",")), params = as.list(ids))
  # excluded units are never exported, so they do not count
  on <- !(cur$problematic %in% "yes") & (cur$lock == 1 |
    unit_key(cur$ID, cur$path, cur$scaffold) %in% unit_key(u$ID, u$path, u$scaffold))
  bad <- names(Filter(function(p) length(unique(p)) > 1, split(cur$path[on], cur$ID[on])))
  if (length(bad) == 0) return(TRUE)
  mp_alert(
    title = "Only one assembly path can be locked per sample",
    text = paste0(mp_n(length(bad), "sample"), " would have more than one assembly path locked: ",
      mp_id_list(bad), ". Paths are alternative reconstructions of the same rDNA repeat, ",
      "and a sample exports one. Lock the units of one path and leave the other paths unlocked."),
    type = "warning", session = session
  )
  FALSE
}
