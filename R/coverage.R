# Coverage: parse `samtools depth -a` output (reads mapped back to one assembly
# unit's scaffold sequence) into a per-unit coverage row.

#' Write per-position read depth for one assembly unit (ID, path, scaffold)
#'
#' @param con Open project DB connection.
#' @param id,path,scaffold Assembly unit.
#' @param depth_file `samtools depth -a` output (ref, pos, depth), tab-separated.
#' @return (invisibly) the mean depth.
#' @export
write_coverage <- function(con, id, path = 1L, scaffold = 1L, depth_file) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))

  dep <- integer(0)
  if (file.exists(depth_file) && file.info(depth_file)$size > 0) {
    d <- utils::read.table(depth_file, sep = "\t", header = FALSE,
      col.names = c("ref", "pos", "depth"), stringsAsFactors = FALSE)
    dep <- d$depth[order(d$pos)]
  }
  if (length(dep) == 0) {
    DBI::dbAppendTable(con, "coverage", data.frame(ID = id, path = path, scaffold = scaffold,
      mean_depth = 0, min_depth = 0L, max_depth = 0L, depths = "", time_stamp = ts,
      stringsAsFactors = FALSE))
    return(invisible(0))
  }
  DBI::dbAppendTable(con, "coverage", data.frame(ID = id, path = path, scaffold = scaffold,
    mean_depth = round(mean(dep), 2), min_depth = as.integer(min(dep)),
    max_depth = as.integer(max(dep)), depths = paste(dep, collapse = ","),
    time_stamp = ts, stringsAsFactors = FALSE))
  invisible(mean(dep))
}

#' Split a combined `samtools depth -a` (all scaffolds mapped at once) into
#' per-unit coverage rows
#'
#' The reference names are the scaffold fasta headers `path.scaffold`, so each
#' distinct ref is written as one coverage unit.
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param depth_file `samtools depth -a` output (ref, pos, depth), ref = path.scaffold.
#' @return (invisibly) the number of units written.
#' @export
write_coverage_units <- function(con, id, depth_file) {
  DBI::dbExecute(con, "DELETE FROM coverage WHERE ID = ?", params = list(id))
  if (!file.exists(depth_file) || file.info(depth_file)$size == 0) {
    return(invisible(0L))
  }
  d <- utils::read.table(depth_file, sep = "\t", header = FALSE,
    col.names = c("ref", "pos", "depth"), stringsAsFactors = FALSE)
  ts <- as.integer(Sys.time())
  n <- 0L
  for (r in unique(d$ref)) {
    ps <- strsplit(as.character(r), "\\.")[[1]]
    if (length(ps) != 2) next
    sub <- d[d$ref == r, ]
    dep <- sub$depth[order(sub$pos)]
    if (length(dep) == 0) next
    DBI::dbAppendTable(con, "coverage", data.frame(ID = id,
      path = as.integer(ps[1]), scaffold = as.integer(ps[2]),
      mean_depth = round(mean(dep), 2), min_depth = as.integer(min(dep)),
      max_depth = as.integer(max(dep)), depths = paste(dep, collapse = ","),
      time_stamp = ts, stringsAsFactors = FALSE))
    n <- n + 1L
  }
  invisible(n)
}
