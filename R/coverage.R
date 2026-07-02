# WF1 coverage: parse `samtools depth -a` output (reads mapped back to the
# stitched nrDNA sequence) into a per-sample coverage row.

#' Write per-position read depth for one sample's stitched sequence
#'
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param depth_file `samtools depth -a` output (ref, pos, depth), tab-separated.
#' @return (invisibly) the mean depth.
#' @export
write_coverage <- function(con, id, depth_file) {
  ts <- as.integer(Sys.time())
  DBI::dbExecute(con, "DELETE FROM coverage WHERE ID = ?", params = list(id))

  dep <- integer(0)
  if (file.exists(depth_file) && file.info(depth_file)$size > 0) {
    d <- utils::read.table(depth_file, sep = "\t", header = FALSE,
      col.names = c("ref", "pos", "depth"), stringsAsFactors = FALSE)
    dep <- d$depth[order(d$pos)]
  }
  if (length(dep) == 0) {
    DBI::dbAppendTable(con, "coverage", data.frame(ID = id, mean_depth = 0,
      min_depth = 0L, max_depth = 0L, depths = "", time_stamp = ts,
      stringsAsFactors = FALSE))
    return(invisible(0))
  }
  DBI::dbAppendTable(con, "coverage", data.frame(ID = id,
    mean_depth = round(mean(dep), 2), min_depth = as.integer(min(dep)),
    max_depth = as.integer(max(dep)), depths = paste(dep, collapse = ","),
    time_stamp = ts, stringsAsFactors = FALSE))
  invisible(mean(dep))
}
