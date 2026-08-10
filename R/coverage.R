# Coverage: parse `samtools depth -a` output (reads mapped back to one assembly
# unit's scaffold sequence) into a per-unit coverage row.
#
# Both a mean and a MEDIAN per-position depth are stored. Coverage along an rDNA
# scaffold is not uniform: the conserved 18S/28S cores recruit multi-mapping
# reads from every copy of the array (and from any other organism's rDNA in the
# sample) while the spacers do not, so a mean is pulled up by those regions --
# and pulled up unevenly, because how much a scaffold is inflated depends on how
# much conserved core it happens to contain. That biases exactly the comparison
# the app makes with it (this scaffold's depth against its siblings'), which is
# supposed to reflect abundance rather than element composition. On real
# RiboPilot samples mean/median runs 1.01-1.26; the short conserved-core-heavy
# fragments where it should be worst are the ones a curator most needs to judge.

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
      mean_depth = 0, median_depth = 0, min_depth = 0L, max_depth = 0L, depths = "",
      time_stamp = ts, stringsAsFactors = FALSE))
    return(invisible(0))
  }
  DBI::dbAppendTable(con, "coverage", data.frame(ID = id, path = path, scaffold = scaffold,
    mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2),
    min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)),
    depths = paste(dep, collapse = ","), time_stamp = ts, stringsAsFactors = FALSE))
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
      mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2),
      min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)),
      depths = paste(dep, collapse = ","), time_stamp = ts, stringsAsFactors = FALSE))
    n <- n + 1L
  }
  invisible(n)
}

#' Reduce a per-position depth track to a bounded number of points
#'
#' Takes the value at evenly spaced positions rather than a mean, so a dropout
#' stays visible instead of being averaged away.
#' @param x Numeric depths, in position order.
#' @param n Maximum points to keep.
#' @return Numeric vector of length at most `n`, first and last preserved.
#' @noRd
downsample_depths <- function(x, n = 100) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) <= n) return(x)
  x[unique(round(seq(1, length(x), length.out = n)))]
}

#' Parse a stored `coverage.depths` blob into a numeric vector
#' @param depths_blob The stored `coverage.depths` value.
#' @return Numeric vector in position order, empty when there is no usable track.
#' @noRd
parse_depths <- function(depths_blob) {
  if (length(depths_blob) == 0 || is.na(depths_blob[1]) || !nzchar(depths_blob[1])) {
    return(numeric(0))
  }
  v <- suppressWarnings(as.numeric(strsplit(depths_blob[1], "[,;[:space:]]+")[[1]]))
  v[!is.na(v)]
}

#' Downsampled depth track as a comma-separated string for a table cell
#' @param depths_blob The stored `coverage.depths` value.
#' @param n Maximum points to keep.
#' @return Comma-separated string, or NA when there is no usable track.
#' @noRd
depth_track_string <- function(depths_blob, n = 100) {
  v <- parse_depths(depths_blob)
  if (length(v) == 0) return(NA_character_)
  paste(round(downsample_depths(v, n), 1), collapse = ",")
}
