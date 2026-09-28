# Coverage: parse `samtools depth -a` output (reads mapped back to one assembly
# unit's scaffold sequence) into a per-unit coverage row. The pipeline maps with
# a MAPQ >= 10 filter, so reads that map equally well to several of the sample's
# scaffolds are excluded and a depth reflects reads placed on that scaffold
# specifically. Conserved-core stretches can therefore read as low or zero depth.
#
# Both a mean and a MEDIAN per-position depth are stored. Coverage along an rDNA
# scaffold is not uniform: the conserved 18S/28S cores and the spacers do not
# behave the same way, so a mean is skewed by those regions -- and skewed
# unevenly, because how far a scaffold is off depends on how
# much conserved core it happens to contain. That biases exactly the comparison
# the app makes with it (this scaffold's depth against its siblings'), which is
# supposed to reflect abundance rather than element composition. On real
# RiboPilot samples mean/median runs 1.01-1.26; the short conserved-core-heavy
# fragments where it should be worst are the ones a curator most needs to judge.

#' Split a combined `samtools depth -a` (one path mapped at a time) into
#' per-unit coverage rows
#'
#' The reference names are the scaffold fasta headers `path.scaffold`, so each
#' distinct ref is written as one coverage unit.
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param depth_file `samtools depth -a` output (ref, pos, depth), ref = path.scaffold.
#' @param clips_file Optional per-position soft-clip counts (ref, pos, n), tab-separated.
#' @return (invisibly) the number of units written.
#' @export
write_coverage_units <- function(con, id, depth_file, clips_file = NULL) {
  # a dropped scaffold is not re-measured, so its track is kept for a put-back
  DBI::dbExecute(con,
    "DELETE FROM coverage WHERE ID = ? AND NOT EXISTS (SELECT 1 FROM assemblies a
       WHERE a.ID = coverage.ID AND a.path = coverage.path AND a.scaffold = coverage.scaffold
         AND COALESCE(a.ignore, 0) = 1)", params = list(id))
  if (!file.exists(depth_file) || file.info(depth_file)$size == 0) {
    return(invisible(0L))
  }
  d <- utils::read.table(depth_file, sep = "\t", header = FALSE,
    col.names = c("ref", "pos", "depth"), colClasses = c("character", "integer", "integer"))
  clips <- NULL
  if (!is.null(clips_file) && file.exists(clips_file) && file.size(clips_file) > 0) {
    clips <- utils::read.table(clips_file, sep = "\t", col.names = c("ref", "pos", "n"),
      colClasses = c("character", "integer", "integer"), quote = "", comment.char = "")
  }
  clip_track <- function(ref, len) {
    v <- integer(len)
    if (!is.null(clips)) {
      c1 <- clips[clips$ref == ref & clips$pos >= 1 & clips$pos <= len, , drop = FALSE]
      if (nrow(c1)) v[c1$pos] <- v[c1$pos] + c1$n
    }
    paste(v, collapse = ",")
  }
  ts <- as.integer(Sys.time())
  n <- 0L
  for (r in unique(d$ref)) {
    ps <- strsplit(as.character(r), "\\.")[[1]]
    if (length(ps) != 2) next
    sub <- d[d$ref == r, ]
    dep <- sub$depth[order(sub$pos)]
    if (length(dep) == 0) next
    DBI::dbExecute(con, "DELETE FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
      params = list(id, as.integer(ps[1]), as.integer(ps[2])))
    DBI::dbAppendTable(con, "coverage", data.frame(ID = id,
      path = as.integer(ps[1]), scaffold = as.integer(ps[2]),
      mean_depth = round(mean(dep), 2), median_depth = round(stats::median(dep), 2),
      min_depth = as.integer(min(dep)), max_depth = as.integer(max(dep)),
      depths = paste(dep, collapse = ","), clips = clip_track(as.character(r), length(dep)),
      time_stamp = ts, stringsAsFactors = FALSE))
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
