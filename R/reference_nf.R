# Launching the reference workflow, and getting its results back into R.
#
# Counts come back as one small file per clade rather than appended rows: the
# tasks run in parallel and would otherwise race on a single cache file.

#' Write the clade work list the reference workflow reads
#' @noRd
write_clade_manifest <- function(path, clades) {
  cols <- c("ref_id", "taxid", "clade_name", "query")
  d <- clades[, cols, drop = FALSE]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(d, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  path
}

#' Read the per-clade count files a COUNT run produced
#' @noRd
read_count_files <- function(dir) {
  empty <- data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE)
  files <- list.files(dir, pattern = "\\.count$", full.names = TRUE)
  if (length(files) == 0) return(empty)
  rows <- lapply(files, function(f) {
    x <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
    x <- trimws(x[nzchar(trimws(x))])
    if (length(x) == 0) return(NULL)
    p <- strsplit(x[1], "\t")[[1]]
    if (length(p) < 2) return(NULL)
    tx <- suppressWarnings(as.integer(p[1])); n <- suppressWarnings(as.integer(p[2]))
    if (is.na(tx) || is.na(n)) return(NULL)
    data.frame(taxid = tx, n_records = n, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(empty)
  do.call(rbind, rows)
}

#' Merge fresh counts into a reference folder's cache
#' @noRd
merge_counts <- function(ref_dir, counts) {
  cache <- read_counts_cache(ref_dir)
  if (nrow(counts) > 0) {
    cache <- cache[!(cache$taxid %in% counts$taxid), , drop = FALSE]
    cache <- rbind(cache, counts[, c("taxid", "n_records"), drop = FALSE])
  }
  dir.create(path.expand(ref_dir), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(cache, ref_counts_path(ref_dir), sep = "\t",
    row.names = FALSE, quote = FALSE, na = "")
  cache
}
