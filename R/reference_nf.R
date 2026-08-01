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

#' Container bind option for the reference folder, if the config uses a container
#'
#' The folder lives outside the project, so a container task cannot see it. The
#' mount is added per run rather than written into the project config, so a user
#' who moves the folder needs no project edits.
#' @noRd
ref_bind_option <- function(cfg_lines, ref_dir) {
  txt <- paste(cfg_lines, collapse = "\n")
  d <- path.expand(ref_dir)
  if (grepl("docker\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions", paste0("-v ", d, ":", d)))
  }
  if (grepl("(singularity|apptainer)\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions", paste0("--bind ", d)))
  }
  character(0)
}

#' Nextflow arguments for the reference workflow
#' @noRd
nf_reference_command <- function(project, entry, manifest, ref_dir, out_dir,
                                 api_key = "", refseq = FALSE, resume = TRUE,
                                 work = file.path(project, "work_ref")) {
  entry <- match.arg(entry, c("COUNT", "BUILD"))
  main_nf <- app_sys("nextflow", "main_reference.nf")
  cfg <- file.path(project, ".config")
  cfg_lines <- if (file.exists(cfg)) readLines(cfg) else character(0)
  c("run", main_nf, "-entry", entry,
    if (file.exists(cfg)) c("-c", cfg),
    "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    "--clades", path.expand(manifest),
    "--refDir", path.expand(ref_dir),
    "--outDir", path.expand(out_dir),
    "--ribopilot_src", ribopilot_src_default(),
    if (isTRUE(refseq)) c("--refseq", "true"),
    if (nzchar(api_key %||% "")) c("--ncbi_api_key", api_key),
    ref_bind_option(cfg_lines, ref_dir))
}
