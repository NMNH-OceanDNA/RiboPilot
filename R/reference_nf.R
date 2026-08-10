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
  container_bind_opts(cfg_lines, ref_dir)
}

#' Bind options mounting each directory into the container the config uses
#'
#' @param cfg_lines Lines of the project `.config`.
#' @param dirs Directories the tasks must be able to read.
#' @return character(0), or c("-process.containerOptions", "<value>").
#' @noRd
container_bind_opts <- function(cfg_lines, dirs) {
  dirs <- unique(path.expand(dirs[!is.na(dirs) & nzchar(dirs)]))
  if (length(dirs) == 0) return(character(0))
  txt <- paste(cfg_lines, collapse = "\n")
  if (grepl("docker\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions",
      paste(paste0("-v ", dirs, ":", dirs), collapse = " ")))
  }
  if (grepl("(singularity|apptainer)\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) {
    return(c("-process.containerOptions", paste(paste0("--bind ", dirs), collapse = " ")))
  }
  character(0)
}

#' Directories the project config's own runOptions already bind-mount
#'
#' Mounting them again would be harmless but noisy, and hides what this run
#' actually had to add.
#' @noRd
config_mounted_dirs <- function(cfg_lines) {
  txt <- paste(cfg_lines, collapse = "\n")
  # perl = TRUE throughout: TRE mishandles \\s inside a bracket expression
  ro <- regmatches(txt, gregexpr('runOptions\\s*=\\s*"[^"]*"', txt, perl = TRUE))[[1]]
  if (length(ro) == 0) return(character(0))
  v <- unlist(regmatches(ro, gregexpr("-v\\s+[^\\s:]+", ro, perl = TRUE)))
  b <- unlist(regmatches(ro, gregexpr("(?:-B|--bind)\\s+[^\\s,]+", ro, perl = TRUE)))
  d <- c(sub("^-v\\s+", "", v, perl = TRUE), sub("^(?:-B|--bind)\\s+", "", b, perl = TRUE))
  d <- d[!grepl("${", d, fixed = TRUE)]   # ${params.projectDir} is passed separately
  unique(path.expand(d))
}

#' Directories holding the reference databases a project's presets point at
#'
#' The Reference tab writes absolute paths that can live anywhere, and the
#' project config only mounts what was known at project creation, so the
#' pipeline has to mount them itself or the assembler cannot read its seeds.
#'
#' @param project Project directory.
#' @param mounted Directories the config already mounts; anything inside one of
#'   them is dropped.
#' @return character vector of directories, possibly empty.
#' @noRd
reference_mount_dirs <- function(project, mounted = character(0)) {
  con <- tryCatch(DBI::dbConnect(RSQLite::SQLite(), file.path(project, ".sqlite")),
    error = function(e) NULL)
  if (is.null(con)) return(character(0))
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  dbs <- tryCatch(DBI::dbGetQuery(con,
    "SELECT seeds_db, labels_db FROM assemble_opts"), error = function(e) NULL)
  paths <- if (is.null(dbs)) character(0) else c(dbs$seeds_db, dbs$labels_db)
  paths <- paths[!is.na(paths) & nzchar(paths)]
  # only what a preset actually names: the configured reference folder may hold
  # databases this project never uses, and mounting it would say otherwise
  dirs <- unique(dirname(path.expand(paths)))
  dirs <- dirs[dir.exists(dirs)]
  dirs <- drop_covered(dirs, path.expand(mounted[nzchar(mounted %||% "")]))
  # a folder and its own seeds/ and labels/ subfolders would be three mounts of
  # the same tree; the parent covers them
  unname(drop_covered(dirs, dirs))
}

#' Drop directories already covered by one of `by` (a parent, or the same path)
#'
#' A directory is never treated as covering itself, so passing the same vector
#' as both arguments collapses children into their parent.
#' @noRd
drop_covered <- function(dirs, by) {
  if (length(dirs) == 0 || length(by) == 0) return(dirs)
  inside <- function(d, m) !identical(d, m) && startsWith(paste0(d, "/"), paste0(m, "/"))
  keep <- vapply(dirs, function(d) {
    !any(vapply(by, function(m) identical(d, m) && !identical(dirs, by) || inside(d, m),
      logical(1)))
  }, logical(1))
  unname(dirs[keep])
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
  # every container executor template interpolates params.projectDir and
  # params.rawDir into docker/singularity runOptions, so a run that omits them
  # cannot even parse the project config (see nf_command())
  proj_dir <- path.expand(project)
  raw_dir <- config_param(cfg_lines, "rawDir") %||% ""
  if (length(raw_dir) == 0 || is.na(raw_dir) || !nzchar(raw_dir)) raw_dir <- proj_dir
  c("run", main_nf, "-entry", entry,
    if (file.exists(cfg)) c("-c", cfg),
    "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    "--projectDir", proj_dir,
    "--rawDir", path.expand(raw_dir),
    "--clades", path.expand(manifest),
    "--refDir", path.expand(ref_dir),
    "--outDir", path.expand(out_dir),
    "--ribopilot_src", ribopilot_src_default(),
    if (isTRUE(refseq)) c("--refseq", "true"),
    if (nzchar(api_key %||% "")) c("--ncbi_api_key", api_key),
    ref_bind_option(cfg_lines, ref_dir))
}

#' The NCBI API key a project should use
#'
#' new_project() stores the key in the project config, so the tab should not
#' make the user type it again; the environment variable is the fallback.
#' @noRd
project_ncbi_key <- function(project) {
  cfg <- file.path(project %||% "", ".config")
  v <- if (file.exists(cfg)) {
    config_param(readLines(cfg, warn = FALSE), "ncbi_api_key")
  } else {
    NA_character_
  }
  if (length(v) == 0 || is.na(v) || v %in% c("NA", "NULL") || !nzchar(trimws(v))) {
    return(Sys.getenv("NCBI_API_KEY", ""))
  }
  trimws(v)
}

#' Run a Nextflow command to completion, returning its exit status
#'
#' Split out so tests can substitute a runner that spawns nothing.
#' @param log Optional file the run's output is written to, so a failure has a
#'   diagnostic instead of a discarded stream.
#' @noRd
nf_run_blocking <- function(project, args, log = NULL) {
  if (!nzchar(Sys.which("nextflow"))) {
    stop("Nextflow not found on PATH.", call. = FALSE)
  }
  # unpinned, a Nextflow 26 install rejects -entry with "not supported with the
  # strict parser", so every count run fails on an otherwise fine workflow
  nf_env <- c("current", COLUMNS = "500")
  pin <- nf_pin_version()
  if (!is.na(pin)) nf_env <- c(nf_env, NXF_VER = pin)
  p <- processx::run("nextflow", args, wd = project, error_on_status = FALSE,
    env = nf_env, stdout = "|", stderr = "2>&1")
  if (!is.null(log) && nzchar(log %||% "")) {
    append_run_log(log, strsplit(p$stdout %||% "", "\n", fixed = TRUE)[[1]])
  }
  p$status
}

#' Count GenBank records for a set of clades, waiting for the run
#'
#' Only uncached clades are sent to Nextflow; the walk reads the cache.
#' @noRd
run_clade_counts <- function(project, ref_dir, clades, api_key = "", refseq = FALSE,
                             runner = nf_run_blocking) {
  if (is.null(project) || length(project) != 1 || is.na(project) ||
      !nzchar(project) || !dir.exists(project)) {
    stop("run_clade_counts: 'project' must be an existing directory.", call. = FALSE)
  }
  cache <- read_counts_cache(ref_dir)
  # a cached NA is not an answer: keying on presence alone would skip that clade
  # for good, so only rows with a real count count as cached
  counted <- cache$taxid[!is.na(cache$taxid) & !is.na(cache$n_records)]
  need <- clades[!(as.integer(clades$taxid) %in% counted), , drop = FALSE]
  need$taxid <- as.integer(need$taxid)
  need <- need[!duplicated(need$taxid), , drop = FALSE]
  if (nrow(need) == 0) return(cache)
  rundir <- file.path(project, ".reference")
  out <- file.path(rundir, "counts")
  unlink(out, recursive = TRUE)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  man <- write_clade_manifest(file.path(rundir, "clades_count.tsv"), need)
  log <- file.path(rundir, paste0("count_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
  unlink(log); file.create(log)
  args <- nf_reference_command(project, "COUNT", man, ref_dir, out,
    api_key = api_key, refseq = refseq, resume = FALSE)
  # a non-zero status is not fatal: clades that did report are still worth
  # keeping, and the walk treats a missing count as unknown
  status <- runner(project, args, log)
  if (!identical(status, 0L) && !identical(status, 0)) {
    warning("run_clade_counts: Nextflow COUNT run exited with status ", status,
      call. = FALSE)
  }
  res <- merge_counts(ref_dir, read_count_files(out))
  # carried so the caller's failure message can quote the run that produced it
  attr(res, "log") <- log
  res
}
