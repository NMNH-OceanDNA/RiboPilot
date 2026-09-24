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

#' Bind options mounting each directory into the container the config uses
#'
#' @param cfg_lines Lines of the project `.config`.
#' @param dirs Directories the tasks must be able to read.
#' @return character(0), or c("-process.containerOptions", "<value>").
#' @noRd
container_bind_opts <- function(cfg_lines, dirs) {
  dirs <- unique(path.expand(dirs[!is.na(dirs) & nzchar(dirs)]))
  if (length(dirs) == 0) return(character(0))
  switch(config_container_engine(cfg_lines) %|NA|% "",
    docker = c("-process.containerOptions",
      paste(paste0("-v ", shQuote(dirs), ":", shQuote(dirs)), collapse = " ")),
    singularity = c("-process.containerOptions",
      paste(paste0("--bind ", shQuote(dirs)), collapse = " ")),
    character(0))
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
  con <- tryCatch(open_project_db(file.path(project, ".sqlite")),
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
  # the cross-project reference store also lives under this folder, and the
  # fetch step's priming/push calls need it mounted just like a seeds/labels db.
  # Read reference_config.ref_dir directly, not through ref_config(), whose
  # fallback chain would otherwise mount the running machine's default folder
  # for a project that never set one (see nf_ref_dir()).
  ref <- tryCatch(DBI::dbGetQuery(con, "SELECT ref_dir FROM reference_config WHERE rowid_key = 1")$ref_dir,
    error = function(e) NA_character_)
  ref <- if (length(ref) > 0) ref[1] else NA_character_
  if (!is.na(ref) && nzchar(trimws(ref))) dirs <- c(dirs, path.expand(trimws(ref)))
  # the local BLAST panel is a database PREFIX, not a directory; mount its
  # containing folder, and only when a project has explicitly configured one
  bdb <- tryCatch(DBI::dbGetQuery(con, "SELECT blastdb FROM reference_config WHERE rowid_key = 1")$blastdb,
    error = function(e) NA_character_)
  bdb <- if (length(bdb) > 0) bdb[1] else NA_character_
  if (!is.na(bdb) && nzchar(trimws(bdb))) dirs <- c(dirs, dirname(path.expand(trimws(bdb))))
  dirs <- unique(dirs)
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
  proj_dir <- path.expand(project)
  # ref_dir is deliberately not bind-mounted: the tasks write into their own
  # work directory and Nextflow publishes to the folder from the host, so the
  # mount only gave Singularity a missing path to refuse
  c("run", main_nf, "-entry", entry,
    if (file.exists(cfg)) c("-c", cfg),
    "-work-dir", work,
    if (isTRUE(resume)) "-resume",
    # the container runOptions interpolate params.projectDir
    "--projectDir", proj_dir,
    "--clades", path.expand(manifest),
    "--refDir", path.expand(ref_dir),
    "--outDir", path.expand(out_dir),
    "--ribopilot_src", ribopilot_src_default(),
    if (isTRUE(refseq)) c("--refseq", "true"),
    if (nzchar(api_key %||% "")) c("--ncbi_api_key", api_key))
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

#' Stderr from the tasks that failed during a run
#'
#' An ignored failure leaves only a NOTE in the run log, so the reason the task
#' died stays in its work directory and never reaches the user. The task hash in
#' the NOTE is what connects the two.
#'
#' @param log The run log written by [nf_run_blocking()].
#' @param work The run's work directory.
#' @param n Lines kept from the end of each task's stderr.
#' @return Unique error strings, newest run only; character(0) when none failed.
#' @noRd
failed_task_errors <- function(log, work, n = 2L) {
  if (is.null(log) || !nzchar(log %||% "") || !file.exists(log)) return(character(0))
  x <- readLines(log, warn = FALSE)
  x <- x[grepl("terminated with an error", x, fixed = TRUE)]
  rx <- "\\[([0-9a-f]{2}/[0-9a-f]{6,})\\]"
  x <- x[grepl(rx, x)]
  if (length(x) == 0) return(character(0))
  hashes <- unique(sub(paste0("^.*", rx, ".*$"), "\\1", x))
  errs <- vapply(hashes, function(h) {
    p <- strsplit(h, "/", fixed = TRUE)[[1]]
    d <- list.dirs(file.path(work, p[1]), recursive = FALSE, full.names = TRUE)
    d <- d[startsWith(basename(d), p[2])]
    if (length(d) == 0) return("")
    f <- file.path(d[1], ".command.err")
    if (!file.exists(f)) return("")
    y <- trimws(readLines(f, warn = FALSE))
    y <- y[nzchar(y)]
    if (length(y) == 0) "" else paste(utils::tail(y, n), collapse = " ")
  }, character(1), USE.NAMES = FALSE)
  unique(errs[nzchar(errs)])
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
  work <- file.path(project, "work_ref")
  args <- nf_reference_command(project, "COUNT", man, ref_dir, out,
    api_key = api_key, refseq = refseq, resume = FALSE, work = work)
  # a non-zero status is not fatal: clades that did report are still worth
  # keeping, and the walk treats a missing count as unknown
  status <- runner(project, args, log)
  # ignored failures leave the run status at 0, so the errors are collected
  # whatever the status says
  errs <- tryCatch(failed_task_errors(log, work), error = function(e) character(0))
  if (!identical(status, 0L) && !identical(status, 0)) {
    warning("run_clade_counts: Nextflow COUNT run exited with status ", status,
      if (length(errs) > 0) paste0(": ", errs[1]), call. = FALSE)
  }
  res <- merge_counts(ref_dir, read_count_files(out))
  # carried so the caller's failure message can quote the run that produced it
  attr(res, "log") <- log
  attr(res, "errors") <- errs
  res
}
