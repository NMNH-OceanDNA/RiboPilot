# Initialize a RiboPilot test project: fetch the bundled test samples' reads,
# write the mapping, and initialize the project. Mirrors MitoPilot's
# new_test_project, but fetches reads with an ENA-then-NCBI-SRA fallback (the ENA
# fastq mirror is firewall-blocked on some networks, e.g. Smithsonian, where NCBI
# SRA works). rDNA is present in the full reads, so there is no pre-filtered set.

#' Build the ENA fastq download URL for a run accession + mate
#' @noRd
.ena_fastq_url <- function(acc, mate) {
  pre <- substr(acc, 1, 6)
  n <- nchar(acc)
  sub <- if (n <= 9) "" else {
    d <- n - 9  # 1 extra digit -> 00X/, 2 -> 0XX/, 3 -> XXX/
    paste0(strrep("0", max(0, 3 - d)), substr(acc, n - d + 1, n), "/")
  }
  sprintf("https://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s%s/%s_%d.fastq.gz", pre, sub, acc, acc, mate)
}

#' TRUE if `fn` exists and is a valid (non-empty) gzip file
#' @noRd
.gz_ok <- function(fn) {
  file.exists(fn) && file.info(fn)$size > 0 &&
    isTRUE(tryCatch(
      system2("gzip", c("-t", shQuote(fn)), stdout = FALSE, stderr = FALSE) == 0L,
      error = function(e) FALSE))
}

#' Download one ENA mate to `dest`, retrying until the gzip validates
#' @noRd
.fetch_ena_mate <- function(acc, mate, dest, max_tries = 3) {
  url <- .ena_fastq_url(acc, mate)
  for (i in seq_len(max_tries)) {
    # -f: fail (no HTML error body) on 403/404 so .gz_ok correctly rejects it
    try(system2("curl", c("-fsSL", "--retry", "2", shQuote(url), "-o", shQuote(dest)),
      stdout = FALSE, stderr = FALSE), silent = TRUE)
    if (.gz_ok(dest)) return(TRUE)
  }
  FALSE
}

#' Fetch a run's paired reads via NCBI SRA (fasterq-dump), gzip to R1/R2 names
#' @noRd
.fetch_sra_pair <- function(acc, r1_dest, r2_dest, dir, sra_bin = "", scratch = tempdir()) {
  fqd <- if (nzchar(sra_bin)) file.path(sra_bin, "fasterq-dump") else Sys.which("fasterq-dump")
  if (!nzchar(fqd) || !file.exists(fqd)) return(FALSE)
  dir.create(scratch, showWarnings = FALSE, recursive = TRUE)
  try(system2(fqd, c("--split-files", "-e", "4", "-O", shQuote(dir),
    "-t", shQuote(scratch), shQuote(acc)), stdout = FALSE, stderr = FALSE), silent = TRUE)
  f1 <- file.path(dir, paste0(acc, "_1.fastq"))
  f2 <- file.path(dir, paste0(acc, "_2.fastq"))
  if (!file.exists(f1) || !file.exists(f2)) return(FALSE)
  gz <- if (nzchar(Sys.which("pigz"))) "pigz" else "gzip"
  system2(gz, c("-f", shQuote(f1), shQuote(f2)), stdout = FALSE, stderr = FALSE)
  ok <- file.rename(paste0(f1, ".gz"), r1_dest) && file.rename(paste0(f2, ".gz"), r2_dest)
  isTRUE(ok) && .gz_ok(r1_dest) && .gz_ok(r2_dest)
}

#' Fetch one sample's paired reads: ENA first, NCBI SRA fallback
#' @noRd
.download_test_reads <- function(acc, r1_name, r2_name, dir, sra_bin = "", scratch = tempdir()) {
  r1 <- file.path(dir, r1_name)
  r2 <- file.path(dir, r2_name)
  if (.gz_ok(r1) && .gz_ok(r2)) {
    message("    already present, skipping")
    return(invisible(TRUE))
  }
  if (.fetch_ena_mate(acc, 1, r1) && .fetch_ena_mate(acc, 2, r2)) {
    return(invisible(TRUE))
  }
  message("    ENA unavailable; trying NCBI SRA (fasterq-dump) ...")
  if (.fetch_sra_pair(acc, r1, r2, dir, sra_bin, scratch)) {
    return(invisible(TRUE))
  }
  stop("Could not download reads for ", acc,
    " from ENA or NCBI SRA. If ENA is blocked on your network, ensure sra-tools ",
    "is on PATH or pass sra_bin=.", call. = FALSE)
}

#' Stage each mapping row's reads into `dir`: bundled (Local) rows are copied
#' from the package, the rest downloaded.
#' @noRd
.stage_test_reads <- function(mapping, dir, local_dir = app_sys("test_data"),
                              sra_bin = "", scratch = tempdir()) {
  local <- (mapping$Local %||% rep(FALSE, nrow(mapping))) %in% TRUE
  for (i in seq_len(nrow(mapping))) {
    message(sprintf("  [%d/%d] %s - %s", i, nrow(mapping), mapping$ID[i], mapping$Taxon[i]))
    if (local[i]) {
      for (f in c(mapping$R1[i], mapping$R2[i])) file.copy(file.path(local_dir, f), file.path(dir, f), overwrite = TRUE)
    } else {
      .download_test_reads(mapping$ID[i], mapping$R1[i], mapping$R2[i], dir, sra_bin, scratch)
    }
  }
  invisible(NULL)
}

#' Keep the first `n` downloadable samples; bundled (Local) rows are always kept.
#' @noRd
.select_test_samples <- function(mapping, n = Inf) {
  local <- (mapping$Local %||% rep(FALSE, nrow(mapping))) %in% TRUE
  if (!is.finite(n)) return(mapping)
  rbind(utils::head(mapping[!local, , drop = FALSE], n), mapping[local, , drop = FALSE])
}

#' Initialize a RiboPilot test project
#'
#' Sets up a ready-to-run test project from the bundled 10-sample test set (USNM
#' fish genome skims, study SRP313821): downloads each sample's paired reads,
#' writes the mapping file, and initializes the project database + config.
#'
#' Reads are fetched from ENA when reachable, falling back to NCBI SRA
#' (`fasterq-dump`, which must be on `PATH` or given via `sra_bin`) on networks
#' where the ENA fastq mirror is blocked. The full set is several GB and takes a
#' while; use `n` to fetch fewer samples for a quicker setup. Already-downloaded,
#' valid read files are skipped, so an interrupted run resumes cleanly.
#'
#' @param path Project directory (created if missing; default current dir).
#' @param n Number of downloaded samples to include (default all 10). The
#'   bundled synthetic samples CHIM00 and CHIM01 are always included.
#' @param executor Executor / config profile (default "conda").
#' @param container Container image for the pipeline.
#' @param min_depth Minimum trimmed paired reads to proceed with assembly.
#' @param sra_bin Directory holding `fasterq-dump` (for the SRA fallback);
#'   default searches `PATH`.
#' @param Rproj Initialize an RStudio project (default TRUE; no effect outside
#'   RStudio).
#' @param force Recreate an existing project database / config (default FALSE).
#' @param ... Passed to [new_project()] (e.g. `custom_seeds_db`).
#' @return (invisibly) the project path.
#' @export
new_test_project <- function(path = ".", n = Inf, executor = "conda",
                             container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
                             min_depth = 500, sra_bin = NULL, Rproj = TRUE, force = FALSE, ...) {
  if (!dir.exists(path)) {
    message("Creating project directory: ", path)
    dir.create(path, recursive = TRUE)
  }
  path <- normalizePath(path)

  mapping <- utils::read.csv(app_sys(file.path("test_data", "mapping_test.csv")),
    stringsAsFactors = FALSE)
  mapping <- .select_test_samples(mapping, n)

  data_dir <- file.path(path, "data")
  dir.create(data_dir, showWarnings = FALSE)
  scratch <- file.path(tempdir(), "ribopilot_sra")
  sra_bin <- sra_bin %||% ""

  message("Fetching ", nrow(mapping), " test sample(s). This downloads several GB ",
    "and may take a while (go grab a coffee)...")
  .stage_test_reads(mapping, data_dir, sra_bin = sra_bin, scratch = scratch)
  unlink(scratch, recursive = TRUE)

  mapping$Local <- NULL
  readr::write_csv(mapping, file.path(path, "mapping.csv"), quote = "none", na = "")

  new_project(
    path = path,
    mapping_fn = file.path(path, "mapping.csv"),
    mapping_id = "ID",
    data_path = file.path(path, "data", ""),
    min_depth = min_depth,
    container = container,
    executor = executor,
    Rproj = Rproj,
    force = force,
    ...
  )
  invisible(path)
}
