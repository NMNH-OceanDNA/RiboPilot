# Build a per-clade GetOrganelle nuclear-rDNA reference database from a GenBank
# taxonomic query. Wraps the shipped builder (inst/ref_dbs/getOrganelle/
# GenBankDownload/download_rrna.sh + parse_rrna_fasta.py). Databases are written
# to a user-chosen directory OUTSIDE any project so one clade DB serves every
# project for that clade (mirrors MitoPilot's per-clade reference builder).

#' Default directory for user-built reference databases
#'
#' A per-user, cross-project location (override with the `RIBOPILOT_REF_DIR`
#' environment variable). Kept outside project directories so a clade DB built
#' once is reused by every project. Uses [tools::R_user_dir()].
#' @export
ribopilot_ref_dir <- function() {
  d <- Sys.getenv("RIBOPILOT_REF_DIR", "")
  if (!nzchar(d)) d <- file.path(tools::R_user_dir("RiboPilot", which = "data"), "ref_dbs")
  path.expand(d)
}

#' Build a clade nuclear-rDNA reference database (GetOrganelle seeds + labels)
#'
#' Downloads 18S / 5.8S / 28S / ITS records for a GenBank taxonomic query and
#' builds the GetOrganelle seed + label FASTAs used by [new_project()]
#' (`custom_seeds_db` / `custom_labels_db`). Run this once per clade, before
#' running the pipeline; point every project for that clade at the result.
#'
#' Requires Entrez Direct + biopython (in the RiboPilot conda env). Pass
#' `conda_env` to put them on PATH.
#'
#' This is the CLI entry point for building a database outside the app. The
#' RiboPilot app's Reference tab builds databases through its own Nextflow
#' workflow instead and does not call this function; use it directly for
#' scripted or headless use, or to build a database ahead of time.
#'
#' @param query Advanced GenBank search term, e.g. `'"Gastropoda"[Organism]'`
#'   or `'"Percidae"[Organism]'`.
#' @param clade Short name for the output files, e.g. `"gastropoda"`. Files are
#'   written as `<out_dir>/seeds/<clade>_rrna_seeds.fasta` and
#'   `<out_dir>/labels/<clade>_rrna_labels.fasta`.
#' @param out_dir Directory to store databases (default [ribopilot_ref_dir()],
#'   outside any project so it is reused across projects).
#' @param api_key Optional NCBI API key (raises the request-rate limit).
#' @param refseq Restrict to RefSeq records.
#' @param count_only Only print the matching-record count and return (no
#'   download).
#' @param conda_env Optional conda env (name or path) providing Entrez Direct +
#'   biopython; its `bin/` is prepended to PATH for the run.
#' @param overwrite Rebuild even if the seed/label FASTAs already exist. When
#'   `FALSE` (default) an existing DB is reused (no download) and its paths
#'   returned.
#' @return (invisibly) a list with `seeds` and `labels` file paths (both `NA`
#'   for `count_only`).
#' @export
build_reference_db <- function(query, clade, out_dir = ribopilot_ref_dir(),
                               api_key = NULL, refseq = FALSE, count_only = FALSE,
                               conda_env = NULL, overwrite = FALSE) {
  stopifnot(is.character(query), nzchar(query), is.character(clade), nzchar(clade))
  script <- app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
  if (!nzchar(script) || !file.exists(script)) {
    stop("Reference-DB builder script not found in the installed package.", call. = FALSE)
  }

  if (!is.null(conda_env)) {
    bin <- resolve_conda_bin(conda_env)
    old_path <- Sys.getenv("PATH")
    Sys.setenv(PATH = paste(bin, old_path, sep = .Platform$path.sep))
    on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  }

  out_dir <- path.expand(out_dir)
  paths <- reference_db_paths(out_dir, clade)
  prefix <- paths$prefix
  seeds <- paths$seeds
  labels <- paths$labels

  if (isTRUE(count_only)) {
    args <- c(shQuote(script), "--count-only")
    if (isTRUE(refseq)) args <- c(args, "--refseq")
    if (!is.null(api_key) && nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))
    args <- c(args, shQuote(query), shQuote(prefix))
    system2("bash", args, stdin = "/dev/null")
    return(invisible(list(seeds = NA_character_, labels = NA_character_)))
  }

  if (!overwrite && file.exists(seeds) && file.exists(labels)) {
    message("Reference DB for '", clade, "' already exists in ", out_dir,
      " (use overwrite = TRUE to rebuild). Reusing.")
    return(invisible(list(seeds = seeds, labels = labels)))
  }

  dir.create(file.path(out_dir, "seeds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(out_dir, "labels"), recursive = TRUE, showWarnings = FALSE)

  # Build into a temp staging dir (the script writes flat <prefix>_{seeds,labels}
  # + <prefix>.raw.fasta), then place into the seeds/ and labels/ subdirs.
  stage <- file.path(tempdir(), paste0("rp_refdb_", clade, "_", as.integer(Sys.time())))
  dir.create(stage, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  args <- c(shQuote(script), "--outdir", shQuote(stage))
  if (isTRUE(refseq)) args <- c(args, "--refseq")
  if (!is.null(api_key) && nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))
  args <- c(args, shQuote(query), shQuote(prefix))
  status <- system2("bash", args, stdin = "/dev/null")
  if (status != 0) {
    stop("Reference-DB build failed (Entrez Direct / biopython on PATH? try conda_env=).",
      call. = FALSE)
  }

  st_seeds <- file.path(stage, paste0(prefix, "_seeds.fasta"))
  st_labels <- file.path(stage, paste0(prefix, "_labels.fasta"))
  if (!file.exists(st_seeds) || !file.exists(st_labels)) {
    stop("Builder finished but seed/label FASTAs were not produced.", call. = FALSE)
  }
  file.copy(st_seeds, seeds, overwrite = TRUE)
  file.copy(st_labels, labels, overwrite = TRUE)

  message("Built reference DB for '", clade, "':")
  message("  seeds:  ", seeds)
  message("  labels: ", labels)
  message("Pass these to new_project(custom_seeds_db=, custom_labels_db=).")
  invisible(list(seeds = seeds, labels = labels))
}
