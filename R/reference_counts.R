# How many nuclear rRNA records GenBank holds for a clade.
#
# This is the evidence the clade suggestion walks on: a clade nobody has
# sequenced makes a useless seed database, and one covering half of Metazoa
# makes an unspecific one. The count comes from the same query the builder will
# run, so the number shown is the number that will be downloaded.

#' GenBank organism term for a clade name
#' @noRd
clade_query <- function(name) sprintf("\"%s\"[Organism]", trimws(name))

#' Pull the record count out of the builder script's output
#' @param lines Character vector of script stdout.
#' @return integer, or NA when the script did not report a count.
#' @noRd
parse_record_count <- function(lines) {
  hit <- grep("^Found [0-9]+ matching", lines, value = TRUE)
  if (length(hit) == 0) return(NA_integer_)
  suppressWarnings(as.integer(sub("^Found ([0-9]+) .*$", "\\1", hit[1])))
}

#' Count GenBank nuclear rRNA records for a query
#'
#' Runs the shipped builder script in `--count-only` mode, so the count matches
#' exactly what a build of the same query would fetch. Requires Entrez Direct on
#' PATH (pass `conda_env` when the app is not running inside that environment).
#'
#' @param query Advanced GenBank search term, e.g. `'"Asteroidea"[Organism]'`.
#' @param api_key Optional NCBI API key.
#' @param refseq Restrict to RefSeq records.
#' @param conda_env Optional conda env (name or path) providing Entrez Direct.
#' @return integer record count, or `NA_integer_` if the count could not be read.
#' @noRd
genbank_rrna_count <- function(query, api_key = NULL, refseq = FALSE, conda_env = NULL) {
  script <- app_sys("ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh")
  if (!nzchar(script) || !file.exists(script)) return(NA_integer_)
  if (!is.null(conda_env) && nzchar(conda_env)) {
    bin <- tryCatch(resolve_conda_bin(conda_env), error = function(e) NULL)
    if (!is.null(bin)) {
      old <- Sys.getenv("PATH")
      Sys.setenv(PATH = paste(bin, old, sep = .Platform$path.sep))
      on.exit(Sys.setenv(PATH = old), add = TRUE)
    }
  }
  args <- c(shQuote(script), "--count-only")
  if (isTRUE(refseq)) args <- c(args, "--refseq")
  if (!is.null(api_key) && nzchar(api_key)) args <- c(args, "--api-key", shQuote(api_key))
  args <- c(args, shQuote(query), "count_probe")
  out <- suppressWarnings(tryCatch(system2("bash", args, stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)))
  parse_record_count(out)
}

#' Path to a reference folder's record-count cache
#' @noRd
ref_counts_path <- function(dir) file.path(path.expand(dir), "counts.tsv")

#' Record-count lookup with an on-disk cache
#'
#' Counts change slowly and the walk revisits the same high ranks for every
#' project, so they are cached per reference folder. Failures are not cached: a
#' missing Entrez Direct today should not poison the folder forever.
#'
#' @param dir Reference folder holding the cache.
#' @param count_fn Function(query, api_key, refseq, conda_env) -> integer. Injected
#'   so tests never reach the network.
#' @param ... Passed to `count_fn` (api_key, refseq, conda_env).
#' @return function(taxid, name) -> integer count (NA when unavailable).
#' @noRd
cached_count_fn <- function(dir, count_fn = genbank_rrna_count, ...) {
  extra <- list(...)
  p <- ref_counts_path(dir)
  read_cache <- function() {
    if (!file.exists(p)) {
      return(data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE))
    }
    d <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || !all(c("taxid", "n_records") %in% names(d))) {
      return(data.frame(taxid = integer(0), n_records = integer(0), stringsAsFactors = FALSE))
    }
    d$taxid <- suppressWarnings(as.integer(d$taxid))
    d$n_records <- suppressWarnings(as.integer(d$n_records))
    d
  }
  function(taxid, name) {
    taxid <- suppressWarnings(as.integer(taxid))
    cache <- read_cache()
    i <- match(taxid, cache$taxid)
    if (!is.na(i) && !is.na(cache$n_records[i])) return(cache$n_records[i])
    n <- do.call(count_fn, c(list(clade_query(name)), extra))
    n <- suppressWarnings(as.integer(n))
    if (is.na(n)) return(NA_integer_)
    dir.create(path.expand(dir), recursive = TRUE, showWarnings = FALSE)
    cache <- cache[cache$taxid != taxid | is.na(cache$taxid), , drop = FALSE]
    cache <- rbind(cache, data.frame(taxid = taxid, n_records = n, stringsAsFactors = FALSE))
    utils::write.table(cache, p, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    n
  }
}
