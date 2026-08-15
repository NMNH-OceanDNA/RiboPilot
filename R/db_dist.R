# Getting a reference panel onto disk, and knowing which one it is.

#' Pull the digest out of one line of sha256sum/shasum/openssl output
#'
#' Position matters: sha256sum/shasum print the hex first ("hex  path"), so
#' that shape is anchored to the line start; openssl prints the path first
#' ("SHA2-256(path)= hex"), so that shape is anchored to the "= " right
#' before the hex. An unanchored search would match a hex-looking path
#' component instead of the real digest.
#' @return 64-char lowercase hex, or NA_character_ when unrecognised.
#' @noRd
parse_sha256_line <- function(line) {
  if (length(line) != 1 || is.na(line)) return(NA_character_)
  m <- regmatches(line, regexpr("^[0-9a-fA-F]{64}", line))
  if (length(m) == 1) return(tolower(m))
  m <- regmatches(line, regexpr("=\\s*[0-9a-fA-F]{64}\\s*$", line))
  if (length(m) == 1) {
    hex <- regmatches(m, regexpr("[0-9a-fA-F]{64}", m))
    return(tolower(hex))
  }
  NA_character_
}

#' SHA-256 of a file, as lowercase hex
#'
#' tools::sha256sum() arrived in R 4.5.0 and the package floor is 4.4.0, so it
#' is used only when present and a command-line digest is the fallback.
#' @return 64-char lowercase hex, or NA_character_ when unobtainable.
#' @noRd
rp_sha256 <- function(path) {
  if (length(path) != 1 || is.na(path) || !file.exists(path)) return(NA_character_)
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) {
    v <- tryCatch(unname(get("sha256sum", envir = asNamespace("tools"))(path)),
      error = function(e) NA_character_)
    if (length(v) == 1 && !is.na(v)) return(tolower(v))
  }
  cand <- list(c("sha256sum", path), c("shasum", "-a", "256", path),
    c("openssl", "dgst", "-sha256", path))
  for (cc in cand) {
    if (!nzchar(Sys.which(cc[1]))) next
    out <- tryCatch(suppressWarnings(
      system2(cc[1], shQuote(cc[-1]), stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0))
    if (length(out) < 1) next
    hit <- parse_sha256_line(out[1])
    if (!is.na(hit)) return(hit)
  }
  NA_character_
}

#' Read a panel directory's manifest.json
#'
#' Distinct from read_ref_manifest(), which indexes built clade databases.
#' @return Named list, or NULL when there is no usable manifest.
#' @noRd
read_panel_manifest <- function(dir) {
  if (length(dir) != 1 || is.na(dir)) return(NULL)
  p <- file.path(path.expand(dir), "manifest.json")
  if (!file.exists(p)) return(NULL)
  m <- tryCatch(jsonlite::fromJSON(p, simplifyVector = TRUE), error = function(e) NULL)
  if (!is.list(m)) return(NULL)
  ok <- function(x) length(x) == 1 && !is.na(x) && nzchar(as.character(x))
  if (!ok(m$db_version) || !ok(m$db_name)) return(NULL)
  m
}

#' Where downloaded BLAST panels live
#'
#' Honours RIBOPILOT_DB_DIR so a cluster user can put a multi-hundred-MB panel
#' off their home quota, then falls back to a subfolder of the reference dir.
#' @noRd
ribopilot_db_dir <- function() {
  d <- Sys.getenv("RIBOPILOT_DB_DIR", "")
  if (!nzchar(d)) d <- file.path(ribopilot_ref_dir(), "blastdb")
  path.expand(d)
}
