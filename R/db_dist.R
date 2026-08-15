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

#' Base URL of the panel release assets
#' @noRd
RP_DB_RELEASE_BASE <- "https://github.com/dmacguigan/RiboPilot/releases/download"

#' Default panel version
#' @noRd
RP_DB_VERSION <- "2026.08.1"

#' Download one asset and check it against its published .sha256
#' @noRd
fetch_verified <- function(base_url, tag, asset, dir) {
  url <- paste(base_url, tag, asset, sep = "/")
  dest <- file.path(dir, asset)
  if (utils::download.file(url, dest, mode = "wb", quiet = TRUE) != 0 || !file.exists(dest)) {
    stop("Could not download ", url, call. = FALSE)
  }
  sum_dest <- paste0(dest, ".sha256")
  if (utils::download.file(paste0(url, ".sha256"), sum_dest, mode = "wb", quiet = TRUE) != 0) {
    stop("Could not download the checksum for ", asset, call. = FALSE)
  }
  line <- readLines(sum_dest, warn = FALSE)[1]
  want <- regmatches(line, regexpr("[0-9a-fA-F]{64}", line))
  got <- rp_sha256(dest)
  if (length(want) != 1 || is.na(got)) {
    stop("Could not verify the checksum of ", asset, call. = FALSE)
  }
  if (!identical(tolower(want), got)) {
    stop("Checksum mismatch for ", asset, ": the download is corrupt or truncated.",
      call. = FALSE)
  }
  dest
}

#' Download and install a RiboPilot BLAST panel
#'
#' Extraction is staged in a sibling directory and moved into place only after
#' every asset verifies, so an interrupted or corrupt download never leaves a
#' half-panel that blast_db_exists() would call usable.
#'
#' @param version Panel version, matching release tag `db-<version>`.
#' @param dest Parent directory; the panel lands in `<dest>/<version>`.
#' @param decoys Also fetch the NCBI decoy databases (adds roughly 110 MB).
#' @param base_url Release download base; overridable for testing.
#' @param force Re-download even when the panel is already present.
#' @return (invisibly) the panel prefix.
#' @export
install_ribopilot_db <- function(version = RP_DB_VERSION,
                                 dest = ribopilot_db_dir(),
                                 decoys = TRUE,
                                 base_url = RP_DB_RELEASE_BASE,
                                 force = FALSE) {
  dest <- path.expand(dest)
  final <- file.path(dest, version)
  prefix <- file.path(final, "ribopilot_rrna")
  if (!force && blast_db_exists(prefix)) {
    message("Panel ", version, " is already installed at ", final)
    return(invisible(prefix))
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (unname(file.access(dest, 2)) != 0) {
    stop(dest, " is not writable. Set RIBOPILOT_DB_DIR to a writable location.",
      call. = FALSE)
  }
  stage <- file.path(dest, paste0(".stage-", version))
  unlink(stage, recursive = TRUE)
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  tag <- paste0("db-", version)
  assets <- sprintf("ribopilot_rrna-%s.tar.gz", version)
  if (isTRUE(decoys)) assets <- c(assets, sprintf("decoys-%s.tar.gz", version))
  for (a in assets) {
    message("Downloading ", a, " ...")
    tgz <- fetch_verified(base_url, tag, a, stage)
    utils::untar(tgz, exdir = stage)
    unlink(c(tgz, paste0(tgz, ".sha256")))
  }
  if (!blast_db_exists(file.path(stage, "ribopilot_rrna"))) {
    stop("The downloaded archive did not contain a usable panel.", call. = FALSE)
  }
  unlink(final, recursive = TRUE)
  if (!file.rename(stage, final)) {
    stop("Could not move the panel into ", final, call. = FALSE)
  }
  message("Panel ", version, " installed at ", final)
  invisible(prefix)
}
