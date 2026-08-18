# Getting a reference panel onto disk, and knowing which one it is.

#' Pull the digest out of one line of sha256sum/shasum/openssl output
#'
#' Position matters: sha256sum/shasum print the hex first ("hex  path"), so
#' that shape is anchored to the line start and bounded by a separator or
#' end-of-string (a 65-hex run must reject, not truncate to a spurious 64);
#' openssl prints the path first ("SHA2-256(path)= hex"), so that shape is
#' anchored to the "= " right before the hex. An unanchored search would
#' match a hex-looking path component instead of the real digest.
#' @return 64-char lowercase hex, or NA_character_ when unrecognised.
#' @noRd
parse_sha256_line <- function(line) {
  if (length(line) != 1 || is.na(line)) return(NA_character_)
  line <- sub("^\\s+", "", line)
  m <- regmatches(line, regexpr("^[0-9a-fA-F]{64}(\\s|$)", line))
  if (length(m) == 1) return(tolower(substr(m, 1, 64)))
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
#'
#' `YYYY.MM.DD`, the date the panel was harvested. Dotted numbers only:
#' `numeric_version()` returns NA for anything else, and discover_panel() then
#' cannot see the panel at all. The date is the version because the source grows
#' by roughly 370 records a day, so two builds are never the same content; a
#' hand-bumped counter just gets forgotten, which is how two different panels
#' both came to be called 2026.08.1.
#' @noRd
RP_DB_VERSION <- "2026.08.18"

#' Fetch one URL, or stop naming the asset and the reason
#'
#' download.file() signals a condition rather than returning nonzero, so a
#' failure has to be caught to be reported at all. options(timeout) is a cap on
#' the whole transfer rather than on idle time, and its 60s default cannot carry
#' a 76 MB asset on any link slower than 1.3 MB/s, so it is raised here.
#' @noRd
download_or_stop <- function(url, dest, what, quiet = TRUE) {
  old <- options(timeout = max(getOption("timeout"), 3600))
  on.exit(options(old), add = TRUE)
  why <- ""
  note <- function(c) {
    why <<- conditionMessage(c)
    FALSE
  }
  ok <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = quiet) == 0,
    error = note, warning = note)
  if (!isTRUE(ok) || !file.exists(dest)) {
    stop("Could not download ", what, if (nzchar(why)) paste0(": ", why),
      ". A slow or interrupted connection is the usual cause; try again.",
      call. = FALSE)
  }
  dest
}

#' Download one asset and check it against its published .sha256
#' @noRd
fetch_verified <- function(base_url, tag, asset, dir) {
  url <- paste(base_url, tag, asset, sep = "/")
  # Progress on the asset itself: it is the part worth waiting on.
  dest <- download_or_stop(url, file.path(dir, asset), url, quiet = FALSE)
  sum_dest <- download_or_stop(paste0(url, ".sha256"), paste0(dest, ".sha256"),
    paste0("the checksum for ", asset))
  want <- parse_sha256_line(readLines(sum_dest, warn = FALSE)[1])
  got <- rp_sha256(dest)
  if (is.na(want) || is.na(got)) {
    stop("Could not verify the checksum of ", asset, call. = FALSE)
  }
  if (!identical(want, got)) {
    stop("Checksum mismatch for ", asset, ": the download is corrupt or truncated.",
      call. = FALSE)
  }
  dest
}

#' Which pieces of a whole panel are absent from a directory
#'
#' blast_db_exists() is satisfied by an index file alone, so a partial unpack
#' would otherwise install and then be searched as if it were whole. Naming the
#' absentees keeps the refusal diagnosable, for instance when a future panel
#' ships as .nal volumes rather than the single-volume layout expected here.
#' The extension set is PANEL_REQUIRED_EXT, the same one discovery uses, so an
#' installed panel cannot be one discovery then refuses to see.
#' @return character(0) when the directory holds a whole panel.
#' @noRd
panel_missing_parts <- function(dir) {
  want <- paste0("ribopilot_rrna.", PANEL_REQUIRED_EXT)
  miss <- want[!file.exists(file.path(dir, want))]
  if (is.null(read_panel_manifest(dir))) miss <- c(miss, "manifest.json")
  miss
}

#' Reclaim staging and set-aside directories left by interrupted installs
#'
#' Each abandoned one can hold a few hundred MB. Age stands in for a liveness
#' check so a concurrent install in another process is never swept out from
#' under itself.
#' ponytail: age heuristic; use a real liveness check if installs ever run long.
#' @return (invisibly) the paths removed.
#' @noRd
sweep_stale_installs <- function(dest, version, max_age_days = 1) {
  hits <- list.files(dest, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  hits <- hits[dir.exists(hits)]
  nm <- basename(hits)
  hits <- hits[startsWith(nm, paste0(".stage-", version, "-")) |
    startsWith(nm, paste0(version, ".old-"))]
  age <- difftime(Sys.time(), file.mtime(hits), units = "days")
  hits <- hits[!is.na(age) & age > max_age_days]
  if (length(hits)) {
    message("Reclaiming space from an interrupted install: ",
      paste(basename(hits), collapse = ", "))
    unlink(hits, recursive = TRUE)
  }
  invisible(hits)
}

#' Download and install a RiboPilot BLAST panel
#'
#' Extraction is staged in a sibling directory and moved into place only after
#' every asset verifies, so an interrupted or corrupt download never leaves a
#' half-panel that blast_db_exists() would call usable. An already-installed
#' panel is only deleted once the replacement is in place.
#'
#' @param version Panel version, matching release tag `db-<version>`.
#' @param dest Parent directory; the panel lands in `<dest>/<version>`.
#' @param decoys Also fetch the NCBI decoy databases (adds roughly 110 MB).
#' @param base_url Release download base; overridable for testing.
#' @param force Re-download even when the panel is already present.
#' @return (invisibly) the panel prefix.
#' @noRd
install_ribopilot_db <- function(version = RP_DB_VERSION,
                                 dest = ribopilot_db_dir(),
                                 decoys = TRUE,
                                 base_url = RP_DB_RELEASE_BASE,
                                 force = FALSE) {
  dest <- path.expand(dest)
  final <- file.path(dest, version)
  prefix <- file.path(final, "ribopilot_rrna")
  # Refused up front, not after the download: blastn cannot express a
  # whitespace path, so validate_blastdb_prefix() would reject the result.
  if (grepl("\\s", prefix)) {
    stop(dest, " contains whitespace; blastn cannot read a panel there. ",
      "Set RIBOPILOT_DB_DIR to a path without spaces.", call. = FALSE)
  }
  if (!force && blast_db_exists(prefix)) {
    message("Panel ", version, " is already installed at ", final)
    return(invisible(prefix))
  }
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  if (unname(file.access(dest, 2)) != 0) {
    stop(dest, " is not writable. Set RIBOPILOT_DB_DIR to a writable location.",
      call. = FALSE)
  }
  sweep_stale_installs(dest, version)
  # The pid keeps two concurrent installs of one version off each other's stage.
  stage <- file.path(dest, paste0(".stage-", version, "-", Sys.getpid()))
  unlink(stage, recursive = TRUE)
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  tag <- paste0("db-", version)
  assets <- sprintf("ribopilot_rrna-%s.tar.gz", version)
  if (isTRUE(decoys)) assets <- c(assets, sprintf("decoys-%s.tar.gz", version))
  for (a in assets) {
    message("Downloading ", a, " ...")
    tgz <- fetch_verified(base_url, tag, a, stage)
    st <- tryCatch(utils::untar(tgz, exdir = stage),
      error = function(e) 1L, warning = function(w) 1L)
    if (!isTRUE(st == 0)) {
      stop("Could not unpack ", a, ": the archive is damaged, or the disk is full.",
        call. = FALSE)
    }
    unlink(c(tgz, paste0(tgz, ".sha256")))
  }
  miss <- panel_missing_parts(stage)
  if (length(miss)) {
    stop("The downloaded archive did not contain a usable panel; it is missing ",
      paste(miss, collapse = ", "), ".", call. = FALSE)
  }

  # Set the existing panel aside rather than delete it, so a failed rename cannot
  # leave a gutted directory that blast_db_exists() still calls usable. A name
  # nothing else holds: a leftover .old- may be a crashed run's only copy, and
  # pids get recycled.
  old <- tempfile(paste0(version, ".old-"), tmpdir = dest)
  had_old <- file.exists(final)
  if (had_old && !file.rename(final, old)) {
    stop("Could not set the panel at ", final, " aside; it is unchanged.", call. = FALSE)
  }
  if (!file.rename(stage, final)) {
    if (had_old && !file.rename(old, final)) {
      stop("Could not move the panel into ", final, ", and the previous panel is now ",
        "at ", old, "; move it back by hand.", call. = FALSE)
    }
    stop("Could not move the panel into ", final, call. = FALSE)
  }
  if (had_old && unlink(old, recursive = TRUE) != 0) {
    warning("The replaced panel is still at ", old, "; delete it to reclaim the space.",
      call. = FALSE)
  }
  message("Panel ", version, " installed at ", final)
  message("Reference database path to use in RiboPilot: ", prefix)
  invisible(prefix)
}

#' Record which panel a run searched against
#'
#' Written when a run launches so a result can be traced to the database that
#' produced it. Silent and non-fatal on a project predating the column.
#' @return (invisibly) the recorded version, or "" when there is nothing to record.
#' @noRd
record_blastdb_provenance <- function(con, prefix) {
  prefix <- trimws(prefix %||% "")
  if (!nzchar(prefix)) return(invisible(""))
  m <- read_panel_manifest(dirname(path.expand(prefix)))
  if (is.null(m)) return(invisible(""))
  v <- paste(m$db_name, m$db_version)
  ok <- tryCatch({
    DBI::dbExecute(con, "UPDATE reference_config SET blastdb_version = ? WHERE rowid_key = 1",
      params = list(v))
    TRUE
  }, error = function(e) FALSE)
  # The panel usually lives inside the container, where the app cannot read its
  # manifest, so the run writes the display line out for the options modal.
  tryCatch(DBI::dbExecute(con,
    "UPDATE reference_config SET blastdb_info = ? WHERE rowid_key = 1",
    params = list(panel_summary(m))), error = function(e) NULL)
  invisible(if (ok) as.character(m$db_version) else "")
}

#' Thousands-separated count, or "" for anything unusable
#' @noRd
fmt_count <- function(x) {
  if (length(x) != 1) return("")
  n <- suppressWarnings(as.numeric(x))
  if (is.na(n)) return("")
  formatC(n, format = "d", big.mark = ",")
}

#' One-line description of a panel, from its manifest
#'
#' Every part after the name and version is optional: an older manifest carries
#' no counts block, and a note that reads "NA sequences" is worse than a short one.
#' @return character(1); "" when there is no manifest.
#' @noRd
panel_summary <- function(m) {
  if (!is.list(m)) return("")
  out <- trimws(paste(m$db_name, m$db_version))
  built <- substr(m$built_utc %||% "", 1, 10)
  if (nzchar(built)) out <- paste0(out, ", built ", built)
  n <- fmt_count(m$counts$records_after_dedup)
  if (nzchar(n)) out <- paste0(out, ", ", n, " sequences")
  tx <- fmt_count(m$counts$distinct_taxids)
  if (nzchar(tx)) out <- paste0(out, " from ", tx, " taxa")
  mb <- suppressWarnings(as.numeric(m$counts$total_bases))
  if (length(mb) == 1 && !is.na(mb)) {
    out <- paste0(out, ", ", fmt_count(round(mb / 1e6)), " Mbp")
  }
  out
}

#' The modal's panel line: the manifest if the app can see it, else the last run's
#' @param m Panel manifest, or NULL when the panel is not visible here.
#' @param recorded `reference_config.blastdb_info`, written by the last local run.
#' @noRd
panel_note_text <- function(m, recorded = "") {
  s <- panel_summary(m)
  if (!nzchar(s)) s <- trimws(recorded %||% "")
  if (!nzchar(s)) {
    return(paste0("Bundled BLAST reference database: which database this ",
      "project searches is recorded the first time a local search runs."))
  }
  paste0("Bundled BLAST reference database ", s, ".")
}

#' The image named by a project config, or "" when it names none
#' @noRd
config_container_image <- function(cfg_lines) {
  hit <- grep("^\\s*container\\s*=", cfg_lines, value = TRUE, perl = TRUE)
  if (length(hit) < 1) return("")
  v <- sub("^\\s*container\\s*=\\s*", "", hit[1], perl = TRUE)
  # one capture, not two alternations: R's engine prefers the longest match, so
  # a "strip a leading quote OR a trailing quote and the rest" pattern ate the
  # whole string
  trimws(sub("^['\"]?([^'\"]*).*$", "\\1", trimws(v)))
}

#' Ask the project's container image to describe the database it carries
#'
#' The database lives inside the image, where the app cannot read its manifest,
#' so on the container path this is the only way to describe it before a run has
#' happened. One short call, cached by the caller.
#' @return character(1); "" when there is no image, no engine, or no answer.
#' @noRd
image_panel_summary <- function(project = ".") {
  cfg <- file.path(project, ".config")
  if (!file.exists(cfg)) return("")
  lines <- readLines(cfg, warn = FALSE)
  engine <- config_container_engine(lines)
  img <- config_container_image(lines)
  if (is.na(engine) || !nzchar(img)) return("")
  bin <- if (identical(engine, "docker")) "docker" else {
    w <- Sys.which(c("singularity", "apptainer"))
    w <- w[nzchar(w)]
    if (length(w) < 1) return("") else names(w)[1]
  }
  if (!nzchar(Sys.which(bin))) return("")
  code <- paste0("p <- RiboPilot:::discover_panel(); ",
    "if (!is.null(p)) cat(RiboPilot:::panel_summary(",
    "RiboPilot:::read_panel_manifest(dirname(p))))")
  args <- if (identical(engine, "docker")) {
    c("run", "--rm", img, "Rscript", "-e", shQuote(code))
  } else {
    c("exec", img, "Rscript", "-e", shQuote(code))
  }
  # timeout: a missing image makes docker try to pull one, which on a slow link
  # would hang the modal open
  out <- tryCatch(suppressWarnings(
    system2("timeout", c("60", bin, args), stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(out) < 1) "" else trimws(paste(out, collapse = " "))
}

#' What the last local run recorded about the panel; "" when nothing has run
#' @noRd
blastdb_info <- function(con) {
  v <- tryCatch(DBI::dbGetQuery(con,
    "SELECT blastdb_info FROM reference_config WHERE rowid_key = 1")$blastdb_info,
    error = function(e) NULL)
  if (length(v) < 1 || is.na(v[1])) "" else trimws(v[1])
}

#' Panel line for a project, preferring a panel this session can actually read
#'
#' `prefix` is the project's configured panel; empty means "whatever is
#' discovered", which on the container path is nothing the app can see.
#' @noRd
project_panel_note <- function(con, prefix = "", project = ".",
                               probe = image_panel_summary) {
  prefix <- trimws(prefix %||% "")
  if (!nzchar(prefix)) prefix <- discover_panel() %||% ""
  m <- if (nzchar(prefix)) read_panel_manifest(dirname(path.expand(prefix))) else NULL
  recorded <- blastdb_info(con)
  # nothing here and nothing recorded: the database is inside the image, so ask
  # it once and keep the answer, rather than showing a placeholder until the
  # first run finishes
  if (is.null(m) && !nzchar(recorded)) {
    recorded <- tryCatch(probe(project), error = function(e) "")
    if (nzchar(recorded)) {
      tryCatch(DBI::dbExecute(con,
        "UPDATE reference_config SET blastdb_info = ? WHERE rowid_key = 1",
        params = list(recorded)), error = function(e) NULL)
    }
  }
  panel_note_text(m, recorded)
}

#' The recorded panel, for display; "" when none
#' @noRd
blastdb_provenance <- function(con) {
  v <- tryCatch(DBI::dbGetQuery(con,
    "SELECT blastdb_version FROM reference_config WHERE rowid_key = 1")$blastdb_version,
    error = function(e) NULL)
  if (length(v) < 1 || is.na(v[1]) || !nzchar(trimws(v[1]))) "" else trimws(v[1])
}
