#' If x is `NULL` (or length 0), return y, otherwise return x
#' @noRd
"%||%" <- function(x, y) {
  if (length(x) == 0L) y else x
}

#' If x is `NA`, return y, otherwise return x
#' @noRd
"%|NA|%" <- function(x, y) {
  if (is.na(x)) y else x
}


#' Container engine a Nextflow config actually turns on
#'
#' A config can declare both a `docker {}` and a `singularity {}` block with only
#' one `enabled = true`, so the block's presence says nothing; the flag does.
#' @param cfg_lines Lines of a Nextflow config.
#' @return "docker", "singularity", or `NA_character_`.
#' @noRd
config_container_engine <- function(cfg_lines) {
  txt <- paste(cfg_lines, collapse = "\n")
  if (grepl("docker\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) return("docker")
  if (grepl("(singularity|apptainer)\\s*\\{[^}]*enabled\\s*=\\s*true", txt)) return("singularity")
  NA_character_
}

#' Classify a URL as "ok", "missing" (definitively absent), or "unknown"
#'
#' "unknown" covers the offline / transient-failure case, which must not be
#' treated as absence.
#' @noRd
url_status <- function(u) {
  msg <- NULL
  ok <- tryCatch({
    con <- url(u, open = "rb")
    on.exit(close(con), add = TRUE)
    length(readBin(con, "raw", 1L)) > 0
  },
  error = function(e) {
    msg <<- conditionMessage(e)
    FALSE
  },
  warning = function(w) {
    msg <<- conditionMessage(w)
    FALSE
  })
  if (isTRUE(ok)) return("ok")
  # only a reported HTTP status counts: a proxy or DNS message mentioning "404"
  # is a transient failure, not an absent file
  if (!is.null(msg) && grepl("HTTP status was '40[34]\\b", msg)) {
    return("missing")
  }
  "unknown"
}

#' Warn when a container executor is pointed at a Docker Hub image that is not
#' published
#'
#' Only Docker Hub images are checkable; anything with an explicit registry
#' host, or a locally built image, is left alone. Warns rather than errors
#' because a locally built image is legitimate and invisible to this check.
#' @param container e.g. "macguigand/ribopilot:0.0.0.9000".
#' @noRd
check_container_image <- function(container) {
  if (is.null(container) || is.na(container) || !nzchar(container)) return(invisible(container))
  # Skip explicit registries (host contains a dot or a port) and local images.
  parts <- strsplit(container, "/", fixed = TRUE)[[1]]
  if (length(parts) > 2 || (length(parts) == 2 && grepl("[.:]", parts[1]))) {
    return(invisible(container))
  }
  ref <- strsplit(utils::tail(parts, 1), ":", fixed = TRUE)[[1]]
  if (length(ref) != 2) return(invisible(container))
  # A one-part name (no namespace) is also what a locally built image looks
  # like; skip the Docker Hub check when docker already has it.
  docker <- Sys.which("docker")
  if (length(parts) == 1 && nzchar(docker) &&
      identical(system2(docker, c("image", "inspect", container),
                         stdout = FALSE, stderr = FALSE), 0L)) {
    return(invisible(container))
  }
  ns <- if (length(parts) == 2) parts[1] else "library"
  u <- sprintf("https://hub.docker.com/v2/repositories/%s/%s/tags/%s", ns, ref[1], ref[2])
  if (identical(url_status(u), "missing")) {
    warning(sprintf(paste0(
      "Container image '%s' is not published on Docker Hub, so container ",
      "executors will fail at image pull.\nBuild and push it with ",
      "docker/deploy-dockerhub.sh, or pass an image you can pull via ",
      "new_project(container = ...)."), container), call. = FALSE)
  }
  invisible(container)
}

#' Check a GetOrganelle reference DB is resolvable, warning when it is not
#'
#' A missing seeds/labels DB otherwise surfaces hours into a run as an opaque
#' GetOrganelle failure. Warns rather than errors: the database can be fixed
#' from the app's reference database step, so init should not die on it.
#' @param x Path or URL. @param what "seeds" or "labels".
#' @noRd
check_reference_db <- function(x, what) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(invisible(x))
  hint <- paste0(
    "Set up a database for your samples' group in the app's Reference step, or build one\n",
    "with build_reference_db() and pass db$seeds and db$labels to new_project() as\n",
    "custom_seeds_db and custom_labels_db.")
  if (!grepl("^https?://", x)) {
    if (!file.exists(x)) {
      warning(sprintf("%s reference database not found:\n  %s\n\n%s", what, x, hint),
        call. = FALSE)
    }
    return(invisible(x))
  }
  st <- url_status(x)
  if (identical(st, "missing")) {
    warning(sprintf("%s reference database is not reachable (HTTP 404/403):\n  %s\n\n%s",
      what, x, hint), call. = FALSE)
  }
  if (identical(st, "unknown")) {
    warning(sprintf("Could not verify %s reference database (offline?):\n  %s", what, x),
      call. = FALSE)
  }
  invisible(x)
}

.lst <- function(x, max = 10L) {
  x <- as.character(x)
  if (length(x) <= max) return(paste(x, collapse = ", "))
  paste0(paste(x[seq_len(max)], collapse = ", "), ", ... and ", length(x) - max, " more")
}
