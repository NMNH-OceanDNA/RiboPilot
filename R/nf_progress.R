# Parse Nextflow's ANSI progress board (run with -ansi-log). Each redraw is a
# full snapshot, so key each process row by its name and overwrite in place;
# keep the latest executor line; banner lines before the first process are the
# header, and anything after (warnings, completion) is footer.

nf_proc_rx <- "^\\[(?:[0-9a-f]{2}/[0-9a-f]+|-\\s*)\\]\\s+(\\S+).*$"

#' Most footer lines kept; only the tail is ever rendered, and a long run would
#' otherwise accumulate every warning line for its whole life.
#' @noRd
NF_FOOTER_MAX <- 200L

#' An empty Nextflow progress state
#' @noRd
nf_progress_empty <- function() {
  list(header = character(0), executor = "", proc = list(), more = "", footer = character(0))
}

#' Fold raw Nextflow output lines into a progress state
#' @noRd
nf_progress_parse <- function(state, lines) {
  lines <- cli::ansi_strip(lines)
  for (ln in lines) {
    if (grepl("^executor\\s*>", ln)) {
      state$executor <- ln
      state$more <- ""
    } else if (grepl("^Plus \\d+ more processes", ln)) {
      state$more <- ln
    } else if (grepl(nf_proc_rx, ln, perl = TRUE)) {
      key <- sub(nf_proc_rx, "\\1", ln, perl = TRUE)
      state$proc[[key]] <- ln
    } else if (!nzchar(trimws(ln))) {
      next
    } else if (length(state$proc) == 0 && !nzchar(state$executor)) {
      state$header <- unique(c(state$header, ln))
    } else {
      state$footer <- utils::tail(c(state$footer, ln), NF_FOOTER_MAX)
    }
  }
  state
}

#' Render a progress state as one block
#'
#' Skips empty sections so there's no blank header/footer space and no
#' inter-element whitespace inside the pre-wrapped code block.
#' @noRd
nf_progress_text <- function(state, footer_n = 40) {
  parts <- c(
    paste(state$header, collapse = "\n"),
    state$executor,
    paste(unlist(state$proc), collapse = "\n"),
    state$more %||% "",
    paste(utils::tail(state$footer, footer_n), collapse = "\n")
  )
  paste(parts[nzchar(parts)], collapse = "\n")
}
