#!/usr/bin/env Rscript
# Print the composed -db argument for a local search, or exit 1 if unusable.
# Args: [panel_prefix] [decoy_dir]
# An empty or absent prefix means "use whatever panel is installed", which is
# how a containerised run finds the panel baked into the image.
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
local({
  panel <- if (length(a) >= 1 && !is.na(a[1]) && nzchar(trimws(a[1]))) trimws(a[1]) else ""
  if (!nzchar(panel)) {
    found <- RiboPilot:::discover_panel()
    if (is.null(found)) quit(status = 1L)
    panel <- found
  }
  # decoys ship beside the panel, so default to its own directory; the optional
  # second argument exists only to let a site keep them somewhere else
  dd <- if (length(a) >= 2 && nzchar(a[2])) a[2] else dirname(panel)
  arg <- RiboPilot:::local_blast_db_arg(panel, dd)
  if (is.null(arg)) quit(status = 1L)
  cat(arg, "\n", sep = "")
})
