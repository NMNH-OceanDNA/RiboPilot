#!/usr/bin/env Rscript
# Print the composed -db argument for a local search, or exit 1 if unusable.
# Args: <panel_prefix> [decoy_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 1, "rp_dbarg.R <panel_prefix> [decoy_dir]")
local({
  # decoys ship beside the panel, so default to its own directory; the optional
  # second argument exists only to let a site keep them somewhere else
  dd <- if (length(a) >= 2 && nzchar(a[2])) a[2] else dirname(a[1])
  arg <- suppressWarnings(RiboPilot:::local_blast_db_arg(a[1], dd))
  if (is.null(arg)) quit(status = 1L)
  cat(arg, "\n", sep = "")
})
