#!/usr/bin/env Rscript
# Re-derive reference candidates + screen from stored blast_hits.
# Args: <db> [id ...]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 1, "rp_rederive.R <db> [id ...]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  ids <- if (length(a) > 1) a[-1] else NULL
  n <- rederive_blast_products(con, ids = ids)
  cat("units re-derived:", n, "\n")
})
