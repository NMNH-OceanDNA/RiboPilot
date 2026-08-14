#!/usr/bin/env Rscript
# Copy stored references into the project before fetching.
# Args: <db> <ref_dir> [id ...]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_prime.R <db> <ref_dir> [id ...]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  ids <- if (length(a) > 2) a[-(1:2)] else NULL
  n <- prime_from_ref_store(con, a[2], ids = ids)
  cat("primed from store:", n, "\n")
})
