#!/usr/bin/env Rscript
# Curate one unit's reconciled feature model (structural + reference flags).
# Args: <db> <id> <path> <scaffold>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_curate.R <db> <id> <path> <scaffold>")
res <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  res <- curate_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
  res
})
cat("curated:", a[2], a[3], a[4], "-", res$n_flags, "flags\n")
