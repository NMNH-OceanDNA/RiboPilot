#!/usr/bin/env Rscript
# Per-feature sample-vs-reference alignment. Args: <db> <id> <path> <scaffold>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_ref_align.R <db> <id> <path> <scaffold>")
n <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  compute_ref_alignment(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
})
cat("ref alignments:", a[2], a[3], a[4], n, "\n")
