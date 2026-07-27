#!/usr/bin/env Rscript
# Per-feature sample-vs-reference alignment. Args: <db> <id> <path> <scaffold>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- compute_ref_alignment(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
cat("ref alignments:", a[2], a[3], a[4], n, "\n")
