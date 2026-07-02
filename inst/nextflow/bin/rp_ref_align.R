#!/usr/bin/env Rscript
# Per-feature sample-vs-reference alignment. Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- compute_ref_alignment(con, a[2])
cat("ref alignments:", a[2], n, "\n")
