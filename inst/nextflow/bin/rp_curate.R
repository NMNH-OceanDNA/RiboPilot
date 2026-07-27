#!/usr/bin/env Rscript
# Curate one unit's reconciled feature model (structural + reference flags).
# Args: <db> <id> <path> <scaffold>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
res <- curate_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]))
cat("curated:", a[2], a[3], a[4], "-", res$n_flags, "flags\n")
