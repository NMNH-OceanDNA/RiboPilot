#!/usr/bin/env Rscript
# Curate one sample's reconciled feature model (structural + reference flags).
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
res <- curate_one(con, a[2])
cat("curated:", a[2], "-", res$n_flags, "flags\n")
