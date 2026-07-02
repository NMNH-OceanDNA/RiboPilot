#!/usr/bin/env Rscript
# Reconcile barrnap+ITSx into the feature model. Args: <db> <id> <gff> <positions>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
reconcile_one(con, a[2], gff = a[3], positions = a[4])
cat("reconciled:", a[2], "\n")
