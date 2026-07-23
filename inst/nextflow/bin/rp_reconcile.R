#!/usr/bin/env Rscript
# Reconcile barrnap+ITSx into the feature model for one unit.
# Args: <db> <id> <path> <scaffold> <gff> <positions>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
reconcile_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]),
  gff = a[5], positions = a[6])
cat("reconciled:", a[2], a[3], a[4], "\n")
