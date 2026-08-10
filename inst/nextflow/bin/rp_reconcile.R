#!/usr/bin/env Rscript
# Reconcile barrnap+ITSx into the feature model for one unit.
# Args: <db> <id> <path> <scaffold> <gff> <positions>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_reconcile.R <db> <id> <path> <scaffold> <gff> <positions>")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  reconcile_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]),
    gff = a[5], positions = a[6])
})
cat("reconciled:", a[2], a[3], a[4], "\n")
