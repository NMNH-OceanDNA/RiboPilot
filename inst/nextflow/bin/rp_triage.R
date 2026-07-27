#!/usr/bin/env Rscript
# Propose keep/drop/undecided for every unit of one sample, once its evidence
# (annotation, GenBank screen, coverage) is in place.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- auto_triage_sample(con, a[2])
cat("auto-triage:", a[2], "-", n, "unit(s) updated\n")
