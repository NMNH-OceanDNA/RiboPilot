#!/usr/bin/env Rscript
# Propose keep/drop/undecided for every unit of one sample, once its evidence
# (annotation, GenBank screen, coverage) is in place.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
# 0 updated means two different things; say which one this run was. Resolved
# before the write so a rename of an internal cannot fail the process with the
# database already updated.
off <- RiboPilot:::triage_all_off(
  RiboPilot:::normalize_triage_actions(RiboPilot:::triage_opts_for(con, a[2])))
n <- auto_triage_sample(con, a[2])
cat("auto-triage: ", a[2], " - ", n, " unit(s) updated",
  if (off) " (every rule is off)" else "", "\n", sep = "")
