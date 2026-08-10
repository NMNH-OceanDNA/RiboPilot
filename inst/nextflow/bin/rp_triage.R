#!/usr/bin/env Rscript
# Propose keep/drop/undecided for every unit of one sample, once its evidence
# (annotation, GenBank screen, coverage) is in place.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_triage.R <db> <id>")
res <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  # 0 updated means two different things; say which one this run was. Resolved
  # before the write so a rename of an internal cannot fail the process with the
  # database already updated.
  off <- RiboPilot:::triage_all_off(
    RiboPilot:::normalize_triage_actions(RiboPilot:::triage_opts_for(con, a[2])))
  list(off = off, n = auto_triage_sample(con, a[2]))
})
cat("auto-triage: ", a[2], " - ", res$n, " unit(s) updated",
  if (res$off) " (every rule is off)" else "", "\n", sep = "")
