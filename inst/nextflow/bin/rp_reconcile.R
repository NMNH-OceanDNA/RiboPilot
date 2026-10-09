#!/usr/bin/env Rscript
# Reconcile barrnap+ITSx into the feature model for one unit.
# Args: <db> <id> <path> <scaffold> <gff> <positions> [itsx_check_prefix]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_reconcile.R <db> <id> <path> <scaffold> <gff> <positions> [itsx_check_prefix]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  reconcile_one(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]),
    gff = a[5], positions = a[6])
  # stock-ITSx warnings: advisory, must never sink the unit's annotation
  if (length(a) >= 7) tryCatch(
    itsx_check_one(con, a[2], as.integer(a[3]), as.integer(a[4]), a[7]),
    error = function(e) {
      message("ITSx check failed: ", conditionMessage(e))
      tryCatch(RiboPilot:::itsx_store_not_run(con, a[2], as.integer(a[3]), as.integer(a[4]),
        paste("ITSx check failed:", conditionMessage(e))),
        error = function(e2) message("ITSx check result not stored: ", conditionMessage(e2)))
    })
  RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
})
cat("reconciled:", a[2], a[3], a[4], "\n")
