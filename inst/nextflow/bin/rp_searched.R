#!/usr/bin/env Rscript
# Record that a unit was searched. Args: <db> <id> <path> <scaffold> <md5> <n_hits> [source]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_searched.R <db> <id> <path> <scaffold> <md5> <n_hits> [source]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  src <- if (length(a) >= 7 && nzchar(a[7])) a[7] else "remote"
  mark_unit_searched(con, a[2], as.integer(a[3]), as.integer(a[4]), a[5], as.integer(a[6]),
    source = src)
})
cat("searched:", a[2], a[3], a[4], "\n")
