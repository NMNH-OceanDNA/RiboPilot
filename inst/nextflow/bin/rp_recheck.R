#!/usr/bin/env Rscript
# Args: <db> <id> <path> <scaffold> <rev> <mode> <tools_dir> [dbarg] [taxdump_dir] [label_db]
# Exits 0 when the recheck ran, whether the copy ended done, stale, or failed:
# the database holds the outcome.
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 7, "rp_recheck.R <db> <id> <path> <scaffold> <rev> <mode> <tools_dir> [dbarg] [taxdump_dir] [label_db]")
opt <- function(i) if (length(a) >= i && nzchar(a[i])) a[i] else NULL
r <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  RiboPilot:::recheck_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), rev = a[5], mode = a[6],
    tools_dir = a[7], dbarg = opt(8), taxdump_dir = opt(9), label_db = opt(10),
    threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "recheck_work"))
})
cat("recheck:", a[2], a[3], a[4], "-", r, "\n")
