#!/usr/bin/env Rscript
# Args: <db> <id> <path> <scaffold> [dbarg] [taxdump_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_feature_blast.R <db> <id> <path> <scaffold> [dbarg] [taxdump_dir]")
opt <- function(i) if (length(a) >= i && nzchar(a[i])) a[i] else NULL
r <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  r <- feature_blast_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), dbarg = opt(5),
    taxdump_dir = opt(6), threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "feature_blast_work"))
  if (!is.null(r)) RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
  r
})
cat("feature_blast:", a[2], a[3], a[4], "-", NROW(r), "features\n")
