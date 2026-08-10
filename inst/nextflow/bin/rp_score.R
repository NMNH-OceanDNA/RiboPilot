#!/usr/bin/env Rscript
# Load a sample's assemblies + score candidate paths. Args: <db> <id> <out_dir> [taxdump_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 3, "rp_score.R <db> <id> <out_dir> [taxdump_dir]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  load_assemblies(con, a[3], a[2])
  score_one(con, a[2], a[3], taxdump_dir = if (length(a) >= 4 && nzchar(a[4])) a[4] else NULL)
})
cat("scored:", a[2], "\n")
