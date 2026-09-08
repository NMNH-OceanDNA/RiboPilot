#!/usr/bin/env Rscript
# Split a combined per-scaffold `samtools depth -a` into per-unit coverage rows.
# Args: <db> <id> <depth_file> [clips_file]  (depth ref column = path.scaffold)
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 3, "rp_coverage.R <db> <id> <depth_file> [clips_file]")
n <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  write_coverage_units(con, a[2], a[3], if (length(a) >= 4) a[4] else NULL)
})
cat("coverage:", a[2], "-", n, "unit(s)\n")
