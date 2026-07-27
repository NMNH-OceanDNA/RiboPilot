#!/usr/bin/env Rscript
# Split a combined per-scaffold `samtools depth -a` into per-unit coverage rows.
# Args: <db> <id> <depth_file>  (depth ref column = path.scaffold)
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- write_coverage_units(con, a[2], a[3])
cat("coverage:", a[2], "-", n, "unit(s)\n")
