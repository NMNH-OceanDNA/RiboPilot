#!/usr/bin/env Rscript
# Parse samtools depth output into the coverage table. Args: <db> <id> <depth_file>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
write_coverage(con, a[2], a[3])
cat("coverage written:", a[2], "\n")
