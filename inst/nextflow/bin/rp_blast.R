#!/usr/bin/env Rscript
# Parse remote-GenBank blastn hits into ranked reference candidates + the
# contamination screen.
# Args: <db> <id> <path> <scaffold> <hits_file> [taxdump_dir]
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
taxdump_dir <- if (length(a) >= 6 && nzchar(a[6])) a[6] else NULL
parse_blast_candidates(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]),
  hits_file = a[5], taxdump_dir = taxdump_dir)
cat("blast candidates written:", a[2], a[3], a[4], "\n")
