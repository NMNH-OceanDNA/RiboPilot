#!/usr/bin/env Rscript
# Parse remote-GenBank blastn hits into ranked reference candidates + the
# contamination screen.
# Args: <db> <id> <hits_file> [taxdump_dir]
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
taxdump_dir <- if (length(a) >= 4 && nzchar(a[4])) a[4] else NULL
parse_blast_candidates(con, a[2], a[3], taxdump_dir = taxdump_dir)
cat("blast candidates written:", a[2], "\n")
