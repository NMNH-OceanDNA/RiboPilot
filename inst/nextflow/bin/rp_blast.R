#!/usr/bin/env Rscript
# Parse remote-GenBank blastn hits into ranked reference candidates.
# Args: <db> <id> <hits_file>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
parse_blast_candidates(con, a[2], a[3])
cat("blast candidates written:", a[2], "\n")
