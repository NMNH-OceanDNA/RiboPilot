#!/usr/bin/env Rscript
# Stitch a sample's selected path. Args: <db> <id> <fasta_out>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
stitch_one(con, a[2], fasta_out = a[3])
cat("stitched:", a[2], "\n")
