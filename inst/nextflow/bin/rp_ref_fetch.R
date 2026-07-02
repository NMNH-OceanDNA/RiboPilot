#!/usr/bin/env Rscript
# Store a fetched NCBI reference. Args: <db> <accession> <fasta_file> <ft_file>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
write_blast_ref(con, a[2], a[3], a[4])
cat("ref stored:", a[2], "\n")
