#!/usr/bin/env Rscript
# Store a fetched NCBI reference. Args: <db> <accession> <fasta_file> <ft_file>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_ref_fetch.R <db> <accession> <fasta_file> <ft_file>")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  write_blast_ref(con, a[2], a[3], a[4])
})
cat("ref stored:", a[2], "\n")
