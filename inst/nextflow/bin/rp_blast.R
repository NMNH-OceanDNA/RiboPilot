#!/usr/bin/env Rscript
# Parse remote-GenBank blastn hits into ranked reference candidates + the
# contamination screen.
# Args: <db> <id> <path> <scaffold> <hits_file> [taxdump_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 5, "rp_blast.R <db> <id> <path> <scaffold> <hits_file> [taxdump_dir]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  taxdump_dir <- if (length(a) >= 6 && nzchar(a[6])) a[6] else NULL
  parse_blast_candidates(con, a[2], path = as.integer(a[3]), scaffold = as.integer(a[4]),
    hits_file = a[5], taxdump_dir = taxdump_dir)
})
cat("blast candidates written:", a[2], a[3], a[4], "\n")
