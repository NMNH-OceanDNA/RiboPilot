#!/usr/bin/env Rscript
# Annotate a multi-record reference fasta and store it.
# Args: <db> <fasta> [kingdom] [cpu]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_store.R <db> <fasta> [kingdom] [cpu]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  kingdom <- if (length(a) >= 3 && nzchar(a[3])) a[3] else "euk"
  cpu <- if (length(a) >= 4 && nzchar(a[4])) as.integer(a[4]) else 2L
  ann <- annotate_records(a[2], kingdom = kingdom, cpu = cpu)
  n <- write_blast_refs(con, a[2], ann)
  cat("references stored:", n, "\n")
})
