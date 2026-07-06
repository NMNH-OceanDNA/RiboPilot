#!/usr/bin/env Rscript
# Assembly fallback for one sample (Option B2 baiting / Option A de-novo).
# Writes <id>_assembly_<n>.fasta into cwd (published by the process).
# Args: <db> <id> <r1> <r2> <method> <threads>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
opt <- DBI::dbGetQuery(con,
  "SELECT ao.seeds_db, ao.labels_db FROM assemble a
     JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
  params = list(a[2]))
fallback_assemble_one(a[2], a[3], a[4],
  reference = opt$seeds_db[1], labels_db = opt$labels_db[1],
  out_dir = ".", method = a[5], threads = as.integer(a[6]))
cat("fallback done:", a[2], "\n")
