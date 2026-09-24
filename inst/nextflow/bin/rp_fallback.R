#!/usr/bin/env Rscript
# Assembly fallback for one sample (Option B2 baiting / Option A de-novo).
# Writes <id>_assembly_<n>.fasta into cwd (published by the process).
# Args: <db> <id> <r1> <r2> <method> <threads> <attempt>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 7, "rp_fallback.R <db> <id> <r1> <r2> <method> <threads> <attempt>")
opt <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  DBI::dbGetQuery(con,
    "SELECT ao.seeds_db, ao.labels_db,
            COALESCE(ao.bowtie2, '--local') AS bowtie2,
            COALESCE(ao.spades, '--isolate') AS spades
       FROM assemble a
       JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts WHERE a.ID = ?",
    params = list(a[2]))
})
n <- fallback_assemble_one(a[2], a[3], a[4],
  reference = opt$seeds_db[1], labels_db = opt$labels_db[1],
  out_dir = ".", method = a[5], threads = as.integer(a[6]),
  bowtie2_opts = opt$bowtie2[1], spades_opts = opt$spades[1])
st <- attr(n, "retry_status")
if (!is.null(st) && as.integer(a[7]) <= 3L) quit(status = st, save = "no")
cat("fallback done:", a[2], "\n")
