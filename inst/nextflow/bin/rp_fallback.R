#!/usr/bin/env Rscript
# Assembly fallback for one sample (Option B2 baiting / Option A de-novo).
# Writes <id>_assembly_<n>.fasta into cwd (published by the process).
# Args: <id> <r1> <r2> <method> <threads> <attempt> [seeds] [labels] [bowtie2] [spades] [seeds_url] [labels_url]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_fallback.R <id> <r1> <r2> <method> <threads> <attempt> [seeds] [labels] [bowtie2] [spades] [seeds_url] [labels_url]")
db <- function(i) if (!is.na(a[i]) && nzchar(a[i])) a[i] else a[i + 4]
n <- fallback_assemble_one(a[1], a[2], a[3],
  reference = db(7), labels_db = db(8),
  out_dir = ".", method = a[4], threads = as.integer(a[5]),
  bowtie2_opts = a[9], spades_opts = a[10])
st <- attr(n, "retry_status")
if (!is.null(st) && as.integer(a[6]) <= 3L) quit(status = st, save = "no")
cat("fallback done:", a[1], "\n")
