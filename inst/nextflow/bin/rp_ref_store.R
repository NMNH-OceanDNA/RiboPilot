#!/usr/bin/env Rscript
# Annotate a multi-record reference fasta and store it.
# Args: <db> <fasta> [kingdom] [cpu] [taxon] [ref_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_store.R <db> <fasta> [kingdom] [cpu] [taxon] [ref_dir]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  kingdom <- if (length(a) >= 3 && nzchar(a[3])) a[3] else "euk"
  cpu <- if (length(a) >= 4 && nzchar(a[4])) as.integer(a[4]) else 2L
  taxon <- if (length(a) >= 5 && nzchar(a[5])) a[5] else "all"
  ref_dir <- if (length(a) >= 6 && nzchar(a[6])) a[6] else ""
  ann <- annotate_records(a[2], kingdom = kingdom, cpu = cpu, taxon = taxon)
  n <- write_blast_refs(con, a[2], ann)
  stop_if_nothing_stored(a[2], n)
  cat("references stored:", n, "\n")
  if (nzchar(ref_dir)) {
    # the project-side store above has already succeeded; the shared store is
    # an optimisation and must never fail an otherwise successful task
    pushed <- tryCatch(
      push_to_ref_store(ref_dir, a[2], ann, kingdom, taxon, "ribopilot-annotate"),
      error = function(e) {
        warning("push to reference store failed: ", conditionMessage(e), call. = FALSE)
        0L
      })
    cat("references pushed to store:", pushed, "\n")
  }
})
