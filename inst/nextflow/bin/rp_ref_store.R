#!/usr/bin/env Rscript
# Annotate a multi-record reference fasta and store it.
# Args: <db> <fasta> [cpu] [ref_dir]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_ref_store.R <db> <fasta> [cpu] [ref_dir]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  cpu <- if (length(a) >= 3 && nzchar(a[3])) as.integer(a[3]) else 2L
  ref_dir <- if (length(a) >= 4 && nzchar(a[4])) a[4] else ""
  RiboPilot:::normalize_pdb_deflines(a[2])
  # reference annotation is pinned to the metazoan constants, not a per-run setting
  ann <- annotate_records(a[2], kingdom = RiboPilot:::REF_ANNOT_KINGDOM, cpu = cpu,
    taxon = RiboPilot:::REF_ANNOT_TAXON)
  n <- write_blast_refs(con, a[2], ann)
  stop_if_nothing_stored(a[2], n)
  cat("references stored:", n, "\n")
  if (nzchar(ref_dir)) {
    # the project-side store above has already succeeded; the shared store is
    # an optimisation and must never fail an otherwise successful task
    pushed <- tryCatch(
      push_to_ref_store(ref_dir, a[2], ann, "ribopilot-annotate"),
      error = function(e) {
        warning("push to reference store failed: ", conditionMessage(e), call. = FALSE)
        0L
      })
    cat("references pushed to store:", pushed, "\n")
  }
})
