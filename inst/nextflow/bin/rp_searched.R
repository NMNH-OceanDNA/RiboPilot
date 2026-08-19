#!/usr/bin/env Rscript
# Record that a unit was searched. Args: <db> <id> <path> <scaffold> <md5> <n_hits> [source] [panel]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 6, "rp_searched.R <db> <id> <path> <scaffold> <md5> <n_hits> [source] [panel]")
local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  src <- if (length(a) >= 7 && nzchar(a[7])) a[7] else "remote"
  mark_unit_searched(con, a[2], as.integer(a[3]), as.integer(a[4]), a[5], as.integer(a[6]),
    source = src)
  # the blast_searched row now exists, so re-derive from it: on a unit's first
  # local search, parse_blast_candidates ran before this row existed and
  # defaulted the label to core_nt.
  derive_blast_products(con, a[2], as.integer(a[3]), as.integer(a[4]))
  # panel is the prefix DBARG actually resolved to for this unit's local
  # search (blast.nf passes "" when the unit did not search locally); a panel
  # with no manifest beside it records nothing, defensively, not an error.
  if (length(a) >= 8 && nzchar(a[8])) {
    tryCatch(record_blastdb_provenance(con, a[8]), error = function(e) NULL)
  }
})
cat("searched:", a[2], a[3], a[4], "\n")
