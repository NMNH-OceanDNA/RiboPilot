#!/usr/bin/env Rscript
# Args: <db> <id> <path> <scaffold> [dbarg] [taxdump_dir] [label_db]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_chimera.R <db> <id> <path> <scaffold> [dbarg] [taxdump_dir] [label_db]")
opt <- function(i) if (length(a) >= i && nzchar(a[i])) a[i] else NULL
r <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  r <- chimera_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), dbarg = opt(5),
    taxdump_dir = opt(6), label_db = opt(7), threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "chimera_work"))
  if (!is.null(r)) RiboPilot:::review_follow_original(con, a[2], as.integer(a[3]), as.integer(a[4]))
  r
})
cat("chimera:", a[2], a[3], a[4], "-", if (is.null(r)) "skipped" else r$label, "\n")
