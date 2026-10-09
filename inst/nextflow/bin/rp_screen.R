#!/usr/bin/env Rscript
# Score the contamination screen for every unit of one sample, from the BLAST
# hit distribution already stored by blast_genbank.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 2, "rp_screen.R <db> <id>")
n <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  n <- screen_sample(con, a[2], variants = FALSE)
  w <- attr(n, "written")
  for (i in seq_len(nrow(w))) RiboPilot:::review_follow_original(con, a[2], w$path[i], w$scaffold[i])
  n
})
cat("screen: ", a[2], " - ", n, " unit(s)\n", sep = "")
