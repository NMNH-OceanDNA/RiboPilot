#!/usr/bin/env Rscript
# Score the contamination screen for every unit of one sample, from the BLAST
# hit distribution already stored by blast_genbank.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- screen_sample(con, a[2])
cat("screen: ", a[2], " - ", n, " unit(s)\n", sep = "")
