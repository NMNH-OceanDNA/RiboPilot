# Shared setup for chimera-detection tests: builds a minimal project DB via
# new_db() (which needs a mapping CSV and seed/label FASTA, not just a bare
# sqlite file) and returns an open connection to it.
chimera_test_con <- function(id = "S1") {
  d <- withr::local_tempdir(.local_envir = parent.frame())
  map <- file.path(d, "mapping.csv")
  utils::write.csv(data.frame(ID = id, Taxon = "Danio rerio",
    R1 = paste0(id, "_R1.fq.gz"), R2 = paste0(id, "_R2.fq.gz")), map, row.names = FALSE)
  db <- file.path(d, ".sqlite")
  fa <- file.path(d, "x.fasta")
  writeLines(c(">18S", "ACGT"), fa)
  RiboPilot::new_db(db_path = db, mapping_fn = map, check_reads = FALSE,
    seeds_db = fa, labels_db = fa)
  RiboPilot:::open_project_db(db)
}
