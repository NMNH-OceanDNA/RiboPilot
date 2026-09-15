# Tables and columns the app reads that the hand-built fixtures below do not
# declare themselves. Kept in one place so both app fixtures agree with new_db().

fixture_fill_schema <- function(con) {
  for (ddl in list(RiboPilot:::BLAST_HITS_DDL, RiboPilot:::TRIAGE_OPTS_DDL,
                   RiboPilot:::SCREEN_DDL, RiboPilot:::ASSEMBLY_BACKUP_DDL,
                   RiboPilot:::REFERENCE_DDL, RiboPilot:::REFERENCE_CONFIG_DDL,
                   RiboPilot:::CHIMERA_DDL, RiboPilot:::CHIMERA_WINDOWS_DDL)) {
    DBI::dbExecute(con, ddl)
  }
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con,
    "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
     VALUES (1, NULL, 1000, 25000)")
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS curate_opts (curate_opts TEXT PRIMARY KEY,
    ruleset TEXT, min_anchor_pid REAL, check_contaminant INTEGER)")
  DBI::dbExecute(con, "INSERT INTO curate_opts VALUES ('default','metazoa_generic',85,1)")
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS blast_opts (blast_opts TEXT PRIMARY KEY,
    run_blast INTEGER, entrez_query TEXT, max_target_seqs INTEGER, extra_opts TEXT)")
  DBI::dbExecute(con, "INSERT OR IGNORE INTO blast_opts VALUES ('default',1,NULL,50,'')")
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS annotate_opts (annotate_opts TEXT PRIMARY KEY,
    cpus INTEGER, memory INTEGER, rrna_tool TEXT, rrna_kingdom TEXT, itsx_taxon TEXT, itsx_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO annotate_opts VALUES ('default',4,8,'pybarrnap','euk','all','')")
  if (!("assemble_opts" %in% DBI::dbListTables(con))) {
    DBI::dbExecute(con, "CREATE TABLE assemble_opts (assemble_opts TEXT PRIMARY KEY, cpus INTEGER,
      memory INTEGER, getOrganelle TEXT, seeds_db TEXT, labels_db TEXT, target TEXT, max_paths INTEGER,
      max_scaffolds INTEGER, min_assembly_length INTEGER, assembler TEXT, bowtie2 TEXT, spades TEXT)")
    DBI::dbExecute(con, "INSERT INTO assemble_opts VALUES
      ('default',4,16,'','s.fa','l.fa',NULL,10,10,1000,'getorganelle','--local','--isolate')")
  }
  RiboPilot:::ensure_export_opts(con)

  add_col <- function(tbl, col, decl) {
    if (tbl %in% DBI::dbListTables(con) && !(col %in% DBI::dbListFields(con, tbl))) {
      DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", tbl, col, decl))
    }
  }
  add_col("assemblies", "triage", "TEXT")
  add_col("assemblies", "triage_reason", "TEXT")
  add_col("assemblies", "triage_source", "TEXT")
  add_col("assemblies", "triage_rule", "TEXT")
  add_col("assemblies", "source_contig", "TEXT")
  add_col("assemble", "triage_opts", "TEXT DEFAULT 'default'")
  add_col("annotate", "chimeric", "INTEGER")
  add_col("annotations", "model_coverage", "REAL")
  add_col("coverage", "median_depth", "REAL")
  add_col("samples", "reference_support", "INTEGER")
  add_col("samples", "ref_id", "TEXT")
  invisible()
}
