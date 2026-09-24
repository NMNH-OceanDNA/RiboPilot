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
  add_col("annotate", "curate_flags", "TEXT")
  add_col("annotations", "model_coverage", "REAL")
  add_col("coverage", "median_depth", "REAL")
  add_col("samples", "reference_support", "INTEGER")
  add_col("samples", "ref_id", "TEXT")
  add_col("samples", "tax_alternatives", "TEXT")
  invisible()
}

mk_app_project <- function(env = parent.frame()) {
  dir <- withr::local_tempdir(.local_envir = env)
  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(dir, ".sqlite"))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, Taxon TEXT, sci_name TEXT, taxid INT, rank TEXT,
    match_type TEXT, lineage TEXT, lineage_taxids TEXT, taxonomy_switch INT, taxonomy_lock INT,
    export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, assemble_lock INT, assemble_switch INT,
    assemble_opts TEXT, blast_opts TEXT, paths INT, scaffolds INT, length TEXT, assemble_notes TEXT, find_rrna_notes TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE preprocess(ID TEXT, pre_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, topology TEXT, length INT,
    sequence TEXT, ignore INT, blast_species TEXT, blast_accession TEXT, blast_pident REAL,
    blast_qcovs REAL, blast_source TEXT, source_contig TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, problematic TEXT,
    partial TEXT, n_features INT, missing TEXT, n_flags INT, annotate_notes TEXT, itsx_warnings TEXT, export_group TEXT, time_stamp INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT, start INT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL,
    min_depth INT, max_depth INT, depths TEXT)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, accession TEXT, title TEXT,
    species TEXT, staxid INT, lineage TEXT, lineage_taxids TEXT, rank_distance INT, pident REAL,
    qcovs REAL, evalue REAL, bitscore REAL, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates(ID TEXT, path INT, scaffold INT, rank INT,
    accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL, bitscore REAL, evalue REAL)")
  DBI::dbExecute(con, "CREATE TABLE path_scores(ID TEXT, path INT, length INT, depth REAL, anchors TEXT,
    blast_pident REAL, hit_id TEXT, hit_species TEXT, score REAL, rank INT, auto_select INT)")
  DBI::dbExecute(con, "INSERT INTO path_scores (ID,path,score,rank,auto_select) VALUES('S1',1,0.912,1,1)")
  DBI::dbExecute(con, "INSERT INTO path_scores (ID,path,score,rank,auto_select) VALUES('S1',2,0.431,2,0)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio','Danio rerio',7955,'species','exact',NULL,NULL,1,1,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',1,2,'default','default',2,4,'100',NULL,NULL,1)")
  DBI::dbExecute(con, "INSERT INTO preprocess VALUES('S1','default')")
  # four units; the SHORTEST is p2.s2, so a length-sorted view reorders the table
  for (u in list(c(1,1,6120), c(1,2,1442), c(2,1,5980), c(2,2,890))) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',?,?,'linear',?,'ACGT',0,'Danio rerio','A1',99,88,'core_nt',NULL)",
      params = list(u[1], u[2], u[3]))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',?,?,1,0,'default','default','no','no','no',3,'',0,NULL,NULL,NULL,1)",
      params = list(u[1], u[2]))
    DBI::dbExecute(con, "INSERT INTO coverage (ID,path,scaffold,mean_depth) VALUES('S1',?,?,100)",
      params = list(u[1], u[2]))
  }
  fixture_fill_schema(con)
  DBI::dbDisconnect(con)
  dir
}
