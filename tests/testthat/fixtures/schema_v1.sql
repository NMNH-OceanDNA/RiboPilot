CREATE TABLE samples ('ID', 'Taxon', 'R1', 'R2', 'Assembly', 'taxid', 'sci_name', 'rank', 'match_type', 'tax_alternatives', 'lineage', 'lineage_taxids', 'ref_id', 'user_group', 'ref_source', 'taxonomy_switch', 'taxonomy_lock', 'time_stamp', PRIMARY KEY (ID))
;
CREATE TABLE preprocess (
      ID TEXT NOT NULL,
      R1 TEXT,
      R2 TEXT,
      pre_opts TEXT NOT NULL,
      reads INTEGER,
      trimmed_reads INTEGER,
      mean_length INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID)
    )
;
CREATE TABLE pre_opts (
      pre_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      fastp TEXT,
      PRIMARY KEY (pre_opts)
    )
;
CREATE TABLE assemble (
      ID TEXT NOT NULL,
      length TEXT,
      topology TEXT,
      paths INTEGER,
      scaffolds INTEGER,
      assemble_notes TEXT,
      find_rrna_notes TEXT,
      user_notes TEXT,
      assemble_switch INTEGER,
      assemble_lock INTEGER,
      assemble_opts TEXT,
      blast_opts TEXT,
      triage_opts TEXT DEFAULT 'default',
      time_stamp INTEGER,
      PRIMARY KEY (ID)
    )
;
CREATE TABLE assemble_opts (
      assemble_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      getOrganelle TEXT,
      seeds_db TEXT,
      labels_db TEXT,
      max_paths INTEGER,
      max_scaffolds INTEGER,
      min_assembly_length INTEGER,
      assembler TEXT DEFAULT 'getorganelle',
      find_rrna INTEGER DEFAULT 1,
      bowtie2 TEXT DEFAULT '--local',
      spades TEXT DEFAULT '--isolate',
      PRIMARY KEY (assemble_opts)
    )
;
CREATE TABLE blast_opts (
      blast_opts TEXT NOT NULL,
      run_blast INTEGER,
      entrez_query TEXT,
      remote_blast INTEGER,
      remote_fallback INTEGER,
      extra_opts TEXT,
      max_target_seqs INTEGER,
      PRIMARY KEY (blast_opts)
    )
;
CREATE TABLE assemblies (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL,
      scaffold INTEGER NOT NULL,
      topology TEXT,
      length INTEGER,
      sequence TEXT,
      -- the contig name this unit had in the file it came from: the original
      -- GetOrganelle/SPAdes header, or the user's own contig name when the
      -- assembly was supplied. Only the latter is shown in the app.
      source_contig TEXT,
      ignore INTEGER,
      edited INTEGER,
      blast_accession TEXT,
      blast_species TEXT,
      blast_pident REAL,
      blast_qcovs REAL,
      blast_evalue REAL,
      blast_lineage TEXT,
      blast_source TEXT,
      triage TEXT,
      triage_reason TEXT,
      triage_source TEXT,
      triage_rule TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    )
;
CREATE TABLE blast_hits (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  hit_rank INTEGER NOT NULL,
  hsp INTEGER NOT NULL DEFAULT 1,
  accession TEXT,
  title TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  pident REAL,
  qcovs REAL,
  evalue REAL,
  bitscore REAL,
  aln_len INTEGER,
  qstart INTEGER,
  qend INTEGER,
  sstart INTEGER,
  send INTEGER,
  slen INTEGER,
  is_nuclear_rrna INTEGER,
  hit_class TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, hit_rank, hsp)
)
;
CREATE TABLE blast_searched (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  query_md5 TEXT NOT NULL,
  n_hits INTEGER,
  source TEXT,
  db TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
)
;
CREATE TABLE blast_ref_candidates (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      rank INTEGER NOT NULL,
      accession TEXT,
      species TEXT,
      title TEXT,
      pident REAL,
      qcovs REAL,
      evalue REAL,
      bitscore REAL,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, rank)
    )
;
CREATE TABLE blast_screen (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      accession TEXT,
      title TEXT,
      species TEXT,
      staxid INTEGER,
      lineage TEXT,
      lineage_taxids TEXT,
      rank_distance INTEGER,
      pident REAL,
      qcovs REAL,
      evalue REAL,
      bitscore REAL,
      is_nuclear_rrna INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    )
;
CREATE TABLE blast_ref_sequences (
      accession TEXT NOT NULL,
      species TEXT,
      title TEXT,
      taxid TEXT,
      lineage TEXT,
      sequence TEXT,
      length INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (accession)
    )
;
CREATE TABLE blast_ref_annotations (
      accession TEXT NOT NULL,
      element TEXT NOT NULL,
      type TEXT,
      start INTEGER,
      stop INTEGER,
      strand TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (accession, element)
    )
;
CREATE TABLE blast_ref_alignment (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      accession TEXT NOT NULL,
      element TEXT NOT NULL,
      pid REAL,
      aln_cov REAL,
      aln_len INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, accession, element)
    )
;
CREATE TABLE coverage (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      mean_depth REAL,
      median_depth REAL,
      min_depth INTEGER,
      max_depth INTEGER,
      depths TEXT,
      clips TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    )
;
CREATE TABLE assembly_backup (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  sequence TEXT,
  length INTEGER,
  annotations TEXT,
  coverage TEXT,
  chimera TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
)
;
CREATE TABLE annotations (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      element TEXT NOT NULL,
      type TEXT,
      start INTEGER,
      stop INTEGER,
      length INTEGER,
      strand TEXT,
      model_coverage REAL,
      partial_start INTEGER,
      partial_stop INTEGER,
      tool TEXT,
      status TEXT,
      flags TEXT,
      edited INTEGER,
      notes TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, element)
    )
;
CREATE TABLE annotate (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      annotate_switch INTEGER,
      annotate_lock INTEGER,
      n_features INTEGER,
      missing TEXT,
      chimeric INTEGER,
      annotate_opts TEXT,
      curate_opts TEXT,
      n_flags INTEGER,
      annotate_notes TEXT,
      curate_flags TEXT,
      itsx_warnings TEXT,
      reviewed TEXT,
      problematic TEXT,
      partial TEXT DEFAULT 'no',
      export_group TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    )
;
CREATE TABLE annotate_opts (
      annotate_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      rrna_tool TEXT,
      rrna_kingdom TEXT,
      itsx_taxon TEXT,
      itsx_opts TEXT,
      PRIMARY KEY (annotate_opts)
    )
;
CREATE TABLE curate_opts (
      curate_opts TEXT NOT NULL,
      ruleset TEXT,
      PRIMARY KEY (curate_opts)
    )
;
CREATE TABLE triage_opts (
  triage_opts            TEXT NOT NULL,
  rule_not_annotated     TEXT,
  rule_no_screen         TEXT,
  rule_no_feat_non_rdna  TEXT,
  rule_no_feat_rdna_hit  TEXT,
  rule_feat_non_rdna_hit TEXT,
  rule_odd_one_out       TEXT,
  rule_screen_off_target TEXT,
  rule_chimeric          TEXT,
  PRIMARY KEY (triage_opts)
)
;
CREATE TABLE rrna_candidates (
  ID TEXT NOT NULL,
  contig TEXT NOT NULL,
  length INTEGER,
  rrna_features TEXT,
  rrna_bases INTEGER,
  rank INTEGER,
  selected INTEGER,
  reason TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, contig)
)
;
CREATE TABLE screen (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  screen_verdict    TEXT,
  screen_class      TEXT,
  screen_confidence REAL,
  screen_reason     TEXT,
  screen_version    INTEGER,
  lca_taxid         INTEGER,
  lca_taxon         TEXT,
  lca_rank          TEXT,
  lca_n_hits        INTEGER,
  host_concordant   INTEGER,
  margin            REAL,
  variant_class     TEXT,
  variant_group     TEXT,
  variant_identity  REAL,
  time_stamp        INTEGER,
  PRIMARY KEY (ID, path, scaffold)
)
;
CREATE TABLE chimera (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  score REAL,
  label TEXT,
  junctions TEXT,
  reason TEXT,
  segments TEXT,
  gene_taxa TEXT,
  channels TEXT,
  ribotyper_codes TEXT,
  db_source TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
)
;
CREATE TABLE chimera_windows (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  start INTEGER NOT NULL,
  stop INTEGER,
  accession TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  hit_class TEXT,
  pident REAL,
  bitscore REAL,
  depth_median REAL,
  clip_frac REAL,
  db_source TEXT,
  PRIMARY KEY (ID, path, scaffold, start)
)
;
CREATE TABLE path_scores (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL,
      length INTEGER,
      depth REAL,
      anchors TEXT,
      blast_pident REAL,
      hit_id TEXT,
      hit_taxid TEXT,
      hit_species TEXT,
      blast_lineage TEXT,
      s_tax REAL,
      s_identity REAL,
      s_anchors REAL,
      s_depth REAL,
      s_length REAL,
      score REAL,
      rank INTEGER,
      auto_select INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path)
    )
;
CREATE TABLE reference (
  ref_id           TEXT NOT NULL,
  clade_taxid      INTEGER,
  clade_name       TEXT,
  clade_rank       TEXT,
  query            TEXT,
  n_samples        INTEGER,
  n_records        INTEGER,
  reference_switch INTEGER,
  source           TEXT,
  seeds_db         TEXT,
  labels_db        TEXT,
  log_file         TEXT,
  pid              INTEGER,
  preset_name      TEXT,
  fail_reason      TEXT,
  time_stamp       INTEGER,
  PRIMARY KEY (ref_id)
)
;
CREATE TABLE reference_config (
  rowid_key   INTEGER NOT NULL,
  ref_dir     TEXT,
  min_records INTEGER,
  max_records INTEGER,
  blastdb     TEXT,
  blastdb_version TEXT,
  blastdb_info TEXT,
  PRIMARY KEY (rowid_key)
)
;
CREATE TABLE export_opts (
  export_opts TEXT NOT NULL,
  fasta_header TEXT,
  fasta_header_marker TEXT,
  PRIMARY KEY (export_opts)
)
;
CREATE TABLE schema_history (
  event TEXT NOT NULL,
  schema_version INTEGER NOT NULL,
  ribopilot_version TEXT,
  time_stamp INTEGER,
  note TEXT
)
;
