#' Initialize a new RiboPilot project database
#'
#' Creates the SQLite tracking database. Mirrors MitoPilot's `new_db()` spine
#' (samples / preprocess / pre_opts / assemble / assemble_opts / blast_opts /
#' assemblies) but with nuclear ribosomal (nrDNA) assembly options. Annotation,
#' curation, and export tables are added in later build increments.
#'
#' @param db_path Path to the new database file.
#' @param mapping_fn Path to the mapping CSV (columns ID, Taxon, R1, R2).
#' @param mapping_id Column name to use as the primary key (default "ID").
#' @param mapping_taxon Column name with a taxonomic identifier (default "Taxon").
#' @param assemble_cpus Default number of CPUs for assembly.
#' @param assemble_memory Default memory (GB) for assembly.
#' @param seeds_db Path/URL to the GetOrganelle nr seeds database (built per
#'   clade, see the reference-DB builder). May be NULL until built.
#' @param labels_db Path/URL to the GetOrganelle nr labels database.
#' @param getOrganelle Default GetOrganelle command-line options.
#' @param target Read-recruitment target passed to GetOrganelle (default
#'   "anonym", using the custom seed/label DBs).
#' @param max_paths Max assembly paths for a sample to continue (default 10).
#' @param max_scaffolds Max scaffolds within a path to continue (default 10).
#' @param min_assembly_length Minimum contig length (bp) to retain (default 500).
#' @param ... Reserved for forward compatibility.
#'
#' @export
new_db <- function(
    db_path = "./.sqlite",
    mapping_fn = NULL,
    mapping_id = "ID",
    mapping_taxon = "Taxon",
    assemble_cpus = 6,
    assemble_memory = 24,
    seeds_db = "https://raw.githubusercontent.com/dmacguigan/RiboPilot/main/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
    labels_db = "https://raw.githubusercontent.com/dmacguigan/RiboPilot/main/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
    getOrganelle = paste(
      "-F anonym",
      "-R 10 -k 21,45,65,85,105,115",
      "--larger-auto-ws"
    ),
    target = "anonym",
    max_paths = 10,
    max_scaffolds = 10,
    min_assembly_length = 500,
    ...) {
  # Read mapping file ----
  if (is.null(mapping_fn)) {
    mapping_fn <- "./mapping.csv"
    if (!file.exists(mapping_fn)) {
      stop("Mapping file not found")
    }
  }
  mapping <- utils::read.csv(mapping_fn)
  mapping[[mapping_id]] <- as.character(mapping[[mapping_id]])

  # Validate ID column ----
  if (any(duplicated(mapping[[mapping_id]]))) {
    bad <- unique(mapping[[mapping_id]][duplicated(mapping[[mapping_id]])])
    message("problematic IDs:\n", paste(bad, collapse = ", "))
    stop("Duplicate IDs found in mapping file")
  }
  if (any(nchar(mapping[[mapping_id]]) > 18)) {
    bad <- mapping[[mapping_id]][nchar(mapping[[mapping_id]]) > 18]
    message("problematic IDs:\n", paste(bad, collapse = ", "))
    stop("IDs must be no more than 18 characters")
  }
  if (any(!grepl("^[a-zA-Z0-9_:-]+$", mapping[[mapping_id]]))) {
    bad <- mapping[[mapping_id]][!grepl("^[a-zA-Z0-9_:-]+$", mapping[[mapping_id]])]
    message("problematic IDs:\n", paste(bad, collapse = ", "))
    stop("IDs must contain only alphanumeric characters, dashes, underscores, and colons")
  }

  # Default GetOrganelle nr databases (broad Actinopterygii) when unset ----
  seeds_db <- seeds_db %||% DEFAULT_SEEDS_DB
  labels_db <- labels_db %||% DEFAULT_LABELS_DB
  # expand '~' in local paths so Nextflow/GetOrganelle can find them
  if (!grepl("^https?://", seeds_db)) seeds_db <- path.expand(seeds_db)
  if (!grepl("^https?://", labels_db)) labels_db <- path.expand(labels_db)
  # Fail here, not four hours into a run as an opaque GetOrganelle error.
  check_reference_db(seeds_db, "seeds")
  check_reference_db(labels_db, "labels")

  # Connect ----
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = db_path)
  on.exit(DBI::dbDisconnect(con))

  # samples ----
  mapping <- mapping |>
    dplyr::mutate(
      ID = .data[[mapping_id]],
      Taxon = .data[[mapping_taxon]],
      taxid = NA_integer_,
      sci_name = NA_character_,
      rank = NA_character_,
      match_type = NA_character_,
      lineage = NA_character_,
      lineage_taxids = NA_character_,
      taxonomy_switch = 1L,
      taxonomy_lock = 0L,
      export_group = NA_character_,
      reference_support = NA_integer_,
      time_stamp = NA_integer_
    )
  glue::glue_sql(
    "CREATE TABLE samples ({cols*}, PRIMARY KEY (ID))",
    cols = colnames(mapping),
    .con = con
  ) |> DBI::dbExecute(con, statement = _)
  dplyr::tbl(con, "samples") |>
    dplyr::rows_upsert(mapping, in_place = TRUE, copy = TRUE, by = "ID")

  # preprocess ----
  DBI::dbExecute(
    con,
    "CREATE TABLE preprocess (
      ID TEXT NOT NULL,
      R1 TEXT,
      R2 TEXT,
      pre_opts TEXT NOT NULL,
      reads INTEGER,
      trimmed_reads INTEGER,
      mean_length INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (ID)
    );"
  )
  dplyr::tbl(con, "preprocess") |>
    dplyr::rows_upsert(
      mapping |>
        dplyr::select(ID, R1, R2) |>
        dplyr::mutate(
          pre_opts = "default",
          reads = NA_real_,
          trimmed_reads = NA_real_,
          mean_length = NA_real_,
          time_stamp = NA_integer_
        ),
      in_place = TRUE, copy = TRUE, by = "ID"
    )

  # pre_opts ----
  DBI::dbExecute(
    con,
    "CREATE TABLE pre_opts (
      pre_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      fastp TEXT,
      PRIMARY KEY (pre_opts)
    );"
  )
  dplyr::tbl(con, "pre_opts") |>
    dplyr::rows_upsert(
      data.frame(
        pre_opts = "default",
        cpus = 4,
        memory = 16,
        fastp = "--trim_poly_g --correction --detect_adapter_for_pe --dont_eval_duplication"
      ),
      in_place = TRUE, copy = TRUE, by = "pre_opts"
    )

  # assemble ----
  DBI::dbExecute(
    con,
    "CREATE TABLE assemble (
      ID TEXT NOT NULL,
      length TEXT,
      topology TEXT,
      paths INTEGER,
      scaffolds INTEGER,
      assemble_notes TEXT,
      assemble_switch INTEGER,
      assemble_lock INTEGER,
      hide_switch INTEGER,
      assemble_opts TEXT,
      blast_opts TEXT,
      triage_opts TEXT DEFAULT 'default',
      blast_accession TEXT,
      blast_species TEXT,
      blast_pident REAL,
      blast_qcovs REAL,
      blast_evalue REAL,
      blast_lineage TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID)
    );"
  )
  dplyr::tbl(con, "assemble") |>
    dplyr::rows_upsert(
      mapping |>
        dplyr::select(ID) |>
        dplyr::mutate(
          length = NA_character_,
          topology = NA_character_,
          paths = NA_integer_,
          scaffolds = NA_integer_,
          assemble_notes = NA_character_,
          assemble_switch = 1,
          assemble_lock = 0,
          hide_switch = 0,
          assemble_opts = "default",
          blast_opts = "default",
          time_stamp = NA_integer_
        ),
      in_place = TRUE, copy = TRUE, by = "ID"
    )

  # assemble_opts (nrDNA / GetOrganelle) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE assemble_opts (
      assemble_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      getOrganelle TEXT,
      seeds_db TEXT,
      labels_db TEXT,
      target TEXT,
      max_paths INTEGER,
      max_scaffolds INTEGER,
      min_assembly_length INTEGER,
      assembler TEXT DEFAULT 'getorganelle',
      bowtie2 TEXT DEFAULT '--local',
      spades TEXT DEFAULT '--isolate',
      PRIMARY KEY (assemble_opts)
    );"
  )
  dplyr::tbl(con, "assemble_opts") |>
    dplyr::rows_upsert(
      data.frame(
        assemble_opts = "default",
        cpus = assemble_cpus,
        memory = assemble_memory,
        getOrganelle = getOrganelle,
        seeds_db = seeds_db %||% NA_character_,
        labels_db = labels_db %||% NA_character_,
        target = target,
        max_paths = max_paths,
        max_scaffolds = max_scaffolds,
        min_assembly_length = min_assembly_length,
        assembler = "getorganelle",
        bowtie2 = "--local",
        spades = "--isolate"
      ),
      in_place = TRUE, copy = TRUE, by = "assemble_opts"
    )

  # blast_opts (against nrDNA reference DB) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_opts (
      blast_opts TEXT NOT NULL,
      run_blast INTEGER,
      entrez_query TEXT,
      extra_opts TEXT,
      max_target_seqs INTEGER,
      PRIMARY KEY (blast_opts)
    );"
  )
  dplyr::tbl(con, "blast_opts") |>
    dplyr::rows_upsert(
      data.frame(
        blast_opts = "default",
        run_blast = 1L,
        # Unrestricted by default. Restricting the search to nuclear rRNA titles
        # also restricts what a contaminant is allowed to look like, so a
        # symbiont / bacterial 16S / mitochondrial assembly could only ever
        # match a distant nuclear rDNA record and appeared plausible. The
        # nuclear-rDNA filter is applied client-side in
        # parse_blast_candidates() instead. Setting this blinds the
        # contamination screen to anything outside the query.
        entrez_query = "",
        extra_opts = "",
        # Needs enough depth that real nuclear rDNA references still make the
        # panel after the client-side filter drops non-rDNA hits.
        max_target_seqs = 50L
      ),
      in_place = TRUE, copy = TRUE, by = "blast_opts"
    )

  # assemblies (per path/scaffold sequences) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE assemblies (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL,
      scaffold INTEGER NOT NULL,
      topology TEXT,
      length INTEGER,
      length_raw INTEGER,
      sequence TEXT,
      depth TEXT,
      gc TEXT,
      errors TEXT,
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
      edit_positions TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    );"
  )

  # blast_hits (full retained hit distribution per unit; see BLAST_HITS_DDL) ----
  DBI::dbExecute(con, BLAST_HITS_DDL)

  # assembly_blast (per-path BLAST hits) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE assembly_blast (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL,
      blast_opts TEXT,
      blast_accession TEXT,
      blast_species TEXT,
      blast_pident REAL,
      blast_qcovs REAL,
      blast_evalue REAL,
      blast_lineage TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path)
    );"
  )

  # blast_ref_candidates (ranked remote-GenBank hits per assembly unit) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_ref_candidates (
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
    );"
  )

  # blast_screen (contamination screen: top remote hit over ALL of core_nt) ----
  # Distinct from blast_ref_candidates, which is filtered to nuclear rDNA. This
  # row is the unfiltered best hit and is what the contaminant check reads: it
  # is the only place that can say "this is not the host, and not even nuclear
  # rDNA".
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_screen (
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
    );"
  )

  # blast_ref_sequences (fetched NCBI reference records, deduped by accession) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_ref_sequences (
      accession TEXT NOT NULL,
      species TEXT,
      title TEXT,
      taxid TEXT,
      lineage TEXT,
      sequence TEXT,
      length INTEGER,
      time_stamp INTEGER,
      PRIMARY KEY (accession)
    );"
  )

  # blast_ref_annotations (rRNA/ITS features parsed from a reference record) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_ref_annotations (
      accession TEXT NOT NULL,
      element TEXT NOT NULL,
      type TEXT,
      start INTEGER,
      stop INTEGER,
      strand TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (accession, element)
    );"
  )

  # blast_ref_alignment (per-feature sample-vs-reference nucleotide alignment) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE blast_ref_alignment (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      accession TEXT NOT NULL,
      element TEXT NOT NULL,
      pid REAL,
      aln_cov REAL,
      aln_len INTEGER,
      sample_aln TEXT,
      ref_aln TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, accession, element)
    );"
  )

  # coverage (per-position read depth per assembly unit) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE coverage (
      ID TEXT NOT NULL,
      path INTEGER NOT NULL DEFAULT 1,
      scaffold INTEGER NOT NULL DEFAULT 1,
      mean_depth REAL,
      median_depth REAL,
      min_depth INTEGER,
      max_depth INTEGER,
      depths TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    );"
  )

  # assembly_backup (original pre-trim state so a trim can be undone) ----
  DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)

  # annotations (nrDNA feature model per assembly unit) ----
  # status/flags are filled by curation (curate_one); tool/edited by annotation.
  DBI::dbExecute(
    con,
    "CREATE TABLE annotations (
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
      tool TEXT,
      status TEXT,
      flags TEXT,
      edited INTEGER,
      notes TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold, element)
    );"
  )

  # annotate (per-unit annotation + curation summary) ----
  # One row per assembly unit (ID, path, scaffold); seeded per unit by
  # load_assemblies() after assembly. n_flags / annotate_notes filled by curation;
  # curate_opts selects the per-clade ruleset (see curate_opts table).
  DBI::dbExecute(
    con,
    "CREATE TABLE annotate (
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
      reviewed TEXT,
      problematic TEXT,
      partial TEXT,
      export_group TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    );"
  )
  # Seed a placeholder unit (1,1) per sample; load_assemblies() fills the real
  # per-scaffold units once assembly has produced them.
  dplyr::tbl(con, "annotate") |>
    dplyr::rows_upsert(
      data.frame(ID = mapping$ID, path = 1L, scaffold = 1L,
        annotate_switch = 1L, annotate_lock = 0L,
        annotate_opts = "default", curate_opts = "default", reviewed = "no"),
      in_place = TRUE, copy = TRUE, by = c("ID", "path", "scaffold")
    )

  # annotate_opts (barrnap + ITSx settings) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE annotate_opts (
      annotate_opts TEXT NOT NULL,
      cpus INTEGER,
      memory INTEGER,
      rrna_tool TEXT,
      rrna_kingdom TEXT,
      itsx_taxon TEXT,
      itsx_opts TEXT,
      PRIMARY KEY (annotate_opts)
    );"
  )
  # rrna_tool: 'pybarrnap' (default; euk-capable Python reimplementation) or
  # 'barrnap'. rrna_kingdom values differ by tool (pybarrnap: bac/arc/euk/all;
  # barrnap: bac/arc/fun).
  dplyr::tbl(con, "annotate_opts") |>
    dplyr::rows_upsert(
      data.frame(annotate_opts = "default", cpus = 4L, memory = 8L,
        rrna_tool = "pybarrnap", rrna_kingdom = "euk",
        itsx_taxon = "all", itsx_opts = ""),
      in_place = TRUE, copy = TRUE, by = "annotate_opts"
    )

  # curate_opts (per-clade curation ruleset selection) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE curate_opts (
      curate_opts TEXT NOT NULL,
      ruleset TEXT,
      min_anchor_pid REAL,
      check_contaminant INTEGER,
      PRIMARY KEY (curate_opts)
    );"
  )
  dplyr::tbl(con, "curate_opts") |>
    dplyr::rows_upsert(
      data.frame(curate_opts = "default", ruleset = "metazoa_generic",
        min_anchor_pid = 85, check_contaminant = 1L),
      in_place = TRUE, copy = TRUE, by = "curate_opts"
    )

  # triage_opts (per-rule triage actions; see TRIAGE_RULE_DEFS) ----
  DBI::dbExecute(con, TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, TRIAGE_OPTS_SEED)

  # screen (contamination verdict per unit; see SCREEN_DDL) ----
  DBI::dbExecute(con, SCREEN_DDL)

  # path_scores (taxonomy-aware candidate scoring + auto-select) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE path_scores (
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
    );"
  )

  # export_opts (named FASTA-header configurations) ----
  ensure_export_opts(con)

  # stamp the current schema version so migrate_db() is a no-op on fresh DBs
  DBI::dbExecute(con, sprintf("PRAGMA user_version = %d", RIBOPILOT_SCHEMA_VERSION))

  invisible(return())
}
