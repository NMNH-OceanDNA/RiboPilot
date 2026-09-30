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
#' @param data_path Directory holding the raw read files named in the mapping
#'   file. Read files are resolved as `<data_path>/<R1>`, the same way the
#'   Nextflow pipeline does it. When NULL, read files cannot be located and
#'   only the paired-end check is applied.
#' @param assembly_path Directory holding the assembly files named in the
#'   mapping file's optional `Assembly` column. A sample with an assembly skips
#'   read preprocessing and assembly: its nuclear rDNA contigs are found inside
#'   the supplied file instead. Its `R1`/`R2` may be blank.
#' @param check_reads (logical) Check that the read files named in the mapping
#'   file exist (default TRUE). Set FALSE when the data will be staged later.
#' @param check_assemblies (logical) Check that the assembly files named in the
#'   mapping file exist under `assembly_path` (default TRUE).
#' @param assemble_cpus Default number of CPUs for assembly.
#' @param assemble_memory Default memory (GB) for assembly.
#' @param seeds_db Path or URL to the GetOrganelle seeds database. NULL
#'   (default) names an Actinopterygii database that RiboPilot does not ship;
#'   the Reference step replaces it.
#' @param labels_db Path or URL to the GetOrganelle labels database. NULL
#'   (default) names an Actinopterygii database that RiboPilot does not ship;
#'   the Reference step replaces it.
#' @param getOrganelle Default GetOrganelle command-line options.
#' @param max_paths Max assembly paths for a sample to continue (default 10).
#' @param max_scaffolds Max contigs within a path to continue (default 10).
#' @param min_assembly_length Minimum contig length (bp) to retain (default 500).
#' @param ref_dir Folder for built reference databases (default: the
#'   `RIBOPILOT_REF_DIR` environment variable, else [ribopilot_ref_dir()]).
#' @param mapping_geome Name of the mapping-file column holding GEOME BCIDs
#'   (optional). Stored as `GEOME_BCID`. See `vignette("Sample-Metadata")`.
#' @param fetch_geome Fetch GEOME metadata for samples with a BCID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_geome()] later.
#' @param mapping_gbif Name of the mapping-file column holding GBIF occurrence
#'   IDs (optional). Stored as `GBIF_ID`. See `vignette("Sample-Metadata")`.
#' @param fetch_gbif Fetch GBIF metadata for samples with a GBIF ID during setup
#'   (default TRUE). Set FALSE when offline and run [fetch_gbif()] later.
#' @param mapping_biosample Name of the mapping-file column holding NCBI
#'   BioSample or SRA accessions (optional). Must not be the sample ID column.
#'   Stored as `BioSample`. See `vignette("Sample-Metadata")`.
#' @param fetch_biosample Fetch NCBI metadata for samples with a BioSample
#'   value during setup (default TRUE). Set FALSE when offline and run
#'   [fetch_biosample()] later.
#' @param link_sources Follow links between GEOME, GBIF, and NCBI records to
#'   fill in IDs a sample does not have yet (default FALSE). Saved as the project
#'   setting. See `vignette("Sample-Metadata")`.
#'
#' @export
new_db <- function(
    db_path = "./.sqlite",
    mapping_fn = NULL,
    mapping_id = "ID",
    mapping_taxon = "Taxon",
    data_path = NULL,
    assembly_path = NULL,
    check_reads = TRUE,
    check_assemblies = TRUE,
    assemble_cpus = 6,
    assemble_memory = 24,
    seeds_db = NULL,
    labels_db = NULL,
    getOrganelle = paste(
      "-F anonym",
      "-R 10 -k 21,45,65,85,105,115",
      "--larger-auto-ws"
    ),
    max_paths = 10,
    max_scaffolds = 10,
    min_assembly_length = 500,
    ref_dir = NULL,
    mapping_geome = "GEOME_BCID",
    fetch_geome = TRUE,
    mapping_gbif = "GBIF_ID",
    fetch_gbif = TRUE,
    mapping_biosample = "BioSample",
    fetch_biosample = TRUE,
    link_sources = FALSE) {
  # Read and check the mapping file ----
  # Fail here, not per-task inside a work/ directory as an opaque tool error.
  mapping <- read_mapping(mapping_fn %||% "./mapping.csv", mapping_id, mapping_taxon) |>
    check_mapping(mapping_id, data_path, assembly_path, check_reads, check_assemblies,
      mapping_geome, mapping_gbif, mapping_biosample)

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
  # dbplyr is never called by name but is required: dplyr::tbl() on a DBI
  # connection loads it to translate, so dropping it from Imports breaks every
  # table write below on a machine that does not happen to have it installed.
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
      tax_alternatives = NA_character_,
      lineage = NA_character_,
      lineage_taxids = NA_character_,
      ref_id = NA_character_,
      # user_group: samples sharing a label must land in one reference group.
      # ref_source: 'auto' (proposer) or 'manual' (hand-assigned).
      user_group = NA_character_,
      ref_source = NA_character_,
      taxonomy_switch = 1L,
      taxonomy_lock = 0L,
      time_stamp = NA_integer_
    )
  mapping <- .meta_take_cols(mapping,
    c(GEOME = mapping_geome, GBIF = mapping_gbif, NCBI = mapping_biosample),
    keep = c(mapping_id, mapping_taxon))
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
        # Adapter handling plus 3' quality trimming and a length floor. Without
        # --cut_tail, degraded read ends reach GetOrganelle / SPAdes and, at rDNA
        # copy number, become spurious variant groups in screen_variants().
        # --dont_eval_duplication stays: for rDNA, duplicate reads are real copy
        # number, so they are never deduplicated.
        fastp = paste("--trim_poly_g --correction --detect_adapter_for_pe",
          "--dont_eval_duplication --cut_tail --cut_mean_quality 20 -l 50")
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
      find_rrna_notes TEXT,
      user_notes TEXT,
      assemble_switch INTEGER,
      assemble_lock INTEGER,
      assemble_opts TEXT,
      blast_opts TEXT,
      triage_opts TEXT DEFAULT 'default',
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
          find_rrna_notes = NA_character_,
          assemble_switch = 1,
          assemble_lock = 0,
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
      max_paths INTEGER,
      max_scaffolds INTEGER,
      min_assembly_length INTEGER,
      assembler TEXT DEFAULT 'getorganelle',
      find_rrna INTEGER DEFAULT 1,
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
        seeds_db = seeds_db,
        labels_db = labels_db,
        max_paths = max_paths,
        max_scaffolds = max_scaffolds,
        min_assembly_length = min_assembly_length,
        assembler = "getorganelle",
        find_rrna = 1L,
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
      remote_blast INTEGER,
      remote_fallback INTEGER,
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
        # Local-first. The shipped panel answers most scaffolds in milliseconds,
        # so going to NCBI by default would trade minutes per sample for hits the
        # panel already has. Set this to search NCBI instead of the panel.
        remote_blast = 0L,
        # ...but the panel is metazoan nuclear rDNA only, so a contaminant or an
        # off-target scaffold can be unnameable locally and still perfectly
        # nameable at NCBI. Those scaffolds escalate unless this is turned off.
        remote_fallback = 1L,
        # Significance floor for the remote search, shown and editable as "Extra
        # blastn options". blastn's own default (-evalue 10) keeps hits expected
        # by chance. Deliberately NO -perc_identity: this one search feeds both
        # the reference panel and the contamination screen, and the screen exists
        # to see the distant cross-kingdom hits an identity floor would remove.
        # The reference panel gets its own floor in parse_blast_candidates(),
        # after the screen has been written.
        extra_opts = "-evalue 1e-20",
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
    );"
  )

  # blast_hits (full retained hit distribution per unit; see BLAST_HITS_DDL) ----
  DBI::dbExecute(con, BLAST_HITS_DDL)

  # blast_searched (resume marker keyed on query hash; see BLAST_SEARCHED_DDL) ----
  DBI::dbExecute(con, BLAST_SEARCHED_DDL)

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
      clips TEXT,
      time_stamp INTEGER,
      PRIMARY KEY (ID, path, scaffold)
    );"
  )

  # assembly_backup (original pre-trim state so a trim can be undone) ----
  DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)

  # annotations (nrDNA feature model per assembly unit) ----
  # status/flags are filled by curation (curate_one); tool/edited by annotation.
  # partial_start/partial_stop: 5'/3' end partial in the feature's own orientation.
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
      partial_start INTEGER,
      partial_stop INTEGER,
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
  # load_assemblies() after assembly. n_flags / curate_flags filled by curation,
  # annotate_notes by the curator; curate_opts selects the per-clade ruleset.
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
      curate_flags TEXT,
      itsx_warnings TEXT,
      reviewed TEXT,
      problematic TEXT,
      partial TEXT DEFAULT 'no',
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
      barrnap_opts TEXT,
      pybarrnap_opts TEXT,
      PRIMARY KEY (annotate_opts)
    );"
  )
  # rrna_tool: 'pybarrnap' (default; euk-capable Python reimplementation) or
  # 'barrnap'. rrna_kingdom values differ by tool (pybarrnap: bac/arc/euk/all;
  # barrnap: bac/arc/fun). *_opts: extra CLI flags for that tool only.
  dplyr::tbl(con, "annotate_opts") |>
    dplyr::rows_upsert(
      data.frame(annotate_opts = "default", cpus = 4L, memory = 8L,
        rrna_tool = "pybarrnap", rrna_kingdom = "euk",
        itsx_taxon = "all", itsx_opts = "", barrnap_opts = "", pybarrnap_opts = ""),
      in_place = TRUE, copy = TRUE, by = "annotate_opts"
    )

  # curate_opts (per-clade curation ruleset selection) ----
  DBI::dbExecute(
    con,
    "CREATE TABLE curate_opts (
      curate_opts TEXT NOT NULL,
      ruleset TEXT,
      PRIMARY KEY (curate_opts)
    );"
  )
  dplyr::tbl(con, "curate_opts") |>
    dplyr::rows_upsert(
      data.frame(curate_opts = "default", ruleset = "metazoa_generic"),
      in_place = TRUE, copy = TRUE, by = "curate_opts"
    )

  # triage_opts (per-rule triage actions; see TRIAGE_RULE_DEFS) ----
  DBI::dbExecute(con, TRIAGE_OPTS_DDL)
  write_triage_opts(con, triage_actions_default())

  # rrna_candidates (user-assembly rDNA contig search; see RRNA_CANDIDATES_DDL) ----
  DBI::dbExecute(con, RRNA_CANDIDATES_DDL)

  # screen (contamination verdict per unit; see SCREEN_DDL) ----
  DBI::dbExecute(con, SCREEN_DDL)
  DBI::dbExecute(con, CHIMERA_DDL)
  DBI::dbExecute(con, CHIMERA_WINDOWS_DDL)

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

  # reference (clade DB groups) + reference_config ----
  DBI::dbExecute(con, REFERENCE_DDL)
  DBI::dbExecute(con, REFERENCE_CONFIG_DDL)
  DBI::dbExecute(con,
    "INSERT INTO reference_config (rowid_key, ref_dir, min_records, max_records)
     VALUES (1, ?, 1000, 25000)",
    params = list(ref_dir %||% NA_character_))

  # export_opts (named FASTA-header configurations) ----
  ensure_export_opts(con)
  meta_apply_schema(con)

  if (isTRUE(link_sources)) .meta_set_link_enabled(con, TRUE)
  .meta_fetch_new(con, mapping,
    list(GEOME = fetch_geome, GBIF = fetch_gbif, NCBI = fetch_biosample),
    link = isTRUE(link_sources))
  record_schema_event(con, "created", current_schema_version())

  invisible(return())
}

#' Sample columns new_db() adds itself; a mapping column may not reuse one
#'
#' Compared ignoring case, as SQLite compares column names.
#' @noRd
SAMPLE_BUILTIN_COLS <- c("taxid", "sci_name", "rank", "match_type", "tax_alternatives", "lineage",
  "lineage_taxids", "ref_id", "user_group", "ref_source", "taxonomy_switch",
  "taxonomy_lock", "time_stamp")

#' Read a mapping CSV as typed, or stop naming the problem and the fix
#'
#' Every column is read as text, so IDs like 001 stay as typed. Values are
#' trimmed, blank cells become NA, and fully blank rows are dropped. Column
#' names are matched to ID, Taxon, R1, R2, Assembly, GEOME_BCID, GBIF_ID, and
#' BioSample ignoring case.
#'
#' @param mapping_fn Path to the mapping CSV.
#' @param mapping_id,mapping_taxon ID and taxon column names.
#' @return the mapping data.frame.
#' @noRd
read_mapping <- function(mapping_fn, mapping_id = "ID", mapping_taxon = "Taxon") {
  if (length(mapping_fn) != 1 || is.na(mapping_fn) || !file.exists(mapping_fn) ||
      dir.exists(mapping_fn)) {
    stop("Mapping file not found: ", mapping_fn %||% "(none given)", call. = FALSE)
  }
  excel <- "In Excel, use File > Save As > CSV UTF-8 (Comma delimited)."
  if (any(readBin(mapping_fn, "raw", 512L) == as.raw(0))) {
    stop("The mapping file is not a plain-text CSV (it looks like an Excel workbook ",
      "or another binary file): ", mapping_fn, ". ", excel, call. = FALSE)
  }
  hdr <- readLines(mapping_fn, n = 1L, warn = FALSE)
  if (length(hdr) == 0 || !nzchar(trimws(hdr))) {
    stop("The mapping file is empty: ", mapping_fn, call. = FALSE)
  }
  if (!grepl(",", hdr) && grepl("[;\t]", hdr)) {
    stop("The mapping file must be comma-separated, but its header uses ",
      if (grepl(";", hdr)) "semicolons" else "tabs", ": ", mapping_fn, ". ", excel,
      call. = FALSE)
  }
  mapping <- utils::read.csv(mapping_fn, colClasses = "character", na.strings = c("", "NA"))
  mapping[] <- lapply(mapping, function(x) {
    x <- trimws(x)
    x[which(x == "")] <- NA_character_
    x
  })

  for (w in setdiff(c(mapping_id, mapping_taxon, "R1", "R2", "Assembly", META_SAMPLE_COLS), names(mapping))) {
    i <- which(tolower(names(mapping)) == tolower(w))
    if (length(i) == 1) {
      message("Using mapping column '", names(mapping)[i], "' as '", w, "'.")
      names(mapping)[i] <- w
    }
  }
  if ("Assembly" %in% names(mapping)) {
    for (w in setdiff(c("R1", "R2"), names(mapping))) mapping[[w]] <- NA_character_
  }
  miss <- setdiff(c(mapping_id, mapping_taxon, "R1", "R2"), names(mapping))
  if (length(miss)) {
    stop("The mapping file is missing column(s): ", paste(miss, collapse = ", "), ".\n",
      "Columns found: ", paste(names(mapping), collapse = ", "), ".\n",
      "It needs ", mapping_id, " and ", mapping_taxon, ", plus R1 and R2 (paired reads) ",
      "or Assembly, in a comma-separated file.", call. = FALSE)
  }
  low <- tolower(names(mapping))
  twins <- names(mapping)[low %in% low[duplicated(low)]]
  if (length(twins)) {
    stop("Mapping column names must differ by more than upper/lower case: ",
      paste(twins, collapse = ", "), ".", call. = FALSE)
  }
  reserved <- c(SAMPLE_BUILTIN_COLS, if (mapping_id != "ID") "ID",
    if (mapping_taxon != "Taxon") "Taxon")
  hit <- names(mapping)[low %in% tolower(reserved)]
  if (length(hit)) {
    stop("The mapping file has column(s) RiboPilot fills in itself: ",
      paste(hit, collapse = ", "), ". Rename them (for example my_", tolower(hit[1]), ").",
      call. = FALSE)
  }

  blank <- rowSums(!is.na(mapping)) == 0
  no_id <- is.na(mapping[[mapping_id]]) & !blank
  if (any(no_id)) {
    stop("Some rows have no ", mapping_id, " (blank or NA): line(s) ",
      paste(truncate_ids(which(no_id) + 1L), collapse = ", "), " of the mapping file.",
      call. = FALSE)
  }
  mapping <- mapping[!blank, , drop = FALSE]
  rownames(mapping) <- NULL
  if (nrow(mapping) == 0) {
    stop("The mapping file has a header row but no samples: ", mapping_fn, call. = FALSE)
  }

  ids <- mapping[[mapping_id]]
  k <- tolower(ids)
  dup <- unique(ids[k %in% k[duplicated(k)]])
  if (length(dup)) {
    stop("Each ", mapping_id, " must be unique, ignoring upper/lower case (sample folders ",
      "would clash on macOS and Windows). Repeated: ",
      paste(truncate_ids(dup), collapse = ", "), call. = FALSE)
  }
  bad <- ids[!grepl("^[A-Za-z0-9_:-]+$", ids)]
  if (length(bad)) {
    stop("IDs may use only letters, digits, dash (-), underscore (_), and colon (:), with ",
      "no spaces, dots, or slashes. Change: ",
      paste(truncate_ids(sQuote(bad, FALSE)), collapse = ", "), call. = FALSE)
  }
  long <- ids[nchar(ids) > 18]
  if (length(long)) {
    stop("IDs must be 18 characters or fewer. Shorten: ",
      paste(truncate_ids(long), collapse = ", "), call. = FALSE)
  }
  mapping
}

#' Check a mapping's read and assembly columns against the files on disk
#'
#' @return the mapping, with a normalized `Assembly` column.
#' @noRd
check_mapping <- function(mapping, mapping_id = "ID", data_path = NULL, assembly_path = NULL,
                          check_reads = TRUE, check_assemblies = TRUE,
                          mapping_geome = "GEOME_BCID", mapping_gbif = "GBIF_ID",
                          mapping_biosample = "BioSample") {
  ids <- mapping[[mapping_id]]
  src_cols <- c(GEOME = mapping_geome, GBIF = mapping_gbif, NCBI = mapping_biosample)
  for (src in names(src_cols)[src_cols == mapping_id]) {
    stop("The sample ID column '", mapping_id, "' can't also be the ", src,
      " ID column; copy the IDs into their own column.", call. = FALSE)
  }
  for (src in names(src_cols)) {
    col <- src_cols[[src]]
    std <- META_SOURCES[[src]]$col
    if (col != std && std %in% names(mapping)) {
      stop("The mapping file has both a ", col, " column and a ", std,
        " column for the same source. Keep one of them.", call. = FALSE)
    }
    if (!col %in% names(mapping)) {
      if (col != std) {
        stop("The ", META_SOURCES[[src]]$label, " column '", col, "' was not found in the mapping file.",
          call. = FALSE)
      }
      next
    }
    raw <- trimws(.meta_chr(mapping[[col]]))
    raw[is.na(raw)] <- ""
    bad <- nzchar(raw) & is.na(META_SOURCES[[src]]$normalize(raw))
    if (any(bad)) {
      warning(structure(class = c("ribopilot_meta_id_warning", "warning", "condition"),
        list(message = paste0("mapping ", META_SOURCES[[src]]$label, " ", META_SOURCES[[src]]$id_label,
          ": not a valid ID for ", .lst(ids[bad]), "; these samples will show a failed ",
          META_SOURCES[[src]]$label, " fetch"), call = NULL)))
    }
  }
  stop_unsafe_paths(c(
    stats::setNames(mapping$R1, paste(ids, "R1")),
    stats::setNames(mapping$R2, paste(ids, "R2")),
    if ("Assembly" %in% names(mapping)) stats::setNames(mapping$Assembly, paste(ids, "Assembly"))),
    mount = FALSE)
  mapping <- normalize_assembly_column(mapping, mapping_id, assembly_path)
  if (isTRUE(check_reads)) check_mapping_reads(mapping, mapping_id, data_path)
  if (isTRUE(check_assemblies)) check_mapping_assemblies(mapping, mapping_id, assembly_path)
  mapping
}

#' Stop when a path cannot hold a new project
#'
#' @return (invisibly) the path, made absolute.
#' @noRd
check_project_path <- function(path, force = FALSE) {
  path <- normalizePath(path, mustWork = FALSE)
  if (!startsWith(path, "/")) path <- file.path(getwd(), path)
  if (file.exists(path) && !dir.exists(path)) {
    stop("The project path is a file, not a folder: ", path, call. = FALSE)
  }
  if (file.exists(file.path(path, ".sqlite")) && !force) {
    stop("A RiboPilot project already exists in ", path, ". Use force = TRUE to ",
      "replace it (its results will be lost).", call. = FALSE)
  }
  top <- path
  while (!dir.exists(top) && dirname(top) != top) top <- dirname(top)
  if (unname(file.access(top, 2)) != 0) {
    stop("Cannot write to ", top, ": permission denied.", call. = FALSE)
  }
  stop_unsafe_paths(c(path = path))
  invisible(path)
}

#' Stop on paths the pipeline passes unquoted to Docker, Singularity, and bash
#'
#' @param x Named character vector; names say what each path is. NA, "" and
#'   URLs are skipped.
#' @param mount Also refuse "," and ":", which split a container mount.
#' @noRd
stop_unsafe_paths <- function(x, mount = TRUE) {
  x <- x[!is.na(x) & nzchar(x) & !grepl("^[a-z]+://", x)]
  re <- if (mount) {
    "[\\s'\"`$\\\\;&|<>()*?#!{}\\[\\],:]"
  } else {
    "[\\s'\"`$\\\\;&|<>()*?#!{}\\[\\]]"
  }
  hit <- regmatches(x, gregexpr(re, x, perl = TRUE))
  bad <- which(lengths(hit) > 0)
  if (length(bad) == 0) return(invisible())
  lines <- vapply(bad, function(i) {
    ch <- unique(hit[[i]])
    what <- unique(ifelse(grepl("\\s", ch, perl = TRUE), "space", sQuote(ch, FALSE)))
    sprintf("  %s: '%s' (remove: %s)", names(x)[i], x[i], paste(what, collapse = ", "))
  }, character(1))
  stop("These paths have characters the pipeline cannot pass to Docker, Singularity, ",
    "or Nextflow:\n", paste(truncate_ids(lines), collapse = "\n"),
    "\nRename them to use only letters, digits, and . _ - / characters.", call. = FALSE)
}

#' Check that the read files named in a mapping table exist
#'
#' @param mapping Mapping data frame (needs `R1` and `R2` columns).
#' @param mapping_id Column name holding the sample ids.
#' @param data_path Raw data directory the read names are relative to. When
#'   NULL, only the paired-end check runs.
#' @noRd
check_mapping_reads <- function(mapping, mapping_id = "ID", data_path = NULL) {
  ids <- as.character(mapping[[mapping_id]])
  blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))

  # A sample supplying its own assembly may name no reads at all. It may still
  # name a pair, which is used for coverage only; half a pair is an error either
  # way, since single-end reads are not supported.
  has_asmb <- if ("Assembly" %in% colnames(mapping)) !blank(mapping$Assembly) else rep(FALSE, nrow(mapping))
  no_reads <- blank(mapping$R1) & blank(mapping$R2)

  missing_input <- no_reads & !has_asmb
  if (any(missing_input)) {
    stop(
      "Every sample needs either paired reads or an assembly, but these have neither:\n",
      paste(truncate_ids(ids[missing_input]), collapse = "\n")
    )
  }

  unpaired <- (blank(mapping$R1) | blank(mapping$R2)) & !no_reads
  if (any(unpaired)) {
    stop(
      "Every sample with reads needs both an R1 and an R2 file, but these are missing one:\n",
      paste(truncate_ids(ids[unpaired]), collapse = "\n"),
      "\nSingle-end reads are not supported by the RiboPilot pipeline."
    )
  }
  same <- !no_reads & basename(mapping$R1) == basename(mapping$R2)
  if (any(same)) {
    stop(
      "R1 and R2 must name two different read files, but these use one file name for both:\n",
      paste(truncate_ids(paste0(ids[same], ": ", mapping$R1[same])), collapse = "\n")
    )
  }

  mapping <- mapping[!no_reads, , drop = FALSE]
  ids <- ids[!no_reads]
  if (nrow(mapping) == 0) {
    return(invisible())
  }

  if (length(data_path) != 1 || is.na(data_path) || !nzchar(data_path)) {
    return(invisible())
  }
  if (!dir.exists(data_path)) {
    stop("data_path does not exist or is not a folder: ", data_path)
  }
  files <- as.vector(rbind(mapping$R1, mapping$R2))
  full <- grepl("^[/~]", files)
  if (any(full)) {
    stop(
      "R1 and R2 take file names relative to data_path (", data_path, "), not full paths:\n",
      paste(truncate_ids(paste0(rep(ids, each = 2)[full], ": ", files[full])), collapse = "\n")
    )
  }
  reads <- data.frame(
    ID = rep(ids, each = 2),
    path = file.path(data_path, files)
  )
  bad <- reads[!file.exists(reads$path), , drop = FALSE]
  if (nrow(bad)) {
    stop(
      "Read files named in the mapping file were not found in ", data_path, ":\n",
      paste(truncate_ids(paste0(bad$ID, ": ", bad$path)), collapse = "\n")
    )
  }
  invisible()
}

#' Trim a list of problem entries to a readable length
#' @noRd
truncate_ids <- function(x, n = 10) {
  if (length(x) <= n) {
    return(x)
  }
  c(x[seq_len(n)], sprintf("... and %d more", length(x) - n))
}

#' Ensure the mapping has an Assembly column, and that it is usable
#'
#' The column is optional in the user's file but always present in the samples
#' table, so the pipeline's sample query can name it unconditionally.
#'
#' @param mapping Mapping data.frame.
#' @param mapping_id ID column name.
#' @param assembly_path Directory the assembly files live in, or NULL.
#' @return the mapping, with a normalized `Assembly` column.
#' @noRd
normalize_assembly_column <- function(mapping, mapping_id = "ID", assembly_path = NULL) {
  if (!("Assembly" %in% colnames(mapping))) {
    mapping$Assembly <- NA_character_
    return(mapping)
  }
  mapping$Assembly <- trimws(as.character(mapping$Assembly))
  mapping$Assembly[!nzchar(mapping$Assembly)] <- NA_character_
  named <- !is.na(mapping$Assembly)
  if (any(named) && (length(assembly_path) != 1 || is.na(assembly_path) || !nzchar(assembly_path))) {
    stop(
      "The mapping file names assembly files, but no assembly_path was given.\n",
      "Pass assembly_path = the directory holding them."
    )
  }
  mapping
}

#' Check that the assembly files named in the mapping file exist
#'
#' Fail here, not per-task inside a work/ directory as a dangling symlink.
#'
#' @param mapping Mapping data.frame (already normalized).
#' @param mapping_id ID column name.
#' @param assembly_path Directory the assembly files live in, or NULL.
#' @noRd
check_mapping_assemblies <- function(mapping, mapping_id = "ID", assembly_path = NULL) {
  if (!("Assembly" %in% colnames(mapping))) {
    return(invisible())
  }
  named <- !is.na(mapping$Assembly) & nzchar(mapping$Assembly)
  if (!any(named) || length(assembly_path) != 1 || !nzchar(assembly_path)) {
    return(invisible())
  }
  if (!dir.exists(assembly_path)) {
    stop("assembly_path does not exist or is not a folder: ", assembly_path)
  }
  ids <- as.character(mapping[[mapping_id]])[named]
  files <- mapping$Assembly[named]
  full <- grepl("^[/~]", files)
  if (any(full)) {
    stop(
      "Assembly takes file names relative to assembly_path (", assembly_path, "), not full paths:\n",
      paste(truncate_ids(paste0(ids[full], ": ", files[full])), collapse = "\n")
    )
  }
  paths <- file.path(assembly_path, files)
  bad <- !file.exists(paths)
  if (any(bad)) {
    stop(
      "Assembly files named in the mapping file were not found in ", assembly_path, ":\n",
      paste(truncate_ids(paste0(ids[bad], ": ", paths[bad])), collapse = "\n")
    )
  }
  gz <- grepl("\\.gz$", files, ignore.case = TRUE)
  if (any(gz) && !all(gz)) {
    stop(
      "Assembly files must be all gzipped (.gz) or all plain in one project; the pipeline ",
      "cannot mix them. Gzipped: ", paste(truncate_ids(ids[gz]), collapse = ", "),
      ". Plain: ", paste(truncate_ids(ids[!gz]), collapse = ", "), "."
    )
  }
  first <- vapply(paths, function(p) {
    con <- gzfile(p, "rt")
    on.exit(close(con))
    l <- suppressWarnings(readLines(con, n = 1L, warn = FALSE))
    if (length(l)) l else ""
  }, character(1), USE.NAMES = FALSE)
  notfa <- !startsWith(first, ">")
  if (any(notfa)) {
    stop(
      "Assembly files must be FASTA (first line starts with '>'), but these are not:\n",
      paste(truncate_ids(paste0(ids[notfa], ": ", paths[notfa])), collapse = "\n")
    )
  }
  invisible()
}
