# Backwards-compatible DB migration. Brings a project database created by an
# older RiboPilot up to the current schema. Migrations are additive and
# idempotent (guarded ALTER/CREATE), tracked by SQLite's `PRAGMA user_version`.
# new_db() stamps the current version; migrate_db() upgrades older projects.

#' Current project-DB schema version. Bump when adding a migration step below.
#' @noRd
RIBOPILOT_SCHEMA_VERSION <- 15L

#' DDL for blast_hits: the full retained BLAST hit distribution per unit.
#'
#' One row per (subject, HSP). The pipeline already pays for a remote core_nt
#' search and used to keep a single row of up to 50 subjects; everything a
#' contamination call needs -- the hit distribution, the on/off-target bitscore
#' margin, and the positional structure that separates a chimera from a clean
#' scaffold -- was discarded. Shared by new_db() and the v9 migration.
#' @noRd
BLAST_HITS_DDL <- "CREATE TABLE IF NOT EXISTS blast_hits (
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
  is_nuclear_rrna INTEGER,
  hit_class TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold, hit_rank, hsp)
);"

#' DDL for the triage_opts table (per-rule action presets).
#' Shared by new_db() and the v12 migration so both agree on the schema.
#' @noRd
TRIAGE_OPTS_DDL <- "CREATE TABLE IF NOT EXISTS triage_opts (
  triage_opts            TEXT NOT NULL,
  rule_not_annotated     TEXT,
  rule_no_screen         TEXT,
  rule_no_feat_non_rdna  TEXT,
  rule_no_feat_rdna_hit  TEXT,
  rule_feat_non_rdna_hit TEXT,
  rule_odd_one_out       TEXT,
  rule_screen_off_target TEXT,
  PRIMARY KEY (triage_opts)
);"

#' Seed row for triage_opts, matching TRIAGE_RULE_DEFS defaults.
#' @noRd
TRIAGE_OPTS_SEED <- "INSERT INTO triage_opts
  (triage_opts, rule_not_annotated, rule_no_screen, rule_no_feat_non_rdna,
   rule_no_feat_rdna_hit, rule_feat_non_rdna_hit, rule_odd_one_out,
   rule_screen_off_target)
  VALUES ('default', 'flag', 'flag', 'drop', 'flag', 'flag', 'flag', 'flag');"

#' DDL for the screen table (machine contamination verdict, one row per unit).
#' Shared by new_db() and the v14 migration so both agree on the schema.
#' @noRd
SCREEN_DDL <- "CREATE TABLE IF NOT EXISTS screen (
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
);"

#' DDL for the assembly_backup table (original pre-trim state, for undo/restore).
#' Shared by new_db() and the v6 migration so both agree on the schema.
#' @noRd
ASSEMBLY_BACKUP_DDL <- "CREATE TABLE IF NOT EXISTS assembly_backup (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  sequence TEXT,
  length INTEGER,
  annotations TEXT,
  coverage TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' The pre-v4 restrictive remote-BLAST entrez_query.
#'
#' Retired because it restricted the search space, not just the ranking, so the
#' contamination screen could only ever see nuclear rDNA. Matched literally so
#' the v4 migration only clears it where the user never customized it.
#' @noRd
LEGACY_ENTREZ_QUERY <- paste(
  "(internal transcribed spacer[Title] OR 18S ribosomal RNA[Title]",
  "OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title]",
  "OR small subunit ribosomal RNA[Title]",
  "OR large subunit ribosomal RNA[Title])",
  "NOT (mitochondrion[Title] OR chloroplast[Title] OR plastid[Title])")

#' Bring a project database up to the current schema (idempotent, additive)
#'
#' Adds any tables/columns introduced after a project was created, without
#' touching existing data. Safe to call repeatedly. New columns are only ever
#' added, never dropped or retyped.
#'
#' @param con Open DBI connection to a project database.
#' @return (invisibly) the schema version after migration.
#' @export
migrate_db <- function(con) {
  ver <- tryCatch(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], error = function(e) 0L)
  has_tbl <- function(tbl) tbl %in% DBI::dbListTables(con)
  has_col <- function(tbl, col) has_tbl(tbl) && col %in% DBI::dbListFields(con, tbl)
  add_col <- function(tbl, col, decl) {
    if (has_tbl(tbl) && !has_col(tbl, col)) {
      DBI::dbExecute(con, sprintf("ALTER TABLE %s ADD COLUMN %s %s", tbl, col, decl))
    }
  }

  # --- schema v1: reference-guided curation (curate_opts + status/flags) ---
  if (!has_tbl("curate_opts")) {
    DBI::dbExecute(con,
      "CREATE TABLE curate_opts (
        curate_opts TEXT NOT NULL,
        ruleset TEXT,
        min_anchor_pid REAL,
        check_contaminant INTEGER,
        PRIMARY KEY (curate_opts)
      );")
  }
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM curate_opts")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO curate_opts (curate_opts, ruleset, min_anchor_pid, check_contaminant)
       VALUES ('default', 'metazoa_generic', 85, 1)")
  }
  add_col("annotate", "curate_opts", "TEXT DEFAULT 'default'")
  add_col("annotate", "n_flags", "INTEGER")
  add_col("annotations", "status", "TEXT")
  add_col("annotations", "flags", "TEXT")

  # --- schema v2: editable annotation options (annotate_opts) ---
  if (!has_tbl("annotate_opts")) {
    DBI::dbExecute(con,
      "CREATE TABLE annotate_opts (
        annotate_opts TEXT NOT NULL,
        cpus INTEGER,
        memory INTEGER,
        rrna_tool TEXT,
        rrna_kingdom TEXT,
        itsx_taxon TEXT,
        itsx_opts TEXT,
        PRIMARY KEY (annotate_opts)
      );")
  }
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM annotate_opts")$n == 0) {
    DBI::dbExecute(con,
      "INSERT INTO annotate_opts (annotate_opts, cpus, memory, rrna_tool,
         rrna_kingdom, itsx_taxon, itsx_opts)
       VALUES ('default', 4, 8, 'pybarrnap', 'euk', 'all', '')")
  }
  add_col("annotate", "annotate_opts", "TEXT DEFAULT 'default'")

  # --- schema v3: in-app taxonomy stage (display + lock/state on samples) ---
  had_tax <- has_col("samples", "taxonomy_lock")
  add_col("samples", "sci_name", "TEXT")
  add_col("samples", "rank", "TEXT")
  add_col("samples", "match_type", "TEXT")
  add_col("samples", "taxonomy_switch", "INTEGER DEFAULT 1")
  add_col("samples", "taxonomy_lock", "INTEGER DEFAULT 0")
  if (!had_tax && has_tbl("samples")) {
    # Grandfather samples from projects that predate the taxonomy stage so they
    # are not hidden from Assemble; mark them assigned + locked.
    DBI::dbExecute(con, "UPDATE samples SET taxonomy_lock = 1, taxonomy_switch = 2")
  }

  # --- schema v4: contamination screen over the unrestricted remote BLAST ---
  if (!has_tbl("blast_screen")) {
    DBI::dbExecute(con,
      "CREATE TABLE blast_screen (
        ID TEXT NOT NULL,
        accession TEXT,
        title TEXT,
        species TEXT,
        staxid INTEGER,
        lineage TEXT,
        pident REAL,
        qcovs REAL,
        evalue REAL,
        is_nuclear_rrna INTEGER,
        time_stamp INTEGER,
        PRIMARY KEY (ID)
      );")
  }
  # ITSx chimera flag surfaced through curation.
  add_col("annotate", "chimeric", "INTEGER")
  # Fraction of the rRNA model an anchor covers (truncation signal).
  add_col("annotations", "model_coverage", "REAL")
  # Fraction of a sample element covered by its reference alignment.
  add_col("blast_ref_alignment", "aln_cov", "REAL")
  # Curator flags beyond `reviewed`: problematic (exclude from export) + partial.
  add_col("annotate", "problematic", "TEXT")
  add_col("annotate", "partial", "TEXT")
  # Retire the restrictive entrez_query, but only where it is still verbatim the
  # old default -- a query the user chose is theirs to keep.
  if (has_tbl("blast_opts")) {
    DBI::dbExecute(con,
      "UPDATE blast_opts SET entrez_query = '' WHERE TRIM(COALESCE(entrez_query,'')) = TRIM(?)",
      params = list(LEGACY_ENTREZ_QUERY))
    # 5 was tuned for a pre-filtered search; the client-side filter needs depth.
    DBI::dbExecute(con,
      "UPDATE blast_opts SET max_target_seqs = 50 WHERE COALESCE(max_target_seqs,5) = 5")
  }

  # --- schema v5: taxonomy-stage timestamp on samples ("Last Updated" column) ---
  add_col("samples", "time_stamp", "INTEGER")

  # --- schema v6: assembly_backup (original pre-trim state, for restore) ---
  if (!has_tbl("assembly_backup")) DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)

  # --- schema v7: per-(ID,path,scaffold) assembly units (drop the stitch model) ---
  # annotate/annotations/coverage move from per-sample to per-unit. Existing
  # stitch-based annotations + coverage are invalid under the per-scaffold model
  # (their coords are on the stitched sequence, not a single scaffold), so they
  # are cleared and annotation is re-queued (annotate_switch = 1).
  if (has_tbl("annotate") && !has_col("annotate", "path")) {
    DBI::dbExecute(con, "CREATE TABLE annotate_new (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        annotate_switch INTEGER, annotate_lock INTEGER, n_features INTEGER, missing TEXT,
        chimeric INTEGER, annotate_opts TEXT, curate_opts TEXT, n_flags INTEGER,
        annotate_notes TEXT, reviewed TEXT, problematic TEXT, partial TEXT,
        export_group TEXT, time_stamp INTEGER, PRIMARY KEY (ID, path, scaffold));")
    if (has_tbl("assemblies")) {
      DBI::dbExecute(con,
        "INSERT OR IGNORE INTO annotate_new
           (ID, path, scaffold, annotate_switch, annotate_lock, annotate_opts, curate_opts, reviewed)
         SELECT a.ID, a.path, a.scaffold, 1, 0,
                COALESCE(o.annotate_opts,'default'), COALESCE(o.curate_opts,'default'), 'no'
         FROM assemblies a LEFT JOIN annotate o ON o.ID = a.ID
         WHERE COALESCE(a.ignore, 0) = 0;")
    }
    # samples with no (non-ignored) assemblies keep a placeholder (1,1) unit
    DBI::dbExecute(con,
      "INSERT OR IGNORE INTO annotate_new
         (ID, path, scaffold, annotate_switch, annotate_lock, annotate_opts, curate_opts, reviewed)
       SELECT o.ID, 1, 1, 1, 0,
              COALESCE(o.annotate_opts,'default'), COALESCE(o.curate_opts,'default'), 'no'
       FROM annotate o;")
    if (has_col("samples", "export_group")) {
      DBI::dbExecute(con,
        "UPDATE annotate_new SET export_group =
           (SELECT s.export_group FROM samples s WHERE s.ID = annotate_new.ID)
         WHERE EXISTS (SELECT 1 FROM samples s WHERE s.ID = annotate_new.ID
                       AND TRIM(COALESCE(s.export_group,'')) != '');")
    }
    DBI::dbExecute(con, "DROP TABLE annotate;")
    DBI::dbExecute(con, "ALTER TABLE annotate_new RENAME TO annotate;")
  }
  if (has_tbl("annotations") && !has_col("annotations", "path")) {
    DBI::dbExecute(con, "DROP TABLE annotations;")
    DBI::dbExecute(con, "CREATE TABLE annotations (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        element TEXT NOT NULL, type TEXT, start INTEGER, stop INTEGER, length INTEGER, strand TEXT,
        model_coverage REAL, tool TEXT, status TEXT, flags TEXT, edited INTEGER, notes TEXT,
        time_stamp INTEGER, PRIMARY KEY (ID, path, scaffold, element));")
    if (has_tbl("annotate")) DBI::dbExecute(con, "UPDATE annotate SET annotate_switch = 1;")
    # The GenBank screen + reference alignment were computed on the (now invalid)
    # stitched sequence and keyed per sample; clear them so curation does not
    # mis-attribute one sample-level hit to every scaffold unit. They are
    # repopulated per unit by the (re-added) per-scaffold BLAST step.
    for (t in c("blast_screen", "blast_ref_alignment", "blast_ref_candidates")) {
      if (has_tbl(t)) DBI::dbExecute(con, sprintf("DELETE FROM %s;", t))
    }
  }
  if (has_tbl("coverage") && !has_col("coverage", "path")) {
    DBI::dbExecute(con, "DROP TABLE coverage;")
    DBI::dbExecute(con, "CREATE TABLE coverage (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        mean_depth REAL, min_depth INTEGER, max_depth INTEGER, depths TEXT, time_stamp INTEGER,
        PRIMARY KEY (ID, path, scaffold));")
  }
  # assembly_backup (trim undo state) becomes per-unit; it is transient, so just
  # rebuild it empty when it predates the unit key.
  if (has_tbl("assembly_backup") && !has_col("assembly_backup", "path")) {
    DBI::dbExecute(con, "DROP TABLE assembly_backup;")
    DBI::dbExecute(con, ASSEMBLY_BACKUP_DDL)
  }

  # --- schema v8: per-unit GenBank BLAST reference tables (rebuild empty) ---
  # blast_ref_candidates / blast_screen / blast_ref_alignment move to
  # per-(ID,path,scaffold). v7 already cleared them, so recreate with the new keys.
  if (has_tbl("blast_ref_candidates") && !has_col("blast_ref_candidates", "path")) {
    DBI::dbExecute(con, "DROP TABLE blast_ref_candidates;")
    DBI::dbExecute(con, "CREATE TABLE blast_ref_candidates (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        rank INTEGER NOT NULL, accession TEXT, species TEXT, title TEXT, pident REAL, qcovs REAL,
        evalue REAL, time_stamp INTEGER, PRIMARY KEY (ID, path, scaffold, rank));")
  }
  if (has_tbl("blast_screen") && !has_col("blast_screen", "path")) {
    DBI::dbExecute(con, "DROP TABLE blast_screen;")
    DBI::dbExecute(con, "CREATE TABLE blast_screen (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        accession TEXT, title TEXT, species TEXT, staxid INTEGER, lineage TEXT, pident REAL, qcovs REAL,
        evalue REAL, is_nuclear_rrna INTEGER, time_stamp INTEGER, PRIMARY KEY (ID, path, scaffold));")
  }
  if (has_tbl("blast_ref_alignment") && !has_col("blast_ref_alignment", "path")) {
    DBI::dbExecute(con, "DROP TABLE blast_ref_alignment;")
    DBI::dbExecute(con, "CREATE TABLE blast_ref_alignment (
        ID TEXT NOT NULL, path INTEGER NOT NULL DEFAULT 1, scaffold INTEGER NOT NULL DEFAULT 1,
        accession TEXT NOT NULL, element TEXT NOT NULL, pid REAL, aln_cov REAL, aln_len INTEGER,
        sample_aln TEXT, ref_aln TEXT, time_stamp INTEGER,
        PRIMARY KEY (ID, path, scaffold, accession, element));")
  }

  # --- schema v9: retained BLAST hit distribution + rank-aligned lineages ---
  DBI::dbExecute(con, BLAST_HITS_DDL)
  # Which search wrote assemblies.blast_*: the local clade seed DB (host-clade
  # by construction, so its identity carries no contamination signal) or the
  # unrestricted remote core_nt search. Previously indistinguishable.
  add_col("assemblies", "blast_source", "TEXT")
  # Rank-aligned taxid path for the frozen sample lineage, so rank_distance()
  # can compare on taxids rather than name tokens.
  add_col("samples", "lineage_taxids", "TEXT")
  add_col("blast_screen", "lineage_taxids", "TEXT")
  add_col("blast_screen", "rank_distance", "INTEGER")
  add_col("blast_screen", "bitscore", "REAL")
  add_col("blast_ref_candidates", "bitscore", "REAL")
  # blast_hits is the ranking source of truth from v9 on; drop any partial rows
  # written by an interrupted pre-v9 run so a re-BLAST repopulates cleanly.
  if (ver < 9L && has_tbl("blast_hits")) DBI::dbExecute(con, "DELETE FROM blast_hits;")

  # --- schema v10: median per-position depth ---
  # The relative-depth comparison in the app is meant to reflect abundance, but a
  # mean over an rDNA scaffold is inflated by the conserved cores in proportion
  # to how much core the scaffold contains. The per-position track is already
  # stored, so existing projects can be back-filled without re-running coverage.
  add_col("coverage", "median_depth", "REAL")

  # --- schema v11: tri-state per-scaffold triage ---
  # `ignore` is a binary with no "not looked at yet" state and no record of why,
  # so a fresh sample and a fully reviewed one are indistinguishable and a
  # machine-made drop is indistinguishable from the curator's. triage carries the
  # decision, triage_reason carries the justification in plain English, and
  # triage_source says who made it so an auto re-run never overwrites a human.
  # `ignore` stays the single downstream gate (WF2, the app and export all filter
  # on it) and becomes derived from triage rather than hand-written.
  add_col("assemblies", "triage", "TEXT")
  add_col("assemblies", "triage_reason", "TEXT")
  add_col("assemblies", "triage_source", "TEXT")
  # Existing projects: an ignored scaffold was a deliberate human decision, so
  # adopt it as such rather than leaving it undecided and re-triaging over it.
  if (ver < 11L && has_col("assemblies", "triage")) {
    DBI::dbExecute(con,
      "UPDATE assemblies SET triage = 'drop', triage_source = 'user',
         triage_reason = 'dropped before triage was recorded'
       WHERE COALESCE(ignore,0) = 1 AND triage IS NULL")
  }

  # --- schema v12: named, individually controllable triage rules ---
  # The rules were an anonymous if/else chain with no per-rule setting, so the
  # WF1 triage process and the app could not be made to agree on anything the
  # curator changed. triage_opts holds one action per rule; triage_rule records
  # which rule produced a unit's current verdict.
  if (!has_tbl("triage_opts")) DBI::dbExecute(con, TRIAGE_OPTS_DDL)
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM triage_opts")$n == 0) {
    DBI::dbExecute(con, TRIAGE_OPTS_SEED)
  }
  add_col("assemble", "triage_opts", "TEXT DEFAULT 'default'")
  add_col("assemblies", "triage_rule", "TEXT")
  # unit_evidence() hard-selects coverage.depths for the sparkline, so its
  # absence would break the whole scaffold query rather than one cell
  add_col("coverage", "depths", "TEXT")
  if (ver < 12L && has_col("assemble", "triage_opts")) {
    DBI::dbExecute(con, "UPDATE assemble SET triage_opts = 'default' WHERE triage_opts IS NULL")
  }

  # --- schema v13: three-state hit classification ---
  # is_nuclear_rrna was a literal-phrase regex over ONE hit's title, so a
  # chromosome-level host record read as "not nuclear rDNA". Classification is a
  # pure function of the stored title, so history is repairable offline: the
  # backfill costs no remote BLAST, which matters because a real run takes hours.
  add_col("blast_hits", "hit_class", "TEXT")
  # gated on the version, not just on the table existing: parse_blast_candidates()
  # calls migrate_db() once per unit, so an ungated backfill would rescan the
  # whole project on every unit, O(units^2) across the run this repair exists to
  # avoid repeating. Rows written after v13 get their class from the parse path.
  if (ver < 13L && has_tbl("blast_hits")) backfill_hit_class(con)

  # --- schema v14: the screen verdict layer ---
  # The machine's opinion about a scaffold lived only as a substring inside
  # annotate_notes, which curate_one() overwrites wholesale and the modal
  # textarea also writes to. Created empty on purpose: scoring needs a resolved
  # sample lineage and a taxdump, which a migration cannot assume.
  if (!has_tbl("screen")) DBI::dbExecute(con, SCREEN_DDL)
  add_col("samples", "reference_support", "INTEGER")
  add_col("triage_opts", "rule_screen_off_target", "TEXT")
  if (ver < 14L && has_col("triage_opts", "rule_screen_off_target")) {
    DBI::dbExecute(con,
      "UPDATE triage_opts SET rule_screen_off_target = 'flag'
       WHERE rule_screen_off_target IS NULL")
  }

  # --- schema v15: intragenomic isomers and paralogs ---
  # Orthogonal to the verdict: a paralogous variant is emphatically on_target,
  # so the finding gets its own columns rather than a new verdict string. On the
  # screen table because the finding is per unit and the recompute path already
  # exists; a second table would need its own version gate kept in step with
  # this one. No backfill: SCREEN_VERSION at 3L rewrites every stored row on the
  # next run, which is what fills these.
  add_col("screen", "variant_class", "TEXT")
  add_col("screen", "variant_group", "TEXT")
  add_col("screen", "variant_identity", "REAL")

  if (ver < RIBOPILOT_SCHEMA_VERSION) {
    DBI::dbExecute(con, sprintf("PRAGMA user_version = %d", RIBOPILOT_SCHEMA_VERSION))
  }
  invisible(RIBOPILOT_SCHEMA_VERSION)
}

#' Back-fill `coverage.median_depth` from the stored per-position tracks
#'
#' Computed lazily rather than inside [migrate_db()]: a depth track is one
#' comma-joined integer per base, so a large project holds hundreds of MB of
#' them and parsing every row would make opening the app crawl. Rows written
#' before schema v10 keep a NULL median and fall back to the mean until this
#' runs.
#'
#' @param con Open project DB connection.
#' @param limit Maximum rows to fill in one call (NULL = all).
#' @return (invisibly) the number of rows filled.
#' @export
backfill_coverage_median <- function(con, limit = NULL) {
  if (!("median_depth" %in% DBI::dbListFields(con, "coverage"))) {
    return(invisible(0L))
  }
  q <- "SELECT ID, path, scaffold, depths FROM coverage
        WHERE median_depth IS NULL AND depths IS NOT NULL AND depths != ''"
  if (!is.null(limit)) q <- paste(q, "LIMIT", as.integer(limit))
  d <- tryCatch(DBI::dbGetQuery(con, q), error = function(e) data.frame())
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  upd <- DBI::dbSendStatement(con,
    "UPDATE coverage SET median_depth = ? WHERE ID = ? AND path = ? AND scaffold = ?")
  on.exit(DBI::dbClearResult(upd), add = TRUE)
  n <- 0L
  for (i in seq_len(nrow(d))) {
    v <- suppressWarnings(as.numeric(strsplit(d$depths[i], ",", fixed = TRUE)[[1]]))
    v <- v[!is.na(v)]
    if (!length(v)) next
    DBI::dbBind(upd, list(round(stats::median(v), 2), d$ID[i], d$path[i], d$scaffold[i]))
    n <- n + 1L
  }
  invisible(n)
}

#' Classify stored BLAST hit titles and recompute the per-unit screen flag
#'
#' Runs offline over titles already in the database, so an existing project is
#' repaired without repeating its remote BLAST. Units with no `blast_hits` rows
#' are left exactly as they are: absence of hits is not evidence.
#'
#' @param con Open project DB connection.
#' @return (invisibly) the number of `blast_hits` rows classified.
#' @noRd
backfill_hit_class <- function(con) {
  if (!("hit_class" %in% DBI::dbListFields(con, "blast_hits"))) {
    return(invisible(0L))
  }
  todo <- DBI::dbGetQuery(con,
    "SELECT rowid AS rid, title FROM blast_hits WHERE hit_class IS NULL")
  n <- nrow(todo)
  if (n > 0) {
    todo$hit_class <- classify_hit_title(todo$title)
    DBI::dbWithTransaction(con, {
      up <- DBI::dbSendStatement(con, "UPDATE blast_hits SET hit_class = ? WHERE rowid = ?")
      # not on.exit(): inside dbWithTransaction the block is forced in the
      # caller's frame, so on.exit would register on this function rather than
      # the block. finally still releases the handle if dbBind throws.
      tryCatch(DBI::dbBind(up, list(todo$hit_class, todo$rid)),
        finally = DBI::dbClearResult(up))
    })
  }
  # recompute the screen flag for every unit that HAS hits
  units <- DBI::dbGetQuery(con,
    "SELECT DISTINCT ID, path, scaffold FROM blast_hits")
  if (nrow(units) > 0 && "blast_screen" %in% DBI::dbListTables(con)) {
    DBI::dbWithTransaction(con, {
      for (i in seq_len(nrow(units))) {
        h <- DBI::dbGetQuery(con,
          "SELECT hit_class, bitscore FROM blast_hits
           WHERE ID = ? AND path = ? AND scaffold = ?",
          params = list(units$ID[i], units$path[i], units$scaffold[i]))
        DBI::dbExecute(con,
          "UPDATE blast_screen SET is_nuclear_rrna = ?
           WHERE ID = ? AND path = ? AND scaffold = ?",
          params = list(consensus_hit_class(h$hit_class, h$bitscore),
            units$ID[i], units$path[i], units$scaffold[i]))
      }
    })
  }
  invisible(n)
}

#' Copy a project's `.sqlite` into `<path>/.old_sqlite_dbs/` before an upgrade
#'
#' Migrations `ALTER` tables in place, which makes the DB unopenable by an older
#' RiboPilot the user may still have on a cluster; keep a copy to fall back to.
#' @noRd
backup_project_db <- function(path) {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) return(invisible(NULL))
  bdir <- file.path(path, ".old_sqlite_dbs")
  dir.create(bdir, showWarnings = FALSE)
  dest <- file.path(bdir, sprintf("sqlite.%d", as.integer(Sys.time())))
  tryCatch(file.copy(db, dest, overwrite = FALSE),
    error = function(e) warning("Could not back up project DB: ", conditionMessage(e)))
  invisible(dest)
}

#' Migrate the project database at a path to the current schema
#'
#' @param path Project directory containing `.sqlite`.
#' @return (invisibly) the schema version after migration.
#' @export
migrate_project <- function(path = ".") {
  db <- file.path(path, ".sqlite")
  if (!file.exists(db)) {
    stop("No project database found at ", db, call. = FALSE)
  }
  # Back up first when an actual upgrade is pending (older schema version).
  ver <- tryCatch({
    c0 <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(c0), add = TRUE)
    DBI::dbGetQuery(c0, "PRAGMA user_version")[[1]]
  }, error = function(e) 0L)
  if (isTRUE(ver < RIBOPILOT_SCHEMA_VERSION)) backup_project_db(path)

  con <- DBI::dbConnect(RSQLite::SQLite(), db)
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
}
