# Package index

## All functions

- [`RiboPilot()`](https://dmacguigan.github.io/RiboPilot/reference/RiboPilot.md)
  : Open the RiboPilot GUI
- [`annotate_project()`](https://dmacguigan.github.io/RiboPilot/reference/annotate_project.md)
  : Annotate assembled samples for a project
- [`annotate_seq()`](https://dmacguigan.github.io/RiboPilot/reference/annotate_seq.md)
  : Annotate the nuclear rDNA feature model on a single sequence (all in
  R)
- [`apply_project_taxonomy()`](https://dmacguigan.github.io/RiboPilot/reference/apply_project_taxonomy.md)
  : Freeze confirmed taxids + lineages from the review CSV into
  \`samples\`
- [`assign_export_group()`](https://dmacguigan.github.io/RiboPilot/reference/assign_export_group.md)
  : Assign selected samples to an export group
- [`build_reference_db()`](https://dmacguigan.github.io/RiboPilot/reference/build_reference_db.md)
  : Build a clade nuclear-rDNA reference database (GetOrganelle seeds +
  labels)
- [`clear_export_group()`](https://dmacguigan.github.io/RiboPilot/reference/clear_export_group.md)
  : Clear the export group of selected samples
- [`compute_ref_alignment()`](https://dmacguigan.github.io/RiboPilot/reference/compute_ref_alignment.md)
  : Compute per-feature sample-vs-reference alignments for one sample
- [`curate_actinopterygii_rrna()`](https://dmacguigan.github.io/RiboPilot/reference/curate_actinopterygii_rrna.md)
  : Thin actinopterygii curation wrapper.
- [`curate_one()`](https://dmacguigan.github.io/RiboPilot/reference/curate_one.md)
  : Curate one sample's annotations and write results back to the DB
- [`curate_rrna_core()`](https://dmacguigan.github.io/RiboPilot/reference/curate_rrna_core.md)
  : Curate one reconciled nrDNA feature model against a clade ruleset
  (pure)
- [`ensure_curate_schema()`](https://dmacguigan.github.io/RiboPilot/reference/ensure_curate_schema.md)
  : Ensure the curation columns/tables exist (idempotent)
- [`ensure_taxdump()`](https://dmacguigan.github.io/RiboPilot/reference/ensure_taxdump.md)
  : Download a pinned NCBI taxdump tarball into \`dir\` (once per
  project)
- [`export_bundle()`](https://dmacguigan.github.io/RiboPilot/reference/export_bundle.md)
  : Build a GenBank submission bundle for a set of samples (one export
  group)
- [`export_project()`](https://dmacguigan.github.io/RiboPilot/reference/export_project.md)
  : Export every assigned export group in a project
- [`extract_rdna_contigs()`](https://dmacguigan.github.io/RiboPilot/reference/extract_rdna_contigs.md)
  : Keep the rDNA contigs (BLAST hit vs the label DB) and write assembly
  fastas
- [`fallback_assemble_one()`](https://dmacguigan.github.io/RiboPilot/reference/fallback_assemble_one.md)
  : Run an assembly fallback for one sample (Option B2 baiting, or
  Option A)
- [`generate_config()`](https://dmacguigan.github.io/RiboPilot/reference/generate_config.md)
  : Generate and save a reusable cluster config profile
- [`list_configs()`](https://dmacguigan.github.io/RiboPilot/reference/list_configs.md)
  : List available executor configs
- [`load_assemblies()`](https://dmacguigan.github.io/RiboPilot/reference/load_assemblies.md)
  : Load one sample's published assembly fastas into the database
- [`migrate_db()`](https://dmacguigan.github.io/RiboPilot/reference/migrate_db.md)
  : Bring a project database up to the current schema (idempotent,
  additive)
- [`migrate_project()`](https://dmacguigan.github.io/RiboPilot/reference/migrate_project.md)
  : Migrate the project database at a path to the current schema
- [`new_db()`](https://dmacguigan.github.io/RiboPilot/reference/new_db.md)
  : Initialize a new RiboPilot project database
- [`new_project()`](https://dmacguigan.github.io/RiboPilot/reference/new_project.md)
  : Initialize a new RiboPilot project
- [`new_test_project()`](https://dmacguigan.github.io/RiboPilot/reference/new_test_project.md)
  : Initialize a RiboPilot test project
- [`open_project_db()`](https://dmacguigan.github.io/RiboPilot/reference/open_project_db.md)
  : Open a project DB connection with WAL + busy timeout (safe for
  concurrent writes from parallel Nextflow processes)
- [`parse_blast_candidates()`](https://dmacguigan.github.io/RiboPilot/reference/parse_blast_candidates.md)
  : Parse remote blastn hits into ranked reference candidates
- [`parse_feature_table()`](https://dmacguigan.github.io/RiboPilot/reference/parse_feature_table.md)
  : Parse an NCBI feature table (efetch -format ft) into rRNA/ITS
  features
- [`reconcile_one()`](https://dmacguigan.github.io/RiboPilot/reference/reconcile_one.md)
  : Reconcile barrnap + ITSx outputs for one sample; write annotations
- [`resolve_project_taxonomy()`](https://dmacguigan.github.io/RiboPilot/reference/resolve_project_taxonomy.md)
  : Resolve all sample Taxon strings and write a review table for
  confirmation
- [`resolve_taxonomy()`](https://dmacguigan.github.io/RiboPilot/reference/resolve_taxonomy.md)
  : Resolve sample Taxon strings to NCBI taxids + lineages (for
  confirmation)
- [`ribopilot_config_dir()`](https://dmacguigan.github.io/RiboPilot/reference/ribopilot_config_dir.md)
  : RiboPilot user config directory
- [`ribopilot_ref_dir()`](https://dmacguigan.github.io/RiboPilot/reference/ribopilot_ref_dir.md)
  : Default directory for user-built reference databases
- [`run_annotate()`](https://dmacguigan.github.io/RiboPilot/reference/run_annotate.md)
  : Run the annotation workflow (WF2) for a project
- [`run_app()`](https://dmacguigan.github.io/RiboPilot/reference/run_app.md)
  : Run the RiboPilot Shiny application
- [`run_assembly()`](https://dmacguigan.github.io/RiboPilot/reference/run_assembly.md)
  : Run the assembly workflow (WF1) for a project
- [`run_fallback()`](https://dmacguigan.github.io/RiboPilot/reference/run_fallback.md)
  : Run an assembly fallback (Option B2 baiting / Option A de-novo)
- [`score_assembly()`](https://dmacguigan.github.io/RiboPilot/reference/score_assembly.md)
  : Score GetOrganelle candidate paths from real assembly output
- [`score_candidates()`](https://dmacguigan.github.io/RiboPilot/reference/score_candidates.md)
  : Score and rank candidate paths, auto-selecting the best
- [`score_one()`](https://dmacguigan.github.io/RiboPilot/reference/score_one.md)
  : Score one sample's candidate paths and write path_scores (+
  auto_select)
- [`score_project()`](https://dmacguigan.github.io/RiboPilot/reference/score_project.md)
  : Score assembled candidate paths for a project and write
  \`path_scores\`
- [`stitch_one()`](https://dmacguigan.github.io/RiboPilot/reference/stitch_one.md)
  : Stitch one sample's selected path; write the \`stitched\` table +
  return fasta
- [`stitch_scaffolds()`](https://dmacguigan.github.io/RiboPilot/reference/stitch_scaffolds.md)
  : Order + orient a path's scaffolds along a reference and stitch them
- [`write_blast_ref()`](https://dmacguigan.github.io/RiboPilot/reference/write_blast_ref.md)
  : Store one fetched NCBI reference (sequence + features) in the
  project DB
- [`write_coverage()`](https://dmacguigan.github.io/RiboPilot/reference/write_coverage.md)
  : Write per-position read depth for one sample's stitched sequence
- [`write_path_scores()`](https://dmacguigan.github.io/RiboPilot/reference/write_path_scores.md)
  : Write scored candidate paths for a sample into the \`path_scores\`
  table
