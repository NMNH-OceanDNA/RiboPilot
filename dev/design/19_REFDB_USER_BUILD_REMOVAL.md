## 1. WHAT A CLADE DATABASE IS FOR

A clade database is not a BLAST database. It is a pair of plain FASTA files: a GetOrganelle seeds file (`>ACC_Genus_species`) and a labels file (`>18S`/`>5.8S`/`>28S`). The shipped panel is a BLAST v5 index. Different file types, different tools, different pipeline stages. They are never alternatives to each other.

Seeds/labels are consumed in three places, none of which is identification:
- `inst/nextflow/modules/assemble.nf:24` passes them as GetOrganelle `-s` and `--genes`. GetOrganelle runs in `-F anonym` mode (`R/init_db.R:46`), the mode that requires user-supplied seeds. There is no built-in-database GetOrganelle path in this codebase.
- `R/assemble_fallback.R:177,188` uses seeds as the bowtie2 bait index and labels as the BLAST filter that decides which SPAdes contigs are rDNA.
- `R/pipeline_steps.R:21` builds a temporary BLAST db from seeds and feeds it to path scoring. That supplies the taxonomy and identity terms, 0.50 of 1.00 total weight (`R/scaffold_scoring.R:12`), plus the per-scaffold "seed file" hint in the UI.

The panel's only consumer is `inst/nextflow/modules/blast.nf:50,60,109` via `discover_panel()`, where it substitutes for a remote core_nt search. Nothing in the assemble or score chain reads `params.blastdb`; nothing in blast.nf reads seeds or labels.

So: not redundant. The panel cannot backfill a clade database, and a clade database cannot backfill the panel.

Coverage: the clade builder's query is the sample's own taxid with no kingdom clause (`R/reference_counts.R:17-23`), so it covers fungi, plants, protists, anything in GenBank. The panel is `txid33208` metazoa only (`dev/refdb/build_refdb.sh:26-30`) and `R/blast_local.R:5-9` says so. Today a non-metazoan user builds their own seeds in the Reference tab and assembles normally. After removal they have no seeds at all unless they supply their own files. That is the real cost.

## 2. SAFE TO DELETE

Whole files:
- `inst/nextflow/main_reference.nf` (72). Only referrer is `R/reference_nf.R:152`. Holds `params.clades/refDir/outDir/refseq`, processes `count_clade`/`build_ref_db`, workflows COUNT/BUILD. Nothing in `main.nf:26-32` includes it. No `withName` in nextflow.config names either process.
- `inst/nextflow/bin/rp_ref_db.R` (35). Only callers `main_reference.nf:33,56`. Kills the RP_REFSEQ channel with it.
- `inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh` (165) and `parse_rrna_fasta.py` (138). Callers are `R/build_reference_db.R:57` and `rp_ref_db.R:18` only. Keep the sibling `bin/curl`, see section 3.
- `R/reference_groups.R` (240) and `R/reference_counts.R` (60). Sole consumers are the Propose-groups and ladder observers in `R/app_reference.R`.
- Tests: `test-reference-modal.R`, `test-reference-groups.R`, `test-reference-nf.R`, `test-reference-workflow.R`, `test-reference-counts.R`.
- `vignettes/Reference-Databases.Rmd`, minus the "reference folder" section at :96-114 which must be relocated.

By line range, never by file:
- `R/build_reference_db.R`: delete `build_reference_db()` (:53-122). Move `ribopilot_ref_dir()` (:13-17) elsewhere first.
- `R/reference_nf.R`: delete :1-56 and :146-295. Keep :57-145.
- `R/reference_dir.R`: delete the manifest block :102-150.
- `R/reference_build.R`: delete everything except `usable_fasta()` :329-335, `reference_db_ok()` :343-346, `missing_assembly_dbs()` :357-369, which move out.
- `R/app_reference.R`: delete ~740 of 802. Keep `reference_settings_modal()` :139-177, its opener :779-785, its apply observer :786-795.
- `R/app_modals.R`: `reference_run_modal()` :395-447. `R/app_run.R`: `build_reference_submit_script()` :130-163.
- `R/app_server.R`: :101 mode_meta, :118 table call, :264-285 colDefs, :296-298 nodata, :488-495 build poll, :666-675 pending warning. Keep :653-665.
- `R/app_ui.R`: :51 mode entry, :60-62 control bar.
- Schema: `reference` table (`R/db_schema.R:120-138`), `samples.ref_id`, `reference_config.min_records/max_records`. Every reader dies in the same change, so no migration is needed; old project DBs keep the columns inert.
- `NAMESPACE:8` (`build_reference_db`), `R/constants.R:22` stale link, `R/utils.R:98-101` error hint.
- Test edits, not deletes: `test-reference-dir.R` (3 of 15 tests), `test-reference-app.R`, `test-reference-schema.R:22`, `test-write-atomicity.R:80-107`, `test-nextflow.R:116,346,376,479`, `test-nf-bin-contract.R:68`, `test-nf-command-quoting.R:111,121`, `helper-app-fixture.R:12,46`, `test-app-opts.R:193-210`.

## 3. MUST SURVIVE

- **Cross-project reference store** (`R/ref_store.R`). Sits under the same folder as clade artifacts but shares no file, no reader, no function. Minimal change: keep the `reference_config.ref_dir` column, `nf_ref_dir()` (`R/nf_command.R:71-79`), the ref_dir branch of `reference_mount_dirs()` (`R/reference_nf.R:108-116`), and `new_project(ref_dir=)`. Everything fails closed (`R/ref_store.R:215,265`), so losing ref_dir turns the cache off silently, with no error, on the very branch that exists to build it.
- **Per-run GenBank fetch and align**. `blast_ref_fetch` -> `rp_ref_store.R` -> `rp_ref_align.R` -> `curate`. Reads no clade artifact anywhere. No change needed. Note this chain still hits NCBI live on every run; the shipped panel removes the search round trip, not the sequence fetch, so runs are not offline either way.
- **`ribopilot_ref_dir()`**, currently inside `R/build_reference_db.R:13-17`. It is the fallback for `ribopilot_db_dir()` (`R/db_dist.R:76`), which is the default arg of `discover_panel()`. Docker sets `RIBOPILOT_DB_DIR` so the image is safe, but a conda or local install falls through to it. Delete the file whole and every non-Docker run silently drops to remote BLAST. Move the 5 lines into `R/db_dist.R`.
- **`R/reference_nf.R:57-145`**: `container_bind_opts()`, `config_mounted_dirs()`, `reference_mount_dirs()`, `drop_covered()`. `R/nf_command.R:57-58` calls them on every pipeline run to mount the store, the panel folder, and local seed FASTAs.
- **`usable_fasta()`, `reference_db_ok()`, `missing_assembly_dbs()`** from `R/reference_build.R`. Consumers are `R/app_server.R:657` (the Run-button gate) and `R/app_opts.R:457`. That gate is wrapped in tryCatch, so deleting them does not crash: it degrades to "nothing dead" and the user submits a cluster job that exits 0 with zero assemblies. Keep `test-reference-build.R:277-300` with them; it is the only coverage.
- **`reference_config.blastdb` / `blastdb_version`**, `validate_blastdb_prefix()`, `blastdb_provenance()`. Panel concerns, not clade concerns. The settings modal is the only writer of `blastdb` anywhere; there is no non-GUI setter.
- **`assemble_opts.seeds_db/labels_db`**, `live_db()`, the `ao_seeds`/`ao_labels` fields (`R/app_opts.R:138-141`), and `new_project(custom_seeds_db=, custom_labels_db=)`. These become the only remaining route to seeds.
- **`inst/ref_dbs/getOrganelle/GenBankDownload/bin/curl`**. `dev/refdb/build_refdb.sh:46-47` puts it on PATH and :256-257 hard-fails if curl does not resolve there. Deleting the parent tree breaks the maintainer's own panel rebuild. Edit that directory's README, do not delete it.
- `ncbi_api_key` in all ten config templates and the `<<EXTRA_MOUNTS>>` token. Both serve the main pipeline.

## 4. THE UI AFTER REMOVAL

The Reference tab loses its table, Propose groups, the clade ladder, the row-action modal, the run modal, the build poller, and the log tail. What is left is one gear button opening a modal whose three sections split as: reference folder (keep, it is the store), record window (delete), BLAST panel prefix plus "Last run searched against" (keep, panel only).

That is not a tab. Kill the Reference mode from `R/app_ui.R:51` and rehouse the gear as a mode-independent Settings button. Drop the record-window section. Also reword `R/app_server.R:662`, which currently tells users to "Rebuild them in the Reference tab".

Worth flagging: the provenance line only fills in when the user explicitly set `reference_config.blastdb` (`R/nf_command.R:105-106` returns early otherwise), so in the shipped-panel-only world this change is heading toward, it is already blank.

## 5. WHAT THE USER LOSES

- **Any non-fish user loses the only in-product way to get assembly seeds.** The default is a hardcoded Actinopterygii pair fetched from raw.githubusercontent.com (`R/init_db.R:43-44`). The image ships no seeds or labels; `docker/Dockerfile:100-102` adds only the panel and the decoys. `dev/refdb/build_refdb.sh` has no FASTA export stage. Nothing replaces this today.
- **The degradation is silent, not an error.** `main.nf:79-86` routes a seedless sample to the bait fallback and rewrites its method; the run exits 0. With labels blank, `R/assemble_fallback.R:118-120` skips the rDNA filter entirely and keeps every contig over 300 bp, which either trips the max_scaffolds gate or pushes junk into annotation.
- **Fish seeds on a sea star is a wrong-result path, not just weaker baiting.** Path auto-selection scores against fish rDNA for 50% of its weight, and every scaffold card in the app shows a fish species as its closest hit.
- Non-metazoan users (fungi, plants, protists, diatoms) go from "build your own clade database in the app" to "produce your own FASTAs or accept fish".
- The guided "use a database I already have" flow goes. The Assemble options text fields do the same job, minus the hard FASTA validation and minus the manifest list of previously built databases.
- Per-clade GenBank record counts and automatic clade grouping go entirely.

No workflow becomes strictly impossible, provided the user can produce seeds and labels themselves. Everything above is a downgrade, and most of it is quiet.

## 6. SIZE

Roughly 45 to 50 files, about 4,900 lines. R sources ~2,150 (the largest bucket), tests ~1,970, shipped scripts ~486, docs ~250 hand-written plus regenerated man pages and pkgdown output.

Five files must be cut by line range, never deleted: `R/build_reference_db.R`, `R/reference_nf.R`, `R/reference_dir.R`, `R/reference_build.R`, `R/app_reference.R`. Plus one directory carve-out, `inst/ref_dbs/getOrganelle/GenBankDownload/bin/`.

Fix the test fixtures in the same commit. `helper-app-fixture.R:12,46`, `test-nextflow.R:346`, and `test-nf-command-quoting.R:111` all INSERT the dropped columns; miss them and a third of the suite fails at helper load, which will look far worse than it is. Three tests hard-error on a missing `main_reference.nf` (`test-nextflow.R:376,479`, `test-nf-bin-contract.R:69`) and one file aborts at source time (`test-reference-workflow.R:1-3`).

`_pkgdown.yml:42` must lose the article or the site build fails. Do not touch `_pkgdown.yml:27`; that is the function reference index, unrelated.

## 7. OPEN QUESTIONS FOR THE MAINTAINER

1. **Will the image ship GetOrganelle seeds and labels alongside the panel?** This is the only real blocker. `dev/refdb` builds none today, and it needs a new export stage. If the answer is no, then GetOrganelle becomes user-supplied-FASTA-only, `main.nf` should say so loudly instead of silently falling back, and `README.Rmd:54` must be rewritten since it currently promises RiboPilot "can build the reference data it needs automatically".
2. **Does the "use a database I already have" modal survive in some form?** The Assemble options fields cover the mechanics, but a user with lab-shared FASTAs loses the guided path and the hard validation.
3. **Where do the surviving settings live?** The reference folder field and the BLAST panel prefix field both need a home once the Reference mode is gone, and `ref_dir` in particular currently has no other in-app setter.