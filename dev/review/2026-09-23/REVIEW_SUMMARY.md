# RiboPilot review summary

Date: 2026-09-23. Code reviewed: `main` at 3362308 (version 0.0.0.9018), read on branch `review/2026-09-full-review`.

This page pulls together the three track reports in this folder:

- [CODE_REVIEW.md](CODE_REVIEW.md): the R package, Shiny app, Nextflow pipeline, containers, tests, and packaging.
- [DOCS_REVIEW.md](DOCS_REVIEW.md): the unpublished pkgdown site (README, three vignettes, function reference), `docker/README.md`, and in-app help.
- [USABILITY_REVIEW.md](USABILITY_REVIEW.md): hands-on use of the app on real projects.

## Verdict

RiboPilot works end to end today, and a lot of it works well: two operators took real projects from raw reads to a GenBank bundle using only the app's buttons, a third swept common user mistakes, the off-target screen caught both real problems in the test data, the test suite passes (3,952 checks, 0 failures), and reviewers praised the evidence views, the Reference tab, and the Export window. It is not ready for outside users or for publishing the docs yet, for three reasons. First, it can write wrong GenBank records with no warning: a sea anemone contaminant went out under a snail's name behind a green all-clear, truncated 18S and 28S pieces are labelled "complete sequence", and assembling a sample a second time can leave old gene coordinates, verdicts, and scaffolds attached to the new sequence. Second, a new user cannot get started: the install command, the first project command outside RStudio, the default executor, and (while the repo is private) the built-in fish reference all fail, and the Get started guide describes an older, console-driven, three-tab app. Third, the app often says things went well when they did not: failed runs end with a green "Done", and closing the browser tab kills a run. None of this needs a redesign. Most fixes are small and local, and a handful of focused changes (an export gate, partial markers from gene coverage, clearing a unit's evidence when its sequence changes, four small setup fixes, and a rewritten Get started) remove the worst of it. No code finding is critical; 20 are high.

## Fixed on this branch after the review

Five batches of fixes were made on `review/2026-09-full-review` after the review was written. The first batch (2026-09-23 and 2026-09-24) covered actions 1 to 4, in the table below. The second batch (2026-09-24) covered actions 6 to 15, name and path validation, new functions for adding samples, code simplification, and the docs rewrite; it is described under [Second batch](#second-batch-actions-6-to-15). Action 5 was deferred by the author. A third batch (2026-09-26) fixed the re-run, caching, BLAST, and visibility items; see [Third batch](#third-batch-re-runs-caching-blast-and-visibility). A fourth batch (2026-09-26) fixed working at scale, reproducibility and packaging, the docs items, and the Hydra memory setting; see [Fourth batch](#fourth-batch-working-at-scale-reproducibility-and-docs). A fifth batch (2026-09-28) added per-gene partial ends, batch Trim and Restore, the test and CI gaps, and the remaining docs wording; see [Fifth batch](#fifth-batch-partial-ends-batch-trim-tests-and-docs). The full test suite passes (1,302 tests, 0 failures).

### First batch: actions 1 to 4

| Problem | Items | What changed | Files | Tests |
|---|---|---|---|---|
| Contaminated records exported under the sample's name | UX-31, action 1 | `submission_risk` now carries `off_target(<taxon>)` from the contamination screen and `contaminant` from curation, so the summary CSV, HTML report, done dialog, console warning, and bundle README all show it. The Export window now lists every off_target or contaminant unit not marked Reviewed, and offers "Leave them out" (default), "Export anyway", or "Cancel". | `R/export.R`, `R/app_modals.R`, `R/app_server.R` | `test-export.R`, `test-app-export.R` |
| Truncated genes exported as "complete sequence" | DB-H2, UX-32, action 2 | A feature the curation step flagged `truncated` is now partial at both ends in the `.tbl` (`<start >stop`) and in the definition line ("partial sequence"). | `R/export.R`, `R/export_opts.R` | `test-export.R` |
| Old gene positions and screen results kept after re-assembly | DB-H1, SCI-H1 (sequence part), action 3 | When a re-assembly gives a surviving unit a different sequence, `load_assemblies()` deletes that unit's annotations, coverage, screen verdict, chimera call, BLAST results, reference alignments, and trim backup. It also clears its lock, Reviewed, Partial, Problematic, and export group, then re-queues it for annotation. Units whose sequence did not change keep everything. | `R/pipeline_steps.R` | `test-load-assemblies.R` |
| First commands fail on a new machine | PKG-H1, PKG-H2, PIPE-M11, PIPE-H2, UX-01, UX-04, action 4 | DESCRIPTION gains an empty `biocViews:` line, so `install_github` finds the Bioconductor packages. The RStudio step runs only inside RStudio, so `new_project()` and `new_test_project()` work from a terminal, Rscript, or HPC session. Both default to `executor = "docker"`, and conda warns. The local, docker, and conda configs cap every task at this machine's CPUs and memory (`resourceLimits`, filled by `new_project()`), so an oversized request no longer ends the run. The README and vignette lines saying conda is the default are updated. | `DESCRIPTION`, `R/init_project.R`, `R/init_test_project.R`, `R/generate_config.R`, `inst/config.{docker,local,conda}`, `man/`, README, vignettes | `test-init-reads.R` |
| The docs described an older, console-driven app, and key pages were missing | DOC-28, UX-49, DOC-20 and UX-09 (docs side), DOC-37, DOC-38, DOC-40, UX-37, DOC-30, DOC-35 (docs half), DOC-36, DOC-41, DOC-42 (docs half), DOC-43, DOC-44, DOC-45, DOC-46, DOC-47, DOC-49, DOC-56, DOC-57, DOC-58, DOC-59, DOC-61, DOC-01 to DOC-06, DOC-08 to DOC-11, DOC-13, DOC-16, DOC-21 to DOC-24, DOC-27, DOC-52, DOC-55 (rule 8), DOC-64, DOC-66, DOC-68, DOC-90, DOC-94; action 6 | Get started now walks the practice project through the app's five steps with 18 screenshots from the end-to-end check, and new pages cover installation, your own project, GenBank submission, scripting, clusters, caches and offline use, the unsupported conda route, and troubleshooting. The Reference and Curation guides are updated (run lifecycle, rule 8, chimera levels, depth, variants), the README is ASCII, and the function reference is grouped with internals hidden. | `vignettes/*.Rmd`, `vignettes/figures/`, `_pkgdown.yml`, `README.Rmd`, `README.md`, `R/init_test_project.R` and `R/run_assembly.R` (roxygen only), `man/new_test_project.Rd`, `man/run_assembly.Rd` | `test-docs-site.R` |

Checked on real data: Tret_01 now reads `off_target(Actiniaria);contaminant`, and the Export window asks about it ([screenshot](screenshots/fix_01_contamination_prompt.png)). "Leave them out" exported only the five *Euspira* records. SRR33849683_p1_s2 now exports as `<155 >1488` with "28S ribosomal RNA gene, partial sequence."

### Second batch: actions 6 to 15

Specs were written per work package, built in parallel in separate git worktrees, combined, simplified, checked end to end in the app, reviewed once, and given one final fix pass. Author decisions taken along the way are marked "(author)".

| Action | Items | What changed | Main files |
|---|---|---|---|
| 6. Get started follows the app | DOC-28, UX-49, DOC-20, UX-09, PIPE-H5, DOC-37, UX-16, DOC-40, UX-37, DOC-30, DOC-35, DOC-38 | The pipeline and the Update count only take samples whose taxonomy is locked. `run_fallback()` is deleted; the rescue is the Assembler option in Assembly Options (author). `export_project()` exports every exportable unit as one "export" bundle when no groups exist, and uses the same exportable-unit rule as the app (both locks). The docs rewrite is in the docs row of the first-batch table. | `inst/nextflow/modules/pipeline_r.nf`, `R/nf_command.R`, `R/export.R`, `R/run_fallback.R` (deleted) |
| 7. Keep the taxon right | REF-H1, REF-H2, REF-H3, UX-08, REF-M1, APP-L16, SCI-H1 (taxon half), DOC-24 | Changing a taxon writes `lineage_taxids` in the same update and re-runs the screen for that sample. Homonyms are never picked automatically: scientific-name matches win, and the alternatives are stored with their kingdom (`samples.tax_alternatives`) and shown in the table and search. The Assign popup counts samples by confidence. Locking a sample with no taxon asks first (author), and a taxonomy-locked sample refuses Assign on the server. Propose groups names the samples it left out. The taxdump is downloaded to a temporary file with a long timeout and NCBI's md5, the old index is deleted, a failed untar stops, and an index with too few nodes is refused. Clade record counts and builds query by taxid. | `R/taxonomy.R`, `R/tax_index.R`, `R/tax_search.R`, `R/reference_groups.R`, `R/reference_build.R`, `R/app_server.R`, `R/init_db.R` |
| 8. Tell the truth about runs | APP-M2, UX-42, DOC-88, UX-15, DOC-08, UX-45, APP-M3, PIPE-M6, PIPE-M9 | A shared run window (`R/nf_run.R`, `R/app_run.R`) reads the exit status and the Nextflow log. The end-of-run alert says "Run finished with no failures", or names what failed. Each failed task gets a line with the RAW exit code (for example 140), sample ID, pipeline step, likely cause, and attempts (author). The report has a one-click copy block with project, versions, run name, and log path, for sending to the developers (author), and is saved as `nextflow_logs/<run>_report.txt`. Update with nothing queued says "Nothing to update". A missing Nextflow, or an error while building the report, no longer ends the session. Each failed sample gets a plain note (gate limit, below minDepth with its read count, SPAdes exit code and error line), shown as a red banner in the details window with the next step. | `R/nf_run.R` (new), `R/app_run.R`, `R/app_server.R`, `R/app_reference.R`, `R/app_modals.R`, `R/assemble_fallback.R`, `inst/nextflow/main.nf`, `inst/nextflow/modules/preprocess.nf`, `inst/nextflow/modules/fallback.nf` |
| 9. App bugs that end the session or overwrite work | APP-H1, APP-H2, APP-M4, APP-M8, UX-18, APP-M14, DB-M4 | The options cells work on samples with several units. Picking an existing option set reloads its values and unticks Edit; a new name creates the set from the values shown (author). Curation writes its flags to a new `annotate.curate_flags` column, so the curator's Notes are never overwritten. The app opens the project database with a busy timeout everywhere, and a busy database gives a message instead of ending the session. A second trim reads the unit's current span, and the range is checked at Apply. | `R/app_opts.R`, `R/app_server.R`, `R/curate_rrna.R`, `R/init_db.R`, `R/export.R`, `R/utils.R` |
| 10. Scaffold and path bookkeeping | SCI-H3, SCI-M8, SCI-H2, SCI-M6, PIPE-H4, PIPE-M5 | Path numbers come from the file name, not the file listing order. Scaffold ids are read as text, so 1.10 is not 1.1. Circular units are rotated to start in the middle of the IGS before they are stored (author: no app Rotate button for now). The 18S and 28S anchors are the copies nearest the ITS region, on the majority strand. Supplied assemblies with more than 5,000 contigs no longer collide on chunk file names. Max paths 0 means no limit, and GetOrganelle exports one path more than the limit so the gate can see overflow (author: default stays 10). | `R/select_paths.R`, `R/pipeline_steps.R`, `R/annotate.R`, `R/coverage.R`, `inst/nextflow/modules/find_rrna.nf`, `inst/nextflow/modules/assemble.nf` |
| 11. The rest of each GenBank record | DB-M1, DOC-45, UX-36, DOC-44, DB-M5, UX-84, DOC-39, DB-L11, DOC-49, DB-L14 | `[organism=]` uses the name confirmed in the Taxonomy step for exact matches and manual species picks, and a `taxon_differs(<name>)` flag marks records whose Taxon label differs. The made-up isolates (`ID_p1_s1` ...) stay, because they let several scaffolds of one sample go into one GenBank submission (author); a template may carry its own `[isolate=]` unless the group has such a sample. ITS spacers are snapped to the rRNA genes when the gap is 5 bp or less, and larger gaps get a "gap" curation flag. Assembly gaps use "unspecified" linkage evidence. Empty source modifiers give a warning. Unannotated ETS and IGS flanks are kept, because GenBank does its own annotation (author). | `R/export_opts.R`, `R/export.R`, `R/annotate.R`, `R/curate_rrna.R` |
| 12. The cluster route | PIPE-H3, DOC-29, PIPE-M14, APP-M10, UX-54, PIPE-L2, PIPE-M10, DOC-71, DOC-69, UX-79 | Generated container profiles no longer mount `rawDir`, so assembly-only projects parse; the reads folder goes in the project's own mounts. Head-job scripts get a header for SLURM, SGE, PBS, or LSF from the profile (queue, account, walltime), with a single `%` in SLURM log names, and LSF jobs are submitted on stdin so `#BSUB` lines are read. Work Dirs finds cluster runs through `nextflow_logs/`. Singularity profiles get a pull timeout and a shared image cache under the reference folder. awsbatch is removed until it is supported. Cluster mode (job script instead of Run from App) now depends only on a scheduler on PATH or RStudio Server, not on `launch.browser = FALSE` (author). | `R/generate_config.R`, `R/app_run.R`, `R/hydra.R`, `R/app_utils.R`, `R/reference_nf.R`, `inst/config.*` |
| 13. Stop runs cleanly | APP-M1, UX-43, UX-44, DOC-36, APP-M9, UX-46, APP-L11, REF-L17, PIPE-M2 | Runs do NOT survive the browser (author): pressing Stop or closing the tab stops the run gracefully. The app sends SIGINT, then SIGTERM after 10 s (for an app started with SIGINT ignored, for example with nohup), then kills the process tree after 60 s. This project's leftover task containers are removed. One project-level run lock (`.ribopilot_run.lock`) covers pipeline runs, reference builds, and reference counts across browser windows and R sessions (author), and is held until the head process and this project's task containers are gone. Reference runs launch from their own `.reference` folder, so `-resume` of the pipeline is not disturbed. | `R/nf_run.R`, `R/app_run.R`, `R/app_server.R`, `R/app_reference.R`, `R/reference_nf.R`, `inst/app/www/custom.js`, `inst/nextflow/nextflow.config` |
| 14. Locks lock, marks say what they do (MitoPilot model) | APP-M7, UX-20, DOC-90, APP-M15, UX-25, UX-33, DOC-41, APP-L10, UX-21, UX-23, TEST-L4 | A locked unit is read only: the window shows a lock banner and disables Edit, Trim, Restore, and the marks, and the server refuses them too. Export needs both the assembly lock and the annotation lock. A second assembly path is refused everywhere, as in MitoPilot (author); the Lock anyway prompt is gone. Lock & Close only locks; Reviewed stays the curator's own sign-off (author). Problematic is renamed Exclude from export, and the Export window names the units it leaves out and why. Toolbar toggles say what they will do and skip locked units. The options guard names the right lock and tab. | `R/app_server.R`, `R/app_modals.R`, `R/app_opts.R`, `R/export.R`, `R/constants.R` |
| 15. Chimera levels | SCI-M1, SCI-M2, UX-83, DOC-61, SCI-M3, SCI-M4, UX-82, DOC-57, DOC-56, APP-M11 | Labels are now clean, minor, and severe (author). Severe is a lineage switch at phylum or kingdom rank with at least 90 percent identity on both sides (for example CHIM01, fish plus *Penicillium*); it goes to `submission_risk` as chimeric and to the Export prompt. Minor is a weaker disagreement, including an ancestor pair where one side resolves only to phylum or above (author: the Chordata against Pleuronectiformes case); it is a warning in the app and report only. Other ancestor pairs count as agreement, ties at the top bitscore use their common ancestor, soft clips come from a local alignment, depth steps skip scaffold ends and need a minimum depth, and junctions move with a trim. ITSx's own call stays as `chimeric(ITSx)`. The chimera panel shows the unit's level and reason. After a user report (2026-09-25), a trim that removes every junction sets the label to clean, and Restore brings back the chimera result and windows saved with the trim backup (new `assembly_backup.chimera` column). | `R/chimera.R`, `inst/nextflow/modules/coverage.nf`, `R/export.R`, `R/app_modals.R`, `R/pipeline_steps.R` |
| Names and paths | DB-M2, DB-M3, PIPE-M13, UX-35, SEC-M1 | One set of checks runs at `new_project()`, `add_samples()`, and `update_sample_metadata()`: required columns in any case, trimmed values, duplicate IDs, IDs unsafe for file names or shells (with a suggested safe ID), missing read or assembly files (listed), R1 equal to R2, BOM, CRLF, semicolon CSVs, xlsx given as CSV, empty files, and mixed gzipped and plain assemblies. Export group names must be folder-safe, with a suggestion. Header templates are filled by a plain column lookup and are never run as R code. | `R/init_db.R`, `R/init_project.R`, `R/export_opts.R`, `R/app_server.R` |
| New: add samples and metadata | (author request) | `add_samples()` and `update_sample_metadata()`, ported from MitoPilot. New samples are seeded exactly like `new_db()`, must use the project's own read and assembly folders, and duplicate IDs are refused. Metadata updates never touch pipeline state or locks, except that a changed Taxon resets that sample's taxonomy step. | `R/add_samples.R`, `R/update_sample_metadata.R`, `NAMESPACE` |
| Simplification | P2, P5 to P9, P12 (part), P13, P14 (part), P21, P28, P30a, P35 | About 775 fewer lines of R. Removed: dead files (`resolve_project_taxonomy.R`, `score_project.R`, `rp_ref_fetch.R`), dead helpers and JavaScript handlers, the stringr, crayon, and golem dependencies, dead tables and columns (for example `assembly_blast`, `assemblies.depth/gc/errors/length_raw`, `samples.reference_support`), and the unused Curation options. `triage_opts` is seeded from one source. | many |

How it was checked:
- **Tests:** 1,227 tests, 0 failures, on a PATH without conda.
- **End to end:** a fresh practice project at a neutral path, run from the app in headless Chrome on a local image built from the combined code (`ribopilot:review-fixes`). 25 behaviors were confirmed. Examples: Stop and tab close stop the run and release the lock after cleanup; a second window cannot start a run; Resume after Stop reused cached tasks; the failed sample shows a banner with the next step; CHIM01 is severe and one fish unit minor; the Export window asks about CHIM01 and names left-out units; SRR33849683_p1_s2 exports as `<155 >1488`; locked units are read only.
- **Review pass:** 5 reviewers found 24 issues. The one high finding (the header check still ran template code through glue) was confirmed and fixed. A final pass fixed 18 regressions in the new code, for example: a run failing before any task ends the session; export_project() ignoring the assembly lock; case and Taxon handling in update_sample_metadata(); the folder check in add_samples(); severe chimera calls depending on pair order; pkgdown topics; LSF submission. The end-to-end gaps were fixed in the same passes: Stop when SIGINT is ignored, a report when the tab closes, the taxonomy lock, the README wording, the chimera panel, the BLAST note for failed samples, the empty NCBI key warning, and the reference build's failure reason.

### Third batch: re-runs, caching, BLAST, and visibility

Built on 2026-09-26 the same way: specs, four packages in parallel worktrees, combined, an end-to-end run of the practice project from the app on a local image built from the combined code, one review pass with skeptic checks, and one bounded fix pass. The package version was deliberately left at 0.0.0.9018 (author).

| Area | Items | What changed | Main files |
|---|---|---|---|
| Pipeline caching | PIPE-H1, PIPE-M4, PIPE-M7, PIPE-M8 | Scoring and loading read the assembly files and graph staged into their own task, never the publish folder, so a re-run cannot pick up old path files. The fallback's seeds, labels, bowtie2, and SPAdes options are task inputs, so `-resume` re-runs a fallback when they change; local databases are staged into separate folders and URLs are passed as values, so one bad URL fails one sample, not the run. `find_rrna_pick` and `stage_assembly` always run their database writes. An incomplete rRNA scan of a supplied assembly fails the sample with "rRNA search incomplete: N of M part(s) ... could not be scanned" (author: fail, do not keep a partial result). | `inst/nextflow/main.nf`, `modules/*.nf`, `bin/rp_score.R`, `bin/rp_fallback.R`, `bin/rp_find_rrna.R`, `R/pipeline_steps.R` |
| Re-run state | SCI-M5, UX-17, PIPE-M3 | The curator's own triage decisions (source user) follow their sequence through a re-assembly, including on the gate path, and identical copies keep at least one kept copy (author: only the curator's decisions are carried; a curator's keep beats the minimum length). A rescue through the Assembler option annotates its own new scaffolds in the same run. Changing annotation options only re-queues annotation: no re-assembly, and trims stay. The gate path now also clears the evidence of scaffolds whose sequence changed. | `R/pipeline_steps.R`, `R/app_opts.R` |
| BLAST and science | PIPE-M1, REF-M5, SCI-M7, PIPE-M16, SCI-M9 | A unit is searched again when its BLAST settings change (the settings are part of the search key), and a new search clears the old screen verdict so it is recomputed. PDB reference hits are stored under the key BLAST uses, so the top reference is aligned. The isomer and paralog classes are computed after annotation finishes. The bait fallback refuses to run without a usable seed or label database instead of assembling every read (author). An unreadable seed database no longer stops a sample: its paths are scored without it and a note says so (author). A scaffold put back by the curator after a re-assembly is now BLASTed and screened. | `inst/nextflow/modules/blast.nf`, `R/blast_genbank.R`, `bin/rp_ref_store.R`, `R/screen.R`, `R/screen_variants.R`, `R/assemble_fallback.R`, `R/select_paths.R` |
| Seeing problems | UX-14, UX-22, UX-34, action 1 | A Checks column right after the ID on Assemble and Annotation Review sums up each sample or unit with colour and text: red for off target, contaminant, severe chimera, or ITSx chimeric; amber for minor chimera, an uncertain screen, ITSx or curation warnings, and feature flags; green "ok". Every Lock (Assemble toolbar, Annotation Review toolbar, and the unit window's Lock and Lock & Close) lists flagged units not marked Reviewed and asks before locking (author: both tabs ask). The Export complete dialog names the flagged records in a bold amber line, and the ITSx line is never green next to a risk. The old "Checks" column group is now "Check details" (author). | `R/export.R`, `R/utils_reactable.R`, `R/app_server.R`, `R/app_modals.R`, `R/constants.R`, vignettes |

How it was checked:
- **Tests:** 1,264 tests, 0 failures (`NOT_CRAN=true`, so the tests that run the real bin scripts also ran).
- **End to end:** a fresh practice project at a neutral path, run from the app on a local image built from the combined code. Confirmed: the Checks column on both tabs (CHIM01 red, SRR33849683 amber, clean fish green); the Lock prompt on all four Lock paths, with Cancel locking nothing; a curator's drop surviving a re-assembly; a trim surviving an annotation-option change, with only annotation re-run; SRR22396893 rescued with baited SPAdes in one Update; a re-search after a BLAST setting change, and no re-search without one; the Export complete dialog naming CHIM01; no process, container, or run lock left after any run.
- **Review pass:** 3 reviewers found 15 issues. The 2 high ones were confirmed and fixed: a re-search left the old screen verdict in place, and the gate path kept stale evidence on changed scaffolds. The batch's own medium regressions were fixed too: the fallback URL and file-name staging, a curator's keep on identical isomers, and put-back scaffolds never searched. The e2e mismatch between the Checks "warnings" and the Warnings column was fixed by showing feature flags separately.

### Fourth batch: working at scale, reproducibility, and docs

Built on 2026-09-26 the same way: specs, three code packages in parallel worktrees (scale, reproducibility, testing), combined, a local image built from the new pins, then the docs package on the combined code, one review pass with skeptic checks, one bounded fix pass, and a final end-to-end run with the screenshots retaken. The version stays 0.0.0.9018 (author).

| Area | Items | What changed | Main files |
|---|---|---|---|
| Working at scale | APP-M5, APP-M6, APP-M12, UX-57, UX-58, B3-app-2 | The Lock and State filters now take rows out of the table's data, so a filtered table pages and counts only what passes (no more blank pages). A refresh keeps the page and the ticked rows, matched by sample or unit. The browser reports which rows the search box and column filters leave, and bulk actions (Lock, State, marks, Triage) skip hidden rows; the header tick box ticks only rows shown. The busy loop in `custom.js` is gone: an idle table at 504 samples now uses under 1 percent of the browser's time (was 10 percent at 50 rows per page and about 50 percent at All). Triage takes the ticked, unlocked samples when its window opens and skips assembly-locked ones with a note. Lock anyway after a background refresh now locks what was ticked. | `R/app_server.R`, `R/app_ui.R`, `R/app_ui_helpers.R`, `inst/app/www/custom.js` |
| Reproducibility | PIPE-M12, PIPE-M15, PIPE-M17 | Every tool in the image is pinned to the version actually installed, `environment.yml` carries the same pins, the ribovore installer is checked against a sha256, and seven unused R packages are dropped. When the app opens a project whose `.config` names the RiboPilot image for another version, it says so ("Container version does not match RiboPilot") and how to switch, then carries on; custom images are not compared. This follows MitoPilot (author, 2026-09-27: warn, do not block); a first version that stopped every pipeline R step on a mismatch was replaced. The batch-3 guard in `rp_score.R` is gone. Which database answered each unit's BLAST search (the bundled rDNA database and its release, or NCBI core_nt) is recorded and shown in the unit window. | `docker/Dockerfile`, `environment.yml`, `R/db_dist.R`, `R/app_server.R`, `inst/nextflow/bin/rp_searched.R`, `R/db_schema.R`, `R/app_modals.R` |
| Hydra memory | B3-pipeline-4 | The Hydra profile clears the memory setting for every step, so no job asks the grid engine for `h_rss` or `mem_free`, which Hydra rejects. The rRNA search steps get their own memory request (the per-sample memory for the scan, 4 GB for the pick and staging steps). Checked with a stand-in grid engine; no real Hydra run was tried. | `inst/config.NMNH_Hydra` |
| Testing and packaging | TEST-M1, action 4 (CI), DOC-95, PKG-M1 | New tests run the pipeline's R scripts in order and the real `main.nf` end to end on the practice samples, with stand-in tools for fastp, GetOrganelle, barrnap, ITSx, BLAST, and the rest. A new GitHub Actions job installs RiboPilot on a clean runner the way the docs say and runs them. The site URL is in DESCRIPTION. Roxygen markdown is on, and every help page is regenerated with no warnings. | `tests/testthat/helper-fake-tools.R`, `test-pipeline-drivers.R`, `test-nf-main-run.R`, `.github/workflows/clean-install.yaml`, `DESCRIPTION`, `man/` |
| Docs | DOC-07, DOC-12, DOC-15, DOC-53, DOC-60, B3-app-3 | Installation has an Upgrading section (recreate projects when the database changes; match the container image; pull again or build it yourself). NEWS.md, `inst/CITATION`, and a Citing page with a tool-version table that matches the pins (13 references checked against Crossref). The docs and the app use the same names ("step", "BLAST screen", "bundled rDNA database", "BLAST evidence", Hit Source "bundled", "NCBI", "seed file"), and the evidence line names the database that answered. The warnings table lists every code and says what happens when a lookup is missing. The Checks sentence is fixed and 8 screenshots are retaken. | `vignettes/*.Rmd`, `NEWS.md`, `inst/CITATION`, `_pkgdown.yml`, `README.Rmd`, in-app strings in `R/` |

How it was checked:
- **Tests:** 1,281 tests, 0 failures (`NOT_CRAN=true`, so the new pipeline tests and the stand-in grid engine test ran).
- **Image:** a local image built from the new pins; every pinned package resolved to its pin.
- **Scale:** a 504-sample copy of the practice project driven in headless Chrome: 17 checks, all as expected (idle cost, paging under State = Failed, refresh keeping page 3 and 2 ticked rows, search and column filters dropping hidden rows from the selection, the header tick box, the details window order, Unlock during a background refresh, and the three Triage cases). Rechecked after the fix pass.
- **End to end:** a fresh practice project at a neutral path, run from the app on a local image built from the final code (`ribopilot:batch4final`). The first run finished with only the expected failure (SRR22396893, every piece under the minimum length); the rescue with baited SPAdes finished with no failures; exporting practice-getorganelle with Leave them out wrote 11 records, and the final ITSx check flagged one. No "Version mismatch" appeared, and every launch command carried `--ribopilot_version 0.0.0.9018`. All 8 screenshots were retaken from this run, and the Get started text was matched to them (the fourth Set up database choice shows only for a single group; the export summary shows a risk line only when a flagged record is exported). No process, container, or run lock was left.
- **Review pass:** 6 reviewers (3 on the code, 3 on the docs) found 41 issues. The one high finding was confirmed and fixed: the docs said a version mismatch stops a run before any work, but each sample fails at its first R step, after read cleaning and assembly. Batch regressions were fixed too: a refresh right after switching steps could tick rows the user never ticked; the list of shown rows could go stale after a data update; the ID column could be cut off; the version check was skipped when a source folder was visible but could not be loaded (the old guard caught that case); and the docs had dropped the note that `low_identity` is skipped silently without the reference download. Also fixed: the mismatch message sent same-version users to the image they already use and blamed a container when there was none; the Citing page's Nextflow version and run-report pointer; "GenBank hits" wording left in the app; three screenshots taken in the wrong state.

### Fifth batch: partial ends, batch trim, tests, and docs

Built on 2026-09-28 the same way: four packages in parallel worktrees (partial ends and data, tests and CI, docs, and batch trim, which the author asked for during the batch), combined, one review pass with skeptic checks, and one bounded fix pass. The version stays 0.0.0.9018 (author). The partial design copies MitoPilot (author).

| Area | Items | What changed | Main files |
|---|---|---|---|
| Partial ends per gene | action 2, DB-L8, DB-L10, DOC-42 (5' or 3' part), review A-1, A-4 | Each gene has its own 5' and 3' partial marks (`annotations.partial_start`, `partial_stop`, in the gene's own direction), as in MitoPilot. The unit window's Edit controls have 5' and 3' buttons that show what will be exported; Save keeps a change, Reset drops it, and locked units refuse it. The curator's choice always wins (author). With no choice made, an end is partial when curation found the gene truncated and that end meets no spacer, when it touches an untrimmed scaffold end, or when a trim cut through the gene (author). The unit Partial button now changes only `{completeness}` to "partial sequence". Minus-strand genes are written stop first, and `<` marks the gene's 5' end (checked with table2asn). | `R/export.R`, `R/export_opts.R`, `R/curate_rrna.R`, `R/pipeline_steps.R`, `R/app_server.R`, `R/app_modals.R`, `R/init_db.R`, `R/constants.R`, vignettes |
| Batch Trim and Restore (new, author request) | none | Trim and Restore buttons on the Assemble toolbar (every kept scaffold of the ticked samples) and the Annotation Review toolbar (ticked units); the Assemble buttons were later removed (author: the Assemble details window has no trim). Trim offers "Match gene ends" (18S to 28S) or "Gene ends + flank". Only units with all five of 18S, ITS1, 5.8S, ITS2, and 28S are trimmed; the others are named, with the reason, in a small scrolling list before anything happens. Restore undoes the trim on ticked trimmed units, and the unit window's own Restore works on batch-trimmed units. A trim that lands on gene ends marks nothing partial. The Assemble toolbar skips samples with a locked assembly, and while a run holds the project, units queued in it are skipped on both toolbars (E-6, author). | `R/pipeline_steps.R`, `R/app_server.R`, `R/app_modals.R`, `R/app_ui.R`, `R/constants.R`, `vignettes/RiboPilot.Rmd` |
| Lock prompt and re-run state | B3-app-1, DB-L1, SCI-L10, put-back coverage (part) | Lock no longer asks about units marked Exclude from export. The re-assembly clean-up now removes leftover rows for vanished scaffolds from every per-unit evidence table. An already-searched unit's GenBank best hit is re-derived after the seeds pass, so it is no longer replaced by the seed file's hit. A scaffold that stays dropped through a re-assembly keeps its coverage track when put back. | `R/app_server.R`, `R/pipeline_steps.R`, `R/coverage.R` |
| Tests and CI | T1 to T6 | The `main.nf` test checks step order (curation after the reference alignment, the chimera check after coverage and BLAST) and was shown to fail when two steps are swapped. New runs cover the depth gate, a supplied assembly, a re-queued unit, and the fallback's routing. Leftover processes are killed on timeout, and the test skips offline with an unsupported Nextflow. The clean-install job fails when the install fails and installs rstudioapi so the RStudio guard runs. | `tests/testthat/test-nf-main-run.R`, `.github/workflows/clean-install.yaml` |
| Docs | ACC-02, ACC-05 to ACC-08, ACC-10, R6 to R10, R12, R13 | A Singularity route for a self-built image; a local image no longer triggers the "not published on Docker Hub" warning, and devtools is listed for `deploy-local.sh`; the smaller wording fixes. | vignettes, `NEWS.md`, `R/utils.R`, `DESCRIPTION` |

How it was checked:
- **Tests:** 1,302 tests, 0 failures (`NOT_CRAN=true`, conda off PATH).
- **table2asn:** a minus-strand gene written `<260 100` is read as 5' partial at 260.
- **Review pass:** 5 reviewers. The docs review found nothing. The two high findings were confirmed and fixed: after a batch trim, reopening a unit showed the old coordinates, and a nudge then wrote wrong ones; and a trim on one side turned off the scaffold-end rule on both sides. Also fixed: a Restore then a new trim reused the first saved copy, so a second Restore lost edits made in between; a spacer overlapping a truncated gene by a few bases no longer counted as its neighbour; a blank flank trimmed to gene ends without a word; the gene-ends span now counts only the five genes; a curator's explicit choice could not clear a scaffold-end mark (A-1).
- **End to end:** a fresh practice project at a neutral path, run from the app on a local image built from the final code (`ribopilot:batch5`, removed afterwards). The first run failed only SRR22396893, and the baited SPAdes rescue finished with no failures. In the browser: the Partial column and the 5' and 3' buttons (Reset, Save, a locked unit refusing a forced click); the export wrote `<194 2028` and "partial sequence" for a toggled 18S and `<155 >1488` for SRR33849683_p1_s2; batch Trim on both tabs trimmed only five-feature units, listed the rest with reasons, set no partial marks, and showed current coordinates on reopening; batch Restore and the unit window's Restore brought back the original lengths; a blank flank disabled Trim; Lock did not ask about a unit marked Exclude from export. No process, container, or run lock was left.

### After batch 5: app polish from hands-on testing (2026-09-28 and 2026-09-29)

Author requests from testing the batch 5 image, each checked in the app on a copy of the author's test project; no separate review pass.

- **Sequence viewer:** MitoPilot's canvas viewer, ported without its circular, amino-acid, and read-lane code, replaces Read coverage in the Assemble and Annotation Review details windows: depth, clipped-read share (orange at the chimera check's 40% cutoff), genes with dashed partial ends, and bases when zoomed in. Clicking a gene makes it the active feature; Copy sequence copies the bases in view as FASTA. `R/seqview.R`, `inst/app/www/seqviewer.js`.
- **Chimera check:** the plot keeps the pre-trim check at full length with the cut ends shaded and genes following boundary edits; a Possible junctions table names the checks behind each junction and why; Possible junctions and All windows sit in one nested box; tied windows name the clade the tied hits share ("Chordata (50 tied hits)") instead of one arbitrary hit, in a Best Match column; the header is one call badge plus a short paragraph, with the full reason on hover.
- **Exclude & Lock** replaces Exclude from export: one click excludes and locks, and excluded units no longer count against one locked path per sample (author).
- **Tables and windows:** User Notes reach the table when a details window closes and sit before Last Updated; Partial "no" is a green pill; the Assemble toolbar loses batch Trim and Restore (its details window has no trim); Edit active feature boundaries and Reference alignment start collapsed; flagged features show orange in the Assemble rRNA features column, with the flags on hover.
- **Fixes:** the Nextflow progress window no longer repeats "Plus N more processes"; an app-wide empty-table style no longer leaves a blank strip under every table; app CSS and JS load with a file-time tag, so browsers drop stale copies after an update.
- **Tests:** 1,323 tests, 0 failures.

Still open (action 5, shipping the fish reference, was deferred by the author on 2026-09-24):
- **Container image, IMPORTANT:** the version stays 0.0.0.9018 (author), so the default container tag still points at the published 0.0.0.9018 image, which has old code. The app cannot tell (the tag matches), and container runs can fail in their R steps until that tag is rebuilt and pushed from this branch (or the version is bumped). Until then, build the image locally with `docker/deploy-local.sh` (see Upgrading in the Installation guide). Recreate existing projects as well: the schema changed again (`blast_searched.db` in batch 4, `annotations.partial_start` and `partial_stop` in batch 5). A project curated before batch 5 has explicit 0 marks that now override the scaffold-end rule (E-7).
- **Action 3:** a scaffold dropped at its very first load has no coverage track when put back (needs the reads in an annotation-only run). A kept track was measured on the earlier reads. After a rescue, the publish folder keeps the old GetOrganelle graph and log (cosmetic). A gated re-assembly skips the leftover-row clean-up (A-2). After an update, the best BLAST reference can change while its old alignment is kept (A-3).
- **Action 4:** the clean-install CI job has not run on GitHub yet (it runs on the next push).
- **Multi-user and cluster:** the shared reference folder is not safe for several users (REF-M2), concurrent taxonomy-index builds (REF-M6), clade downloads not checked (REF-M3), rebuilding a group re-queues finished samples (REF-M4), the run lock records no user (R-runs-3), a second "Submit as job" is not refused (APP-L11), and no real cluster run was tried (the Hydra fix was checked with a stand-in grid engine).
- **Reproducibility:** the version check compares the image tag only, so an image rebuilt with different code under the same version is not noticed (repro-4). Database migrations (release blocker).
- **Tests:** the fallback test checks routing only, not a finished baited assembly; triage has no timestamp, so its order after the screen is not checked; the SCI-L10 test mocks the derive step (A-5).
- **Docs:** the "stale shown rows" fix from batch 4 was checked by reading and by the scale drive, but a live case that forces it could not be staged.
- **Low findings not fixed:** the option-set toast still says "set to Ready to run" for annotation and curation sets; B3-state-5 (gate failure can still lose curator decisions on flipped or rotated scaffolds), B3-state-6, B3-state-7 (an existing but unusable label database is not refused), B3-state-8; R-data-6; R-data-8; R-app-2.
- **Not checked by hand in batch 3:** PIPE-M4, PIPE-M7, PIPE-M8 (the practice project has no supplied assemblies), SCI-M9, PIPE-M16, REF-M5, SCI-M7, and the Taxonomy lock prompt after its move to the shared confirm. All have unit tests.
- **Simplification:** groups 4 and 5 of the ponytail list, and dropping dplyr.

## What this review ran

### Reviewers and personas

All reviewers, personas, operators, and verifiers were review agents working from written roles. No outside users took part.

| Track | Who | What they did |
|---|---|---|
| Code | 9 reviewers: code-shiny (Shiny app), code-nextflow (pipeline and containers), code-science (science and algorithms), code-data (database and export), code-refdata (reference data and taxonomy), code-security, code-tests, code-ponytail (over-engineering), and code-pkg (packaging) | Read the code and reproduced problems on copies of the repo and of real projects: `shiny::testServer`, small Nextflow and Docker runs, headless Chrome, the full test suite, `R CMD check`, roxygen on a copy, lintr, and function-level coverage tracing. |
| Docs | 7 reviewers: docs-novice (field biologist, never used a terminal), docs-intermediate (grad student, 80 gastropods on Slurm), docs-hpc (research computing engineer), docs-domain (rDNA systematist and GenBank submitter), docs-editor (technical writer), docs-api (R reference auditor), and docs-appclaims (QA analyst checking claims against the app) | Followed the README and vignettes, built the pkgdown site on a copy, checked links, and checked 72 claims against the code. docs-hpc ran the practice project from the app under Singularity. |
| Operators (shared by docs and usability) | op-practice (first-time user, practice project), op-novel (grad student, gastropod genome skims), and op-edge (support engineer reproducing common mistakes) | Followed the docs literally on real projects and drove the app to a finished export or a full error sweep. |
| Usability | 5 personas: ux-novice (field biologist), ux-mitopilot (MitoPilot veteran, about 200 samples), ux-power (core facility, 500 samples per batch), ux-design (UI and WCAG 2.2 AA accessibility), and ux-curator (GenBank curator) | Worked hands-on in copies of the operators' projects and gave SUS-style estimates (average 54, range 45 to 62). |
| Verification | Separate verifier agents | Re-checked findings against the code, often re-ran them, and could confirm, refute, mark uncertain, or change severity. Every docs and usability finding and every medium or higher code finding was verified. Low and nit code findings were mostly left unverified to save usage (117 are marked unverified-low). |

### Projects and data

- **Practice project**, made with `new_test_project(executor = "docker")`: 12 samples. 10 are fish SRA read subsets (for example SRR14433899 *Brosme brosme*, SRR33849683 *Ancylopsetta cycloidea*, SRR22396893 *Sagamichthys abei*), plus the synthetic CHIM00 (clean control) and CHIM01 (fish 18S joined to a *Penicillium* 28S).
- **Gastropod project**: PRJEB115569 North Sea genome skims, first 3 million read pairs each of ERR17620304 (*Euspira nitida*) and ERR17620292 (*Tritia reticulata*), with two order-level reference databases (Littorinimorpha, Neogastropoda) built from GenBank in the app.
- **Mistake sweep** (op-edge): the bundled CHIM00 and CHIM01 reads, about 60 hand-made mapping files, synthetic assemblies (rDNA in a scaffold, no rRNA, CRLF, gzip, and a FASTQ named as an assembly), and projects with spaces in the path, a bad container tag, and the conda default.
- **Scale**: an inflated 504-sample copy of the practice project database, used to time the app (not to run the pipeline).
- **Exports** were checked with NCBI table2asn 1.29 and against NCBI's help pages. Nothing was submitted to NCBI.

### Setup

Linux workstation; Docker image `macguigand/ribopilot:0.0.0.9018`; Nextflow launcher 26.04.6, pinned by RiboPilot to 25.10.6 through `NXF_VER`; the development code loaded with `pkgload::load_all()`. The app was driven through headless Chrome with real clicks and key presses. Every pipeline run was started from the app's Update and Run from App buttons; none were scripted.

### Time

- Practice project: the first run aborted after about 1 minute (fish reference 404); the second run took about 16 minutes (169 tasks); a baited SPAdes rescue took about 1 minute; each export about 25 seconds.
- Gastropod project: taxonomy about 1 minute, reference build 1 minute 20 seconds, first assembly about 3 minutes, rescues about 1 to 2 minutes.
- Singularity practice run (docs-hpc): 44 tasks in about 5.5 minutes for 3 samples.
- Test suite: 230 seconds; `R CMD check` with CI flags: 217 seconds.
- Reported hands-on time: op-practice about 150 minutes over two sessions, op-novel about 130 minutes over three, op-edge about 30 minutes in its last session plus earlier attempts.
- File times in the review workspace fall between about 12:00 and 19:20 on 2026-09-23, with a pause of several hours in the afternoon.

### Resumed once

The first attempt was stopped by a usage limit, and the review was resumed once. Reviewers reused sound earlier work: the docs-novice and docs-intermediate reviews come from the first attempt (verified then); some operators and personas worked in two or three sittings; the code lead re-checked the earlier draft against the merged data (all 333 ids, statuses, and severities); the docs lead re-mapped every finding onto current ids; and the docs screenshots were checked against the operators' originals. To save usage, low and nit code findings were left unverified (unverified-low).

### Counts

| Track | Source findings | Confirmed | Uncertain | Unverified-low | Refuted | Merged items |
|---|---|---|---|---|---|---|
| Code | 333 | 211 | 2 | 117 | 3 | 217: 17 high, 55 medium, 105 low, 6 nit, plus 34 ponytail (low or nit) |
| Docs | 247 (227 from docs reviewers, 20 from operators) | 246 | 0 | 0 | 1 | 99: 5 blocker, 43 major, 40 minor, 11 polish |
| Usability | 194 (78 from operators, 116 from personas) | 192 | 0 | 0 | 2 | 87: 3 blocker, 30 major, 44 minor, 10 polish |
| **Distinct findings** | **754** (the 20 operator findings in docs also appear in usability) | **629** | **2** | **117** | **6** | 403 items across tracks, many describing the same problem (see Appendix B) |

No finding in any track ended as disputed. No code finding is critical; four came in as critical and verifiers lowered them to high. Severity scales differ by track: code uses critical, high, medium, low, and nit; docs and usability use blocker, major, minor, and polish.

How to read the status labels: **confirmed** means a verifier re-checked it and agrees; **uncertain** means the code behaves as described but harm was not shown in practice; **unverified-low** means a low or nit code finding that no verifier checked (a strong lead, not a proven bug); **refuted** findings are listed only in Appendix A.

## Top 15 actions

Ranked by: wrong data reaching GenBank or the database first, then what stops new users, then how many reviewers hit it. Each action names every merged item and source id; "Read more" links go to the full write-ups. Status is confirmed unless an id says otherwise.

### 1. Stop contaminated and chimeric records from leaving without a check

- **Ids:** UX-31 (U-curator-1, O-novel-1, U-curator-7); UX-14 (U-novice-1, U-mitopilot-7, U-curator-8, O-practice-12, U-novice-13, U-curator-19); UX-22 (U-design-2, U-curator-9, U-curator-12); UX-34 (U-novice-4, U-design-10, U-mitopilot-21, U-curator-18 in part); UX-86 (U-curator-16).
- **Severity:** blocker (UX-31); major (UX-14, UX-22, UX-34); minor (UX-86). All confirmed.
- **Where:** `risk_flags()` and `submission_risk` in `R/export.R`; `export_check_summary_ui()` in `R/app_modals.R:676-679`; the Assemble table (no Screen or Chimera column) and the Annotation Review table (those columns are 17 to 23 of 26, off screen); the Lock handlers.
- **Problem:** Tret_01, a whelk sample whose rDNA is 99.6 percent sea anemone and screened off_target, passed both locks and the Export window with no prompt, and went out as `>Tret_01 [organism=Tritia reticulata] ... complete sequence` under a green "Final ITSx check: nothing flagged". `risk_flags()` never reads the screen verdict, and the tables show the same green check for CHIM01 (fish plus *Penicillium*) as for clean samples.
- **Why it matters:** this is the most direct route to a wrong public record, and the app already knows the answer.
- **Fix:** add the screen verdict and the contaminant flag to `submission_risk`, the summary CSV, and the report. Before lock and before export, list units that are off_target, contaminant, or chimera suspect or likely and not yet marked Reviewed, and have the user confirm or exclude each one (exclude by default). Add one "Checks" column right after the ID in Assemble and Annotation Review, name each flagged record in the export dialog, and never show a green all-clear while such a flag is open. Show the agreement rank on the on_target badge.
- **Read more:** [UX-31](USABILITY_REVIEW.md#ux-31), [UX-14](USABILITY_REVIEW.md#ux-14), [UX-22](USABILITY_REVIEW.md#ux-22), [UX-34](USABILITY_REVIEW.md#ux-34), [UX-86](USABILITY_REVIEW.md#ux-86).
- **Status:** fixed on this branch. Screen and contaminant flags reach `submission_risk` and the Export window asks (first batch); the Checks column, the Lock prompt, and the named records in the Export complete dialog came in the [third batch](#third-batch-re-runs-caching-blast-and-visibility).

### 2. Mark truncated genes as partial on export

- **Ids:** DB-H2 (C-data-2, C-science-7); UX-32 (U-curator-2, O-practice-3); DOC-42 (D-domain-1); DB-L8 (C-data-21).
- **Severity:** high (DB-H2); blocker (UX-32); major (DOC-42); low (DB-L8). All confirmed.
- **Where:** `R/export_opts.R:134-141` (definition line), `R/export.R:144-153` (feature table `<` and `>`), `R/export.R:200-222` (`risk_flags`), `R/curate_rrna.R:115-120`; the Partial tooltip at `R/app_server.R:1933`.
- **Problem:** a feature is called partial only when it touches a scaffold end or the curator pressed Partial. Model coverage and the curation "truncated" flag are ignored, so SRR33849683_p1_s2 (28S at 46 percent of the model) and CHIM01's 18S (27 percent) were written as "complete sequence" with no `<` or `>`. When Partial is used, it marks both ends even if only one is incomplete.
- **Why it matters:** public records claim complete genes that are not; table2asn flags SHORT_RRNA on three records, and the export gives no warning.
- **Fix:** treat an end as partial when its model coverage is below `min_anchor_coverage` or the feature is flagged truncated, and use the same rule for the definition line and the feature table. Let Partial choose 5', 3', or both. Add a "partial not marked" submission risk when the curator's setting disagrees. Explain the rule in the export docs.
- **Read more:** [CODE 4.4](CODE_REVIEW.md#44-database-and-export) (DB-H2, DB-L8), [UX-32](USABILITY_REVIEW.md#ux-32), [DOC-42](DOCS_REVIEW.md#doc-42).
- **Status:** DB-H2 and UX-32 are fixed on this branch: truncated features export as partial in both the `.tbl` and the definition line. The 5' and 3' marks per gene, DOC-42, and DB-L8 came in the [fifth batch](#fifth-batch-partial-ends-batch-trim-tests-and-docs).

### 3. Clear a unit's old evidence when its sequence changes

- **Ids:** DB-H1 (C-data-1, C-tests-3); PIPE-H1 (C-nextflow-1, C-data-5); SCI-H1 (C-science-3, C-tests-4); SCI-M5 (C-science-8, C-tests-2); UX-17 (O-practice-9, U-novice-6, U-power-1); PIPE-M1 (C-nextflow-5, C-refdata-14); DB-L1 (C-tests-5); SCI-L10 (C-science-24).
- **Severity:** high (DB-H1, PIPE-H1, SCI-H1); medium (SCI-M5, PIPE-M1); major (UX-17); low (DB-L1, SCI-L10). Confirmed, except C-science-24 (unverified-low).
- **Where:** `R/pipeline_steps.R:384-447` (`load_assemblies()`: carry-over by path and scaffold number, orphan cleanup) and `:200-206` (`gate_failed()`); `inst/nextflow/bin/rp_score.R:11` and `inst/nextflow/modules/pipeline_r.nf:123` (score reads the publish folder); `R/screen.R:329, 367`; `inst/nextflow/modules/blast.nf:104-107`.
- **Problem:** when a sample is assembled a second time (after an option change, a rescue, or an unlock), the pipeline can load old path files left in `out/<ID>/assemble`, keep a locked unit's old gene coordinates and export them onto the new sequence, keep the old contamination verdict, skip BLAST because the sequence md5 matches even though the settings changed, bring a curator's dropped scaffold back as "keep", and give new rescue scaffolds the old "too short" drop flags so the sample stays "In progress".
- **Why it matters:** silent wrong data reaches export, curator decisions are lost, and the documented rescue does not recover.
- **Fix:** (1) have scoring and loading read the files staged into the task, not the publish folder (stopgap: delete `out/<ID>/assemble/<ID>_assembly_*` at the start of each assembly). (2) In `load_assemblies()`, compare each surviving unit's new sequence with the old one; when it changed, delete its annotations, BLAST, screen, chimera, and coverage rows, unlock and un-review it, and re-queue it. The simpler guard is to refuse to re-assemble a sample that has an annotation-locked unit. (3) Carry all four triage columns keyed on a sequence hash, and skip the carry-over when the assembler changes. (4) Store what the screen verdict and the "already searched" mark were computed from (sequence md5, taxid, BLAST settings) and recompute when they differ. Add tests that go through the functions the app calls (TEST-M1 explains why the current tests missed this).
- **Read more:** [CODE 4.4](CODE_REVIEW.md#44-database-and-export) (DB-H1, DB-L1), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-H1, PIPE-M1), [CODE 4.3](CODE_REVIEW.md#43-science-and-algorithms) (SCI-H1, SCI-M5, SCI-L10), [UX-17](USABILITY_REVIEW.md#ux-17).
- **Status:** fixed on this branch except DB-L1 and SCI-L10. Fix (2) came in the first batch (DB-H1, SCI-H1 re-assembly half), the taxon half of SCI-H1 in the second, and fixes (1), (3), and (4) in the [third batch](#third-batch-re-runs-caching-blast-and-visibility) (PIPE-H1, SCI-M5, UX-17, PIPE-M1).

### 4. Make the first commands work on a new machine

- **Ids:** PKG-H1 (C-pkg-1) and DOC-17 (D-novice-2, D-intermediate-3, D-hpc-2, D-editor-3, D-api-2); PKG-H2 (C-data-3, C-pkg-2), DOC-18 (D-novice-3, D-intermediate-2, D-editor-2, D-api-1), and UX-01 (O-practice-1, O-edge-1, O-novel-2); PIPE-M11 (C-nextflow-20, C-pkg-30), DOC-16 (D-novice-9, D-hpc-10, D-novice-34, D-editor-20, D-hpc-25), and UX-04 (O-edge-8); PIPE-H2 (C-nextflow-7) and DOC-03 (D-novice-4, D-hpc-12, D-intermediate-10, D-editor-6 in part).
- **Severity:** high (PKG-H1, PKG-H2, PIPE-H2); blocker (DOC-17, DOC-18, UX-01); medium (PIPE-M11); major (DOC-16, UX-04, DOC-03). Confirmed, except C-pkg-30 (unverified-low).
- **Where:** `DESCRIPTION` (no `biocViews`); `R/init_project.R:98-105` and `R/init_test_project.R:158, 166, 199` (RStudio check); `R/init_project.R:52` and `new_test_project()` (default `executor = "conda"`); `inst/config.local`, `inst/config.docker`, and `inst/config.conda` (no `resourceLimits`); `inst/nextflow/nextflow.config:63-66` (memory doubles on retry).
- **Problem:** `remotes::install_github()` cannot find Biostrings, DECIPHER, or pwalign. `new_project()` and `new_test_project()` stop with "RStudio not running" in a terminal, Rscript, VS Code, Positron, or an HPC session, because `!isFALSE(Sys.getenv("RSTUDIO", FALSE))` is always TRUE; the test project fails only after downloading every read. The default executor is conda, which the docs call unsupported and whose `environment.yml` lacks jq, sqlite3, and R, so every preprocess task fails and the run still says Done. On the local executor a task that asks for more memory than the machine has ends the whole run, and the defaults ask 24 GB per assembly (then 48 and 72 GB on retry).
- **Why it matters:** every new user hits at least one of these in the first hour, and terminal and HPC users are the main audience for genome skims.
- **Fix:** add an empty `biocViews:` line; test `requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()` and fix the help text; default `executor` to "docker" and warn once on conda; add `process.resourceLimits` with the host's CPUs and memory, filled in by `new_project()`. Add a clean-machine check to CI (see theme 4).
- **Read more:** [CODE 4.7](CODE_REVIEW.md#47-packaging) (PKG-H1, PKG-H2), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-H2, PIPE-M11), [DOC-17](DOCS_REVIEW.md#doc-17), [DOC-18](DOCS_REVIEW.md#doc-18), [DOC-16](DOCS_REVIEW.md#doc-16), [DOC-03](DOCS_REVIEW.md#doc-03), [UX-01](USABILITY_REVIEW.md#ux-01), [UX-04](USABILITY_REVIEW.md#ux-04).
- **Status:** fixed on this branch (see [Fixed on this branch](#fixed-on-this-branch-after-the-review)). The clean-machine CI check was added in the [fourth batch](#fourth-batch-working-at-scale-reproducibility-and-docs); it has not run on GitHub yet.

### 5. Ship the built-in fish reference, and treat "repo public" as a release gate

- **Ids:** DOC-19 (D-hpc-3, D-novice-19, D-editor-7, D-api-22, D-intermediate-31 in part, O-novel-7 in part); UX-11 (O-practice-2, O-edge-14, O-novel-7, O-practice-22, O-practice-23); REF-L4 (C-refdata-15, C-security-12, C-data-29, C-pkg-20); SCI-M9 (C-science-18); REF-L9 (C-refdata-17); UX-56 (O-edge-24, O-edge-23).
- **Severity:** blocker (DOC-19); major (UX-11); medium (SCI-M9); low (REF-L4, REF-L9); minor (UX-56). Confirmed, except C-pkg-20 and C-refdata-17 (unverified-low).
- **Where:** `R/constants.R:20-29` and `R/init_db.R:51-52, 102` (raw `main`-branch URLs); `R/score_project.R:11-17` (`local_ref()`); `R/db_dist.R:142`; `README.Rmd`; `docker/README.md:107-109`.
- **Problem:** the default fish seeds and labels are raw.githubusercontent.com URLs on `main`, which return 404 while the repo is private; no token helps, because Nextflow and `download.file()` fetch them without auth. The practice project's first run aborts (or scoring fails for every sample) and the app still says Done. Install and the practice reads need `GITHUB_PAT`, which no doc mentions, and the error first blames a slow connection. After release, the unpinned URL means any push to `main` silently changes the seeds for existing projects.
- **Why it matters:** every pre-release tester stops here; after release, results can drift with no version change.
- **Fix:** ship the fish seeds and labels inside the image, or as checksummed release assets pinned to the package version, and have `new_db()` use the constants. Wrap the seed download in `tryCatch` and return NA. Check URL databases before launch. Until the repo is public, add a short private-preview note, and report a 404 as "not found or private; set GITHUB_PAT".
- **Read more:** [DOC-19](DOCS_REVIEW.md#doc-19), [UX-11](USABILITY_REVIEW.md#ux-11), [UX-56](USABILITY_REVIEW.md#ux-56), [CODE 4.5](CODE_REVIEW.md#45-reference-data-and-taxonomy) (REF-L4, REF-L9), [CODE 4.3](CODE_REVIEW.md#43-science-and-algorithms) (SCI-M9).
- **Status:** deferred by the author (2026-09-24). Most users will build their own reference databases, so shipping the built-in fish reference is not a priority or a release gate.

### 6. Rewrite Get started around the app's five steps

- **Ids:** DOC-28 (D-novice-1, D-intermediate-1, D-editor-1, D-appclaims-1, D-novice-10, D-intermediate-5, D-editor-4, D-api-5, O-practice-5, O-practice-6, O-novel-6, O-edge-27, D-hpc-6 in part) and UX-49 (O-practice-5, O-practice-6, O-novel-6, O-edge-27, U-novice-7, U-mitopilot-13, U-power-17); DOC-20 (D-intermediate-4, D-editor-5, D-appclaims-2, D-api-5 in part) and UX-09 (O-edge-11, O-practice-7, O-practice-6 in part); DOC-37 (D-appclaims-8, O-practice-10, D-novice-10 in part, D-intermediate-16 in part), UX-16 (O-practice-10, U-power-14, U-novice-12, U-novice-6 in part), and PIPE-H5 (C-pkg-5, C-nextflow-22, C-ponytail-3); DOC-38 (D-appclaims-5); DOC-40 (D-novice-21, D-intermediate-8, D-api-6, D-appclaims-3, O-practice-14, O-practice-26) and UX-37 (O-practice-14); DOC-30 (D-novice-12, D-intermediate-20, D-editor-9, D-api-9, D-novice-41); DOC-35 (D-novice-13, D-appclaims-18, D-domain-22, O-practice-19, D-api-26, D-intermediate-19 in part, D-editor-21 in part).
- **Severity:** blocker (DOC-28); high (PIPE-H5); major (UX-49, DOC-20, UX-09, DOC-37, UX-16, DOC-38, DOC-40, UX-37, DOC-30, DOC-35). All confirmed.
- **Where:** `vignettes/RiboPilot.Rmd:51-57, 140-173, 245-310, 294-296, 342-347, 385-394`; `R/nf_command.R:163-176` and `inst/nextflow/modules/pipeline_r.nf:21-40` (no taxonomy-lock filter); `R/run_fallback.R:13-17`.
- **Problem:** the guide describes three tabs and runs `run_assembly()`, `run_fallback()`, and `export_project()` from the console right after `new_project()`. The app has five steps (Taxonomy, Reference, Assemble, Annotation Review, Export), opens on Taxonomy, and runs everything from Update and Run from App. Because the pipeline ignores the taxonomy lock, following the guide assembles samples with no lineage and then shows an empty Assemble table. `run_fallback()` also re-assembles good samples that have one fragment scaffold, `export_project()` stops unless groups were made in the app, and `RiboPilot("~/my_project")` fails without `library()`.
- **Why it matters:** 9 of the 10 people in the docs track and 6 of the 8 in the usability track hit this. It is the main reason nobody finished a project from the docs, and it works against the app-first design.
- **Fix:** rewrite Steps 2 to 9 as an app walkthrough, following the [proposed outline](DOCS_REVIEW.md#proposed-outline-for-the-published-site), with one screenshot per step (15 candidates are in `screenshots/docs_*.png`; retake them at a neutral path after the fixes), a "What you should see" box for the practice project, and a "keep the tab open" warning. Move the console functions to an "Advanced: scripting" page. In code, require `taxonomy_lock = 1` in the pipeline's sample query and in the Update count, and delete `run_fallback()` (the Assembler option already does the rescue) or restrict it to samples with no complete kept unit.
- **Read more:** [DOC-28](DOCS_REVIEW.md#doc-28), [UX-49](USABILITY_REVIEW.md#ux-49), [DOC-20](DOCS_REVIEW.md#doc-20), [UX-09](USABILITY_REVIEW.md#ux-09), [DOC-37](DOCS_REVIEW.md#doc-37), [UX-16](USABILITY_REVIEW.md#ux-16), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-H5), [DOC-38](DOCS_REVIEW.md#doc-38), [DOC-40](DOCS_REVIEW.md#doc-40), [DOC-30](DOCS_REVIEW.md#doc-30), [DOC-35](DOCS_REVIEW.md#doc-35).
- **Status:** fixed on this branch (second batch): taxonomy-lock gate, `run_fallback()` deleted, `export_project()` without groups, and the Get started rewrite. See [Second batch](#second-batch-actions-6-to-15).
- **Status:** docs rewritten on this branch (see [Fixed on this branch](#fixed-on-this-branch-after-the-review)). The "Next: open the app" message at the end of `new_project()` (DOC-35 code half) is still open.

### 7. Keep the taxon right from the Taxonomy step to the screen

- **Ids:** REF-H1 (C-refdata-1); REF-H3 (C-refdata-3) and UX-08 (O-edge-16, O-edge-17); REF-H2 (C-refdata-2, C-pkg-4, C-security-18) and DOC-24 (D-hpc-15, D-novice-29, D-editor-36); REF-M1 (C-refdata-4, C-security-9, C-ponytail-34); APP-L16 (C-shiny-31).
- **Severity:** high (REF-H1, REF-H2, REF-H3); medium (REF-M1); major (UX-08); minor (DOC-24); low (APP-L16). Confirmed, except C-shiny-31 (unverified-low).
- **Where:** `R/tax_search.R:61, 80` and `R/taxonomy_ranks.R:219-230` (lineage taxids); `R/tax_index.R:294` and `R/taxonomy.R:122, 245-265` (homonyms); `R/taxonomy.R:54-93` and `R/tax_index.R:137, 177` (taxdump download and index); `R/reference_counts.R:17-22` and its callers (clade queries by name).
- **Problem:** correcting a sample's taxon leaves `lineage_taxids` pointing at the old organism, so the screen scores true host hits as contaminants. A name that matches several taxa takes the lowest taxid ("Morus" becomes the mulberry, not the gannet) with only a grey badge, and the alert counts unresolved samples as assigned. A slow taxdump download is kept half-finished, and the index built from it is trusted forever. Clade databases are fetched by name, so the wheatear genus *Oenanthe* pulls in about 91 percent plant rDNA.
- **Why it matters:** the screen, reference groups, clade databases, and scoring all run against the wrong organism, silently.
- **Fix:** write `lineage_taxids` in the same update as `taxid` (or set it to NULL so it is rebuilt), and clear the sample's screen rows. Leave homonyms unresolved, store and show the alternatives with their kingdom, and report counts by confidence. Download the taxdump to a temp file with a long timeout and NCBI's md5, stop when untar fails, and require a node-count floor before trusting the index. Pass the stored `clade_taxid` to every clade query and to `download_rrna.sh --taxid`. Warn before locking unresolved or ambiguous samples.
- **Read more:** [CODE 4.5](CODE_REVIEW.md#45-reference-data-and-taxonomy) (REF-H1, REF-H2, REF-H3, REF-M1), [UX-08](USABILITY_REVIEW.md#ux-08), [DOC-24](DOCS_REVIEW.md#doc-24), [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-L16).
- **Status:** fixed on this branch (second batch). See [Second batch](#second-batch-actions-6-to-15).

### 8. Tell the truth at the end of a run and when a sample fails

- **Ids:** APP-M2 (C-shiny-3, C-nextflow-12), UX-42 (O-edge-2, U-power-6, U-design-9, O-practice-4, O-novel-5), and DOC-88 (D-appclaims-15, D-appclaims-21, in part); UX-15 (O-novel-4, U-novice-5, O-practice-11, U-power-8, U-curator-14) and DOC-08 (D-intermediate-12, D-novice-45, D-intermediate-33, O-edge-27 in part, D-editor-17 in part); UX-45 (U-power-10, U-mitopilot-10); APP-M3 (C-shiny-2); PIPE-M6 (C-nextflow-14); PIPE-M9 (C-nextflow-17).
- **Severity:** medium (APP-M2, APP-M3, PIPE-M6, PIPE-M9); major (UX-42, UX-15, DOC-08); minor (UX-45, DOC-88). All confirmed.
- **Where:** `R/app_server.R:1009-1018` and `R/app_reference.R:962-968` (the "Workflow run finished." alert); `R/app_server.R:982-998` (`nf_start`); the details window (`R/app_modals.R:748`); `R/assemble_fallback.R:64-91`; `inst/nextflow/main.nf:91`.
- **Problem:** every run ends with a green "Done. Workflow run finished.", including aborted runs, runs where every task failed, and runs with nothing queued, because the exit status is never read. Failed samples show "failed assembly" or nothing, the real cause sits in the work folder, and the details window says "the remote BLAST is still running" after the run has ended. Samples below minDepth are skipped with no note, and a SPAdes crash exits 0. With Nextflow missing, Run from App ends the session.
- **Why it matters:** this was the most repeated complaint (5 of 8 usability reviewers). Users trust the popup and move on.
- **Fix:** read `p$get_exit_status()` and the failed and ignored task counts. On failure, show an error alert with the first "Caused by" line and the log path; on success, list per-sample results ("11 succeeded, 1 failed: SRR22396893"). Refuse to start with nothing queued. Wrap the launch in `tryCatch`. Write a plain note for each failure (gate, minDepth, SPAdes error lines), show it as a banner in the details window with the next step ("Try the baited SPAdes assembler in Assembly Options"), and add a Troubleshooting page.
- **Read more:** [UX-42](USABILITY_REVIEW.md#ux-42), [UX-15](USABILITY_REVIEW.md#ux-15), [UX-45](USABILITY_REVIEW.md#ux-45), [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-M2, APP-M3), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-M6, PIPE-M9), [DOC-08](DOCS_REVIEW.md#doc-08), [DOC-88](DOCS_REVIEW.md#doc-88), [error message tables](USABILITY_REVIEW.md#error-messages-from-the-user-mistake-sweep).
- **Status:** fixed on this branch (second batch), with raw exit codes and a copyable report as the author asked. The failed-task lines with a real exit code were checked by unit tests, not by hand. See [Second batch](#second-batch-actions-6-to-15).

### 9. Fix the app bugs that end the session or overwrite work

- **Ids:** APP-H1 (C-shiny-1); APP-H2 (C-shiny-5); APP-M4 (C-shiny-4); APP-M8 (C-shiny-13) and UX-18 (O-novel-3); APP-M14 (C-ponytail-1); DB-M4 (C-data-13, C-shiny-9).
- **Severity:** high (APP-H1, APP-H2); medium (APP-M4, APP-M8, APP-M14, DB-M4); major (UX-18). All confirmed.
- **Where:** `R/app_opts.R:384` (options lookup); `R/app_server.R:2158, 2217-2220` (trim span); `R/app_opts.R:434-442, 487, 568-586` (option sets); `R/curate_rrna.R:351-353` (notes); `R/export.R:25, 36, 76, 539, 742` and `R/nf_command.R:163-165` (database opened with no busy timeout).
- **Problem:** clicking the Annotation or Curation options cell on any sample with more than one unit ends the Shiny session (3 of 12 practice samples). A second trim reuses the old span, so "Annotated + flank" can delete 18S and ITS1. Picking another existing set while Edit is ticked overwrites that set. Saving a new set name with Edit unticked points samples at a set that was never written. Curation overwrites the curator's Notes. Assign Group and export fail with "database is locked" (and can end the session) while the pipeline writes.
- **Why it matters:** crashes, lost settings, and lost notes during normal use.
- **Fix:** add `ORDER BY x.path, x.scaffold LIMIT 1`; read `edit_tick()` in `trim_span` and check the range at Apply; reload a known set's values and untick Edit; create or refuse an unknown set; give curation its own notes column; open the database with `open_project_db()` and wrap app writes in `tryCatch`.
- **Read more:** [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-H1, APP-H2, APP-M4, APP-M8, APP-M14), [CODE 4.4](CODE_REVIEW.md#44-database-and-export) (DB-M4), [UX-18](USABILITY_REVIEW.md#ux-18).
- **Status:** fixed on this branch (second batch). See [Second batch](#second-batch-actions-6-to-15).

### 10. Fix scaffold and path bookkeeping

- **Ids:** SCI-H3 (C-science-10); SCI-M8 (C-science-9); SCI-H2 (C-science-5); SCI-M6 (C-science-6, C-tests-6); PIPE-H4 (C-nextflow-2); PIPE-M5 (C-nextflow-10).
- **Severity:** high (SCI-H3, SCI-H2, PIPE-H4); medium (SCI-M8, SCI-M6, PIPE-M5). All confirmed.
- **Where:** `R/select_paths.R:123-124` and `R/pipeline_steps.R:24` (path numbers from `Sys.glob()` order); `R/coverage.R:65`; `R/pipeline_steps.R:376` (no rotation); `R/annotate.R:105` (`pick()`); `inst/nextflow/modules/find_rrna.nf:59`; `inst/nextflow/modules/pipeline_r.nf:26`.
- **Problem:** with 10 or more paths (the default cap is 10, and one practice sample has exactly 10), scores land on the wrong path and differ between the host and the container. Scaffold "1.10" is read as the number 1.1, so two scaffolds share one depth track. Circular units are never rotated, so a start point inside 18S or 28S gives a broken model that trim cannot fix. The feature model takes the first hit per gene, so a partial leading copy replaces the real anchor and ITS2 is lost. A supplied assembly with more than 5,000 contigs aborts the whole run. Max paths 0 (documented as "no limit") crashes GetOrganelle.
- **Why it matters:** wrong scores, depth, and gene models on normal data, with no error.
- **Fix:** read the path number from the file name; `read.table(..., colClasses = c("character", "integer", "integer"))`; rotate circular units to start in the IGS before annotation (and add a rotate action); pick the nearest 18S upstream and 28S downstream of the chosen 5.8S or ITS; `path(gffs, stageAs: 'chunk_*.gff')`; omit the max-paths flag for 0, and let GetOrganelle export max_paths + 1 so the gate can see overflow.
- **Read more:** [CODE 4.3](CODE_REVIEW.md#43-science-and-algorithms) (SCI-H2, SCI-H3, SCI-M6, SCI-M8), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-H4, PIPE-M5).
- **Status:** fixed on this branch (second batch); an app Rotate button was skipped (author). See [Second batch](#second-batch-actions-6-to-15).

### 11. Get the rest of each GenBank record right

- **Ids:** DB-M1 (C-data-4) and DOC-45 (D-domain-5, D-intermediate-7 in part); UX-36 (U-curator-4, O-practice-15) and DOC-44 (D-domain-13, D-appclaims-13, D-intermediate-7); DB-M5 (C-data-17, C-science-16), UX-84 (U-curator-3, O-novel-12), and DOC-39 (D-domain-4); DOC-43 (D-domain-12); DB-L11 (C-data-24) and DOC-49 (D-domain-16); DB-L14 (C-data-28); DOC-46 (D-novice-38, D-domain-23).
- **Severity:** medium (DB-M1, DB-M5); major (DOC-45, UX-36, DOC-44, UX-84, DOC-43); minor (DOC-39, DOC-49, DOC-46); low (DB-L11, DB-L14). All confirmed.
- **Where:** `R/export_opts.R:18, 118-122, 164, 219-254`; `R/annotate.R:97-100, 167-180` (`clip_spacer`); `R/export.R:110-158, 136-141, 623-627`.
- **Problem:** `[organism=]` uses the raw mapping Taxon, not the name confirmed in the Taxonomy step. The template refuses any "isolate", yet two scaffolds of one fish get made-up isolates (`ID_p1_s1`, `ID_p1_s2`), and NCBI's required source modifiers are listed nowhere (the guide spells them with underscores). ITS spacers stop 1 to 4 bp short of the rRNA genes on almost every record, which gave 19 ITSdoesNotAbutRRNA warnings in the curator's table2asn checks. Exports carry unannotated ETS and IGS flanks, every N-run claims "align genus" linkage evidence, and empty modifiers are written with no warning.
- **Why it matters:** GenBank staff will ask about each of these, some records may be removed at the Portal, and some source data in public records would be wrong.
- **Fix:** default the organism to `sci_name` when the match was exact, and flag units where Taxon and `sci_name` differ (keep labels like "Upeneus sp. A"). Allow a curator isolate or voucher, and tell variant records apart with a note, not an isolate. Snap the 18S end and 28S start to the ITSx spacers (or close gaps of 5 bp or less and report it), and add a gap flag. Use "unspecified" linkage evidence unless it is known. Warn on empty modifiers. Write a "Preparing a GenBank submission" page: modifiers in NCBI's spelling, which files to upload, trimming flanks, and the partial rule.
- **Read more:** [CODE 4.4](CODE_REVIEW.md#44-database-and-export) (DB-M1, DB-M5, DB-L11, DB-L14), [UX-36](USABILITY_REVIEW.md#ux-36), [UX-84](USABILITY_REVIEW.md#ux-84), [DOC-39](DOCS_REVIEW.md#doc-39), [DOC-43](DOCS_REVIEW.md#doc-43), [DOC-44](DOCS_REVIEW.md#doc-44), [DOC-45](DOCS_REVIEW.md#doc-45), [DOC-46](DOCS_REVIEW.md#doc-46), [DOC-49](DOCS_REVIEW.md#doc-49).
- **Status:** fixed on this branch (second batch), keeping the made-up isolates and the unannotated flanks (author). See [Second batch](#second-batch-actions-6-to-15).

### 12. Make the cluster route work, and write it down

- **Ids:** PIPE-H3 (C-nextflow-3, C-data-10, C-pkg-9) and DOC-29 (D-hpc-1); PIPE-M14 (C-nextflow-23, C-security-13); APP-M10 (C-shiny-16, C-pkg-10) and UX-54 (U-mitopilot-15, U-power-9); PIPE-L2 (C-nextflow-24); PIPE-M10 (C-nextflow-18) and DOC-71 (D-hpc-8); DOC-04 (D-hpc-6, D-intermediate-9, D-hpc-19, D-hpc-21); DOC-05 (D-hpc-4, D-novice-23, D-intermediate-27); DOC-06 (D-hpc-5); DOC-69 (D-api-7, D-novice-26, D-intermediate-26) and UX-79 (O-novel-15, O-edge-29, U-mitopilot-20, U-power-22, U-design-30).
- **Severity:** high (PIPE-H3); blocker (DOC-29); medium (PIPE-M14, APP-M10, PIPE-M10); major (DOC-04, DOC-05, DOC-06, DOC-69); minor (UX-54, DOC-71); low (PIPE-L2); polish (UX-79). All confirmed.
- **Where:** `R/generate_config.R:32-38` (rawDir bind); `R/app_run.R:27, 82-100` and `R/hydra.R:21-24` (job scripts); `R/app_utils.R:26-29` (Work Dirs); `inst/config.awsbatch`; `man/RiboPilot.Rd` and `man/run_app.Rd`; `R/app_run.R:194`.
- **Problem:** every Singularity or Docker profile made by `generate_config()` mounts `rawDir`, which is NA for assembly-only projects, so Nextflow cannot even parse the config. Head-job scripts fit only Hydra-style SGE: PBS ignores them and kills the head job at the default walltime, and SLURM log names print a literal `%j`. Work Dirs finds nothing for cluster runs. Generated profiles have no pull timeout and no shared image cache. awsbatch cannot run as written. `launch.browser = FALSE`, which the SSH-tunnel route needs, silently turns on headless mode and hides Run from App. No page covers tunnels, caches, network needs, or shared storage.
- **Why it matters:** cluster users are a main audience, and a research computing engineer had to read the code to get a working run.
- **Fix:** drop the rawDir bind; write scheduler-specific job headers from the profile (partition, account, walltime), export `NXF_OPTS`, and use a single `%`; read `nextflow_logs/*.nextflow.log` and match any `*TaskHandler`; add `pullTimeout = '60 min'` and a shared `cacheDir`; remove awsbatch until it is supported; fix the headless help text and console hint; add "Running on a cluster" and "Disk space, caches, and offline use" pages.
- **Read more:** [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-H3, PIPE-M10, PIPE-M14, PIPE-L2), [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-M10), [DOC-29](DOCS_REVIEW.md#doc-29), [DOC-04](DOCS_REVIEW.md#doc-04), [DOC-05](DOCS_REVIEW.md#doc-05), [DOC-06](DOCS_REVIEW.md#doc-06), [DOC-69](DOCS_REVIEW.md#doc-69), [DOC-71](DOCS_REVIEW.md#doc-71), [UX-54](USABILITY_REVIEW.md#ux-54), [UX-79](USABILITY_REVIEW.md#ux-79).
- **Status:** fixed on this branch (second batch), but no real cluster run was tried. See [Second batch](#second-batch-actions-6-to-15).

### 13. Keep runs alive when the browser goes away, and stop them cleanly

- **Ids:** APP-M1 (C-shiny-6, C-nextflow-8, C-security-6, C-pkg-14), UX-43 (O-edge-10, U-power-4), UX-44 (U-power-5), and DOC-36 (D-appclaims-6, D-novice-11); APP-M9 (C-shiny-15) and UX-46 (O-edge-9); APP-L11 (C-shiny-26); REF-L17 (C-refdata-27); PIPE-M2 (C-nextflow-9, C-shiny-32).
- **Severity:** medium (APP-M1, APP-M9, PIPE-M2); major (UX-43, UX-44, DOC-36); minor (UX-46); low (APP-L11, REF-L17). Confirmed, except C-shiny-26, C-refdata-27, and C-shiny-32 (unverified-low).
- **Where:** `R/app_server.R:1022-1033` (`kill_nf`, `onSessionEnded`); `R/app_reference.R:973-983`; `R/app_lifecycle.R`; `R/app_server.R:982` (`nf_start`); `R/reference_nf.R:211`.
- **Problem:** Stop, closing the tab, or reloading the page sends SIGKILL to the Nextflow JVM only. Nextflow cannot clean up, so GetOrganelle, BLAST, and R-step containers keep running and writing to the project database, and the reloaded page says nothing. Nothing stops two launches at once (a double click, or two browser windows). The docs say you can "come back later". Reference count and build runs share the project launch folder, so they break `-resume` for the main pipeline.
- **Why it matters:** a 16-minute run or a 500-sample batch is lost to one accidental reload, and two sets of writers can race on the same samples.
- **Fix:** send SIGINT first (`p$interrupt()`), wait, then `kill_tree()`, and do not kill from `onSessionEnded` inside the lifecycle grace period. Better, launch Nextflow detached and re-attach on load with a "Run in progress" banner. Add a project-level run lock and a `beforeunload` warning. Launch reference runs from their own folder. Do this once in a shared run helper (ponytail P1) so the Update and Reference windows do not drift apart. Until then, the docs should say to keep the tab open.
- **Read more:** [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-M1, APP-M9, APP-L11), [CODE 4.2](CODE_REVIEW.md#42-pipeline-and-containers) (PIPE-M2), [CODE 4.5](CODE_REVIEW.md#45-reference-data-and-taxonomy) (REF-L17), [UX-43](USABILITY_REVIEW.md#ux-43), [UX-44](USABILITY_REVIEW.md#ux-44), [UX-46](USABILITY_REVIEW.md#ux-46), [DOC-36](DOCS_REVIEW.md#doc-36).
- **Status:** fixed on this branch (second batch) the author's way: runs stop cleanly when the browser goes away or Stop is pressed, and one run per project. See [Second batch](#second-batch-actions-6-to-15).

### 14. Make locks lock, and make marks say what they do

- **Ids:** APP-M7 (C-shiny-10, C-security-10), UX-20 (U-mitopilot-1, U-mitopilot-3, U-curator-10), and DOC-90 (D-appclaims-16); APP-M15 (C-ponytail-39) and UX-25 (U-mitopilot-8); UX-33 (U-novice-3, U-mitopilot-4, U-design-14, U-curator-10 in part, O-practice-17) and DOC-41 (D-intermediate-22, D-appclaims-4, D-editor-16, D-novice-28, O-practice-17); APP-L10 (C-shiny-23); UX-21 (U-mitopilot-2); UX-23 (U-design-1, U-power-11, U-novice-14); TEST-L4 (C-tests-11).
- **Severity:** medium (APP-M7, APP-M15); major (UX-20, UX-33, DOC-41, UX-21); minor (DOC-90, UX-25, UX-23); low (APP-L10, TEST-L4). Confirmed, except C-shiny-23 (unverified-low).
- **Where:** `R/app_server.R:2089-2289` (edit, trim, restore), `:700-728` and `:2342-2350` (two lock rules), `:1931` (Problematic tooltip), `:751` (`toolbar_flag`); `R/export.R:542-549`; `R/app_opts.R:333-337, 419-423`.
- **Problem:** an annotation-locked, already exported unit can still have its boundaries edited, be trimmed, and be restored, and an unlocked assembly stays in the Export table. Toolbar Lock offers "Lock anyway" for a second path while Lock & Close refuses the same state. "Problematic" (tooltip: "Flag this annotation for another look") silently removes the unit from every export. The options guard names the wrong lock and sends users round a loop across two tabs. Toolbar toggles on a mixed selection clear every row.
- **Why it matters:** signed-off records can change, records can vanish from a submission, and users cannot predict what a button will do.
- **Fix:** make a locked unit read only with a lock banner, and refuse edits on the server with `need_unlocked()`; require both locks for export; pick one rule for a second path; rename Problematic "Exclude from export" (or say so) and name excluded units in the Export window; skip locked units in `toolbar_flag()`; name the right lock and tab in the guard; relabel toggles to show what they will do. Add export tests for Problematic and unlocked units.
- **Read more:** [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-M7, APP-M15, APP-L10), [CODE 6](CODE_REVIEW.md#6-test-suite-results) (TEST-L4), [UX-20](USABILITY_REVIEW.md#ux-20), [UX-21](USABILITY_REVIEW.md#ux-21), [UX-23](USABILITY_REVIEW.md#ux-23), [UX-25](USABILITY_REVIEW.md#ux-25), [UX-33](USABILITY_REVIEW.md#ux-33), [DOC-41](DOCS_REVIEW.md#doc-41), [DOC-90](DOCS_REVIEW.md#doc-90).
- **Status:** fixed on this branch (second batch), following MitoPilot: a second path is refused. See [Second batch](#second-batch-actions-6-to-15).

### 15. Make the chimera label mean something before it goes to GenBank

- **Ids:** SCI-M1 (C-science-1); SCI-M2 (C-science-2), UX-83 (O-practice-28, U-curator-11), and DOC-61 (D-domain-18); SCI-M3 (C-science-11); SCI-M4 (C-science-12), UX-82 (O-practice-13, U-curator-6), and DOC-57 (D-domain-11); DOC-56 (D-domain-10, D-editor-13 in part, D-novice-15 in part, D-intermediate-19 in part); APP-M11 (C-shiny-18).
- **Severity:** medium (SCI-M1 to SCI-M4, APP-M11); major (UX-82, DOC-56, DOC-57); minor (UX-83, DOC-61). All confirmed.
- **Where:** `R/chimera.R:68, 99-100, 219, 237`; `inst/nextflow/modules/coverage.nf:47-48`; `submission_risk` in `R/export.R`; the panel text at `R/app_modals.R:104-108`.
- **Problem:** 15 of 22 clean fish units are "suspect", while the planted fish plus *Penicillium* chimera is also only "suspect". The soft-clip channel can never fire (bowtie2 runs end to end), depth dips from IGS repeats and scaffold ends look like junctions, an ancestor lineage (Chordata against Pleuronectiformes) counts as a conflict, and ties in conserved 18S are broken arbitrarily, so a flounder was called fish plus dog. "Suspect" goes into `submission_risk` as "chimeric", and the bundle README warns that NCBI may remove chimeric sequences. After a trim, the junction hint keeps the old coordinates.
- **Why it matters:** it pushes users to drop good records and teaches them to ignore the flag.
- **Fix:** count clips from a `--local` alignment; compute depth steps without the mapping-quality filter (or mask repeats), skip scaffold ends, and require a minimum depth; treat an ancestor as agreement; take the LCA of hits tied at the top bitscore; re-tune the cut-offs on CHIM01; write "chimeric" to `submission_risk` only when there is a lineage switch; shift junctions on trim; document the four channels and their known false positives.
- **Read more:** [CODE 4.3](CODE_REVIEW.md#43-science-and-algorithms) (SCI-M1 to SCI-M4), [CODE 4.1](CODE_REVIEW.md#41-shiny-app) (APP-M11), [UX-82](USABILITY_REVIEW.md#ux-82), [UX-83](USABILITY_REVIEW.md#ux-83), [DOC-56](DOCS_REVIEW.md#doc-56), [DOC-57](DOCS_REVIEW.md#doc-57), [DOC-61](DOCS_REVIEW.md#doc-61).
- **Status:** fixed on this branch (second batch) with clean, minor, and severe levels. See [Second batch](#second-batch-actions-6-to-15).

Just below the top 15: validating names and paths (DB-M2, DB-M3, PIPE-M13, UX-35; most of it is in the quick wins), the variant check running before annotation finishes (SCI-M7, UX-85), the out-of-date Curation guide (DOC-53 to DOC-60), the function reference (DOC-66, DOC-72, PKG-M1), bulk actions reaching hidden rows and blank filtered pages at scale (UX-58, UX-57, APP-M12), and header templates run as R code (SEC-M1).

## Quick wins

Each should take under an hour (the reviewers' estimate; check the linked write-up before starting). Status is confirmed unless noted.

| # | Ids | What to change | Where |
|---|---|---|---|
| 1 | PKG-H1 (C-pkg-1); DOC-17 | Add an empty `biocViews:` line so `install_github()` finds the Bioconductor packages. | `DESCRIPTION` |
| 2 | PKG-H2 (C-data-3, C-pkg-2); DOC-18; UX-01 | Replace the RStudio test with `requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()`, and fix the Rproj help text. | `R/init_project.R:98-105`, `R/init_test_project.R` |
| 3 | PIPE-M11 (C-nextflow-20; C-pkg-30 unverified-low); DOC-16; UX-04 | Default `executor` to "docker" in `new_project()` and `new_test_project()`, and warn once on conda. | `R/init_project.R:52`, `R/init_test_project.R:164` |
| 4 | PIPE-H3 (C-nextflow-3, C-data-10, C-pkg-9); DOC-29 | Drop the `${params.rawDir}` bind from generated profiles, and treat 'NA' as empty in the reference command. | `R/generate_config.R:32-38`, `R/reference_nf.R:159-160` |
| 5 | PIPE-H2 (C-nextflow-7); DOC-03 | Add `process.resourceLimits` (host CPUs and memory) to the local, docker, and conda templates. | `inst/config.local`, `inst/config.docker`, `inst/config.conda` |
| 6 | DOC-20; UX-09 | Add `taxonomy_lock = 1` to the pipeline's sample query and to the Update count. | `inst/nextflow/modules/pipeline_r.nf:21-40`, `R/nf_command.R:163-176` |
| 7 | APP-H1 (C-shiny-1) | Add `ORDER BY x.path, x.scaffold LIMIT 1` to the options lookup. | `R/app_opts.R:384` |
| 8 | SCI-H3 (C-science-10) | Read the path number from the file name, not the glob order. | `R/select_paths.R:123-124` |
| 9 | SCI-M8 (C-science-9) | Read coverage with `colClasses = c("character", "integer", "integer")`. | `R/coverage.R:65` |
| 10 | DB-L5 (C-data-15); UX-03 (in part) | Read the mapping file with `colClasses = "character"` so IDs like 001 stay as typed. | `R/init_db.R:71` |
| 11 | PIPE-H4 (C-nextflow-2) | Stage scan outputs as `path(gffs, stageAs: 'chunk_*.gff')`. | `inst/nextflow/modules/find_rrna.nf:59` |
| 12 | PIPE-M5 (C-nextflow-10) | Omit `--max-paths-num` when max paths is 0 or less. | `inst/nextflow/modules/pipeline_r.nf:26` |
| 13 | PIPE-M17 (C-tests-1) | Call `RiboPilot:::record_blastdb_provenance` and log the error instead of swallowing it. This clears the one `R CMD check` error. | `inst/nextflow/bin/rp_searched.R:22` |
| 14 | APP-M2 (C-shiny-3, C-nextflow-12); UX-42 | First step: read `p$get_exit_status()` and show an error alert when it is not 0. (Per-sample counts can follow.) | `R/app_server.R:1009-1018`, `R/app_reference.R:962-968` |
| 15 | UX-45 (U-power-10, U-mitopilot-10) | Refuse to start a run with nothing queued. | Update handler, `R/app_server.R` |
| 16 | APP-M5 (C-shiny-7) | Pass `page = getReactableState("tbl", "page")` and re-select rows after each refresh. | `R/app_server.R:421-423` |
| 17 | APP-M6 (C-shiny-8); UX-59 | Observe only the `#tbl` subtree and ignore the probe's own mutations. | `inst/app/www/custom.js:12-36` |
| 18 | DB-M3 (C-data-9, C-security-4) | `shQuote()` the docker and singularity mount and working-folder arguments. | `R/itsx_check.R:65, 71-72` |
| 19 | DB-M2 (C-data-7, C-security-1, C-pkg-8); UX-35 (in part) | Validate export group names when assigned (letters, digits, dot, dash, underscore; never `.` or `..`), and check the folder is inside `out/export` before deleting it. | `R/app_server.R:778-781`, `R/export.R:571, 610-611` |
| 20 | PIPE-M13 (C-security-5); UX-35 (in part) | Refuse whitespace and quotes in project, data, assembly, and reference paths with a clear message. | `R/init_project.R`, `validate_ref_dir()` |
| 21 | DB-L2 (C-data-11); DB-L3 (C-data-12) | `unlink()` the old zip before zipping, normalize the path before `setwd()`, and check the zip status. | `R/export.R:518, 726-731` |
| 22 | REF-M5 (C-refdata-10) | Convert `pdb\|X\|Y` to `X_Y` when parsing deflines. | `R/blast_ref_fetch.R:53-55` |
| 23 | REF-M1 (C-refdata-4, C-security-9, C-ponytail-34) | Pass the stored `clade_taxid` to `clade_query()` and `--taxid` to `download_rrna.sh` (the plumbing exists). About an hour. | `R/reference_counts.R:17-22` and callers |
| 24 | UX-07 (O-practice-21); DOC-31 | Default `ncbi_api_key` to `Sys.getenv("NCBI_API_KEY")`, and comment out the placeholder key in the guide. | `R/init_project.R:181`, `vignettes/RiboPilot.Rmd:195` |
| 25 | UX-48 (O-practice-31, U-novice-21, U-mitopilot-19, U-power-21 in part) | Set `NXF_DISABLE_CHECK_LATEST=true` for app runs so the progress box does not invite users to upgrade past the pin. | run environment |
| 26 | UX-52; DOC-89 | Point the Help link at the docs site; say "samples" in the Reference Samples view (also UX-78, APP-L8, C-shiny-21 unverified-low). | `R/app_server.R:84-85, 91` |
| 27 | UX-79; DOC-69 (in part) | Print the tunnel hint only when headless; otherwise print the local URL. | `R/app_run.R:194` |
| 28 | SEC-L3 (C-security-11, C-ponytail-36; C-pkg-15 unverified-low) | Add `.claude`, `.superpowers`, `.nextflow.log*`, `update-version.log`, and `RiboPilot_test*.R` to `.Rbuildignore`. (Keeping them out of the image layer is a slightly bigger change.) | `.Rbuildignore` |
| 29 | PKG-L3 (C-pkg-18, unverified-low); PKG-L11 (C-tests-19, C-pkg-16, unverified-low) | Add rlang to Suggests; use `.data$y` and drop five stale `globalVariables` names. | `DESCRIPTION`, `R/app_modals.R:853`, `R/zzz.R:57-60` |
| 30 | DOC-27, DOC-22; PKG-L8 (C-pkg-25, unverified-low) | Knit README with `md_extensions: -smart`, and move the README.md-only paragraph into README.Rmd. | `README.Rmd` |
| 31 | DOC-30, DOC-33, DOC-52 | Use `RiboPilot::RiboPilot("~/my_project")` and define `project`; say "ray-finned fishes"; say the Reference step is optional only for ray-finned fishes. | `vignettes/RiboPilot.Rmd:143-144, 291, 391`, `vignettes/Reference-Databases.Rmd:91-94` |
| 32 | DOC-95; DOC-84 and UX-81 | Add the site URL to DESCRIPTION and re-run `build_favicons()`; reword the `.config` comment so it does not contain the `<<EXTRA_MOUNTS>>` token. | `DESCRIPTION:12`, `inst/config.docker:11` |
| 33 | APP-L3 (C-shiny-14, C-security-15) | Rebuild the cluster job script when the Resume box changes. | `R/app_server.R:888, 898-899` |
| 34 | REF-L11 (C-refdata-21, unverified-low); PIPE-L5 (C-security-23) | Treat only 1, true, or yes as opting out of local BLAST; send bowtie2 baiting SAM output to `/dev/null`. | `R/blast_local.R:88`, `R/assemble_fallback.R:32` |

## Release blockers

### Before publishing the docs site

1. **Get started describes a different app.** DOC-28 and UX-49 (action 6). The site's main walkthrough must match the five-step app.
2. **The first commands on the site must work.** DOC-17 / PKG-H1, DOC-18 / PKG-H2 / UX-01, and DOC-16 / PIPE-M11 (quick wins 1 to 3). Fixed on this branch (action 4).
3. **The repo must be public and the fish reference must load.** DOC-19 and UX-11 (action 5). The fish-reference half is deferred by the author; the repo-public half still applies. Until then the site URL, every Source link, BugReports, the in-app Help link, and the practice project all fail. Note that GitHub Pages from a private repo needs a paid plan and would be public anyway (DOC-96).
4. **"Already have an assembly?" must work on generated cluster profiles.** DOC-29 / PIPE-H3 (quick win 4).
5. **Do not publish wrong science or GenBank guidance.** Either fix the code first or document what ships: the partial rule (DOC-42, DB-H2, UX-32), isomers that ignore the ITS (DOC-58, SCI-L5), "cannot be contamination" (DOC-59), the "GenBank screen" that is really a local panel by default (DOC-60), source modifiers and the organism name (DOC-44, DOC-45), and the out-of-date Curation guide (DOC-53 to DOC-56: missing flags, 2 of 5 rulesets, 7 of 8 triage rules, no chimera section).
   DOC-60 and DOC-53 are fixed on this branch (fourth batch).
6. **Write the missing essentials.** Install prerequisites, supported computers, and hardware (DOC-01, DOC-02, DOC-03); Troubleshooting (DOC-08); mapping-file rules (DOC-09); the "keep the tab open" warning (DOC-36); and move all conda material off the main path (DOC-16). The rest of the missing-content list (DOC-04 to DOC-07, DOC-10 to DOC-14) can follow soon after.
7. **Fix the function reference.** Group it and hide internals (DOC-66, P10), turn on roxygen markdown (DOC-72, PKG-M1), and correct the wrong pages (DOC-67, DOC-68, DOC-69).
   Roxygen markdown (PKG-M1) is on as of the fourth batch.
8. **Make the site build and deploy cleanly.** The DESCRIPTION URL and favicons (DOC-95), the deploy trigger that never fires on a release (DOC-96, TEST-L1), and curly quotes that break copy and paste (DOC-94, DOC-27).
   The DESCRIPTION URL (DOC-95) is fixed on this branch; the deploy trigger is open.
9. **Retake the screenshots.** Use a fresh practice project at a neutral path after the fixes, and never publish an image that shows the NCBI API key (DOC-91, UX-47, SEC-L1). Some operator screenshots in the review workspace show a real key; they were kept out of this folder and will go when the workspace is deleted.
   Retaken in the fourth batch from a fresh practice project at a neutral path, with no key set.

### Before making the repo public or taking outside users

1. **No database schema migration path (known gap).** There is no way to move an existing project's `.sqlite` to a newer schema. This was deferred on purpose while RiboPilot is pre-release, and that is fine for now, but it becomes a blocker as soon as anyone outside the lab has projects: any change to a table or to what a column means will then break their projects, possibly silently. PIPE-M15 shows the risk is already real: after a package upgrade, the app writes with the new code while every pipeline task runs the older image code against the same database, and nothing compares the two versions (DOC-07 notes the docs are silent on upgrades too). Before release: store a schema version in each project database; check it whenever the app or a pipeline task opens a project, and stop with a clear message on a mismatch; write forward migrations from the first public version on; check the image version against the package version (PIPE-M15); and document the upgrade rule (DOC-07). Do the schema clean-up first (P13 dead schema, P14 back-compat shims, the new notes column in APP-M14, and the gate column suggested for UX-17), so the first public schema is one you want to keep.
   The image tag check (PIPE-M15, warn only) and the upgrade docs (DOC-07) are done on this branch (fourth batch); the schema version and migrations are open.
2. **Fix the wrong-record problems.** Actions 1, 2, 3, 7, and 11 (UX-31; DB-H2 and UX-32; DB-H1, PIPE-H1, SCI-H1, and SCI-M5; REF-H1, REF-H2, and REF-H3; DB-M1, DB-M5, and UX-36). Outside users will submit what RiboPilot writes.
   UX-31 (at export), DB-H2, UX-32, DB-H1, and the re-assembly half of SCI-H1 are fixed on this branch; the rest of this item is open.
3. **Close the security and privacy gaps.** SEC-M1 (header templates are run as R code, including templates stored in the project database, and the headless app has no login on shared nodes); SEC-L1 and SEC-L2 (the NCBI API key in process lists, logs, `.command.run` files, job scripts, and the run window); SEC-L4 (GenBank titles inserted as raw HTML). SEC-L3 matters already: the Docker image on Docker Hub is built from the working tree, so it ships `.claude` and `.superpowers` notes and local paths today. Also move the BUG and HANDOFF notes and personal driver scripts off the repo root (DOC-97).
4. **Make results reproducible and citable.** Pin every tool in the Dockerfile and match `environment.yml` (PIPE-M12, DOC-26); add NEWS, CITATION, and a per-tag tool-version table (DOC-12); pin the default seeds with a checksum (REF-L4).
   PIPE-M12 and DOC-12 are fixed on this branch (fourth batch); REF-L4 is open.
5. **Let CI catch what only shows up off the author's machine.** Add a CI job with Nextflow and the pipeline tools (TEST-L1), so bugs like PIPE-M17 fail CI, plus a clean-machine install and `new_test_project()` from Rscript (theme 4).
   The fourth batch adds a clean-install job with Nextflow that runs `main.nf` on the practice samples with stand-in tools; it has not run on GitHub yet, and a job with the real pipeline tools is still open.
6. **Offer only what works.** Remove awsbatch until it runs (PIPE-M10), stop defaulting to conda (PIPE-M11), remove or warn on the barrnap choice that loses the truncation flag (SCI-L2), and delete or wire the Curation options that do nothing (APP-M13, DOC-87).
7. **Remove the private-preview workarounds** (`GITHUB_PAT` notes and local seed paths) once the repo and assets are public.

## Cross-track themes

These are root problems that show up in more than one track. Fixing the root usually closes items in all three.

### 1. The app has the evidence but does not act on it or say it

- **Code:** `risk_flags()` ignores the screen verdict; export ignores model coverage and the truncated flag (DB-H2); the run's exit status is never read (APP-M2); two Curation options are saved but never read (APP-M13).
- **Usability:** contaminated and chimeric units show a green check (UX-14, UX-22), pass lock and export with no prompt (UX-31), and are not named in the export dialog (UX-34); failed runs say Done (UX-42); failed samples give no reason (UX-15); Problematic removes units silently (UX-33).
- **Docs:** flags users see are missing from the Curation guide (DOC-53); the Problematic rule and the partial rule are undocumented (DOC-41, DOC-42).
- **Root and fix:** signals are computed and stored, but the decision points (tables, lock, export, end of run) do not read them. One per-unit "checks" summary, used by the tables, both locks, the export gate, and the run summary, would close most of these.

### 2. The second run is less safe than the first

- **Code:** old path files, old coordinates, old verdicts, lost triage decisions, and skipped searches after re-assembly (PIPE-H1, DB-H1, SCI-H1, SCI-M5, PIPE-M1, DB-L1, SCI-L10); option changes that re-assemble and undo trims (PIPE-M3); `-resume` reusing stale fallback or pick results (PIPE-M4, PIPE-M7); a finished reference build re-queuing samples (REF-M4); tests that set database state directly and so miss all of this (TEST-M1).
- **Usability:** a rescue skips its own new scaffolds and sticks "In progress" (UX-17); Save points samples at a set that was never written (UX-18); the run window does not say what will re-run (UX-24).
- **Docs:** the guide's advice on when to re-run is wrong (DOC-38), and the documented rescue re-assembles good samples (DOC-37).
- **Root and fix:** per-unit evidence is keyed on path and scaffold numbers, not on the sequence, and nothing is invalidated when inputs change. Key carry-over on a sequence hash, record what each derived result was computed from, and recompute or clear when it changes (action 3).

### 3. The console-first past and the app-first present disagree

- **Docs:** Get started runs the pipeline from the console and describes three tabs (DOC-28); console rescue and export are the documented routes (DOC-37, DOC-40, DOC-32); console messages point to `.config`, conda, and maintainer scripts (DOC-79); help pages describe older behavior (DOC-67, DOC-68).
- **Code:** `run_fallback()` re-assembles good samples (PIPE-H5); dead console-era workflows and wrappers (P2, P4, P5, P7, P20, P34); about 40 exports that exist only for scripts (P10); the pipeline ignores the taxonomy lock that only the app's table enforces (the code side of DOC-20 and UX-09); conda hints in reference help (REF-L18).
- **Usability:** the guide describes a different app (UX-49); the rescue is hard to find in the app because the docs show only the console route (UX-16); `export_project()` errors without app-made groups (UX-37).
- **Root and fix:** the scripted route predates the five-step app and was never retired. Put the guards in the pipeline (not only in the app), move scripting to one clearly labelled advanced page, and delete what the app already covers.

### 4. It works on the author's machine

- **Code:** CI installs with pak, which hides the missing `biocViews` (PKG-H1); tests and driver scripts pass `Rproj = FALSE`, which hides the RStudio crash (PKG-H2); the author's system has jq, sqlite3, and R, which hides the conda gaps (PIPE-M11); `load_all()` exposes internal functions and CI has no Nextflow, which hides PIPE-M17 (TEST-L1); the author has local seed files, which hides the 404 (SCI-M9, REF-L4).
- **Usability:** the first command fails for all three operators (UX-01); the default executor fails silently (UX-04); the practice project's first run aborts (UX-11); paths with spaces break every task (UX-35).
- **Docs:** `GITHUB_PAT` is never mentioned (DOC-19); install prerequisites, supported computers, and hardware are missing (DOC-01, DOC-02, DOC-03).
- **Root and fix:** nothing tests a fresh machine. A small CI job (fresh R, `remotes::install_github()`, `new_test_project()` from Rscript, and a `-stub-run` of `main.nf`) would have caught most of these before a user did.

### 5. Free text flows unchecked into folder names, shell commands, and R code

- **Code:** group names become deleted folders (DB-M2); unquoted mounts break on spaces (DB-M3, PIPE-M13); header templates are run as R code (SEC-M1); the API key is pasted into every task's shell (SEC-L2); config values go into Groovy unescaped (PIPE-L3, PIPE-L4); IDs are read as numbers (DB-L5); reserved column names clash (DB-L6); titles are inserted as raw HTML (SEC-L4).
- **Usability:** spaces and slashes break runs and export with Docker errors (UX-35); mapping mistakes give R internals (UX-55, UX-03).
- **Docs:** the mapping-file rules are not written down (DOC-09).
- **Root and fix:** validate at the point of entry (names, paths, IDs, keys), quote at the point of use, and never evaluate user templates.

### 6. Locks, marks, and names mean different things in different places

- **Code:** locked units can be edited (APP-M7); two opposite rules for a second path (APP-M15); toolbar marks change locked units (APP-L10).
- **Usability:** locks do not protect (UX-20); the options guard names the wrong lock (UX-21); toggles flip the whole selection (UX-23); Lock & Close marks Reviewed but toolbar Lock does not (UX-25); Problematic excludes from export (UX-33); one scaffold has five names (UX-77).
- **Docs:** "Locking a sample signs off on it" is not true (DOC-90); the docs and app use different names for the same things (DOC-15); tooltips describe old behavior (DOC-88).
- **Root and fix:** one rule per concept (a locked unit is read only; one rule for a second path; Problematic says it excludes), and one vocabulary taken from the toolbar, used in the app and the docs.

### 7. The taxon is not carried through consistently

- **Code:** a corrected taxon leaves the old lineage (REF-H1); homonyms pick the lowest taxid (REF-H3); clade databases are fetched by name (REF-M1, REF-L2); a partial taxdump is trusted (REF-H2); export uses the raw Taxon (DB-M1); locked and unresolved taxa can be changed or locked (APP-L1, APP-L16).
- **Usability:** unresolved samples are counted as assigned, and "Morus" goes to a plant (UX-08); on_target is green even at phylum level (UX-86).
- **Docs:** `[organism=]` uses the raw label (DOC-45); "cannot be contamination" overstates a phylum-level screen (DOC-59); the taxonomy download size and time are understated (DOC-24).
- **Root and fix:** names and taxids are used interchangeably, and a correction updates some columns but not others. Carry the taxid (with its lineage) as the one identity from the Taxonomy step to the screen, the reference groups, and export.

### 8. Runs belong to a browser tab

- **Code:** Stop, tab close, and reload SIGKILL Nextflow and leave containers running (APP-M1); double launches are possible (APP-M9, APP-L11); reference builds are not detached (REF-L17); reference runs break `-resume` (PIPE-M2).
- **Usability:** a reload silently kills a run (UX-43); Stop leaves containers running (UX-44); two sessions can run at once (UX-46).
- **Docs:** the docs say you can keep working or come back later (DOC-36); cluster users who tunnel lose the run when the tunnel drops (DOC-04).
- **Root and fix:** the pipeline process is owned by the Shiny session. Launch it detached with a project-level run lock, and let the app re-attach (action 13).

### 9. Science signals are noisy or run at the wrong time

- **Code:** the chimera channels misfire (SCI-M1 to SCI-M4); the variant check runs before annotation finishes (SCI-M7); isomers are judged on the anchors only (SCI-L5); ITS spacers do not touch the genes (DB-M5); path and scaffold numbering errors (SCI-H3, SCI-M8).
- **Usability:** a flounder is called fish plus dog (UX-82); depth and scaffold ends raise chimera flags (UX-83); the variant check is usually blank (UX-85); 1 to 4 bp gaps give table2asn warnings (UX-84).
- **Docs:** the chimera false positive (DOC-57); chimera detection is undocumented (DOC-56); isomers are said not to matter (DOC-58); the depth explanation does not match the code (DOC-61); spacers are said to fill the gaps (DOC-39).
- **Root and fix:** the signals were designed well, but were not tuned against the clean practice data and the known chimera. Re-tune them on CHIM00, CHIM01, and the clean fish, and run the variant check after annotation (action 15).

## Limits of this review

- **No real users.** Every reviewer, persona, operator, and verifier was a review agent working from a written role. The SUS scores are the personas' estimates, not questionnaire results.
- **One platform.** Everything ran on one Linux workstation. macOS (Intel or Apple Silicon) and Windows or WSL2 were not tested, and the app was launched from R scripts and driven through headless Chrome rather than from RStudio.
- **No real cluster.** No SLURM, SGE, PBS, or LSF job was submitted. Singularity was tested only on the local executor, with a hand-edited profile (docs-hpc). Generated profiles and job scripts were checked by reading, config parsing, and small Nextflow runs. The NMNH Hydra and NOAA SEDNA profiles and AWS Batch were not run.
- **Supported route only, one Nextflow version.** Runs used the Docker image and the pinned Nextflow 25.10.6. The conda route was only seen failing (op-edge) and through the test suite's tool-dependent tests. Nextflow 24.10.x was not tried.
- **Small, fish-heavy data.** Reads were filtered subsets (the practice project) or the first 3 million pairs (two gastropods). The invertebrate rulesets (asteroidea, scyphozoa, octocorallia) were not exercised on real data. User-supplied assemblies were synthetic. No full-depth skims and no large pipeline batch were run; the 504-sample project was an inflated database used only to time the app. Remote NCBI BLAST was barely exercised: no practice unit was searched remotely.
- **No submission.** Exports were checked with table2asn 1.29 and against NCBI's help pages. Nothing was sent to NCBI, so the Portal's own ITSx and length checks were not seen.
- **No upgrade path.** No project from an older RiboPilot was opened or migrated (by the no-migrations rule). The package and image version mismatch (PIPE-M15) was found by reading.
- **Security by reading.** Security findings come from reading and small local reproductions. There was no penetration test, and exposure on a shared login node was reasoned about, not tested with a second user.
- **Accessibility by script.** WCAG 2.2 AA checks used headless Chrome (the accessibility tree, real key events, and contrast math). No screen reader or other assistive technology was used.
- **Partial verification in code.** 117 low and nit code findings were not verified, to save usage. The ponytail line counts are estimates. Coverage was measured by tracing functions (covr was not installed), and lintr ran against the installed older RiboPilot, so its "no visible function" hits were ignored.
- **Docs not deployed.** The site was built on a copy. Links into the private repo return 404, so they could not be checked as public pages. The docs-novice and docs-intermediate reviews come from the first attempt and were not re-run. The screenshots come from review folders and show known bugs.
- **Interrupted run.** The review was stopped once by a usage limit and resumed. Some reviewers worked in two or three sittings, and the leads reused earlier drafts after re-checking them.
- **MitoPilot comparisons** were made by reading MitoPilot's source, not by running it side by side.

## Completeness check: gaps and caveats

A final checker read all four reports after they were written. Its notes are below. They were not fed back into a new review round.

### Correct these before acting on them

- **SCI-H3 (path scores on the wrong path) was shown on an older project.** Its real-data proof ("SRR14433904 has exactly 10 paths") comes from the author's `RiboPilot_test/2026_09_02` project. In this review's 0.0.0.9018 practice run, every sample has 1 path, so action 10's "one practice sample has exactly 10" does not hold for the current code. SCI-H2 (circular units not rotated) was shown only on a hand-rotated unit. Re-check both on a sample that gives several paths.
- **Chimera noise rates come from the older project.** SCI-M1 to SCI-M4, action 15, and theme 9 quote "15 of 22 clean fish units are suspect". The 0.0.0.9018 practice project gives 3 of 14 fish units suspect (plus CHIM01). The problem is real, but it is about a third as bad as stated. Other statistics taken from `2026_09_02` may be old too.
- **APP-H1 (options cell ends the session) was reproduced only with `shiny::testServer`.** Nobody hit it in a browser, even though the practice project has 2 multi-unit samples (SRR33849683 and SRR22396893). Click those cells in the app before you fix it.
- **Some top actions and quick wins rest on unverified findings.** Quick wins 29, 30, and 34 are unverified-low. Actions 3, 5, 7, 13, and 14 cite unverified ids: C-science-24, C-pkg-20, C-refdata-17, C-shiny-31, C-shiny-26, C-refdata-27, C-shiny-32, and C-shiny-23. The ribotyper part of PIPE-M11 (C-pkg-30) is also unverified. 24 of the 57 ponytail ids and the 2,600-line estimate were not checked. Unexporting functions (P10) would break bin scripts that call `RiboPilot::` unless those scripts change too.
- **Export findings were judged with table2asn.** NCBI's Submission Portal takes rRNA-ITS as FASTA only and runs its own ITSx. So `.tbl`-only findings (DB-M5 gaps, partial markers in the feature table, linkage evidence) matter only for the table2asn and BankIt route. Defline, flank, organism, and isolate findings matter for both. Tag each export finding by route and re-rank.
- **The author's known-issue list was not matched.** `dev/plans/2026-09-15-ui-review-findings.md` (37 items), the BUG and HANDOFF notes, `dev/design/HANDOFF.md`, and `dev/design/05_ROADMAP.md` were not cross-referenced. Some items were found again (for example, prior item 4 is now UX-23), and others may have been missed or already fixed.

### Journeys not exercised

- Choosing a Selected Path, locking a second path, and exporting several paths (no sample had more than 1 path on current code).
- Running from an installed package (`install_github` or a built tarball) instead of `pkgload::load_all()`.
- The headless cluster route in the app (Save job script, Submit as job, Work Dirs after a cluster run, RStudio Server detection), and parsing the sge, pbs, and lsf profiles made by `generate_config()`.
- Remote NCBI BLAST, including NCBI throttling and losing the connection partway through a run.
- The de novo SPAdes fallback (every rescue used baited SPAdes).
- Real user-supplied genome assemblies (`.fna.gz` with thousands of contigs and tandem rDNA arrays). Only synthetic ones were tried.
- A large reference build (for example a class over the 25,000-record cap), reuse of the shared reference store by a second project, and two builds at once.
- The invertebrate rulesets on real data, a mixed-clade project, and non-animal samples.
- Full-depth reads and a concurrent batch of 30 to 50 samples. The hardware advice and the database-lock risk have no measurements behind them.
- The fish reference served by URL. Every successful run used the local repo copy, so the "repo public" release gate is untested.
- Nextflow 24.10.x, the low end of the stated range.
- Accessibility and 500-sample checks on the Taxonomy and Reference steps, and moving a project folder to a new location.

## Appendix A: refuted findings

Listed only for the record. The track reports give the reasons.

| Id | Track | Title |
|---|---|---|
| C-data-8 | Code | An ITSx check that could not run overwrites the pipeline's stored ITSx warnings with `not_run(...)` (intended and tested) |
| C-data-16 | Code | Zero-feature and undecided-triage units are exported with no submission-risk flag (the final ITSx check flags them) |
| C-tests-12 | Code | The scientific-defaults test checks its own copy of the reference filter (another test runs the real one) |
| D-appclaims-20 | Docs | Re-annotating a unit is undocumented and needs the assembly unlocked first (the State button can re-queue it; the vague hint is kept as UX-28) |
| D-domain-15 (part) | Docs | README.txt describes the ITSx "chimeric" flag wrongly (ITSx 1.1.3 does mark it "Chimeric!"; the rest of the finding stands in DOC-48) |
| U-mitopilot-9 | Usability | Toolbar flag buttons do not say what they will do, unlike MitoPilot (MitoPilot's toolbar labels are static too; the direction problem is kept as UX-23) |
| U-mitopilot-22 | Usability | The Triage dialog puts Apply left of Cancel, unlike other dialogs (all destructive dialogs share that MitoPilot convention; the rest is kept as UX-64) |

## Appendix B: every item, by track

Every finding that was not refuted is listed here, merged into the items the track reports use, with all of its source ids. Nothing was dropped. Each item's location, why it matters, and the fix are in the linked track report. The last column links items in different tracks that describe the same problem, so a fix can close several at once. Status is **confirmed** for every docs and usability item. For code items, the status column names any id whose status differs from the rest (unverified-low or uncertain).

### B.1 Code review items

Severity scale: high, medium, low, nit. Status is per finding; where the ids in one item differ, the odd ones are named. Links go to the area section of CODE_REVIEW.md, where each item has its location, why it matters, and the fix.


**Shiny app** ([CODE_REVIEW.md section 4.1](CODE_REVIEW.md#41-shiny-app))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| APP-H1 | high | confirmed | C-shiny-1 | The Annotation or Curation options cell crashes the session for any multi-unit sample | - |
| APP-H2 | high | confirmed | C-shiny-5 | A second trim reuses the old span and length, so a flank trim can delete 18S and ITS1 | - |
| APP-M1 | medium | confirmed | C-shiny-6, C-nextflow-8, C-security-6, C-pkg-14 | Stop, tab close, and page refresh SIGKILL Nextflow, leaving task containers running | DOC-36, UX-43, UX-44 |
| APP-M2 | medium | confirmed | C-shiny-3, C-nextflow-12 | Failed Nextflow runs are announced as a success | DOC-88, UX-42 |
| APP-M3 | medium | confirmed | C-shiny-2 | "Run from App" crashes the session when Nextflow is missing or the command cannot be built | - |
| APP-M4 | medium | confirmed | C-shiny-4 | Picking another existing parameter set while Edit is ticked overwrites that set with the values on screen | - |
| APP-M5 | medium | confirmed | C-shiny-7 | Every refresh resets the main table to page 1 and clears the selection | - |
| APP-M6 | medium | confirmed | C-shiny-8 | `custom.js` re-triggers itself about 20 times a second while idle | UX-59 |
| APP-M7 | medium | confirmed | C-shiny-10, C-security-10 | Annotation-locked units can still be edited, trimmed, and restored | DOC-90, UX-20 |
| APP-M8 | medium | confirmed | C-shiny-13 | A new set name saved with Edit unticked points samples at a set that does not exist | UX-18 |
| APP-M9 | medium | confirmed | C-shiny-15 | Nothing stops Nextflow from being launched twice | UX-46 |
| APP-M10 | medium | confirmed | C-shiny-16, C-pkg-10 | The Work Dirs browser misses every cluster run | UX-54 |
| APP-M11 | medium | confirmed | C-shiny-18 | The chimera junction hint and track keep pre-trim coordinates | - |
| APP-M12 | medium | confirmed | C-shiny-19 | Bulk triage can silently switch to the whole table, including locked samples | UX-58 |
| APP-M13 | medium | confirmed | C-science-17, C-ponytail-2, C-pkg-3 | The Curation Options "Min anchor identity" and "Flag contaminants" controls do nothing | DOC-87 |
| APP-M14 | medium | confirmed | C-ponytail-1 | Curation overwrites the curator's typed Notes | - |
| APP-M15 | medium | confirmed | C-ponytail-39 | Two opposite rules for locking a second assembly path | DOC-90, UX-25 |
| APP-L1 | low | confirmed | C-shiny-11 | Clicking a locked sample in Taxonomy mode opens the review and lets the user assign a new taxon; the taxonomy lock is not checked | - |
| APP-L2 | low | confirmed | C-shiny-12 | Closing the unsaved-edits prompt with its X leaves the edit dirty; re-ticking Edit then overwrites the backup with the nudged coordinates | - |
| APP-L3 | low | confirmed | C-shiny-14, C-security-15 | The cluster job script is built once with `resume = TRUE` and never rebuilt when the "Resume previous run?" box changes; the preview does follow the box | - |
| APP-L4 | low | confirmed | C-shiny-17 | Once the Export modal has been opened, every refresh re-runs about 4 queries per exportable unit (including full sequences) to validate the header | - |
| APP-L5 | low | confirmed | C-ponytail-4 | Opening BLAST Options can run a synchronous `docker run` (60 s timeout) to describe the panel, and an empty answer is never cached | DOC-04 |
| APP-L6 | low | confirmed | C-ponytail-40, C-data-27 | `set_export_opts()` turns a failed save into a warning, and the app then always says "Saved header template" | - |
| APP-L7 | low | unverified-low | C-shiny-20 | The curation warnings summary depends only on the selected unit, so it is stale after Save or Trim | - |
| APP-L8 | low | unverified-low | C-shiny-21 | The status line says "groups" in the Reference Samples view | DOC-89, UX-78 |
| APP-L9 | low | unverified-low | C-shiny-22 | Pressing Cancel on "This export already exists" rebuilds the Export modal from saved values | - |
| APP-L10 | low | unverified-low | C-shiny-23 | The toolbar Reviewed, Problematic, and Partial buttons change annotation-locked units, while the modal disables the same toggles | UX-20 |
| APP-L11 | low | unverified-low | C-shiny-26 | "Submit as job" has no duplicate-submission guard (the reference flow has one) | UX-46 |
| APP-L12 | low | unverified-low | C-shiny-27 | `unit_by_uid` ignores the sample part of the uid | - |
| APP-L13 | low | unverified-low | C-shiny-28 | Reference settings accept a minimum larger than the maximum | - |
| APP-L14 | low | unverified-low | C-shiny-29 | Each K or D click fully re-renders the scaffold table | - |
| APP-L15 | low | unverified-low | C-shiny-30 | Triage Apply always writes the project-wide "default" preset, while the modal says it is for N samples | - |
| APP-L16 | low | unverified-low | C-shiny-31 | Taxonomy Lock accepts samples with no taxid | UX-08 |
| APP-L17 | low | unverified-low | C-shiny-33 | The progress observer reads `nf_state()` without `isolate()`, so each write re-runs it immediately and bypasses the 250 ms timer | - |

**Pipeline and containers** ([CODE_REVIEW.md section 4.2](CODE_REVIEW.md#42-pipeline-and-containers))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| PIPE-H1 | high | confirmed | C-nextflow-1, C-data-5 | Re-assembly loads old path files left in the publish folder | - |
| PIPE-H2 | high | confirmed | C-nextflow-7 | On the local executor, a memory request larger than the machine ends the whole run | DOC-03 |
| PIPE-H3 | high | confirmed | C-nextflow-3, C-data-10, C-pkg-9 | Cluster profiles made by `generate_config()` cannot run assembly-only projects | DOC-29 |
| PIPE-H4 | high | confirmed | C-nextflow-2 | A supplied assembly with more than 5,000 contigs aborts the whole run | - |
| PIPE-H5 | high | confirmed | C-pkg-5, C-nextflow-22, C-ponytail-3 | `run_fallback()` re-assembles good samples and chains preset names | DOC-37, UX-16 |
| PIPE-M1 | medium | confirmed | C-nextflow-5, C-refdata-14 | Changed BLAST settings never cause a new search | - |
| PIPE-M2 | medium | confirmed; C-shiny-32 unverified-low | C-nextflow-9, C-shiny-32 | Reference count and build runs break `-resume` for the main pipeline | - |
| PIPE-M3 | medium | confirmed | C-nextflow-4 | Changing annotation options re-assembles the sample and undoes curator trims | - |
| PIPE-M4 | medium | confirmed | C-nextflow-6 | `-resume` reuses an old fallback assembly after its options change | - |
| PIPE-M5 | medium | confirmed | C-nextflow-10 | Max paths 0 crashes GetOrganelle, and the max-paths gate can never fire | - |
| PIPE-M6 | medium | confirmed | C-nextflow-14 | A fallback SPAdes crash (including out of memory) exits 0, so the memory retry never runs | UX-15 |
| PIPE-M7 | medium | confirmed | C-nextflow-15 | A cached `find_rrna_pick` skips its database updates, so the sample stays queued | - |
| PIPE-M8 | medium | confirmed | C-nextflow-16 | A failed `find_rrna_scan` chunk silently shrinks the search | - |
| PIPE-M9 | medium | confirmed | C-nextflow-17 | Samples below minDepth are silently skipped and stay queued; the depth unit does not match the docs | UX-15 |
| PIPE-M10 | medium | confirmed | C-nextflow-18 | The awsbatch config cannot run as written | DOC-71 |
| PIPE-M11 | medium | confirmed; C-pkg-30 unverified-low | C-nextflow-20, C-pkg-30 | The default executor is conda, and `environment.yml` lacks tools the pipeline calls | DOC-16, UX-04 |
| PIPE-M12 | medium | confirmed | C-nextflow-19, C-security-22 | Tool versions are not pinned, and the conda and container versions differ | DOC-12, DOC-26 |
| PIPE-M13 | medium | confirmed | C-security-5 | A space in a project, read, or reference path breaks every container task | UX-35 |
| PIPE-M14 | medium | confirmed | C-nextflow-23, C-security-13 | Cluster head-job scripts only fit Hydra-style SGE | DOC-04 |
| PIPE-M15 | medium | confirmed | C-nextflow-25 | Container and package version mismatch goes unnoticed, and every task prints a developer warning | DOC-07, DOC-85 |
| PIPE-M16 | medium | confirmed | C-nextflow-13, C-tests-8 | With no usable reference, the bait fallback assembles the whole genome and keeps every contig | - |
| PIPE-M17 | medium | confirmed | C-tests-1 | BLAST database provenance is never recorded in container runs | - |
| PIPE-L1 | low | uncertain | C-nextflow-11 | The per-unit database writers (reconcile, curate, chimera, triage, screen) use `errorStrategy 'ignore'` with no retry, and the busy timeout is 10 s | - |
| PIPE-L2 | low | confirmed | C-nextflow-24 | Generated singularity profiles set no `pullTimeout` and no shared `cacheDir`, and every small R step inherits the 16 GB per-task memory default | DOC-04, DOC-05 |
| PIPE-L3 | low | confirmed | C-security-8 | Values put into Groovy config strings (queue, walltime, cluster options, paths) are not escaped | - |
| PIPE-L4 | low | confirmed | C-security-19 | The config reader cuts everything after `//` and strips all quotes | - |
| PIPE-L5 | low | confirmed | C-security-23 | Bowtie2 baiting writes a full SAM file to the R temp folder that nothing reads | - |
| PIPE-L6 | low | confirmed; C-nextflow-42 unverified-low | C-nextflow-42, C-security-26 | `nf_pending_count()` re-implements the pipeline's sample queries (and already differs from them), and `nf_pin_version()` starts a JVM (`nextflow -version`) on every call although the version is cached per session | DOC-20, UX-24, UX-45 |
| PIPE-L7 | low | unverified-low | C-nextflow-27 | The command-line `-process.containerOptions` for reference mounts replaces any `containerOptions` the user set in `.config` | - |
| PIPE-L8 | low | unverified-low | C-nextflow-28 | Score reads `projectDir/out`, while the assemble steps publish to `launchDir/out` | - |
| PIPE-L9 | low | unverified-low | C-nextflow-29 | `stage_assembly` runs `awk` directly on the supplied file, so a gzipped assembly becomes binary garbage | UX-05 |
| PIPE-L10 | low | unverified-low | C-nextflow-30 | Exit code 140 (a walltime kill on SGE and LSF) is in the out-of-memory retry list, and retries only raise memory | - |

**Science and algorithms** ([CODE_REVIEW.md section 4.3](CODE_REVIEW.md#43-science-and-algorithms))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| SCI-H1 | high | confirmed | C-science-3, C-tests-4 | The contamination screen verdict is never recomputed when hits or the sample's taxon change | - |
| SCI-H2 | high | confirmed | C-science-5 | Circular units are never rotated, so a start point inside 28S or 18S gives a broken model | - |
| SCI-H3 | high | confirmed | C-science-10 | Path scores are assigned by file-listing order, so with 10 or more paths they land on the wrong path | - |
| SCI-M1 | medium | confirmed | C-science-1 | The soft-clip chimera channel can never fire | - |
| SCI-M2 | medium | confirmed | C-science-2 | The depth-step channel fires on normal IGS and ETS repeats | DOC-61, UX-83 |
| SCI-M3 | medium | confirmed | C-science-11 | The gene channel counts a shallow (unresolved) element as a lineage conflict | DOC-57 |
| SCI-M4 | medium | confirmed | C-science-12 | The window channel picks an arbitrary best hit among ties in conserved regions | DOC-57, UX-82 |
| SCI-M5 | medium | confirmed | C-science-8, C-tests-2 | Re-assembly wipes triage decisions, so a curator's hand drop comes back as "keep" | UX-17 |
| SCI-M6 | medium | confirmed | C-science-6, C-tests-6 | The feature model takes the first hit per gene, so a partial leading copy replaces the real anchor | - |
| SCI-M7 | medium | confirmed | C-science-4 | The isomer and paralog classes are computed before annotation finishes | DOC-58, UX-85 |
| SCI-M8 | medium | confirmed | C-science-9 | Coverage names like "1.10" are read as the number 1.1, merging two scaffolds' depth tracks | - |
| SCI-M9 | medium | confirmed | C-science-18 | A failed seeds download stops everything after scoring for that sample | DOC-19, UX-11 |
| SCI-L1 | low | confirmed | C-science-13 | Pybarrnap reports "partial" only when the hit is under 80 percent of a fixed 28S length (2,912 bp), and otherwise `parse_barrnap` sets coverage to 100 | - |
| SCI-L2 | low | confirmed | C-science-14 | Perl barrnap 1.10.5 with kingdom `fun` has no 5.8S model and writes no partial notes, so every hit gets coverage 100 | - |
| SCI-L3 | low | confirmed | C-science-15 | ITSx runs with `--complement T` (sense-frame coordinates) while barrnap coordinates are in the input frame, and the two are combined unmapped | - |
| SCI-L4 | low | uncertain | C-science-19 | Hits titled like genome records ("chromosome, complete genome") are left out of the taxonomy LCA as well as the rDNA vote | - |
| SCI-L5 | low | confirmed | C-science-20 | The isomer test compares only the 18S and 28S anchors, but true isomers share identical ITS, while real paralogs differ mainly in ITS | DOC-58 |
| SCI-L6 | low | confirmed; C-science-27 unverified-low | C-science-27, C-pkg-12 | An explicit `metazoa_generic` ruleset is treated as "not set" and replaced by lineage auto-selection | - |
| SCI-L7 | low | unverified-low | C-science-21 | `annotate_records` treats a reference record as minus-strand if any single anchor is on minus | - |
| SCI-L8 | low | unverified-low | C-science-22 | `consensus_hit_class` votes per HSP, not per subject | - |
| SCI-L9 | low | unverified-low | C-science-23 | `scaffold_depth` drops the first and last segment of circular GetOrganelle headers | - |
| SCI-L10 | low | unverified-low | C-science-24 | On a rerun with an unchanged sequence, the seeds-database best hit overwrites the GenBank best hit in `assemblies.blast_*`, because BLAST is skipped as already searched | - |
| SCI-L11 | low | unverified-low | C-science-25 | The max-paths and max-scaffolds gate does not stop the scoring pass | - |
| SCI-L12 | low | unverified-low | C-science-26 | Find_rrna caps kept contigs at `max_scaffolds` (10) and uses the "default" annotate options, while doc 20 says a cap of 20 and the sample's own options | - |
| SCI-L13 | low | unverified-low | C-science-28 | Unanimous phylum-level discordance with 3 or more hits gives `off_target` at high confidence | - |
| SCI-L14 | low | unverified-low | C-science-29 | "odd one out" fires at species rank with contamination wording | - |
| SCI-L15 | low | unverified-low | C-science-30 | Ribotyper scoring keeps only fatal codes (so the documented LowScoreDifference weight never applies) and uses the default model set | - |
| SCI-L16 | low | unverified-low | C-science-32 | The gene-channel reason reports whole-HSP identity as the element identity | - |
| SCI-L17 | low | unverified-low | C-science-33 | When a supplied contig is reverse-complemented, no flag is stored | - |
| SCI-L18 | low | unverified-low | C-science-34 | The contaminant flag uses the single top hit while the screen uses an LCA | - |
| SCI-N1 | nit | unverified-low | C-science-31 | The chimeric triage reason includes "depth step at ...", while doc 09 says depth never appears in a reason | - |
| SCI-N2 | nit | unverified-low | C-science-35 | `run_itsx` leaves ITSx output files in the temp folder | - |

**Database and export** ([CODE_REVIEW.md section 4.4](CODE_REVIEW.md#44-database-and-export))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| DB-H1 | high | confirmed | C-data-1, C-tests-3 | Re-assembly keeps a locked unit's old gene coordinates, and export writes them onto the new sequence | UX-17 |
| DB-H2 | high | confirmed | C-data-2, C-science-7 | Truncated genes are exported as "complete sequence" with no partial markers | DOC-42, UX-32 |
| DB-M1 | medium | confirmed | C-data-4 | `[organism=]` uses the raw mapping label, not the organism confirmed in the Taxonomy tab | DOC-45 |
| DB-M2 | medium | confirmed | C-data-7, C-security-1, C-pkg-8 | Export group names are used as folder names with no checks | UX-35 |
| DB-M3 | medium | confirmed | C-data-9, C-security-4 | The export ITSx check breaks on a space in the group or project path | UX-35 |
| DB-M4 | medium | confirmed | C-data-13, C-shiny-9 | Export and pending-count code opens the database with no busy timeout | - |
| DB-M5 | medium | confirmed | C-data-17, C-science-16 | ITS spacers do not touch the rRNA genes, so table2asn warns on every ITS record | DOC-39, UX-84 |
| DB-L1 | low | confirmed | C-tests-5 | The re-assembly cleanup removes orphan rows only from `annotate`, `annotations`, and `coverage`; `screen`, `blast_hits`, `blast_screen`, `blast_searched`, `chimera`, and `chimera_windows` keep rows for scaffolds that no longer exist | - |
| DB-L2 | low | confirmed | C-data-11 | `utils::zip(flags = "-rq")` adds to an existing zip, so records removed from a group stay in it | DOC-47, UX-38 |
| DB-L3 | low | confirmed | C-data-12 | With a relative project path (the default "."), the zip path is built before `setwd()`, so zip fails, and its status is not checked | DOC-47, UX-38 |
| DB-L4 | low | confirmed | C-data-14, C-data-18, C-pkg-7 | `new_project()` creates `.sqlite` before it finishes checking the mapping and the config, and copies `mapping.csv` over the existing one before the "already exists" check | DOC-70, UX-02 |
| DB-L5 | low | confirmed | C-data-15 | `read.csv()` without `colClasses` turns all-digit IDs like `00123` into `123` | DOC-09, UX-03 |
| DB-L6 | low | confirmed | C-data-19 | Mapping columns that clash with built-in sample columns (case-insensitive in SQLite) crash init (`TaxID`) or are silently overwritten (`lineage`, `rank`) | DOC-09 |
| DB-L7 | low | confirmed | C-data-20 | Trim, restore, and several delete-then-write functions commit in several steps with no transaction | - |
| DB-L8 | low | confirmed | C-data-21 | `partial = 'yes'` marks both the first and last features partial, wherever they sit | DOC-42, UX-32 |
| DB-L9 | low | confirmed | C-data-22 | The pipeline runs ITSx with `\|\| true` and then creates an empty positions file, so a crash is recorded as `itsx_not_detected` | - |
| DB-L10 | low | confirmed | C-data-23 | The feature table ignores `annotations.strand` and always writes start then stop | - |
| DB-L11 | low | confirmed | C-data-24 | Every N-run is written as `assembly_gap` with `linkage_evidence align genus` | DOC-49 |
| DB-L12 | low | confirmed | C-data-25 | The SeqID is the sample ID when a sample has one exportable unit, else `ID_p<path>_s<scaffold>` | - |
| DB-L13 | low | confirmed | C-data-26, C-security-20 | The list of network filesystems that must not use WAL misses `ceph`, `9p`/`v9fs` (WSL2 drives), and `wekafs`; unknown types get WAL | - |
| DB-L14 | low | confirmed | C-data-28 | Blank tokens render empty source modifiers such as `[specimen_voucher=]`, and validation only dry-runs the first row | - |
| DB-L15 | low | confirmed | C-data-31 | An all-assembly mapping is still required to have R1 and R2 columns | DOC-09, UX-06 |
| DB-L16 | low | confirmed | C-pkg-11 | `export_bundle()` and `export_project()` pointed at a non-project folder create an empty `.sqlite`, and `export_project()` leaks its connection on the error | - |

**Reference data and taxonomy** ([CODE_REVIEW.md section 4.5](CODE_REVIEW.md#45-reference-data-and-taxonomy))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| REF-H1 | high | confirmed | C-refdata-1 | Changing a sample's taxon never updates `samples.lineage_taxids` | - |
| REF-H2 | high | confirmed | C-refdata-2, C-pkg-4, C-security-18 | A partial taxdump download is kept, and the index built from it is trusted forever | DOC-24 |
| REF-H3 | high | confirmed | C-refdata-3 | Homonym taxon names silently resolve to the lowest taxid, usually the wrong kingdom | UX-08 |
| REF-M1 | medium | confirmed | C-refdata-4, C-security-9, C-ponytail-34 | Clade databases are counted and downloaded by name, so homonym clades pull in another organism's rDNA | - |
| REF-M2 | medium | confirmed; C-shiny-25 unverified-low | C-refdata-5, C-shiny-25 | The shared reference folder is not safe for several users, and one error fails every running build | - |
| REF-M3 | medium | confirmed | C-refdata-6 | Clade downloads are never checked for completeness, with no retry and no limit on parallel requests | - |
| REF-M4 | medium | confirmed | C-refdata-9 | Build on an already-built group re-queues finished samples and flattens their presets | DOC-51 |
| REF-M5 | medium | confirmed | C-refdata-10 | PDB hits are stored under the wrong key, so the top reference is never aligned | - |
| REF-M6 | medium | confirmed | C-refdata-12 | Concurrent taxonomy-index builds share one temp file and can install a broken index | - |
| REF-L1 | low | confirmed | C-refdata-7 | When a record count fails, the clade walk skips that rank (`if (is.na(n)) next`) | - |
| REF-L2 | low | confirmed | C-refdata-8 | An existing clade database is reused by its file-name slug alone; the manifest's taxid and RefSeq flag are never compared | - |
| REF-L3 | low | confirmed | C-refdata-11 | The shared reference store keeps records with no annotations, and never compares the `ribopilot_version` it writes | - |
| REF-L4 | low | confirmed; C-pkg-20 unverified-low | C-refdata-15, C-security-12, C-data-29, C-pkg-20 | The default fish seeds and labels come from an unpinned `main` branch URL (currently 404 while the repo is private), with no checksum, and `new_db()` repeats the same URLs as literals | DOC-19, UX-11 |
| REF-L5 | low | confirmed; C-refdata-18 unverified-low | C-nextflow-21, C-refdata-18 | A cluster-submitted build that fails leaves no FASTA, and the poller treats "no fresh FASTA" as still running, with no timeout; the poller also accepts a FASTA that is still being copied | - |
| REF-L6 | low | confirmed | C-security-17 | Build polling trusts a bare pid on the current host | - |
| REF-L7 | low | confirmed; C-refdata-19 unverified-low, C-pkg-27 unverified-low | C-refdata-19, C-security-16, C-pkg-27 | The manifest is written without quotes and read back with quote handling, so `"Mullidae"[Organism]` loses its quotes on the next append; the manifest has no lock, and builds publish straight onto the shared seeds and labels names | - |
| REF-L8 | low | unverified-low | C-refdata-16 | In the container, efetch fails once on every call (curl 56) and succeeds on retry; the HTTP/1.1 shim is on PATH only in the download scripts, and pacing is per task while three tasks run | - |
| REF-L9 | low | unverified-low | C-refdata-17 | `GITHUB_PAT` headers are ignored when a site sets the download method to curl or wget, and the error then blames a slow connection | DOC-19, UX-56 |
| REF-L10 | low | unverified-low | C-refdata-20 | Hybrid names keep only the first two words, so "Citrus maxima x Citrus reticulata" resolves as an exact match to one parent | - |
| REF-L11 | low | unverified-low | C-refdata-21 | `RIBOPILOT_NO_LOCAL_BLAST=0` turns local BLAST off, because any non-empty value counts | DOC-60 |
| REF-L12 | low | unverified-low | C-refdata-22 | "Leave these samples on the built-in fish database" does not revert a group that was already wired to a clade database | DOC-51 |
| REF-L13 | low | unverified-low | C-refdata-23 | `extract_species()` returns prefixes such as "PREDICTED: Danio" or "Chain A5," as the species | - |
| REF-L14 | low | unverified-low | C-refdata-24 | The awk accession scan uses `NR==FNR`, which stays true for stdin when the wanted list is empty | - |
| REF-L15 | low | unverified-low | C-refdata-25 | Each project keeps its own 400 MB taxonomy index plus the 80 MB tarball | DOC-05, DOC-24 |
| REF-L16 | low | unverified-low | C-refdata-26 | The record-count cache never expires and ignores the RefSeq flag | - |
| REF-L17 | low | unverified-low | C-refdata-27 | The file header says in-app builds are detached and survive a closed app, but the app launches them with `cleanup = TRUE` | DOC-36, UX-43 |
| REF-L18 | low | unverified-low | C-refdata-28 | The README ("From R (recommended)" with `conda_env`), the `build_reference_db` help, and the `check_reference_db` hint push conda and scripted builds | DOC-79 |

**Security** ([CODE_REVIEW.md section 4.6](CODE_REVIEW.md#46-security))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| SEC-M1 | medium | confirmed | C-security-3, C-data-6, C-pkg-13 | FASTA header templates are run as R code, including templates stored in the project database | - |
| SEC-L1 | low | confirmed; C-nextflow-26 unverified-low, C-shiny-24 unverified-low | C-security-2, C-refdata-13, C-nextflow-26, C-shiny-24 | The NCBI API key is protected in `.config` (mode 0600) but then passed as `--ncbi_api_key`, so it lands in `ps` output, `.nextflow.log`, `.nextflow/history`, every task's `.command.run` (mode 664), job scripts saved as 0755, the command preview, and the Reference modal's plain text box | DOC-91, UX-47 |
| SEC-L2 | low | confirmed | C-security-7 | The API key is written unescaped into every task's shell as `export NCBI_API_KEY="..."` | DOC-31 |
| SEC-L3 | low | confirmed; C-pkg-15 unverified-low | C-security-11, C-ponytail-36, C-pkg-15 | `.Rbuildignore` misses `.claude`, `.superpowers` (39 files of agent reports and diffs, about 1.9 MB), `.nextflow.log`, `update-version.log`, and `RiboPilot_test_userAsmb.R`; `devtools::build` does not read `.gitignore`, so all of them go into the tarball, and the Dockerfile copies that tarball into its own layer of the public image (the later `rm` does not remove a layer) | DOC-97 |
| SEC-L4 | low | confirmed | C-security-14 | GenBank titles, species, clade names, and `fail_reason` are inserted as raw HTML, and the Taxon is pasted into R Markdown source that is knitted | - |
| SEC-L5 | low | confirmed | C-security-21 | A test-read download that fails its checksum is left in place, and the next run accepts any file that passes a gzip test | - |

**Packaging** ([CODE_REVIEW.md section 4.7](CODE_REVIEW.md#47-packaging))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| PKG-H1 | high | confirmed | C-pkg-1 | A clean install fails: DESCRIPTION has no `biocViews` field | DOC-17 |
| PKG-H2 | high | confirmed | C-data-3, C-pkg-2 | `new_project()` and `new_test_project()` crash outside RStudio with the default `Rproj = TRUE` | DOC-18, UX-01 |
| PKG-M1 | medium | confirmed | C-pkg-6 | Roxygen markdown is off, so every help page shows raw markdown | DOC-72 |
| PKG-L1 | low | confirmed; C-pkg-23 unverified-low | C-data-30, C-pkg-23 | `.onAttach` runs `nextflow -version` on every attach and calls `stop()` when Nextflow is too old | DOC-25 |
| PKG-L2 | low | unverified-low | C-pkg-17 | NAMESPACE is marked "Hand-maintained", so roxygen skips it, although it is byte-identical to roxygen's output | DOC-82 |
| PKG-L3 | low | unverified-low | C-pkg-18 | Tests call `rlang::current_env()`, but rlang is not in Suggests | - |
| PKG-L4 | low | confirmed; C-pkg-19 unverified-low | C-pkg-19, C-ponytail-41 | Several roxygen blocks sit on the wrong object or are orphaned, roxygen 8.1 reports 4 "@noRd must not be followed by any text" errors, and an inline `@noRd` fails so `man/check_container_image.Rd` documents an internal function | DOC-66, DOC-81 |
| PKG-L5 | low | unverified-low | C-pkg-21 | 7 exported pages have no `\value`, only 1 of 53 has examples, and there is no package-level help page (`?RiboPilot` opens the app launcher) | DOC-66, DOC-82 |
| PKG-L6 | low | unverified-low | C-pkg-22 | `run_app()` sets `options(RiboPilot.headless = ...)` and never restores it | - |
| PKG-L7 | low | unverified-low | C-pkg-24, C-pkg-38 | Temp folders are named with `Sys.time()` (and `sample.int()`) | - |
| PKG-L8 | low | unverified-low | C-pkg-25 | `README.md` has curly quotes from pandoc, and a paragraph about filtered reads and `full_size = TRUE` exists only in the `.md`, not in `README.Rmd` | DOC-22, DOC-27 |
| PKG-L9 | low | unverified-low | C-pkg-28 | The Dockerfile's R package list has drifted from Imports (it adds r-config, r-gggenes, r-purrr, and r-svglite, and says Biostrings is the only Bioconductor dependency), and the package is installed with `dependencies = FALSE` | - |
| PKG-L10 | low | unverified-low | C-nextflow-32 | The script commits the version bump before the image build and push, which can fail (for example `docker login` with no TTY) | - |
| PKG-L11 | low | unverified-low | C-tests-19, C-pkg-16 | `ggplot2::aes(y = y)` in `plot_chimera` gives a "no visible binding for global variable y" NOTE, and five `globalVariables` entries (start, label, taxon, pident, src) are no longer needed | - |
| PKG-N1 | nit | unverified-low | C-pkg-26 | A comment says "constants.R collates last", but with no Collate field files load alphabetically and it loads 13th | - |
| PKG-N2 | nit | unverified-low | C-pkg-29 | `LazyData: true` with no `data/`, a Title not in title case, a Description starting "An R package for", and CI building with `--no-build-vignettes` although every vignette chunk is `eval = FALSE` | DOC-82, DOC-99 |

**Tests** ([CODE_REVIEW.md section 6](CODE_REVIEW.md#6-test-suite-results))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| TEST-M1 | medium | confirmed | C-tests-7 | The per-unit pipeline drivers and the WF1 wiring are never exercised | - |
| TEST-L1 | low | confirmed; C-nextflow-31 unverified-low | C-tests-14, C-nextflow-31 | CI runs `R CMD check` on Ubuntu with no Nextflow, sqlite3, or BLAST, so 132 tests skip, including every `blast.nf` end-to-end test, and no job parses `main.nf` with the pinned Nextflow | DOC-96 |
| TEST-L2 | low | confirmed | C-tests-9 | The fake-tool tests prepend their folder to PATH, so a real pybarrnap on PATH still wins | - |
| TEST-L3 | low | confirmed | C-tests-10 | The app's 5 s shutdown timer runs after every `testServer()` session | - |
| TEST-L4 | low | confirmed | C-tests-11 | The test named "ignored and problematic units are excluded from the bundle" only sets `ignore = 1`; no test exports with `problematic = 'yes'` or `annotate_lock = 0` | DOC-41, UX-33 |
| TEST-L5 | low | unverified-low | C-tests-13 | The ITSx roll-up test runs a pasted copy of the app's SQL, not the app | - |
| TEST-L6 | low | unverified-low | C-tests-15 | Many tests grep source text instead of checking behaviour, and one regex stops at the first nested `)` | - |
| TEST-L7 | low | unverified-low | C-tests-16 | The only real-data chimera and ITSx checks need `RIBOPILOT_TEST_PANEL`, `RIBOPILOT_TEST_TAXDUMP`, and `RIBOPILOT_TEST_PROJECT`, which are documented nowhere | - |
| TEST-L8 | low | unverified-low | C-tests-17 | Fixed sleeps, elapsed-time limits under concurrency, and a vacuous `/tmp` check | - |
| TEST-L9 | low | unverified-low | C-tests-18 | Tests expect one warning where the code emits several, so 9 warnings leak on every run, plus raw `tar` error text | - |
| TEST-L10 | low | unverified-low | C-tests-26 | Four `blast.nf` end-to-end tests rebuild the same fixture instead of using `mk_dbarg_project()` | - |
| TEST-N1 | nit | unverified-low | C-tests-20 | A comment says `reference_handlers()` is not mounted by `app_server`, but `R/app_server.R:625` mounts it | - |
| TEST-N2 | nit | unverified-low | C-tests-27 | The app DDL, BLAST table schema, and chimera fakes are copied across files, and positional INSERTs into `reference_config` break on any new column | - |

**Ponytail (delete or simplify)** ([CODE_REVIEW.md section 5](CODE_REVIEW.md#5-ponytail-what-to-delete-or-simplify))

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| P1 | low | confirmed; C-shiny-34 unverified-low | C-shiny-34, C-security-25, C-ponytail-5 | One Nextflow run surface instead of two | - |
| P2 | low | confirmed; C-refdata-30 unverified-low, C-tests-21 unverified-low, C-pkg-32 unverified-low | C-refdata-30, C-tests-21, C-pkg-32, C-ponytail-6 | Delete the dead CSV taxonomy workflow | - |
| P3 | low | confirmed; C-refdata-32 unverified-low | C-refdata-32, C-ponytail-7 | Delete the network `datasets` taxonomy fallback | - |
| P4 | low | confirmed; C-refdata-31 unverified-low, C-pkg-33 unverified-low | C-refdata-31, C-pkg-33, C-ponytail-8 | Delete `install_ribopilot_db()` and helpers | - |
| P5 | low | confirmed; C-nextflow-33 unverified-low, C-refdata-29 unverified-low | C-nextflow-33, C-refdata-29, C-security-24, C-ponytail-9 | Delete the dead `rp_ref_fetch.R` chain | - |
| P6 | low | confirmed; C-shiny-35 unverified-low | C-shiny-35, C-ponytail-10 | Delete dead JS handlers | - |
| P7 | low | confirmed; C-science-36 unverified-low, C-tests-22 unverified-low, C-pkg-31 unverified-low | C-science-36, C-tests-22, C-pkg-31, C-ponytail-11, C-data-33 | Delete `score_project()` | - |
| P8 | low | confirmed; C-shiny-36 unverified-low, C-shiny-37 unverified-low, C-science-39 unverified-low, C-science-40 unverified-low, C-science-42 unverified-low, C-science-43 unverified-low, C-science-44 unverified-low, C-tests-24 unverified-low, C-tests-25 unverified-low, C-pkg-35 unverified-low | C-shiny-36, C-shiny-37, C-science-39, C-science-40, C-science-42, C-science-43, C-science-44, C-tests-24, C-tests-25, C-pkg-35, C-ponytail-29, C-security-29 | Delete small dead functions, arguments, and constants | - |
| P9 | low | confirmed; C-science-37 unverified-low, C-tests-23 unverified-low, C-pkg-34 unverified-low | C-science-37, C-tests-23, C-pkg-34, C-data-32, C-ponytail-24 | One coverage writer and one depth parser | - |
| P10 | low | confirmed; C-pkg-37 unverified-low | C-ponytail-15, C-pkg-37 | Unexport the pipeline internals | DOC-66 |
| P11 | low | confirmed; C-nextflow-40 unverified-low | C-nextflow-40, C-security-30, C-ponytail-16 | One bin-script dispatcher | - |
| P12 | low | confirmed; C-pkg-36 unverified-low | C-data-35, C-ponytail-12, C-ponytail-13, C-ponytail-14, C-pkg-36 | Drop five dependencies | - |
| P13 | low | confirmed | C-ponytail-17 | Drop dead schema | UX-40 |
| P14 | low | confirmed | C-ponytail-18, C-data-39 | Remove back-compat shims | - |
| P15 | low | confirmed | C-ponytail-19 | One busy-alert helper | - |
| P16 | low | confirmed; C-nextflow-35 unverified-low | C-nextflow-35, C-ponytail-20 | One deploy script | - |
| P17 | low | confirmed; C-nextflow-34 unverified-low | C-nextflow-34, C-ponytail-21 | Remove unused image tools and fix the smoke test | - |
| P18 | nit | confirmed | C-ponytail-22 | One pill renderer | - |
| P19 | low | confirmed; C-nextflow-37 unverified-low | C-nextflow-37, C-ponytail-23 | Simpler `blast.nf` | - |
| P20 | nit | confirmed | C-ponytail-25, C-data-34, C-data-37 | Trim `export_bundle()` and the export wrappers | DOC-40 |
| P21 | nit | confirmed | C-ponytail-26 | One triage settings row | - |
| P22 | low | confirmed; C-nextflow-39 unverified-low | C-nextflow-39, C-security-27, C-ponytail-27 | One reference-config reader | - |
| P23 | nit | confirmed | C-ponytail-28 | Build the variant clause at display time | - |
| P24 | nit | confirmed | C-ponytail-33 | Drop unused scoring options | - |
| P25 | nit | confirmed; C-shiny-39 unverified-low | C-ponytail-30, C-shiny-39 | One alert and one toast | - |
| P26 | low | confirmed; C-shiny-38 unverified-low | C-shiny-38, C-data-36, C-ponytail-31 | One exportable-units query | - |
| P27 | low | confirmed; C-nextflow-36 unverified-low | C-nextflow-36, C-ponytail-32 | Merge `config.local` into `config.docker` | - |
| P28 | nit | confirmed | C-ponytail-35 | `RiboPilot <- run_app` | - |
| P29 | nit | confirmed | C-ponytail-37 | Move the `dev/` script tests out of the package suite | - |
| P30 | low | confirmed; C-science-38 unverified-low | C-science-38, C-ponytail-38 | One anchor caller | - |
| P31 | low | unverified-low | C-nextflow-38 | Declare Nextflow params once | - |
| P32 | nit | unverified-low | C-science-41, C-nextflow-41, C-refdata-33 | Small repeats | - |
| P33 | nit | confirmed | C-data-38 | Drop the hand-built HTML report fallback | - |
| P34 | nit | confirmed | C-security-28 | Launch Nextflow with processx in `run_assembly()` | - |

### B.2 Docs review items

Severity scale: blocker, major, minor, polish. All confirmed. Each item links to its write-up in DOCS_REVIEW.md.

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| [DOC-01](DOCS_REVIEW.md#doc-01) | major | confirmed | D-novice-8, D-intermediate-15, D-hpc-11, D-editor-6 | The install list leaves out Docker or Singularity, Java 17+, and a Nextflow version | - |
| [DOC-02](DOCS_REVIEW.md#doc-02) | major | confirmed | D-novice-5 | Supported operating systems are never stated | - |
| [DOC-03](DOCS_REVIEW.md#doc-03) | major | confirmed | D-novice-4, D-hpc-12, D-intermediate-10, D-editor-6 (part) | No hardware, disk, or run-time guidance, and the default resource requests stop small machines | PIPE-H2 |
| [DOC-04](DOCS_REVIEW.md#doc-04) | major | confirmed | D-hpc-6, D-intermediate-9, D-hpc-19, D-hpc-21 | No "Running on a cluster" page: opening the app on a login node, the head-job script, and shared storage | APP-L5, PIPE-L2, PIPE-M14, UX-79 |
| [DOC-05](DOCS_REVIEW.md#doc-05) | major | confirmed | D-hpc-4, D-novice-23, D-intermediate-27 | Nothing says where caches and big files go, how big they get, or how to move them | PIPE-L2, REF-L15 |
| [DOC-06](DOCS_REVIEW.md#doc-06) | major | confirmed | D-hpc-5 | No guidance on which steps need internet, or on compute nodes without it | - |
| [DOC-07](DOCS_REVIEW.md#doc-07) | major | confirmed | D-hpc-7 | Nothing explains what an upgrade does to existing projects | PIPE-M15 |
| [DOC-08](DOCS_REVIEW.md#doc-08) | major | confirmed | D-intermediate-12, D-novice-45, D-intermediate-33, O-edge-27 (part), D-editor-17 (part) | No troubleshooting page; the only troubleshooting text is a conda-only BLAST box | UX-15, UX-42, UX-54, UX-55 |
| [DOC-09](DOCS_REVIEW.md#doc-09) | major | confirmed | D-intermediate-6, D-novice-22, O-novel-8, O-edge-26, O-edge-27 (part), D-editor-17 (part) | The mapping-file rules are not written down | DB-L5, DB-L6, DB-L15, UX-03, UX-06, UX-55 |
| [DOC-10](DOCS_REVIEW.md#doc-10) | major | confirmed | D-intermediate-11 | Where outputs land, and what can be deleted, is not documented | - |
| [DOC-11](DOCS_REVIEW.md#doc-11) | major | confirmed | D-intermediate-13 | There is no way to add samples to an existing project, and the docs do not say so | - |
| [DOC-12](DOCS_REVIEW.md#doc-12) | major | confirmed | D-intermediate-14, D-domain-21, D-editor-17 | No NEWS, no citation, and no tool-version list for methods sections | PIPE-M12 |
| [DOC-13](DOCS_REVIEW.md#doc-13) | major | confirmed | D-novice-18, D-appclaims-19, D-editor-17 (part) | No screenshots, and many app features are never described | UX-51, UX-53 |
| [DOC-14](DOCS_REVIEW.md#doc-14) | minor | confirmed | D-novice-36, D-intermediate-23 | No glossary; sample, path, scaffold, unit, and record are never tied together | - |
| [DOC-15](DOCS_REVIEW.md#doc-15) | major | confirmed | D-editor-10, D-novice-17, D-appclaims-24, O-practice-16 | The docs and the app use different names for the same things | UX-50, UX-77 |
| [DOC-16](DOCS_REVIEW.md#doc-16) | major | confirmed | D-novice-9, D-hpc-10, D-novice-34, D-editor-20, D-hpc-25 | The unsupported conda route is the main install recipe, and the code still defaults to it | PIPE-M11, UX-04 |
| [DOC-17](DOCS_REVIEW.md#doc-17) | blocker | confirmed | D-novice-2, D-intermediate-3, D-hpc-2, D-editor-3, D-api-2 | The install command fails on a fresh R: the Bioconductor dependencies are not found | PKG-H1 |
| [DOC-18](DOCS_REVIEW.md#doc-18) | blocker | confirmed | D-novice-3, D-intermediate-2, D-editor-2, D-api-1 | The documented new_test_project() and new_project() calls crash outside RStudio | PKG-H2, UX-01 |
| [DOC-19](DOCS_REVIEW.md#doc-19) | blocker | confirmed | D-hpc-3, D-novice-19, D-editor-7, D-api-22, D-intermediate-31 (part), O-novel-7 (part) | While the repo is private, three documented steps need a GitHub token, and the default fish reference fails even with one | REF-L4, REF-L9, SCI-M9, UX-11 |
| [DOC-20](DOCS_REVIEW.md#doc-20) | major | confirmed | D-intermediate-4, D-editor-5, D-appclaims-2, D-api-5 (part) | The README says the taxonomy lock gates assembly; the pipeline does not check it | PIPE-L6, UX-09 |
| [DOC-21](DOCS_REVIEW.md#doc-21) | minor | confirmed | D-novice-32, D-editor-25, D-hpc-17, D-api-25 | The Hydra section is out of date: the app runs hydra_setup() itself | - |
| [DOC-22](DOCS_REVIEW.md#doc-22) | minor | confirmed | D-hpc-13, D-novice-33, D-editor-19, D-api-27 (part) | README.md has a paragraph that README.Rmd lacks | PKG-L8 |
| [DOC-23](DOCS_REVIEW.md#doc-23) | minor | confirmed | D-domain-17, D-editor-23, D-novice-39 | The repeat diagram puts LSU under 5.8S, draws ETS as a gene, and has no key | - |
| [DOC-24](DOCS_REVIEW.md#doc-24) | minor | confirmed | D-hpc-15, D-novice-29, D-editor-36 | The taxonomy snapshot's size and time are understated, and its download can time out | REF-H2, REF-L15 |
| [DOC-25](DOCS_REVIEW.md#doc-25) | minor | confirmed | D-novice-31, D-intermediate-32 | The Nextflow version claim and hint point the wrong way | PKG-L1 |
| [DOC-26](DOCS_REVIEW.md#doc-26) | minor | confirmed | D-hpc-18 | environment.yml says the conda and container paths agree; they do not | PIPE-M12 |
| [DOC-27](DOCS_REVIEW.md#doc-27) | polish | confirmed | D-hpc-26, D-editor-42 | The knitted README.md is not ASCII | PKG-L8 |
| [DOC-28](DOCS_REVIEW.md#doc-28) | blocker | confirmed | D-novice-1, D-intermediate-1, D-editor-1, D-appclaims-1, D-novice-10, D-intermediate-5, D-editor-4, D-api-5, O-practice-5, O-practice-6, O-novel-6, O-edge-27, D-hpc-6 (part) | Get started describes an older, console-driven app with three tabs | UX-09, UX-49 |
| [DOC-29](DOCS_REVIEW.md#doc-29) | blocker | confirmed | D-hpc-1 | The documented "Already have an assembly?" route fails on every generate_config() Singularity profile | PIPE-H3 |
| [DOC-30](DOCS_REVIEW.md#doc-30) | major | confirmed | D-novice-12, D-intermediate-20, D-editor-9, D-api-9, D-novice-41 | RiboPilot("~/my_project") fails as written, and the summary block uses an undefined variable | - |
| [DOC-31](DOCS_REVIEW.md#doc-31) | major | confirmed | D-novice-6, D-editor-22, O-practice-21 | The NCBI API key: the placeholder is sent to NCBI if left in, and the environment variable is ignored | SEC-L2, UX-07 |
| [DOC-32](DOCS_REVIEW.md#doc-32) | major | confirmed | D-novice-7, D-intermediate-16, O-practice-20, D-editor-21 (part) | Step 2 sends non-fish users to an unsupported console build before the project exists | UX-50 |
| [DOC-33](DOCS_REVIEW.md#doc-33) | polish | confirmed | D-domain-24 | "If your samples are fish" should say ray-finned fishes | - |
| [DOC-34](DOCS_REVIEW.md#doc-34) | minor | confirmed | D-novice-37, D-editor-21 | Code fragments need a comma the reader is not told about, and rely on an in-memory object | - |
| [DOC-35](DOCS_REVIEW.md#doc-35) | major | confirmed | D-novice-13, D-appclaims-18, D-domain-22, O-practice-19, D-api-26, D-intermediate-19 (part), D-editor-21 (part) | The practice project leaves the reader without a next step, and its samples are unexplained | UX-11 |
| [DOC-36](DOCS_REVIEW.md#doc-36) | major | confirmed | D-appclaims-6, D-novice-11 | Closing a run window or the browser tab stops the run, but the docs say you can keep working or come back later | APP-M1, REF-L17, UX-43 |
| [DOC-37](DOCS_REVIEW.md#doc-37) | major | confirmed | D-appclaims-8, O-practice-10, D-novice-10 (part), D-intermediate-16 (part) | Step 5 rescue: the app route is undocumented, and run_fallback() re-assembles good samples | PIPE-H5, UX-15, UX-16, UX-17 |
| [DOC-38](DOCS_REVIEW.md#doc-38) | major | confirmed | D-appclaims-5 | The Assemble-step instructions describe controls and re-runs that do not exist | - |
| [DOC-39](DOCS_REVIEW.md#doc-39) | minor | confirmed | D-domain-4 | "The spacers are simply the stretches between those landmarks" is not what ships: 1 bp and 4 bp are left unannotated | DB-M5, UX-84 |
| [DOC-40](DOCS_REVIEW.md#doc-40) | major | confirmed | D-novice-21, D-intermediate-8, D-api-6, D-appclaims-3, O-practice-14, O-practice-26 | export_project() does not "export everything" | P20, UX-37 |
| [DOC-41](DOCS_REVIEW.md#doc-41) | major | confirmed | D-intermediate-22, D-appclaims-4, D-editor-16, D-novice-28, O-practice-17 | Step 9 leaves out behaviors that change what is submitted; Problematic silently drops units | TEST-L4, UX-25, UX-33, UX-34 |
| [DOC-42](DOCS_REVIEW.md#doc-42) | major | confirmed | D-domain-1 | Partial or complete is decided by scaffold position, not gene completeness, and no doc explains it | DB-H2, DB-L8, UX-32 |
| [DOC-43](DOCS_REVIEW.md#doc-43) | major | confirmed | D-domain-12 | Exports carry unannotated ETS and IGS flanks; trimming and Portal length limits are never mentioned | - |
| [DOC-44](DOCS_REVIEW.md#doc-44) | major | confirmed | D-domain-13, D-appclaims-13, D-intermediate-7 | Required source modifiers are not stated, and isolate is banned in templates but added automatically | UX-36, UX-40 |
| [DOC-45](DOCS_REVIEW.md#doc-45) | major | confirmed | D-domain-5, D-intermediate-7 (part) | The default [organism=] uses the raw Taxon label, not the name confirmed in the Taxonomy step | DB-M1 |
| [DOC-46](DOCS_REVIEW.md#doc-46) | minor | confirmed | D-novice-38, D-domain-23 | Step 9 does not say which files to upload, and README.txt names a table2asn / BankIt route without steps | - |
| [DOC-47](DOCS_REVIEW.md#doc-47) | minor | confirmed | D-appclaims-12, D-intermediate-29, D-novice-40, D-api-15, O-practice-25, O-practice-27, O-practice-33 | The export layout block and the export_bundle help have small errors | DB-L2, DB-L3, UX-38, UX-39 |
| [DOC-48](DOCS_REVIEW.md#doc-48) | minor | confirmed | D-domain-14, D-domain-15 | The ITSx notes in README.txt and Step 9 overstate: "ITS too long" is called misassembly, and NCBI is said to remove records "silently" | - |
| [DOC-49](DOCS_REVIEW.md#doc-49) | minor | confirmed | D-domain-16 | Every N-run in the .tbl claims linkage_evidence "align genus" | DB-L11 |
| [DOC-50](DOCS_REVIEW.md#doc-50) | minor | confirmed | D-novice-24, D-domain-20, D-editor-24 | The 25,000-record ceiling is described as a hard limit | UX-12 |
| [DOC-51](DOCS_REVIEW.md#doc-51) | minor | confirmed | D-novice-25, D-appclaims-17, O-novel-10, D-editor-24 (part) | Reference tab details differ from the app | REF-L12, REF-M4, UX-50 |
| [DOC-52](DOCS_REVIEW.md#doc-52) | minor | confirmed | D-novice-20 | The Reference tab is called optional, which non-fish users will take literally | - |
| [DOC-53](DOCS_REVIEW.md#doc-53) | major | confirmed | D-novice-14, D-domain-9, D-editor-11, D-appclaims-7, D-intermediate-19, O-novel-9 (part) | The flag table misses flags users see, and misstates what happens without a lookup | UX-51 |
| [DOC-54](DOCS_REVIEW.md#doc-54) | major | confirmed | D-novice-16, D-domain-7, D-editor-12, D-appclaims-10, O-novel-9, D-intermediate-19 (part) | Five rulesets ship and are picked from the lineage; the docs describe two | UX-51 |
| [DOC-55](DOCS_REVIEW.md#doc-55) | major | confirmed | D-novice-15, D-editor-13, D-domain-19, D-appclaims-11 | The triage section counts seven rules; eight ship | - |
| [DOC-56](DOCS_REVIEW.md#doc-56) | major | confirmed | D-domain-10, D-editor-13 (part), D-novice-15 (part), D-intermediate-19 (part) | Chimera detection is not documented, and the in-app description is incomplete | - |
| [DOC-57](DOCS_REVIEW.md#doc-57) | major | confirmed | D-domain-11 | The chimera check scores a clean flatfish the same as the planted chimera, so a good record is exported as chimeric | SCI-M3, SCI-M4, UX-82, UX-83 |
| [DOC-58](DOCS_REVIEW.md#doc-58) | major | confirmed | D-domain-2 | The isomer test ignores the ITS, so real ITS variants are told the pick "does not matter" | SCI-L5, SCI-M7 |
| [DOC-59](DOCS_REVIEW.md#doc-59) | major | confirmed | D-domain-3 | "Cannot be contamination" contradicts the phylum-level screen described on the same page | UX-86 |
| [DOC-60](DOCS_REVIEW.md#doc-60) | major | confirmed | D-domain-6, D-appclaims-14, D-intermediate-24, D-editor-27 | The "GenBank screen" is by default a local panel, but the docs describe an NCBI search | REF-L11 |
| [DOC-61](DOCS_REVIEW.md#doc-61) | minor | confirmed | D-domain-18 | The depth explanation does not match how depth is computed | SCI-M2, UX-83 |
| [DOC-62](DOCS_REVIEW.md#doc-62) | minor | confirmed | D-editor-26 | The rule 7 note is about the wrong thing and contradicts the tooltip | - |
| [DOC-63](DOCS_REVIEW.md#doc-63) | minor | confirmed | D-novice-27 | The page explains what flags mean but not what to do, and the triage prose is written for experts | - |
| [DOC-64](DOCS_REVIEW.md#doc-64) | polish | confirmed | D-novice-42, D-editor-34 | Spaced hyphens are used as dashes | - |
| [DOC-65](DOCS_REVIEW.md#doc-65) | polish | confirmed | D-editor-35 | The curation tables use equal column widths | - |
| [DOC-66](DOCS_REVIEW.md#doc-66) | major | confirmed | D-editor-14, D-api-8, D-api-28, D-api-11, D-novice-35, D-intermediate-30 | The reference index is one flat list of 54 topics, most of them pipeline internals | P10, PKG-L4, PKG-L5 |
| [DOC-67](DOCS_REVIEW.md#doc-67) | major | confirmed | D-hpc-9, D-intermediate-18, D-api-3, D-editor-15 (part) | new_project() help says NMNH_Hydra and NOAA_SEDNA need generate_config() and calls conda ready to use | - |
| [DOC-68](DOCS_REVIEW.md#doc-68) | major | confirmed | D-api-4, D-intermediate-21, D-appclaims-9, D-editor-15, D-novice-30 | run_assembly() help says annotation waits for a second run after locking | - |
| [DOC-69](DOCS_REVIEW.md#doc-69) | major | confirmed | D-api-7, D-novice-26, D-intermediate-26 | The headless (job-script) rule is described wrongly, and launch.browser = FALSE silently removes Run from App | UX-79 |
| [DOC-70](DOCS_REVIEW.md#doc-70) | major | confirmed | D-intermediate-17, D-api-14 | A failed new_project() leaves a half-made project, and the retry silently skips the config | DB-L4, UX-02 |
| [DOC-71](DOCS_REVIEW.md#doc-71) | minor | confirmed | D-hpc-8 | generate_config() ignores container_engine for "local" and queue for "awsbatch" | PIPE-M10 |
| [DOC-72](DOCS_REVIEW.md#doc-72) | minor | confirmed | D-editor-18, D-api-10, D-intermediate-30 (part), D-hpc-16 (part) | Roxygen markdown is off: raw backticks and [fn()] brackets on every page, and no cross-links | PKG-M1 |
| [DOC-73](DOCS_REVIEW.md#doc-73) | minor | confirmed | D-hpc-16, D-api-17 | generate_config() help has no example and does not say where profiles are saved or how to share one | - |
| [DOC-74](DOCS_REVIEW.md#doc-74) | minor | confirmed | D-api-13, D-novice-30 (part) | new_project() help omits ID rules, the return value, and the data loss from force, and typos pass silently | - |
| [DOC-75](DOCS_REVIEW.md#doc-75) | minor | confirmed | D-api-12 | The RiboPilot() and run_app() pages lack a return value and explain defaults poorly | - |
| [DOC-76](DOCS_REVIEW.md#doc-76) | minor | confirmed | D-api-16, D-intermediate-28 | build_reference_db(): count_only reports success when tools are missing, and a retry printed the API key | UX-47 |
| [DOC-77](DOCS_REVIEW.md#doc-77) | minor | confirmed | D-api-19 | Pipeline jargon in user-facing titles and messages | - |
| [DOC-78](DOCS_REVIEW.md#doc-78) | minor | confirmed | D-api-20 | new_db() is exported with a stale description and misleading defaults | - |
| [DOC-79](DOCS_REVIEW.md#doc-79) | minor | confirmed | D-api-21, D-editor-28, D-intermediate-31, O-novel-7, D-novice-13 (part) | Console messages point to .config, conda, or maintainer scripts instead of the app | REF-L18, UX-56 |
| [DOC-80](DOCS_REVIEW.md#doc-80) | minor | confirmed | D-api-18 | The ribopilot_ref_dir() page lacks a return value, an example, and the real default path | - |
| [DOC-81](DOCS_REVIEW.md#doc-81) | polish | confirmed | D-editor-40, D-api-27 (part) | roxygen2 warnings from glued comment blocks | PKG-L4 |
| [DOC-82](DOCS_REVIEW.md#doc-82) | polish | confirmed | D-api-27 | Small help-page and metadata issues | PKG-L2, PKG-L5, PKG-N2 |
| [DOC-83](DOCS_REVIEW.md#doc-83) | minor | confirmed | D-hpc-14, D-api-24, D-editor-38, D-intermediate-35, D-novice-44, D-appclaims-25 | docker/README.md examples use tags that do not exist, trigger a false warning, and give stale mount advice | - |
| [DOC-84](DOCS_REVIEW.md#doc-84) | polish | confirmed | D-novice-43, D-intermediate-34, D-api-24 (part) | The generated .config has the mount list pasted into a comment | UX-81 |
| [DOC-85](DOCS_REVIEW.md#doc-85) | minor | confirmed | D-hpc-20 | Every Singularity task log starts with duplicate-bind and RIBOPILOT_SRC noise | PIPE-M15 |
| [DOC-86](DOCS_REVIEW.md#doc-86) | polish | confirmed | D-hpc-22, D-hpc-23 | The Docker Hub page is blank, and config.local has a garbled sentence | - |
| [DOC-87](DOCS_REVIEW.md#doc-87) | major | confirmed | D-domain-8, D-intermediate-25 | The Curation options "Min anchor identity" and "Flag contaminants" controls do nothing | APP-M13 |
| [DOC-88](DOCS_REVIEW.md#doc-88) | minor | confirmed | D-appclaims-15, D-appclaims-21 | Several tooltips and option-window texts describe behavior the app no longer has | APP-M2, UX-42 |
| [DOC-89](DOCS_REVIEW.md#doc-89) | minor | confirmed | D-editor-29, D-appclaims-25 (part) | Help opens the GitHub repo, the Samples view says "groups", and the Triage window uses jargon | APP-L8, UX-52, UX-78 |
| [DOC-90](DOCS_REVIEW.md#doc-90) | minor | confirmed | D-appclaims-16 | Locking a unit in Annotation Review does not freeze it, and the one-path rule depends on which button you press | APP-M7, APP-M15, UX-20, UX-25 |
| [DOC-91](DOCS_REVIEW.md#doc-91) | minor | confirmed | D-appclaims-22 | The run windows show the NCBI API key in plain text | SEC-L1, UX-47 |
| [DOC-92](DOCS_REVIEW.md#doc-92) | polish | confirmed | D-appclaims-23, D-editor-43 | Jargon and developer-facing text in app messages | UX-76, UX-80 |
| [DOC-93](DOCS_REVIEW.md#doc-93) | polish | confirmed | D-domain-25 | Small naming and wording mismatches in app and guide text | - |
| [DOC-94](DOCS_REVIEW.md#doc-94) | minor | confirmed | D-editor-8, D-editor-41 | Code in the alert boxes renders with curly quotes and breaks copy and paste, and the alert CSS is copied into each vignette | - |
| [DOC-95](DOCS_REVIEW.md#doc-95) | minor | confirmed | D-editor-32, D-api-23, D-intermediate-30 (part) | check_pkgdown fails on the DESCRIPTION URL, the favicons are stale, and the theme upper-cases every heading | - |
| [DOC-96](DOCS_REVIEW.md#doc-96) | minor | confirmed | D-editor-30 | The deploy workflow only builds when a commit message contains "BUILD SITE" | TEST-L1 |
| [DOC-97](DOCS_REVIEW.md#doc-97) | minor | confirmed | D-editor-31 | A stale Hydra driver script and internal notes sit at the repo root | SEC-L3 |
| [DOC-98](DOCS_REVIEW.md#doc-98) | minor | confirmed | D-editor-39, D-editor-33 | Article group names, titles, and link text do not match, and Get started is too long | - |
| [DOC-99](DOCS_REVIEW.md#doc-99) | polish | confirmed | D-editor-37, D-hpc-24, D-api-27 (part) | The DESCRIPTION Title and Description undersell the package and never mention containers | PKG-N2 |

### B.3 Usability review items

Severity scale: blocker, major, minor, polish. All confirmed. Each item links to its write-up in USABILITY_REVIEW.md.

| Item | Severity | Status | Ids | Problem | Same issue in other tracks |
|---|---|---|---|---|---|
| [UX-01](USABILITY_REVIEW.md#ux-01) | blocker | confirmed | O-practice-1, O-edge-1, O-novel-2 | The guide's first command fails outside RStudio | PKG-H2, DOC-18 |
| [UX-02](USABILITY_REVIEW.md#ux-02) | major | confirmed | O-edge-5, O-edge-6 | A failed setup leaves a half-made project, and running it again quietly does nothing | DB-L4, DOC-70 |
| [UX-03](USABILITY_REVIEW.md#ux-03) | major | confirmed | O-edge-12, O-edge-13, O-edge-3, O-edge-28 | The mapping file is accepted or changed silently | DB-L5, DOC-09 |
| [UX-04](USABILITY_REVIEW.md#ux-04) | major | confirmed | O-edge-8 | The default executor is conda, which the docs call unsupported | PIPE-M11, DOC-16 |
| [UX-05](USABILITY_REVIEW.md#ux-05) | major | confirmed | O-edge-15, O-edge-25 | A gzipped assembly mixed with plain ones aborts the run; a FASTQ given as an assembly fails silently | PIPE-L9 |
| [UX-06](USABILITY_REVIEW.md#ux-06) | minor | confirmed | O-novel-8, O-edge-26 | The mapping-file rules are not written down, and assembly-only samples need R1/R2 columns anyway | DB-L15, DOC-09 |
| [UX-07](USABILITY_REVIEW.md#ux-07) | minor | confirmed | O-practice-21 | `new_project()` ignores `NCBI_API_KEY` in the environment | DOC-31 |
| [UX-08](USABILITY_REVIEW.md#ux-08) | major | confirmed | O-edge-16, O-edge-17 | Assign Taxonomy counts unresolved samples as assigned, and sends an ambiguous genus to a plant | APP-L16, REF-H3 |
| [UX-09](USABILITY_REVIEW.md#ux-09) | major | confirmed | O-edge-11, O-practice-7, O-practice-6 | Update runs samples whose taxonomy is not locked | DOC-20, DOC-28 |
| [UX-10](USABILITY_REVIEW.md#ux-10) | minor | confirmed | U-design-20 | A step that has not run yet is shown in error red | - |
| [UX-11](USABILITY_REVIEW.md#ux-11) | major | confirmed | O-practice-2, O-edge-14, O-novel-7, O-practice-22, O-practice-23 | The built-in fish reference is unreachable, so the practice project's first run aborts | REF-L4, SCI-M9, DOC-19, DOC-35 |
| [UX-12](USABILITY_REVIEW.md#ux-12) | minor | confirmed | O-practice-8, O-practice-29 | Fish users are told to skip the Reference step, but Propose groups offers an over-limit clade | DOC-50 |
| [UX-13](USABILITY_REVIEW.md#ux-13) | minor | confirmed | O-edge-18 | The existing-database picker pre-selects an empty database | - |
| [UX-14](USABILITY_REVIEW.md#ux-14) | major | confirmed | U-novice-1, U-mitopilot-7, U-curator-8, O-practice-12, U-novice-13, U-curator-19 | The Assemble table hides screen and chimera verdicts, so contaminated samples look clean | - |
| [UX-15](USABILITY_REVIEW.md#ux-15) | major | confirmed | O-novel-4, U-novice-5, O-practice-11, U-power-8, U-curator-14 | Failed samples give no plain reason and no next step | PIPE-M6, PIPE-M9, DOC-08, DOC-37 |
| [UX-16](USABILITY_REVIEW.md#ux-16) | major | confirmed | O-practice-10, U-power-14, U-novice-12, U-novice-6 | The rescue assembler is hidden and not documented for the app | PIPE-H5, DOC-37 |
| [UX-17](USABILITY_REVIEW.md#ux-17) | major | confirmed | O-practice-9, U-novice-6, U-power-1 | A rescue run skips its own new scaffolds and leaves the sample "In progress" | DB-H1, SCI-M5, DOC-37 |
| [UX-18](USABILITY_REVIEW.md#ux-18) | major | confirmed | O-novel-3 | Save can point samples at a parameter set that was never written | APP-M8 |
| [UX-19](USABILITY_REVIEW.md#ux-19) | minor | confirmed | U-power-7 | The Length column sorts as text | - |
| [UX-20](USABILITY_REVIEW.md#ux-20) | major | confirmed | U-mitopilot-1, U-mitopilot-3, U-curator-10 | Locks do not protect what they lock | APP-L10, APP-M7, DOC-90 |
| [UX-21](USABILITY_REVIEW.md#ux-21) | major | confirmed | U-mitopilot-2 | The annotation options guard sends the user round in a loop | - |
| [UX-22](USABILITY_REVIEW.md#ux-22) | major | confirmed | U-design-2, U-curator-9, U-curator-12 | Annotation Review hides the contamination and chimera flags | - |
| [UX-23](USABILITY_REVIEW.md#ux-23) | minor | confirmed | U-design-1, U-power-11, U-novice-14 | Toolbar toggles and Lock flip the whole selection and never show which way | - |
| [UX-24](USABILITY_REVIEW.md#ux-24) | minor | confirmed | U-mitopilot-6 | Annotation Review has no Update button, and the run window does not say what will run | PIPE-L6 |
| [UX-25](USABILITY_REVIEW.md#ux-25) | minor | confirmed | U-mitopilot-8 | Lock & Close marks Reviewed, but toolbar Lock does not | APP-M15, DOC-41, DOC-90 |
| [UX-26](USABILITY_REVIEW.md#ux-26) | minor | confirmed | U-design-21, U-mitopilot-16 | The notes box has no label, and deleting its text does not delete the note | - |
| [UX-27](USABILITY_REVIEW.md#ux-27) | minor | confirmed | U-curator-13, U-mitopilot-17 | The boundary editor lacks context, and there is no sequence viewer | - |
| [UX-28](USABILITY_REVIEW.md#ux-28) | minor | confirmed | O-practice-18, U-curator-17 | The reference alignment panel wrongly advises running Update again | - |
| [UX-29](USABILITY_REVIEW.md#ux-29) | polish | confirmed | U-power-16 | Lock & Close does not move on to the next unit | - |
| [UX-30](USABILITY_REVIEW.md#ux-30) | polish | confirmed | U-curator-20 | Small cosmetic issues in the review windows | - |
| [UX-31](USABILITY_REVIEW.md#ux-31) | blocker | confirmed | U-curator-1, O-novel-1, U-curator-7 | Contaminant records export under the sample's species name with no gate, and can get a green all-clear | - |
| [UX-32](USABILITY_REVIEW.md#ux-32) | blocker | confirmed | U-curator-2, O-practice-3 | Truncated 18S and 28S pieces are exported as "complete sequence" | DB-H2, DB-L8, DOC-42 |
| [UX-33](USABILITY_REVIEW.md#ux-33) | major | confirmed | U-novice-3, U-mitopilot-4, U-design-14, U-curator-10, O-practice-17 | "Problematic" silently drops a unit from export | TEST-L4, DOC-41 |
| [UX-34](USABILITY_REVIEW.md#ux-34) | major | confirmed | U-novice-4, U-design-10, U-mitopilot-21, U-curator-18 | The Export complete dialog does not name the flagged records, so a suspected chimera ships | DOC-41 |
| [UX-35](USABILITY_REVIEW.md#ux-35) | major | confirmed | U-novice-2, U-mitopilot-5, O-edge-7 | Spaces or slashes in names and paths break runs and export | DB-M2, DB-M3, PIPE-M13 |
| [UX-36](USABILITY_REVIEW.md#ux-36) | major | confirmed | U-curator-4, O-practice-15 | The header template blocks any isolate, and pieces of one specimen get made-up isolates | DOC-44 |
| [UX-37](USABILITY_REVIEW.md#ux-37) | major | confirmed | O-practice-14 | `export_project()` errors unless groups were set in the app | DOC-40 |
| [UX-38](USABILITY_REVIEW.md#ux-38) | minor | confirmed | O-practice-26, O-practice-27, O-practice-33 | The export layout in the guide and the bundle README is out of date | DB-L2, DB-L3, DOC-47 |
| [UX-39](USABILITY_REVIEW.md#ux-39) | minor | confirmed | U-curator-18, O-practice-25, O-practice-32, U-novice-20 | The summary CSV and HTML report are hard to read and partly wrong | DOC-47 |
| [UX-40](USABILITY_REVIEW.md#ux-40) | minor | confirmed | U-novice-15, U-curator-4 | The header template's column list shows internal fields and has no preview | P13, DOC-44 |
| [UX-41](USABILITY_REVIEW.md#ux-41) | polish | confirmed | U-design-28 | The Assign Group dialog is awkward from the keyboard | - |
| [UX-42](USABILITY_REVIEW.md#ux-42) | major | confirmed | O-edge-2, U-power-6, U-design-9, O-practice-4, O-novel-5 | Every run ends with a green "Done", even when it aborted or failed | APP-M2, DOC-08, DOC-88 |
| [UX-43](USABILITY_REVIEW.md#ux-43) | major | confirmed | O-edge-10, U-power-4 | Reloading or closing the tab silently kills a running pipeline | APP-M1, REF-L17, DOC-36 |
| [UX-44](USABILITY_REVIEW.md#ux-44) | major | confirmed | U-power-5 | Stop leaves task containers running and gives no feedback | APP-M1 |
| [UX-45](USABILITY_REVIEW.md#ux-45) | minor | confirmed | U-power-10, U-mitopilot-10 | The run window does not list what will run, and runs with nothing queued | PIPE-L6 |
| [UX-46](USABILITY_REVIEW.md#ux-46) | minor | confirmed | O-edge-9 | Two browser sessions can run at once and overwrite each other's locks | APP-L11, APP-M9 |
| [UX-47](USABILITY_REVIEW.md#ux-47) | minor | confirmed | O-novel-11, U-power-18 | The NCBI API key is shown in plain text and saved with the project | SEC-L1, DOC-76, DOC-91 |
| [UX-48](USABILITY_REVIEW.md#ux-48) | polish | confirmed | O-practice-31, U-novice-21, U-mitopilot-19, U-power-21 | The progress box starts with a Nextflow update nag and a stray escape fragment | - |
| [UX-49](USABILITY_REVIEW.md#ux-49) | major | confirmed | O-practice-5, O-practice-6, O-novel-6, O-edge-27, U-novice-7, U-mitopilot-13, U-power-17 | The Get started guide describes a different app | DOC-28 |
| [UX-50](USABILITY_REVIEW.md#ux-50) | minor | confirmed | O-practice-16, O-practice-19, O-practice-20, O-novel-10 | Smaller mismatches between the guides and the app | DOC-15, DOC-32, DOC-51 |
| [UX-51](USABILITY_REVIEW.md#ux-51) | minor | confirmed | O-novel-9, U-curator-15, O-practice-17, U-power-20 | Several flags, rulesets, review marks, and shortcuts are not documented | DOC-13, DOC-53, DOC-54 |
| [UX-52](USABILITY_REVIEW.md#ux-52) | minor | confirmed | O-practice-24, U-novice-8, U-mitopilot-12, U-design-24 | The Help link opens the private GitHub repo, not the docs | DOC-89 |
| [UX-53](USABILITY_REVIEW.md#ux-53) | minor | confirmed | U-novice-9, U-mitopilot-11, U-design-23 | The app always opens on step 1, with no "what next" hint | DOC-13 |
| [UX-54](USABILITY_REVIEW.md#ux-54) | minor | confirmed | U-mitopilot-15, U-power-9 | The "output" button opens Nextflow work folders, with no log viewer | APP-M10, DOC-08 |
| [UX-55](USABILITY_REVIEW.md#ux-55) | major | confirmed | O-edge-4, O-edge-19, O-edge-20, O-edge-21, O-edge-22 | Common mapping-file mistakes give R internals instead of plain messages | DOC-08, DOC-09 |
| [UX-56](USABILITY_REVIEW.md#ux-56) | minor | confirmed | O-edge-24, O-edge-23 | Some hints are aimed at developers or blame the wrong cause | REF-L9, DOC-79 |
| [UX-57](USABILITY_REVIEW.md#ux-57) | major | confirmed | U-power-3 | Lock and State filters hide rows with CSS, so pages come up blank at scale | - |
| [UX-58](USABILITY_REVIEW.md#ux-58) | major | confirmed | U-power-2, U-power-13 | Bulk actions include rows hidden by search or column filters | APP-M12 |
| [UX-59](USABILITY_REVIEW.md#ux-59) | minor | confirmed | U-power-15 | The browser stays busy while the page is idle | APP-M6 |
| [UX-60](USABILITY_REVIEW.md#ux-60) | minor | confirmed | U-power-12 | CSV downloads are raw database dumps | - |
| [UX-61](USABILITY_REVIEW.md#ux-61) | minor | confirmed | U-novice-10, U-mitopilot-18, U-design-26 | Result columns sit off screen, behind five option columns, and headers are clipped | - |
| [UX-62](USABILITY_REVIEW.md#ux-62) | minor | confirmed | U-design-18, U-design-19 | Layout problems at 1024x768 and 320 px, and Save below the fold | - |
| [UX-63](USABILITY_REVIEW.md#ux-63) | minor | confirmed | U-design-13, U-novice-16 | Empty tables show a misleading filter message, and the real reason is clipped | - |
| [UX-64](USABILITY_REVIEW.md#ux-64) | polish | confirmed | U-design-27, U-power-23 | Button order and control styles are inconsistent | - |
| [UX-65](USABILITY_REVIEW.md#ux-65) | minor | confirmed | U-design-3 | Focus drops to the top of the page when any dialog closes | - |
| [UX-66](USABILITY_REVIEW.md#ux-66) | minor | confirmed | U-design-4 | The step switcher shows no keyboard focus | - |
| [UX-67](USABILITY_REVIEW.md#ux-67) | minor | confirmed | U-design-5 | Some radios and checkboxes have no accessible name | - |
| [UX-68](USABILITY_REVIEW.md#ux-68) | minor | confirmed | U-design-6 | Toasts are not announced and can cover the main button | - |
| [UX-69](USABILITY_REVIEW.md#ux-69) | minor | confirmed | U-design-7 | Text on the dark sequence and evidence boxes has low contrast | - |
| [UX-70](USABILITY_REVIEW.md#ux-70) | minor | confirmed | U-design-8 | Keep/Drop controls are tiny letters, with state shown by fill only | - |
| [UX-71](USABILITY_REVIEW.md#ux-71) | minor | confirmed | U-design-11 | Tables need many Tab presses per row | - |
| [UX-72](USABILITY_REVIEW.md#ux-72) | minor | confirmed | U-design-12 | Tooltip-only text is not reachable by keyboard | - |
| [UX-73](USABILITY_REVIEW.md#ux-73) | minor | confirmed | U-design-15, U-design-16 | Dialogs lack dialog semantics, and Escape and the close X vary | - |
| [UX-74](USABILITY_REVIEW.md#ux-74) | minor | confirmed | U-design-17 | Plots have the alt text "Plot object" | - |
| [UX-75](USABILITY_REVIEW.md#ux-75) | minor | confirmed | U-design-22 | Hard-coded grey text and badges fall below 4.5:1 | - |
| [UX-76](USABILITY_REVIEW.md#ux-76) | minor | confirmed | U-novice-11, U-design-29, U-curator-15 | Internal codes and jargon are shown to users | DOC-92 |
| [UX-77](USABILITY_REVIEW.md#ux-77) | minor | confirmed | U-mitopilot-14, U-novice-19 | One scaffold has several different names | DOC-15 |
| [UX-78](USABILITY_REVIEW.md#ux-78) | polish | confirmed | O-novel-14, U-novice-17, U-power-19, U-design-25 | The Reference tab's Samples view says "groups" | APP-L8, DOC-89 |
| [UX-79](USABILITY_REVIEW.md#ux-79) | polish | confirmed | O-novel-15, O-edge-29, U-mitopilot-20, U-power-22, U-design-30 | The console says "RiboPilot is headless" when it is not, and other launcher wording slips | DOC-04, DOC-69 |
| [UX-80](USABILITY_REVIEW.md#ux-80) | polish | confirmed | U-novice-18, U-power-21 | Toasts and confirms use programmer wording | DOC-92 |
| [UX-81](USABILITY_REVIEW.md#ux-81) | polish | confirmed | O-practice-30, O-practice-34 | Small text glitches in a dialog and the config file | DOC-84 |
| [UX-82](USABILITY_REVIEW.md#ux-82) | major | confirmed | O-practice-13, U-curator-6 | The chimera check breaks bitscore ties arbitrarily, giving false calls | SCI-M4, DOC-57 |
| [UX-83](USABILITY_REVIEW.md#ux-83) | minor | confirmed | O-practice-28, U-curator-11 | Read and depth signals alone raise chimera flags, including at scaffold ends | SCI-M2, DOC-57, DOC-61 |
| [UX-84](USABILITY_REVIEW.md#ux-84) | major | confirmed | U-curator-3, O-novel-12 | Adjacent features do not touch, and the gaps are never flagged | DB-M5, DOC-39 |
| [UX-85](USABILITY_REVIEW.md#ux-85) | major | confirmed | U-curator-5 | The variant check runs before annotation finishes, so its result depends on timing and is usually blank | SCI-M7 |
| [UX-86](USABILITY_REVIEW.md#ux-86) | minor | confirmed | U-curator-16 | "on_target" is green even when hits agree only at class or phylum | DOC-59 |
| [UX-87](USABILITY_REVIEW.md#ux-87) | polish | confirmed | O-novel-13 | ITSx origin calls on tiny fragments are noise | - |
