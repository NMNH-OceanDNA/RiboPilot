# RiboPilot code review

Date: 2026-09-23. Code reviewed: `main` at 3362308 (version 0.0.0.9018), read on branch `review/2026-09-full-review`.

Nine reviewers looked at the code: the Shiny app (code-shiny), the Nextflow pipeline and containers (code-nextflow), the science and algorithms (code-science), the database and export (code-data), reference data and taxonomy (code-refdata), security (code-security), the tests (code-tests), over-engineering (code-ponytail), and packaging (code-pkg). Every medium or higher finding was then checked by a separate verifier, who could confirm it, refute it, mark it uncertain, or change its severity. Low and nit findings were not all sent to a verifier, to save usage; those are marked **unverified-low**.

How to read the status labels:

- **confirmed**: a verifier re-checked the code (and often re-ran a reproduction) and agrees it is real.
- **uncertain**: the code behaves as described, but the verifier could not show it causes harm in practice.
- **unverified-low**: a low or nit finding that no verifier checked. Treat it as a strong lead, not a proven bug.
- **refuted**: the verifier showed it is not a bug. These are listed only in the appendix.

No finding ended as "disputed". Where two reviewers found the same problem, the items are merged below and every id is kept. Where verifiers disagreed on severity, the item uses the higher one and says so.

Reviewer notes, reproduction scripts, and logs are in `/media/dmacguig/Mockra/DJMCompBio/temp/RiboPilot_review/` (reports in `reports/code-*.md`, scripts in each reviewer's folder). That folder is temporary.

---

## 1. Summary

RiboPilot is in good shape for its stage. The test suite passes (3,952 checks, 0 failures on a clean PATH), `R CMD check` gives no errors or warnings at the default level (2 NOTEs), roxygen output matches the committed NAMESPACE and man pages, and the reviewers found no SQL injection (every query value is a bound parameter) and no schema drift. Most of what follows is about what happens on the second run, on a new user's machine, or on a cluster, not the first run on the author's machine.

The main themes, most important first:

1. **Running an assembly again mixes old and new results.** This is the biggest cluster. When a sample is assembled a second time, the pipeline can load old path files left in the output folder, keep a locked unit's old gene coordinates on the new sequence, keep the old contamination-screen verdict, skip BLAST because the sequence did not change, and forget the curator's keep/drop decisions. Each one is a separate finding, but they share one root cause: nothing clears the evidence tied to a unit when its sequence changes. One focused change (read the task's own input files instead of the publish folder, and clear or recompute everything derived from a unit when its sequence changes) would fix most of them.
2. **Export can send wrong details to GenBank.** A 28S gene that barrnap says is only 46 percent present is written as "complete sequence" with no partial markers. ITS spacers stop 1 to 4 bp short of the rRNA genes on almost every record, which table2asn flags. The `[organism=]` field uses the raw mapping label, not the name the user confirmed in the Taxonomy tab.
3. **A new user's first hour can fail in several places.** `remotes::install_github()` cannot find the three Bioconductor packages. `new_project()` and `new_test_project()` crash outside RStudio with the default `Rproj = TRUE`. The default executor is conda, which is unsupported and whose `environment.yml` lacks `jq`, `sqlite3`, and R. The default seed URLs 404 while the repo is private. On a laptop, one out-of-memory retry can end the whole local run.
4. **Taxonomy can quietly pick the wrong organism.** Homonyms (for example "Morus", a bird genus and a plant genus) resolve to the lowest taxid, usually the wrong kingdom. Correcting a sample's taxon does not update the taxid lineage the screen uses. A slow taxdump download leaves a partial file that is then trusted forever. Clade databases are searched by name, so homonym clades pull in another organism's rDNA.
5. **The app has a few sharp edges.** Clicking the Annotation or Curation options cell on any multi-unit sample crashes the session. Stop, closing the tab, or refreshing the page sends SIGKILL to Nextflow, which leaves Docker task containers running and writing to the database. Failed runs get a green "Done". Every refresh resets the table to page 1, and `custom.js` re-triggers itself about 20 times a second.
6. **The chimera detector mostly reports noise on clean data.** 15 of 22 clean fish units are labelled "suspect", while the made-up fish plus fungus chimera is also only "suspect". The soft-clip channel can never fire (bowtie2 runs end-to-end), depth dips from repeats look like junctions, and a conserved 18S that resolves only to phylum is counted as a lineage conflict.
7. **About 2,600 lines of code can go**, plus about 3,000 lines of tests that only exercise `dev/` scripts, and 5 dependencies (see the ponytail section).

---

## 2. Counts

Merge script output for this run: `applied 105 Counter({'confirmed': 211, 'unverified-low': 117, 'refuted': 3, 'uncertain': 2})`.

Findings by final severity and status (333 findings from 9 reviewers):

| Severity | Confirmed | Uncertain | Unverified-low | Refuted | Total |
|---|---|---|---|---|---|
| Critical | 0 | 0 | 0 | 0 | 0 |
| High | 20 | 0 | 0 | 0 | 20 |
| Medium | 76 | 0 | 0 | 2 | 78 |
| Low | 90 | 2 | 94 | 0 | 186 |
| Nit | 25 | 0 | 23 | 1 | 49 |
| **Total** | **211** | **2** | **117** | **3** | **333** |

Verifiers lowered 72 severities and raised 1. Four findings came in as critical and were lowered to high (C-nextflow-1, C-data-1, C-refdata-1, C-refdata-2), mostly because the user has to take a deliberate step, or because the first failure does show an error. C-ponytail-39 was raised from low to medium.

After merging duplicates, the 330 findings that are not refuted become 217 distinct items: 17 high, 55 medium, 105 low, and 6 nit items in sections 4 and 6, plus 34 ponytail items in section 5 (all low or nit).

---

## 3. Fix first

The item codes point to the full write-ups: APP, PIPE, SCI, DB, REF, SEC, and PKG codes are in section 4, P codes in section 5, and TEST codes in section 6.

### High (all confirmed)

1. **PKG-H1. A clean install fails.** No `biocViews` field in DESCRIPTION, so `remotes::install_github()` cannot find Biostrings, DECIPHER, or pwalign (C-pkg-1). One line to fix.
2. **PKG-H2. `new_project()` and `new_test_project()` crash outside RStudio** with the default `Rproj = TRUE`, including the README quick start (C-data-3, C-pkg-2). One line to fix.
3. **APP-H1. The Annotation or Curation options cell crashes the session** for any sample with more than one unit (C-shiny-1).
4. **PIPE-H1. Re-assembly loads old path files** left in `out/<ID>/assemble` as if they were new (C-nextflow-1, C-data-5).
5. **DB-H1. Re-assembly keeps a locked unit's old gene coordinates** on the new sequence, and export writes them (C-data-1, C-tests-3).
6. **SCI-H1. The contamination screen verdict is never recomputed** when BLAST hits or the sample's taxon change (C-science-3, C-tests-4).
7. **REF-H1. Changing a sample's taxon does not update `lineage_taxids`**, so the screen and grouping keep using the old organism (C-refdata-1).
8. **DB-H2. Truncated genes are exported as "complete sequence"** with no partial markers (C-data-2, C-science-7).
9. **SCI-H3. Path scores land on the wrong path** when a sample has 10 or more paths, because paths are numbered by file-listing order (C-science-10).
10. **SCI-H2. Circular units are never rotated**, so a start point inside 28S or 18S gives a broken model that trim cannot fix (C-science-5).
11. **REF-H2. A partial taxdump download is kept**, and the index built from it is trusted forever (C-refdata-2, C-pkg-4, C-security-18).
12. **REF-H3. Homonym taxon names resolve to the lowest taxid**, usually the wrong kingdom, with only a grey badge as a hint (C-refdata-3).
13. **PIPE-H2. On the local executor, a memory request larger than the machine ends the whole run**, and the out-of-memory retry makes that likely (C-nextflow-7).
14. **PIPE-H3. Cluster profiles made by `generate_config()` cannot run assembly-only projects** because they always mount `rawDir`, which is `NA` (C-nextflow-3, C-data-10, C-pkg-9).
15. **PIPE-H4. A supplied assembly with more than 5,000 contigs aborts the whole run** with a file-name collision (C-nextflow-2).
16. **PIPE-H5. `run_fallback()` re-assembles good samples** that have one fragment scaffold, and the vignette tells users to run it (C-pkg-5, C-nextflow-22, C-ponytail-3).
17. **APP-H2. A second trim reuses the old trim span**, so "Annotated + flank" can delete 18S and ITS1 (C-shiny-5).

### Medium items with the widest reach

1. **APP-M1. Stop, tab close, and page refresh SIGKILL Nextflow** and leave task containers writing to the database (C-shiny-6, C-nextflow-8, C-security-6, C-pkg-14). Hits every user who stops a run.
2. **APP-M2. Failed runs get a green "Done: Workflow run finished."** (C-shiny-3, C-nextflow-12).
3. **PIPE-M1. Changed BLAST settings never take effect** because the "already searched" check looks only at the sequence (C-nextflow-5, C-refdata-14).
4. **PIPE-M2. Reference count and build runs break `-resume` for the main pipeline** because all runs share one launch folder (C-nextflow-9, C-shiny-32).
5. **PIPE-M3. Changing annotation options re-assembles the sample** and undoes curator trims (C-nextflow-4).
6. **SCI-M5. Re-assembly wipes triage decisions**, so a curator's hand drop comes back as "keep" (C-science-8, C-tests-2).
7. **DB-M4. Export and pending-count code opens the database with no busy timeout**, so it fails, and can end the app session, while the pipeline writes (C-shiny-9, C-data-13).
8. **PIPE-M13. A space in a project, read, or reference path breaks every container task** (C-security-5).
9. **DB-M2. Export group names are used as folder names unchecked**: `..` deletes `out/`, `/` crashes halfway (C-data-7, C-security-1, C-pkg-8).
10. **PIPE-M11. The default executor is conda, and `environment.yml` lacks `jq`, `sqlite3`, R, and ribotyper** (C-nextflow-20, C-pkg-30).
11. **APP-M5 and APP-M6. Every refresh resets the main table to page 1**, and `custom.js` rescans the table about 20 times a second while idle (C-shiny-7, C-shiny-8).
12. **SCI-M1 to SCI-M4. Chimera labels are mostly noise on clean data** (C-science-1, C-science-2, C-science-11, C-science-12).
13. **DB-M5. ITS spacers do not touch the rRNA genes** on almost every record; table2asn warns on each one (C-data-17, C-science-16).
14. **SEC-M1. FASTA header templates are run as R code**, including templates stored in the project database (C-security-3, C-data-6, C-pkg-13).
15. **APP-M14. Curation overwrites the curator's typed Notes** (C-ponytail-1).

---

## 4. Findings by area

Each item gives: code, title, ids, severity, status, where, the problem, why it matters, and the fix. Low and nit items use a shorter one-paragraph form. Where an item says "see", the other item holds the full write-up.

### 4.1 Shiny app

#### High

**APP-H1. The Annotation or Curation options cell crashes the session for any multi-unit sample**
- Ids: C-shiny-1. Severity: high. Status: confirmed (two verifiers).
- Where: `R/app_opts.R:384` (lookup), `:269` and `:301` (modal builders), `R/utils.R:9` (`%|NA|%`).
- Problem: the lookup joins `annotate` (one row per unit) to the options table with no `LIMIT`, so a sample with 3 units returns 3 rows. The modal then runs `cur$rrna_tool %|NA|% ...` on a length-3 vector, and `%|NA|%` is a bare `if (is.na(x))`, which errors.
- Why it matters: the error happens inside an observer with no `tryCatch`, so Shiny closes the session (grey screen). Multi-unit samples are common: 3 of 12 in the test project.
- Fix: add `ORDER BY x.path, x.scaffold LIMIT 1` to the lookup (the Assemble table subquery already does this), or take `cur[1, ]`.
- Status: fixed on `review/2026-09-full-review`. `new_project()` and `new_test_project()` default to `docker`, and `executor = "conda"` warns that it is not supported. The C-pkg-30 ribotyper part is not addressed. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- Status: fixed on `review/2026-09-full-review`. `config.docker`, `config.local`, and `config.conda` carry `resourceLimits = [cpus: <<HOST_CPUS>>, memory: <<HOST_MEMORY>>.GB]`, which `new_project()` fills from `ps`. With Nextflow 25.10.6, a 1,000 GB, 64-CPU request now runs capped at the limit instead of stopping the run. Cluster profiles are unchanged. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- Status: fixed on `review/2026-09-full-review`. The RStudio step runs only when `rstudioapi::isAvailable()`; the Rproj help text says it has no effect outside RStudio. Test: `test-init-reads.R`. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- Status: fixed on `review/2026-09-full-review`. DESCRIPTION has an empty `biocViews:` line, and `remotes::dev_package_deps()` now resolves Biostrings, DECIPHER, and pwalign from Bioconductor (it returned NA before). See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).

**APP-H2. A second trim reuses the old span and length, so a flank trim can delete 18S and ITS1**
- Ids: C-shiny-5. Severity: high. Status: confirmed.
- Where: `R/app_server.R:2158` (`trim_span`), `:2217-2220` (`trim_apply`), `R/pipeline_steps.R:487, 513-517` (`trim_assembly`).
- Problem: `trim_span` depends only on `sel_unit()`, which a trim does not change, so after one trim it still offers the old coordinates and length. Reproduced: trim 6,310 bp down to 4,161 bp, reopen Trim, pick "Annotated + flank", and the offered range is 2050 to 6310 on a 4,161 bp sequence. `trim_assembly()` clamps the end and drops every feature that starts before 2050.
- Why it matters: a curator who trims tightly and then tries to add back some flank deletes 18S, ITS1, and 5.8S instead. Nothing warns. Restore recovers the sequence but throws away any boundary edits.
- Fix: read `edit_tick()` inside `trim_span` (and anywhere else that caches unit length), and check the range against the current length when Apply is pressed.

#### Medium

**APP-M1. Stop, tab close, and page refresh SIGKILL Nextflow, leaving task containers running**
- Ids: C-shiny-6, C-nextflow-8, C-security-6, C-pkg-14. Severity: medium. Status: confirmed (all four).
- Where: `R/app_server.R:1022-1033` (`kill_nf` and `onSessionEnded`), `R/app_reference.R:973-983` (same for reference builds).
- Problem: all of these call processx `p$kill()`, which is SIGKILL to the Nextflow JVM only. Nextflow never gets to run its shutdown, so it cannot stop its Docker containers or cancel scheduler jobs. Three reviewers reproduced it: after `p$kill()` the task containers were still "Up". A page refresh ends the Shiny session, so it triggers the same kill, even though `app_lifecycle.R` is written to survive refreshes.
- Why it matters: after the app says the run stopped, GetOrganelle, BLAST, and R-step containers keep running and keep writing to `.sqlite`. If the user starts a new run, two sets of writers race on the same samples and work folders.
- Fix: send SIGINT first (`p$interrupt()`), wait a few seconds for Nextflow to exit, then fall back to `p$kill_tree()`. Do not kill from `onSessionEnded` before the lifecycle grace period has passed. Fix it once in a shared helper (see ponytail item P1) so the main run and the reference build do not drift.

**APP-M2. Failed Nextflow runs are announced as a success**
- Ids: C-shiny-3 (medium), C-nextflow-12 (low). Severity: medium. Status: confirmed.
- Where: `R/app_server.R:1010-1017`, `R/app_reference.R:962-968`, `R/app_notify.R:17`.
- Problem: when the process ends, the app always shows `rp_alert("Workflow run finished.", type = "message")`, which maps to a green "Done". It never reads `p$get_exit_status()`, and it never counts the tasks that failed and were ignored (`failed_task_errors()` is only used for count runs).
- Why it matters: a config error, a missing container, or every assembly running out of memory all end with a green "Done", while samples stay queued. The progress panel underneath does show "failed: N" and Nextflow's error text, which is why this is medium, not high.
- Fix: check the exit status. After each run, parse the log for "terminated with an error" notes and show the count and the first few errors in the alert.

**APP-M3. "Run from App" crashes the session when Nextflow is missing or the command cannot be built**
- Ids: C-shiny-2. Severity: medium. Status: confirmed.
- Where: `R/app_server.R:982-998` (`nf_start`).
- Problem: `processx::process$new("nextflow", ...)` and `nf_command()` run with no `tryCatch`. `run_app()` only warns at startup when Nextflow is missing. The reference build already does this right (`nf_spawn()` plus `tryCatch`, `R/reference_build.R:378`, `R/app_reference.R:927-937`).
- Why it matters: RStudio opened from the dock often has no conda PATH. The user clicks Run from App and the whole session dies with no message.
- Fix: wrap the launch in `tryCatch`, show an alert with the reason, restore the Start button, and reuse `nf_spawn()` for the PATH check.

**APP-M4. Picking another existing parameter set while Edit is ticked overwrites that set with the values on screen**
- Ids: C-shiny-4. Severity: medium. Status: confirmed.
- Where: `R/app_opts.R:434-442` (name observer), `:568-579` (`opts_apply`).
- Problem: the name box only reacts to unknown names. Choosing a known set does not reload its values, so the fields still show the first set. Save then upserts the chosen name with those values. Reproduced: opened on "actinopteri", ticked Edit, picked "default", saved; "default" now had actinopteri's seeds path and 64 GB memory.
- Why it matters: every sample on the overwritten set gets different GetOrganelle flags, memory, and seed and label databases on the next run. The confirm step only fires when other samples share the set.
- Fix: on a known name, reload every field from that set and untick Edit (the MitoPilot pattern). On an unknown name, keep the values and tick Edit.

**APP-M5. Every refresh resets the main table to page 1 and clears the selection**
- Ids: C-shiny-7. Severity: medium. Status: confirmed.
- Where: `R/app_server.R:421-423`.
- Problem: the refresh observer calls `updateReactable("tbl", data = ...)` with data only. When data is given, reactable resets page, selection, and expanded rows (checked in the installed source). The comment above the code says the page survives; it does not.
- Why it matters: a reviewer on page 3 who clicks Mark reviewed, Lock and Close, K or D, assigns a group, or saves options lands back on page 1 with nothing selected. Background reference polls also bump refresh, so this can happen with no click at all.
- Fix: pass `page = getReactableState("tbl", "page")` and re-select by key, as `set_flag` already does.

**APP-M6. `custom.js` re-triggers itself about 20 times a second while idle**
- Ids: C-shiny-8. Severity: medium. Status: confirmed.
- Where: `inst/app/www/custom.js:12-17` (probe), `:35-36` (observer).
- Problem: `rpReportVisible()` adds a probe element inside `#tbl`, sets its class, and removes it. A body-wide MutationObserver watches exactly those changes and queues another report 50 ms later, so the loop never stops. Each pass calls `getComputedStyle` on every row. Headless Chrome counted 56 calls in 2.8 s of idle time after one DOM change.
- Why it matters: any open project burns CPU continuously, worse with hundreds of rows and on laptops.
- Fix: observe only the `#tbl` subtree and skip mutations on the probe, or call `rpReportVisible` only from reactable's state-change hook and the filter event.

**APP-M7. Annotation-locked units can still be edited, trimmed, and restored**
- Ids: C-shiny-10 (medium), C-security-10 (low). Severity: medium. Status: confirmed.
- Where: `R/app_server.R:2142` (`trim_btn`), `:2217` (`trim_apply`), `:2240-2262` (restore), `:2097` and `:2289` (edit, `nudge_boundary`); `R/app_modals.R:349-385` (controls always shown).
- Problem: only the Reviewed, Problematic, and Partial toggles are disabled for locked units (`:1925`). Trim, Restore, and boundary Edit never check `annotate_lock`, and neither do `trim_assembly()` or `restore_assembly()`. `need_unlocked()` exists in `R/app_guards.R` but is not used here. Reproduced: a locked unit's 18S start moved from 2150 to 2149 and saved.
- Why it matters: a unit signed off for export can change after sign-off with no unlock step, and the next export ships the changed sequence and feature table. MitoPilot makes a locked assembly read-only in this window.
- Fix: hide or disable Edit, Trim, and Restore for locked units, and refuse them on the server with `need_unlocked()`.

**APP-M8. A new set name saved with Edit unticked points samples at a set that does not exist**
- Ids: C-shiny-13. Severity: medium. Status: confirmed.
- Where: `R/app_opts.R:487` (`opts_edit_check`), `:578` and `:586` (`opts_apply`).
- Problem: `opts_apply` only creates the set when Edit is ticked, but it always links the samples to the typed name. Reproduced: typed "remote_only" (Edit auto-ticked), unticked Edit, saved; the samples pointed at "remote_only" and no such row existed.
- Why it matters: `blast.nf` (lines 36-45) joins on the options table, gets nothing back, and silently falls back to shell defaults. The same applies to the other options tables.
- Fix: when the name is not an existing set and Edit is off, refuse with a message, or create the set from the values shown.

**APP-M9. Nothing stops Nextflow from being launched twice**
- Ids: C-shiny-15. Severity: medium. Status: confirmed (by reading; not reproduced live).
- Where: `R/app_server.R:982` (`nf_start`), `R/app_reference.R:915` (`ref_nf_start`).
- Problem: neither handler checks for a live process, and the Start button is hidden only after a server round trip, so a fast double click queues two launches. `nf_process(p)` overwrites the first handle.
- Why it matters: the first run is then never read (its output pipe can fill and stall it) and never killed by Stop or session end. Or the second run dies on Nextflow's lock while the app tracks the dead one and reports "finished".
- Fix: return early when `nf_process()` is alive, and disable the button in the browser on click.

**APP-M10. The Work Dirs browser misses every cluster run**
- Ids: C-shiny-16, C-pkg-10. Severity: medium. Status: confirmed (both).
- Where: `R/app_utils.R:26` and `:29` (`find_workdirs`), `R/app_run.R:114-124` (job scripts).
- Problem: `find_workdirs()` only reads `<project>/.nextflow.log*`, but submitted job scripts write their logs to `nextflow_logs/<base>.nextflow.log`. Also, the regex `Task completed > TaskHandler\[` does not match `GridTaskHandler` lines, so SLURM and SGE runs would show nothing even if that folder were read. (The claim that desktop runs rotate out of the log was judged overstated: Nextflow keeps 9 rotated logs.)
- Why it matters: on HPC, the headless path the app recommends, a curator cannot find a failed task's `.command.err` from the app; it always says "No work directories found".
- Fix: also read `nextflow_logs/*.nextflow.log`, and widen the regex to any `*TaskHandler`.

**APP-M11. The chimera junction hint and track keep pre-trim coordinates**
- Ids: C-shiny-18. Severity: medium. Status: confirmed.
- Where: `R/pipeline_steps.R:479-538` (`trim_assembly`), `R/app_server.R:2198` (`trim_chimera_hint`), `R/app_modals.R:840-845` (`plot_chimera`).
- Problem: `trim_assembly()` shifts annotations and coverage but not `chimera.junctions` or `chimera_windows`, and the hint and the plot read them as stored.
- Why it matters: after trimming 2,149 bp off the 5' end, the Trim modal still suggests a cut at 4400, but the real junction is now at 2251. A curator who cuts at the hint cuts in the wrong place. Design doc 21 says junctions are meant to feed the trim UI.
- Fix: shift the chimera junctions and windows in `trim_assembly()` and `restore_assembly()`, or clear them with a "rerun chimera check" note.

**APP-M12. Bulk triage can silently switch to the whole table, including locked samples**
- Ids: C-shiny-19. Severity: medium. Status: confirmed.
- Where: `R/app_server.R:1676` (`triage_target_ids`), `:1701` (`triage_pv`), `:1768` (`triage_apply`).
- Problem: the target list falls back to every row when nothing is selected, and it is re-read live. A refresh clears the selection (see APP-M5), and refreshes fire while the modal is open when a Nextflow run or a reference build finishes. `auto_triage_sample()` skips annotation-locked units but not assemble-locked ones, and rows hidden by filters are included.
- Why it matters: the user selects 2 samples and opens Triage; a background refresh fires; Apply now changes every sample, while the title still says "2 samples".
- Fix: snapshot the target ids when the modal opens (as the options and taxonomy modals do), and exclude assemble-locked samples.

**APP-M13. The Curation Options "Min anchor identity" and "Flag contaminants" controls do nothing**
- Ids: C-science-17, C-ponytail-2, C-pkg-3. Severity: medium. Status: confirmed (all three).
- Where: `R/app_opts.R:305-308` (controls), `:541-542` (save), `R/curate_rrna.R:255` and `:303` (`resolve_ruleset`, `curate_one`), `R/params_rrna.R:84-87` and `:106`.
- Problem: both values are saved to `curate_opts`, but `curate_one()` builds its settings only from `rrna_params(ruleset)`, where `check_contaminant = TRUE` is hard-coded and the identity floors are per element (fish: 97, 95, and 85). Nothing reads the saved columns. The code comments and `dev/design/08` already list this as a known gap.
- Why it matters: a curator who turns off contaminant flags, or lowers the identity floor for a divergent clade, sees no change, and the modal's default of 85 misstates the real floors. No data is lost, but the UI says a policy is in force that is not.
- Fix: the ponytail choice is to delete both controls, both columns, the seed in `new_db`, and the scalar branch of `anchor_pid_floor` (`R/params_rrna.R:144-155`). The alternative is to merge them into the ruleset as overrides.

**APP-M14. Curation overwrites the curator's typed Notes**
- Ids: C-ponytail-1. Severity: medium. Status: confirmed (two verifiers).
- Where: `R/curate_rrna.R:351-353` (`curate_one`), `R/app_server.R:1801-1802` and `:2265-2284` (Notes box), `R/app_modals.R:325-327`.
- Problem: `annotate.annotate_notes` has two owners. The Notes box reads and writes it, and `curate_one()` replaces it with the joined flag list after every Save, Trim, Restore, or pipeline run. Reproduced: the note "check 3prime end with PCR" became "missing_5.8S; missing_28S; ...".
- Why it matters: the curator's own notes are lost, and when the modal reopens the Notes box shows curation codes as if the user typed them. "Clear notes" also blanks the warnings tooltip while `n_flags` stays above zero. Export does not read notes, so no GenBank data is affected.
- Fix: give curation its own column (for example `curate_notes`), or derive the flag text from `annotations.flags` when needed, and leave `annotate_notes` to the user. The design docs (08 and 11) already propose this split.

**APP-M15. Two opposite rules for locking a second assembly path**
- Ids: C-ponytail-39. Severity: medium (raised from low). Status: confirmed.
- Where: `R/app_server.R:700-728` (bulk Lock), `R/app_modals.R:200-213` (`multi_path_lock_modal`), `R/app_server.R:2342-2350` (Lock and Close).
- Problem: the bulk Lock button offers "Lock anyway" and says both paths will be annotated and exported. Lock and Close refuses the same state with "Only one assembly path can be locked per sample".
- Why it matters: the curator is told two different things about what export will contain, depending only on which button they pressed.
- Fix: pick one rule and use it in both places. Delete the unused `LOCK_UNIT_SQL` at the same time.

#### Low

- **APP-L1. C-shiny-11** (low, confirmed; was medium). `R/app_server.R:435, 622`, `R/tax_search.R:73`. Problem: clicking a locked sample in Taxonomy mode opens the review and lets the user assign a new taxon; the taxonomy lock is not checked. Why: taxid, `sci_name`, and lineage change on a reviewed sample, and anything using them goes stale (the default header uses `{Taxon}`, so `[organism=]` only changes if the template uses `sci_name`). Fix: filter locked ids out of the review targets with `need_unlocked()`.
- **APP-L2. C-shiny-12** (low, confirmed; was medium). `R/app_alerts.R:48`, `R/app_server.R:2103`. Problem: closing the unsaved-edits prompt with its X leaves the edit dirty; re-ticking Edit then overwrites the backup with the nudged coordinates. Why: Reset can no longer undo the change, and `curate_one()` never re-runs on it. Fix: only take a new backup when none exists, or keep Save and Reset visible while edits are dirty.
- **APP-L3. C-shiny-14, C-security-15** (low, confirmed). `R/app_server.R:888, 898-899`. Problem: the cluster job script is built once with `resume = TRUE` and never rebuilt when the "Resume previous run?" box changes; the preview does follow the box. Why: a user who unticks Resume to force a clean run still submits a resuming job. Fix: rebuild the script text when `nf_resume` changes, as the Reference modal already does (`R/app_reference.R:823-829`).
- **APP-L4. C-shiny-17** (low, confirmed; was medium). `R/app_server.R:1088, 1167`, `R/export_opts.R:117`. Problem: once the Export modal has been opened, every refresh re-runs about 4 queries per exportable unit (including full sequences) to validate the header. Why: in projects with hundreds of locked units every click gets slower. Fix: isolate the Save-template toggle from `refresh()`. The verifier notes that the suggested `head(u, 1)` shortcut would break the duplicate-id check (`R/export_opts.R:228`), so keep the full data for that check.
- **APP-L5. C-ponytail-4** (low, confirmed; was medium). `R/db_dist.R:434, 477, 487-492`. Problem: opening BLAST Options can run a synchronous `docker run` (60 s timeout) to describe the panel, and an empty answer is never cached. Why: before the image is pulled, every click on a BLAST options cell can freeze the app for up to a minute. Fix: cache the empty result too, or drop the image probe and show the recorded panel or "Recorded after the first local search".
- **APP-L6. C-ponytail-40, C-data-27** (low, confirmed). `R/export_opts.R:71`, `R/app_server.R:1188-1192`. Problem: `set_export_opts()` turns a failed save into a warning, and the app then always says "Saved header template". Why: if the database is locked during a run, the template is not saved but the user is told it was. Fix: let the error reach the caller and show an alert on failure.
- **APP-L7. C-shiny-20** (low, unverified-low). `R/app_server.R:2047`. Problem: the curation warnings summary depends only on the selected unit, so it is stale after Save or Trim. Why: the flag table shows 0 warnings while the line above still says 2. Fix: add `edit_tick()` to `output$curate_notes`.
- **APP-L8. C-shiny-21** (low, unverified-low). `R/app_server.R:91`. Problem: the status line says "groups" in the Reference Samples view. Why: "Showing 12 of 12 groups" for a list of samples. Fix: use "sample" for that view.
- **APP-L9. C-shiny-22** (low, unverified-low). `R/app_server.R:1245`. Problem: pressing Cancel on "This export already exists" rebuilds the Export modal from saved values. Why: the chosen group and unsaved header are lost, and pressing Export then writes the first group. Fix: pass the current inputs back into `export_modal()`.
- **APP-L10. C-shiny-23** (low, unverified-low). `R/app_server.R:751`. Problem: the toolbar Reviewed, Problematic, and Partial buttons change annotation-locked units, while the modal disables the same toggles. Why: marking a locked unit Problematic from the toolbar silently drops it from the next export while it still shows in the Export table. Fix: skip or refuse locked units in `toolbar_flag()`.
- **APP-L11. C-shiny-26** (low, unverified-low). `R/app_server.R:911`. Problem: "Submit as job" has no duplicate-submission guard (the reference flow has one). Why: two jobs on one project race on the Nextflow lock or the database. Fix: track the pending job and confirm before a second submit.
- **APP-L12. C-shiny-27** (low, unverified-low). `R/app_server.R:1303`. Problem: `unit_by_uid` ignores the sample part of the uid. Why: a late K or D event for sample A can act on sample B's scaffold with the same path and scaffold numbers. Fix: require `p[1] == sel_id()`.
- **APP-L13. C-shiny-28** (low, unverified-low). `R/app_reference.R:1033`. Problem: Reference settings accept a minimum larger than the maximum. Why: Propose groups can never find a clade and gives no hint why. Fix: check `min >= 1` and `min <= max` before saving.
- **APP-L14. C-shiny-29** (low, unverified-low). `R/app_server.R:1299, 1409`. Problem: each K or D click fully re-renders the scaffold table. Why: sort, filters, page, and the ticks used for bulk K/D are lost on every click. Fix: update the scaffold table with `updateReactable()` instead of re-rendering.
- **APP-L15. C-shiny-30** (low, unverified-low). `R/app_server.R:1772`. Problem: Triage Apply always writes the project-wide "default" preset, while the modal says it is for N samples. Why: trying rules on 2 samples changes triage for every later run. Fix: say so in the modal, or save only on an explicit "save as project default".
- **APP-L16. C-shiny-31** (low, unverified-low). `R/app_server.R:684`. Problem: Taxonomy Lock accepts samples with no taxid. Why: unresolved samples move on to Assemble and Reference, where grouping and screening need a lineage. Fix: warn or refuse when taxid is NA.
- **APP-L17. C-shiny-33** (low, unverified-low). `R/app_server.R:979`, `R/app_reference.R:955`. Problem: the progress observer reads `nf_state()` without `isolate()`, so each write re-runs it immediately and bypasses the 250 ms timer. Why: chatty runs re-render the progress text many times a second on the single R process. Fix: use `isolate(nf_state())`.

### 4.2 Pipeline and containers

#### High

**PIPE-H1. Re-assembly loads old path files left in the publish folder**
- Ids: C-nextflow-1 (came in as critical), C-data-5. Severity: high. Status: confirmed (three verifiers).
- Where: `inst/nextflow/modules/pipeline_r.nf:123` (score passes the publish folder), `inst/nextflow/bin/rp_score.R:11`, `R/pipeline_steps.R:24-29` (`score_one`) and `:251` (`load_assemblies`); publishDirs in `assemble.nf:7`, `fallback.nf:11`, `find_rrna.nf:55, 117`.
- Problem: scoring and loading glob `out/<ID>/assemble/<ID>_assembly_*.fasta` and `*.gfa`. `publishDir` with `overwrite: true` replaces files with the same name but never deletes the others, and nothing in `R/` or the `.nf` files clears the folder. Reproduced twice: a first run published paths 1 and 2, a rerun published only path 1, and the database ended with the new path 1 plus the old path 2. Nextflow also documents that publishing is asynchronous, so a downstream task should not read the publish folder at all.
- Why it matters: after a rerun that yields fewer paths (new options, a lower max_paths, or a switch to the bait fallback), old paths are loaded, scored, annotated, and BLASTed as current, and can be auto-selected and exported next to the new assembly. Scoring may also read an old graph file. This is silent. (The common "GetOrganelle failed, then fallback" case writes only the `_0` marker, which is skipped, so that case is safe.)
- Fix: have `rp_score.R` and `load_assemblies()` read the files staged into the score task (the process already receives `path(fastas)`; pass the graph through the channel too), not the publish folder. As a stopgap, delete `out/<ID>/assemble/<ID>_assembly_*` at the start of each assembly.

**PIPE-H2. On the local executor, a memory request larger than the machine ends the whole run**
- Ids: C-nextflow-7. Severity: high. Status: confirmed.
- Where: `inst/config.local`, `inst/config.docker`, `inst/config.conda` (no `resourceLimits`); `inst/nextflow/nextflow.config:63-66` (memory doubles on OOM retry); defaults `R/init_db.R:50` (assemble 24 GB) and `:193` (preprocess 16 GB).
- Problem: when a task asks for more memory than the host has, Nextflow's local executor does not cap it; it stops the session. Reproduced with Nextflow 25.10.6: "Process requirement exceeds available memory", then "Killing running tasks", even with `errorStrategy 'ignore'`.
- Why it matters: on a 16 GB laptop the first assemble task already asks for more than the machine has. On a 32 GB workstation, one sample that runs out of memory retries at double, and that retry kills every other sample's running work.
- Fix: add `process.resourceLimits = [cpus: <host>, memory: <host>]` to the local, docker, and conda templates (supported since Nextflow 24.04), so requests are capped instead of fatal. `new_project()` can fill the host values.

**PIPE-H3. Cluster profiles made by `generate_config()` cannot run assembly-only projects**
- Ids: C-nextflow-3 (high), C-data-10 (medium), C-pkg-9 (medium). Severity: high (two verifiers rated it medium because editing the config by hand works around it). Status: confirmed (all three).
- Where: `R/generate_config.R:32-38` (`container_engine_block`), `R/init_project.R:177-178` (writes `rawDir = 'NA'`), `R/nf_command.R:52` (omits `--rawDir`), `R/reference_nf.R:159-160` (passes the literal `NA`).
- Problem: the generated docker and singularity blocks always mount `${params.rawDir}`. A project with no reads has `rawDir = 'NA'`. With singularity profiles Nextflow fails to even read the config ("Unknown config attribute `singularity.params.rawDir`", reproduced on the pinned versions); with docker every task fails with "mount path must be absolute". The shipped `config.docker` was already fixed for this, but `generate_config()` profiles were not.
- Why it matters: on SLURM, SGE, PBS, or LSF, a project where every sample brings its own assembly fails before any task starts, and so do its reference builds.
- Fix: drop the `rawDir` bind from `container_engine_block()` (new_project already adds the data folder through `EXTRA_MOUNTS`), and treat `'NA'` like empty in `nf_reference_command()`.

**PIPE-H4. A supplied assembly with more than 5,000 contigs aborts the whole run**
- Ids: C-nextflow-2. Severity: high. Status: confirmed.
- Where: `inst/nextflow/modules/find_rrna.nf:29` (scan output), `:59` (pick input), `:14` (chunk size).
- Problem: every scan chunk emits a file named `chunk.gff`, and `groupTuple` hands them all to `find_rrna_pick`, which stages them into one folder. Reproduced: "input file name collision ... chunk.gff", then Nextflow kills all running tasks, despite `errorStrategy 'ignore'`.
- Why it matters: draft genomes with many contigs are exactly the case the chunking was written for, and the failure takes every other sample's tasks down with it.
- Fix: declare the input as `path(gffs, stageAs: 'chunk_*.gff')`, or name each scan output after its chunk.

**PIPE-H5. `run_fallback()` re-assembles good samples and chains preset names**
- Ids: C-pkg-5 (high), C-nextflow-22 (medium), C-ponytail-3 (medium). Severity: high. Status: confirmed (all three).
- Where: `R/run_fallback.R:13-17` (`fallback_candidates`), `vignettes/RiboPilot.Rmd:279, 390` (Step 5).
- Problem: the candidate query joins the per-unit `annotate` rows and picks a sample if any one unit is missing 18S or 28S, with no `DISTINCT` and no ignore filter. Since annotation became per scaffold, any sample with one fragment or ITS-only scaffold qualifies. A sample with two matching units is returned twice, and the loop then creates `default_bait` and `default_bait_bait`.
- Why it matters: a user who follows vignette Step 5 re-queues every unlocked sample whose genes sit on separate scaffolds (a normal outcome), switches it to SPAdes baiting, and replaces a good GetOrganelle assembly. Combined with PIPE-H1, old paths then mix in. The lock step comes after Step 5 in the vignette, so good samples are not yet protected.
- Fix: the author wants runs launched from the app, and the Assembly Options assembler picker already switches a sample to bait or de novo. The ponytail fix is to delete `run_fallback()`, `fallback_candidates()`, and `set_sample_assembler()`, plus the man page and vignette step. If kept, pick samples where no non-ignored unit has both anchors (`NOT EXISTS`), use `SELECT DISTINCT`, and strip an existing `_bait` or `_denovo` suffix before deriving the preset name.

#### Medium

**PIPE-M1. Changed BLAST settings never cause a new search**
- Ids: C-nextflow-5, C-refdata-14. Severity: medium. Status: confirmed (both).
- Where: `inst/nextflow/modules/blast.nf:104-107` (skip test), `:181-188` (failed escalation still recorded), `R/app_opts.R:330-331`.
- Problem: a unit is skipped as "already searched" when its sequence md5 matches `blast_searched.query_md5`. The check ignores the source (local or remote), Entrez query, extra options, `max_target_seqs`, and panel version, and no code ever deletes `blast_searched` rows. Changing BLAST options only re-queues assembly, which usually gives back the same sequence.
- Why it matters: switching to remote BLAST, installing a panel, or changing the Entrez query has no effect, with no message. A unit whose remote escalation failed once is never retried. The old hits stay valid, so nothing is corrupted.
- Fix: store a hash of the effective settings with the md5 and compare both, or clear `blast_searched` for samples whose BLAST options changed. A per-sample "search again" action would also help.

**PIPE-M2. Reference count and build runs break `-resume` for the main pipeline**
- Ids: C-nextflow-9 (medium), C-shiny-32 (low, unverified-low). Severity: medium. Status: confirmed (C-shiny-32 unverified-low).
- Where: `R/reference_nf.R:211` (COUNT), `R/app_reference.R:927` (BUILD), `R/app_server.R:995` (main run); all use the project folder as the launch folder.
- Problem: a bare `-resume` picks the last session in `.nextflow/history`, whatever script it ran. Reproduced with the pinned Nextflow: after run B, `A -resume` picked up B's session and re-ran A's finished task.
- Why it matters: after any reference count or build, re-running the pipeline for re-queued samples starts GetOrganelle and preprocessing from scratch (hours per sample) instead of using the cache.
- Fix: launch reference workflows from `<project>/.reference/` (their own history), or record the main workflow's session id and pass it to `-resume`.

**PIPE-M3. Changing annotation options re-assembles the sample and undoes curator trims**
- Ids: C-nextflow-4. Severity: medium (was high). Status: confirmed.
- Where: `inst/nextflow/modules/pipeline_r.nf:40` (`query_samples` takes states 1 and 4), `R/app_opts.R:589-592`, `R/pipeline_steps.R:224` (`refresh_sample_state`).
- Problem: state 4 means "annotation pending", but `query_samples` treats 4 as "assemble". Changing annotate or curate options sets units back to pending, the sample moves to 4, and the next run re-assembles it while the unit query picks nothing. Reproduced: a trimmed, edited unit came back as the raw sequence.
- Why it matters: trims and undo snapshots are lost, and any ignored annotation failure also keeps a sample at 4, so it is re-assembled on every run. It needs the assembly to be unlocked; the verifier found the locked-unit part of the claim is not reachable from the app.
- Fix: give "annotation pending" its own state (or make `query_samples` select only 1), and have `refresh_sample_state()` ignore locked units.

**PIPE-M4. `-resume` reuses an old fallback assembly after its options change**
- Ids: C-nextflow-6. Severity: medium (was high). Status: confirmed.
- Where: `inst/nextflow/modules/fallback.nf:15`, `inst/nextflow/bin/rp_fallback.R:13-19`.
- Problem: the task's inputs are only id, reads, cpus, memory, and method. The seeds and labels databases and the bowtie2 and SPAdes options are read from the database at run time, so they are not part of Nextflow's cache key.
- Why it matters: after changing SPAdes options or the reference database for a bait-assembled sample, Nextflow reports the task as cached and reloads the old assembly. Unticking Resume works around it.
- Fix: pass seeds and labels as `path` inputs and the tool options as `val` inputs, as the GetOrganelle arm already does.

**PIPE-M5. Max paths 0 crashes GetOrganelle, and the max-paths gate can never fire**
- Ids: C-nextflow-10. Severity: medium. Status: confirmed.
- Where: `inst/nextflow/modules/pipeline_r.nf:26`, `R/app_opts.R:158-163`, `R/pipeline_steps.R:342`.
- Problem: the app allows 0 and its help says "0 means no limit", but the value goes straight to `--max-paths-num 0`, which GetOrganelle rejects with an assertion. Separately, GetOrganelle already cuts its output to max-paths, so `load_assemblies()` never sees more paths than the limit. In the test project, SRR14433904's log says "Only exporting 10 out of all 24 possible paths", and the sample is marked done with no note.
- Why it matters: (a) a user who sets 0 gets a sample that stays queued with no message; (b) a sample with 24 candidate paths is silently cut to 10 instead of being failed by the gate.
- Fix: pass a large number (or omit the flag) when max paths is 0 or less, and let GetOrganelle export max_paths + 1 so the gate can detect overflow, or read the "Only exporting" warning.

**PIPE-M6. A fallback SPAdes crash (including out of memory) exits 0, so the memory retry never runs**
- Ids: C-nextflow-14. Severity: medium. Status: confirmed.
- Where: `R/assemble_fallback.R:64-91` (`run_spades`, `bait_reads`), `:186` (temp folder), `:191-197`.
- Problem: SPAdes and bowtie2 run with their output and exit status thrown away; the wrapper only checks whether `contigs.fasta` exists, then writes the `_0` marker and returns normally. SPAdes also gets no `-m` limit, and its working data goes to the R temp folder.
- Why it matters: on SLURM with memory limits, SPAdes is killed, the task exits 0, it is never retried with more memory, the sample becomes "failed assembly", and `.command.err` has no SPAdes output to explain why. On nodes with a small `/tmp`, SPAdes runs out of space.
- Fix: exit with the tool's non-zero status, send tool output to stderr, pass `-m` from the task memory, and put SPAdes and baited reads under the task folder.

**PIPE-M7. A cached `find_rrna_pick` skips its database updates, so the sample stays queued**
- Ids: C-nextflow-15. Severity: medium. Status: confirmed.
- Where: `inst/nextflow/modules/find_rrna.nf:50-74`, `inst/nextflow/bin/rp_find_rrna.R`.
- Problem: this step writes the sample state, notes, and candidates to the database, and reads `max_scaffolds` from the database, but it is not `cache false` and its inputs do not change on a re-queue.
- Why it matters: re-queue a supplied assembly that found no rDNA (for example after raising max_scaffolds), and with `-resume` the step is cached, nothing updates the row, and the sample shows "queued" forever.
- Fix: add `cache false` to `find_rrna_pick` and `stage_assembly` (both are cheap and write to the database), or pass `max_scaffolds` as an input.

**PIPE-M8. A failed `find_rrna_scan` chunk silently shrinks the search**
- Ids: C-nextflow-16. Severity: medium. Status: confirmed.
- Where: `inst/nextflow/modules/find_rrna.nf:23` (ignore), `:67` (concatenates), `:89` (`groupTuple` with no size).
- Problem: a scan chunk that fails is ignored, and the pick step just combines whatever chunk results arrived.
- Why it matters: if the chunk that held the rDNA contig fails (memory, pre-emption), the sample is recorded as "no rRNA genes found in the supplied assembly", with no warning. Rare, but silent and wrong.
- Fix: retry or stop on scan failure, and have the pick step check it got every chunk (`groupTuple(size:)` or `groupKey`).

**PIPE-M9. Samples below minDepth are silently skipped and stay queued; the depth unit does not match the docs**
- Ids: C-nextflow-17. Severity: medium. Status: confirmed.
- Where: `inst/nextflow/main.nf:91`, `inst/nextflow/modules/preprocess.nf:22`, `R/init_project.R:19`, `R/nf_command.R:167`.
- Problem: the filter drops shallow samples with no note or state change. The count it compares is fastp's total reads (both mates), while the docs say "trimmed paired reads".
- Why it matters: a shallow sample never assembles, shows "queued" after every run, keeps counting as pending, and is re-preprocessed each time with no explanation. A `min_depth` of 2,000,000 actually lets through a sample with 1,000,000 pairs.
- Fix: send low samples to a small step that sets them to failed with a note such as "below minDepth (N reads)", and use one unit (pairs or reads) in both the code and the docs.

**PIPE-M10. The awsbatch config cannot run as written**
- Ids: C-nextflow-18. Severity: medium. Status: confirmed.
- Where: `inst/config.awsbatch:4-21`, `R/nf_command.R:50`, `R/constants.R:36`, `R/generate_config.R:171`.
- Problem: `nf_command()` always passes a local `-work-dir`, and Nextflow's AWS Batch executor needs an S3 work folder. AWS Batch also ignores `docker.runOptions`, so the database mount in the template never happens. The template's own header says pure S3 cannot work.
- Why it matters: awsbatch is still offered as an executor, but no combination of settings runs.
- Fix: remove awsbatch from the offered executors until it is supported (simplest), or template `aws.batch.volumes` and let `nf_command()` honour a configured work folder.

**PIPE-M11. The default executor is conda, and `environment.yml` lacks tools the pipeline calls**
- Ids: C-nextflow-20 (medium), C-pkg-30 (low, unverified-low). Severity: medium. Status: confirmed (C-pkg-30 unverified-low).
- Where: `environment.yml:11` onward, `inst/nextflow/modules/preprocess.nf:22` (`jq`), `R/init_project.R:52` and `new_test_project()` (default `executor = 'conda'`), `R/chimera.R:317`.
- Problem: `environment.yml` has no `jq`, no `sqlite3` command, no R or R packages, and no ribotyper, but the pipeline calls all of them. The author's own env works only because the system has `jq`, `sqlite3`, and R. Yet conda is the default executor, and the README says conda is not supported.
- Why it matters: a user who leaves the executor at its default gets a config for an unsupported route. Without system `jq`, every preprocess task fails, is ignored, and the run still says "finished". A missing ribotyper quietly weakens chimera checks.
- Fix: default `new_project()` and `new_test_project()` to `docker` and drop the README workaround text. If conda stays, add `jq`, `sqlite`, `r-base` plus the R dependencies, and ribovore to `environment.yml`.

**PIPE-M12. Tool versions are not pinned, and the conda and container versions differ**
- Ids: C-nextflow-19 (medium), C-security-22 (low). Severity: medium. Status: confirmed (both).
- Where: `docker/Dockerfile:32-60`, `environment.yml:38`, `inst/nextflow/modules/pipeline_r.nf:164`, `find_rrna.nf:34`.
- Problem: the Dockerfile leaves sqlite, fastqc, BLAST (`>=2.12` only), seqtk, seqkit, pybarrnap, barrnap, ITSx, hmmer, entrez-direct, biopython, and datasets unpinned, and runs the ribovore `install.sh` from a moving tag with no checksum. `environment.yml` pins barrnap 0.9, which has no `fun` kingdom, while the pipeline translates `euk` to `fun` for barrnap.
- Why it matters: rebuilding the same image tag later can change annotation results with no version change. In an env built from `environment.yml`, choosing barrnap fails every annotation task (pybarrnap, the default, is not affected).
- Fix: pin every tool in the Dockerfile to the versions that were validated (the image has barrnap 1.10.5), use the same pins in `environment.yml`, and check the install script's sha256 before running it.

**PIPE-M13. A space in a project, read, or reference path breaks every container task**
- Ids: C-security-5. Severity: medium. Status: confirmed.
- Where: `R/generate_config.R:32-34`, `inst/config.docker:20`, `config.local:18`, `config.awsbatch:21`, `R/init_project.R:161-166` (extra mounts), `R/reference_nf.R:61-64` (`container_bind_opts`).
- Problem: mount paths are pasted into `runOptions` with no quotes, and Nextflow passes `runOptions` to docker as-is, so a space splits the `-v` argument and the next word is read as the image name ("Unable to find image 'test/my:latest'"). Only `validate_blastdb_prefix()` refuses whitespace.
- Why it matters: a Mac user with a project under "Google Drive" or a reference folder on "/Volumes/Lab Share" gets every task failing with a confusing image-name error. See also DB-M3 for the export ITSx check.
- Fix: refuse whitespace and quotes in project, data, assembly, and reference paths in `new_project()` and `validate_ref_dir()`, with a clear message. Quoting every mount is the harder alternative.

**PIPE-M14. Cluster head-job scripts only fit Hydra-style SGE**
- Ids: C-nextflow-23, C-security-13. Severity: medium. Status: confirmed (both).
- Where: `R/app_run.R:27` (`detect_scheduler`), `:82-100` (`job_script`), `R/hydra.R:21-24`.
- Problem: any `qsub` system gets SGE `#$` lines with a hard-coded `-pe mthread 2` (a Hydra name) and `h_vmem` with no `NXF_OPTS`, so PBS ignores them all. The SLURM and LSF log names render as `%%j` and `%%J`, so the scheduler writes a literal `%j`. `is_hydra_cluster()` matches "hydra" anywhere in `/etc/hosts`, so another cluster named Hydra gets the NMNH script.
- Why it matters: on PBS the Nextflow head job runs with the queue's default walltime (often 1 hour) and is killed mid-run. An SGE site without an `mthread` environment rejects the job. The script is editable before submit, which only partly helps.
- Fix: tell PBS and SGE apart and write `#PBS` lines for PBS, take the parallel environment from the saved profile, export `NXF_OPTS='-Xms512m -Xmx4g'` in every script, match Hydra by host name, use single `%` in the log patterns, and quote the log path.

**PIPE-M15. Container and package version mismatch goes unnoticed, and every task prints a developer warning**
- Ids: C-nextflow-25. Severity: medium. Status: confirmed.
- Where: `inst/nextflow/bin/rp_env.R:53-59`, `R/nf_command.R:150-153`, `R/init_project.R:53`.
- Problem: for normal installs, `nf_command()` passes the installed package folder as `--ribopilot_src`, which `rp_env.R` rejects, so every task logs "Uncommitted package edits are NOT in this run ... rebuild the image". The project's container tag is fixed when the project is made, and nothing compares the image's RiboPilot version with the package's.
- Why it matters: after a package upgrade, the app writes new-format data while every task runs the older image code against the same database (and there are no migrations pre-release). Meanwhile every end user's logs carry a misleading warning.
- Fix: pass the package version to the tasks and stop (or warn loudly) in `rp_env.R` when the image version differs. Pass `--ribopilot_src` only for a real source checkout.

**PIPE-M16. With no usable reference, the bait fallback assembles the whole genome and keeps every contig**
- Ids: C-nextflow-13 (low), C-tests-8 (medium). Severity: medium. Status: confirmed (both).
- Where: `R/assemble_fallback.R:118-120` (`extract_rdna_contigs`), `:179-181` (`bait_reads`), `inst/nextflow/main.nf:53, 100-110`.
- Problem: sending seedless samples to the bait method is intended. But when the labels database is blank or unreachable, `extract_rdna_contigs()` keeps every contig, and `bait_reads()` falls back to all reads. None of the fallback assembler functions are tested.
- Why it matters: SPAdes runs on the full genome-skim read set for hours and every contig of 300 bp or more is written as an "rDNA" scaffold. The default `max_scaffolds = 10` gate usually fails the sample afterward; with the gate at 0, thousands of contigs go on to BLAST and annotation. The app blocks missing file databases, so this needs a blank or dead database.
- Fix: when the rDNA filter cannot run, keep nothing and write a clear note (fail closed). Add tests for `fallback_assemble_one()` and `extract_rdna_contigs()` with fake `bowtie2`, `spades.py`, and `blastn` on PATH.

**PIPE-M17. BLAST database provenance is never recorded in container runs**
- Ids: C-tests-1. Severity: medium (was high). Status: confirmed.
- Where: `inst/nextflow/bin/rp_searched.R:22`.
- Problem: the script calls `record_blastdb_provenance()`, which is not exported, without `RiboPilot:::`, inside `tryCatch(error = function(e) NULL)`. From a source checkout `load_all()` exposes it, so it works for the author. In the Docker or Singularity image `library(RiboPilot)` is used, so the call fails and the error is swallowed. `R CMD check` fails on exactly this (test-nextflow.R:1051), but CI skips that test because its runner has no Nextflow.
- Why it matters: for discovered or baked-in local panels, the project never records which panel version it searched. Panels configured at launch are still recorded, so the loss is display-only provenance, but the bug class (bare calls to internals from bin scripts) can hide worse failures.
- Fix: call `RiboPilot:::record_blastdb_provenance`, and log instead of swallowing the error. Extend the existing check in `test-reference-workflow.R:293` to scan every `inst/nextflow/bin/*.R` for bare calls to unexported functions.

#### Low

- **PIPE-L1. C-nextflow-11** (low, uncertain; was medium). `inst/nextflow/modules/pipeline_r.nf:188-259`, `R/pipeline_steps.R:171`. Problem: the per-unit database writers (reconcile, curate, chimera, triage, screen) use `errorStrategy 'ignore'` with no retry, and the busy timeout is 10 s. Why: on a cluster with many units on NFS, a "database is locked" error would leave a unit silently unannotated; the verifier found no such error in the author's test project, so harm is unproven. Fix: retry these steps with backoff, raise the busy timeout to 60 to 120 s, and cap concurrent database writers (`maxForks`, or run the tiny R steps on the local executor).
- **PIPE-L2. C-nextflow-24** (low, confirmed; was medium). `R/generate_config.R:37-46`, `inst/config.slurm:27`. Problem: generated singularity profiles set no `pullTimeout` and no shared `cacheDir`, and every small R step inherits the 16 GB per-task memory default. Why: the 7.6 GB image can hit the 20 minute pull timeout on first use (only `config.NMNH_Hydra` sets 60 min), each project re-pulls it, and hundreds of 16 GB jobs are submitted for sub-GB R scripts. Fix: write `pullTimeout = '60 min'` and a shared `cacheDir`, and give the R-script steps a small-memory label.
- **PIPE-L3. C-security-8** (low, confirmed; was medium). `R/generate_config.R:~60, 234`, `R/init_project.R:170-183`. Problem: values put into Groovy config strings (queue, walltime, cluster options, paths) are not escaped. Why: a cluster option with a quote, or a path like `/data/O'Neil_lab/`, writes a config Nextflow cannot parse; the failure is loud but cryptic. Fix: escape quotes and backslashes (and `$` in run options), or refuse those characters with a clear message; optionally run `nextflow config` once on the written file.
- **PIPE-L4. C-security-19** (low, confirmed). `R/run_assembly.R:35` (`config_param`). Problem: the config reader cuts everything after `//` and strips all quotes. Why: `rawDir = 's3://bucket/raw'` is read as `s3:`, and an apostrophe in a path is lost; `nf_command()` passes the result on as `--rawDir`. Fix: take the quoted value first, then strip a trailing comment only outside quotes.
- **PIPE-L5. C-security-23** (low, confirmed). `R/assemble_fallback.R:32` (`bt2_args`). Problem: bowtie2 baiting writes a full SAM file to the R temp folder that nothing reads. Why: on a deep skim this can fill a small `/tmp` on a compute node. Fix: use `-S /dev/null`.
- **PIPE-L6. C-nextflow-42 (unverified-low), C-security-26 (confirmed)** (low). `R/nf_command.R:163`, `R/nextflow_version.R:141`, `R/app_reference.R:823-829`. Problem: `nf_pending_count()` re-implements the pipeline's sample queries (and already differs from them), and `nf_pin_version()` starts a JVM (`nextflow -version`) on every call although the version is cached per session. Why: the Update modal can count samples the pipeline will not run, and each launch or Reference settings change pays 0.5 to 3 s (Shiny's text debounce softens the per-keystroke case). Fix: keep the query SQL in one place, and use the session-cached version in `nf_pin_version()`.
- **PIPE-L7. C-nextflow-27** (low, unverified-low). `R/nf_command.R:62`. Problem: the command-line `-process.containerOptions` for reference mounts replaces any `containerOptions` the user set in `.config`. Why: options like `--shm-size` or extra binds silently vanish once a project uses an outside reference folder. Fix: add the reference mounts through a `-c` overlay on `runOptions`, or merge with the existing value.
- **PIPE-L8. C-nextflow-28** (low, unverified-low). `inst/nextflow/modules/pipeline_r.nf:124`. Problem: score reads `projectDir/out`, while the assemble steps publish to `launchDir/out`. Why: launching from any other folder, or setting an absolute `publishDir`, makes score find no files and fail every sample. Fix: use one base for both; this goes away if score reads staged inputs (PIPE-H1).
- **PIPE-L9. C-nextflow-29** (low, unverified-low). `inst/nextflow/modules/find_rrna.nf:128`. Problem: `stage_assembly` runs `awk` directly on the supplied file, so a gzipped assembly becomes binary garbage. Why: a sample with `contig.fa.gz` and find_rrna off fails. Fix: read through `zcat -f` (or seqkit) first.
- **PIPE-L10. C-nextflow-30** (low, unverified-low). `inst/nextflow/nextflow.config:65`. Problem: exit code 140 (a walltime kill on SGE and LSF) is in the out-of-memory retry list, and retries only raise memory. Why: an assembly that hits walltime is resubmitted up to 3 times with the same time limit. Fix: drop 140 from the list, or also scale time with the attempt.

### 4.3 Science and algorithms

#### High

**SCI-H1. The contamination screen verdict is never recomputed when hits or the sample's taxon change**
- Ids: C-science-3, C-tests-4. Severity: high. Status: confirmed (both).
- Where: `R/screen.R:329` (skip test), `:367` (the only `DELETE FROM screen`), `inst/nextflow/bin/rp_screen.R:13`.
- Problem: `screen_one()` returns the stored verdict whenever it is at the current `SCREEN_VERSION` (unless it was `no_reference`), and the pipeline never passes `force`. Nothing deletes screen rows when BLAST hits are replaced or a sequence changes. Reproduced twice: fish hits gave `on_target`; after replacing them with fungal hits (what a re-search stores) the verdict stayed `on_target`, and only `force = TRUE` gave `off_target`. Changing the sample's lineage to a starfish also left both units `on_target`.
- Why it matters: after a re-assembly turns scaffold 1.1 into a different contig, or after the user corrects a mis-identified taxon, the screen, the `screen_off_target` triage rule, and export all act on the old verdict. REF-H1 makes the taxon case worse, because the lineage the screen reads is not updated either.
- Fix: store what the verdict was computed from (for example the unit's `query_md5` plus the sample taxid) and recompute when it differs. The simpler option is to delete the unit's screen row whenever its hits are replaced or its sequence changes, and when the sample's taxon changes. Add a regression test that does not pass `force`.
- Status: partly fixed on `review/2026-09-full-review`. A re-assembly that changes a unit's sequence now deletes its screen row (`reset_unit_evidence()` in `R/pipeline_steps.R`), so the verdict is recomputed. A taxon change still leaves the old verdict. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).

**SCI-H2. Circular units are never rotated, so a start point inside 28S or 18S gives a broken model**
- Ids: C-science-5. Severity: high. Status: confirmed.
- Where: `R/pipeline_steps.R:376` (orientation only), no rotate function anywhere; design `dev/design/02_ASSEMBLY.md:111`, `01_ARCHITECTURE.md:124`.
- Problem: GetOrganelle labels some units circular and picks the start point by sorted vertex name, which has nothing to do with 18S. The `topology` value is stored but never used, and the planned rotate-to-18S step was never built. The reviewer rotated a real unit (SRR14433901 1.1) to start inside 28S and ran the real tools: the model became a 3' fragment of 28S, then 18S, ITS1, and 5.8S, with ITS2 and the 5' half of 28S lost.
- Why it matters: for any circular unit whose start falls inside 18S, ITS, or 28S (roughly half the unit length), the annotation is split and wrong. Trim can only remove sequence, so the complete unit can never be exported correctly. The current test units happen to start in the IGS.
- Fix: for circular units, rotate before annotation so the sequence starts in the IGS (for example a fixed offset upstream of the 18S found by pybarrnap), and add a rotate action in the app as MitoPilot has.

**SCI-H3. Path scores are assigned by file-listing order, so with 10 or more paths they land on the wrong path**
- Ids: C-science-10. Severity: high. Status: confirmed.
- Where: `R/select_paths.R:123-124` (`path = i`), `R/pipeline_steps.R:24` (`Sys.glob`), `R/app_utils.R:191` (join on path).
- Problem: scoring numbers paths by their position in `Sys.glob()` output, not by the number in the file name. In the host locale the order is `_10, _1, _2, ...`; in the C locale it is `_1, _10, _2, ...`. Real sample SRR14433904 has exactly 10 paths (the default cap), all tied at 0.948, and file `_10` was scored as path 1 and auto-selected. `format_assembly.sh` already fixed the same bug for file names but not here.
- Why it matters: path score, depth, anchors, best hit, and "Selected Path" are shown against the wrong path, and the result differs between the host and the container.
- Fix: read the path number from the file name (`sub('.*_assembly_([0-9]+)\\.fasta$', '\\1', f)`) or from the FASTA header.

#### Medium

**SCI-M1. The soft-clip chimera channel can never fire**
- Ids: C-science-1. Severity: medium (was high). Status: confirmed.
- Where: `inst/nextflow/modules/coverage.nf:47`, `R/chimera.R:237`.
- Problem: coverage maps reads with bowtie2 in end-to-end mode, which never soft-clips. All 24 coverage rows in the test project (including the made-up chimera CHIM01) have zero clips. A direct test gave 0 soft-clipped reads end-to-end and 64 with `--local`.
- Why it matters: design doc 21 describes this channel as the only one that can reach "likely" on its own, and it is the one that sees same-host misjoins the taxonomy channels cannot. The other channels still run, so this is a dead signal, not a crash.
- Fix: count clips from a `--local` alignment (a second bowtie2 pass, or `--local` for the coverage map with `-q 10` kept), then re-check the 0.3 threshold on CHIM01-style fixtures.

**SCI-M2. The depth-step channel fires on normal IGS and ETS repeats**
- Ids: C-science-2. Severity: medium (was high). Status: confirmed.
- Where: `R/chimera.R:219`, `inst/nextflow/modules/coverage.nf:48` (`samtools view -q 10`).
- Problem: 15 of the 22 clean fish units are "suspect" from depth steps alone. The dips sit on short repeats inside the IGS and ETS: reads there map to two places, get a mapping quality under 10, and are dropped, which looks like a junction. The real cross-kingdom chimera CHIM01 is also only "suspect", so the label does not separate the two.
- Why it matters: "suspect" only feeds a warn-only export note (the triage flag needs "likely"), so impact is limited, but a depth step near another channel's junction can push a clean unit to "likely", and the export risk flags already treat "suspect" as chimeric.
- Fix: compute steps on a depth track without the mapping-quality filter, or mask positions that occur twice on the scaffold before calling steps, and require the step to persist rather than be a single-window dip.

**SCI-M3. The gene channel counts a shallow (unresolved) element as a lineage conflict**
- Ids: C-science-11. Severity: medium. Status: confirmed.
- Where: `R/chimera.R:68` (`chimera_channel_genes`).
- Problem: each element's lineage is cut below its LCA rank, and `rank_distance()` then reports the depth of the shallower one as a conflict. Reproduced on real data: 18S resolves only to Chordata while 5.8S resolves to Pleuronectiformes, both with the same top species, and the unit scores 0.6 "18S resolves to Chordata but 5.8S to Pleuronectiformes".
- Why it matters: conserved 18S or 5.8S often resolves only to phylum or class while ITS or 28S reaches genus, so clean samples get a 0.6 "suspect" call with a misleading reason.
- Fix: compare only ranks where both elements have a taxid, and fire only when they disagree at class or above at a rank both resolve. An ancestor is not a conflict.

**SCI-M4. The window channel picks an arbitrary best hit among ties in conserved regions**
- Ids: C-science-12. Severity: medium. Status: confirmed.
- Where: `R/chimera.R:99-100`.
- Problem: the best hit per window is the first row after sorting by bitscore, so ties are broken by BLAST output order. Real example: a 500 bp window inside fish 18S came back as "Mammalia (PREDICTED: Canis, 100.0%)" next to Actinopteri windows, a class switch that scored 0.5 on both units.
- Why it matters: 500 bp windows in conserved 18S and 28S cores match many classes at 100 percent, so the class is arbitrary and clean scaffolds get a window-channel chimera signal.
- Fix: per window, take the LCA of all hits tied at the top bitscore (or within a small margin), and treat windows whose LCA is above class as uninformative.

**SCI-M5. Re-assembly wipes triage decisions, so a curator's hand drop comes back as "keep"**
- Ids: C-science-8, C-tests-2. Severity: medium (both were high). Status: confirmed (both).
- Where: `R/pipeline_steps.R:384-402` (carry-over), `:348-369` (gate path), `R/triage.R:129`.
- Problem: `load_assemblies()` rewrites the sample's `assemblies` rows and carries over only `ignore`, matched on the path and scaffold numbers. The triage, reason, source, and rule columns come back empty. `auto_triage_sample()` skips only rows whose source is "user", so the carried-over drop now looks undecided and is re-decided. Reproduced: a user drop went to "ignore 1, no source" after re-assembly, then to "keep, auto, ignore 0" after triage. Matching by number also means that if contigs come back in a different order, the drop lands on the wrong scaffold. The code comment promises to carry curator decisions.
- Why it matters: a scaffold the curator rejected goes back into Annotation Review with no record it was ever rejected, and can be locked and exported. Existing tests miss it because they set `ignore = 1` directly instead of going through `set_unit_triage()`.
- Fix: carry all four triage columns with `ignore`, and do not discard user decisions on the gate path. Better, key the carry-over on a hash of the sequence (see DB-H1). Add a test: user drop, then `load_assemblies()`, then `auto_triage_sample()`.

**SCI-M6. The feature model takes the first hit per gene, so a partial leading copy replaces the real anchor**
- Ids: C-science-6, C-tests-6. Severity: medium (C-science-6 was high). Status: confirmed (both).
- Where: `R/annotate.R:105` (`pick()`), `:170-180` (`clip_spacer`).
- Problem: `pick()` returns the first row per element, in GFF start order. Reproduced twice: a scaffold with a 28S tail at the start and a full 28S later gave a model of 28S (the tail), 18S, ITS1, and 5.8S; the full 28S and ITSx's ITS2 were dropped.
- Why it matters: user contigs with several rDNA copies, rotated circles, and scaffolds that start with a partial repeat get a model built from the wrong copy, losing ITS2 (the main barcode). Curation flags `missing_ITS2` and `out_of_order`, so it is not silent, but the curator is never offered the right model.
- Fix: pick anchors relative to the chosen 5.8S or ITS: the nearest 18S upstream and the nearest 28S downstream on the same strand. Add a test with two partial 28S copies. Rotation (SCI-H2) removes the common circular case.

**SCI-M7. The isomer and paralog classes are computed before annotation finishes**
- Ids: C-science-4. Severity: medium (was high). Status: confirmed.
- Where: `R/screen.R:395` (`screen_variants`), `inst/nextflow/main.nf:170` (screen waits only for BLAST).
- Problem: the screen step, which also runs `screen_variants()`, waits only for BLAST, which runs in parallel with annotation. In the test project, SRR14433904 has 10 identical paths; the screen ran 2 to 8 seconds before 9 of the 10 were annotated, so one unit is stored as "singleton" and nine as NA. Rerunning now gives "isomer" for all 10.
- Why it matters: the curator reads "singleton, the sample's only comparable unit" for a sample with 10 copies. These columns are information only, with no gating.
- Fix: call `screen_variants()` from `auto_triage_sample()` (which already waits for curation and chimera), or make the screen wait for reconcile for that sample.

**SCI-M8. Coverage names like "1.10" are read as the number 1.1, merging two scaffolds' depth tracks**
- Ids: C-science-9. Severity: medium (was high). Status: confirmed.
- Where: `R/coverage.R:65`.
- Problem: `read.table()` with no `colClasses` turns the reference column into numbers. Reproduced: with scaffolds 1.1, 1.10, and 1.2, only 2 rows were written, scaffold 1.1 got the two tracks interleaved, and 1.10 got nothing.
- Why it matters: any path with 10 or more scaffolds (fragmented assemblies, user contigs) gets a corrupted depth track and median depth for scaffold 1 (which also feeds the chimera read channel) and no coverage for scaffold 10. The default gate allows exactly 10.
- Fix: `read.table(..., colClasses = c('character', 'integer', 'integer'))`.

**SCI-M9. A failed seeds download stops everything after scoring for that sample**
- Ids: C-science-18. Severity: medium. Status: confirmed.
- Where: `R/score_project.R:17` (`local_ref`), `R/pipeline_steps.R:21` (`score_one`), `R/init_db.R:102`.
- Problem: `download.file()` throws on an HTTP error, but `local_ref()` only checks a return code, and `score_one()` calls it with no `tryCatch`. Reproduced with the default seeds URL (currently 404).
- Why it matters: the default seeds database is a URL. On a compute node with no internet, or with a dead URL, a successfully assembled sample gets no BLAST, coverage, annotation, or triage, and the only error is in the Nextflow log.
- Fix: wrap the download in `tryCatch` and return `NA` (the callers already handle `NA`), or reuse the seeds file Nextflow already staged for the assembler.

#### Low

- **SCI-L1. C-science-13** (low, confirmed; was medium). `R/annotate.R:44`, `R/curate_rrna.R:117`. Problem: pybarrnap reports "partial" only when the hit is under 80 percent of a fixed 28S length (2,912 bp), and otherwise `parse_barrnap` sets coverage to 100. Why: the 90 percent floor behaves like 80, and a fish 28S (about 3,900 bp) missing up to about 40 percent can pass as complete when it does not touch a scaffold end. Fix: judge completeness from model start and end or a clade expected length, or pass `--lencutoff` to match `min_anchor_coverage`.
- **SCI-L2. C-science-14** (low, confirmed; was medium). `R/annotate.R:202`, `R/app_opts.R:280`. Problem: Perl barrnap 1.10.5 with kingdom `fun` has no 5.8S model and writes no partial notes, so every hit gets coverage 100. Why: a user who picks barrnap loses the truncated flag for every anchor and can get an 18S that runs into an adjacent chimeric 28S. Fix: remove barrnap as a choice or warn on selection, and store coverage as `NA` (not 100) when the tool gives none.
- **SCI-L3. C-science-15** (low, confirmed; was medium). `inst/nextflow/modules/pipeline_r.nf:171`, `R/pipeline_steps.R:95, 377`. Problem: ITSx runs with `--complement T` (sense-frame coordinates) while barrnap coordinates are in the input frame, and the two are combined unmapped. Why: a reverse-oriented unit that was not flipped (a tie, or no barrnap anchors while ITSx still finds spacers) gets ITS1, 5.8S, and ITS2 stored at mirrored positions. Rare. Fix: run ITSx with `--complement F` on the oriented unit, as `annotate_records` does.
- **SCI-L4. C-science-19** (low, uncertain; was medium). `R/screen.R:213-218`. Problem: hits titled like genome records ("chromosome, complete genome") are left out of the taxonomy LCA as well as the rDNA vote. Why: a contaminant that only hits complete genomes could come out as `no_reference`; the verifier notes this exclusion is deliberate (commit 122f417) and the case was not shown on real data. Fix: keep genome-record hits out of the "is it rDNA" vote but include them in the taxonomy LCA.
- **SCI-L5. C-science-20** (low, confirmed; was medium). `R/screen_variants.R:16-19, 72`. Problem: the isomer test compares only the 18S and 28S anchors, but true isomers share identical ITS, while real paralogs differ mainly in ITS. Why: two real rDNA variants that differ only in ITS1 are labelled "isomer, pick either"; information only, and doc 12 calls the threshold an unvalidated prior. Fix: also compare ITS1 plus ITS2, and call units isomers only when both match.
- **SCI-L6. C-science-27 (unverified-low), C-pkg-12 (confirmed)** (low). `R/curate_rrna.R:260-261`, `R/init_db.R:624`, `R/app_opts.R:301`. Problem: an explicit `metazoa_generic` ruleset is treated as "not set" and replaced by lineage auto-selection. Why: a curator cannot choose the wide generic ranges for an unusual fish; affects flags only. Already listed as a known defect in `dev/design/08`. Fix: store `auto` (or NULL) for "not chosen" and honour an explicit `metazoa_generic`.
- **SCI-L7. C-science-21** (low, unverified-low). `R/annotate.R:466`. Problem: `annotate_records` treats a reference record as minus-strand if any single anchor is on minus. Why: a record with one inverted partial copy has all its plus anchors mirrored. Fix: use the same majority rule as `orient_rrna_sense`.
- **SCI-L8. C-science-22** (low, unverified-low). `R/hit_class.R:67`. Problem: `consensus_hit_class` votes per HSP, not per subject. Why: one multi-copy clone with several HSPs can flip a unit to "not nuclear rDNA". Fix: call it on `screen_one_vote_per_subject(h)`.
- **SCI-L9. C-science-23** (low, unverified-low). `R/scaffold_scoring.R:54`. Problem: `scaffold_depth` drops the first and last segment of circular GetOrganelle headers. Why: path depth for circular paths omits segments and warns every run. Fix: take the last whitespace token and strip `(circular)` before splitting.
- **SCI-L10. C-science-24** (low, unverified-low). `R/pipeline_steps.R:40`. Problem: on a rerun with an unchanged sequence, the seeds-database best hit overwrites the GenBank best hit in `assemblies.blast_*`, because BLAST is skipped as already searched. Why: the review pane shows the seeds species instead of the GenBank hit. Fix: call `derive_blast_products` for already-searched units.
- **SCI-L11. C-science-25** (low, unverified-low). `inst/nextflow/bin/rp_score.R:12`. Problem: the max-paths and max-scaffolds gate does not stop the scoring pass. Why: a failed 47-scaffold sample still runs every anchor call and BLAST. Fix: skip `score_one` when the gate failed.
- **SCI-L12. C-science-26** (low, unverified-low). `inst/nextflow/bin/rp_find_rrna.R:26`, `pipeline_r.nf:31-34`. Problem: find_rrna caps kept contigs at `max_scaffolds` (10) and uses the "default" annotate options, while doc 20 says a cap of 20 and the sample's own options. Why: per-sample annotate options are ignored by the contig search, and genomes with more rDNA contigs lose the rest. Fix: update the doc, or read the sample's own options.
- **SCI-L13. C-science-28** (low, unverified-low). `R/screen.R:276`. Problem: unanimous phylum-level discordance with 3 or more hits gives `off_target` at high confidence. Why: a sample from a poorly sampled phylum is flagged as contamination. Fix: require high identity (for example 97 percent) for a discordance-only `off_target`, or report it as uncertain.
- **SCI-L14. C-science-29** (low, unverified-low). `R/triage_rules.R:96`, `R/app_utils.R:253`. Problem: "odd one out" fires at species rank with contamination wording. Why: a congener among three top hits is called a likely contaminant. Fix: evaluate at genus or family, or scale the wording to the rank.
- **SCI-L15. C-science-30** (low, unverified-low). `R/chimera.R:191`. Problem: ribotyper scoring keeps only fatal codes (so the documented LowScoreDifference weight never applies) and uses the default model set. Why: the documented weak signal never counts, and non-eukaryote models can win. Fix: parse non-fatal codes too (or update doc 21) and pass the eukaryote model set.
- **SCI-L16. C-science-32** (low, unverified-low). `R/chimera.R:46`. Problem: the gene-channel reason reports whole-HSP identity as the element identity. Why: "5.8S 82.2%" comes from HSPs spanning ITS1 to ITS2 and suggests a contaminant that is not there. Fix: label it HSP identity, or compute identity over the element.
- **SCI-L17. C-science-33** (low, unverified-low). `R/pipeline_steps.R:378`. Problem: when a supplied contig is reverse-complemented, no flag is stored. Why: users mapping coordinates back to their own contig get mirrored positions. Fix: store an "oriented" flag and show it next to `source_contig`.
- **SCI-L18. C-science-34** (low, unverified-low). `R/curate_rrna.R:323`. Problem: the contaminant flag uses the single top hit while the screen uses an LCA. Why: one mislabelled top record raises "contaminant" while the screen says on target. Fix: base the flag on the screen's LCA or host concordance.

#### Nit

- **SCI-N1. C-science-31** (nit, unverified-low). `R/triage_rules.R:110`. Problem: the chimeric triage reason includes "depth step at ...", while doc 09 says depth never appears in a reason. Why: design and behaviour disagree. Fix: update doc 09 or strip the read-channel text.
- **SCI-N2. C-science-35** (nit, unverified-low). `R/annotate.R:285`. Problem: `run_itsx` leaves ITSx output files in the temp folder. Why: long reference-store runs fill the session temp folder. Fix: `on.exit(unlink(Sys.glob(paste0(pre, '*'))), add = TRUE)`.

### 4.4 Database and export

#### High

**DB-H1. Re-assembly keeps a locked unit's old gene coordinates, and export writes them onto the new sequence**
- Ids: C-data-1 (came in as critical), C-tests-3 (medium). Severity: high. Status: confirmed (three verifiers).
- Where: `R/pipeline_steps.R:398-401` (rows replaced), `:427-435` (`INSERT OR IGNORE` into `annotate`), `:438-447` (orphan cleanup), `inst/nextflow/modules/pipeline_r.nf:103` (only unlocked units are re-annotated), `R/export.R:544-549`.
- Problem: when a sample is assembled again, a unit whose (ID, path, scaffold) key survives keeps its `annotate_lock = 1`, `reviewed = 'yes'`, export group, and old annotation rows. Locked units are never re-annotated, and export picks units by `annotate_lock` alone. Reproduced by two reviewers and a verifier: a locked 3,000 bp unit with 18S at 100-1900 was re-assembled to 1,200 bp, and `export_bundle()` wrote `100 >1900 rRNA` against the 1,200 bp sequence. If the new sequence was reverse-complemented, every stored coordinate is also mirrored.
- Why it matters: the feature table, deflines, and marker slices describe a different sequence. When lengths match, nothing flags it; when they differ, table2asn rejects out-of-range features. It takes deliberate steps (lock annotation, then unlock assembly and re-run), and the Portal route uploads only the FASTA, which is why it is high and not critical. Nothing in the app blocks those steps.
- Fix: in `load_assemblies()`, compare each surviving unit's new sequence with the old one (md5). When it changed, delete its annotations and all other per-unit evidence (BLAST, screen, chimera, coverage), set `annotate_lock = 0` and `reviewed = 'no'`, clear partial and export group, and re-queue it. The simpler guard is to refuse to queue assembly for a sample that has any annotation-locked unit. This one change also covers much of SCI-H1, SCI-M5, PIPE-M1, and DB-L1.
- Status: fixed on `review/2026-09-full-review`. `load_assemblies()` compares each surviving unit's new and old sequence. When they differ, `reset_unit_evidence()` deletes the unit's rows in 11 per-unit tables, clears its lock, Reviewed, Partial, Problematic, and export group, and re-queues it. Test: `test-load-assemblies.R`. Container runs get it after the image is rebuilt. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).

**DB-H2. Truncated genes are exported as "complete sequence" with no partial markers**
- Ids: C-data-2 (high), C-science-7 (medium). Severity: high (one verifier rated it medium because the curator can see the flag in the app). Status: confirmed (three verifiers).
- Where: `R/export_opts.R:134-141` (defline), `R/export.R:144-153` (feature table `<` and `>`), `:628-631`, `:200-222` (`risk_flags`), `R/curate_rrna.R:115-120`.
- Problem: a feature is called partial only when it touches a scaffold end or when the curator pressed the Partial button. Nothing reads `annotations.model_coverage` or the curation "truncated" flag, even though the curation code's own comment says truncated anchors are "flagged rather than exported as a complete gene". In the real export, SRR33849683_p1_s2 (28S at 155-1488 on a 1,694 bp scaffold, 45.81 percent model coverage, flags "length,truncated") got the defline "28S ribosomal RNA gene, complete sequence." and `155 1488 rRNA` with no `<` or `>`. CHIM01's 18S (27 percent coverage) was also called complete.
- Why it matters: GenBank receives records claiming complete genes that are not, and no export warning mentions it.
- Fix: treat a feature as partial when its model coverage is below `min_anchor_coverage` or its flags include "truncated". Use the same rule for the defline and the feature-table markers, and add a "partial not marked" submission risk when they disagree with the curator's setting.
- Status: fixed on `review/2026-09-full-review`. A feature flagged `truncated` now gets `<` and `>` in the `.tbl` and "partial sequence" in the definition line (`is_truncated()` in `R/export.R`, used by `feature_table_text_con()` and `export_header_data()`). The flag already carries the ruleset's `min_anchor_coverage`. The "partial not marked" risk was not added, because the two outputs now agree. Test: `test-export.R`. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).

#### Medium

**DB-M1. `[organism=]` uses the raw mapping label, not the organism confirmed in the Taxonomy tab**
- Ids: C-data-4. Severity: medium (was high). Status: confirmed (two verifiers).
- Where: `R/export_opts.R:18, 21` (default headers use `{Taxon}`), `R/export.R:308` (summary), `R/tax_search.R:61, 80`, `R/export.R:339-341`.
- Problem: the Taxonomy tab writes the confirmed name to `samples.sci_name` and never changes `Taxon`, but the default header and the summary CSV use `Taxon`. The `missing_organism` check only tests for an empty value.
- Why it matters: a typo such as "Upeneous parvus" that the user corrected in the app still goes out as `[organism=Upeneous parvus]`, even though the README promises users can verify or correct the organism. A user can type `{sci_name}` in the template, and NCBI checks names on submission, so it is not silent all the way.
- Fix: the verifier warns that switching blindly to `sci_name` would lose valid labels like "Upeneus sp. A" (fuzzy genus matches set `sci_name` to "Upeneus"). Better: default to `sci_name` when the match was exact, keep `Taxon` otherwise, and flag every unit where the two differ so the curator sees which name is exported.

**DB-M2. Export group names are used as folder names with no checks**
- Ids: C-data-7, C-security-1 (was high), C-pkg-8. Severity: medium. Status: confirmed (four verifiers).
- Where: `R/export.R:571, 610-611` (`label <- group`, then `unlink(recursive = TRUE)`), `R/export.R:652, 663, 748`, `R/app_server.R:778-781` (only trims and checks non-empty), `R/export.R:46` (`set_export_group`).
- Problem: the group text becomes `out/export/<label>`, which is deleted recursively before writing. Reproduced: `.` deleted every other bundle in `out/export`; `..` deleted all of `out/` (every published pipeline output); `../S1` deleted that sample's output; `18S/ITS` or `2026/01` failed halfway with "cannot open file" after `samples/` was written. `../..` points at the project root, which holds `.sqlite`.
- Why it matters: a slash in a group name is a plausible choice and leaves a half-written bundle. The dot names are unlikely, and the overwrite dialog shows the path first, but they destroy data silently once confirmed.
- Fix: validate group names when they are assigned (for example `^[A-Za-z0-9][A-Za-z0-9 ._-]*$`, never `.` or `..`, no `/` or `\`), and in `export_bundle()` check that the normalized folder is inside `out/export` before deleting it.

**DB-M3. The export ITSx check breaks on a space in the group or project path**
- Ids: C-data-9, C-security-4 (was high). Severity: medium. Status: confirmed (three verifiers).
- Where: `R/itsx_check.R:65, 71-72` (mount arguments), `:101` (`system2`), `:302-305`, `:422-423`.
- Problem: the docker `-v wd:wd -w wd` and singularity `--bind wd --pwd wd` arguments are not quoted, and `system2()` joins them into one shell command. Reproduced: an export group named "Batch 1" gave "docker: invalid reference format", while "Batch1" worked. After the failure every unit is marked `not_run(...)`, which replaces the flags stored by the pipeline's ITSx step.
- Why it matters: the final ITSx check silently does not run for any group or project folder with a space, and earlier ITSx warnings (such as `itsx_chimeric`) disappear from Annotation Review until the unit is annotated again. Note: overwriting the warnings on a failed check is by design and tested (see refuted C-data-8), so the fix here is the quoting.
- Fix: `shQuote()` the mount and working-folder arguments (or call `processx::run()` with an argument vector). Rejecting spaces at `new_project()` (PIPE-M13) and in group names (DB-M2) also prevents it.

**DB-M4. Export and pending-count code opens the database with no busy timeout**
- Ids: C-data-13, C-shiny-9. Severity: medium. Status: confirmed (both).
- Where: `R/export.R:25, 36, 76, 539, 742` (`assign_export_group`, `clear_export_group`, `sample_assembler`, `export_bundle`, `export_project`), `R/nf_command.R:163-165` (`nf_pending_count`); compare `open_project_db()` at `R/pipeline_steps.R:170-175` (10 s).
- Problem: these functions call `DBI::dbConnect()` directly, so the busy timeout is 0. Reproduced: while another process held a write lock, Assign Group failed in under 0.1 s with "database is locked", and the Shiny session closed, because `do_assign_group` (`R/app_server.R:768`) has no `tryCatch`. `itsx_store_warnings()` swallows the same error, so export-time ITSx flags are lost quietly.
- Why it matters: curators assign groups and export while the pipeline writes other samples, which is the normal workflow. On NFS projects even the reads can fail.
- Fix: open the database through `open_project_db()` (or pass in `app_con()`), and wrap app-started writes in `tryCatch` with an alert.

**DB-M5. ITS spacers do not touch the rRNA genes, so table2asn warns on every ITS record**
- Ids: C-data-17 (medium), C-science-16 (low). Severity: medium. Status: confirmed (both).
- Where: `R/annotate.R:97-100` (promises a gap-free model), `:167-180` (`clip_spacer` only trims overlaps), `R/export.R:110-158`.
- Problem: when ITSx reports a spacer that stops short of barrnap's 18S or 28S, the gap is kept. In the test project, 21 of 22 units have 1 unannotated base between 18S and ITS1, and 20 have 4 between ITS2 and 28S. table2asn gave `SEQ_FEAT.ITSdoesNotAbutRRNA` on all 6 ITS-bearing records of the exported fish bundle. ITSx's 5.8S matches the vertebrate length better than barrnap's, so the 18S 3' end and 28S 5' end are systematically short.
- Why it matters: every table2asn or BankIt submission carries validator warnings that GenBank staff will ask about. The Portal route sends a single `misc_RNA`, so it is less affected.
- Fix: when ITSx reports the neighbouring spacer, snap the anchor boundary to it (18S end = ITS1 start - 1, 28S start = ITS2 end + 1), or close gaps of 5 bp or less at export and report the change.

#### Low

- **DB-L1. C-tests-5** (low, confirmed; was medium). `R/pipeline_steps.R:441`, `R/screen.R:414-419`. Problem: the re-assembly cleanup removes orphan rows only from `annotate`, `annotations`, and `coverage`; `screen`, `blast_hits`, `blast_screen`, `blast_searched`, `chimera`, and `chimera_windows` keep rows for scaffolds that no longer exist. Why: a later assembly that reuses the key inherits stale evidence; the claimed effect on `samples.reference_support` is real but that column is never read. Fix: add those tables to the orphan cleanup (part of the DB-H1 fix).
- **DB-L2. C-data-11** (low, confirmed; was medium). `R/export.R:611, 726-731`. Problem: `utils::zip(flags = "-rq")` adds to an existing zip, so records removed from a group stay in it. Why: the zip sent to a collaborator or used for table2asn still holds a dropped contaminant. The app writes no zip, so only the R export route is affected. Fix: `unlink(zf)` before zipping.
- **DB-L3. C-data-12** (low, confirmed; was medium). `R/export.R:518, 727-730`, `man/export_bundle.Rd:50`. Problem: with a relative project path (the default "."), the zip path is built before `setwd()`, so zip fails, and its status is not checked. The docs also give the wrong zip location. Why: `export_project()` from the project folder says "Exported" but writes no zip. Fix: normalize the path first, check the zip status, and fix the documented location (`out/export/<label>.zip`).
- **DB-L4. C-data-14, C-data-18, C-pkg-7** (low, confirmed; all were medium or low). `R/init_db.R:115-123`, `R/init_project.R:84-86, 109-111, 118, 187`. Problem: `new_project()` creates `.sqlite` before it finishes checking the mapping and the config, and copies `mapping.csv` over the existing one before the "already exists" check. Why: a mapping without a Taxon column, or `executor = 'slurm'`, leaves an empty database and no `.config`; the corrected re-run only prints "Database already exists. Use force = TRUE" and writes nothing. Re-running on an existing project silently replaces `mapping.csv`. Fix: check required columns and build and validate the config before creating the database; delete the file on failure (or build into a temp file and rename); copy the mapping after the existence check.
- **DB-L5. C-data-15** (low, confirmed; was medium). `R/init_db.R:71`. Problem: `read.csv()` without `colClasses` turns all-digit IDs like `00123` into `123`. Why: catalog or plate-style IDs are renamed in the database, output folders, and GenBank records (only when every ID is numeric). Fix: `read.csv(mapping_fn, colClasses = 'character')`.
- **DB-L6. C-data-19** (low, confirmed). `R/init_db.R:120-140`. Problem: mapping columns that clash with built-in sample columns (case-insensitive in SQLite) crash init (`TaxID`) or are silently overwritten (`lineage`, `rank`). Why: mappings from other pipelines cannot be loaded or lose values. Fix: check reserved names case-insensitively and stop with a clear message.
- **DB-L7. C-data-20** (low, confirmed). `R/pipeline_steps.R:479, 622-636`, `R/coverage.R:61`, `R/blast_ref_align.R:85`, `R/find_rrna.R:83`. Problem: trim, restore, and several delete-then-write functions commit in several steps with no transaction. Why: a failure halfway leaves a trimmed sequence with old feature coordinates. Fix: wrap each in `DBI::dbWithTransaction()` and compute results before deleting.
- **DB-L8. C-data-21** (low, confirmed). `R/export.R:144-153`. Problem: `partial = 'yes'` marks both the first and last features partial, wherever they sit. Why: a complete 18S at an internal position is declared 5' partial, which table2asn reports. Fix: apply partial per feature (from truncation or coverage).
- **DB-L9. C-data-22** (low, confirmed). `inst/nextflow/modules/pipeline_r.nf:177-179`, `R/itsx_check.R:387`. Problem: the pipeline runs ITSx with `|| true` and then creates an empty positions file, so a crash is recorded as `itsx_not_detected`. Why: the curator is told the Portal cannot annotate the record when the check never ran. Fix: treat "empty positions and no `_no_detections` file" as not run.
- **DB-L10. C-data-23** (low, confirmed). `R/export.R:154`. Problem: the feature table ignores `annotations.strand` and always writes start then stop. Why: a reverse unit that the curator locks anyway is exported with plus-strand features. Fix: write stop then start for minus features, or block export when any feature is minus.
- **DB-L11. C-data-24** (low, confirmed). `R/export.R:137-141`. Problem: every N-run is written as `assembly_gap` with `linkage_evidence align genus`. Why: records from user assemblies claim an evidence type that may not apply. Fix: choose the evidence by assembly source, or use "unspecified".
- **DB-L12. C-data-25** (low, confirmed). `R/export.R:566-570`. Problem: the SeqID is the sample ID when a sample has one exportable unit, else `ID_p<path>_s<scaffold>`. Why: locking a second unit later renames the first record, and an ID like `A_p1_s2` can collide with unit p1 s2 of sample `A`. Fix: always add the unit suffix, or check for collisions before writing.
- **DB-L13. C-data-26, C-security-20** (low, confirmed). `R/pipeline_steps.R:122` (`NETWORK_FS`). Problem: the list of network filesystems that must not use WAL misses `ceph`, `9p`/`v9fs` (WSL2 drives), and `wekafs`; unknown types get WAL. (The verifier notes FUSE mounts are already covered as `fuseblk`.) Why: a project on CephFS or a WSL2 Windows drive shares a WAL database across hosts, which SQLite does not support. Fix: default to the rollback journal unless the filesystem is known to be local (ext4, xfs, btrfs, apfs, zfs, tmpfs).
- **DB-L14. C-data-28** (low, confirmed). `R/export_opts.R:164`, `R/export.R:339-341`. Problem: blank tokens render empty source modifiers such as `[specimen_voucher=]`, and validation only dry-runs the first row. Why: some records go out with empty modifiers and no warning. Fix: warn at export when any bracketed modifier renders empty for any record.
- **DB-L15. C-data-31** (low, confirmed). `R/init_db.R:166, 691`. Problem: an all-assembly mapping is still required to have R1 and R2 columns. Why: `ID,Taxon,Assembly` is refused with a message about read files (the vignette does say to leave them blank). Fix: add missing R1 and R2 columns as NA when an Assembly column is present.
- **DB-L16. C-pkg-11** (low, confirmed; was medium). `R/export.R:539, 742`. Problem: `export_bundle()` and `export_project()` pointed at a non-project folder create an empty `.sqlite`, and `export_project()` leaks its connection on the error. Why: the stray file then passes `run_app()`'s project check and blocks `new_project()` there. Fix: stop when `.sqlite` is missing (or open with `flags = RSQLite::SQLITE_RW`), and add `on.exit(dbDisconnect(...))`.

### 4.5 Reference data and taxonomy

#### High

**REF-H1. Changing a sample's taxon never updates `samples.lineage_taxids`**
- Ids: C-refdata-1 (came in as critical). Severity: high. Status: confirmed.
- Where: `R/tax_search.R:61, 80` (`resolve_samples_taxonomy`, `assign_samples_taxid`), `R/resolve_project_taxonomy.R:91`, `R/taxonomy_ranks.R:219-230` (`ensure_sample_lineage_taxids` only fills empty rows).
- Problem: these updates set `taxid`, `sci_name`, and `lineage` but not `lineage_taxids`, and the backfill only fills rows where it is empty. Reproduced: a sample resolved as plant "Morus sp." (3497) and then corrected to the gannet genus (37577) kept the plant's taxid path, and a real gannet hit then scored rank distance 7 ("shares nothing").
- Why it matters: after the user fixes a wrong taxon (via the review modal, or a second Update after the Reference tab has filled the column), the contamination screen scores true host hits as possible contaminants and hits from the old organism as the host, and reference groups are built on the old clade.
- Fix: write `lineage_taxids` in the same update as `taxid` (the lookups already return it), or set it to NULL whenever the taxid changes so the backfill recomputes it. Also clear the sample's screen rows (SCI-H1).

**REF-H2. A partial taxdump download is kept, and the index built from it is trusted forever**
- Ids: C-refdata-2 (came in as critical), C-pkg-4 (medium), C-security-18 (low). Severity: high. Status: confirmed (all three).
- Where: `R/taxonomy.R:54` (`have_taxdump`), `:68` (untar status ignored), `:84-93` (`ensure_taxdump`), `R/tax_index.R:137, 177`, `R/taxonomy.R:292` (falls back to the datasets CLI).
- Problem: the roughly 80 MB taxdump downloads straight to its final name with R's default 60 s timeout and no md5 check, and a later call skips the download whenever the file exists. The untar exit status is ignored, and the index is accepted if it has more than 0 rows. Reproduced on a simulated 1 MB/s link: 63 of 80 MB kept, an index of 99,032 nodes instead of about 2.7 million marked valid, and "Danio rerio" then matched "exact" with a placeholder lineage. With a truncated tarball in another test, taxonomy failed with the misleading "NCBI datasets CLI not found".
- Why it matters: on home Wi-Fi, a VPN, or a busy login node, every later lookup uses the broken index: green "exact" badges with missing lineages, reference ladders with no family or order, and a screen with no usable rank distances. Replacing the tarball by hand does not help, because the bad index still counts as valid. The first failure does show an error, which is why this is high and not critical.
- Fix: download to a temp file with `options(timeout = 3600)` (as `download_or_stop()` in `db_dist.R` already does), check it against NCBI's `taxdump.tar.gz.md5`, then rename. Stop when untar fails. Require a node-count floor (for example over 2 million) before writing `tax_meta`, and store the tarball hash so a new tarball rebuilds the index. When the index cannot be built, say so instead of falling back to the datasets CLI (see ponytail P3).

**REF-H3. Homonym taxon names silently resolve to the lowest taxid, usually the wrong kingdom**
- Ids: C-refdata-3. Severity: high. Status: confirmed.
- Where: `R/tax_index.R:294` (`td_name2taxid` orders by taxid), `R/taxonomy.R:122` (`tk_name2taxid` takes the first), `:245-265` (fuzzy pass), `R/tax_search.R:61`, `R/utils_reactable.R:550`.
- Problem: when a name matches several taxa, the lowest taxid wins. Exact homonyms are marked "ambiguous" but the alternatives are thrown away, and the badge for "ambiguous" is plain grey. The fuzzy pass also takes the lowest taxid with no ambiguity flag. Tested: "Morus sp." becomes the mulberry (3497), the typo "Morus bassanu" becomes the plant genus, and "Oenanthe" (a bird genus and a plant genus) becomes the plant. The index has 1,134 genus-level homonyms. Search results show two identical labels "Morus [genus]" that differ only by taxid.
- Why it matters: a gannet labelled "Morus sp." gets a plant lineage, and reference groups, clade databases, the screen, and scoring all run against plants, with only a grey badge as a hint. Nothing blocks the pipeline.
- Fix: prefer scientific-name matches; when more than one scientific name matches, leave the sample unresolved and ask for a manual pick. In the fuzzy pass try every candidate genus. Store and show the alternatives, and add the kingdom (or NCBI's unique name, for example "Morus <birds>") to search labels.

#### Medium

**REF-M1. Clade databases are counted and downloaded by name, so homonym clades pull in another organism's rDNA**
- Ids: C-refdata-4 (was high), C-security-9 (low), C-ponytail-34 (low). Severity: medium. Status: confirmed (all three).
- Where: `R/reference_counts.R:17-22` (`clade_query` supports a taxid), callers at `R/reference_groups.R:155`, `R/reference_build.R:54, 241, 252`, `R/app_reference.R:654`; `inst/nextflow/bin/rp_ref_db.R`, `download_rrna.sh:49` (has `--taxid`, never passed).
- Problem: every caller passes only the clade name. Checked live: `"Oenanthe"[Organism]` returns 163 records, 149 of them the plant; `"Proboscidea"[Organism]` mixes elephants with a plant genus. `dev/design/19` says the builder queries by taxid.
- Why it matters: a user who picks the wheatear genus gets a seed and label database that is about 91 percent plant rDNA, and GetOrganelle and the scorer then recruit and credit plant sequences. It needs a homonym clade (mostly at genus rank; family-level cases exist but the foreign share is small).
- Fix: pass `clade_taxid` (already stored in the reference table) to every `clade_query()` call, and pass `--taxid` through `rp_ref_db.R` to `download_rrna.sh`. The plumbing already exists.

**REF-M2. The shared reference folder is not safe for several users, and one error fails every running build**
- Ids: C-refdata-5, C-shiny-25 (low, unverified-low). Severity: medium. Status: confirmed (C-shiny-25 unverified-low).
- Where: `R/reference_build.R:609` (`append_ref_manifest` unguarded), `R/app_reference.R:377-392` (`poll_now` error handler), `R/reference_dir.R:70, 148` (`validate_ref_dir`), `R/reference_nf.R:46` (`merge_counts`).
- Problem: `counts.tsv`, `manifest.tsv`, `refstore/`, and the seeds and labels files get the first user's file permissions, and `validate_ref_dir()` only checks that the folder is writable. Reproduced: with read-only files, `merge_counts` and `append_ref_manifest` fail with "Permission denied". When any poll step errors, the handler sets every building group (state 4) to failed, not just the one that errored, and those groups are never wired when their jobs finish. The manifest and counts are rewritten in place, so two projects finishing together can drop rows.
- Why it matters: the vignette advertises lab-shared reference folders. A second user's Propose groups fails, and healthy cluster builds are marked failed with another group's error.
- Fix: catch errors per row inside the poll loop; wrap the manifest and count-cache writes in `tryCatch` (they are only caches); write through a temp file plus rename; check that existing files are writable in `validate_ref_dir()`; and document a group-writable umask for shared folders.

**REF-M3. Clade downloads are never checked for completeness, with no retry and no limit on parallel requests**
- Ids: C-refdata-6. Severity: medium. Status: confirmed.
- Where: `inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh:146-147`, `inst/nextflow/main_reference.nf:37-40`, `inst/nextflow/nextflow.config:49`.
- Problem: `esearch | efetch > raw.fasta` is never compared with the expected record count, and the author's own `build_refdb.sh:753` notes that efetch exits 0 on failure. `build_ref_db` ignores errors with no retry, and neither it nor `count_clade` has `maxForks`. A live sequential test already hit HTTP 429.
- Why it matters: during a multi-clade build, one failed batch gives seed and label files built from a partial download, the group is marked "built", and the manifest records the full count. Assembly then quietly recruits from a smaller seed set.
- Fix: in `download_rrna.sh`, compare the record count in `raw.fasta` with the expected count (allow a small shortfall) and exit 1 otherwise. Give `build_ref_db` a retry with backoff, and set `maxForks` to 2 or 3 on `count_clade` and `build_ref_db`.

**REF-M4. Build on an already-built group re-queues finished samples and flattens their presets**
- Ids: C-refdata-9. Severity: medium. Status: confirmed.
- Where: `R/app_reference.R:307-336, 779-786` (`do_build`), `R/reference_build.R:311-313, 337` (`wire_reference_db`).
- Problem: when the database already exists, Build still calls `wire_reference_db()`, which sets every unlocked sample in the group back to "ready to run", even when the seeds and labels have not changed. The base preset is taken with `LIMIT 1` and no `ORDER BY`. The app never passes `overwrite = TRUE`, so a stale or partial database cannot be rebuilt from the app.
- Why it matters: clicking "Download sequences from GenBank" on a finished group sends every unlocked sample back to the queue, and the next run re-assembles them for nothing. Samples on custom presets silently get another sample's settings.
- Fix: skip the rewire and re-queue when the preset already points at the same files; clone per existing preset instead of `LIMIT 1`; add an explicit Rebuild action that passes `overwrite = TRUE`.

**REF-M5. PDB hits are stored under the wrong key, so the top reference is never aligned**
- Ids: C-refdata-10. Severity: medium. Status: confirmed.
- Where: `R/blast_ref_fetch.R:53-55` (`parse_defline_ids`), `R/blast_ref_align.R:89-99`, `R/blast_genbank.R:364` (`extract_species`).
- Problem: efetch returns PDB records as `>pdb|4V6X|A5 Chain A5, 28S ribosomal RNA`, and the record is stored as `pdb|4V6X|A5` while the candidate table uses `4V6X_A5`. The title still classifies as rDNA, and the species becomes "Chain A5,". `build_refdb.sh` already converts `pdb|X|Y` to `X_Y`, and the local panel holds 298 PDB records.
- Why it matters: for hosts with ribosome structures (human, rabbit, yeast, Drosophila, parasites), a full-length near-identical PDB chain ranks first, and reference-guided alignment is silently skipped. The record is fetched again every run.
- Fix: apply the same `pdb|X|Y` to `X_Y` rule in `parse_defline_ids()` (and `push_to_ref_store()`), or exclude PDB accessions as reference candidates.

**REF-M6. Concurrent taxonomy-index builds share one temp file and can install a broken index**
- Ids: C-refdata-12. Severity: medium. Status: confirmed.
- Where: `R/tax_index.R:92-93, 153`, `inst/nextflow/nextflow.config:49-50` (three BLAST forks), `inst/nextflow/bin/rp_blast.R`.
- Problem: every builder writes `<index>.tmp`, deletes it first, and renames it at the end. Reproduced with 3 concurrent builds: one got "database is locked" and 0 lineages, another "cannot rename", and the final index had no `tax_meta` table.
- Why it matters: after a `TAX_INDEX_VERSION` bump, or when the app's taxonomy step was skipped, the first run's three BLAST tasks rebuild the index together, some units get NA lineages for every hit, and the index on disk is left broken.
- Fix: use `tempfile(tmpdir = taxdump_dir)` plus a lock file, or build and check the index in the app before launching Nextflow and let tasks open it read-only.

#### Low

- **REF-L1. C-refdata-7** (low, confirmed; was medium). `R/reference_groups.R:179`, `R/app_reference.R:426`. Problem: when a record count fails, the clade walk skips that rank (`if (is.na(n)) next`). Why: the proposal lands on a broader clade (for example the order instead of the family) with no hint that a count failed. Fix: stop the walk and flag the unit, or re-run failed counts first, and list failed clades in the proposal.
- **REF-L2. C-refdata-8** (low, confirmed; was medium). `R/reference_build.R:460, 512-517`. Problem: an existing clade database is reused by its file-name slug alone; the manifest's taxid and RefSeq flag are never compared. Why: two projects sharing a reference folder can silently reuse a homonym clade's database or a RefSeq-only build. Rare. Fix: require a matching taxid and RefSeq flag before reuse, or put the taxid in the file name.
- **REF-L3. C-refdata-11** (low, confirmed; was medium). `R/ref_store.R:103-110, 271-292`. Problem: the shared reference store keeps records with no annotations, and never compares the `ribopilot_version` it writes. Why: one run where barrnap or ITSx failed leaves unannotated references for every project in that folder, and annotation-code fixes never reach cached rows. Fix: do not push records with zero annotations, add an annotation version to the filter, and reject sequences whose length differs from `blast_hits.slen`.
- **REF-L4. C-refdata-15, C-security-12, C-data-29 (confirmed), C-pkg-20 (unverified-low)** (low). `R/constants.R:20-29`, `R/init_db.R:51-52, 102`. Problem: the default fish seeds and labels come from an unpinned `main` branch URL (currently 404 while the repo is private), with no checksum, and `new_db()` repeats the same URLs as literals. Why: before release, every default project warns and cannot stage seeds; after release, any push to `main` silently changes the seeds for existing projects (the file already changed once, commit 38c0aca). The duplicate literal also causes the `new_db.Rd` line-width NOTE. Fix: pin to a release tag or commit with a sha256 (or ship the defaults in the image), and default `new_db()`'s arguments to NULL so the constants are the single source.
- **REF-L5. C-nextflow-21 (confirmed), C-refdata-18 (unverified-low)** (low). `R/reference_build.R:579, 611-614`, `inst/nextflow/main_reference.nf`. Problem: a cluster-submitted build that fails leaves no FASTA, and the poller treats "no fresh FASTA" as still running, with no timeout; the poller also accepts a FASTA that is still being copied. Why: the group shows "building" forever, or a half-copied seeds file is wired. Fix: have the build write a completion or failure marker, wire only after it exists, and add an age timeout.
- **REF-L6. C-security-17** (low, confirmed). `R/reference_build.R:385, 597`. Problem: build polling trusts a bare pid on the current host. Why: after a reboot or on another login node, an unrelated process keeps a group "building", or a live build is marked failed. Fix: store the host name and process start time with the pid.
- **REF-L7. C-refdata-19 (unverified-low), C-security-16 (confirmed), C-pkg-27 (unverified-low)** (low). `R/reference_dir.R:117-133`, `inst/nextflow/main_reference.nf:43-44`. Problem: the manifest is written without quotes and read back with quote handling, so `"Mullidae"[Organism]` loses its quotes on the next append; the manifest has no lock, and builds publish straight onto the shared seeds and labels names. Why: the recorded query for a multi-word clade changes meaning, and two projects building the same clade overwrite each other. Fix: read with `quote = ""` (or write with quotes), lock the manifest, and publish to a temp name then rename.
- **REF-L8. C-refdata-16** (low, unverified-low). `inst/nextflow/modules/blast_ref.nf:44`. Problem: in the container, efetch fails once on every call (curl 56) and succeeds on retry; the HTTP/1.1 shim is on PATH only in the download scripts, and pacing is per task while three tasks run. Why: doubled requests push past NCBI's rate limit, and each failure writes the full URL (with the API key) to `.command.err`. Fix: put the shim on the image PATH and divide the pacing by `maxForks`.
- **REF-L9. C-refdata-17** (low, unverified-low). `R/db_dist.R:142`. Problem: `GITHUB_PAT` headers are ignored when a site sets the download method to curl or wget, and the error then blames a slow connection. Why: test-read and panel downloads fail on such sites while the repo is private. Fix: pass `method = 'libcurl'` when headers are set, and map a 404 to "asset not found or private; set GITHUB_PAT".
- **REF-L10. C-refdata-20** (low, unverified-low). `R/taxonomy.R:34`. Problem: hybrid names keep only the first two words, so "Citrus maxima x Citrus reticulata" resolves as an exact match to one parent. Why: hybrid samples are recorded as a pure species with no flag. Fix: detect " x " and return unresolved or fuzzy.
- **REF-L11. C-refdata-21** (low, unverified-low). `R/blast_local.R:88`. Problem: `RIBOPILOT_NO_LOCAL_BLAST=0` turns local BLAST off, because any non-empty value counts. Why: a user who sets 0 or "false" meaning "keep local" silently gets slow remote searches. Fix: treat only 1, true, or yes as opt-out.
- **REF-L12. C-refdata-22** (low, unverified-low). `R/reference_build.R:263`. Problem: "Leave these samples on the built-in fish database" does not revert a group that was already wired to a clade database. Why: the samples keep assembling against the clade database. Fix: revert the samples to the default preset, or reword the button for wired groups.
- **REF-L13. C-refdata-23** (low, unverified-low). `R/blast_genbank.R:364`. Problem: `extract_species()` returns prefixes such as "PREDICTED: Danio" or "Chain A5," as the species. Why: the Species column, reference species, and chimera reasons show the wrong name. Fix: strip those prefixes, or take the species from the staxid lineage.
- **REF-L14. C-refdata-24** (low, unverified-low). `dev/refdb/build_refdb.sh:635`. Problem: the awk accession scan uses `NR==FNR`, which stays true for stdin when the wanted list is empty. Why: the scan loads a multi-GB dump into memory for nothing. Fix: skip the scan when the list is empty, or use `FILENAME == ARGV[1]`.
- **REF-L15. C-refdata-25** (low, unverified-low). `R/tax_index.R:23`. Problem: each project keeps its own 400 MB taxonomy index plus the 80 MB tarball. Why: 20 projects spend about 10 GB of HPC quota on identical copies. Fix: share the index in `ref_dir`, keyed by the tarball hash (the author's rule is that shared caches go under `ref_dir`).
- **REF-L16. C-refdata-26** (low, unverified-low). `R/reference_counts.R:41`. Problem: the record-count cache never expires and ignores the RefSeq flag. Why: old counts drive the walk after GenBank grows, and RefSeq-only databases list the full GenBank count. Fix: store the count date and RefSeq flag, and recount after a set age.
- **REF-L17. C-refdata-27** (low, unverified-low). `R/reference_build.R:4, 368-381`. Problem: the file header says in-app builds are detached and survive a closed app, but the app launches them with `cleanup = TRUE`. Why: a user who starts a build and closes the tab returns to a failed build. Fix: launch detached with output to a log file, or correct the header and warn in the UI.
- **REF-L18. C-refdata-28** (low, unverified-low). `inst/ref_dbs/getOrganelle/GenBankDownload/README.md:25`, `R/utils.R:97-100`. Problem: the README ("From R (recommended)" with `conda_env`), the `build_reference_db` help, and the `check_reference_db` hint push conda and scripted builds. Why: new users follow the unsupported conda path instead of the app's Reference tab. Fix: point to the Reference tab, keep a short headless Docker note, and drop the conda instructions.

### 4.6 Security

No reviewer found SQL injection: every query value is a bound parameter. The app has no login, which is normal for a local Shiny app but matters on shared login nodes.

#### Medium

**SEC-M1. FASTA header templates are run as R code, including templates stored in the project database**
- Ids: C-security-3 (was high), C-data-6 (low), C-pkg-13 (low). Severity: medium (the C-data-6 and C-pkg-13 verifiers rated it low, because an attacker needs write access to the project, which already lets them plant a Nextflow config). Status: confirmed (four verifiers).
- Where: `R/export_opts.R:214, 250` (`stringr::str_glue_data`), `R/export_opts.R:47-56` (templates loaded from the database), `R/app_server.R:1149-1160` (validated on open and on each debounced keystroke), `R/export.R:594-602, 646`, `R/app_run.R:192-200`.
- Problem: glue evaluates any R expression inside braces. Reproduced: `validate_fasta_header('{seqid} {nchar(system("id -un", intern=TRUE))}', d)` ran the command and returned "ok". A saved template runs as soon as the Export modal opens. Headless mode binds `127.0.0.1` with no authentication.
- Why it matters: another user on a shared login node who can reach the port, or anyone who can write a shared project `.sqlite`, can run commands as the owner. The editable job script and Submit button already give similar access, so glue adds one more route rather than a new class of risk, but it is also the easiest to close.
- Fix: replace glue with plain token substitution: match `\{([A-Za-z0-9_.]+)\}` and look the name up in the data columns, rejecting anything else (this also removes the last real use of stringr). For headless use, document that the app must not be exposed to other users, and consider a random per-session token in the URL.

#### Low

- **SEC-L1. C-security-2 (confirmed; was high), C-refdata-13 (confirmed; was medium), C-nextflow-26 (unverified-low), C-shiny-24 (unverified-low)** (low). `R/nf_command.R:31, 58`, `R/reference_nf.R:175`, `inst/nextflow/nextflow.config:37`, `R/app_server.R:851-854, 904`, `R/app_reference.R:837-841, 858`, `R/app_modals.R:483, 523`, `dev/refdb/lib.sh:47`. Problem: the NCBI API key is protected in `.config` (mode 0600) but then passed as `--ncbi_api_key`, so it lands in `ps` output, `.nextflow.log`, `.nextflow/history`, every task's `.command.run` (mode 664), job scripts saved as 0755, the command preview, and the Reference modal's plain text box. The maintainer's log scrubber misses the `-api_key KEY` form that EDirect prints on errors. Why: on a shared node or a group-readable project other users can read the key, and docs screenshots of the run modal will show it. Severity is low because an E-utilities key only raises rate limits (the worst case is a used-up quota or a revoked key); `dev/design/17` already lists this as deferred. Fix: pass the key only through the launched process's environment (and `docker.envWhitelist`), not argv or `params`; mask it in the preview and do not prefill the text box; write job scripts as 0700; change `RP_SCRUB_RE` to match `-?api[_-]?key[=: ]+`; redact keys in `log_tail` and `failed_task_errors`. Check any docs screenshots before publishing.
- **SEC-L2. C-security-7** (low, confirmed; was medium). `inst/nextflow/nextflow.config:37`, `R/app_modals.R:523`, `R/init_project.R:181`. Problem: the API key is written unescaped into every task's shell as `export NCBI_API_KEY="..."`. Reproduced: a value with a quote and a command ran that command in every task. Why: this is self-injection with no privilege gain unless the app is exposed; the realistic harm is that a pasted key with a stray quote or `$` breaks every task with an unclear bash error. Fix: validate the key where it enters (NCBI keys are 36 lowercase hex characters) and stop routing it through `params` (SEC-L1).
- **SEC-L3. C-security-11 (confirmed; was medium), C-ponytail-36 (confirmed), C-pkg-15 (unverified-low)** (low). `.Rbuildignore`, `docker/Dockerfile:163-166`, `docker/deploy-dockerhub.sh:71`, `docker/deploy-local.sh`. Problem: `.Rbuildignore` misses `.claude`, `.superpowers` (39 files of agent reports and diffs, about 1.9 MB), `.nextflow.log`, `update-version.log`, and `RiboPilot_test_userAsmb.R`; `devtools::build` does not read `.gitignore`, so all of them go into the tarball, and the Dockerfile copies that tarball into its own layer of the public image (the later `rm` does not remove a layer). The repo root also holds personal driver scripts, a committed build log, BUG and HANDOFF notes, a broken `dev/run_test_assembly.sh` (passes `--samplesheet`, which `main.nf` no longer reads), `LazyData: true` with no `data/`, and a stale GenBankDownload README. Why: the public image ships local paths, review notes, and a hook config today; if the author ever runs Nextflow by hand from the repo root with the key set, `.nextflow.log` would carry the key into the image. The tarball is 1.43 MB instead of 0.86 MB, and checks raise hidden-file and top-level-file NOTEs. Fix: add `^\.claude$`, `^\.superpowers$`, `^\.nextflow\.log.*$`, `^update-version\.log$`, and `^RiboPilot_test.*\.R$` to `.Rbuildignore` (or build from `git archive`); install the tarball in the image with `RUN --mount=type=bind` or a multi-stage build; delete or move the root clutter to `dev/`; delete `dev/run_test_assembly.sh`; drop `LazyData`; update or delete the GenBankDownload README.
- **SEC-L4. C-security-14** (low, confirmed). `R/app_server.R:1496`, `R/utils_reactable.R:235-262, 326-341`, `R/export.R:474-499`, `inst/export_report.Rmd:49-55`. Problem: GenBank titles, species, clade names, and `fail_reason` are inserted as raw HTML, and the Taxon is pasted into R Markdown source that is knitted. Why: a title with quotes or angle brackets breaks table cells; a crafted top hit could run script in the app (which can set inputs such as the header template), and a Taxon containing inline R code runs while knitting the report. Fix: escape with `htmltools::htmlEscape()` in R and the existing replace chain in the JS renderers, and pass the Taxon to the report as data, not Markdown source.
- **SEC-L5. C-security-21** (low, confirmed). `R/db_dist.R:173`, `R/init_test_project.R:87-91`. Problem: a test-read download that fails its checksum is left in place, and the next run accepts any file that passes a gzip test. Why: a wrong but valid gzip is used without re-checking (unlikely). Fix: download to a temp name and rename only after the checksum matches, or delete on mismatch.

### 4.7 Packaging

#### High

**PKG-H1. A clean install fails: DESCRIPTION has no `biocViews` field**
- Ids: C-pkg-1. Severity: high. Status: confirmed.
- Where: `DESCRIPTION:16` (Imports), `README.Rmd:91`, `vignettes/RiboPilot.Rmd:102`.
- Problem: Imports lists three Bioconductor packages (Biostrings, DECIPHER, pwalign), but `remotes` only looks in Bioconductor when a `biocViews` field is present. Reproduced: `remotes::dev_package_deps()` with CRAN repos shows all three as not available; after adding an empty `biocViews:` line it found them. CI passes only because it installs with pak, which searches Bioconductor on its own.
- Why it matters: a new user on a fresh R install who follows the README runs `remotes::install_github('dmacguigan/RiboPilot')` and gets "dependencies Biostrings, DECIPHER, pwalign are not available".
- Fix: add `biocViews:` (it can be empty) to DESCRIPTION. Optionally also show `BiocManager::install()` in the README.

**PKG-H2. `new_project()` and `new_test_project()` crash outside RStudio with the default `Rproj = TRUE`**
- Ids: C-data-3, C-pkg-2. Severity: high. Status: confirmed (three verifiers).
- Where: `R/init_project.R:98-105`, `R/init_test_project.R:158, 166, 199`, `README.md:166`, `vignettes/RiboPilot.Rmd:128, 190`.
- Problem: the guard is `Rproj && !isFALSE(Sys.getenv("RSTUDIO", FALSE))`. `Sys.getenv()` returns the string "FALSE", and `isFALSE("FALSE")` is FALSE, so the test is always TRUE. `rstudioapi` is in Imports, so `rstudioapi::initializeProject()` is then called and stops with "RStudio not running". Reproduced from Rscript: the folder is left with only `mapping.csv` (or the staged reads), no `.sqlite`, and no `.config`. The tests and the author's own driver scripts always pass `Rproj = FALSE`, which is why it stayed hidden, and the `new_test_project` help wrongly says Rproj has "no effect outside RStudio".
- Why it matters: anyone following the README quick start from a terminal, Positron or VS Code, an HPC login node, or Rscript fails at step one, and nothing tells them about `Rproj = FALSE`.
- Fix: `if (Rproj && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())`, and register the open-project step only after the project is fully built. Consider defaulting `Rproj` to FALSE.

#### Medium

**PKG-M1. Roxygen markdown is off, so every help page shows raw markdown**
- Ids: C-pkg-6. Severity: medium. Status: confirmed.
- Where: `DESCRIPTION:57` (no `Roxygen: list(markdown = TRUE)`), for example `man/run_app.Rd:28`.
- Problem: 44 of the roughly 54 Rd files contain literal backticks, 25 contain literal `[fn()]` brackets, and there are no `\link` entries at all. Turning markdown on (in a copy) produced 16 links from exported pages to internal, undocumented topics (`detect_headless`, `rrna_params`, `ref_config`, `ribopilot_config_dir`, and others).
- Why it matters: `?new_project` and the unpublished pkgdown reference show raw markdown with no links. If markdown is simply turned on, `R CMD check` gains "Missing link" warnings.
- Fix: add `Roxygen: list(markdown = TRUE)`, re-document, and replace links to internal topics in exported docs with plain code text.

#### Low

- **PKG-L1. C-data-30 (confirmed), C-pkg-23 (unverified-low)** (low). `R/zzz.R:31`. Problem: `.onAttach` runs `nextflow -version` on every attach and calls `stop()` when Nextflow is too old. Why: a user with an old Nextflow module cannot even `library(RiboPilot)` to open the app, review, or export, none of which need Nextflow, and every conda-route task pays about 0.5 s of JVM start. `run_app`, `run_assembly`, and `hydra_setup` already check the version. Fix: keep at most a startup message, and drop the `stop()` and `in_r_build_context()`.
- **PKG-L2. C-pkg-17** (low, unverified-low). `NAMESPACE:1`. Problem: NAMESPACE is marked "Hand-maintained", so roxygen skips it, although it is byte-identical to roxygen's output. Why: a new `@export` tag is silently ignored, and a bin script calling that function after `library(RiboPilot)` fails in the container (the same bug class as PIPE-M17). Fix: restore the roxygen header and let roxygen own NAMESPACE.
- **PKG-L3. C-pkg-18** (low, unverified-low). `tests/testthat/test-app-utils.R:83, 97, 111`. Problem: tests call `rlang::current_env()`, but rlang is not in Suggests. Why: `--as-cran` WARNING, and the tests fail on a machine that installs only declared packages. Fix: add rlang to Suggests, or use base `environment()`.
- **PKG-L4. C-pkg-19 (unverified-low), C-ponytail-41 (nit, confirmed)** (low). `R/app_modals.R:722-723, 904-906`, `R/app_run.R:64-68`, `R/generate_config.R:49-55`, `R/screen.R:91-101`, `R/utils.R:66`. Problem: several roxygen blocks sit on the wrong object or are orphaned, roxygen 8.1 reports 4 "@noRd must not be followed by any text" errors, and an inline `@noRd` fails so `man/check_container_image.Rd` documents an internal function. Some comments are stale (coverage per path vs all at once, the curl shim naming `download_mitogenomes.sh`, `blast_genbank` claiming it sets `assemble.blast_*`). Why: the pkgdown reference will list an internal helper, and docs attach to the wrong objects. Fix: move each block onto its own object, put `@noRd` on its own line, and fix or delete the stale comments.
- **PKG-L5. C-pkg-21** (low, unverified-low). `man/*.Rd`. Problem: 7 exported pages have no `\value`, only 1 of 53 has examples, and there is no package-level help page (`?RiboPilot` opens the app launcher). Why: CRAN or Bioconductor would reject it, and users get no return-value descriptions. Fix: add `@return`, a `_PACKAGE` doc, and `\dontrun` examples for the user-facing entry points.
- **PKG-L6. C-pkg-22** (low, unverified-low). `R/app_run.R:191`. Problem: `run_app()` sets `options(RiboPilot.headless = ...)` and never restores it. Why: later code in the same session sees a stale value, and CRAN policy disallows it. Fix: pass headless through the app options instead of a global option.
- **PKG-L7. C-pkg-24 (unverified-low), C-pkg-38 (unverified-low, nit)** (low). `R/export_report.R:118-119`, `R/build_reference_db.R:95`, `R/assemble_fallback.R:186`. Problem: temp folders are named with `Sys.time()` (and `sample.int()`). Why: `sample.int()` changes the user's random-number state after `set.seed()`, and two runs in the same second share and delete one folder. Fix: use `tempfile('ribopilot_report_')` and similar.
- **PKG-L8. C-pkg-25** (low, unverified-low). `README.md:29, 51, 114, 134, 172-175`. Problem: `README.md` has curly quotes from pandoc, and a paragraph about filtered reads and `full_size = TRUE` exists only in the `.md`, not in `README.Rmd`. Why: the next render deletes that paragraph, and the ASCII-only rule is broken. Fix: move the paragraph into `README.Rmd` and turn off pandoc's smart quotes (`md_extensions: -smart`).
- **PKG-L9. C-pkg-28** (low, unverified-low). `docker/Dockerfile`, `DESCRIPTION:16`. Problem: the Dockerfile's R package list has drifted from Imports (it adds r-config, r-gggenes, r-purrr, and r-svglite, and says Biostrings is the only Bioconductor dependency), and the package is installed with `dependencies = FALSE`. Why: adding an Import without editing the Dockerfile gives an image where every pipeline task fails at `library(RiboPilot)`. Fix: generate the install list from DESCRIPTION at build time, or install with dependencies.
- **PKG-L10. C-nextflow-32** (low, unverified-low). `docker/update-version.sh:29`. Problem: the script commits the version bump before the image build and push, which can fail (for example `docker login` with no TTY). Why: DESCRIPTION is bumped and committed with no image for that tag, and new projects then point at a tag that was never published. Fix: build and smoke-test first, then commit, then push.
- **PKG-L11. C-tests-19 (unverified-low, nit), C-pkg-16 (unverified-low)** (low). `R/app_modals.R:853`, `R/zzz.R:57-60`. Problem: `ggplot2::aes(y = y)` in `plot_chimera` gives a "no visible binding for global variable y" NOTE, and five `globalVariables` entries (start, label, taxon, pident, src) are no longer needed. Why: a standing NOTE hides new real ones, and stale entries could hide a real typo. Fix: use `.data$y` (or add `y`), and drop the five stale names.

#### Nit

- **PKG-N1. C-pkg-26** (nit, unverified-low). `R/constants.R:194`. Problem: a comment says "constants.R collates last", but with no Collate field files load alphabetically and it loads 13th. Why: a developer may move a dependent constant into an earlier file and break installation. Fix: fix the comment, or add a Collate field if order matters.
- **PKG-N2. C-pkg-29** (nit, unverified-low). `DESCRIPTION:55`, `.github/workflows/R-CMD-check.yaml`. Problem: `LazyData: true` with no `data/`, a Title not in title case, a Description starting "An R package for", and CI building with `--no-build-vignettes` although every vignette chunk is `eval = FALSE`. Why: check NOTEs, and vignette breakage goes unnoticed in CI. Fix: drop LazyData, use title case, start the Description with a verb phrase, and let CI build the vignettes.

---

## 5. Ponytail: what to delete or simplify

The code-ponytail reviewer estimates the total at **about 2,600 lines of R, JS, shell, and Nextflow, plus about 3,000 test lines that only exercise `dev/` scripts, about 1,200 generated Rd lines, and 5 dependencies** (dplyr, dbplyr, golem, stringr, and crayon; readr and ps are optional extras). Other reviewers' dead-code counts (about 250 lines from code-data, 330 from code-refdata, 350 from code-tests, 400 from code-pkg) overlap with this list. Line counts below are the reviewers' rough estimates. The "why" for every item is the same unless noted: less code to read, test, and keep in step, and fewer places for two copies to drift apart.

Biggest wins first:

1. **P1. One Nextflow run surface instead of two.** Ids: C-shiny-34 (unverified-low), C-security-25, C-ponytail-5 (low, confirmed). Where: `R/app_modals.R:426, 482`, `R/app_server.R:896-1052`, `R/app_reference.R:743-1000`, `R/reference_nf.R:208`. Problem: the run modal, output streaming, kill, session end, close confirm, script save and submit, and the Nextflow environment setup are written twice (the environment up to four times). Why beyond upkeep: they have already drifted (the Resume bug in APP-L3, the PATH check in APP-M3), and APP-M1, APP-M2, and APP-M9 would each need fixing twice. Fix: one `nf_run_modal(prefix, ...)`, one server helper, and one `nf_env()`. About 150 lines (one verifier thinks about 80 is more realistic).
2. **P2. Delete the dead CSV taxonomy workflow.** Ids: C-refdata-30, C-tests-21, C-pkg-32 (unverified-low), C-ponytail-6 (low, confirmed). Where: `R/resolve_project_taxonomy.R:1-106`. Problem: `resolve_project_taxonomy()` and `apply_project_taxonomy()` have no callers, and they write taxonomy without `sci_name` or `lineage_taxids` (a trap if revived). Fix: delete both and move the live `write_path_scores()` next to `score_one()`. About 105 lines.
3. **P3. Delete the network `datasets` taxonomy fallback.** Ids: C-refdata-32 (unverified-low), C-ponytail-7 (low, confirmed). Where: `R/taxonomy.R:279-384`, `:18`. Problem: it only runs when the local index fails, then silently switches to slower per-taxon network calls with different homonym rules (and gives the misleading error in REF-H2). Fix: delete it, stop with a clear message when the index cannot be built, and drop ncbi-datasets-cli from the image and env. About 110 lines.
4. **P4. Delete `install_ribopilot_db()` and helpers.** Ids: C-refdata-31, C-pkg-33 (unverified-low), C-ponytail-8 (low, confirmed). Where: `R/db_dist.R:188-316`, `tests/testthat/test-db-dist.R`. Problem: no production caller; the panel is baked into the image. Fix: delete it with `sweep_stale_installs()` and `panel_missing_parts()` and their tests; keep `fetch_verified()`, `download_or_stop()`, and `rp_sha256()`. About 130 lines plus about 400 test lines.
5. **P5. Delete the dead `rp_ref_fetch.R` chain.** Ids: C-nextflow-33, C-refdata-29 (unverified-low), C-security-24, C-ponytail-9 (low, confirmed). Where: `inst/nextflow/bin/rp_ref_fetch.R`, `R/blast_ref_fetch.R:1-116`. Problem: no `.nf` file runs the script (a test asserts this), yet `write_blast_ref()` is still exported with a man page. Fix: delete the script, `write_blast_ref`, `parse_feature_table`, `ref_product_to_element`, `read_fasta_seq`, and their tests; keep `parse_defline_ids` and `is_fasta_file`. About 100 lines.
6. **P6. Delete dead JS handlers.** Ids: C-shiny-35 (unverified-low), C-ponytail-10 (low, confirmed). Where: `inst/app/www/custom.js:40-123`. Problem: `hScroll`, `toggleDetails`, `copy_to_clipboard`, `rightScroll`, `msaScrollToCol`, `scrollProgress`, the popover init, and the `.mp-sample-table` selector have no sender or target (MitoPilot leftovers). Fix: delete them and set the table selector to `#tbl`. About 85 lines.
7. **P7. Delete `score_project()`.** Ids: C-science-36, C-tests-22, C-pkg-31 (unverified-low), C-ponytail-11, C-data-33 (low, confirmed). Where: `R/score_project.R:37-102`. Problem: a pre-Nextflow driver with no callers that copies `score_one()` (including the path-numbering bug in SCI-H3). C-data-33 also lists `resolve_project_taxonomy`, `apply_project_taxonomy` (P2), and `LOCK_UNIT_SQL` (P8). Fix: delete it; move `local_ref()` and `make_blast_db()` next to `score_one()`. About 75 lines.
8. **P8. Delete small dead functions, arguments, and constants.** Ids: C-shiny-36, C-shiny-37, C-science-39, C-science-40, C-science-42, C-science-43, C-science-44, C-tests-24, C-tests-25, C-pkg-35 (all unverified-low), C-ponytail-29, C-security-29 (nit, confirmed). Where: `rt_bool_bttn` (`R/utils_reactable.R:606-662`), `row_in_selection` (`R/app_guards.R:38`), `annotate_seq` (`R/annotate.R:374`), `curate_actinopterygii_rrna` and `curate_rrna_core(chimeric=)` (`R/curate_rrna.R:44, 166`), `screen_anchor_string` (`R/screen_variants.R:253`), `shared_taxon` (`R/taxonomy_ranks.R:137`), `LOCK_UNIT_SQL` (`R/constants.R:171`), `CHIMERA_LABELS`, `HIT_CLASSES`, `SCREEN_VERDICTS`, `SCREEN_CLASSES`, `SCREEN_NO_REF_CODES`, `SCREEN_WEIGHTS`, `chimera_channel_ribotyper(spans=)` and the duplicate `track()` depth parser (`R/chimera.R:200, 363`), the scalar branch of `anchor_pid_floor` (`R/params_rrna.R:144`), `rp_sample_nav(pfx)`, and the unreachable `processx_launcher` default (`R/reference_build.R:368`). Problem: each has no production caller; several are kept alive only by their own tests. Fix: delete them with their tests and Rd files (or use the vocabulary constants in `stopifnot()` checks on write). About 140 lines.
9. **P9. One coverage writer and one depth parser.** Ids: C-science-37, C-tests-23, C-pkg-34 (unverified-low), C-data-32, C-ponytail-24 (confirmed). Where: `R/coverage.R:17-97`, `R/db_schema.R:200-234`, `R/pipeline_steps.R:544-549`, `R/chimera.R:363`, `R/app_modals.R:754`. Problem: `write_coverage()` and `backfill_coverage_median()` are exported with no callers (the latter's help text is false), and stats are computed three ways and depths parsed four ways. Fix: delete both functions, their Rd files, and the comment at `R/app_utils.R:207`; use one `depth_stats()` and `parse_depths()` everywhere. About 95 lines.
10. **P10. Unexport the pipeline internals.** Ids: C-ponytail-15 (low, confirmed), C-pkg-37 (unverified-low). Where: `NAMESPACE` (53 exports), `inst/nextflow/bin/*.R`. Problem: about 40 exports exist only for bin scripts or internal use, and `_pkgdown.yml` has no `reference:` section, so the site will list about 50 flat topics users should never call, each a public API to keep stable. Fix: keep about 10 user-facing exports (new_project, new_test_project, generate_config, hydra_setup, RiboPilot or run_app, ribopilot_ref_dir, export_bundle, build_reference_db); call the rest with `RiboPilot:::` from bin scripts (several already do). Do this after PIPE-M17's scan test exists. About 1,200 generated Rd lines.
11. **P11. One bin-script dispatcher.** Ids: C-nextflow-40 (unverified-low), C-security-30, C-ponytail-16 (nit, confirmed). Where: `inst/nextflow/bin/rp_*.R` (18 scripts). Problem: each repeats the same 4-line bootstrap and connect block. Fix: one `rp.R <step> <db> args...` dispatcher mapping step names to functions. About 100 lines.
12. **P12. Drop five dependencies.** Ids: C-data-35, C-ponytail-12, C-ponytail-13, C-ponytail-14 (low, confirmed), C-pkg-36 (unverified-low). Where: `R/init_db.R` (9 `rows_upsert` calls into tables just created, so empty), `R/app_opts.R:579`, `R/zzz.R:10, 36-41`, `R/app_utils.R:6-10, 104`, `R/app_run.R:202`, `R/generate_config.R:72`, `R/export_opts.R:214, 250`, `R/nextflow_version.R:68`. Problem: dplyr and dbplyr are used only for those upserts, golem for 5 calls base Shiny covers, stringr for 4 or 5 calls, and crayon only for the startup banner. `--as-cran` notes "Imports includes 25 non-default packages". Fix: `DBI::dbAppendTable()` and one `INSERT ... ON CONFLICT`; `shiny::addResourcePath`, a favicon tag, and `shinyOptions(project=)`; `gsub(fixed = TRUE)`, `regmatches`, and plain token substitution (SEC-M1); plain text or cli for the banner. Remove them from the Dockerfile too. About 55 lines, 5 fewer dependencies.
13. **P13. Drop dead schema.** Ids: C-ponytail-17 (low, confirmed). Where: `R/init_db.R` (several), `R/screen.R:392, 412-430`, `R/init_project.R:447`, `R/generate_config.R:58`. Problem: the `assembly_blast` table, `assemble.hide_switch`, six `assemble.blast_*` columns, `assemblies.depth/gc/errors/edit_positions/length_raw`, `samples.export_group`, the write-only `samples.reference_support` and `update_reference_support()`, `assemble_opts.target` and `new_db(target=)`, the `GENETIC_CODE` token, and `new_db(...)` are never used. Why beyond upkeep: users see these always-blank columns among the header template's "RiboPilot columns". Fix: drop them; pre-release, so no migration is needed. About 75 lines.
14. **P14. Remove back-compat shims.** Ids: C-ponytail-18 (low, confirmed), C-data-39 (nit, confirmed). Where: `R/pipeline_steps.R:312-321, 492, 557, 623`, `R/blast_genbank.R:202-212, 229-275, 293`, `R/curate_rrna.R:295-298`, `R/screen.R:346`, `R/taxonomy_ranks.R:216`, `R/export_opts.R:37-57`, `R/itsx_check.R:360-368`, `R/app_utils.R:168-169`, `R/utils_reactable.R:481-484`, `R/app_server.R:1127`, `inst/nextflow/modules/blast.nf:42-45`, `R/utils.R:15` (`%nin%`, one use). Problem: code for "older projects" that the no-migrations rule says are not supported. Why beyond upkeep: the `tryCatch(error = NULL)` wrappers around tables `new_db` always creates also hide real errors, such as a locked database, as silent defaults. Fix: delete them (the verifier notes `new_db` itself uses `ensure_export_opts` to create its table, so keep that DDL). About 110 lines.
15. **P15. One busy-alert helper.** Ids: C-ponytail-19 (low, confirmed). Where: `R/app_server.R:922, 949, 957-960`, `R/app_reference.R:403, 618, 887`, `R/tax_search.R:141, 148-153`. Problem: the busy alert, `later()`, `tryCatch`, result alert pattern is copied six times, and two copies close the database only on success. Why beyond upkeep: if taxonomy resolve or assign errors, a connection stays open for the life of the session. Fix: one `with_busy_alert()` that owns the connection with `on.exit`. About 50 lines.
16. **P16. One deploy script.** Ids: C-nextflow-35 (unverified-low), C-ponytail-20 (nit, confirmed). Where: `docker/deploy-local.sh:32-118`, `docker/deploy-dockerhub.sh:64-160`. Problem: about 58 to 65 lines (tarball build, asset sweep, fetch, checksum) are copied and already differ. Fix: source one `docker/lib.sh`, or have the Docker Hub script call the local one and then tag and push. About 65 lines.
17. **P17. Remove unused image tools and fix the smoke test.** Ids: C-nextflow-34 (unverified-low), C-ponytail-21 (low, confirmed). Where: `docker/Dockerfile:15-51, 78-106`, `environment.yml:18-31`, `docker/deploy-dockerhub.sh`. Problem: fastqc (plus a JRE), minimap2, seqtk, seqkit, GNU parallel, build-essential, r-gggenes, r-purrr, and r-svglite are never called, and the smoke test checks seqkit but not pybarrnap, ribotyper, or jq. Why beyond upkeep: a 7.66 GB image slows pulls (PIPE-L2), and an image missing pybarrnap would pass the smoke test and then annotate nothing. Fix: remove them and smoke-test the tools the pipeline calls. About 15 lines and a smaller image.
18. **P18. One pill renderer.** Ids: C-ponytail-22 (nit, confirmed). Where: `R/utils_reactable.R:121, 172, 326, 394, 411, 545`, `R/app_modals.R:9`, `R/app_server.R:1906`. Problem: seven badge renderers hand-roll pill HTML beside `rt_pill`, with different escaping rules; `status_badge()` has unused arguments. Fix: express the badges through `rt_pill(map, labels)` and reuse `status_badge()` for the modal pills. About 60 lines.
19. **P19. Simpler `blast.nf`.** Ids: C-nextflow-37 (unverified-low), C-ponytail-23 (nit, confirmed). Where: `inst/nextflow/modules/blast.nf:36-45, 119-194`. Problem: six `sqlite3` calls read one options row (with `RB` reused for two meanings), and the `blastn` command is written five times. Fix: one `read ... <<< "$(sqlite3 ...)"` and one shell `run_blastn()` function. About 40 lines.
20. **P20. Trim `export_bundle()` and the export wrappers.** Ids: C-ponytail-25, C-data-34, C-data-37 (nit, confirmed). Where: `R/export.R:24-67, 161-191, 532-538, 702-721, 741-750`. Problem: options no caller passes (`template`, `run_table2asn`, `export_opts`, `label`, `ids`, `itsx_cpu`), `submission_risks()` (tests only), one-caller wrappers that each open their own connection (see DB-M4), and a manual begin and commit. `export_project()` is used only by the vignette. Fix: drop the unused options and wrappers, use `dbWithTransaction()`, and replace the vignette's scripted export with the app's Export tab. About 80 lines.
21. **P21. One triage settings row.** Ids: C-ponytail-26 (nit, confirmed). Where: `R/db_schema.R:60-79`, `R/triage_opts.R`, `R/init_db.R:221, 631`. Problem: named triage presets never go beyond "default", and the seed repeats the rule defaults, so a new rule needs edits in three places. Fix: one row seeded from `triage_actions_default()`. About 25 lines.
22. **P22. One reference-config reader.** Ids: C-nextflow-39 (unverified-low), C-security-27, C-ponytail-27 (nit, confirmed). Where: `R/nf_command.R:97-121`, `R/reference_nf.R:96`. Problem: the same row is read three ways with copied NA handling. Fix: one `project_config_value(path, col)`. About 25 lines.
23. **P23. Build the variant clause at display time.** Ids: C-ponytail-28 (nit, confirmed). Where: `R/screen_variants.R:32, 281-289, 350-364`. Problem: the variant sentence is spliced into `screen_reason` and parsed back out, although `variant_*` columns exist. Fix: store only the columns and render the sentence in the UI. About 25 lines.
24. **P24. Drop unused scoring options.** Ids: C-ponytail-33 (nit, confirmed). Where: `R/select_paths.R:43-75, 100-152`. Problem: `with_staxids = TRUE`, `expected_length`, and `weights` are never set by a caller. Fix: drop the branch and arguments. About 20 lines.
25. **P25. One alert and one toast.** Ids: C-ponytail-30 (nit, confirmed), C-shiny-39 (unverified-low). Where: `R/app_notify.R:39`, `R/app_alerts.R`, `R/app_server.R:683`, `R/constants.R:361`. Problem: `mp_*`, `rp_*`, `notify()`, and 14 direct `sendSweetAlert` calls under a "one alert system" comment (`rp_toast` even uses a different duration); `next_step` duplicates `MP_NEXT`, and `lock_next$u` is stored but never read. Fix: keep one alert and one toast, use `MP_NEXT`. About 20 lines.
26. **P26. One exportable-units query.** Ids: C-shiny-38 (unverified-low), C-data-36, C-ponytail-31 (low, confirmed). Where: `R/export.R:544-570`, `R/app_server.R:1061-1086`. Problem: the exportable-units SQL is written three times and the SeqID rule twice. Why beyond upkeep: if one copy changes, the header preview checks different SeqIDs than export writes. Fix: `exportable_units(con)` and `unit_seqids()` in `export.R`, used by both. About 15 lines.
27. **P27. Merge `config.local` into `config.docker`.** Ids: C-nextflow-36 (unverified-low), C-ponytail-32 (nit, confirmed). Where: `inst/config.local`, `inst/config.conda` (stale comment). Problem: the two files differ only in comments. Fix: keep one template and alias `local` to it. About 32 lines.
28. **P28. `RiboPilot <- run_app`.** Ids: C-ponytail-35 (nit, confirmed). Where: `R/app_run.R:159-209`. Problem: `RiboPilot()` only forwards to `run_app()`, which carries unused golem arguments. Fix: alias one to the other and drop the unused arguments. About 15 lines.
29. **P29. Move the `dev/` script tests out of the package suite.** Ids: C-ponytail-37 (nit, confirmed). Where: `tests/testthat/test-refdb-build.R` (2,666 lines), `test-db-release.R` (833 lines). Problem: they shell out to `dev/refdb/*.sh`, which `.Rbuildignore` removes, so they all skip under check, and many only grep a script for a string. Fix: keep a few real packaging tests next to `dev/refdb` (run with bats or a small R runner). About 3,000 test lines.
30. **P30. One anchor caller.** Ids: C-science-38 (unverified-low), C-ponytail-38 (nit, confirmed). Where: `R/annotate.R:214, 318, 374`, `R/blast_genbank.R:13, 91-92`. Problem: four ways to call the rRNA tool with different fallback rules (so orientation and annotation can pick different tools for one scaffold), and titles are classified twice per hit. Fix: build everything on `rrna_caller()` plus `run_rrna_tool()` plus `parse_barrnap()`; derive `is_nuclear_rrna` from `hit_class`. About 15 lines (more with `annotate_seq`, counted in P8).
31. **P31. Declare Nextflow params once.** Ids: C-nextflow-38 (unverified-low). Where: `inst/nextflow/main.nf:35`, `nextflow.config`, `blast.nf`, `pipeline_r.nf`, `blast_ref.nf`. Problem: `blastdb`, `ref_dir`, and `taxdump_dir` defaults are declared in up to four places. Fix: declare all params once in `nextflow.config`. Not estimated.
32. **P32. Small repeats.** Ids: C-science-41 (unverified-low, the "failed assembly" UPDATE written twice at `R/pipeline_steps.R:256-260, 298-302`, about 8 lines), C-nextflow-41 (unverified-low, a redundant first `while` loop in the clip awk at `inst/nextflow/modules/coverage.nf:69`), C-refdata-33 (unverified-low, one-caller wrappers such as `tk_genus_species` and `tax_index_con`, an exported `shared_taxon`, a misnamed `bitscore_weighted_lca`, a counts cache re-read on every rung, and unused script options). Fix: inline or delete each. Not estimated.
33. **P33. Drop the hand-built HTML report fallback.** Ids: C-data-38 (nit, confirmed). Where: `R/export.R:474`, `R/export_report.R:106`. Problem: a second report used only when rmarkdown or pandoc fails, and it does not escape HTML (SEC-L4). Fix: skip the report with a message when rmarkdown is missing. Not estimated.
34. **P34. Launch Nextflow with processx in `run_assembly()`.** Ids: C-security-28 (nit, confirmed). Where: `R/run_assembly.R:75`, `R/build_reference_db.R:78, 101`, `inst/nextflow/bin/rp_ref_db.R:22`. Problem: `bash -c` with nested `shQuote`, and explicit `--api-key` arguments EDirect does not need. Fix: `processx::run("nextflow", args, wd = path, env = ...)`, and pass the key through the environment. A few lines.
35. **Repo clutter** (part of SEC-L3, C-ponytail-36): root driver scripts, `update-version.log`, BUG and HANDOFF notes, and `dev/run_test_assembly.sh`. About 700 lines.

---

## 6. Test suite results

All runs were on a copy of the repo (by code-tests).

| Run | Setup | Result | Time |
|---|---|---|---|
| Full suite | `testthat::test_local()`, clean PATH (no conda) | FAIL 0, WARN 9, SKIP 23, PASS 3,952 | 230 s |
| Tool-dependent files | conda `ribopilot` env first on PATH | 2 tests (4 expectations) FAIL in `test-rrna-kingdom.R`; the 20 BLAST tests that skipped above now run and pass | about 2 min |
| `R CMD check` with CI flags | `--no-manual --no-build-vignettes`, installed package | **1 ERROR** (`test-nextflow.R:1051`), SKIP 132, 2 NOTEs | 217 s |
| Function coverage | every namespace function traced (covr not installed) | 496 of 570 functions (87 percent) called at least once | about 4 min |

What this means:

- **No single environment runs the whole suite green.** Without conda on PATH, 20 BLAST tests skip; with it, 2 tests fail (TEST-L2).
- **The one `R CMD check` error is a real bug**, not a test problem: BLAST database provenance is lost in container runs (PIPE-M17). CI cannot see it because its runner has no Nextflow.
- **The 9 warnings are test hygiene** (TEST-L9), not code failures.
- **Slow files:** `test-nextflow.R` 87 s (six real `nextflow run` calls), `test-app-scaffold-uid.R` 68 s, of which about 45 s is the app's 5 s shutdown timer (TEST-L3).
- **Network and Docker are handled well.** No test touches the internet (downloads use `file://` or mocks) and Docker and Singularity are never called; command builders are tested with a fake `which`.
- **Biggest gaps:** the per-unit database drivers (`curate_one`, `reconcile_one`, `score_one`, `blast_scaffolds_local`), the fallback assembler, and the real WF1 channel wiring are never run by any test (TEST-M1, PIPE-M16). Several of the confirmed re-assembly bugs (DB-H1, SCI-H1, SCI-M5, DB-L1) were missed because tests set database state directly instead of going through the functions the app calls.

#### Medium

**TEST-M1. The per-unit pipeline drivers and the WF1 wiring are never exercised**
- Ids: C-tests-7. Severity: medium. Status: confirmed.
- Where: `R/curate_rrna.R:285` (`curate_one`), `R/pipeline_steps.R:16, 90` (`score_one`, `reconcile_one`), `inst/nextflow/bin/rp_*.R`, `tests/testthat/test-nextflow.R:37-155, 379-562`.
- Problem: tracing the whole suite shows `curate_one`, `reconcile_one`, `score_one`, and `blast_scaffolds_local` are never called, and no test runs `rp_curate.R`, `rp_reconcile.R`, `rp_score.R`, `rp_triage.R`, `rp_chimera.R`, `rp_coverage.R`, `rp_fallback.R`, `rp_ref_align.R`, or `rp_rederive.R`. No `.nf` process has a `stub:` block, so a stub run is impossible, and `main.nf` wiring is checked only by regex over its text. The Nextflow run tests use small synthetic harnesses, not the real `main.nf`.
- Why it matters: rename a column that `curate_one()` writes, and in the pipeline the curate step fails on every unit (and is ignored), curation silently stops, and the suite stays green. The same is true of a channel-shape change in `main.nf` that still matches the grep patterns.
- Fix: add database-level tests for `curate_one()`, `reconcile_one()`, and `score_one()` on a `new_db()` project. Add `stub:` blocks and one test that runs `nextflow run main.nf -stub-run` on a small project using the bundled CHIM00 and CHIM01 reads, checking that the database ends with assemblies, annotations, screen, and triage rows.

#### Low

- **TEST-L1. C-tests-14 (confirmed; was medium), C-nextflow-31 (unverified-low)** (low). `.github/workflows/R-CMD-check.yaml:27-34`, `.github/workflows/pkgdown.yaml`. Problem: CI runs `R CMD check` on Ubuntu with no Nextflow, sqlite3, or BLAST, so 132 tests skip, including every `blast.nf` end-to-end test, and no job parses `main.nf` with the pinned Nextflow. The pkgdown job's `BUILD SITE` commit-message condition can never be true for manual, release, or pull-request triggers, and a step name contains a non-ASCII emoji. Why: a regression in `blast.nf`, the local and remote fallback, `rp_env.R` loading, or `main.nf` syntax merges with green CI (as PIPE-M17 did). Fix: add a CI job that runs `testthat::test_local()` from the source tree with Nextflow and the tools installed (the ribopilot image works), plus `nextflow run inst/nextflow/main.nf -preview`; relax the pkgdown condition for non-push events and remove the emoji.
- **TEST-L2. C-tests-9** (low, confirmed; was medium). `tests/testthat/test-rrna-kingdom.R:70, 100, 111, 122, 136`. Problem: the fake-tool tests prepend their folder to PATH, so a real pybarrnap on PATH still wins. Why: with the conda env active (needed for the BLAST tests), 2 tests fail with no code change. Test-only; the author already knows. Fix: `withr::local_path(c(dir, '/usr/bin', '/bin'), action = 'replace')`, or mock `rrna_caller()`.
- **TEST-L3. C-tests-10** (low, confirmed; was medium). `R/app_lifecycle.R:21`, `tests/testthat/helper-app-fixture.R`. Problem: the app's 5 s shutdown timer runs after every `testServer()` session. Why: it adds about 45 s to the suite (12 cases take exactly 5 s), and pending `stopApp()` calls fire inside unrelated later tests. Fix: read the grace period from an option and set it to 0 in the test helper, or stub `register_app_lifecycle()` there.
- **TEST-L4. C-tests-11** (low, confirmed; was medium). `tests/testthat/test-export-units.R:61`, `R/export.R:548-549`. Problem: the test named "ignored and problematic units are excluded from the bundle" only sets `ignore = 1`; no test exports with `problematic = 'yes'` or `annotate_lock = 0`. Why: if a refactor drops either filter, problematic or unlocked units reach the GenBank bundle and the suite stays green. The filter is correct today. Fix: mark one unit problematic and another unlocked, and assert neither appears in the FASTA files or the summary.
- **TEST-L5. C-tests-13** (low, unverified-low). `tests/testthat/test-app-scaffold-uid.R:697`. Problem: the ITSx roll-up test runs a pasted copy of the app's SQL, not the app. Why: if the app's roll-up breaks, the test still passes. Fix: drive it through `shiny::testServer(app_server)` as the other tests in the file do.
- **TEST-L6. C-tests-15** (low, unverified-low). `tests/testthat/test-app-opts.R:106, 129, 137`, `test-app-scaffold-uid.R:431`, `test-screen-variants.R:470`, `test-screen.R:270`, `test-blast-opts-panel.R:106`, `test-nextflow.R:37-155, 379-562`, `test-nf-pending.R:50`. Problem: many tests grep source text instead of checking behaviour, and one regex stops at the first nested `)`. Why: a harmless reformat fails them, while a behaviour change that keeps the strings passes. Fix: use `testServer()` plus output checks, a stub run for Nextflow wiring (TEST-M1), and `getParseData()` instead of regex for call checks.
- **TEST-L7. C-tests-16** (low, unverified-low). `tests/testthat/test-chimera.R:352`, `test-itsx-check.R:203`. Problem: the only real-data chimera and ITSx checks need `RIBOPILOT_TEST_PANEL`, `RIBOPILOT_TEST_TAXDUMP`, and `RIBOPILOT_TEST_PROJECT`, which are documented nowhere. Why: they skipped in every run, so regressions against real data go unnoticed. Fix: document the variables (for example in `tests/README.md`) and add a small dev script that sets them, and run it before releases.
- **TEST-L8. C-tests-17** (low, unverified-low). `tests/testthat/test-app-utils.R:42, 80-106`, `test-refdb-build.R:348, 2051-2131`. Problem: fixed sleeps, elapsed-time limits under concurrency, and a vacuous `/tmp` check. Why: flaky on a loaded runner, and the negative sleep checks can hide a real regression. Fix: poll with a bound for positive cases, mock the launcher for negative ones, and give time limits more headroom.
- **TEST-L9. C-tests-18** (low, unverified-low). `tests/testthat/test-blast-local.R:29, 78, 93, 177`, `test-ref-store.R:217`. Problem: tests expect one warning where the code emits several, so 9 warnings leak on every run, plus raw `tar` error text. Why: a genuinely new warning is easy to miss in the standing noise. Fix: one `expect_warning` per expected warning or `suppressWarnings()` where the warning is not the point; consider one combined warning in `local_blast_db_arg()`.
- **TEST-L10. C-tests-26** (low, unverified-low). `tests/testthat/test-nextflow.R:735-776, 993-1034, 1082-1111, 1194-1235`. Problem: four `blast.nf` end-to-end tests rebuild the same fixture instead of using `mk_dbarg_project()`. Why: a schema change needs four edits, and a missed copy tests a stale schema. Fix: give `mk_dbarg_project()` a stub-body argument and reuse it.

#### Nit

- **TEST-N1. C-tests-20** (nit, unverified-low). `tests/testthat/test-reference-modal.R:143`. Problem: a comment says `reference_handlers()` is not mounted by `app_server`, but `R/app_server.R:625` mounts it. Why: misleads the next maintainer. Fix: update or delete the comment.
- **TEST-N2. C-tests-27** (nit, unverified-low). `tests/testthat/test-app-export.R:4-46`, `test-app-scaffold-uid.R:6-48`, `test-blast-parse.R`, `test-db-dist.R`, `test-write-atomicity.R`, `test-chimera.R`, `test-reference-build.R:15`, `test-reference-modal.R:15, 161`. Problem: the app DDL, BLAST table schema, and chimera fakes are copied across files, and positional INSERTs into `reference_config` break on any new column. Why: fixture copies drift apart and tests end up checking different schemas. Fix: move them into shared helper files (or build from `new_db()`) and use named-column INSERTs.

---

## 7. R CMD check results

From code-pkg (build and check) and code-tests (check with tests), both on copies of the repo:

| Check | Result |
|---|---|
| `R CMD check --no-manual --no-tests` | 0 errors, 0 warnings, 2 NOTEs: hidden files in the tarball (`.nextflow.log`, `.claude`, `.superpowers`), and no visible binding for global `y` in `plot_chimera` |
| same, plus `--as-cran` | 2 WARNINGs and 6 NOTEs in total. Package items: tests use `rlang::` but rlang is not in Suggests; non-standard top-level files (`RiboPilot_test_userAsmb.R`, `update-version.log`); Title not in title case; `new_db.Rd` usage lines wider than 90 characters; "Imports includes 25 non-default packages". Two more items come from the check machine, not the package (no qpdf, no clock check). |
| CI flags with tests (`--no-manual --no-build-vignettes`) | **1 ERROR** at `test-nextflow.R:1051` (PIPE-M17), SKIP 132, 2 NOTEs |
| roxygen 8.1.0 on a copy | NAMESPACE and `man/` match the committed files exactly (no drift), apart from the "Hand-maintained" header; 4 "@noRd must not be followed by any text" errors from glued doc blocks |
| Non-ASCII scan | none in `R/`, `inst/`, or `man/`; `README.md` has pandoc curly quotes |
| Tarball size | 1.43 MB; 0.86 MB with the missing `.Rbuildignore` entries |
| lintr | only real hit is the same global `y` |

Notes:

- The hidden-files NOTE appears only when building from the author's working tree; a fresh git checkout (as in CI) has no `.claude` or `.nextflow.log`. But the Docker image is built from the working tree, so those files do ship there (SEC-L3).
- lintr's `object_usage_linter` checks against the installed RiboPilot 0.0.0.9016 and reports about 150 false "no visible function" hits. Run it with the dev version installed.
- Related items: SEC-L3 (`.Rbuildignore`), PKG-L3 (rlang), PKG-L4 (glued doc blocks), PKG-L8 (curly quotes), PKG-L11 (global `y`), PKG-N2 (DESCRIPTION metadata), REF-L4 (`new_db.Rd` width), P12 (Imports count), PIPE-M17 (the ERROR).

---

## Appendix: refuted findings

- **C-data-8.** An ITSx check that could not run overwrites the pipeline's stored ITSx warnings with `not_run(...)`. (Intended and tested: `tests/testthat/test-itsx-check.R:290-310`, and the column help says export updates it.)
- **C-data-16.** Zero-feature and undecided-triage units are exported with no submission-risk flag. (The app's default final ITSx check flags them and lists missing elements.)
- **C-tests-12.** The scientific-defaults test checks its own copy of the reference filter. (Redundant, not a gap: `test-blast-parse.R` runs the real `ref_candidate_ok()`.)
