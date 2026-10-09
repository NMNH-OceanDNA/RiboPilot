# Contig recheck workflow and review copies

Date: 2026-10-08
Branch: feat/contig-recheck (off devel at a25fd25)
Status: approved design, not yet built

## 1. Goal

- A second Nextflow workflow reruns the checks on contigs the user trimmed or
  edited, one contig at a time.
- An Update button on the Annotation Review step launches it, the same way the
  Assemble step's Update launches the main pipeline.
- A trim or a gene boundary edit queues that contig automatically.
- The current tables always hold the original contig and its check results.
  The Assemble step shows only these.
- New review tables hold the working copy. A contig nobody edited reads its
  original values. An edit marks the copy's checks "not run".
- Restore returns a contig to its original values.
- Export reads contig content and check results only from the review layer.
- Export no longer runs ITSx; it reports the stored ITSx warnings.

## 2. Today, and why it changes

- Trim (`trim_assembly`, R/pipeline_steps.R:572) and Restore
  (`restore_assembly`, :728) rewrite `assemblies`, `annotations`, `coverage`,
  `chimera`, and `chimera_windows` in place. The original survives only in
  `assembly_backup`, snapshotted at the first trim. Boundary nudges and partial
  toggles write `annotations` directly (R/app_server.R:3063-3065, 3077-3079)
  with no backup at all.
- After an edit only `curate_one` reruns, inside the app. ITSx, chimera,
  reference alignment, BLAST, and screen results go stale, and nothing is
  queued. Comments at inst/nextflow/main.nf:20-23 and
  inst/nextflow/modules/pipeline_r.nf:47-50 claim a trim requeues the contig;
  no code does.
- Because the tables are shared, the Assemble step already shows edited state
  (Checks column, contig table, sequence viewer).
- Export reruns stock ITSx and writes the result over
  `annotate.itsx_warnings` (R/export.R:835-849, R/itsx_check.R:445-462).
- History: dev/design/15_WF2_DECOMMISSION.md (2026-07-31) folded the old
  re-annotation workflow (WF2) into main.nf and removed the Annotate Update
  button, because two entry points confused users and the copy drifted. The
  workflow added here has a different job. It never re-annotates; it rechecks
  the review copy through the same R check functions and the same tool shell
  script the main pipeline uses, so there is no second copy to drift.

## 3. Decisions (user, 2026-10-08)

| # | Decision |
|---|---|
| D1 | The recheck reruns every check that works offline, plus a fresh local BLAST: feature counts, barrnap model coverage, ITSx chimera flag, ITSx check, local BLAST, contamination screen, reference alignment, curation, chimera. No reads, no NCBI. Variant grouping and triage are not rechecked; they compare all contigs of a sample and belong to the Assemble step. |
| D2 | Lock is refused while a contig's checks are not done. Exclude & Lock still works. Export leaves out contigs whose checks are not done, as a backstop. |
| D3 | Review copies are made on first touch, and database views fall back to the originals for contigs without a copy (approach A). |
| D4 | Export drops its own ITSx run and reports the stored ITSx warnings. |
| D5 | Restore deletes the copy: trims, boundary edits, and feature-level 5'/3' choices go with it. The unit Partial mark, Reviewed, Exclude, notes, group, and lock stay. |
| D6 | Curation no longer reruns instantly after an edit; flags show "not run" until Update. |
| D7 | State on Annotation Review controls the recheck queue and offers only Ready to run and On hold. |

## 4. Data model

### 4.1 Originals

`assemblies`, `annotations`, `annotate` (check columns), `coverage`, `screen`,
`chimera`, `chimera_windows`, `blast_hits`, `blast_searched`, `blast_screen`,
`blast_ref_candidates`, `blast_ref_alignment`, and `path_scores` keep their
current layout. The main pipeline writes them, plus the existing in-app
exceptions (taxonomy rescreen, Keep/Drop and Triage, option sets). Annotation
Review edits, the recheck workflow, and export never write them.

`annotate` also keeps the user's review state: `annotate_switch` (the main
pipeline's re-annotation queue), `annotate_lock`, `reviewed`, `problematic`,
`partial` (unit mark), `annotate_notes`, `export_group`,
`export_time_stamp`, `annotate_opts`, and `curate_opts`. These are review
decisions, not contig data, so they stay where they are and export keeps
reading them for gating and grouping.

New `annotate` columns, written by the stock ITSx check:
- `itsx_origin TEXT`: ITSx's one-letter origin code as reported in its
  extraction results, `''` when ITSx reports none. The display maps it
  through `ITSX_ORIGINS`.
- `itsx_scores TEXT`: ITSx's group-score string, verbatim.
NULL in either means never stored.

Removed: the `assembly_backup` table and `ensure_assembly_backup()`.
Originals never change now, so there is nothing to back up.

### 4.2 Review tables

A contig (unit `(ID, path, scaffold)`) has a review copy exactly when it has a
`review_units` row.

`review_units`, primary key `(ID, path, scaffold)`:

| Column | Meaning |
|---|---|
| sequence, length | the working sequence |
| trim_start, trim_end | kept window in original coordinates (1 and the original length when untrimmed) |
| edited | 1 once the copy was trimmed, nudged, or given a 5'/3' choice; 0 for a copy made only by State |
| check_switch | 2 done, 1 not run (queued), 4 running, 3 failed, 0 on hold |
| check_note | why the last recheck failed or did not finish |
| rev | starts at the current time in milliseconds when the copy is made and goes up by one on every invalidation, so a value never repeats for a unit |
| checked_at | unix seconds of the last finished recheck |
| n_features, missing, chimeric, n_flags, curate_flags, itsx_warnings, itsx_origin, itsx_scores | per-contig results, same meaning as in `annotate` |
| blast_accession, blast_species, blast_pident, blast_qcovs, blast_evalue, blast_source | top-hit summary, same meaning as in `assemblies` |
| time_stamp | last write |

Mirror tables, with the same columns and keys as their originals:

- `review_annotations`, plus `partial_start_user INTEGER` and
  `partial_stop_user INTEGER` (section 4.4)
- `review_blast_hits`
- `review_blast_screen`
- `review_screen`
- `review_chimera`, `review_chimera_windows`
- `review_ref_alignment`

Not copied:
- Coverage: read from the original and cut to the kept window (section 4.5).
- Reference candidates, sequences, and annotations, and `blast_searched`:
  shared, read-only for the recheck; it uses the original best reference,
  which is already cached.
- Triage, path scores, and rDNA candidates: Assemble decisions.

### 4.3 Views

Each view returns the review rows for contigs that have a copy, and the
original rows for contigs that do not. The fallback keys on `review_units`, so
a copy whose features were all trimmed away returns no features rather than
the original ones.

- `review_units_v`: `review_units`; otherwise `assemblies` joined with
  `annotate`, with `trim_start = 1`, `trim_end = length`, `edited = 0`, and a
  derived `check_switch`: the original's `annotate_switch`, except that 2
  without a stored ITSx value (`itsx_warnings IS NULL`) reads as 1, not run.
  Both kinds also expose `has_copy` (0 or 1), `topology`, and the original
  length `orig_length`.
- `review_annotations_v`: `review_annotations`; otherwise `annotations` with
  NULL user partial columns.
- `review_screen_v`, `review_chimera_v`, `review_chimera_windows_v`: the same
  pattern.

Annotation Review and Export read only these views. The recheck reads and
writes the review tables directly, since a queued contig always has a copy.

### 4.4 Partial ends

Today curation, trims, and the 5'/3' toggles all write
`annotations.partial_start/stop`, and curation fills only NULLs
(R/curate_rrna.R:350-358), so a stale curation-set partial can never be
cleared.

On review copies:
- Trim cut ends and the 5'/3' toggles write `partial_start_user` and
  `partial_stop_user`. The user's choice always wins.
- Curation writes its own verdict into `partial_start/stop`, replacing it on
  every recheck.
- The effective value is `COALESCE(partial_*_user, partial_*)`. Every reader
  on the review layer uses it: export `partial_ends` (R/export.R:213), the
  feature table and toggles (R/app_server.R:2465-2497), the sequence viewer
  (R/seqview.R:74-90), and the header tokens (R/export_opts.R:146-168).

Originals keep today's behavior.

### 4.5 Coverage and BLAST hits on a trimmed copy

A trim only removes sequence from the ends, so original per-position data
maps onto the copy with an offset of `trim_start - 1`:
- Coverage: the original `coverage.depths` and `coverage.clips` cut to the
  window. Clips are no longer lost on a trim.
- Mapped hits (used when the recheck does no fresh search, section 7.3): the
  original HSPs that overlap the window, clipped to it and shifted. A clipped
  HSP gets `aln_len` set to its clipped span, and `sstart` and `send` set to
  NULL. `qcovs` is recomputed per subject as the union of its clipped spans
  over the new length, rounded; that formula reproduces the stored `qcovs` on
  four test databases. Clipped HSPs keep their pident and bitscore, which is
  approximate.

## 5. Lifecycle of a review copy

### 5.1 First touch

These Annotation Review actions make a copy when the contig has none: trim
(single or batch), boundary nudge, 5'/3' toggle, and State. The marks (Lock,
Reviewed, Exclude, unit Partial), notes, and group write `annotate` only and
never make a copy.

`review_copy(con, unit)` does nothing when a copy exists. Otherwise it first
deletes any rows the unit still has in the seven mirror tables (leftovers of
a discarded recheck, section 5.4), then copies the unit's original rows into
all eight review tables. The copy starts with `trim_start = 1`,
`trim_end = length`, `edited = 0`, `rev` = the current time in milliseconds,
and the check state the view showed for the original. Trim, nudge, and the
toggles then set `edited = 1`.

`review_copy`, `review_invalidate`, and `review_drop` use a savepoint rather
than a new transaction, so they also work inside an open transaction (the
upgrade step, the State handler).

### 5.2 Annotation Review writes

| Action | Writes | Check state |
|---|---|---|
| Trim, single or batch | narrows the window; `sequence` cut from the original by the window; features shifted, clamped, or dropped as today; cut ends get `partial_*_user = 1` | invalidated |
| Boundary nudge (Edit mode) | `review_annotations` start, stop, length | invalidated |
| Save (Edit mode) | ends the edit session; no inline curation | unchanged |
| Reset or Discard (Edit mode) | writes the session's starting coordinates and user partials back | invalidated again when a coordinate changes |
| 5'/3' toggles | `partial_*_user` | unchanged; no check reads partial ends |
| Restore assembly, or batch Restore | `review_drop`: deletes the unit's rows from every review table | the view shows the original again |
| State: Ready to run | makes a copy if needed, then invalidates | invalidated |
| State: On hold | `check_switch = 0` | on hold |
| Lock, Reviewed, Exclude, unit Partial, notes, group | `annotate`, as today | unchanged |

State refuses locked contigs ("unlock first"), like the other edits. State:
Ready to run is also refused for a contig without a copy whose original is
queued for the main pipeline (`annotate_switch` 1 or 4) or was never
annotated (`n_features IS NULL`), with "press Update on the Assemble step
first"; a copy made then would freeze features the main pipeline is about to
replace.

`review_invalidate(con, unit)` sets `check_switch = 1`, adds one to `rev`,
clears the copy's check columns, and deletes its rows in
`review_blast_hits`, `review_blast_screen`, `review_screen`,
`review_chimera`, `review_chimera_windows`, and `review_ref_alignment`. It
also clears the curation-owned `status`, `flags`, and `partial_*` on the
copy's features. `n_features` and `missing` are recomputed at once from the
copy's features rather than cleared, since they are plain counts.

Batch Restore offers every contig with a copy, not only trimmed ones.

### 5.3 Events from outside Annotation Review

The main pipeline's per-unit bin scripts that write a unit's check results
(`rp_reconcile.R`, `rp_blast.R`, `rp_searched.R`, `rp_curate.R`,
`rp_chimera.R`, `rp_screen.R`) call `review_follow_original(con, unit)` for
each unit they write:
- an unlocked copy with `edited = 0` is deleted, so the contig follows its
  original again;
- an unlocked edited copy is invalidated; its edits are kept;
- a locked copy is left alone. A locked contig is finished and the recheck
  skips it, so invalidating it would strand it.

The shared check functions do not call it, so the recheck and the in-app
refreshes below never trigger it. Edited copies do not pick up new features
from a re-annotation; Restore adopts the new original.

| Event | Original | Review copy |
|---|---|---|
| Re-assembly changes the sequence (`reset_unit_evidence`) | reset as today, and also clears `itsx_warnings`, `itsx_origin`, `itsx_scores`, `chimeric`, `n_flags`, and `curate_flags` | deleted: the review tables join `UNIT_EVIDENCE_TABLES` (R/pipeline_steps.R:239-241), which the orphan sweep also uses |
| Re-annotation (`reconcile_one`, any route) | annotations rewritten and `chimeric` written from the new positions file, as today; in the same transaction `itsx_warnings`, `itsx_origin`, and `itsx_scores` are cleared until the stock check stores them | `review_follow_original` |
| BLAST, curation, chimera, or screen rerun by the main pipeline | as today | `review_follow_original` |
| Annotation or curation options change (`opts_apply`) | `annotate_switch = 1`, as today | `review_follow_original` on the sample's units, at once |
| Taxonomy change (`rescreen_samples`) | `screen_sample` as today (screen, then variants), plus `curate_one` and the ITSx origin refresh (section 6.3) for every unit | done copies, locked ones included, refreshed in place afterwards: screen re-scored from the copy's hits (carrying the new variant facts), curation, origin flag. Copies at 1 or 4 get `rev` raised by one so an in-flight recheck result is discarded and the contig reruns with the new lineage |
| Keep/Drop and Triage | as today | unchanged; a dropped contig hides its copy with it |

### 5.4 Concurrent edits and runs

- `query_recheck` emits each contig's `rev` and sets `check_switch = 4` in
  the same sqlite3 call. `recheck_tools` writes `unit.fasta` only while the
  copy still has that `rev`, and emits nothing otherwise. The `recheck`
  process passes the `rev` on to `rp_recheck.R`; `recheck_unit()` stops at
  once if the stored `rev` differs, and never reads `rev` afresh.
- The last write sets `check_switch = 2` only
  `WHERE rev = <passed rev> AND check_switch IN (1, 4)`. An edit during the
  run raises `rev`, On hold sets 0, and Restore or re-assembly deletes the
  copy; in each case the write changes nothing and the result is discarded.
- Intermediate writes are not guarded. Rows a discarded recheck already wrote
  either sit under a not-run state, which readers never show as results
  (section 8), or have no `review_units` row, which the views ignore and the
  next `review_copy` deletes. The next recheck replaces them.
- After an in-app run ends, contigs still at 4 go back to 1 with
  `check_note = 'recheck did not finish'`. After a cluster run they stay at 4
  ("running") until the next Update picks them up again.
- The main pipeline writes originals and edits write copies, so a run no
  longer needs to refuse edits. The batch Trim and Restore "in the active
  run" skip (R/pipeline_steps.R:794-796) is removed.
- The two workflows can still overlap on one contig: the main pipeline's
  annotate-only arm re-annotates queued contigs of assembly-locked samples
  (pipeline_r.nf:57-59, 79-81), and cluster jobs take no run lock. The run
  lock (in-app runs), `review_follow_original`, and the `rev` guard keep their
  writes apart.

## 6. Done rule, lock gate, and ITSx states

### 6.1 Done

A contig's checks are done when `review_units_v.check_switch = 2`:
- With a copy: the recheck set it, and the recheck always stores a non-NULL
  `itsx_warnings`.
- Without a copy: `annotate_switch = 2` and `annotate.itsx_warnings IS NOT
  NULL`.

### 6.2 Lock gate and export backstop

- Annotation Review Lock (toolbar `lock_rows`, row icon, and `lock_unit` in
  the details window) refuses contigs that are not done. Exclude & Lock is
  exempt. The refusal names the next step for each state:

  | State | Message |
  |---|---|
  | no copy, original queued for the main pipeline | "Press Update on the Assemble step first." |
  | not run (1) on a copy | "Press Update first." |
  | not run (1) without a copy | "Set State to Ready to run, then press Update." |
  | running (4) | "Wait for the running Update to finish." |
  | failed (3) or on hold (0) | "Set State to Ready to run, then press Update." |

- The Assemble step's Lock is unchanged and reads originals.
- `export_bundle` drops contigs that are not done from its selection before
  its "No contigs to export" check, and stops with a message about unfinished
  checks when none remain. `export_left_out` lists them with the reason
  "checks not finished", and the console message and Export done dialog
  count them. The filter applies to the records written, not to the sibling
  count that names records (`exportable_units`, R/export.R:714-718), so record
  names never change between exports.

### 6.3 ITSx states

- `''` clean; a flag list; `not_run(<reason>)` means the check ran but could
  not decide; NULL means never checked.
- `not_run` counts as done and is lockable. The Checks column shows it, and
  NULL, as an amber "ITSx not checked" item, never as clean. The ITSx column
  keeps its grey "not run" pill for `not_run` and shows NULL as grey "not
  checked".
- A crash is `not_run`, not "ITSx found nothing". `itsx_check_one` counts a
  run as finished only when its summary file contains "ITSx run finished" and
  its captured stderr has no `ERROR` or `FATAL` line. Otherwise it stores
  `not_run(ITSx did not finish)`. ITSx exits 0 on fatal errors, so its exit
  status is not used.
- When `itsx_check_one` itself throws, `rp_reconcile.R` stores
  `not_run(ITSx check failed: <message>)` instead of leaving NULL. NULL then
  only ever means "never checked".
- Origin flag on a taxonomy change: the origin flag
  (R/itsx_check.R:306-314) moves into a helper,
  `itsx_origin_flag(origin, scores, lineage)`, which takes the one-letter
  code. `rescreen_samples` (called by both lineage writers,
  R/tax_search.R:79-80, 99) recomputes only the origin part of
  `itsx_warnings` from the stored `itsx_origin` and `itsx_scores`, for
  originals and done copies. Rows with no stored origin (checked before this
  change) and `not_run` rows are left alone.

## 7. The recheck workflow

### 7.1 Files

- `inst/nextflow/main_recheck.nf`
- `inst/nextflow/modules/recheck.nf`
- `inst/nextflow/bin/rp_recheck.R`
- `inst/nextflow/bin/rp_rrna_tools.sh`: the barrnap/pybarrnap, tuned ITSx,
  and stock ITSx calls, moved out of `annotate_tools`
  (inst/nextflow/modules/pipeline_r.nf:150-187) so the main pipeline and the
  recheck run identical commands. It records each tool's exit status in a
  status file and captures stock ITSx stderr into `unit.itsx_check.err`,
  which the existing `unit.itsx_check*` glob already passes on.

### 7.2 Processes

All are `cache false`, run per contig, and tag `ID.path.scaffold`.

1. `query_recheck` (local executor, sqlite3, retried on a lock race like the
   other `query_*` processes) selects contigs with
   `review_units.check_switch IN (1,4)`, `annotate_lock = 0`, kept
   (`ignore = 0`), on a sample with `assemble_lock = 1` and
   `taxonomy_lock = 1`, and sets them to 4. The taxonomy lock keeps a recheck
   from running without a lineage. It emits ID, path, scaffold, `rev`, cpus
   and memory from `annotate_opts`, and a BLAST mode: `local` when the
   sample's `blast_opts.run_blast = 1` and the contig's original search was
   local (`blast_searched.source = 'local'`), `map` otherwise.
2. `recheck_tools` writes the copy's sequence to `unit.fasta` (section 5.4)
   and runs `rp_rrna_tools.sh` with the unit's `annotate_opts`. It gets the
   same out-of-memory retry and memory closure as `annotate_tools`
   (nextflow.config:59-62). A tool failure does not fail the task; it is
   recorded in the status file.
3. `recheck` resolves the local panel with `rp_dbarg.R` and the labels
   database as the `chimera` process does (pipeline_r.nf:264-284), then calls
   `rp_recheck.R`, which calls `recheck_unit()`.

### 7.3 `recheck_unit()`

In order:

1. Stop if the copy is gone or its `rev` differs from the passed one.
2. Hits:
   - `local`: run blastn megablast on the copy's sequence against the panel
     with the sample's `max_target_seqs` and `extra_opts`, `-max_hsps 5`, and
     the same output format as inst/nextflow/modules/blast.nf:123; retry with
     dc-megablast when there are no hits or the best `qcovs` is under 25.
     The R constants for these settings are tested against blast.nf so they
     cannot drift. Never calls NCBI. If the local search fails, fall back to
     `map`.
   - `map`: mapped original hits (section 4.5). When the original has no
     hits, the copy keeps the `blast_*` summary it copied at first touch (for
     example a seed-database hit).
   Results go to `review_blast_hits`, then the top screening hit to
   `review_blast_screen` and the summary to the `review_units.blast_*`
   columns. Reference candidates are not re-derived, and `blast_searched` is
   not written.
3. Screen: `screen_one` on the review tables with `force = TRUE`. In review
   mode it carries the original's `variant_*` values into the new row;
   `screen_variants` never runs on review tables.
4. `n_features` and `missing` from the copy's features.
5. Model coverage: barrnap hits from `recheck_tools` are matched to the
   copy's barrnap-sourced features by element, strand, and largest overlap,
   and only `model_coverage` is updated. A feature barrnap no longer finds
   gets 0, so curation's "truncated" rule fires instead of being skipped.
   When barrnap itself failed, the values stay as copied.
6. ITSx chimera flag from the tuned run's positions file
   (`itsx_is_chimeric`); left as copied when the tuned run failed.
7. Stock ITSx check: `itsx_check_flags` on the copy's sequence, features, and
   the sample's lineage; stores `itsx_warnings`, `itsx_origin`, and
   `itsx_scores` under the section 6.3 rules.
8. Reference alignment: `compute_ref_alignment` on the review tables against
   the original best candidate.
9. Curation: `curate_one` on the review tables.
10. Chimera: `chimera_unit` on the review tables, with mapped coverage and
    clips, the copy's hits, and the copy's ITSx chimera flag.
11. Set `check_switch = 2` and `checked_at` under the guard (section 5.4).

A check that cannot run stores `not_run(<reason>)` in its own field and does
not fail the contig. An error that escapes the steps sets `check_switch = 3`
and `check_note`, under the same guard.

The existing check functions gain a `review = FALSE` argument (`curate_one`,
`compute_ref_alignment`, `chimera_unit`, `screen_one`, `itsx_check_units`,
`unit_sequences`, `itsx_store_warnings`, `parse_blast_candidates`,
`derive_blast_products`). In review mode:

| Original | Review |
|---|---|
| `annotations` | `review_annotations` |
| `blast_hits`, `blast_screen`, `screen`, `chimera`, `chimera_windows`, `blast_ref_alignment` | the `review_*` mirror of each |
| `assemblies.sequence`, `length`, `blast_*` | `review_units` |
| `annotate` result columns (`n_features`, `missing`, `chimeric`, `n_flags`, `curate_flags`, `itsx_*`) | `review_units` |
| `annotate` review state (`annotate_opts`, `curate_opts`, locks), `samples`, `blast_ref_candidates`, `blast_ref_sequences`, `blast_ref_annotations`, `blast_searched` | unchanged, read-only |
| `coverage` | read through the window helper (section 4.5) |

`derive_blast_products` skips its `blast_ref_candidates` write in review
mode, and the recheck never calls `mark_unit_searched`. Review-mode writes
raise errors instead of swallowing them (today `itsx_store_warnings` swallows
UPDATE errors, R/itsx_check.R:399-402). Both workflows score with the same
code.

### 7.4 Launch from the app

- `nf_recheck_command(path)`: `main_recheck.nf` with `-c <project>/.config`,
  `-work-dir <project>/work_recheck`, `--db`, `--projectDir`, `--rawDir` when
  set, `--ribopilot_src`, `--taxdump_dir`, `--blastdb`, and the same container
  bind options for the panel and labels folders as `nf_command`. Launched
  from `<project>/.recheck`, like reference builds launch from `.reference`,
  so the main pipeline's `-resume` history is never disturbed. No `-resume`:
  every recheck process is uncached.
- `nf_recheck_pending(path)`: runs the same selection as `query_recheck`; a
  test keeps the two identical. When it is 0, "Nothing to update" also says
  how many contigs wait for the Assemble step's Update.
- Annotation Review toolbar: an Update button (`annot_run_btn`) at the end of
  `annot_ctrls` (R/app_ui.R:90-98), styled like Assemble's.
- Flow, as on Assemble: run lock check; "Nothing to update" when nothing is
  queued; a run window titled "RiboPilot recheck - N contig(s) to update";
  Run from App, or a job script on a cluster host; live progress; end-of-run
  report; table refresh.
- `run_modal` takes a prefix that goes on every input and output ID (`rc_`
  for the recheck) and omits the Resume box when the command has no
  `-resume`. The window uses `nf_run_window` with `pfx = "rc_"`,
  `kind = "recheck"`, `workflow = "recheck"`.
- On a cluster host the recheck gets its own Save and Submit observers, which
  use `job_basename("recheck")` and a job script that changes into
  `<project>/.recheck`, as `build_reference_submit_script` does for
  `.reference`.
- Run reports: base name `recheck_<timestamp>`; Run Reports labels it
  "Recheck"; failures are copies at `check_switch = 3`, or left at 4, in the
  run's time window; the stopped-run hint says to press Update on Annotation
  Review.
- The Assemble step's Update count and run are unchanged (main pipeline only).

## 8. App changes

Annotation Review and the Export step read the views. The Assemble step reads
originals.

Annotation Review:
- Table query (R/app_server.R:220-244): length, feature counts, BLAST
  summary, curation and ITSx columns from `review_units_v`; screen from
  `review_screen_v`; chimera from `review_chimera_v`; the State column from
  `review_units_v.check_switch`, with a contig without a copy whose original
  is queued shown as "waiting for Assemble Update".
- Check display is gated on the state. When a contig is not done, every check
  cell shows the state instead of a value ("not run", "running", "failed",
  "on hold", or "waiting for Assemble Update", info pill), and so does the
  Checks column. This also hides the leftovers of a discarded recheck
  (section 5.4).
- Checks column: `unit_checks(by = "unit", review = TRUE)`, with the ITSx rule
  of section 6.3. Today R/export.R:280 hides `not_run`.
- Details window: features from `review_annotations_v` with effective
  partials; sequence viewer with the copy's sequence and cut coverage; trim
  plot with cut coverage; chimera panel from `review_chimera_v`, or "Not run
  yet. Press Update." when not done; curation warnings gated on the state;
  reference alignment view from the copy's sequence and bounds; batch Trim
  planning from the copy.
- Handlers rewritten on the review functions: trim, Restore, batch Trim and
  Restore, nudges, partial toggles, Save, Reset, Discard, State, Lock.
- The inline `curate_one` calls go (R/app_server.R:2817, 2960, 2989;
  R/pipeline_steps.R:825).
- The lock check reuses the export risk path; on Annotation Review it reads
  the review layer, on the Assemble step the originals
  (`export_summary_df` and `unit_checks` gain a `review` argument).

Export step:
- Its table (R/app_server.R:204-218) reads `review_units_v`,
  `review_screen_v`, and `review_chimera_v`, with the marks from `annotate`.
- Contigs that are not done show "checks not finished, left out" and are
  excluded from the group count, the selected count, and the Export button's
  enablement (R/app_server.R:1213-1223, 1750-1775).
- Header previews (R/app_server.R:1244, 1294) read the review layer.

Assemble step:
- Reads originals, as today. The chimera panel shows the original check; the
  pre-trim juggling (`chimera_original`, R/app_modals.R:1020-1044) goes.
- The Checks column applies the same ITSx rule (section 6.3) to originals.
- The ITSx roll-up (R/app_server.R:186-190) never shows a tick when any kept
  contig is NULL or `not_run`, and it counts only kept contigs.

## 9. Export

- Every per-contig value comes from the views: sequence, length, topology,
  features with effective partials, `n_features`, `missing`, `chimeric`,
  `curate_flags`, `itsx_warnings`, the screen verdict, LCA, and variant facts,
  and the chimera label and reason. That includes the header tokens
  (R/export_opts.R) and the HTML report: inst/export_report.Rmd opens its
  own connection and calls `export_summary_df`, `feature_map_svg`,
  `report_elements`, and `gather_feature_seqs` directly, so those calls pass
  the review layer too.
- `unit_trimmed` (R/export.R:227) reads the window: left trimmed when
  `trim_start > 1`, right trimmed when `trim_end < orig_length`. The sequence
  search against the backup goes.
- The ITSx run is removed: `itsx_export_check`, `run_itsx_check`, and
  `itsx_check_command` (R/itsx_check.R), the `itsx_check` and `itsx_cpu`
  arguments of `export_bundle`, the "Final ITSx check" box and its progress
  text in the app, `<label>_itsx_check.csv`, and the `itsx_check/` folder.
- The summary CSV's `itsx_check` column carries the stored `itsx_warnings`,
  with NA written as `not_run(not checked)`.
- The Export done dialog shows counts of flagged, clean, and not checked
  contigs and lists the flagged ones.
- The README says where the flags come from, explains each `not_run`
  reason and what to do about it, and notes that ITSx calls can shift with
  batch size, so a whole-upload check at NCBI may differ on borderline
  contigs.
- Export writes no check result. It still writes `export_time_stamp` and
  export groups, as today.
- Contigs that are not done are left out (section 6.2).
- Copy FASTA on the Assemble step (`unit_fasta_text`) stays on the original.

## 10. Main pipeline changes

- `annotate_tools` calls `rp_rrna_tools.sh` (same commands as today).
- `reconcile_one`: the annotations rewrite and the `annotate` update run in
  one transaction; `chimeric` is written from the new positions file as
  today, and `itsx_warnings`, `itsx_origin`, and `itsx_scores` are cleared.
- `itsx_check_one`: finished-run detection and `not_run` rules (section 6.3);
  stores `itsx_origin` and `itsx_scores`.
- `rp_reconcile.R`: an error from `itsx_check_one` stores
  `not_run(ITSx check failed: <message>)`.
- `reset_unit_evidence` clears the extra columns listed in section 5.3.
- The per-unit bin scripts call `review_follow_original` (section 5.3).
- Stale comments fixed: inst/nextflow/main.nf:20-23, pipeline_r.nf:47-50,
  R/nf_command.R:6-7, R/app_ui.R:89.

## 11. Schema step 4 (version 4 to 5)

`new_db()` and the new step share the DDL constants. The step runs inside
`backwards_compatibility`'s transaction, so everything it calls uses
savepoints or transaction-free inner functions (`review_copy`,
`curate_one`, and the chimera snapshot write-back are each split into a body
and a transaction wrapper).

1. Create the eight review tables and five views.
2. `ALTER TABLE annotate ADD COLUMN itsx_origin TEXT` and
   `ADD COLUMN itsx_scores TEXT`.
3. Convert old edits. This is a safety net: none of the nine local test
   projects has any. A unit qualifies when its current sequence differs from
   its `assembly_backup` sequence (trimmed), or when it has
   `annotations.edited = 1`. A unit whose current sequence equals its backup
   and has no edited features (trimmed, then restored) gets no copy.
   - It gets a review copy of its current state, with `edited = 1`, the
     window found by locating the current sequence in the backup sequence
     (full length when there is no backup), and the partial values on edited
     features moved to the user partial columns.
   - When it was trimmed, its original is written back from the backup:
     sequence, length, annotations JSON, coverage depths and summary, and the
     chimera snapshot. When the snapshot is NULL, the current chimera rows are
     shifted back as `restore_assembly` does. Curation is rerun on the
     restored original and its `itsx_warnings` set to NULL, since the stored
     value may describe the trimmed sequence.
   - Unlocked copies are queued (`check_switch = 1`, checks cleared). Locked
     copies keep their values with `check_switch = 2`; a NULL `itsx_warnings`
     on one becomes `not_run(not checked before upgrade)`.
4. Each locked unit without a copy whose derived state is not done (for
   example a State set by hand on a locked row in an older version) gets a
   copy marked done, with `not_run(not checked before upgrade)` when its ITSx
   value is NULL. The upgrade message lists these units.
5. `DROP TABLE assembly_backup`.

The parity test (tests/testthat/test-backwards-compat.R:118-128) must also
compare views.

## 12. Testing

New:
- First touch, invalidate, Restore, and `review_follow_original`; views fall
  back for contigs without a copy, return the copy for touched ones, and
  return no features for a copy whose features were all trimmed away; a copy
  made after a discarded recheck starts clean.
- Trim on a copy: window arithmetic, features, user partials, original
  untouched. Nudge, toggle, Reset, Discard. State offers only Ready and On
  hold, refuses locked contigs, and refuses contigs waiting for the main
  pipeline.
- Lock gate on the toolbar and in the details window, with the per-state
  messages; Exclude & Lock exempt; the Export step's counts; export left-out
  reason; record names unchanged between exports.
- The `query_recheck` SQL extracted from modules/recheck.nf and run against a
  real schema, as tests/testthat/test-nextflow.R:163-234 does for
  `query_annotate_units`; `nf_recheck_pending` returns the same count.
- `recheck_unit` with fake tools, step by step; the `rev` guard discards a
  result after an edit, after On hold, and after Restore then a new trim;
  `not_run` paths; mapped hits (shift, clip, `aln_len`, `qcovs`); a `map`
  contig with no original hits keeps its seed summary; BLAST constants match
  blast.nf; leftover 4 goes back to 1 after an in-app run.
- `itsx_check_one`: an unfinished summary or an `ERROR` line gives `not_run`;
  an error gives `not_run(ITSx check failed: ...)`; origin code and scores
  stored; a lineage change recomputes the origin flag for originals and done
  copies and raises `rev` on queued copies.
- `reconcile_one` keeps writing `chimeric`; `reset_unit_evidence` deletes
  copies; locked copies are never touched by `review_follow_original`.
- Export reads the review layer, including inst/export_report.Rmd; trimmed
  ends come from the window; no ITSx run; the `itsx_check` column comes from
  stored values; dialog counts.
- Upgrade: a version 4 project with a trimmed, a trimmed-then-restored, a
  nudged, and a locked hand-stated contig ends as section 11 describes;
  parity holds with views.
- A real Nextflow run of `main_recheck.nf` with fake tools, as
  tests/testthat/test-nf-main-run.R does for main.nf, including a contig in
  `map` mode.
- The real-ITSx test (tests/testthat/test-itsx-check.R:220-268) ported to
  `itsx_check_one` on real stock output.
- A static check that export code (R/export.R, R/export_opts.R,
  R/export_report.R, inst/export_report.Rmd) and the Export step query read
  contig data only through the views.

Updated, because they pin the old design or code this design moves:
- test-app-opts.R:129-135 (no Update on Annotate), test-nextflow.R:124-143
  (the list of .nf files), test-nextflow.R:1576 and 1580-1591 (the OOM
  selector string and the barrnap command lines in pipeline_r.nf),
  test-rrna-kingdom.R:157-170 (kingdom translation inside `annotate_tools`),
  test-nf-pending.R, test-utils.R:58-76 (`open_project_db` call count).
- test-batch-trim.R, test-trim-coverage.R, test-chimera.R:433-456
  (`plot_chimera` after a trim), test-app-export.R (second trim, partial
  toggles, the ITSx input), test-export.R (trimmed ends),
  test-load-assemblies.R (backup reset), test-itsx-check.R (the export run;
  the `itsx_check_one` fixture at 297-314 needs a summary file), and the
  export tests that pass `itsx_check = FALSE`.
- helper-fake-tools.R: the fake ITSx writes a summary containing "ITSx run
  finished". helper-app-fixture.R: `fixture_fill_schema` builds the review
  tables, views, and new columns from the shared DDL and stops using
  `ASSEMBLY_BACKUP_DDL`.

The key new tests are mutation-checked: break the code on purpose and
confirm the test fails.

End to end in the app on the test project: trim, Update on Annotation Review,
Lock, Export. Pipeline R steps run the package installed in the container, so
use the conda executor or rebuild the image first.

## 13. Docs

- vignettes/Test-Project-Annotate.Rmd: an edit marks checks "not run"; Update
  on Annotation Review; the lock gate; what Restore restores.
- vignettes/GenBank-Submission.Rmd: no final ITSx check; where the flags come
  from; the batch-size caveat; the `not_run` legend.
- vignettes/Test-Project-Export.Rmd: no ITSx box.
- Link text: Test-Project-Assemble.Rmd:82-83, Curation-Flags.Rmd:44-45,
  Citing.Rmd:51.
- Tooltip at R/constants.R:149-151 drops "updated again by export".
- dev/design/01_ARCHITECTURE.md: the workflow section (two sample-level
  workflows now) and the claim that re-annotation keeps trims and edits.
- NEWS.md entry. A status note in dev/design/15_WF2_DECOMMISSION.md pointing
  here.
- Screenshots to retake (listed for the user): the Annotation Review toolbar,
  export-window.png, export-done.png.

## 14. Deferred and known limits

- No hint yet that a contig's original changed after it was edited; Restore
  adopts the new original.
- A contig without a copy reads as done as soon as its re-annotation and ITSx
  check finish, while reference alignment, curation, and chimera may still be
  running in the main pipeline. Locking in that window is possible today too.
- A crash of the tuned ITSx run still reads as "not chimeric"; the gap exists
  today.
- Version 4 projects whose taxonomy changed after annotation keep the old
  origin flag; none of the local projects is affected.
- Cluster jobs take no run lock; existing gap, unchanged. Contigs of a
  cluster recheck that died show "running" until the next Update.
- Existing NMNH Hydra projects keep their `.config`, so `recheck_tools` uses
  the default memory request there; new projects get the `annotate_tools`
  closure.
- The chimera window `rank_distance` (display only) is not refreshed on a
  taxonomy change.
- Variant grouping and triage are not rechecked (D1).
- Editing a shared option set does not queue other samples using it; existing
  behavior.

## 15. Risks

- Many readers move to the views. A missed one shows original data on
  Annotation Review or in an export. The static export check (section 12)
  and the not-done display gate catch most of them.
- More writers on one SQLite file: recheck tasks plus the app. The existing
  busy timeout and retries apply, and every recheck write is per contig.
- Container runs need the new bin scripts and R functions in the image.
