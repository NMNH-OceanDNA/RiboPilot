# Run reports and .runs/ layout (port from MitoPilot)

Branch: `feat/run-reports` (off `devel`). Source: MitoPilot commit 12efba75
("Add Nextflow run reports and .runs/ file layout", devel). MitoPilot borrowed
RiboPilot's end-of-run report and extended it; this brings the extensions back.

## What RiboPilot has today

- End-of-run report for in-app runs only: `nf_run_report()` (R/nf_run.R:317),
  shown once by `rp_alert(run_report_ui(rep))` (R/app_run.R:325) and saved as
  `nextflow_logs/<kind>_<ts>_report.txt`.
- Log parsing in `nf_log_facts()` (R/nf_run.R:214), `task_err_tail()`,
  `nf_exit_cause()`, `run_failed_samples()`, `failure_hint()`. Same regexes
  MitoPilot uses.
- Files scattered: `jobs/`, `logs/`, `nextflow_logs/`, root `.nextflow.log*`
  (console `run_assembly()`, no `-log`), `.reference/build_<ts>.log`.
- Cluster jobs (`build_submit_script()`, `NF_EXIT_TRAILER` already echoes
  `[RiboPilot] nextflow exited with status $status`) are never reported on.
- No way to look at past runs.

## What MitoPilot added (to port)

1. `.runs/` layout:
   ```
   .runs/jobs/<base>.sh
   .runs/logs/<base>.log            scheduler log, holds exit trailer
   .runs/nextflow/<base>.nextflow.log
   .runs/reports/<base>_report.txt
   .runs/reports/<base>_report.dcf  one record of RUN_META_FIELDS
   ```
   `<base>` = `<workflow>_%Y-%m-%d_%H-%M-%S` (RiboPilot already has
   `job_basename(kind)` in this shape).
2. Report metadata saved as `.dcf` next to the text, so past runs can be listed.
3. Lazy `sync_run_reports()`: any finished log under `.runs/nextflow/` with no
   valid dcf gets a report. Finished = exit trailer in scheduler log, else
   (if job script has trailer, wait) log tail has Goodbye / Session aborted /
   Launcher ERROR. Covers cluster jobs and saved scripts.
4. `list_run_reports()`: dcf rows plus "unfinished" rows for logs without one,
   newest first.
5. Run Reports button + modal: table (New, Started, Duration, Result pill,
   Samples, Failed samples, Failed tasks, Launched App/Job); row click shows
   report text, Copy report, Open report location.
6. "New run reports available (N)" confirm on startup, refresh, and after an
   in-app run; "Open Run Reports" button opens the modal.
7. Small parser upgrades: `tag_sample()` (tag to sample ID, longest known ID
   wins), `entry` and `samples` in facts, `nf_log_time()` year fix,
   `fmt_duration()`, `read_exit_trailer()`, failed-samples window with an
   upper bound (`until`).
8. Workdir browser reads both old and new log locations.

## RiboPilot adaptations

- Workflows: `pipeline` (main.nf, one workflow covering Assemble + Annotate)
  and `reference` (main_reference.nf BUILD). MitoPilot's per-panel
  assemble/annotate split does not apply. Run Reports modal lists both, with a
  Workflow column (or filter), not per-panel.
- Button placement: next to `workdirs_btn` (R/app_ui.R:111). Show on Assemble
  and Annotate (and Reference for build runs), or always.
- Failed samples: keep RiboPilot's `run_failed_samples()` (assemble table,
  `assemble_switch = 3`, assemble + find_rrna notes) and `failure_hint()`;
  add `until` bound. Add annotate-unit failures only if a switch=3 state exists
  for units (check).
- `time_stamp` in RiboPilot is wall-clock "now" when the row is written
  (preprocess.nf, rp_find_rrna.R), not workflow start. `>= since` and
  `<= ended + 1` window still works; overlapping runs are already impossible
  (one run per project lock).
- Tags: `${id}` and `${id}.${upath}.${uscaf}`; `tag_sample()` default regex
  fits as is.
- Strings: "RiboPilot run report", `ribopilot_version`, trailer already
  `[RiboPilot]`.
- Keep RiboPilot helpers: `rp_alert`, `mp_confirm` (add `cancel_label`),
  `run_report_ui`, `rpCopyEl`, `open_path`, reactable `rt_*` helpers. Check
  for `rt_pill` / `mp-pill` CSS; add if missing.
- Reference builds: `start_reference_build()` also tees to
  `.reference/build_<ts>.log` and sets `reference.fail_reason`. Leave that
  path alone; only move its `-log` into `.runs/nextflow/reference_<ts>...`.
- Console `run_assembly()`: add `-log .runs/nextflow/<base>.nextflow.log` so
  console runs show up too (it is still a supported entry).
- `find_workdirs()` (R/app_utils.R:32): add `.runs/nextflow/`, and fix
  `TaskHandler\\[` to `\\w*TaskHandler\\[` so grid runs match.

## Steps

1. **Paths.** `run_dir()`, `run_basename()`, `run_log_base()`,
   `run_base_time()`. Point `JOB_DIRS` / `job_path()` at
   `.runs/{jobs,logs,nextflow}`; reports to `.runs/reports/`. Tests.
2. **Parser.** Extend `nf_log_facts()` with `entry`, `samples`, `started`,
   `ended`, `readable`, `finished`, `fatal`, `tasks$sample`. Add
   `tag_sample()`, `nf_log_time()`, `read_exit_trailer()`. Port MitoPilot
   fixture logs (rewrite process names to RiboPilot ones) into
   `tests/testthat/fixtures/runs/`.
3. **Report object.** Make `nf_run_report()` also return `meta` (RUN_META_FIELDS
   with `ribopilot_version`, `workflow`, `launch`). Add `write_run_report()`
   (txt + dcf), `dcf_ok()`, `run_is_finished()`, `sync_run_reports()`,
   `list_run_reports()`. Keep the existing alert at end of in-app runs.
4. **Launch paths.** `nf_run_window` (pipeline + reference) and
   `build_submit_script` / `build_reference_submit_script` / save-script use
   the shared `<base>` and `.runs/`. Write report via `write_run_report()` in
   the end-of-run path (finish, Stop, Close, session end). Console
   `run_assembly()` gets `-log`.
5. **Run Reports module.** New `R/app_run_reports.R`: `run_reports_ui()`,
   `run_reports_server()`, `check()` exported via `session$userData`, startup
   check after first flush (skip if a startup alert is showing), check after
   in-app run end. Wire in app_ui.R / app_server.R.
6. **Workdir browser.** Read `.runs/nextflow/`, fix grid handler regex.
7. **Docs.** NEWS "Runs" section; vignette/doc mention of `.runs/` if logs
   folders are documented; .gitignore `.runs/` in project template if any.
8. **Verify.** `devtools::test()`, then drive the app headless on the test
   project (one in-app run, one saved-script run via bash, check modal).

## Decisions (2026-09-30: all recommendations taken)

1. **Old log folders** (`jobs/`, `logs/`, `nextflow_logs/`): leave in place
   and only read old `nextflow_logs/` in the workdir browser (MitoPilot way),
   vs. move them into `.runs/` on open. Rec: leave in place.
2. **Reference builds in Run Reports:** include (rec), or pipeline only.
3. **New-reports popup:** port as is (rec), or only the "New" pill in the
   table with no popup.
4. **Button visibility:** always visible in the toolbar (rec, since one
   pipeline covers several panels), vs. only on Assemble/Annotate.
