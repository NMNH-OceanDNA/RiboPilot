# Plan: decommission WF2 (spec: dev/design/15_WF2_DECOMMISSION.md)

Branch `wf2-decommission`, cut from `main`. Task-by-task, suite green between
tasks, mutation-test every new test.

## 1. Pipeline: the annotate-only arm

- Replace `query_annotate_samples` with `query_annotate_units` in
  `modules/pipeline_r.nf`: units queued for annotation, on non-ignored
  scaffolds, EXCLUDING samples this run will assemble.
- `main.nf`: run the new query, mix its units into the existing `annotate_tools`
  input, leave the rest of the chain alone.
- Delete `inst/nextflow/main_annotate.nf`.

## 2. R API

- `nf_command()`: drop the `workflow` argument, always build the WF1 command.
- Delete `run_annotate()`, its Rd, its NAMESPACE line.
- `nf_pending_count()`: count assembly work + annotate-only work, one number.
- `run_assembly()`: docs mention WF2; update.

## 3. App

- Assemble tab Options group gains `annotate_opts` and `curate_opts`; the
  Assemble-mode query must select them.
- Annotate tab: remove the Update button; keep review, marks, lock, trim.
- `run_modal()`: one mode.

## 4. Tests

- `query_annotate_units` excludes a sample queued for assembly (the race).
- It includes a queued unit on an unlocked sample (the dropped `assemble_lock`
  gate).
- No `main_annotate.nf`, no `query_annotate_samples`, no `run_annotate` export.
- `nf_pending_count()` counts both kinds.
- The two pieces of defensive wiring the new arm relies on (`ifEmpty('ready')`
  on the reference sentinel, `remainder: true` on the triage join) are pinned as
  text, with a comment saying why.
- Assemble-mode table exposes the annotation and curation option cells.

## 5. Verify for real

- Run WF1 on a project with a queued unit and no queued assembly: only the
  annotate chain should run.
- Re-run with `-resume` on an already-assembled project and confirm assembly is
  cached, which the spec currently asserts from the process definitions rather
  than from measurement.

## 6. Docs

Doc 15 status, doc 01's workflow section, `HANDOFF.md`, and doc 03 if it names
WF2.
