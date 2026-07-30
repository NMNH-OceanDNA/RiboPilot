# 14 - One workflow, two assemblers

Status: design agreed 2026-07-30. Answers doc 08 open question 1 ("is the
fallback workflow still supported?") and most of open question 2. Schema v16.

## 0. What this changes, in one paragraph

The assembly fallback (Option B2 reference-guided baiting, Option A de-novo)
stops being a separate Nextflow workflow and becomes a per-sample choice of
assembler inside WF1. `main_fallback.nf` is deleted. `main.nf` branches its
assembly step on a new `assemble_opts.assembler` column and merges the two
output channels back together, so everything from `score()` onward is shared
rather than duplicated. `assemble_switch = 5` is retired. `run_fallback()`
survives as a thin wrapper that sets the column and calls the normal assembly
run.

## 1. Why the fallback broke, and why repairing it in place would not hold

`main_fallback.nf` was a copy of the pre-pivot `main.nf` tail. When the
per-scaffold pivot (`b0f1b12`) rebuilt that tail, the copy was not updated.
Nothing detected this because no test reads `main_fallback.nf` and the path has
never been run end to end.

The damage is not one bug. Measured against the current `main.nf`:

1. **Coverage gets the wrong tuple.** `main_fallback.nf:43` joins `stitch.out`
   (`(id, fasta)`, `pipeline_r.nf:128`) with `(id, reads)` and hands the
   resulting 3-tuple to `coverage`, which declares two (`coverage.nf:14-15`).
   `main.nf:52` joins off `score.out`, which is a bare `val id`
   (`pipeline_r.nf:108`), giving a clean 2-tuple. This is the one doc 08
   records, and it is the least of them.
2. **No annotation at all.** `main.nf:60-81` runs `query_sample_units`,
   `annotate_tools`, `reconcile`, `blast_ref_align` and `curate`.
   `main_fallback.nf` includes none of them. A rescued sample reaches the
   Assemble tab with unit rows seeded and empty. `main.nf:5-12` exists to
   explain why that is unacceptable: "nothing found" and "not looked at yet"
   render identically at the keep/drop moment.
3. **No screen, no triage.** `main.nf:89,94-103` run both. The fallback runs
   neither, so a rescued unit gets no contamination verdict and no keep/drop
   proposal.
4. **`stitch` is a hazard in the middle of the chain.** `main_fallback.nf:40-41`
   calls it only to relay an ID. Its output is `optional: true`
   (`pipeline_r.nf:128`), so a sample whose `stitch_selected()` returns NULL
   silently vanishes from BLAST as well.
5. **No depth gate.** `main.nf:39` filters on `params.minDepth`;
   `main_fallback.nf` never declares it, yet `R/nf_command.R:35` passes it.
6. **Samples strand.** `R/run_fallback.R:57-59` sets `assemble_switch = 5`
   before launching and never resets it. Nothing restores it on failure, and
   `assemble_switch = 5` is invisible to both `query_samples`
   (`pipeline_r.nf:29`, `IN (1,4)`) and `nf_pending_count()`
   (`R/nf_command.R:49`), so a stranded sample is never picked up again by an
   ordinary assembly run and never appears in a pending count.

Repairing these in place would mean re-deriving `main.nf:45-103` inside a second
file and then keeping the two in step by hand forever. The drift already
happened once, silently, and the only reason it was caught is that someone read
both files side by side. **The duplication is the defect.** Removing it is the
fix; repairing the copy is not.

## 2. Why folding is cheap

The two assembly processes already have identical contracts. This is not
something the design has to arrange; it is already true.

| | `assemble` (`modules/assemble.nf`) | `fallback_assemble` (`modules/fallback.nf`) |
|---|---|---|
| output | `tuple val(id), path("${id}_assembly_*.fasta"), emit: paths` (`:14`) | the same (`:17`) |
| publishDir | `${launchDir}/${params.publishDir}/${id}/assemble` (`:6`) | the same (`:10`) |
| errorStrategy | `ignore` (`:7`) | the same (`:11`) |
| header format | `>ID.path.scaffold topology` | the same (`R/assemble_fallback.R:110`) |

`modules/fallback.nf:1-4` already calls itself "a drop-in replacement for the
GetOrganelle `assemble` process". So the branch point is the assembly step and
nothing else; the merge point is `score()`.

**Correction to doc 08.** Doc 08 line 1304 says `assemble_fallback.R` "writes
`assemblies` through a path none of this traverses". It writes no database table
at all. It writes FASTA files (`R/assemble_fallback.R:113`), and the shared
`score` process loads them through `load_assemblies()`
(`inst/nextflow/bin/rp_score.R:6`, `R/pipeline_steps.R:216`) exactly as the
GetOrganelle path does. The headers it emits are `>ID.1.k`, which is already the
correct per-scaffold shape. Assembly loading and unit seeding were never broken.

**The fallback does not need stitching.** It emits per-scaffold contigs already.
`stitch` appears in its chain only as an ID relay, a job `score.out` does better
and without the optional-output hazard. Folding therefore removes the last live
caller of the stitch machinery, which is why this document also answers most of
doc 08 open question 2 (see section 7).

## 3. Routing: `assemble_opts.assembler`

New column on `assemble_opts`, values `getorganelle` (default), `bait`,
`denovo`. Nullable, treated as `getorganelle` when NULL so every existing
project keeps its behaviour without being rewritten.

Routing lives in `assemble_opts` rather than in a workflow flag or a switch
value because that is where every other per-sample behavioural choice already
lives, and because it records which assembler produced a given assembly. A
workflow-level `--assembler` flag cannot mix rescued and new samples in one run
and leaves no trace in the database.

**`assemble_switch = 5` is retired.** Rescued samples queue at `1` like
everything else. This is the point of the design rather than a side effect: the
stranding in section 1 item 6 becomes impossible rather than patched, because
there is no longer a queue state that the ordinary queries do not see.

The `bait` / `denovo` distinction is the existing `method` argument
(`R/assemble_fallback.R:139-147`, `modules/fallback.nf:14`), moved from a
workflow parameter to a per-sample column value. Its meaning does not change.

## 4. The workflow change

`query_samples` (`pipeline_r.nf:7-32`) gains `COALESCE(ao.assembler,
'getorganelle') AS assembler` in its select list and its header line.

`main.nf` branches after building `asm_in`, and merges immediately after:

```
asm_branched = asm_in.branch {
    go: it[7] == 'getorganelle'
    fb: true
}
assemble(asm_branched.go.map { id, reads, cpus, mem, go, seeds, labels, _asm ->
    tuple(id, reads, cpus, mem, go, seeds, labels) })
fallback_assemble(asm_branched.fb.map { id, reads, cpus, mem, _go, _seeds, _labels, asm ->
    tuple(id, reads, cpus, mem, asm) })

score(assemble.out.paths.mix(fallback_assemble.out.paths))
```

Everything from `score()` down is untouched. That is the whole point: the
annotate, screen and triage tails that `main_fallback.nf` lacked are not
re-implemented, they are simply reached.

Note the branch is written `fb: true` rather than testing for `bait`/`denovo`
explicitly, so an unrecognised value routes to the fallback assembler rather
than silently vanishing from both branches. A `branch` with no matching
condition drops the item without error, which is exactly the failure mode this
document exists to stop repeating.

`fallback_assemble`'s input tuple already carries `method` as its fifth element
(`modules/fallback.nf:14`), so passing the column value straight through needs
no process change.

## 5. What happens to `run_fallback()`

It stays exported and keeps its signature. No breaking API change, and the
vignette's Step 5 keeps working. It becomes:

1. Resolve candidates as it does now (`fallback_candidates()`,
   `R/run_fallback.R:10-16`) unless `ids` is given.
2. Point those samples at an `assemble_opts` preset whose `assembler` is the
   requested method, creating the preset if absent.
3. Set `assemble_switch = 1` (not 5).
4. Call the ordinary assembly run.

Step 2 is the only part with a real decision in it. Editing the sample's
existing preset in place would silently change the assembler for every other
sample sharing that preset, which is the same class of defect doc 08 records at
`R/app_opts.R:347-351`. So `run_fallback()` derives a new named preset per
(source preset, method) pair instead, copying the source preset's other values.

The completion message at `R/run_fallback.R:70` still says "assemble -> score ->
stitch -> BLAST + coverage". It is corrected to describe the real WF1 tail.

## 6. Schema v16

Additive and idempotent, following the existing pattern
(`R/migrate_db.R`, guarded `ALTER`, shared DDL constants referenced by both
`R/init_db.R` and the migration so a fresh database and a migrated one cannot
diverge).

- `ALTER TABLE assemble_opts ADD COLUMN assembler TEXT DEFAULT 'getorganelle'`,
  guarded on the column not already existing.
- `R/init_db.R:211-223` gains the same column in its `CREATE TABLE`.
- **One data repair, recorded because it is not purely additive.** Any row
  sitting at `assemble_switch = 5` is stranded by construction under the old
  code and unreachable under the new. The migration sets those rows to `1`. This
  is the only migration in the project that rewrites a value rather than adding
  a structure, so it is called out here rather than buried: leaving them at 5
  would strand them permanently once the query that recognised 5 is deleted.

## 7. What this does and does not settle in doc 08 open question 2

**Settled.** After this change nothing in any shipped workflow calls the stitch
machinery. `query_fallback_samples` (`pipeline_r.nf:33-53`) loses its only
caller and is deleted with `main_fallback.nf`.

**Not settled, and deliberately left alone.** `stitch_scaffolds()`,
`stitch_one()`, `stitch_selected()`, `order_scaffolds_by_content()`,
`rp_stitch.R`, the `stitch` process, the `stitched` table, `path_scores`,
`score_project()` and `annotate_project()` all remain. Two of them are exported
public API. Removing them is a separate decision with its own breaking-change
question, and doc 08 open question 2 says explicitly to decide the bundle
together rather than piecemeal. This document removes the last *reason* they are
there without removing them.

Recorded for whoever takes that up: `path_scores` is still written by the live
`score` process on every WF1 run (`R/pipeline_steps.R:37`) and read by the
Shiny app (`R/app_server.R:121-122,830,842-857`), so it is not dead even after
this change. Only its stitch-side reader dies.

## 8. Testing

The existing suite catches this bug class by reading workflow files as text
(`tests/testthat/test-nextflow.R:34-60` does this for `main.nf`). That is the
right instrument here, since the defect was a channel arity mismatch that no R
test could see.

- Assert `main.nf` branches on `assembler` and that both branches reach
  `score` via `mix`.
- Assert `main.nf` is the only workflow file referencing `fallback_assemble`,
  and that no workflow file references `stitch`.
- Assert `query_samples` selects `assembler`.
- Migration test: a v15 database gains the column, defaults to
  `getorganelle`, and any `assemble_switch = 5` row is moved to 1. Idempotent on
  a second run.
- `run_fallback()` unit test against a temporary project: it creates a derived
  preset rather than mutating the shared one, queues at 1, and leaves other
  samples on the source preset untouched.

This is the first test coverage the fallback path will ever have. Note what
tests cannot reach: `bait_reads()` and `run_spades()` shell out to `bowtie2` and
`spades.py` and both degrade to silent no-ops when the tools are absent
(`R/assemble_fallback.R:22-24,47-49`). The assembler itself stays unvalidated by
the suite, and this document does not claim otherwise. Validating it needs a
real run, which is the honest reason to keep the machinery rather than delete
it: the divergent invertebrate clades that block doc 08 section 5's validation
gate are exactly the case B2 baiting exists to serve.

## 9. Out of scope

- Deleting the stitch machinery, `path_scores`, `score_project()` or
  `annotate_project()` (section 7).
- Any change to what the fallback assembler actually does. `bait_reads()` is a
  single bowtie2 pass with no MITObim-style iterate loop, where
  `02_ASSEMBLY.md:118-131` specifies iteration. That gap is real, predates this
  work, and is left as it is.
- Exposing the assembler in the options UI. `R/app_opts.R:52` builds the
  assemble-options modal and would need a new control. Worth doing, but it is a
  UI task and the column is settable from `run_fallback()` without it.
- Validating the fallback assembler end to end on real reads (section 8).

## 10. Documentation to reconcile

- `dev/design/07_SCAFFOLD_SELECTION.md` states that the only remaining stitch
  caller chain runs through `main_fallback.nf`. After this change that is no
  longer true and the sentence needs amending. That text is currently on the
  `design-debt-stitching` branch (PR 9), so this reconciliation is sequenced
  after it merges.
- `dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md` open question 1 is answered
  and should record it; open question 2 is narrowed, not answered. Its line 1304
  claim about `assemble_fallback.R` is wrong and is corrected in section 2.
- `dev/design/02_ASSEMBLY.md:155-159` describes B2 as an "automatic fallback"
  triggered by a QC gate. It is not automatic; `run_fallback()` must be called
  by hand. Record the gap rather than implying it works.
- `vignettes/RiboPilot.Rmd:206-216` Step 5 keeps working unchanged, but its
  prose should stop implying a separate workflow.
