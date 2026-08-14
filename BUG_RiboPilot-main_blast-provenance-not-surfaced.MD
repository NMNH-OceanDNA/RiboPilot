# BUG: BLAST provenance is not surfaced; failed remote BLAST reads as "never run"

Repo: RiboPilot
Branch: main
Verified against commit: f61912b
Date: 2026-08-12
Status: diagnosed, not fixed. No code written for this. No design approved.

---

## One-line summary

RiboPilot assigns species from two different BLAST searches with very different
authority, never labels which one a result came from, and the Assemble details
modal hides the weaker one entirely. A scaffold whose remote BLAST failed is
reported to the user with the same message as a scaffold that never had BLAST
turned on.

---

## Background: there are two BLAST searches, not one

This is the thing to understand first. It is not obvious from the UI.

### 1. Local seed BLAST ("quick hint")

- Function: `blast_scaffolds_local()`, `R/pipeline_steps.R:51`
- Called from: `score_one()`, `R/pipeline_steps.R:40`, wrapped in `tryCatch`,
  best-effort, silently no-ops if blastn or the seed file is missing
- Database: built on the fly from the GetOrganelle **seed** FASTA
  (`R/pipeline_steps.R:21`, `make_blast_db(local_ref(seeds_db))`).
  For the fish test project that is `actinopterygii_rrna_seeds.fasta`,
  3857 sequences, headers of the form `>MW829490.1_Osteoglossum_bicirrhosum`
- Species name comes from parsing that FASTA header, not from any taxonomy lookup
- Writes: `assemblies.blast_accession`, `.blast_species`, `.blast_pident` only
  (`R/pipeline_steps.R:70`). Does NOT write `blast_qcovs` and does NOT write
  `blast_source`
- Runs offline in seconds, right after assembly

**Critical caveat:** these seeds are the same sequences GetOrganelle used to bait
the assembly. So this answers "which seed does my scaffold most resemble", not
"what species is this". It is a hint, not an identification.

### 2. Remote GenBank BLAST ("the real answer")

- Nextflow process `blast_genbank`, megablast against NCBI `core_nt`
- Parsed by `rp_blast.R` -> `parse_blast_candidates()`
- Writes: `blast_hits` (all HSPs), `blast_screen` (best hit + lineage + taxids +
  qcovs + `is_nuclear_rrna`), `blast_ref_candidates`, and overwrites the
  `assemblies.blast_*` columns, this time also setting
  `blast_source = 'core_nt'` (`R/blast_genbank.R:151`)
- Slow, network dependent, and in practice unreliable (see "Why this matters now")

So `blast_source` is the only field that distinguishes the two. It is set to
`'core_nt'` by the remote path and left as empty string by the local path.

---

## The findings

### Finding 1: the details modal hides the local hit entirely

`output$assemble_scaffolds` (`R/app_server.R:1036`) sets
`defaultColDef = reactable::colDef(show = FALSE)`, so only columns with an
explicit `show = TRUE` render.

Every BLAST-ish column that is turned on reads from `blast_screen`:
`screen_species` ("Top GenBank hit"), `screen_accession` ("Accession"),
`screen_pident` ("%ID"), `screen_verdict` ("Screen"), `shares`, `odd_one_out`.

The local seed columns `blast_species`, `blast_accession`, `blast_pident`,
`blast_qcovs`, `blast_source` ARE fetched by the modal's query
(`R/app_utils.R:150`) and are sitting in the data frame, but have no colDef, so
`defaultColDef` hides them.

Result: when remote BLAST has not produced a `blast_screen` row, the modal shows
a blank BLAST section while holding a species call in memory.

### Finding 2: the "no evidence" message misattributes the cause

`output$unit_screen` (`R/app_server.R:1174`) branches on one condition only:

```r
if (is.na(r$screen_accession)) {
  ... "No GenBank screen for this scaffold yet (run WF1 with BLAST enabled)."
}
```

That single message is shown for at least three distinct states:

| Actual state | Message shown | Correct? |
|---|---|---|
| BLAST never run | "run WF1 with BLAST enabled" | yes |
| `blast_opts.run_blast = 0` | "run WF1 with BLAST enabled" | yes |
| Remote BLAST ran and failed / timed out / still queued | "run WF1 with BLAST enabled" | no, misleading |

In the third case the user is told to enable something that is already enabled.
There is nothing in the UI that says the search was attempted and did not come
back.

### Finding 3: the two panels point at each other

`output$unit_candidates` (`R/app_server.R:1211`) renders, when
`blast_ref_candidates` is empty:

> "No nuclear rDNA reference matched this scaffold. See the unrestricted screen above."

But in this failure mode the screen above is itself the empty-state message from
Finding 2. Two panels, each deferring to the other, no evidence in either.

### Finding 4 (the root of 1-3): `blast_source` is never rendered

`blast_source` is selected in `R/app_utils.R:150` and is written only by the
remote path (`R/blast_genbank.R:151`). It appears in no colDef and no UI string
anywhere in `R/`. So nothing in the app ever tells a reviewer whether a species
call came off a 3857-sequence seed file or off core_nt.

Consequence, visible today: the Assemble tab table shows the seed hit under
"BLAST Hit" / "Accession" / "BLAST %ID". Click into that same sample's details
and the modal says there is no GenBank screen. Same scaffold, one view names a
species, the other says nothing is known. Neither view says where its number
came from.

Related minor inconsistency: a `blast_qcovs` colDef exists in the samples table
(`R/app_server.R:213`) but the Annotate-mode query (`R/app_server.R:160`) does
not select `blast_qcovs`, only species/accession/pident.

---

## Why this matters now

Observed in a live full-project test run at
`/home/dmacguig/Documents/DJMCompBio/RiboPilot_test/full_project_docker`
(RiboPilot 0.0.0.9005 image, Docker profile, 10 fish samples).

Remote BLAST is the pipeline's slowest and least reliable stage. In this run,
three `blastn -remote` searches sat in NCBI's queue for over two hours with zero
output. Independently confirmed as an NCBI-side backlog, not a RiboPilot bug: a
420 bp control query submitted directly to the BLAST URL API returned
`RTOE = 30` (a 30 second estimate) and was still `Status=WAITING` 25 minutes
later.

State of that project's database while BLAST was stalled:

```
blast_source   n scaffolds   blast_qcovs NULL
------------   -----------   ----------------
(empty)                 27                 27      <- seed BLAST only
core_nt                  1                  0      <- real screen
```

So 27 of 28 scaffolds carried a confident-looking species name from the seed
file, with no indication anywhere in the app that it was a hint rather than an
identification, and the details modal for all 27 said no GenBank screen exists.

Two calls from that run that a reviewer could easily over-trust:

- `SRR25465041` -> *Stomias affinis* at 94.9% identity, seed BLAST only
- `SRR33849683` -> two scaffolds landing on two different flatfish
  (*Paralichthys squamilentus* 99.5%, *Atheresthes stomias* 98.8%), seed BLAST
  only. Reads as a contamination signal but is just as likely a gap in the seed
  file.

---

## How to reproduce

Easiest path, no waiting on NCBI:

1. Open any project where WF1 has run through the `score` step but where
   `blast_screen` has no row for some scaffold. Setting
   `blast_opts.run_blast = 0` before the run produces this directly.
2. Assemble tab: the sample shows a species under "BLAST Hit" with an accession
   and a %ID.
3. Click the details icon for that sample, then click the scaffold row.
4. Expand "GenBank evidence for the highlighted scaffold".
5. Observed: "No GenBank screen for this scaffold yet (run WF1 with BLAST
   enabled)", and the candidates panel below telling you to consult the empty
   panel above. The seed hit from step 2 appears nowhere in the modal.

To reproduce the specifically misleading variant (BLAST enabled, remote failed),
run WF1 with `run_blast = 1` while offline, or kill the `blast_genbank` task.
The `blast_genbank` process retries and can exhaust its budget, leaving no
`blast_screen` row and the same "run WF1 with BLAST enabled" message.

---

## What is NOT the bug

Do not "fix" these, they are working as designed:

- **`assemble.blast_accession` / `blast_species` / etc. are empty.** Those are
  legacy per-sample columns. Both BLAST paths correctly write the per-scaffold
  `assemblies` table instead, consistent with the per-scaffold annotation pivot.
- **`blast_qcovs` is NULL for seed hits.** `blast_scaffolds_local()` deliberately
  does not compute coverage.
- **The seed BLAST being weak evidence.** That is inherent and fine. The problem
  is that it is not labeled as such.
- **Remote BLAST being slow.** NCBI-side, confirmed independently.
- **`SRR22396893` scaffolds having `ignore = 1`.** They are 469 bp and 379 bp,
  correctly filtered.

---

## Suggested direction (NOT approved, discuss before building)

Repo owner prefers narrow-scope fixes. Smallest change that removes the
misinformation, roughly in priority order:

1. **Make the empty-state message distinguish the states.** `output$unit_screen`
   needs to know whether BLAST was disabled (`blast_opts.run_blast = 0`), never
   run, or attempted and unresolved. `run_blast` is reachable via
   `assemble -> blast_opts`. Minimum viable: two messages instead of one.

2. **Surface the seed hit in the modal as an explicitly labeled fallback**
   rather than hiding it. The data is already in `scaffold_rows()`; this is
   colDefs plus a label, not new plumbing. It must read as "closest seed
   sequence", never as a GenBank identification.

3. **Render `blast_source` wherever a BLAST species is shown**, or derive a
   visible badge from it. This is the root fix. Note the local path currently
   leaves it as empty string rather than writing something like `'seeds'`, so
   this likely needs a one-line addition to the UPDATE at
   `R/pipeline_steps.R:70`.

4. **Fix the dead-end pointer** in `output$unit_candidates` so it does not refer
   the reader to an empty panel.

Open questions for the repo owner:

- Should a seed-only call appear in the Assemble tab's "BLAST Hit" column at all,
  or should that column be reserved for core_nt results with the seed hit living
  only in the modal?
- Is a seed-only species name allowed to feed triage rules or the contamination
  screen? Worth checking what `R/triage_rules.R` and `R/screen.R` consume before
  changing what is displayed, so display and logic do not diverge.
- If `blast_source` starts being written as `'seeds'`, existing projects will
  have empty strings. Repo is pre-release with no migration system by design, so
  the UI probably needs to treat empty as "unknown / legacy" rather than
  backfilling.

---

## Verification status

- All file:line references above were checked against commit `f61912b` on `main`.
- Findings 1-4 were traced by reading the code, NOT by clicking through a running
  Shiny app. The reasoning is straightforward (`defaultColDef(show = FALSE)` plus
  an absent colDef equals a hidden column; a single `is.na()` branch equals a
  single message), but nobody has watched the modal render this state.
  **Confirm in the app before writing a fix.**
- The database state and the NCBI queue behaviour in "Why this matters now" were
  measured directly on the live test project, not inferred.
