# User-supplied assemblies in RiboPilot

Date: 2026-09-01
Branch: `feat/user-assemblies`
Status: approved design, ready for implementation plan

## Problem

RiboPilot only accepts raw reads. A user who already has a genome assembly
(whole genome, or a single nuclear ribosomal contig pulled out of one) has no
way in. MitoPilot solves the equivalent problem with a parallel project type
(`new_project_userAsmb()`, `init_db_userAsmb()`, `find_mito()`), which doubles
the surface area of init, the database, and the app.

RiboPilot takes the other road: user assemblies are an input option on the
existing project type, not a second project type. There is one `new_project()`,
one database schema, one app, one Nextflow workflow.

## Scope

In scope:

- An `Assembly` column in the mapping file plus an `assembly_path` argument.
- Finding nuclear rDNA contigs inside a supplied assembly with barrnap.
- Making raw reads optional per sample.
- The downstream consequences of having no reads.

Out of scope:

- Migrating existing projects. RiboPilot is pre-release; no schema migrations
  are written for old databases.
- Any change to how read-based samples assemble.

## Design

### 1. Input

`new_project()` gains `assembly_path = NULL`. The mapping file may carry an
`Assembly` column naming one FASTA per sample. A row with a non-empty
`Assembly` is a user-assembly sample; its `R1`/`R2` may be blank.

Validation happens at project creation, not four hours into a run:

- `Assembly` values present but `assembly_path` unset is an error.
- Named assembly files must exist under `assembly_path`, gated by a
  `check_assemblies` argument that mirrors the existing `check_reads`.
- `check_mapping_reads()` is relaxed: a row with an assembly and blank reads
  passes; a row with neither reads nor an assembly is still an error.
- A row may carry both an assembly and reads. The assembly is used; the reads
  are used for coverage only.

The `.config` template already carries an `ASMB_DIR` token, currently pinned to
the `"NA"` sentinel. It is filled with `assembly_path`. For container
executors, `assembly_path` joins the bind mounts alongside the project and raw
directories.

Whole genome versus single rDNA contig is not a user-facing distinction.
Everything supplied is scanned; a single rDNA contig passes the scan trivially.
The opt-out is a per-sample `find_rrna` option in `assemble_opts` (default on,
editable in the app like every other option), for the case of an ITS-only
fragment that carries no barrnap-detectable anchor.

### 2. Finding the contigs

A new Nextflow module `find_rrna.nf` runs in place of `assemble` and
`fallback_assemble` for user-assembly rows.

- The supplied FASTA is split with `splitFasta` so a draft genome is scanned in
  parallel and never read whole.
- Each chunk is scanned with the sample's configured rRNA tool, reusing the
  `rrna_tool` and `rrna_kingdom` settings already in `assemble_opts` and the
  same tool-selection shell logic `annotate_tools` uses.
- The merged GFF goes to a new R function `select_rrna_contigs()`.

`select_rrna_contigs(gff, lengths, ...)` keeps a contig when it carries an 18S
or a 28S hit. A 5.8S hit alone is not enough: 5.8S is short, and a lone hit is
as likely to be a degraded pseudogene copy as a real repeat unit. Kept contigs
are ranked by total rRNA bases, best first, and capped (default 20) so a
repeat-rich assembly does not push hundreds of near-identical contigs into the
app. It returns the kept contigs plus one evidence row per contig scored.

The module writes `${id}_assembly_1.fasta` using the exact header contract the
assembler emits (`>id.path.scaffold topology`, produced today by
`format_assembly.sh`). Kept contigs become path 1, scaffolds 1..N. Because the
contract is identical, `score()` and the entire tail of WF1 (BLAST, reference
fetch, annotation, reconcile, curate, screen, triage) need no change.

When nothing is found the sample fails the way a failed assembly already does:
`assemble_switch = 3` with an explanatory note, and its seeded annotate rows
switch off.

### 3. Evidence

A new `rrna_candidates` table records why each contig was kept or dropped:

```
ID, contig, length, rrna_features, rrna_bases, rank, selected, reason, time_stamp
PRIMARY KEY (ID, contig)
```

Written by the module as an upsert keyed on `(ID, contig)`, so a re-run
overwrites rather than duplicating. Surfaced in the app so a user can see that
their assembly was searched and what was found, rather than only seeing a
failure note.

### 4. Reads optional

- `preprocess` and the `min_depth` gate apply only to rows with reads. A
  reads-less row never enters that channel.
- `coverage` is skipped for a reads-less sample, so no `coverage` rows exist
  for it. The app and the export report must render a missing coverage row as
  "not measured" rather than as zero depth. This is the main audit item during
  implementation.
- Path scoring already tolerates absent depth: `score_scaffolds()` sets the
  depth term to 0 when every candidate's depth is NA. No change.
- Triage rules do not read coverage depth (verified: `triage_rules.R` has no
  depth reference), so auto-triage is unaffected.
- The `preprocess` row for a reads-less sample is seeded with blank read names
  and left untouched by the pipeline.

### 5. App

- The finder result appears as an assemble note, which is also what says the
  sample came from a supplied assembly. There is no separate source column: the
  note already carries it.
- `find_rrna` joins the assemble options modal.

Two things the existing assemble details modal cannot show, because both
concern contigs that never became assembly units, are folded into it rather
than given a viewer of their own. Both appear only for a sample that supplied
an assembly:

- A "Your contig" column in the scaffold table, from a new
  `assemblies.source_contig`, so a unit can be traced back to the user's own
  file. `load_assemblies()` fills it from the trailing field of the fasta
  header, which the assembler was already writing and the loader was
  discarding.
- A collapsed "Contigs not carried forward" section listing each dropped
  contig, its length, the rRNA it did carry, and why it was left behind. Hidden
  entirely when a sample dropped nothing.

The rest of `rrna_candidates` is deliberately not surfaced: for a kept contig
the real annotation step re-derives the features properly, and the scaffold
table already shows those alongside length and rank.

## Testing

- `select_rrna_contigs()` against fixture GFFs: whole genome with one rDNA
  contig, a single supplied rDNA contig, nothing found, 5.8S-only, and more
  candidates than the cap.
- `new_project()` / `new_db()` with a mapping carrying `Assembly` and blank
  reads: succeeds, seeds the right rows, and still rejects a row with neither.
- Mapping validation: assembly named but `assembly_path` unset, and an assembly
  file that does not exist.
- Existing tests continue to pass unchanged, which is the check that read-based
  projects are untouched.
