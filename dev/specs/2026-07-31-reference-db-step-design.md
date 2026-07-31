# In-app reference database step

Date: 2026-07-31
Branch: `refdb-picker` (off `main`)

## Problem

GetOrganelle seed and label databases are the single most important input to
assembly, and today they are built entirely outside the app by calling
`build_reference_db()` from the console with a hand-written GenBank query and a
hand-picked clade name. A user has to know that the step exists, know which
clade is appropriate, and know how to wire the result into `new_project()`.
Nothing in the app suggests a clade, warns that the chosen one is too narrow to
be useful or too broad to be specific, or handles a project spanning more than
one distant group.

Every ingredient for doing better is already present: taxonomy is resolved and
frozen per sample before assembly runs, the pinned taxdump supports offline
name search, and `assemble_opts` already stores `seeds_db` / `labels_db` per
option set, so different samples in one project can already use different
databases.

## Goal

Add a reference-database step between taxonomy assignment and the assembly
workflow that:

- proposes one or more clades from the samples' confirmed taxonomy, broad
  enough to have usable reference data but not so broad as to be useless;
- keeps distant groups separate (fish and sea stars in one project get two
  databases, not one animal-wide database);
- lets the user override any proposal with any NCBI taxon;
- builds the databases without blocking the app, and wires the results to the
  right samples automatically;
- can be skipped entirely, and can reuse databases built earlier or by
  someone else;
- stores databases wherever the user wants, including a shared location used
  across projects and users.

Non-goals: changing the Nextflow layer, changing the seed/label file format,
or changing how `build_reference_db()` talks to GenBank.

## Data model

### New project table `reference`

One row per clade group.

| column | meaning |
| --- | --- |
| `ref_id` | group key, e.g. `asteroidea` (also the clade short name used for filenames) |
| `clade_taxid` | chosen NCBI node |
| `clade_name` | scientific name of that node |
| `clade_rank` | rank of that node |
| `query` | GenBank term, `"<clade_name>"[Organism]` |
| `n_samples` | project samples covered |
| `n_records` | GenBank rRNA record count for the clade |
| `state` | 0 hold, 1 ready, 2 done, 3 failed, 4 building, 5 skipped |
| `source` | `built`, `existing`, or `skipped` |
| `seeds_db`, `labels_db` | absolute paths once available |
| `log_file`, `pid` | background build tracking |
| `fail_reason` | one-line summary shown in the table |
| `time_stamp` | last change |

### New project table `reference_config`

Single row, project-level settings.

| column | meaning |
| --- | --- |
| `ref_dir` | folder holding built databases (nullable) |
| `min_records` | suggestion floor, default 1000 |
| `max_records` | suggestion ceiling, default 25000 |

`ref_dir` resolution order: project setting, then the `RIBOPILOT_REF_DIR`
environment variable, then `ribopilot_ref_dir()`. `new_project()` gains a
`ref_dir` argument that seeds this row.

### Change to `samples`

New column `ref_id`, the group a sample belongs to. Null until the step runs.

### Cross-project manifest

`<ref_dir>/manifest.tsv`, appended on every successful build: clade name,
taxid, rank, query, record count, build date, seed file, label file. This is
what makes "you already have a database for Asteroidea" possible in a fresh
project, and what the Use Existing picker lists.

### Unchanged

Seeds and labels continue to reach Nextflow through `assemble_opts`. No
pipeline, module, or samplesheet change.

## Suggestion algorithm

Input is every sample with a confirmed taxid; samples without one are listed
as unresolved rather than grouped.

Sample lineages currently store names only. A helper using
`taxonkit reformat -t` adds the taxid at each rank, giving the
species/genus/family/order/class/phylum/kingdom ladder the picker needs.

For each distinct lineage:

1. Start at family.
2. Query GenBank for that clade's rRNA record count, reusing
   `download_rrna.sh --count-only`.
3. If the count is below `min_records`, step up one rank and repeat.
4. Stop at the first rank meeting the floor; also stop if the next step up
   would exceed `max_records`, or on reaching phylum.

Samples landing on the same node form one group. A group nested inside another
is folded into the parent only if the parent's count is still at or below
`max_records`.

Cost is bounded: only distinct lineages are walked, walks stop early, and
counts are cached by taxid in the manifest. Fish stop at Actinopterygii, sea
stars at Asteroidea, so a mixed project yields two groups without any user
input.

Floor and ceiling are editable on the tab, so a user working on a
record-poor group can lower the floor instead of fighting the suggestion.

## User interface

### New "Reference" mode

Added to the mode picker between Taxonomy and Assemble. Table columns:

`Clade` (link) | `Rank` | `Samples` | `GenBank records` | `Source` | `State` | `Database`

Control bar: Suggest, Edit Clade, Use Existing, Skip, Build, alongside the
existing Refresh and selection controls. The header shows the active reference
folder with an Edit control.

### Clade ladder modal

Opened from the Clade cell. Shows the group's shared lineage as a row of rungs
from species up to phylum. Each rung shows how many of the group's samples it
covers and its GenBank record count; the suggested rung is highlighted. Below
the ladder, an offline taxdump search box accepts any other NCBI taxon.

Two cases get explicit handling:

- A rung covering only part of the group (dropping from Asteroidea to
  Forcipulatida when one of four stars is not a forcipulatid) reports the
  shortfall and offers to split the uncovered samples into their own group.
  This is the primary way a user reshapes groups.
- An off-ladder search result that contains none of the group's samples is
  warned about but allowed, since a sister-group database is sometimes
  deliberate.

### Use Existing

Lists databases found in the reference folder's manifest (clade, record count,
build date) and accepts typed seed and label paths for databases from
elsewhere. Wiring to samples is identical to a successful build; `source` is
recorded as `existing`.

### Skip

Marks selected groups skipped. Their samples keep whatever database their
current option set holds, which for an untouched project is the shipped
Actinopterygii default. Available per group and as Skip All.

### Reference folder

Editable from the tab header: a validated text field (must exist, must be
writable, `~` expanded). Changing it re-scans that folder's manifest so
databases already there appear immediately as reusable. Changing it never
moves or copies existing files, and rows already wired to a database keep
their absolute paths.

### Gating

Assemble is not blocked. If any group is still pending when Update is pressed
on Assemble, the run modal names those samples and states which database they
will actually use.

## Building

Each build spawns a detached `processx` job running `build_reference_db()` in a
fresh R session, logging to `<project>/.reference/<ref_id>.log`, with the pid
and `state = 4` written to the table. Several clades can build concurrently and
the app stays usable.

Polling happens on refresh. Process gone with both FASTAs present means done;
process gone without them means failed. Clicking a building or failed row shows
the log tail.

Entrez Direct and biopython come from the environment running the app, normally
the RiboPilot conda environment. The build modal offers an optional conda
environment field for when it does not, an NCBI API key field defaulting to the
`NCBI_API_KEY` environment variable, and a RefSeq-only toggle.

On success: the FASTAs are placed in the reference folder, the manifest is
appended, an `assemble_opts` set named after the clade is created or updated
with the new seed and label paths (every other setting copied from each
sample's current set), and the group's samples are moved onto it. Users can
still override on the Assemble tab afterwards.

An existing database for the same clade is reused rather than re-downloaded; a
Rebuild action in the build modal forces a fresh download.

### Failure handling

Named failures, each setting `state = 3` with a one-line `fail_reason` and the
full log a click away:

- Entrez Direct or biopython missing from the environment.
- Zero matching GenBank records.
- Download interrupted or network failure.
- Script completed without producing seed or label FASTAs.

Nothing is wired to samples on failure; the previous database assignment
stands.

## Testing

Grouping and ladder construction take the record-count lookup as an injectable
function, so the suite runs offline against stubbed counts:

- fish and sea stars split into two groups;
- a record-poor clade walks further up the ladder;
- the ceiling stops a walk before it over-broadens;
- a nested group merges into its parent only when the parent is under the
  ceiling;
- picking a rung that covers part of a group produces the expected split.

Option-set wiring, folder resolution order, and the manifest round trip are
tested against a temporary SQLite project. No test contacts NCBI.

## Files

New: `R/reference_groups.R` (lineage ladder, record counts, grouping),
`R/reference_build.R` (background build, manifest, option-set wiring),
`R/app_reference.R` (tab table, ladder modal, handlers).

Modified: `R/init_db.R` (two new tables, `samples.ref_id`),
`R/migrate_db.R` (migration for existing projects), `R/init_project.R`
(`ref_dir` argument), `R/app_ui.R` (mode picker and control bar),
`R/app_server.R` (mode wiring), `R/build_reference_db.R` (respect the
project reference folder, expose a count helper).
