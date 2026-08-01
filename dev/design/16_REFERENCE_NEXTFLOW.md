# 16. Reference databases through Nextflow

Status: design, not yet implemented. Supersedes the detached-Rscript build and
the host-side record counting introduced with the Reference tab.

## Why

The Reference tab currently shells out to Entrez Direct on the host for every
record count, and builds each clade database in a detached `Rscript` process.
Both assume the user has `edirect` and biopython installed locally, which is
false for anyone running RiboPilot from the Docker or Singularity image. Builds
also run one at a time, in a process the app has to babysit by pid.

Moving both to Nextflow makes the container the only dependency surface, and
lets several clades count or build in parallel under whatever executor the
project is already configured for (local, slurm, sge, LSF, AWS Batch).

## Scope

In scope: record counting (Suggest and the edit-clade modal) and database
building. Out of scope: how clades are chosen (the ladder walk is unchanged),
the manifest format, preset wiring, and everything downstream of a built
database.

## Workflow

New `inst/nextflow/main_reference.nf`, run with the project's own `.config` so
it inherits the configured executor and container settings.

Two entry points:

- `-entry COUNT` - one task per clade. Runs
  `ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh --count-only` and
  writes `<taxid>.count`, a single line of `taxid<TAB>n_records`, into the run
  output directory.
- `-entry BUILD` - one task per group. Runs the same script in download mode
  plus `parse_rrna_fasta.py`, publishing `<ref_id>_rrna_seeds.fasta` and
  `<ref_id>_rrna_labels.fasta` into the reference folder's `seeds/` and
  `labels/` subdirectories, matching `reference_db_paths()`.

Inputs:

| Param | Meaning |
| --- | --- |
| `--clades` | Path to the work-list TSV the app writes before launching |
| `--refDir` | Reference folder (published output, and the counts cache) |
| `--outDir` | Run output directory for `.count` files |
| `--ncbi_api_key` | Optional, raises the NCBI rate limit |
| `--refseq` | Optional flag, restricts the query to RefSeq |

The work-list TSV has one row per clade with columns `ref_id`, `taxid`,
`clade_name`, `query`. COUNT reads `taxid` and `query`; BUILD reads `ref_id`
and `query`. One format serves both so the writer has one code path.

Every `download_rrna.sh` invocation runs with stdin detached. Entrez Direct
reads stdin whenever it is open and will block forever on an inherited socket.

## Mounting the reference folder

The reference folder lives outside the project (default
`~/.local/share/R/RiboPilot/ref_dbs`), so container tasks cannot see it by
default. The launcher adds the mount for that run only:

- Docker: `-process.containerOptions "-v <refDir>:<refDir>"`
- Singularity: `-process.containerOptions "--bind <refDir>"`

Nothing is written into the project `.config`. A user who changes the reference
folder gets the new mount on the next run with no project edits, and existing
projects need no migration.

## R side

New `R/reference_nf.R` owns launching and reading results. `R/reference_build.R`
keeps the group state machine unchanged.

    write_clade_manifest(path, clades)     # the work-list TSV
    nf_reference_command(project, entry, manifest, ref_dir, out_dir, ...)
    run_clade_counts(con, clades)          # blocking: launch, wait, merge
    read_count_files(dir)                  # <taxid>.count -> data.frame
    merge_counts(ref_dir, counts)          # into counts.tsv

### Suggest

1. Build the ladders (unchanged).
2. Collect every candidate rung: all ranks at or above the walk's start rank,
   across all samples, deduplicated by taxid. Typically 20 to 40 clades.
3. Drop any taxid already in the reference folder's `counts.tsv`.
4. If anything remains, run COUNT and wait, behind the existing
   "Suggesting..." alert. Merge the resulting `.count` files into `counts.tsv`.
5. Walk the ladder offline against the cache. The walk makes no network calls.

The count function passed to `suggest_reference_groups()` becomes cache-only:
a taxid with no cached count returns `NA`, which the walk already skips.

### Edit clade

Applying a clade picked from the taxdump search needs a count the batch run
never produced. It runs COUNT with a single-row work list and waits, behind the
"Applying clade..." alert. One counting code path, no host `edirect`.

### Build

`start_reference_build()` writes one work-list TSV covering every selected
group and launches a single detached BUILD run. Each group row stores that
run's pid and log path, so several rows share both.

`poll_reference_builds()` is unchanged in shape: when the pid is gone, each
group is judged on its own by whether its seed and label FASTAs exist and are
newer than the launch timestamp, then wired or failed individually. One clade
failing inside a run does not fail the others.

The existing reuse shortcut stays: a group whose database is already on disk is
wired without launching anything.

## Errors

| Situation | Behavior |
| --- | --- |
| Nextflow not installed | "Nextflow not found" before anything is launched |
| COUNT run exits non-zero | Merge whatever `.count` files exist; missing counts are `NA` and the walk skips them |
| No counts at all | Suggest fails with the existing "could not count records" message plus the log tail |
| BUILD run exits non-zero | Each group with no fresh FASTAs gets `fail_reason` from the log tail; groups that did produce FASTAs are wired normally |
| Reference folder not writable | Caught before launch by the existing folder check |

## Testing

R-side tests use the existing stub-launcher pattern, so nothing is spawned:

- The work-list TSV has the expected columns and one row per clade.
- The command vector carries the right entry point, params, and mount option
  for Docker and for Singularity.
- `.count` files merge into `counts.tsv` without duplicating existing taxids.
- Suggest walks correctly from a pre-seeded cache with no launcher call, and
  requests counts only for uncached taxids.
- A simulated finished BUILD run wires the groups whose FASTAs are fresh and
  fails the rest, with a shared pid across rows.

The workflow file itself is checked by the existing docker-image smoke test
rather than by unit tests.

## Migration

`build_reference_db()` stays as the CLI entry point for users who want to build
a database outside the app; it keeps its current host-tool behavior. The app no
longer calls it. `genbank_rrna_count()` becomes unused by the app and is kept
only for the CLI path.
