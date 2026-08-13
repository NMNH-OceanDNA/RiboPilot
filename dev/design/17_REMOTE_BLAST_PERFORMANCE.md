# 17. Remote BLAST and reference fetch performance

Status: design, not implemented. No code written.

## Why

Remote BLAST is 94 percent of pipeline CPU. Measured from
`full_project_docker/.nextflow.log` (10 fish samples, 26 active units, submit to
complete, all tasks):

| process | total CPU | median per task | max |
|---|---|---|---|
| blast_genbank | 941.8 min | 7589 s (2.1 h) | 13848 s (3.8 h) |
| assemble | 31.6 min | 190 s | 378 s |
| blast_ref_fetch | 6.9 min | 40.5 s | 116 s |
| blast_ref_align | 1.3 min | 2.3 s | 6.5 s |
| annotate_tools | 1.8 min | 5.0 s | 6.3 s |

One remote megablast against `core_nt` costs roughly two hours of queue time.
That project spent about 17 hours in BLAST at `maxForks 3`. At the 20 to 100
sample scale the pipeline targets, that is days.

`blast_ref_fetch` is 0.7 percent of the run. Optimising it is worth doing for
correctness, not for speed, and this document says so wherever it matters.

Alongside the cost problem there is a silent-yield problem. In the same project,
19 of 26 units end with no reference alignment at all, from two independent
failures that report success:

| stage | units | lost |
|---|---|---|
| active units | 26 | |
| got BLAST hits | 26 | |
| survived the qcovs floor | 16 | 10 |
| got a reference alignment | 7 | 9 |

The first 10 are `BUG_RiboPilot-main_qcovs-floor-rejects-complete-assemblies.MD`.
The other 9 are units whose rank-1 candidate carried no parsable GenBank feature
table, with no fallback to rank 2. Curation then runs without `ref_pid` and
`aln_cov` and does not say so.

## Scope

In scope: how `blast_genbank` and `blast_ref_fetch` search, what they store, how
reference coordinates are derived, and where reference data lives.

Out of scope: assembly, scoring, path selection, the annotation of the sample
itself, and the triage ruleset. The contamination screen's *inputs* are in scope
because narrowing the search narrows them; its *rules* are not.

## Measurements this design rests on

All figures below are measurements on real data, not estimates.

**Local annotation of GenBank reference records.** Running the pipeline's own
`pybarrnap -k euk` and `ITSx -t all --preserve T` over the 517 reference records
already fetched into the three test projects, then feeding the output through
RiboPilot's real `parse_barrnap()`, `parse_itsx()`, and `reconcile_features()`:

| | GenBank feature table | local barrnap + ITSx |
|---|---|---|
| records with >=1 usable element | 215 of 517 (41.6%) | 517 of 517 (100%) |
| records with all 5 elements | 44 (8.5%) | 330 (63.8%) |
| 18S | 44 | 330 |
| ITS1 / 5.8S / ITS2 | 44 each | 346 each |
| 28S | 215 | 517 |

Strict superset: every one of the 371 accession/element pairs already in
`blast_ref_annotations` also got a local call. None were lost.

Local annotation is at ceiling, not merely better. Records that did not get five
elements do not contain five: they are RefSeq `XR_*` 28S-only records around
4 kb, or ITS-only submissions truncated inside ITS1, where ITSx correctly
reports "SSU: Not found".

**Agreement.** Across all 371 shared pairs, no start or stop differs by more
than 20 bp, and the single largest disagreement anywhere is 10 bp. Offsets are
systematic and traceable to Rfam model boundaries plus the existing spacer
clipping: 18S stop -5, 28S start +1 to +9, 28S stop -2, ITS1 start -4, ITS2
stop +5. On a 1840 bp 18S that is 0.3 percent.

**Effect on curation.** Recomputing `compute_ref_alignment()` with local
coordinates and re-running the real `curate_rrna_core()` over 21 rank-1 units
produced zero flag changes. `pid` is provably insensitive: jittering reference
boundaries by 20 bp at both ends moved `pid` by at most 0.55 points, and the
measured worst offset is 10 bp. `aln_cov` moves, and it moves up (ITS1 0.990 to
1.000, ITS2 0.972 to 1.000).

**Annotation cost.** 232 records in 156 s at 8 cores, 130 s at 16. ITSx is 72 to
81 percent of it and scales sublinearly (6.0x at 16 threads). pybarrnap barely
scales past 4 threads. `ITSx -t M` (metazoa profiles) runs 1.6x faster than
`-t all` and is the right default for this package's users.

## Design principle: separate searching from interpreting

`parse_blast_candidates()` currently writes `blast_hits`, `blast_ref_candidates`,
and `blast_screen` in one pass over a hits file, applying the candidate filter at
parse time. Because the filter runs during parsing, changing any threshold costs
a full remote re-search, at two hours per unit.

`blast_hits` already stores everything the filter needs: `pident`, `qcovs`,
`bitscore`, `aln_len`, `qstart`/`qend`, `sstart`/`send`, `hit_class`, `staxid`,
`lineage`, `lineage_taxids`, `rank_distance`.

So candidate selection and the screen become pure, re-runnable functions over
`blast_hits` rather than over the hits file. Consequences:

- The qcovs fix applies retroactively to existing projects at zero NCBI cost.
- Every future threshold change costs seconds, not hours.
- Filter changes become safe to iterate on and could be exposed in the app.

This principle is why Stage 1 comes first: it makes everything after it cheap to
get wrong.

## Stage 1: correctness, no NCBI traffic

**Split the filter out of the parser.** `parse_blast_candidates()` keeps writing
`blast_hits`. A new function derives `blast_ref_candidates` and `blast_screen`
from the stored table, and can be re-run on demand.

**Fix the query-coverage floor.** Adopt Option A from the bug report: accept a
subject on `qcovs >= REF_CANDIDATE_MIN_QCOVS` OR on an absolute aligned-length
minimum. `qcovs` for a complete rDNA assembly measures how much intergenic
spacer the assembler recovered, inverted, so the floor currently rejects the best
assemblies. All 10 affected units in `full_project_docker` have a best HSP of
7208 bp at 100.0 percent identity.

The absolute floor needs sizing against a clade with shorter elements before it
is fixed in code. 18S is about 1.8 kb and 28S about 4 kb, so a value in the 1500
to 2000 range admits a reference covering at least one full subunit, but the fish
numbers must not be assumed to generalise to Asteroidea, Scyphozoa, or
Octocorallia.

**Capture the columns a re-search would otherwise cost.** Each `-outfmt` change
costs one full remote re-run of every project, so add every column that might be
wanted, once, at `inst/nextflow/modules/blast.nf:63`:

- `slen`, which enables the principled subject-coverage filter (Option B) later
  with no second re-run.
- `qseqid`, cheap insurance if per-sample query batching is ever revisited.

`parse_blast_hits_file()` already reads trailing columns through a guarded
accessor, so an older hits file returns NA rather than breaking.

**Mark units that were searched and found nothing.** The resume guard at
`inst/nextflow/modules/blast.nf:51` counts rows in `blast_hits`, so a unit with
no hits is re-searched on every retry and every re-run, forever. Store the
searched state explicitly, keyed on a hash of the query sequence so that a
re-assembled or trimmed unit is correctly re-searched rather than inheriting the
previous contig's hits.

The content hash fixes a second bug in the same change: nothing currently clears
`blast_hits` when an assembly changes, so `blast_screen` and
`blast_ref_candidates` can keep describing a contig that no longer exists.

Note the cost: a re-run will now do real BLAST work where today it wrongly does
none. That is correct, and users will notice.

**Stop swallowing efetch failures.** `inst/nextflow/modules/blast_ref.nf:31-32`
ends both efetch calls with `|| true`, and the error-body check misses NCBI's
letter-spaced `Error: F a i l e d  t o  u n d e r s t a n d  i d` form. In
`e2e_ponytail`, 187 accessions were wanted and 50 fetched; nothing recorded the
other 137, and the run reported success. Reconcile requested against returned,
log the difference once, and stop retrying accessions that are permanently gone.

**Stop leaking the API key.** `export NCBI_API_KEY='${params.ncbi_api_key}'` at
`blast.nf:41` and `blast_ref.nf:23` is Groovy interpolation, written verbatim
into `.command.sh` and `.command.run` in every task work directory. Work
directories on shared HPC filesystems are typically group readable. Use Nextflow
secrets or an `env{}` scope populated from the launching shell. Separately, the
key is dead code in `blast.nf`: it is an E-utilities credential and has no effect
on the BLAST URL API.

## Stage 2: reference yield

**Fetch sequence only.** Nothing consumes what the feature table provides once
annotation is local. `blast_ref_sequences.taxid` and `.lineage` are already
written as NA and have no reader anywhere in the package. `species` comes from
the FASTA defline via `extract_species()`. `hit_class` is computed from the
blastn `stitle` before any fetch. The only real loss is coordinates, which local
annotation replaces as a strict superset.

Dropping `-format ft` halves the request count and makes batching unambiguous: a
concatenated multi-FASTA splits on `>`, whereas a concatenated stream of feature
tables does not.

**Batch the fetch and call R once.** Replace the per-accession loop with a
chunked `efetch -db nuccore -id "a,b,c,..."`, then one R invocation per sample
instead of one per accession. R interpreter startup is 0.87 s with the installed
package and 2.52 s from a source checkout, paid once per accession today, which
is 38 to 65 percent of the step. This is a real fix but it is a fix to 0.7
percent of the run; it must not be sold as a speed win for the pipeline.

Map records back to accessions by the accession parsed from each record, never by
position: efetch silently omits unresolvable ids and EDirect re-sorts. On a
non-empty set difference between requested and returned, retry the missing
individually and record per-accession fetch status. Do not use the
`epost -format acc | efetch` recipe: verified against EDirect 25.5, it returns
unrelated records with exit 0, which would store the wrong sequence under the
right accession.

**Annotate references locally**, with the same `pybarrnap` and `ITSx` path used
for the samples, so the sample-versus-reference comparison is finally
like-for-like.

### Three blockers that must be handled

Each of these fails silently and produces plausible, wrong output.

**Multi-FASTA collapses to one record.** `parse_itsx()` at `R/annotate.R:62` is
`fields <- strsplit(pos_lines[1], "\t")[[1]]`, reading only the first line of the
positions file. `parse_barrnap()` at `R/annotate.R:45-52` returns
`data.frame(element, start, stop, strand, model_cov)`, having parsed the GFF
seqid into a local and never read it, so rows from different records are
indistinguishable. `reconcile_features()`'s `pick()` at `R/annotate.R:98` then
takes `[1, ]`. Verified: real 235-record tool output through the real functions
returns one 5-element model, no error, no warning. `annotate_seq()` documents the
assumption at `R/annotate.R:336` ("a single-sequence fasta").

Either demultiplex tool output by sequence id before calling these, or carry an
id column through both parsers and group inside `reconcile_features()`.

**One record over 100 kb aborts the whole ITSx batch.** ITSx bundles HMMER's
protein pipeline, which hard-fails at `Target sequence length > 100K, over
comparison pipeline limit`. Tested with a 1.4 Mb record: ITSx printed a fatal
error and wrote no positions file at all. Tested with 20 normal records plus one
large one: same abort, all 20 lost. `run_itsx()` at `R/annotate.R:265-271` passes
`stdout=FALSE stderr=FALSE` and returns an empty data frame when the file is
absent, so this is completely silent.

The gate must be per tool, not per accession: skip ITSx above 100 kb, still run
pybarrnap (which handled 1.4 Mb in 44 s), still store the sequence. Dropping the
accession outright would be a regression, since today its GenBank feature table
is parsed. And these records are not rare by accident:
`classify_hit_title()` assigns `uninformative` before `rdna`, so rDNA phrasing
wins, and "chromosome 5, 18S ribosomal RNA gene region" classifies as `rdna`.
Surface the ITSx exit status the way `run_rrna_tool()` already does for barrnap.

**Minus-strand records mix coordinate frames.** barrnap reports the input frame;
ITSx with `--preserve T` silently searches the other strand and reports in the
flipped frame. Verified on reverse-complemented records: for KU240018.1, barrnap
gave 18S at 4883-6722 strand `-` correctly, while ITSx reported the byte-identical
line as the forward run. Reconciling that gives a 5.8S about 1900 bp from its
true position and drops ITS1 and ITS2 entirely via `clip_spacer()`.

The fix is cheap and needs no schema change: take the strand from the batch GFF
already in hand, reverse-complement only the minus-strand records into a second
FASTA, run ITSx on that, and map back with `start' = L - stop + 1`. References do
not need reorienting for alignment, because `compute_ref_alignment()` at
`R/blast_ref_align.R:105-107` and `compute_feature_msa()` at `:30` already
reverse-complement a feature whose stored strand is `-`.

### One guard to add while here

`compute_ref_alignment()` writes a row for every shared element with no floor on
coverage or length. Because `pairwise_nt()` is ends-free, two hypervariable
spacers from different families can collapse to a near-zero-length alignment.
Recomputing with local annotations on `full_project_docker` produced 24 new rows,
of which 8 were degenerate (`pid` NaN, which lands in SQLite as NULL, or
`aln_cov` 0.002). Curation survives, because ITS1 and ITS2 are not in
`RRNA_GENES` and the `min_ref_coverage` gate would reject them, but the rows are
written regardless.

The codebase already learned this lesson: `R/params_rrna.R:61-67` documents a
1277 bp unit scoring 100.00 percent against a 6310 bp one over 0.2 percent of it,
and pairs `min_isomer_anchor_pid` with `min_isomer_anchor_coverage` and
`min_isomer_anchor_bp`. `compute_ref_alignment()` needs the same pair: refuse to
write a row below a coverage and length floor.

### Drop the stored alignment strings

`blast_ref_alignment.sample_aln` and `.ref_aln` are write-only. Verified: the
sole writer is `R/blast_ref_align.R:117` and the sole schema is
`R/init_db.R:460-461`. The only reader of the table at all is
`R/curate_rrna.R:316`, which names `element, pid, aln_cov`. No `SELECT *` runs
against the table, no test references it, and the export report's alignment
section builds sample-versus-sample alignments from `annotations` and
`assemblies` with DECIPHER on the fly (`R/export_report.R:31,70`) rather than
reading this table. The `sample_aln` at `R/blast_ref_align.R:56` is a parameter
name of `alignment_coverage()`, not the column.

Drop both columns and the `r_aln` assignment at `R/blast_ref_align.R:112`, which
feeds storage only. Keep `s_aln` at line 111: both surviving metrics derive from
it, `aln_cov` via `alignment_coverage(s_aln, nchar(s_sub))` and `aln_len` via
`nchar(s_aln)`.

Unusually for this plan, this change does not break existing projects. Both
columns are nullable `TEXT` and neither is in the primary key, so an old table
accepts the narrower `dbAppendTable` frame and writes NULLs.

The motivation is that these are unread columns which invite future misuse, not
disk. Measured storage is 94 KB over 35 rows in `full_project_docker` and 36 KB
over 14 in `e2e_run2`; even at Stage 2's five-fold row growth that stays under a
megabyte.

## Stage 3: the shared reference store

Reference records are project-independent. A second project in the same clade
currently re-fetches and re-annotates everything.

The store holds sequence, locally derived element coordinates, and provenance
(tool versions, kingdom, ITSx taxon, build date), keyed by accession. Projects
copy rows in rather than referencing across files, so a project directory stays
copyable and self-contained.

**Location.** Under the folder resolved by `ribopilot_ref_dir()`
(`R/build_reference_db.R:13`), which already resolves the per-project
`reference_config.ref_dir`, then `RIBOPILOT_REF_DIR`, then
`tools::R_user_dir()`. The app already exposes that folder with a browse and
validate control (`R/app_reference.R:142`, `validate_ref_dir()`). Users on
clusters with home-directory quotas point it at scratch or project space once,
and the clade databases, the reference store, and the BLAST database of Stage 4
all land there together.

Two things this needs that do not exist yet: the path must be passed to the
Nextflow tasks as a param and bind-mounted under Docker in
`R/generate_config.R`, and writes must be safe for a shared filesystem with
concurrent users, unlike the read-modify-write in `append_ref_manifest()`.

**Read-only stores are currently broken.** `open_project_db()`
(`R/pipeline_steps.R:167-171`) runs `PRAGMA journal_mode=wal` at connect. Tested
against a mode 444 file in a mode 555 directory, that fails at connect, so every
read fails too and the wrapped `tryCatch` swallows it. A group-shared store where
one admin populates it and everyone else reads is exactly the recommended HPC
deployment, so the store needs a read-only open path that skips the pragma.

**Cross-project only.** Within a single project, accession overlap between
samples is negligible: 235 distinct across 236 per-sample in `e2e_run2`, and 232
distinct across 243 per-sample (4.7 percent) in `full_project_docker`. The
intra-run dedup win is about 1.00x. This stage is justified by cross-project
reuse and by fixing the check-then-act race on `blast_ref_sequences`, not by
throughput.

## Stage 4: the shipped local rDNA BLAST database

This section is pending the build-and-ship design round and will be completed
before this document is reviewed.

## Migration

Every stage adds columns, and this project deliberately ships no migrations
before release. `R/db_schema.R:147` is the only column-existence check in the
codebase. So each stage breaks projects created before it, and staging is about
sequencing the maintainer's own re-runs, not about backward compatibility.

## Validation

Stages 1 through 3 can be validated on the existing test projects, because their
claims are about yield and about not changing curation, and both are measurable
there.

Stage 4 cannot. Every sample in all three test projects is a ray-finned fish
(`samples.lineage` is `Metazoa;Chordata;Actinopteri;...` for all 10 in
`full_project_docker`), and every pre-existing reference alignment is `pid` 100.0,
so the invertebrate identity floors in `RRNA_ANCHOR_PID_INVERT` are never
reached. Running an invertebrate validation on these projects returns "no change"
and proves nothing. Sensitivity of a local database for a divergent invertebrate
must be tested on the Asteroidea, Scyphozoa, and Octocorallia projects.

One knife-edge worth knowing about: `SRR33849683.1.1`'s 28S scores raw `pid`
84.96, stored as `round(pid, 1)` = 85.0, against a fish floor of 85 tested with
strict `<`. It does not flag, by 0.04 points, because of rounding. Expect a
handful of units parked on a floor like this, ready to flip on any retune.
