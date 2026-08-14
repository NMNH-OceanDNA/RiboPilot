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
directories on shared HPC filesystems are typically group readable. Separately,
the key is dead code in `blast.nf`: it is an E-utilities credential and has no
effect on the BLAST URL API.

Moving the key to an `env{}` scope only gets it out of `.command.sh`. Nextflow
still writes every `env{}` value verbatim into `.command.run` for each task, so
the credential still lands in the same group-readable work directory, just in a
different file. Keeping it out of both files needs Nextflow secrets, or a key
file passed by path and read at runtime from outside the work dir. Deferred to
its own task.

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

**Reference annotation is pinned to kingdom euk and taxon M, not a per-project
setting.** The shipped reference database is metazoan, so there is no per-run
choice to make: barrnap always runs with kingdom "euk" and ITSx always runs
with taxon "M". Taxon M was measured byte-identical to taxon "all" across 82
reference records sampled from the 891k-record Metazoa harvest, spanning
Aedes, Artemia, Aricidea, Donax, Dugesia, Echiniscus, Calopteryx and other
metazoan phyla: 96 rows both ways, 96 identical coordinate pairs, zero
elements lost, zero gained, and about 1.4x faster (12 s versus 17 s). Kingdom
is pinned rather than left to a per-project option because a wrong kingdom
does not merely shift boundaries: on four fish references, barrnap under
kingdom "bac" reports 16S_rRNA/23S_rRNA, which parse_barrnap maps to nothing,
so the anchors silently fall through to ITSx's own SSU/LSU calls and the 28S
extent moves by up to 224 bp. The store's primary key is accession alone;
kingdom and taxon are still recorded on each row as provenance.

## Stage 4: the shipped local rDNA BLAST database

The build is: harvest, byte-exact dereplication, taxid map, `makeblastdb`,
manifest. No annotation pass, no species collapse, no misidentification screen.
The reasoning for each omission is below.

**The escalation path in "When the local panel cannot answer" is a hard
prerequisite for shipping any of this.** Without it, moving the search local
reinstates precisely the restriction that `inst/nextflow/modules/blast.nf:7-13`
was written to remove.

### Scope of the harvest

The GenBank query reuses the title and organelle terms from
`inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh:67,110`: the clade
term, AND the nuclear rRNA title terms (18S, 28S, 5.8S, small subunit, large
subunit, internal transcribed spacer), NOT `mitochondrion[Filter]`, NOT
`plastid[Filter]`.

**It drops that script's `efilter -molecule genomic`.** Measured for
`txid33208[Organism:exp]` (Metazoa):

| query | records |
|---|---|
| title terms, no organelles | 891,382 |
| the same, `-molecule genomic` | 740,941 |
| RefSeq subset of the unfiltered query | 148,209 |

RefSeq `XR_*` and `NR_*` rRNA records carry `MolType=rna` (verified:
`XR_003928005.1`, Biomol rRNA, MolType rna, 3947 bp), so the molecule filter
removes 150,441 records. The intersection of the RefSeq subset with
`-molecule genomic` is **zero**: the filter removes the entire RefSeq rRNA set,
not merely most of it.

Those are exactly the records that dominate real candidate lists. In the
reference sequences the test projects actually fetched from core_nt, RefSeq
`XR_*` 28S-only records around 3.9 to 4.0 kb accounted for 90 of 235 accessions
in `e2e_run2`, 49 in `full_project_docker`, and 32 in `e2e_ponytail`, and they
are why 28S has 215 usable references where every other element has 44. Keeping
the filter would make the local panel systematically worse than the remote search
it replaces, and worse for 28S specifically.

The filter is correct where it currently lives, because `download_rrna.sh`
builds GetOrganelle seeds and labels, which do want genomic DNA. It is wrong for
a BLAST panel. This is a difference between the two builders, not a bug in
either.

891,382 records is roughly 528 MB of FASTA at the 593 bytes per record measured on
the Cnidaria subset.

### The database ships unannotated. Annotation is on demand.

Element coordinates for a reference are derived when a run actually needs that
reference, not when the database is built.

The case is arithmetic. A build-time pass covers all 891,382 records. A project
needs only its own distinct candidate accessions, measured at 235, 232, and 187
in the three test projects. On-demand annotation therefore does about one three
thousandth of the work. At the measured rate (232 records in 156 s at 8 cores) a
build-time pass would take 21 to 34 hours even with ITSx sharded single-threaded
across 16 cores and `-t M` in place of `-t all`, while a project pays a couple of
minutes, parallelised across the per-unit fan-out that already exists, and
Stage 3's store caches the result across projects. The store then warms toward
exactly the records users touch.

This also means Stage 4 adds no annotation machinery at all. Stage 2 already
annotates fetched references at run time; Stage 4 only changes where the sequence
comes from, `blastdbcmd` instead of `efetch`. The ITSx 100 kb ceiling and the
minus-strand frame mismatch are handled once, on the Stage 2 path, rather than in
both places.

### Deduplication is exact-match only

At roughly 528 MB raw the database is small enough that curation is not worth its
risk. The build removes byte-identical sequences and nothing else. No clustering
at a similarity threshold, no collapse to one record per species, no
misidentification screen.

What this avoids is worth stating, because each was a real hazard:

- **Congener collapse.** Many congeneric metazoan 18S sequences are over 99
  percent identical, so clustering at any threshold tight enough to reduce size
  meaningfully would merge species and destroy the species label the app
  displays.
- **Intraspecific ITS variation.** ITS1 and ITS2 vary within a species, and
  between paralogous copies within one genome. A single retained representative
  can be a worse match than the conspecific it displaced.
- **Coverage holes.** Dropping a flagged record that is the only one for its
  taxon removes that taxon from the panel entirely.
- **Circularity.** Screening species labels for error requires trusting species
  labels to group the records being screened.

The aggressive collapse considered earlier was motivated by making a build-time
annotation pass tractable. With annotation moved to run time that motive is gone,
and with it the reason to accept any of the above.

Contamination and misidentification are still caught, just at the point where
there is evidence to catch them: a misidentified reference that becomes a unit's
rank-1 candidate shows up in the contamination screen's lineage comparison and in
triage's odd-one-out check, both of which already exist and both of which see the
sample's own taxonomy. ITSx's chimera flag still fires during run-time
annotation.

### What the build still has to produce

A taxid map. `staxids` is not optional: `parse_blast_candidates()` resolves every
distinct staxid through `tk_ranked_lineage()` to populate `lineage`,
`lineage_taxids`, and `rank_distance`, and the contamination screen's
bitscore-weighted LCA vote and triage's lineage comparison both read them. A
database built without `-taxid_map` would silently strip the taxonomy channel out
of the screen. The map comes from a batched `efetch -format docsum` pass over the
accession list, which is a lookup rather than a curation step.

Key it on `AccessionVersion`, never `Caption`. `Caption` is the unversioned
accession while efetch FASTA deflines carry `ACC.V`, so a map keyed on `Caption`
matches nothing, `staxids` comes back empty, and the failure is silent.

The build must also reconcile the accession set three ways: the `esearch` count,
the records in the FASTA, and the rows in the docsum table. efetch silently omits
unresolvable ids and EDirect re-sorts, which is how `e2e_ponytail` wanted 187
accessions, got 50, and reported success. Any accession in one set and not
another is retried individually and then written to a dropped-records ledger with
a reason. A build that cannot reconcile does not ship.

### When the local panel cannot answer

This is the release blocker, and it is why the local database cannot ship alone.

`inst/nextflow/modules/blast.nf:7-13` records that the search was deliberately
made unrestricted, because restricting it "also restricted what a contaminant
could look like: a symbiont, bacterial 16S, or mitochondrial rRNA assembly could
only ever match some distant nuclear rDNA record and so looked plausible." A
metazoan nuclear-rDNA panel reinstates that restriction in three places:

1. `consensus_hit_class()` (`R/hit_class.R:58`) votes `hit_class` over the top
   hits to set `blast_screen.is_nuclear_rrna`. Against an rDNA-only panel every
   title is rDNA, so the vote returns 1 unconditionally and the flag stops
   carrying information.
2. `screen_class_for()` (`R/screen.R:121`) derives an off-target class from the
   hits' LCA. `SCREEN_CLASSES` is bacterial, archaeal, fungal, algal_plant,
   other_euk, metazoan_other, human. A Metazoa-only panel makes five of seven
   unreachable.
3. A real contaminant returns no hits, which lands in `SCREEN_NO_REF_CODES` as
   `no_hits` and reports `no_reference` rather than `off_target`. Read as
   unchecked rather than clean, which is the safe reading, but the signal is
   gone.

**Superseded 2026-08-14: automatic escalation is shelved.** The maintainer's
decision, and the reasoning is sound: remote search is the thing this whole
design exists to remove, and a pipeline that silently falls back to it can spend
hours without the user choosing to. Instead a unit the local panel cannot answer
is recorded and surfaced, with a suggestion that the user run a remote search
themselves if they want one. The triggers below are retained, not as an
automatic fallback but as the definition of "the panel could not answer this",
which is what gets reported.

Two consequences follow, and the second inverts what this document said before.

**The safe failure is preserved.** Point 3 above is the common contaminant case:
a bacterial or fungal unit returns no hits against a metazoan rDNA panel, lands
in `SCREEN_NO_REF_CODES`, and reports `no_reference`. That reads as unchecked
rather than clean, which is exactly the signal to surface. The pipeline never
calls a contaminant clean; it says it does not know, and now says so to the user.

**Metazoan-versus-metazoan discrimination also survives**, because `staxids`
populate from the local database and the screen's lineage comparison and
`rank_distance` still work. A sea star unit that hits a polychaete reference is
still off-target by lineage.

**What is lost is the specific non-metazoan classification.** Points 1 and 2
stand: `is_nuclear_rrna` becomes constant against an rDNA-only panel, and five of
seven `SCREEN_CLASSES` become unreachable. The app can no longer say "this is
bacterial"; it can only say the local search could not identify it.

**Therefore the decoy panel is promoted from optional to the thing that restores
that signal.** This document previously argued the opposite, on the grounds that
escalation caught what decoys catch. With escalation shelved that argument no
longer holds, and decoys are the only remaining way to distinguish "bacterial
contaminant" from "local panel had nothing".

The triggers, retained as the definition of an unanswerable unit:

| | condition | reading |
|---|---|---|
| E1 | zero local hits | panel cannot answer |
| E2 | `nrow(blast_ref_candidates) == 0` | panel cannot answer |
| E3 | `blast_screen.is_nuclear_rrna == 0` | unit is not nuclear rDNA |
| E4 | `blast_pident < min_anchor_pid["18S"]` | nearest local reference is worse than the identity at which curation would flag the sample's own 18S |
| E5 | `blast_screen.rank_distance >= 5` | best local answer shares nothing below phylum |

Every input already exists. E4 introduces no new constant, reusing the per-clade
floor from `R/params_rrna.R:48-50` (93 invertebrate, 95 metazoan, 97 fish). E5
reuses `rank_distance`, whose docstring rightly warns it must never be read as
contamination on its own; here it is read as "the panel has no near relative",
which is what it measures. Local `dc-megablast` is tried first on the existing
`BESTCOV < 25` trigger, because locally it is free.

E5 will fire often for Asteroidea, Scyphozoa, and Octocorallia. That is correct
behaviour and it is also the honest cost accounting: a clade with thin GenBank
coverage saves less than a fish project does.

Shape: two processes rather than one. `blast_local` with `cpus N`, no network and
no `maxForks` cap; `blast_remote_escalate` with the existing remote body
unchanged at `maxForks 3`. Control via a `blast_opts.remote_fallback` column
(`auto` / `never` / `always`, default `auto`) beside the existing `run_blast` and
`entrez_query`, plus a per-project escalation cap so a badly matched clade cannot
silently become a full remote run.

Two consequences to state plainly. `-entrez_query` is remote-only
(`inst/tool_help/blastn.txt:289`) and is inert on the local path, so log it once
per run and re-apply it on escalation rather than dropping it silently. And
e-values shift by roughly 3,000x between a ~10^12 base `core_nt` and a ~3x10^8
base panel; everything RiboPilot ranks on is bitscore, which is normalised and
unaffected, but a user-supplied `-evalue` floor in `blast_opts.extra_opts` is not,
so warn when one is present and record the effective search space in the manifest.

### The contamination decoys come from NCBI, not from us

The screen's missing classes are restored by searching a decoy panel alongside
our database. `blastn` accepts a space-separated database list in one `-db`
argument, so the panel is composed at query time and never merged into the build.
NCBI already ships the pieces prebuilt, taxid-mapped, and curated:
`16S_ribosomal_RNA` (~72 MB), `SSU_eukaryote_rRNA` (~60 MB),
`LSU_eukaryote_rRNA` (~59 MB), fetched with `update_blastdb.pl`.

These sets are useless as a reference panel, which was established earlier:
frozen at 2022-12-05, and on a real *Acanthaster* assembly the best hit was an
amphibian at 88 percent identity. As decoys the identifier never has to resolve
to a usable reference, only to a kingdom, so that objection does not apply.

**Decided 2026-08-14: the decoys ship.** With automatic escalation shelved they
are the only remaining way to tell a bacterial contaminant from "the local panel
had nothing", so the earlier argument that they were merely a speed optimisation
no longer holds.

Measured sizes and freshness, taken from the FTP headers rather than from an
earlier report that gave a 2022 date:

| database | size | last modified |
|---|---|---|
| `16S_ribosomal_RNA` | 68.6 MB | 2026-08-11 |
| `SSU_eukaryote_rRNA` | 57.0 MB | 2024-04-26 |
| `LSU_eukaryote_rRNA` | 56.6 MB | 2024-04-26 |
| `taxdb` (shared, needed once for staxids) | 61.9 MB | 2026-08-14 |

The 16S set is actively maintained, which is what matters most: bacterial
contamination is the common case. The two eukaryote sets are two years stale,
which is fine for a decoy whose identifier only has to resolve to a kingdom.

Verified working end to end against the built panel: a multi-database search
returns hits from both with `staxids` resolving correctly on each side (81638
*Brosme brosme* from the panel, 562 *Escherichia coli* from a decoy), and a
bacterial query is correctly claimed by the decoy rather than forced onto a
distant metazoan rDNA record.

**Build invariant, discovered by testing rather than assumed:** every database in
a multi-database list must be built with matching `-parse_seqids`. Combining one
that carries accession info with one that does not fails outright with
`BLAST Database error: DB contains no accession info`. NCBI's prebuilt sets and
our panel both carry it, so the composition holds, but any future database added
to the list must too.

### What ships

Layout under `ribopilot_ref_dir()`, beside the clade databases and the Stage 3
store:

```
<ref_dir>/blastdb/<db_version>/
  manifest.json
  ribopilot_rrna.{ndb,nhr,nin,njs,nog,nos,not,nsq,ntf,nto}
  ribopilot_rrna.taxid.map
  ribopilot_rrna.dropped.tsv
```

Size is now measured, not projected. The v1 build (`2026.08.1`, built
2026-08-13) is **342 MB unpacked and 78 MB compressed**, comfortably inside a
GitHub release asset and well under the earlier 120 to 190 MB compressed
estimate. The compression ratio is better than a projection from 2-bit-packed
`.nsq` would suggest, because rDNA is conserved and the retained records are
highly similar to one another.

Measured build, end to end:

| | |
|---|---|
| accessions matching the query | 891,382 |
| records retrieved | 891,371 |
| unavailable | 11, all WGS master records that carry no sequence |
| exact duplicates removed, keyed on (taxid, sequence) | 293,273 (33%) |
| PDB structure records dropped | 110 |
| **records in the database** | **597,988** |
| total bases | 688,449,026 |
| longest sequence | 81,206 bp |
| unpacked / compressed | 342 MB / 78 MB |

Download took about 2 h 20 min in 179 chunks of 5,000. Dereplication and
`makeblastdb` together take under 30 s; `makeblastdb` itself indexed 597,988
sequences in 6.3 s.

PDB structure records must be excluded, and the reason is user-facing rather than
cosmetic. They return deflined `>pdb|28ZU|S2 Chain S2, 18S ribosomal RNA` instead
of `>ACC.V Genus species ...`. `extract_species()` (`R/blast_genbank.R:243`)
takes the first two words of the title, so such a hit would set
`assemblies.blast_species` to "Chain 2," in the app. Their parsed seqid
(`pdb|28ZU|S2`) also fails to match the `28ZU_S2` key in the taxid map, so they
would carry no staxid. They are modelled ribosome chains, overwhelmingly human,
and ordinary GenBank records cover the same sequences.

Every other record in the build carries a plain `ACC.V` defline, so
`classify_hit_title()` and `extract_species()` need no change.

### Measured performance

On the 15,120 bp `SRR14433904` unit, against the full 688 Mbp panel:

| threads | wall | peak RSS |
|---|---|---|
| 1 | 5.46 s | 254 MB |
| 4 | 1.46 s | 256 MB |
| 8 | 0.80 s | 256 MB |

Against a remote `blast_genbank` median of 7,589 s per task, that is roughly
**9,500x** at 8 threads and still about 1,400x single-threaded. Peak RSS is flat
at ~256 MB regardless of thread count, so `maxForks` is bounded by cores rather
than by memory.

Concordance against the answers remote `core_nt` gave for the same units, on
`full_project_docker`:

| unit | remote rank-1 | local rank-1 | pident |
|---|---|---|---|
| SRR14163830.1.1 | MZ188965.1 | MZ188965.1 | 99.807 both |
| SRR14433895.1.1 | MZ286342.1 | MZ286342.1 | 99.986 both |
| SRR14433899.1.1 | MZ286337.1 | MZ286337.1 | 100.000 both |
| SRR14433900.1.1 | MZ286336.1 | MZ286336.1 | 99.957 both |

Exact agreement on accession and on identity to three decimals. `staxids`
populate and `blastdbcmd -entry` retrieves sequence plus taxid plus title, which
is what the on-demand annotation path needs in place of `efetch`.

The panel also returns `MZ286333.1` at 100.000 percent for the 15 kb unit whose
`qcovs` of 46 the candidate floor currently rejects, so the local path and the
Stage 1 filter fix address the same unit from two directions.

**This validates fish only.** Every sample in all three test projects is a
ray-finned fish. The escalation rate on Asteroidea, Scyphozoa, and Octocorallia
remains unmeasured and remains the acceptance gate.

**Not baked into the Docker image.** The image is 5.98 GB and its R-package
layers, which re-push on every version bump, are 537 kB and 1.16 MB. Coupling a
quarterly 220 MB artefact to a daily image cadence is the wrong trade, and
`docker/Dockerfile` already states the intent that reference databases resolve at
run time. It is mounted instead, which needs three things that do not exist
today: `ribopilot_db_dir()` added to the `extra_mounts` built at
`R/init_project.R:137` (without it a containerised run fails with "No alias or
index file found"), `:ro` appended to the reference mounts interpolated at
`R/generate_config.R:34,38`, and `params.blastdb` reaching
`inst/nextflow/modules/blast.nf`, which today hardcodes `-remote -db core_nt` in
both the megablast and the dc-megablast call.

A shared read-only install is the recommended HPC deployment and is blocked today
by `validate_ref_dir()` (`R/reference_dir.R:51`), which fails any directory
without write access. It needs a read-consumer path. Note this function is called
only from the settings-modal Save handler (`R/app_reference.R:770`), so it is not
on the BLAST path.

Distribution is a GitHub release asset: one `.tar.gz` plus a `.sha256`, tagged
`db-<version>`, published with `--latest=false` so it does not steal the badge
from the package release. Verified limits: 2 GiB per file, no total or bandwidth
cap, anonymous download, no LFS. Never re-upload over a published version:
`--clobber` deletes before uploading, so an interrupted upload destroys the
original and a checksum already recorded in someone's project would then point at
different bytes. A botched build gets a new tag. Verification is a shell helper
rather than R, because `tools::sha256sum()` needs R 4.5.0 while `DESCRIPTION`
declares 4.4.0 and the Dockerfile pins `r-base=4.4.*`, and neither `digest` nor
`openssl` is in Imports.

Provenance is recorded per project so a result can be traced to the database that
produced it: database version, checksum, and the effective search space, written
when the pipeline first uses the database.

### Refresh and staleness

Developers rebuild at version bumps. The manifest carries the build date, the
verbatim query, the record counts, the tool versions, and the artefact checksum.
The app surfaces the age of the local copy rather than forcing an upgrade, since
an older database is a reproducibility asset, not a defect: a project that ran
against `2026.08.1` should keep saying so. Old versions stay on disk under their
own version directory.

### Acceptance gate

Revised 2026-08-14. The original gate had three measurements; the maintainer
retired all three, with reasons recorded here so a later reader knows they were
decided rather than forgotten.

**1. Query wall time and peak RSS: answered, and taxon-independent.** Measured on
the built database: 0.80 s at 8 threads, 5.46 s single-threaded, peak RSS 256 MB
flat across thread counts, on a 15,120 bp query against 688 Mbp. Cost is set by
query length and database size, not by which clade the query came from, so the
fish measurement carries. `maxForks` is bounded by cores rather than memory.

**2. Escalation rate: no longer a gate**, because automatic escalation is
shelved. A unit the panel cannot answer is reported to the user rather than
silently re-searched, so there is no hidden remote cost to bound. What matters
instead is the *unanswerable rate*, which is a reporting quality question rather
than a shipping blocker: if most invertebrate units come back unanswered the
feature is unhelpful for those clades, but it is not harmful and it does not cost
hours.

**3. Local-versus-remote concordance: retired.** BLAST ranking is a heuristic,
and the pipeline only ever consumes rank 1 as a comparator, not as a
determination. Exact agreement with a `core_nt` search was never the standard the
remote path itself was held to.

What remains before Stage 4 ships is narrower: confirm on one invertebrate clade
that the unanswerable rate is reported honestly and that the surfaced message is
actionable. That needs one Asteroidea, Scyphozoa, or Octocorallia run, and unlike
the retired gate it does not need a remote baseline to compare against.

### Stage 4a acceptance measurement (2026-08-14)

Ran on copies of all four available validation projects against panel
`2026.08.1` (597,988 sequences, 688 Mbp, no decoy databases downloaded, so
panel-only): Asteroidea (`DRR251133`, *Acanthaster planci*), Scyphozoa
(`ERR3597993`, *Aurelia aurita*), Octocorallia/Keratoisis (`SRR25308289`, 2
units), and Actinopterygii (`SRR14163830`, pufferfish) as a fish control.
**n = 5 units total across 4 projects.** This is a small-n measurement; the
per-clade counts (1, 1, 2, 1) are far too few to estimate a rate, only to
observe individual outcomes.

| unit | clade | top pident | rank distance | anchor floor (18S) | signal | screen verdict |
|---|---|---|---|---|---|---|
| DRR251133.1.1 | Asteroidea | 91.491 | class | 93 | below_anchor_pid_18S | on_target |
| ERR3597993.1.1 | Scyphozoa | 96.751 | genus | 93 | answered | on_target |
| SRR25308289.1.1 | Octocorallia | 99.276 | family | 93 | answered | on_target |
| SRR25308289.1.2 | Octocorallia | 95.533 | family | 93 | answered | on_target |
| SRR14163830.1.1 | Actinopterygii (control) | 100.000 | species | 97 | answered | on_target |

Of the 5 units, 1 was unanswerable (the Asteroidea unit), by the identity
signal only: the top reference candidate (*Asterias amurensis*, `OR630824.1`)
cleared the reference-candidate floors and the taxonomy screen agreed it was
on-target, but its identity (91.491 percent) fell under the invertebrate 18S
anchor floor of 93. No unit in this run returned zero hits or no candidate at
all; every unit saturated `max_target_seqs` (50 distinct subjects each). The
`screen_reason` text a user would see for the unanswerable unit: "hits agree
on Asteroidea (class, 3 voting, mostly nuclear rDNA records); on the sample's
lineage. support 0.30 w 0.30, resolution 0.43 w 0.30, separation 1.00 w 0.40
-> confidence 0.62".

Local query wall time: 0.705-1.437 s per unit (mean 1.14 s), no dc-megablast
fallback triggered on any unit. Peak RSS: 242.8-255.9 MB per unit, consistent
with the flat ~256 MB reported above regardless of query.

Caveat: with n=5, this cannot support a general escalation-rate estimate for
invertebrates. What it shows is that the panel can answer genus- and
family-level invertebrate queries (Scyphozoa, Octocorallia) and can fail a
class-level one on identity alone even while the taxonomy screen still calls
it on-target (Asteroidea): the local path is not merely a fish-only win, but
neither is a single below-floor case grounds to claim it generalizes to
Asteroidea as a clade. More invertebrate projects are needed before an actual
rate can be reported.

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
