# 21 - Chimera detection: per-gene taxonomy, windowed BLAST, ribotyper, read evidence

Status: DESIGN agreed 2026-09-01. Not implemented. Fills the slot doc 08 calls
sub-project 3d (the positional and chimeric channel) and the "second,
taxonomy-free junction detector" in doc 08 section 5. Does NOT touch the screen
verdict enum (doc 11 keeps `chimeric` out of it on purpose).

## 0. The one-paragraph summary

A scaffold assembled from a genome skim can be a join of two things: host 18S
fused to a fungal ITS-28S, host rDNA fused to a symbiont's, two divergent host
array copies spliced through the conserved core, or a tandem misjoin inside the
array. Today the only detector is ITSx's `Chimeric!` token, which is nearly
silent on Metazoa. The pipeline already stores everything a detector needs and
consumes none of it: positional HSPs from BLAST, annotated gene spans, a
per-position depth track, and a local rDNA panel with taxids baked into the
image. This design adds one Nextflow process that computes four independent
signals per unit, stores per-window taxonomy with full lineage text, collapses
it into human-readable lineage segments, scores a junction only when channels
corroborate, and draws the result under the existing depth plot in the annotate
modal. Chimera calls warn. They never block.

## 1. Decisions

- One new process `chimera`, after `reconcile`, gated on `coverage` and
  `blast_genbank` completion with the `collect().map{'ready'}.ifEmpty('ready')`
  sentinel already used for reference alignment in `inst/nextflow/main.nf`.
- Windowed search target is the local panel plus decoys (`ribopilot_rrna`,
  `16S_ribosomal_RNA`, `SSU_eukaryote_rRNA`, `LSU_eukaryote_rRNA`) via
  `rp_dbarg.R`. Windowed remote BLAST stays rejected (doc 08 section 3.4). If no
  panel resolves, windows run against the clade label DB for identity only and
  channel 2 becomes identity-only, which can localize but never score alone.
- Ribovore is used as `ribotyper` only, on the annotated 18S and 28S spans.
  `riboaligner` is not run: pybarrnap and ITSx already give boundaries.
- Read evidence is depth steps from the stored track plus soft-clip pileups from
  one extra samtools pass inside `coverage.nf`. No BAM is published.
- Lineage is surfaced everywhere the score is: per window, per segment, in the
  reason string, in the tooltip, in the plot legend.
- Taxonomy-only evidence caps at `suspect`. Only read evidence or ribotyper can
  raise a unit to `likely`. This is the doc 08 rule about 18S-vs-28S reference
  asymmetry made mechanical.
- Thresholds are constants in `R/chimera.R` for v1. No options UI block.
- Pre-release: new tables are `CREATE TABLE IF NOT EXISTS` DDL constants in
  `R/db_schema.R`, no migration.

## 2. Inputs per unit

| Input | Source | Missing behavior |
|---|---|---|
| sequence | `assemblies.sequence` | process skipped |
| features | `annotations` (18S, 5.8S, 28S, ITS1, ITS2) | channels 1 and 3 unavailable |
| HSPs | `blast_hits` (qstart, qend, bitscore, staxid, lineage, rank_distance) | channel 1 unavailable |
| depth track | `coverage.depths` | channel 4a unavailable |
| clip track | `coverage.clips` (new) | channel 4b unavailable |
| sample lineage | `samples.lineage_taxids` | rank distance undefined; channel 2 identity-only |
| panel | `rp_dbarg.R` discovery | channel 2 identity-only against label DB |

A channel that is unavailable reports `availability = "unavailable"` and
contributes nothing. It is listed as such in `channels` so the curator can see
what the call did not rest on.

## 3. Channels

Each channel returns `score` in [0, 1], `junctions` (integer positions, may be
empty), and `reason` (prose naming taxa).

### 3.1 Per-gene taxonomy from stored HSPs

No new search. For each annotated element, take the HSPs whose query interval
covers at least 50% of the element span. Compute a bitscore-weighted LCA over
one vote per accession (same rule `screen.R` uses). Store per element: taxon
name, rank, taxid, top species, mean pident, n_hits.

Signal: two elements resolve, with at least 3 votes each, to lineages whose LCA
is above class (rank distance between the two element LCAs >= 5 on the
`taxonomy_ranks.R` scale). Score 0.6, junction = midpoint of the gap between the
two elements. An element with no covering HSP is `unresolved` and never votes.
Never fires on a gap.

### 3.2 Sliding-window local BLAST

Tile the scaffold into windows of `CHIMERA_WIN = 500` bp at step
`CHIMERA_STEP = 250`; last window is right-anchored to the scaffold end. All
windows go out as one multi-FASTA in one `blastn -task megablast` call, outfmt
identical to `blast.nf` so `parse_blast_hits_file()` is reused. Windows with no
hit are re-run once as one `dc-megablast` call. Per window keep the best hit by
bitscore: accession, species, taxid, full lineage text, pident, bitscore,
`hit_class` (from `hit_class.R`), rank distance to the sample.

Collapse adjacent windows into **segments**: a run of windows whose best-hit
lineage agrees at the deciding rank. Deciding rank is class for metazoan hits,
kingdom for decoy hits. A segment records start, stop, taxon, rank, top
species, mean pident, n_windows, and `db_source`.

Signal: at least two segments each of >= 2 windows, adjacent, with lineages
differing at class or above, and the shorter one at least 25% of the longer one's
length. Score 0.5. Junction = boundary between the segments, refined to the
window step. Identity-only mode (no panel): segments are runs of pident above
or below `CHIMERA_PID_DROP = 90`, score 0, junctions still reported for
corroboration.

Rationale for the length requirement and the cap: doc 08 section 3.5 (d), a
200 bp sliver must not outvote a 2 kb on-target block, and GenBank 28S coverage
across metazoan phyla is far sparser than 18S.

### 3.3 Ribotyper

Extract the annotated 18S and 28S spans (with 50 bp flanks) to FASTA and run
`ribotyper` with the eukaryotic model set. Parse `*.ribotyper.short.out` for the
`unexpected_features` column.

| Code | Score | Meaning |
|---|---|---|
| InconsistentHits | 0.8 | hits out of order vs model, misjoin |
| DuplicateRegion | 0.8 | same model region twice, tandem misjoin |
| BothStrands | 0.8 | strand flip inside the gene |
| MultipleFamilies | 0.6 | more than one rRNA family in one span |
| LowScoreDifference | 0.2 | weak taxon signal, informational |

Junction = model coordinate of the failing hit boundary mapped back to scaffold
position. The existing ITSx `annotate.chimeric` flag joins this channel at score
0.2 with no junction, so it stops being a separate curation flag source.

### 3.4 Read evidence

**4a Depth steps.** From `coverage.depths`, median depth per 250 bp window.
Flag a boundary where the ratio between adjacent window medians exceeds
`CHIMERA_DEPTH_STEP = 3` and the boundary lies more than 100 bp from any
annotated feature edge. Score 0.4 per step, junction at the boundary.

**4b Soft-clip pileups.** `coverage.nf` adds one pass after `samtools depth`:
`samtools view` piped to awk, counting per position the reads whose CIGAR has a
soft clip >= 20 bp starting or ending there. Written to `clips.txt` in the same
ref/pos/count layout as `depth.txt`, stored in `coverage.clips` as a comma
track by `write_coverage_units()`. Flag a position where clipped reads / depth
>= `CHIMERA_CLIP_FRAC = 0.3` and depth >= 10, excluding 50 bp at each scaffold
end. Score 0.7, junction at the position.

Read evidence is the only channel that sees a host-variant x host-variant join,
which is why it can reach `likely` alone.

## 4. Scoring

```
score = max(channel scores)
if two channels place junctions within 300 bp: score = min(1, score + 0.25)
if only channels 1 and 2 contributed: score = min(score, 0.7)
label = clean (< 0.3) | suspect (0.3 to 0.7) | likely (> 0.7)
```

`reason` is built from segments and junctions, e.g.
`1-2500 Actinopterygii (Salmo salar, 99.1%); 2750-4100 Ascomycota (Penicillium
sp., 97.4%); junction near 2620 (taxonomy, clip pileup 2611)`.

## 5. Schema

`CHIMERA_DDL` and `CHIMERA_WINDOWS_DDL` in `R/db_schema.R`, created by
`new_db()`.

```
chimera (
  ID, path, scaffold,            -- PK
  score REAL, label TEXT,
  junctions TEXT,                -- comma list of positions
  reason TEXT,
  segments TEXT,                 -- JSON array: start, stop, taxon, rank, species, pident, n_windows, db_source
  gene_taxa TEXT,                -- JSON: element -> taxon, rank, taxid, species, pident, n_hits
  channels TEXT,                 -- JSON: channel -> score, availability, junctions, reason
  ribotyper_codes TEXT,
  db_source TEXT,                -- panel prefix or label DB path
  time_stamp INTEGER
)

chimera_windows (
  ID, path, scaffold, start,     -- PK
  stop, accession, species, staxid, lineage, rank_distance, hit_class,
  pident REAL, bitscore REAL, depth_median REAL, clip_frac REAL, db_source
)
```

`coverage` gains `clips TEXT`.

## 6. Downstream surfaces

- Triage rule `chimeric` in `R/triage_rules.R`: actions off/flag, default flag,
  fires on `label == "likely"`, text = `reason`.
- `submission_risks()` in `R/export.R` reads `chimera.label` and reports
  `suspect` and `likely` as warn-only risks with `reason`. No export block.
- `curate_rrna.R`'s `chimeric` flag is retired in favor of the triage rule.
- Scaffold table: `chimera` column, badge by label, tooltip = `reason`. Help
  text in `R/constants.R`.
- Annotate modal, under the depth plot, `chimera_plot` (ggplot2, shared x):
  - segment bar: one rectangle per segment, fill by taxon at deciding rank,
    host-concordant segments grey, label inside when width allows;
  - window tiles: thin row, alpha by pident;
  - junction marks: dashed verticals, color by channel;
  - feature rectangles as on the depth track;
  - legend: taxon name, rank, mean pident per segment.
- Below the plot: reactable of segments; collapsed reactable of all windows
  with full lineage text.
- Junction positions appear as suggested cut points in the existing trim UI.

## 7. Pipeline

- `inst/nextflow/modules/coverage.nf`: add clip pass; `rp_coverage.R` takes an
  optional fourth arg `clips.txt`.
- `inst/nextflow/modules/pipeline_r.nf`: process `chimera`, `cache false`,
  `errorStrategy 'ignore'`, `tag "${id}.${upath}.${uscaf}"`, input
  `tuple val(id), val(upath), val(uscaf)`, script `rp_chimera.R <db> <id>
  <path> <scaffold> <dbarg>`.
- `inst/nextflow/main.nf`: `chim_done` sentinel from `coverage.out` and
  `blast_genbank.out`; `chimera(reconcile.out.combine(chim_done).map{...})`.
  Triage's group barrier must also wait on `chimera.out` so the rule sees the
  label.
- `inst/nextflow/bin/rp_chimera.R`: thin driver, same shape as `rp_screen.R`.
- `R/chimera.R`: tiling, blastn call, segment collapse, ribotyper call and
  parse, channel functions, scoring, writes.
- `docker/Dockerfile`: add pinned `ribovore` from bioconda to the mamba env.

## 8. Testing

- Unit tests in `tests/testthat/test-chimera.R`: tiling edge cases (short
  scaffold, last window anchoring), HSP-to-element intersection with 50%
  coverage rule, segment collapse, length-ratio guard, scoring caps, ribotyper
  short-output parser on a fixture file, clip track parser.
- Synthetic fixtures in `tests/testthat/fixtures/chimera/`: host 18S + fungal
  ITS-28S splice; host-variant join with a 6% identity drop in 28S; tandem 18S
  duplication; clean scaffold. Assert label and that a junction lies within 300
  bp of the splice.
- Property: with only channels 1 and 2 available, score never exceeds 0.7.
- Nextflow smoke on the existing test project: `chimera` rows exist for every
  reconciled unit; run completes when coverage or BLAST failed for a sample.
- Doc 13 validation fixture at line 309 (the known false negative) is re-run
  and its expected outcome updated.

## 9. Out of scope

- Cross-sample and index-hopping detection (doc 08 section 3.6, blocked on
  `samples.batch`).
- Automatic trimming at a junction. The trim UI exists; the call stays human.
- Export blocking on any chimera label.
- Per-position minor-allele fraction (doc 08 section 5).
- A `chimera_opts` block in the options UI.
