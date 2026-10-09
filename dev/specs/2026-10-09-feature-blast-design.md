# Feature BLAST: per-feature BLAST check

Date: 2026-10-09. Branch: `feat/feature-blast` (off `devel`).

## Purpose

A contig can look fine as a whole while one annotated feature (often ITS1 or
ITS2) has no BLAST hits, hits non-rDNA, or hits a distant taxon. That points to
a poor assembly or contamination that the contig-level screen, curation, and
chimera checks can miss: the contig BLAST keeps the top 50 hits, and those are
dominated by 18S and 28S. This feature BLASTs each core feature on its own and
reports problems in Annotation Review only.

## Decisions (user-approved)

- Features searched: 18S, ITS1, 5.8S, ITS2, 28S (`RRNA_GENES` + `RRNA_SPACERS`).
  ETS and IGS are skipped.
- Length rules:
  - under 50 bp: no search; status `too_small`, shown as "no BLAST, feature too small".
  - 50 to 99 bp: `-task blastn` (word 11); result marked low confidence.
  - 100 bp and up: `-task megablast`, retry with `dc-megablast` when there is no hit,
    the same as the contig BLAST.
- Database: the same local bundled panel the contig BLAST uses (`dbarg` from
  `rp_dbarg.R` / `local_blast_db_arg()`). Never NCBI. No panel (or
  `RIBOPILOT_NO_LOCAL_BLAST`): status `not_run` with note "no local BLAST panel".
- Off-target cut: phylum. A feature is off-target when its call shares no taxon
  with the sample at phylum or deeper (`rank_distance >= 6`).
- Shown only in Annotation Review. No triage rule, curation flag, screen vote, or
  export effect.
- Table column name: "Feature BLAST warnings". Problems only; blank when all
  features are fine.

## Thresholds rationale

With Karlin-Altschul statistics an exact match is significant (E < 1e-5) from
about 25 to 30 bp against a panel of 50 to 500 Mbp, so statistics are not the
limit. The limit is taxonomic signal: short pieces match many taxa equally, and
one mismatch costs a large share of the score. NCBI recommends `blastn-short`
under 50 bp, which brings an E threshold of 1000 and more noise; we skip that
range instead. Mini-barcode work (Hajibabaei et al. 2006; Meusnier et al. 2008)
puts reliable identification at about 100 bp and up, hence the 50 to 99 bp
low-confidence band. Citations to be checked before they go in user docs.

## Per-feature calls

Each searched feature gets one `status`:

| status | meaning |
|---|---|
| `too_small` | under 50 bp, not searched |
| `not_run` | no local panel, or the search failed (note says which) |
| `no_hit` | nothing passed `-evalue 1e-5` |
| `non_rdna` | best hit's `classify_hit_title()` is not `rdna` |
| `off_target` | call shares nothing with the sample at phylum or deeper |
| `ok` | on target at phylum or deeper |

Tie collapse (reuse of the chimera window logic, `chimera_window_hits()` in
`R/chimera.R`): among a feature's hits with lineages, if `bitscore_weighted_lca`
with `frac = CHIMERA_TIE_FRAC` (0.02), `min_hits = 1` covers 2 or more hits and
its rank is not species, the call is the shared clade, cut at that rank and
labeled `"<clade> (N tied hits)"`. Otherwise the call is the best hit. The
collapsed lineage, not the single best hit, is compared with the sample.

`match_rank`: `rank_distance_label(rank_distance(call_taxids, sample_taxids))`,
e.g. "species", "genus", "family", ... "none". When the sample has no
`lineage_taxids`, `match_rank` is NA and the feature is never `off_target`
(it stays `ok` with note "sample lineage unknown").

`low_conf`: 1 when the feature is 50 to 99 bp, else 0.

Rank-shared vs `ok`: an `ok` feature that only agrees at phylum or class is still
`ok`; the modal shows its match rank. The table column lists only `no_hit`,
`non_rdna`, `off_target`, and `too_small`. `not_run` is not a warning in the
column; the modal shows it.

## Implementation units

### 1. Core logic: `R/feature_blast.R` (new)

- `feature_blast_targets(con, id, path, scaffold, review)`: core features of the
  unit with start, stop, strand, length, and the sequence cut from
  `assemblies.sequence` or `review_units.sequence` (via `review_annotations_v` /
  `review_units_v` when `review = TRUE`). Reverse-complement minus-strand features.
- `feature_blast_task(len)`: `NA` under 50, `"blastn"` for 50 to 99,
  `"megablast"` for 100 and up. Constants `FEATURE_BLAST_MIN_LEN = 50L`,
  `FEATURE_BLAST_LOWCONF_LEN = 100L`, `FEATURE_BLAST_EVALUE = 1e-5`,
  `FEATURE_BLAST_MAX_TARGETS = 50L`, `FEATURE_BLAST_KEEP_HITS = 10L`.
- `feature_blast_calls(hits, targets, sample_taxids, taxdump_dir)`: pure
  function from a parsed hit frame to one summary row per feature (status,
  call fields, match_rank, low_conf). Holds the tie collapse. Testable without
  BLAST.
- `feature_blast_unit(con, id, path, scaffold, dbarg, taxdump_dir, threads,
  workdir, review = FALSE)`: writes one FASTA per task group, runs blastn with
  the contig BLAST outfmt (reuse `chimera_run_blastn()` or a twin with
  `-evalue` added), dc-megablast retry for megablast features without a hit,
  parses with `parse_blast_hits_file()`, calls `feature_blast_calls()`, and
  replaces the unit's rows in `feature_blast` and `feature_blast_hits` (or the
  review mirrors) in one savepoint. Returns the summary frame.
- `feature_blast_warnings(rows)`: the column text, e.g.
  `"ITS1: no hit; ITS2: off-target (Ascomycota); 5.8S: no BLAST, feature too small"`,
  in `RRNA_ELEMENT_ORDER`. Empty string when nothing to report. Off-target names
  the call's phylum (or deepest named rank when phylum is empty).

### 2. Schema: `R/db_schema.R`, `R/init_db.R`, `R/db_migrate.R`

`FEATURE_BLAST_DDL`, PK (ID, path, scaffold, element):
`ID, path, scaffold, element, start, stop, length, task, low_conf, status, note,
accession, species, staxid, lineage, lineage_taxids, n_tied, match_rank,
rank_distance, pident, qcovs, evalue, bitscore, db_source, time_stamp`.

`FEATURE_BLAST_HITS_DDL`, PK (ID, path, scaffold, element, hit_rank):
`ID, path, scaffold, element, hit_rank, accession, title, staxid, lineage,
pident, qcovs, evalue, bitscore, hit_class` (top `FEATURE_BLAST_KEEP_HITS` per
feature by bitscore, one row per subject).

- `new_db()` creates both.
- `REVIEW_MIRRORS` gains `feature_blast = "review_feature_blast"` and
  `feature_blast_hits = "review_feature_blast_hits"`; `review_apply_schema()`
  creates the mirrors and views `review_feature_blast_v`,
  `review_feature_blast_hits_v` via `review_view_ddl()`.
- `db_migrations()` step 5: create both tables, then drop and recreate review
  views as needed and run `review_apply_schema(con)`. Schema version becomes 6.
- Check that `review_drop()`, `review_invalidate()`, and any copy-creation code
  that iterates `REVIEW_MIRRORS` handle the new mirrors (copying originals into
  the mirror when a copy is made, clearing on drop).

### 3. Main workflow

- New `inst/nextflow/bin/rp_feature_blast.R` (args: db, id, path, scaffold,
  dbarg, taxdump_dir), modeled on `rp_chimera.R`, calling
  `feature_blast_unit()` then `review_follow_original()`.
- New process `feature_blast` in `inst/nextflow/modules/pipeline_r.nf`
  (`errorStrategy 'ignore'`, same DBARG logic as `chimera`).
- `main.nf`: `feature_blast(reconcile.out)`, in parallel with `chimera`. Not
  joined into triage. It needs only annotation and the panel, so no barrier.
- Requeued units (`query_annotate_units` path) go through `reconcile`, so they
  get it too.

### 4. Recheck workflow

- `recheck_unit()` in `R/recheck.R`: after `chimera_unit(... review = TRUE)`,
  call `feature_blast_unit(..., review = TRUE)` with the same `dbarg`. In map
  mode `dbarg` is NULL, so rows are written as `not_run`. A failure inside it is
  caught and stored as `not_run` with the error note; it does not fail the
  recheck.
- `rp_recheck.R` / `recheck.nf` need no new args if `dbarg` and `taxdump_dir`
  already pass through (verify).

### 5. Annotation Review table

- `all_rows` Annotate branch (`R/app_server.R` ~211-245): read
  `review_feature_blast_v` for the units and add `feature_blast_warnings`
  computed by `feature_blast_warnings()`.
- Column def: name "Feature BLAST warnings", help text in `SCAFFOLD_COL_HELP`,
  added to the "Check details" column group near `chimera_label`. Hidden-by-
  default status follows `chimera_label`.
- Stale marker: a feature whose stored start/stop no longer equals its current
  `review_annotations_v` start/stop reads `"<element>: stale, rerun recheck"`.

### 6. Annotate details modal

- New panel "Feature BLAST" in `annotate_modal()` (`R/app_modals.R` ~400),
  after the Chimera check panel.
- Table (reactable), one row per core feature: Feature, Length, Status (badge),
  Call (species or "<clade> (N tied hits)"), Match rank, Identity, Coverage,
  E-value, with "low confidence" note for 50 to 99 bp and the stale note.
- Expandable row detail: the stored top hits (accession, title, identity,
  coverage, bitscore, class), modeled on the assemble modal's BLAST evidence
  panel (`unit_candidates`).
- Reads `review_feature_blast_v` / `review_feature_blast_hits_v`, following the
  same check-state rule the chimera panel uses (`chim_unit` ~2670).

### 7. Docs

- Short section in the Annotation Review article describing the column, the
  panel, and the length rules. Use "contig" in user text.
- NEWS.md entry under the development version.

## Testing

- `test-feature-blast.R`:
  - `feature_blast_task()` boundaries (49, 50, 99, 100).
  - `feature_blast_calls()` on synthetic hit frames with the `helper-taxdump.R`
    taxdump: no hit, non-rDNA title, on-target at genus, off-target (other
    phylum), tied hits collapsing to a clade that contains the sample (ok) and
    one that does not (off_target), sample with no lineage.
  - `feature_blast_warnings()` text and order.
  - `feature_blast_unit()` end to end with a tiny `makeblastdb` panel built on
    the fly (skip if blastn absent), main and `review = TRUE`.
- `test-backwards-compat.R` parity passes with step 5.
- `test-recheck.R`: a recheck writes `review_feature_blast` rows; map mode
  writes `not_run`.
- `test-nextflow.R`: the new process is wired (follow existing wiring checks).
- App check: headless drive (chromote) of a test project, confirm the column and
  the modal panel render.

## Out of scope

Triage, curation flags, screen votes, export blocking, NCBI fallback,
`blastn-short` for features under 50 bp, ETS/IGS.
