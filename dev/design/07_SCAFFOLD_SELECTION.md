# 07 - Scaffold / path selection and taxonomy-aware auto-select

## Status 2026-07-30: stitching was DROPPED, and this document was stale

The body below is the pre-pivot design and is kept for the reasoning, not as a
statement of what the code does. Commit `b0f1b12` ("feat(pipeline): WF1/WF2
per-scaffold workflows (drop stitching)") replaced the two-level scheme with a
per-scaffold model, and no document recorded that until now. Read this section
before the body.

**The decision.** Level 1 (order, orient and stitch a path's scaffolds into one
canonical repeat unit) is not built into the pipeline and will not be. The unit
of work is `(ID, path, scaffold)`. WF1 loads assemblies, seeds one unit per
scaffold, BLASTs each scaffold and splits the depth track per unit. WF2
annotates each non-ignored scaffold independently, pulling that unit's sequence
from the database. There is no stitch process in `main.nf` or
`main_annotate.nf`.

**Why.** Stitching commits to an ordering before the evidence supporting it has
been reviewed, and it destroys the very signal the contamination work needs. A
contaminant scaffold merged into a host repeat unit becomes one sequence with
one taxonomy and one verdict; kept separate, it gets its own BLAST hits, its own
screen verdict and its own keep/drop decision. Everything in
`08_MULTIUNIT_UX_AND_CONTAMINATION.md` onward assumes the per-scaffold grain.

**What this changes in the body below.**

- The "Level 1 - region-label + stitch (per path)" section describes work that
  is not in the pipeline. Its first step, annotating each scaffold to the
  elements it covers, survives as ordinary per-scaffold annotation. Steps 2 and
  3 do not.
- "Implementation status (built)" ends with "Remaining: integrate stitch into
  `score_assembly`". That is the opposite of what happened. It was removed, not
  integrated.
- Recommendation P1 lists "Level-1 stitch" as part of the automation core. It is
  not.
- The three tangled situations in "The problem" are all still real, but only the
  framing changed for the first. Fragmentation is now surfaced to the user as
  several units of one sample rather than silently merged. The other two are
  handled by triage rules (`09_TRIAGE_RULES_AND_DEPTH.md`), the screen verdict
  layer (`11_SCREEN_VERDICT_LAYER.md`) and the isomer/paralog classification
  (`12_INTRAGENOMIC_VARIANTS.md`).

**What remains in the tree, and what it is worth.** The stitching code was not
deleted. `stitch_scaffolds()` and `order_scaffolds_by_content()`
(`R/stitch_scaffolds.R`), `stitch_one()` and `stitch_selected()`
(`R/pipeline_steps.R`), `inst/nextflow/bin/rp_stitch.R` and the `stitch` process
in `modules/pipeline_r.nf` are all still present. None of them is reachable from
the main pipeline. The only remaining caller chain runs through
`main_fallback.nf`, which doc 08 open question 1 records as broken post-pivot,
plus the exported legacy driver `annotate_project()`, which has zero callers
anywhere. `stitch_scaffolds()` and `stitch_one()` are both exported, so they are
public API and cannot simply vanish. No test covers any of it.

Deleting the stitching machinery therefore depends on a decision this document
cannot make on its own: whether the fallback workflow is still supported (doc 08
open question 1) and what happens to `path_scores`, the `stitched` table,
`stitch_one` and `rp_stitch.R` (open question 2, which says explicitly to decide
them together rather than piecemeal). Until that is answered, the machinery
stays and this section is the record of why.

**The `stitched` table is now an orphan.** `R/init_db.R:447` creates it on every
fresh project and no migration drops it. It is written only by the two
unreachable paths above. Its only reader was `compute_feature_alignment()`,
deleted on 2026-07-30 for the reasons doc 08 gives: it read `stitched` and
queried `annotations` with no path or scaffold predicate, so under the
per-scaffold model it mixed features from different units of one sample into a
single alignment. Its successor `compute_feature_msa()` in the same file does
the same job correctly, keying on `assemblies` with full `ID`, `path` and
`scaffold` predicates. So in a normal run `stitched` is now created, never
written and never read. Dropping it is part of open question 2, not a separate
decision.

## The problem

Unlike a mitogenome (one circular molecule that usually resolves to a single
path), the nrDNA tandem repeat graph is messy and GetOrganelle nr routinely
emits **several paths, each with several scaffolds**. Three different things get
tangled together and must be told apart:

1. **Fragmentation** - scaffolds within one path are *different regions of the
   same rDNA repeat unit* (e.g. our fish test: scaffold_3 = 18S + ITS + partial 28S,
   scaffolds 1-2 = 28S pieces). These should be *stitched*, not chosen between.
2. **Alternative paths** - competing graph traversals of the same locus
   (near-duplicates, or genuine intragenomic variants). Pick a representative,
   collapse the rest.
3. **Foreign / paralog contigs** - contaminant rDNA (fungal/algal/microbial) or
   nuclear rDNA pseudogenes. These should be *dropped/flagged*, not selected.

MitoPilot's "rank paths by BLAST + coverage, user locks one" is necessary but not
sufficient here: we also need per-scaffold region labeling + stitching, and a
taxonomy-aware discriminator to separate target from contaminant.

## Signals we already have per contig

- **Coverage / depth** - GetOrganelle's GFA carries `RC` (read count) and `LN`
  per segment, so depth = RC/LN is available for free. The true high-copy nuclear
  rDNA repeat unit should sit among the highest-depth contigs; contaminants and pseudogenes
  are typically lower. Strong discriminator.
- **BLAST identity + query coverage** vs an rDNA reference DB (reuse MitoPilot's
  candidate-BLAST machinery).
- **BLAST top-hit taxonomy / lineage** (the auto-select driver, see below).
- **Anchor annotation** - barrnap/Rfam (`--kingdom fun` = eukaryote Rfam models)
  + ITSx: which of 18S / 5.8S / 28S / ITS a scaffold contains, and in what order.
- **Length vs clade expectation** and **GC** within an expected band.

## Two-level selection

### Level 1 - region-label + stitch (per path) - DROPPED, see the status section

Before choosing between paths, make each path a single comparable object:

1. Annotate each scaffold to the element(s) it covers (barrnap/Rfam + ITSx +
   BLAST region).
2. Order + orient scaffolds along a reference rDNA sequence (a MitoPilot
   `scaffold_join` analog, driven by minimap2/BLAST to a chosen reference) into a
   canonical `18S-ITS1-5.8S-ITS2-28S` layout.
3. Emit one **stitched representative** per path (with gaps flagged where
   scaffolds don't abut), plus a coverage track.

This turns "3 scaffolds x 2 paths" into "2 stitched candidates", far easier to
score and to show the user.

### Level 2 - score + auto-select the best path

Composite score per stitched candidate (weights tunable, stored in curate_opts):

```
score =  w_tax  * taxonomic_match      # top-BLAST-hit lineage vs sample lineage
       + w_id   * blast_identity       # % identity of anchors to reference
       + w_cov  * anchor_completeness  # has 18S + 5.8S + 28S, in order
       + w_dep  * depth_rank           # coverage rank among candidates
       + w_len  * length_plausibility  # total within clade tolerance
       - p_flag * flag_penalties       # contaminant/chimera/paralog penalties
```

Auto-select the top-scoring path (analogous to MitoPilot auto-picking a BLAST
reference), pre-set it as the sample's chosen path, and surface the ranked
alternatives with their sub-scores so the user can override. Never hard-delete
losers - mark them ignored (MitoPilot pattern), keep them reviewable.

## Taxonomy-aware auto-select and the free-text `Taxon` problem

"Choose the scaffold that most closely matches the sample's taxonomy" needs a
**machine-comparable expected lineage** for the sample. Today `Taxon` is free
text with no format contract ("Upeneus parvus", "Mullidae", "reef fish sp.").
Options, roughly in order of robustness:

1. **Resolve `Taxon` -> NCBI taxid -> lineage at project init.** Use Entrez
   `esearch -db taxonomy` or `taxonkit name2taxid`; store the resolved lineage in
   `samples`. Then taxonomic match = lowest shared rank between a scaffold's
   top-BLAST-hit lineage and the sample lineage (species > genus > family > ...).
   Unresolvable/ambiguous names get flagged for the user at init, not silently
   dropped.
2. **Optional structured map-file columns** (backward compatible): `taxid`, or
   `order`/`family`/`genus`. If present, use directly; else fall back to (1).
3. **Graceful rank fallback scoring.** Don't require species match; score by the
   *depth of the shared lineage*. A scaffold whose best hit is congeneric scores
   higher than one matching only at order; one matching only at phylum (or a
   different phylum) is flagged as likely contaminant.
4. **Per-sample reference override.** Let the user name a reference accession (or
   pick one from BLAST hits) to define "the right answer"; auto-select then favors
   the candidate most similar to that reference. Mirrors MitoPilot's
   "Set as reference".
5. **Project-level expected clade.** A single `expected_taxon`/taxid for the whole
   project (e.g. "Actinopterygii") as a coarse contaminant gate when per-sample
   resolution fails.

Recommended: (1) as the backbone (resolve at init, store lineage, flag failures),
(2) as an optional override, (3) for scoring, (4) in the UI for hard cases.

**DECISION (2026-07-02): go with (1) - auto-resolve `Taxon` at init.** No new
required map-file columns. Implications:
- `new_project()` / `new_db()` gains an init step that resolves each `Taxon`
  string to an NCBI taxid + lineage (Entrez `esearch -db taxonomy`, or
  `taxonkit name2taxid | taxid2lineage`; entrez-direct is already in the env).
- `samples` table gains `taxid` and `lineage` columns (nullable).
- Names that fail to resolve or are ambiguous are recorded with an
  `ambiguous-taxon` flag and shown to the user to fix (edit `Taxon`, or supply a
  taxid) rather than silently dropped.
- Scoring uses rank-fallback (option 3): depth of shared lineage between a
  scaffold's top-BLAST-hit lineage and the sample lineage.
- Keep option (2) taxid override and (4) UI "Set as reference" as later add-ons.

## MSA-against-references review UI

For eyeball selection when automation is uncertain:

- For each candidate (or per region: 18S, 28S), fetch the top-N GenBank reference
  rDNA repeat units via BLAST (reuse MitoPilot's ref-fetch), build an MSA
  (MAFFT or DECIPHER), and render it in-app with `msaR` (already a MitoPilot dep,
  used for AA alignments - here nucleotide).
- Columns/rows annotated with element blocks (18S/ITS/5.8S/ITS/28S) so the user
  sees where each candidate aligns and where it diverges/has gaps.
- Per-reference identity + taxonomy shown alongside; user clicks the best
  alignment to select that path, or "Set as reference".
- Complement with a **synteny/dotplot** of candidate vs reference (MitoPilot
  already has synteny plots) to spot chimeras (scaffold regions hitting different
  references/taxa) at a glance.

## Automated flag catalog (surface in the app, MitoPilot-style)

- **contaminant**: top-hit lineage disagrees with sample lineage (shallow shared
  rank) or blobtools/coverage screen disagrees.
- **chimera**: different regions of one scaffold hit different taxa/references.
- **paralog/pseudogene**: extra divergent low-depth copy alongside a dominant one.
- **incomplete**: missing an expected anchor (no 18S, or no 28S), or large
  inter-scaffold gaps after stitching.
- **length-outlier**: total length outside clade tolerance.
- **ambiguous-taxon**: `Taxon` could not be resolved to a taxid at init.
- **low-confidence-autoselect**: top two candidates score within a small margin
  (nudge the user to review the MSA).

## Reducing candidates at the source (cheap wins)

- **Coverage floor**: drop contigs below a fraction of the median rDNA repeat depth
  before selection (kills many pseudogene/contaminant fragments).
- **Dedupe near-identical paths**: collapse paths with >99% pairwise identity to
  one representative, so the UI shows genuine alternatives only.
- **GetOrganelle tuning**: `--reduce-reads-for-coverage`, depth caps, k-list and
  `--max-rounds` adjustments to get cleaner graphs; expose as assemble_opts.

## Implementation status (built)

- `R/taxonomy.R` - typo-tolerant NCBI resolution (per-project pinned taxdump,
  local taxonkit, network fallback); `resolve/apply_project_taxonomy()` confirm +
  freeze into `samples` (taxid/lineage).
- `R/scaffold_scoring.R` - `parse_gfa_depth()`, `taxonomic_match()` (rank
  fallback), `score_candidates()` composite scorer.
- `R/select_paths.R` - `score_assembly()` extracts live depth/anchors/BLAST +
  resolved lineage and scores/ranks/auto-selects.
- `R/stitch_scaffolds.R` - Level-1 stitch: order/orient/dedup scaffolds along a
  reference into one canonical rDNA repeat unit.
- DB: `path_scores` table + `write_path_scores()`.

Remaining: integrate stitch into `score_assembly` (needs the per-path top-hit
reference), the driver step that runs it after assembly, and the Shiny review UI.

**Superseded 2026-07-30.** The stitch was never integrated and will not be; see
the status section at the top. `b0f1b12` removed it from the pipeline.

## Recommendation / phasing

- **P1 (automation core):** per-scaffold BLAST + region label + depth from GFA;
  Level-1 stitch (**dropped, see the status section**); Level-2 composite score
  with the coverage + identity + anchor
  terms; auto-select + ignore-not-delete losers; the flag catalog. Taxonomy term
  uses whatever lineage is available (project-level at first).
- **P2 (taxonomy resolution):** resolve `Taxon` -> taxid -> lineage at init; add
  optional `taxid` map column; per-sample rank-fallback taxonomic scoring; flag
  unresolved.
- **P3 (review UI):** MSA-vs-refs picker + synteny + flag review + "Set as
  reference" override, wired into the Assemble tab.

This keeps the MitoPilot look/flow (ranked candidates, lock/switch/hide/ignore,
synteny, set-as-reference) while adding the two things nrDNA needs that MitoPilot
lacks: **per-scaffold region stitching** and **taxonomy-aware auto-select**.
