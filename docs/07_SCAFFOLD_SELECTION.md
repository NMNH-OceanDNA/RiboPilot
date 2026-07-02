# 07 - Scaffold / path selection and taxonomy-aware auto-select

## The problem

Unlike a mitogenome (one circular molecule that usually resolves to a single
path), the nrDNA tandem repeat graph is messy and GetOrganelle nr routinely
emits **several paths, each with several scaffolds**. Three different things get
tangled together and must be told apart:

1. **Fragmentation** - scaffolds within one path are *different regions of the
   same cistron* (e.g. our fish test: scaffold_3 = 18S + ITS + partial 28S,
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
  cistron should sit among the highest-depth contigs; contaminants and pseudogenes
  are typically lower. Strong discriminator.
- **BLAST identity + query coverage** vs an rDNA reference DB (reuse MitoPilot's
  candidate-BLAST machinery).
- **BLAST top-hit taxonomy / lineage** (the auto-select driver, see below).
- **Anchor annotation** - barrnap/Rfam (`--kingdom fun` = eukaryote Rfam models)
  + ITSx: which of 18S / 5.8S / 28S / ITS a scaffold contains, and in what order.
- **Length vs clade expectation** and **GC** within an expected band.

## Two-level selection

### Level 1 - region-label + stitch (per path)

Before choosing between paths, make each path a single comparable object:

1. Annotate each scaffold to the element(s) it covers (barrnap/Rfam + ITSx +
   BLAST region).
2. Order + orient scaffolds along a reference cistron (a MitoPilot
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

## MSA-against-references review UI

For eyeball selection when automation is uncertain:

- For each candidate (or per region: 18S, 28S), fetch the top-N GenBank reference
  cistrons via BLAST (reuse MitoPilot's ref-fetch), build an MSA
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

- **Coverage floor**: drop contigs below a fraction of the median cistron depth
  before selection (kills many pseudogene/contaminant fragments).
- **Dedupe near-identical paths**: collapse paths with >99% pairwise identity to
  one representative, so the UI shows genuine alternatives only.
- **GetOrganelle tuning**: `--reduce-reads-for-coverage`, depth caps, k-list and
  `--max-rounds` adjustments to get cleaner graphs; expose as assemble_opts.

## Recommendation / phasing

- **P1 (automation core):** per-scaffold BLAST + region label + depth from GFA;
  Level-1 stitch; Level-2 composite score with the coverage + identity + anchor
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
