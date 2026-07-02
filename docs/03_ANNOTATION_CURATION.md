# 03 - Annotation and curation

The whole locus reduces to **five conserved-anchor edges** plus the spacers
between them. If we place the anchors well, the spacers fall out by subtraction.

```
 [====== 18S ======][ ITS1 ][= 5.8S =][ ITS2 ][========= 28S =========]
                   e1       e2        e3      e4
 anchors (conserved, alignable): 18S, 5.8S, 28S core
 spacers (variable, "trim between anchors"): ITS1 = e1..e2, ITS2 = e3..e4
```

## Annotation engines (combine, don't rely on one)

### 1. barrnap - rRNA genes

barrnap predicts rRNA gene locations and outputs GFF3, giving anchored 18S,
5.8S, and 28S.

**Build note (found in Increment 2):** the bioconda `barrnap` that our
`environment.yml` resolves is the **HRGV fork v1.10.5**, whose kingdoms are
`bac arc fun` (no `euk`). Critically, this fork runs **Infernal `cmsearch`
against Rfam eukaryotic covariance models** (SSU_rRNA_eukarya RF01960,
LSU_rRNA_eukarya RF02543, 5_8S RF00002), not nhmmer. Invoke it with
`--kingdom fun`; despite the name it uses the eukaryote Rfam models and
correctly annotated 18S + 28S on the fish test assembly. This means the
structure-aware Rfam anchor calling listed below as an "optional high-accuracy
mode" is effectively built into our barrnap.

Increment 3 decision: standardize on this fork's `--kingdom fun` (Rfam-based),
or pin the classic tseemann `barrnap=0.9` (`--kingdom euk`, nhmmer). Leaning
toward the Rfam fork since it is structure-aware. Either way, ITSx remains the
authority for ITS1/ITS2 boundaries.

### 2. ITSx - spacer boundaries

ITSx uses profile HMMs built from ~20 eukaryote groups (fungi, metazoa, plants,
protists) to detect and extract **ITS1, 5.8S, ITS2** (and flag SSU/LSU ends). It
is HMMER3-based, handles reverse-complement, and is explicitly built for broad
eukaryotes. ITSx is the authority for the ITS1/ITS2 edges; barrnap is the
authority for the 18S/28S extents. We run both and reconcile.

### 3. Reference alignment - anchor + fill (DECIPHER)

Align the assembled unit to a curated reference cistron for the sample's clade
using DECIPHER (already a MitoPilot dependency). The reference provides:
- exact anchor coordinates when barrnap/ITSx disagree,
- 28S D1-D12 expansion-segment landmarks (optional sub-features),
- a similarity track to flag divergent/mis-assembled regions.

### 4. cmsearch / Rfam (optional, high-accuracy anchors)

Infernal `cmsearch` against Rfam covariance models (SSU_rRNA_eukarya RF01960,
LSU_rRNA_eukarya RF02543, 5_8S RF00002) uses secondary-structure-aware models
and is the most robust anchor caller across deep divergence. Heavier dependency;
keep as an optional high-accuracy mode.

## Reconciliation into one feature model

Priority when tools disagree (configurable per clade):

1. **18S 3' end / 28S 5' end** -> barrnap (or cmsearch if enabled), refined by
   reference alignment.
2. **ITS1 / ITS2 extents and 5.8S** -> ITSx primary; barrnap `--accurate` as
   tiebreaker.
3. **Spacers derived by subtraction**: `ITS1 = end(18S)+1 .. start(5.8S)-1`,
   `ITS2 = end(5.8S)+1 .. start(28S)-1`. This guarantees a gap-free, non-overlapping
   model, which is what GenBank wants.

Store every feature with its winning `annotator` and the competing calls in
`notes`, so curation is auditable (mirrors MitoPilot storing provenance).

## Reference-guided curation

Follow MitoPilot's proven design: **a shared core (`curate_rrna_core`) plus thin
per-clade wrappers (`curate_<clade>_rrna`) driven by per-clade params
(`params_<clade>_rrna`)**. Adding a clade = a params file + a 10-line wrapper.

The core does, per assembly:

1. **Anchor check** - confirm 18S, 5.8S, 28S each present, in order, on one
   strand, each within the clade's expected length tolerance. Missing/short
   anchor -> hard flag.
2. **Spacer derivation** - compute ITS1/ITS2 by subtraction; sanity-check lengths
   against clade tolerances (ITS lengths vary a lot, so tolerances are wide and
   clade-specific; the check catches only gross errors like a 3 kb "ITS1").
3. **Boundary refinement** - snap anchor edges to the reference alignment where
   the annotators are fuzzy; record the shift.
4. **Divergence/assembly-quality flags** - reference identity below threshold in
   a conserved anchor -> likely misassembly/pseudogene/contaminant, flag.
5. **Contaminant lineage check** - top BLAST-hit lineage vs `Taxon` (reuse
   MitoPilot's lineage columns); mismatch -> flag.
6. **Emit** feature table + per-feature status (ok / edited / flagged) for the
   Shiny reviewer.

### Why this is harder than MitoPilot, and how we cope

- No codon/frame/translation checks. **Replacement signal = reference-profile
  identity + HMM boundary confidence + structural covariation (cmsearch).**
- Spacers are unalignable across distant taxa. **We never try to base-align
  spacers; we only bound them by anchors and length-check them.** The scientific
  product (ITS1/ITS2 sequence) is delivered as a trimmed span, not a curated
  base-by-base model.
- 28S expansion segments differ by kilobases. **Treat 28S as one feature by
  default; expose D-domain landmarks only as optional annotations for users who
  want them.**
- Intragenomic variation. **Curate the dominant-coverage path; surface variants
  as alternates, not errors.**

## Curation UI (Shiny)

- Linear element map with 4 draggable internal breakpoints; dragging an anchor
  edge live-updates the neighboring spacer span.
- Reference-alignment overlay + identity track under the map to guide edits.
- Per-feature flag list with jump-to-position (reuse MitoPilot's "flag ->
  annotate modal" pattern).
- Lock/switch/hide per sample as in MitoPilot; edits persist to DB via the atomic
  writer.

## Export (GenBank)

Reuse MitoPilot's export pipeline. Feature-table specifics:
- 18S / 5.8S / 28S -> `rRNA` features with `/product="18S ribosomal RNA"` etc.
- ITS1 / ITS2 -> `misc_RNA` (or `misc_feature`) with
  `/note="internal transcribed spacer 1"`.
- Source modifiers from map-file metadata (organism, voucher, country, ...).
- Emit `.fsa` + `.tbl` for tbl2asn/table2asn, plus a combined fasta of each
  extracted marker (18S, 28S, full ITS, ITS1, ITS2) for barcoding submissions.
