# 00 - Overview

## Vision

RiboPilot is an R package + Nextflow pipeline + R-Shiny app for assembling,
annotating, and curating the **nuclear ribosomal repeat (nrDNA)** from genome
skimming data, across a broad taxonomic scope (metazoans first, then a path to
plants/fungi/protists). It is the nuclear-ribosomal counterpart to MitoPilot and
deliberately reuses MitoPilot's architecture so that a lab already running
MitoPilot can point RiboPilot at the same reads and map file.

## Why a separate tool

Genome skims (low-coverage shotgun, ~1-10x nuclear) already recover the
high-copy fraction of the genome for free: the plastid, the mitogenome, and the
**nuclear ribosomal cistron**, which exists in hundreds to thousands of tandem
copies and so is sequenced to very high effective depth. MitoPilot harvests the
mitogenome from these skims; the same libraries contain the full
`18S-ITS1-5.8S-ITS2-28S` array, which yields the most widely used nuclear
phylogenetic and barcoding markers (18S, 28S D1-D3, ITS1, ITS2, full ITS). A
dedicated tool lets us recover both marker sets from one sequencing run.

## Target locus

The eukaryotic 45S ribosomal cistron, transcribed by RNA Pol I as a single
pre-rRNA and processed into the mature rRNAs:

| Element | Type | Notes |
|---|---|---|
| 5' ETS | external transcribed spacer | variable, often not recovered/curated |
| **18S** (SSU) | rRNA gene | ~1.8 kb, highly conserved, deep phylogeny |
| **ITS1** | internal transcribed spacer | fast-evolving, species-level barcode |
| **5.8S** | rRNA gene | ~160 bp, conserved anchor |
| **ITS2** | internal transcribed spacer | fast-evolving, species-level barcode |
| **28S** (LSU) | rRNA gene | ~3.4-5+ kb, conserved core + variable D-domains (D1-D12) |
| 3' ETS / IGS / NTS | spacers | separate copies; usually excluded |

5S rRNA is typically a **separate**, independently transcribed array in most
eukaryotes and is **out of primary scope** (optional stretch goal).

## Core challenge vs MitoPilot

MitoPilot curates protein-coding genes and structural RNAs whose boundaries are
pinned by start/stop codons, reading frames, and tRNA secondary structure.
The nrDNA repeat is:

- **Non-protein-coding** - no codons, no ORFs, no translation check. Boundary
  cues are covariation/secondary structure and alignment to reference profiles.
- **More variable in the spacers** - ITS1/ITS2 and the 28S D-domains vary in
  length and sequence even within genera; expansion segments can differ by
  kilobases across phyla. Reference-guided curation must lean on the conserved
  rRNA cores (18S, 5.8S, 28S core) as anchors and treat spacers as
  "trim between anchors" rather than "align base-for-base."
- **Boundary-defined, not feature-defined** - the interesting products (ITS1,
  ITS2) are literally the gaps between conserved genes, so the annotation
  problem is largely a **boundary-calling** problem. This actually plays to
  profile-HMM tools (ITSx, barrnap) that are built for exactly this.

The curation philosophy therefore shifts from MitoPilot's "validate each gene
model" toward "**place the conserved anchors precisely, derive spacers by
subtraction, flag anything that disagrees with the reference profile.**"

## Design principles carried over from MitoPilot

1. Same input map file (`ID, Taxon, R1, R2` + arbitrary metadata columns).
2. Same project init pattern (`new_project()` -> `.sqlite` DB + `.config`).
3. SQLite database as the single source of truth for state and curation edits.
4. Nextflow for execution; **reuse MitoPilot's executor/config system verbatim**
   (local, slurm, sge, pbs, lsf, awsbatch, NMNH_Hydra, NOAA_SEDNA, saved profiles).
5. R-Shiny app for review/curation with the same three-stage flow
   (assemble -> annotate -> export), adapted to rDNA.
6. Reference-guided curation with per-clade rulesets (thin wrappers over a
   shared core), mirroring MitoPilot's `curate_*` design.
7. GenBank-ready export (features + source modifiers + fasta/tbl/sqn inputs).

See [`01_ARCHITECTURE.md`](01_ARCHITECTURE.md) for how each maps onto code.
