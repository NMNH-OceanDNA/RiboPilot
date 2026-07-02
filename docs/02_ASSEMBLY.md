# 02 - Assembly strategies

The nrDNA repeat is high-copy, so the signal is strong even in shallow skims.
The engineering question is how to pull the `18S-ITS1-5.8S-ITS2-28S` unit out
cleanly without dragging in paralogous/pseudogene copies, intragenomic variants,
or contaminant rDNA (a real risk: fungal/algal/microbial rDNA in the same
library).

Two broad families, in increasing order of preference.

## Option A - Full de novo assembly + contig search (baseline)

Assemble all reads (or a downsampled subset) with a general assembler
(SPAdes/metaSPAdes/MEGAHIT), then find rDNA contigs by homology search against a
reference rRNA database. This mirrors the Smithsonian Ocean-DNA
`BLAST_18s_28s` recipe:

1. `fastp` trim (shared with MitoPilot preprocess).
2. Assemble reads.
3. `blastn` contigs vs an 18S/28S reference DB (SILVA / EUKARYOME / PR2).
4. Keep best-scoring contigs, stitch/extract the repeat region.

**Pros:** simple, robust, no seed tuning, catches unexpected structure.
**Cons:** wasteful (assembles the whole skim), high-copy repeat can fragment or
collapse, harder to get one clean canonical unit, slower. Good as a **fallback
and as a validation cross-check**, not the default.

## Option B - Targeted assembly (preferred default)

Seed the assembly from conserved rRNA anchors and walk out through the spacers,
so only reads belonging to the rDNA repeat unit are used. This is exactly how MitoPilot
uses GetOrganelle for the mitogenome, so it slots into the existing pattern.

### B1. GetOrganelle nr mode (primary)

GetOrganelle has first-class nuclear ribosomal targets: `-F embplant_nr` (plants)
and `-F fungus_nr` (fungi), using a seed DB (conserved rRNA to recruit reads) and
a label DB (to classify assembly-graph contigs). For **animals/broad metazoa**
there is no built-in `nr` preset, so we supply **custom seed + label DBs** via the
same `custom_seeds_db` / `custom_labels_db` arguments MitoPilot already exposes.

Plan:
- **Seed/label DBs built per-clade from a GenBank query, mirroring MitoPilot.**
  MitoPilot ships `ref_dbs/getOrganelle/GenBankDownload/GenBankDownloadUtil.sh`:
  the user runs it once with a taxonomic advanced-search term (e.g.
  `'"Percidae"[Organism]'`), and it `esearch`/`efetch`es all matching records,
  filters, and emits GetOrganelle seed + label FASTAs. RiboPilot mirrors this
  exactly with an rDNA-flavored downloader (see below). This is preferred over
  per-sample auto-fetch (`go_fetch`): the user just grabs everything for a large
  clade once, no reliance on a rank-expansion heuristic, and nrDNA records are
  small so it is cheap in compute and storage. One clade DB serves every sample
  in that clade.
- Run GetOrganelle from reads with those DBs to assemble the repeat unit.

#### Implemented builder (Increment 1)

`ref_dbs/getOrganelle/GenBankDownload/` contains the working builder, adapted
from [GenBankMitoReview](https://github.com/dmacguigan/GenBankMitoReview):
`download_rrna.sh` (Entrez Direct query + fasta download + the HTTP/1.1 curl
shim) and `parse_rrna_fasta.py` (classify records by title into 18S/28S/5.8S/
multi and emit GetOrganelle seed + label FASTAs). A broad **Actinopterygii**
database is built and shipped (3,857 seeds; 1,201 x 18S + 1,030 x 28S labels)
and is the `new_project()` default. See that folder's README.

#### Validation (Increment 2)

Validated end-to-end on a full test sample (`SRR22396758`, *Upeneus parvus*,
4.3M read pairs) with the shipped Actinopterygii DBs:

```
fastp -> get_organelle_from_reads.py -1 R1 -2 R2 -s <seeds> --genes <labels> \
        -F anonym -R 10 -k 21,45,65,85,105,115 --larger-auto-ws -t 8
```

Result: 2 candidate paths, 3 scaffolds each, ~9.77 kb total per path. BLAST vs
the label DB confirms it is genuine rDNA: a full **18S** (98.6% identity to
reference) plus **28S** scaffolds (92-96%). This validates both the assembly
approach and the Increment 1 reference DBs. The two paths are alternative graph
traversals (candidate intragenomic variants) -> emit as separate paths, as
planned.

#### Design of the per-clade nr reference DB (mirror of MitoPilot's util)

A `GenBankDownloadUtil` script (bash + a `parseGB.py`, same shape as MitoPilot's)
that:
1. Takes a taxonomic query term (`'"Gastropoda"[Organism]'`, a BioProject, etc).
2. `esearch | efetch` all nucleotide records for that clade whose products are
   nuclear rRNA / spacers. Filter by product/feature name to
   `18S | 28S | 5.8S | internal transcribed spacer | ITS | small/large subunit
   ribosomal RNA`. These names are **nuclear-specific**, which cleanly excludes
   organellar rRNA (mito 12S/16S = `rrnS`/`rrnL`), so no organelle location
   filter is needed (the inverse of MitoPilot's `-location mitochondrion`).
3. Split records like MitoPilot does:
   - **seeds** = long/multi-element records (whole rDNA repeat or 18S+ITS+28S) used to
     recruit reads,
   - **labels** = single-locus records with the product name in the header
     (`18S ribosomal RNA`, `28S ribosomal RNA`, `5.8S ribosomal RNA`, ITS), used
     to classify assembly-graph contigs,
   - dedup exact sequences.
4. Write `seeds/<clade>_nr_seeds.fasta` + `labels/<clade>_nr_labels.fasta`,
   passed to GetOrganelle via the existing `custom_seeds_db` / `custom_labels_db`
   arguments.

Ship a few common built clade DBs in `ref_dbs/getOrganelle/` (as MitoPilot ships
`fish_mito_*`), and let users build their own for any clade with one command.
Optionally still accept per-sample seed/gene overrides in the map file, and keep
`go_fetch` as an alternative auto-fetch path for one-off odd samples (not the
default).
- Resolve the assembly graph to the monomer repeat unit; report alternative paths
  (intragenomic variants) as separate paths, like MitoPilot's multi-path output.
- Orient and rotate to a canonical start (5' end of 18S) so all samples align.

**Pros:** targeted, fast, reuses MitoPilot's assembler + custom-DB plumbing +
path-scoring UI; produces a clean circular-graph-derived monomer.
**Cons:** needs good seed/label DBs for broad taxa; may need k-mer tuning for
very divergent lineages.

### B2. Reference-guided / iterative baiting (robust for divergent taxa)

When de-novo seeding struggles (very divergent taxon, low depth, contamination),
map/bait reads to a reference rDNA sequence and iteratively extend:

- Bait reads by k-mer/BLAST match to a reference 18S+5.8S+28S (conserved anchors
  only, so it works across large phylogenetic distance).
- Assemble the baited read pool (SPAdes/rnaSPAdes) or do reference-guided
  consensus.
- Iterate (MITObim/aTRAM-style) to fill spacers between anchors.

Candidate engines: `GetOrganelle` (already iterative), `aTRAM`, `MITObim`,
`RAMBLER`/reference-guided consensus, or a simple bait->assemble->extend loop we
control in Nextflow. Prior art `skim2rrna` uses **GetOrganelle for assembly + a
BLAST homology filter + barrnap for annotation** on exactly this data type; it is
the closest existing pipeline and a strong template.

### B3. Region assemblers as cross-check

`rnaSPAdes` / `rnaviralSPAdes` or metaSPAdes on baited reads can validate the
GetOrganelle unit and help resolve expansion segments in 28S.

## Handling multi-copy reality

nrDNA is not perfectly homogenized (incomplete concerted evolution). Expect:
- Multiple near-identical variants -> emit as separate paths (reuse MitoPilot
  path model + scoring); let the user pick the representative in the app.
- Length-variant ITS/28S expansion segments -> keep the dominant-coverage path
  as primary, flag minor variants.
- Contaminant rDNA -> two layers: (a) BLAST-lineage check on each path (reuse
  MitoPilot's candidate-BLAST + lineage columns), and (b) optional
  **blobtools + minimap2** contamination screen as skim2rrna does (map reads back
  to each contig, combine coverage + BLAST taxonomy to spot foreign rDNA). Flag
  paths whose top-hit lineage or blob taxonomy disagrees with `Taxon`.

## Recommended default

**B1 (GetOrganelle nr with custom broad-taxa seed/label DBs) as primary**, with
**B2 reference-guided baiting as automatic fallback** when B1 fails a QC gate
(no recovered 18S+28S, or length outside clade tolerance), and **Option A**
available as a manual whole-assembly cross-check. This matches MitoPilot's
"GetOrganelle-first, alternate assembler fallback" philosophy.

## Open assembly questions

- Exact GenBank product/feature filter that keeps nuclear 18S/5.8S/28S/ITS and
  reliably drops mito 12S/16S and any mislabeled records. Validate on a clade with
  known messy annotations.
- Whether clade seeds should be whole-repeat records only, or also bare 18S/28S,
  to recruit reads without over-recruiting divergent-spacer or contaminant reads.
- Whether to circularize (repeat unit is tandem, so the graph is a cycle through
  the IGS) or emit a linear monomer cut at the 18S start. Leaning linear monomer,
  IGS excluded by default.
- Depth normalization: rDNA can be extremely deep; downsample to avoid
  assembler blow-up (GetOrganelle handles this, but cap `-R`/reads).
