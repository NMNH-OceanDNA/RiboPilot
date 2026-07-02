# 04 - Existing tools, prior art, reference databases

Survey from web research (July 2026). Nothing does the full assemble +
annotate + reference-curate + GenBank-export loop with a review GUI for broad
taxa, which is the RiboPilot gap. But strong components exist to wrap.

## Closest prior-art pipeline

### skim2rrna (o-william-white) - the key precedent

Snakemake pipeline for **batch assembly, annotation, and phylogenetics of
nuclear ribosomal genes from low-coverage genome skims** (museum specimens).
Part of a toolkit with `skim2mito` and `gene2phylo`. Reported ~3 kb rRNA
recovery from 27/27 test samples.

- Assembly: **GetOrganelle**.
- Homology filter: **BLAST vs a SILVA 138-derived DB**.
- Annotation: **barrnap**.
- Published: White et al. 2025, Mol Ecol Resour, doi:10.1111/1755-0998.14036;
  repo: github.com/o-william-white/skim2rrna.

RiboPilot is essentially "skim2rrna's scientific approach, rebuilt on MitoPilot's
architecture (Nextflow + SQLite + Shiny curation + GenBank export + per-clade
reference curation), for broad taxa with a human-in-the-loop GUI." Study this
repo closely before implementing.

## Assembly tools

| Tool | Role | Notes |
|---|---|---|
| **GetOrganelle** | targeted nr assembly | `-F embplant_nr` / `fungus_nr` built in; animals need custom seed/label DBs. Already used by MitoPilot -> plumbing exists. |
| SPAdes / metaSPAdes / MEGAHIT | full de novo (Option A) | contig-search fallback + cross-check |
| rnaSPAdes | region assembly | validate repeat unit / expansion segments |
| aTRAM, MITObim | iterative baiting (Option B2) | reference-guided extension for divergent taxa |
| PhyloHerb | herbarium organelle+nr extraction | wraps GetOrganelle; reference for nr extraction UX |

## Annotation tools

| Tool | Role | Notes |
|---|---|---|
| **barrnap / pyBarrnap** | 18S/5.8S/28S (+5S) genes | nHMMER HMMs, `--kingdom euk`; `--accurate` sharpens 5.8S/ITS2/28S edges; GFF3 out |
| **ITSx** | ITS1 / 5.8S / ITS2 extraction | HMMER3, ~20 eukaryote groups, RC-aware, broad taxa. Primary spacer authority |
| ITSxpress | fast ITS trimming | amplicon-oriented; ITSx more general for assembled units |
| RNAmmer | SSU/LSU prediction | older, HMMER2; barrnap is the modern replacement |
| Infernal cmsearch + Rfam | structure-aware anchors | RF01960 (SSU euk), RF02543 (LSU euk), RF00002 (5.8S); optional high-accuracy mode |
| Metaxa2 | SSU/LSU classification | more for taxonomic sorting than boundary calling |
| DECIPHER | reference alignment/curation | already a MitoPilot dep; use for anchor snapping + identity track |

## Reference databases

| DB | Content | Use in RiboPilot |
|---|---|---|
| **SILVA** (138+) | curated SSU (18S) + LSU (28S) rRNA | BLAST/seed DB for assembly + anchor references; skim2rrna uses this |
| **EUKARYOME** | full-length 18S-ITS-28S for all eukaryotes | best single source for **whole-cistron** references across broad taxa |
| **PR2** | protist-focused SSU (+LSU) | fill protist gaps |
| ITSoneDB | curated ITS1 | ITS1 boundary/QC reference |
| UNITE | fungal ITS | fungal ITS reference (if fungi in scope) |
| Rfam | rRNA covariance models | cmsearch anchors |
| NCBI nt / GenBank | catch-all | candidate-ref BLAST + contaminant lineage check (reuse MitoPilot machinery) |

**EUKARYOME** is the standout: it provides long-read/full-length
18S-ITS-28S references, which is exactly what reference-guided curation of the
whole cistron needs. Plan to build the packaged reference sets primarily from
EUKARYOME (whole cistron) + SILVA (anchor seeds).

## Reference cross-checks / recipes

- Smithsonian NMNH Ocean-DNA `BLAST_18s_28s` recipe (assemble reads, BLAST
  contigs vs 18S/28S) - the Option A baseline the user cited.
- "Skimming for barcodes" (Trevisan et al., PeerJ 2022) - recovering mito +
  nrDNA repeat reference markers from shallow shotgun; motivates the dual
  mito+rDNA-from-one-skim workflow.
- "Decode the complete 18-5.8-28S repeated units by genome skimming" (green
  algae) - confirms full-unit recovery is feasible from skims.

## Sources

- skim2rrna: https://github.com/o-william-white/skim2rrna ; White et al. 2025 https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.14036 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11646300/
- GetOrganelle: https://github.com/Kinggerm/GetOrganelle
- PhyloHerb: https://github.com/lmcai/PhyloHerb
- barrnap: https://github.com/tseemann/barrnap ; pyBarrnap notes https://scienceparkstudygroup.github.io/ibed-bioinformatics-page/source/metagenomics/barrnap.html
- ITSx: https://microbiology.se/software/itsx/ ; Bengtsson-Palme et al. 2013 https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12073
- EUKARYOME: https://pmc.ncbi.nlm.nih.gov/articles/PMC11168333/
- SILVA: https://www.arb-silva.de/
- Trevisan et al. 2022 (skimming for barcodes): https://peerj.com/articles/13790/
- Green algae rRNA repeat by skimming: https://pmc.ncbi.nlm.nih.gov/articles/PMC5713310/
- Smithsonian Ocean-DNA BLAST_18s_28s: https://smithsonian.github.io/NMNH-Ocean-DNA/qmd/hydraScripts/BLAST_18s_28s.html
