# 06 - Test dataset

Ship a small test dataset with the package and provide `new_test_project()`
exactly like MitoPilot (`init_test_project.R`): a `mapping_test.csv` + shippable
read files under `inst/test_data/`, and a `full_size = TRUE` option that fetches
full reads on demand from the ENA FTP mirror.

## Samples

Use the **same test fishes as MitoPilot** (from `inst/test_data/mapping_test.csv`).
MitoPilot ships 13; RiboPilot targets ~10. Proposed 10 (taxonomically spread
across the 13; drop 3 of the most redundant eels/deep-sea forms):

| ID (SRA) | Taxon |
|---|---|
| SRR22396758 | Upeneus parvus |
| SRR22396794 | Psychrolutes paradoxus |
| SRR22396940 | Psenes pellucidus |
| SRR22396740 | Hoplostethus occidentalis |
| SRR21844202 | Fundulus majalis |
| SRR22396640 | Xyrichtys novacula |
| SRR22396627 | Gigantura indica |
| SRR21843972 | Stomias affinis |
| SRR22396843 | Conger oceanicus |
| SRR19434536 | Rhinecanthus rectangulus |

(Exact 10 easy to adjust; keep the same `ID`/`Taxon`/`R1`/`R2` CSV schema.)

## Critical difference from MitoPilot's shipped files

MitoPilot's shipped small files are **NOT** raw reads. `filter_reads.sh` produces
them by: fastp trim -> **bowtie2 map to the mitogenome** -> keep only reads that
map -> `seqtk subseq` the raw reads to those names. So the shipped files contain
essentially **only mitochondrial reads and no nuclear rDNA**. They are useless for
RiboPilot.

Therefore RiboPilot must start from the **full read sets** (`full_size = TRUE`
equivalent: full fastqs pulled from `ftp.sra.ebi.ac.uk/vol1/fastq/...`, same URL
scheme MitoPilot uses), then build its own rDNA-filtered shippable subset.

## Building the shippable subset (mirror `filter_reads.sh`, retargeted to rDNA)

An rDNA analogue of MitoPilot's `filter_reads.sh`:

1. Download full raw reads for each sample (ENA FTP, as in MitoPilot `full_size`).
2. `fastp` adapter/quality trim.
3. **bowtie2 map to an rDNA reference** (an 18S-ITS1-5.8S-ITS2-28S reference; a
   teleost rDNA, since all samples are fish) instead of a mitogenome index.
4. Keep mapped read names; `seqtk subseq` raw reads to those names.
5. Result: small per-sample fastqs containing the rDNA-region reads, safe to ship
   in `inst/test_data/`.

This keeps the package small while guaranteeing the shipped test data actually
contains the target locus. `new_test_project(full_size = TRUE)` still downloads
the complete reads for realistic full runs, identical to MitoPilot.

## Generation is a build task, not shipped-code

Generating the subset requires: (a) multi-GB downloads from ENA, (b) an rDNA
reference to map against, (c) bowtie2/samtools/seqtk/fastp. It should be run once
during Phase 0/1 (needs the reference DB from Phase 1 to exist), the outputs
committed to `inst/test_data/`, and the generator script kept in the repo
(`inst/test_data/filter_reads_rrna.sh`) for reproducibility, exactly as MitoPilot
keeps `filter_reads.sh`.

## Open questions

- Which exact 10 samples (or keep all 13?).
- Which reference rDNA sequence to filter against, and whether to filter loosely (keep
  flanking IGS reads) so assembly has enough overhang to resolve boundaries.
- Target shipped size budget per sample (MitoPilot's mito subsets are ~0.5-1 MB
  each; rDNA subsets will be similar or a bit larger given higher copy number).
