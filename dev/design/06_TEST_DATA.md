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

## Status (2026-09-02): small path implemented, assets pending

`new_test_project(full_size = FALSE)` (the default) now fetches per-sample
`<acc>_R1/_R2.fastq.gz` from the GitHub release tagged `RP_TEST_READS_TAG`
(`R/init_test_project.R`), checksum-verified via `fetch_verified()`. The
subsets are built by `dev/make_test_reads.R`: reads are mapped with minimap2
to the sample's own assembled scaffolds (from a full-size project database)
plus the GetOrganelle seed FASTA, pairs with either mate mapped are kept, and
the pool is capped. Until the repo is public the release assets are not
reachable, so `full_size = TRUE` remains the working path.

## Earlier status (2026-08-10): NOT IMPLEMENTED

The shippable subset above was never built. Current `R/init_test_project.R`
always downloads the **full** read sets (several GB) and there is no
`full_size` argument. `inst/test_data/` holds only `mapping_test.csv`.

What the current implementation does have, and what should be kept when the
small path is added:

- ENA fetch with `curl -f` and a bounded retry count, then an NCBI SRA
  (`fasterq-dump`) fallback for networks where the ENA fastq mirror is blocked
  (e.g. Smithsonian). MitoPilot's version has neither, and loops forever on a
  bad URL.
- Already-present, gzip-valid read files are skipped, so an interrupted setup
  resumes cleanly.
- Sample list has drifted from the table above. Current `mapping_test.csv` is
  10 fish skims from study SRP313821: SRR22396758, SRR14433904, SRR14163830,
  SRR14433900, SRR14433901, SRR14433899, SRR14433895, SRR22396893,
  SRR25465041, SRR33849683. Decide whether to keep these or revert to the
  proposed 10 before generating subsets.

### Notes for whoever implements it

- **Assembly step is compatible.** `assemble.nf` runs GetOrganelle, which
  already does seed-based read recruitment internally, so a pre-baited read
  pool is just a smaller starting pool, not a different kind of input.
  `preprocess.nf` is plain `fastp` with no depth normalization that would care
  about missing background reads.
- **rDNA is high copy**, unlike the mitogenome MitoPilot baits on, so a naive
  bait will pull far more reads than MitoPilot's ~0.5-1 MB per mate. Expect to
  need a cap: subsample to a fixed read count, or bait against only the
  cistron rather than the whole repeat unit.
- **Package size.** MitoPilot ships 17 MB in `inst/test_data/`, well past the
  5 MB CRAN limit. If RiboPilot ever targets CRAN, host the filtered reads as
  a GitHub release asset and download those (~10 MB, seconds, not GB) instead
  of bundling them.
- API shape: add `full_size = FALSE` defaulting to the small path, and keep the
  existing ENA/SRA downloader as the `full_size = TRUE` path.
