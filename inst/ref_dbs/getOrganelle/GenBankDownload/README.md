# nrDNA reference database builder

Builds GetOrganelle nuclear ribosomal (nrDNA) seed + label databases for a
clade by downloading 18S / 5.8S / 28S / ITS records from GenBank. Adapted from
[GenBankMitoReview](https://github.com/dmacguigan/GenBankMitoReview), retargeted
from mitochondrial genomes to the nuclear ribosomal repeat unit.

## Files

- `download_rrna.sh` - query GenBank, download nuclear rRNA records (fasta),
  build seed + label FASTAs.
- `parse_rrna_fasta.py` - classify records (18S / 28S / 5.8S / multi-element
  rDNA repeat unit) by title and write GetOrganelle-format FASTAs.
- `bin/curl` - HTTP/1.1 shim for Entrez Direct (fixes the OpenSSL 3.x
  "curl (56) unexpected eof" failure).

## Usage

### From R (recommended)

`build_reference_db()` wraps this script and writes the seed/label DBs to a
reusable directory outside any project (default `ribopilot_ref_dir()`), so one
clade DB serves every project for that clade:

```r
# count first (no download)
build_reference_db('"Gastropoda"[Organism]', "gastropoda",
                   count_only = TRUE, conda_env = "ribopilot")
# build once, reuse everywhere
db <- build_reference_db('"Gastropoda"[Organism]', "gastropoda",
                         out_dir = "~/ribopilot_refs", conda_env = "ribopilot")
new_project(..., custom_seeds_db = db$seeds, custom_labels_db = db$labels)
```

Re-running with the same `clade`/`out_dir` reuses the existing DB unless
`overwrite = TRUE`.

### Directly (bash)

Run in a conda env with Entrez Direct + biopython (e.g. `mitoreview` or the
RiboPilot env):

```bash
# count first
bash download_rrna.sh --count-only '"Actinopterygii"[Organism]' actinopterygii_rrna
# build
bash download_rrna.sh --outdir build --api-key <KEY> \
  '"Actinopterygii"[Organism]' actinopterygii_rrna
```

Outputs `<prefix>_seeds.fasta` (read-recruiting: 18S/28S/whole rDNA repeat,
`>ACCESSION_Genus_species`) and `<prefix>_labels.fasta` (per-gene: `>18S`/`>28S`,
for contig classification). Place these under `ref_dbs/getOrganelle/seeds/` and
`ref_dbs/getOrganelle/labels/` and point `new_project(custom_seeds_db=,
custom_labels_db=)` at them (or the raw GitHub URLs, which are the `new_db()`
defaults).

## Built databases

- **Actinopterygii** (broad ray-finned fish): 11,753 GenBank records ->
  3,857 seeds (~1,700 species), 1,201 x 18S + 1,030 x 28S labels. Default for
  `new_project()`. Built uncapped (no length-based subsampling).
  Note: standalone 5.8S records are rare (5.8S is embedded in feature-less ITS
  records, which become multi-element seeds), so the label set has no 5.8S; the
  5.8S region is still covered by the rDNA repeat unit seeds.

## Adding a clade

Run `download_rrna.sh` with a different `[Organism]` query (or BioProject, etc),
commit the resulting seed/label FASTAs, and pass them to `new_project()`.
