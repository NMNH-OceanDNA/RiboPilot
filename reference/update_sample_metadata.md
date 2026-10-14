# Update metadata for existing samples

Updates metadata columns for samples already in a project's database
(adding new columns if needed). Never touches R1, R2, Assembly, or any
pipeline state or lock column, except that a sample whose Taxon changes
has its resolved taxonomy cleared, unlocked, and queued for the Taxonomy
step again. Use
[`add_samples()`](https://nmnh-oceandna.github.io/RiboPilot/reference/add_samples.md)
to change what a sample reads or assembles from. Runs the same
read_mapping()/check_mapping() checks as
[`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md).
Backs up the database to `.old_sqlite_dbs/` first.

## Usage

``` r
update_sample_metadata(
  project = ".",
  mapping_fn = NULL,
  mapping_id = "ID",
  mapping_taxon = "Taxon",
  data_path = NULL,
  assembly_path = NULL,
  check_reads = TRUE,
  check_assemblies = TRUE,
  mapping_geome = "GEOME_BCID",
  fetch_geome = TRUE,
  mapping_gbif = "GBIF_ID",
  fetch_gbif = TRUE,
  mapping_ncbi = "BioSample",
  fetch_ncbi = TRUE,
  link_sources = FALSE
)
```

## Arguments

- project:

  Path to the project directory (default = current working directory).

- mapping_fn:

  Path to a comma-separated mapping file. Same shape as
  [`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md)'s
  mapping file. Every `ID` must already be in the project.

- mapping_id:

  Column in the mapping file with unique ids (default "ID").

- mapping_taxon:

  Column in the mapping file with a taxonomic identifier (default
  "Taxon").

- data_path:

  Directory holding the read files named in the mapping file (checked,
  not stored).

- assembly_path:

  Directory holding the assembly files named in the mapping file
  (checked, not stored).

- check_reads:

  (logical) Check that the read files named in the mapping file exist
  under `data_path` (default TRUE).

- check_assemblies:

  (logical) Check that the assembly files named in the mapping file
  exist under `assembly_path` (default TRUE).

- mapping_geome:

  Name of the mapping-file column holding GEOME BCIDs (optional). Stored
  as `GEOME_BCID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_geome:

  Fetch GEOME metadata for samples with a BCID when the identifiers
  change (default TRUE). Set FALSE when offline and run
  [`fetch_geome()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_geome.md)
  later.

- mapping_gbif:

  Name of the mapping-file column holding GBIF occurrence IDs
  (optional). Stored as `GBIF_ID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_gbif:

  Fetch GBIF metadata for samples with a GBIF ID when the identifiers
  change (default TRUE). Set FALSE when offline and run
  [`fetch_gbif()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_gbif.md)
  later.

- mapping_ncbi:

  Name of the mapping-file column holding NCBI BioSample or SRA
  accessions (optional). Must not be the sample ID column. Stored as
  `BioSample`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_ncbi:

  Fetch NCBI metadata for samples with a BioSample value when the
  identifiers change (default TRUE). Set FALSE when offline and run
  [`fetch_ncbi()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_ncbi.md)
  later.

- link_sources:

  Follow links between GEOME, GBIF, and NCBI records for these samples
  (default FALSE). Linking also runs when the project setting is on;
  this argument does not change the setting. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
