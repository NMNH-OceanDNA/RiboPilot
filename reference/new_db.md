# Initialize a new RiboPilot project database

Creates the SQLite tracking database. Mirrors MitoPilot's `new_db()`
spine (samples / preprocess / pre_opts / assemble / assemble_opts /
blast_opts / assemblies) but with nuclear ribosomal (nrDNA) assembly
options. Annotation, curation, and export tables are added in later
build increments.

## Usage

``` r
new_db(
  db_path = "./.sqlite",
  mapping_fn = NULL,
  mapping_id = "ID",
  mapping_taxon = "Taxon",
  data_path = NULL,
  assembly_path = NULL,
  check_reads = TRUE,
  check_assemblies = TRUE,
  assemble_cpus = 6,
  assemble_memory = 24,
  seeds_db = NULL,
  labels_db = NULL,
  getOrganelle = paste("-F anonym", "-R 10 -k 21,45,65,85,105,115", "--larger-auto-ws"),
  max_paths = 10,
  max_scaffolds = 10,
  min_assembly_length = 500,
  ref_dir = NULL,
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

- db_path:

  Path to the new database file.

- mapping_fn:

  Path to the mapping CSV (columns ID, Taxon, R1, R2).

- mapping_id:

  Column name to use as the primary key (default "ID").

- mapping_taxon:

  Column name with a taxonomic identifier (default "Taxon").

- data_path:

  Directory holding the raw read files named in the mapping file. Read
  files are resolved as `<data_path>/<R1>`, the same way the Nextflow
  pipeline does it. When NULL, read files cannot be located and only the
  paired-end check is applied.

- assembly_path:

  Directory holding the assembly files named in the mapping file's
  optional `Assembly` column. A sample with an assembly skips read
  preprocessing and assembly: its nuclear rDNA contigs are found inside
  the supplied file instead. Its `R1`/`R2` may be blank.

- check_reads:

  (logical) Check that the read files named in the mapping file exist
  (default TRUE). Set FALSE when the data will be staged later.

- check_assemblies:

  (logical) Check that the assembly files named in the mapping file
  exist under `assembly_path` (default TRUE).

- assemble_cpus:

  Default number of CPUs for assembly.

- assemble_memory:

  Default memory (GB) for assembly.

- seeds_db:

  Path or URL to the GetOrganelle seeds database. NULL (default) names
  an Actinopterygii database that RiboPilot does not ship; the Reference
  step replaces it.

- labels_db:

  Path or URL to the GetOrganelle labels database. NULL (default) names
  an Actinopterygii database that RiboPilot does not ship; the Reference
  step replaces it.

- getOrganelle:

  Default GetOrganelle command-line options.

- max_paths:

  Max assembly paths for a sample to continue (default 10).

- max_scaffolds:

  Max contigs within a path to continue (default 10).

- min_assembly_length:

  Minimum contig length (bp) to retain (default 500).

- ref_dir:

  Folder for built reference databases (default: the `RIBOPILOT_REF_DIR`
  environment variable, else
  [`ribopilot_ref_dir()`](https://nmnh-oceandna.github.io/RiboPilot/reference/ribopilot_ref_dir.md)).

- mapping_geome:

  Name of the mapping-file column holding GEOME BCIDs (optional). Stored
  as `GEOME_BCID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_geome:

  Fetch GEOME metadata for samples with a BCID during setup (default
  TRUE). Set FALSE when offline and run
  [`fetch_geome()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_geome.md)
  later.

- mapping_gbif:

  Name of the mapping-file column holding GBIF occurrence IDs
  (optional). Stored as `GBIF_ID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_gbif:

  Fetch GBIF metadata for samples with a GBIF ID during setup (default
  TRUE). Set FALSE when offline and run
  [`fetch_gbif()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_gbif.md)
  later.

- mapping_ncbi:

  Name of the mapping-file column holding NCBI BioSample or SRA
  accessions (optional). Must not be the sample ID column. Stored as
  `BioSample`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).

- fetch_ncbi:

  Fetch NCBI metadata for samples with a BioSample value during setup
  (default TRUE). Set FALSE when offline and run
  [`fetch_ncbi()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_ncbi.md)
  later.

- link_sources:

  Follow links between GEOME, GBIF, and NCBI records to fill in IDs a
  sample does not have yet (default FALSE). Saved as the project
  setting. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
