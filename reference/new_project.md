# Initialize a new RiboPilot project

Mirrors MitoPilot's `new_project()`. Creates the project directory,
copies the mapping file, initializes the SQLite tracking database, and
writes a `.config` Nextflow config from the chosen executor template.

## Usage

``` r
new_project(
  path = ".",
  mapping_fn = NULL,
  mapping_id = "ID",
  data_path = NULL,
  assembly_path = NULL,
  min_depth = 2e+06,
  executor = c("docker", "conda", "local", "slurm", "pbs", "sge", "lsf", "NMNH_Hydra",
    "NOAA_SEDNA"),
  container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
  custom_seeds_db = NULL,
  custom_labels_db = NULL,
  ref_dir = NULL,
  config = NULL,
  profile_dir = ribopilot_config_dir(),
  ncbi_api_key = NULL,
  check_reads = TRUE,
  check_assemblies = TRUE,
  mapping_geome = "GEOME_BCID",
  fetch_geome = TRUE,
  mapping_gbif = "GBIF_ID",
  fetch_gbif = TRUE,
  mapping_ncbi = "BioSample",
  fetch_ncbi = TRUE,
  link_sources = FALSE,
  Rproj = TRUE,
  force = FALSE,
  ...
)
```

## Arguments

- path:

  Path to the project directory (default = current working dir).

- mapping_fn:

  Path to a comma-separated mapping file. Minimally an `ID` column
  (unique sample id), a `Taxon` column, and `R1`/`R2` columns naming the
  paired read files relative to `data_path` (or an `Assembly` column).
  Column names are matched ignoring case. Extra metadata columns are
  preserved.

- mapping_id:

  Column in the mapping file with unique ids (default "ID").

- data_path:

  Directory where the raw read data is located. Required when any sample
  lists reads; not needed when every sample supplies its own assembly.

- assembly_path:

  Directory holding the assembly files named in the mapping file's
  optional `Assembly` column. A sample with an assembly is searched for
  its nuclear rDNA contigs instead of being assembled from reads, and
  its `R1`/`R2` may be left blank. See dev/design/20_USER_ASSEMBLIES.md.

- min_depth:

  Minimum number of reads left after trimming, counting both reads of a
  pair (2000000 is 1000000 pairs). A sample below it is not assembled
  and is marked Failed with a note.

- executor:

  Executor for the Nextflow pipeline. Ready to use as-is: "docker"
  (default) and "local" run the RiboPilot container on this machine, and
  "NMNH_Hydra" and "NOAA_SEDNA" are ready-made Singularity profiles for
  those clusters. "conda" uses a self-installed conda env (not
  supported; warns). For any other cluster, build a profile with
  generate_config() and pass its name here. A bare scheduler name
  ("slurm", "sge", "pbs", "lsf") is an error: it has no queue or cluster
  options yet.

- container:

  Container to use for pipeline execution.

- custom_seeds_db:

  Path/URL to a custom GetOrganelle nr seeds database.

- custom_labels_db:

  Path/URL to a custom GetOrganelle nr labels database.

- ref_dir:

  Folder for built reference databases (default: the `RIBOPILOT_REF_DIR`
  environment variable, else
  [`ribopilot_ref_dir()`](https://nmnh-oceandna.github.io/RiboPilot/reference/ribopilot_ref_dir.md)).

- config:

  Optional path to an existing custom Nextflow config file.

- profile_dir:

  Directory searched for saved cluster profiles.

- ncbi_api_key:

  Optional NCBI API key (raises request rate limits).

- check_reads:

  (logical) Check that the read files named in the mapping file exist
  under `data_path` (default TRUE). Set FALSE when the data will be
  staged later.

- check_assemblies:

  (logical) Check that the assembly files named in the mapping file
  exist under `assembly_path` (default TRUE).

- mapping_geome:

  Name of the mapping-file column holding GEOME BCIDs (optional). Stored
  as `GEOME_BCID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
  Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- fetch_geome:

  Fetch GEOME metadata for samples with a BCID during setup (default
  TRUE). Set FALSE when offline and run
  [`fetch_geome()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_geome.md)
  later. Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- mapping_gbif:

  Name of the mapping-file column holding GBIF occurrence IDs
  (optional). Stored as `GBIF_ID`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
  Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- fetch_gbif:

  Fetch GBIF metadata for samples with a GBIF ID during setup (default
  TRUE). Set FALSE when offline and run
  [`fetch_gbif()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_gbif.md)
  later. Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- mapping_ncbi:

  Name of the mapping-file column holding NCBI BioSample or SRA
  accessions (optional). Must not be the sample ID column. Stored as
  `BioSample`. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
  Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- fetch_ncbi:

  Fetch NCBI metadata for samples with a BioSample value during setup
  (default TRUE). Set FALSE when offline and run
  [`fetch_ncbi()`](https://nmnh-oceandna.github.io/RiboPilot/reference/fetch_ncbi.md)
  later. Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- link_sources:

  Follow links between GEOME, GBIF, and NCBI records to fill in IDs a
  sample does not have yet (default FALSE). Saved as the project
  setting. See
  [`vignette("Sample-Metadata")`](https://nmnh-oceandna.github.io/RiboPilot/articles/Sample-Metadata.md).
  Passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).

- Rproj:

  (logical) Initialize/open an RStudio project (default TRUE; no effect
  outside RStudio).

- force:

  (logical) Recreate an existing project database (default FALSE).
  Without it, `new_project()` stops when the folder already holds a
  project.

- ...:

  Additional default processing parameters passed to
  [`new_db()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_db.md).
