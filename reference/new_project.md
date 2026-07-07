# Initialize a new RiboPilot project

Mirrors MitoPilot's \`new_project()\`. Creates the project directory,
copies the mapping file, initializes the SQLite tracking database, and
writes a \`.config\` Nextflow config from the chosen executor template.

## Usage

``` r
new_project(
  path = ".",
  mapping_fn = NULL,
  mapping_id = "ID",
  data_path = NULL,
  min_depth = 2e+06,
  executor = c("conda", "docker", "local", "awsbatch", "slurm", "pbs", "sge", "lsf",
    "NMNH_Hydra", "NOAA_SEDNA"),
  container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
  custom_seeds_db = NULL,
  custom_labels_db = NULL,
  config = NULL,
  profile_dir = ribopilot_config_dir(),
  ncbi_api_key = NULL,
  Rproj = TRUE,
  force = FALSE,
  ...
)
```

## Arguments

- path:

  Path to the project directory (default = current working dir).

- mapping_fn:

  Path to a mapping CSV. Minimally an \`ID\` column (unique sample id),
  a \`Taxon\` column, and \`R1\`/\`R2\` columns naming the paired read
  inputs. Extra metadata columns are preserved.

- mapping_id:

  Column in the mapping file with unique ids (default "ID").

- data_path:

  Directory where the raw read data is located.

- min_depth:

  Minimum trimmed paired reads to proceed with assembly.

- executor:

  Executor for the Nextflow pipeline: a built-in template ("conda"
  (default) for a local conda env, "docker"/"local" for the RiboPilot
  container, "awsbatch", "slurm", "sge", "pbs", "lsf", "NMNH_Hydra",
  "NOAA_SEDNA") or a saved profile from \[generate_config()\].

- container:

  Container to use for pipeline execution.

- custom_seeds_db:

  Path/URL to a custom GetOrganelle nr seeds database.

- custom_labels_db:

  Path/URL to a custom GetOrganelle nr labels database.

- config:

  Optional path to an existing custom Nextflow config file.

- profile_dir:

  Directory searched for saved cluster profiles.

- ncbi_api_key:

  Optional NCBI API key (raises request rate limits).

- Rproj:

  (logical) Initialize/open an RStudio project (default TRUE).

- force:

  (logical) Recreate an existing project database (default FALSE).

- ...:

  Additional default processing parameters passed to \[new_db()\].
