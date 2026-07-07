# Generate and save a reusable cluster config profile

Configure your HPC cluster once and reuse it for every project. Mirrors
MitoPilot's \`generate_config()\`. Builds a Nextflow config from a
generic scheduler template, fills in cluster-level settings, and saves
it as a named profile in \[ribopilot_config_dir()\]. Afterwards
\`new_project(executor = "\<name\>")\` will find and use it
automatically.

## Usage

``` r
generate_config(
  name,
  scheduler = c("slurm", "sge", "pbs", "lsf", "local", "awsbatch"),
  container_engine = c("auto", "singularity", "apptainer", "docker"),
  container_cache = NULL,
  container_run_options = NULL,
  queue = NULL,
  account = NULL,
  cluster_options = NULL,
  penv = "mthread",
  profile_dir = ribopilot_config_dir(),
  interactive = FALSE,
  overwrite = FALSE
)
```

## Arguments

- name:

  Profile name. Saved as \`config.\<name\>\`; pass as \`executor\`.

- scheduler:

  Base template. One of "slurm", "sge", "pbs", "lsf", "local", or
  "awsbatch".

- container_engine:

  "auto" picks docker for local/awsbatch and singularity for HPC
  schedulers; or "singularity"/"apptainer"/"docker".

- container_cache:

  Optional cacheDir for singularity/apptainer.

- container_run_options:

  Optional runOptions for singularity/apptainer.

- queue:

  Partition / queue name. If \`NULL\`, the directive is omitted.

- account:

  Optional accounting string, folded into clusterOptions.

- cluster_options:

  Optional raw clusterOptions string appended verbatim.

- penv:

  SGE parallel environment (default "mthread").

- profile_dir:

  Directory to save the profile in.

- interactive:

  If \`TRUE\` and interactive, prompt for unset values.

- overwrite:

  Overwrite an existing profile (default \`FALSE\`).

## Value

(invisibly) the path to the saved profile.
