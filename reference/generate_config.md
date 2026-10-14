# Generate and save a reusable cluster config profile

Configure your HPC cluster once and reuse it for every project. Mirrors
MitoPilot's `generate_config()`. Builds a Nextflow config from a generic
scheduler template, fills in cluster-level settings, and saves it as a
named profile in `ribopilot_config_dir()`. Afterwards
`new_project(executor = "<name>")` will find and use it automatically.

## Usage

``` r
generate_config(
  name,
  scheduler = c("slurm", "sge", "pbs", "lsf", "local"),
  container_engine = c("auto", "singularity", "apptainer", "docker"),
  container_cache = NULL,
  container_run_options = NULL,
  queue = NULL,
  walltime = NULL,
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

  Profile name. Saved as `config.<name>`; pass as `executor`.

- scheduler:

  Base template: "slurm", "sge", "pbs", "lsf", or "local" (Docker on one
  machine, the same as executor = "local").

- container_engine:

  "auto" (default) picks singularity for the schedulers and docker for
  "local". "singularity", "apptainer", and "docker" apply to slurm, sge,
  pbs, and lsf; "local" always uses Docker.

- container_cache:

  Folder where Singularity or Apptainer keeps the pulled image, shared
  by every project that uses this profile; it must be visible from the
  compute nodes. Default: `NXF_SINGULARITY_CACHEDIR` (or
  `NXF_APPTAINER_CACHEDIR`) when set, otherwise a "containers" folder in
  the reference folder (`RIBOPILOT_REF_DIR`, see ribopilot_ref_dir()).
  Pass "" to leave it to Nextflow (the project's work/singularity
  folder).

- container_run_options:

  Optional extra runOptions for the container engine (any engine),
  appended to the bind mounts.

- queue:

  Partition / queue name. If `NULL`, the directive is omitted.

- walltime:

  Per-job walltime, in a form your scheduler accepts (e.g. "24h",
  "48:00:00"). Many SLURM and PBS sites reject a job that requests no
  walltime instead of applying a default. If `NULL`, the directive is
  omitted and the site default applies.

- account:

  Optional accounting string, folded into clusterOptions.

- cluster_options:

  Optional raw clusterOptions string appended verbatim.

- penv:

  SGE parallel environment (default "mthread").

- profile_dir:

  Directory to save the profile in.

- interactive:

  If `TRUE` and interactive, prompt for unset values.

- overwrite:

  Overwrite an existing profile (default `FALSE`).

## Value

(invisibly) the path to the saved profile.
