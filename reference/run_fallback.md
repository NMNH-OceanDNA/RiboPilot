# Run an assembly fallback (Option B2 baiting / Option A de-novo)

Queues samples for the fallback workflow and launches it. By default it
targets samples that failed the primary GetOrganelle assembly QC gate
(no assembly, or a missing 18S/28S anchor); pass \`ids\` to force a
specific set (e.g. an Option A whole-assembly cross-check of a good
sample).

## Usage

``` r
run_fallback(
  path = ".",
  ids = NULL,
  method = c("bait", "denovo"),
  nextflow = "nextflow",
  work = file.path(path, "work"),
  conda_env = NULL,
  ribopilot_src = NULL
)
```

## Arguments

- path:

  Project directory (contains \`.sqlite\` and \`.config\`).

- ids:

  Sample IDs to run (default: the QC-gate failures).

- method:

  "bait" (B2 reference-guided baiting, default) or "denovo" (Option A,
  assemble all reads).

- nextflow:

  Path to the \`nextflow\` binary.

- work:

  Nextflow work directory (default \`\<path\>/work\`).

- conda_env:

  Optional conda env providing the tools; its \`bin/\` is prepended to
  PATH for the run.

- ribopilot_src:

  Package source for the R-script processes.

## Value

(invisibly) the Nextflow exit status (or 0 if nothing queued).
