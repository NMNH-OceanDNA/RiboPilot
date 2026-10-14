# Run the assembly workflow for a project

Launches the self-contained Nextflow pipeline, the same run as the app's
Update button. For every taxonomy-locked sample that is Ready to run,
one run cleans the reads, assembles, scores, BLASTs, maps coverage,
annotates, curates, screens, checks for chimeras, and triages every
contig, writing everything back to the database. Run it again (or press
Update) only after changing a selection, an option set, or a State.

## Usage

``` r
run_assembly(
  path = ".",
  nextflow = "nextflow",
  work = file.path(path, "work"),
  conda_env = NULL,
  ribopilot_src = NULL
)
```

## Arguments

- path:

  Project directory (contains `.sqlite` and `.config`).

- nextflow:

  Path to the `nextflow` binary.

- work:

  Nextflow work directory (default `<path>/work`).

- conda_env:

  Optional conda env providing the pipeline tools; its `bin/` is
  prepended to PATH for the run.

- ribopilot_src:

  Path to the RiboPilot package source for the R-script processes (only
  needed when the package is not installed; defaults to the loaded
  package location).

## Value

(invisibly) the Nextflow exit status.
