# Run the assembly workflow (WF1) for a project

Launches the self-contained Nextflow pipeline: it reads samples from the
project database, runs fastp + GetOrganelle, then scores and stitches
each assembly and computes the assembly-review aids (reference BLAST +
coverage), writing everything back to the database. Annotation and
curation are handled by the annotation workflow (\[run_annotate()\],
WF2) after you review and lock the assemblies in the app.

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

  Project directory (contains \`.sqlite\` and \`.config\`).

- nextflow:

  Path to the \`nextflow\` binary.

- work:

  Nextflow work directory (default \`\<path\>/work\`).

- conda_env:

  Optional conda env providing the pipeline tools; its \`bin/\` is
  prepended to PATH for the run.

- ribopilot_src:

  Path to the RiboPilot package source for the R-script processes (only
  needed when the package is not installed; defaults to the loaded
  package location).

## Value

(invisibly) the Nextflow exit status.
