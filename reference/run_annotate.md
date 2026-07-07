# Run the annotation workflow (WF2) for a project

Annotates and curates the assemblies you have reviewed and locked in the
app (\`assemble_lock = 1\`): stitch the selected path -\> barrnap + ITSx
-\> reconcile into the feature model -\> reference alignment -\> curate.
Also re-annotates samples whose selected path later changed. Samples are
queued via \`annotate_switch\` (see the Annotate tab).

## Usage

``` r
run_annotate(
  path = ".",
  nextflow = "nextflow",
  work = file.path(path, "work"),
  conda_env = NULL,
  ribopilot_src = NULL
)
```

## Arguments

- path:

  Project directory.

- nextflow:

  Path to the nextflow binary.

- work:

  Nextflow work dir (default \`\<path\>/work\`).

- conda_env:

  Optional conda env providing the tools (see \[run_assembly()\]).

- ribopilot_src:

  Package source for the R-script processes.

## Value

(invisibly) the Nextflow exit status.
