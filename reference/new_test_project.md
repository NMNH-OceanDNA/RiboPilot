# Initialize a RiboPilot test project

Sets up a ready-to-run test project from the bundled 10-sample test set
(USNM fish genome skims, study SRP313821): downloads each sample's
paired reads, writes the mapping file, and initializes the project
database + config.

## Usage

``` r
new_test_project(
  path = ".",
  n = Inf,
  executor = "conda",
  container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
  min_depth = 500,
  sra_bin = NULL,
  Rproj = TRUE,
  force = FALSE,
  ...
)
```

## Arguments

- path:

  Project directory (created if missing; default current dir).

- n:

  Number of samples to include (default all 10).

- executor:

  Executor / config profile (default "conda").

- container:

  Container image for the pipeline.

- min_depth:

  Minimum trimmed paired reads to proceed with assembly.

- sra_bin:

  Directory holding \`fasterq-dump\` (for the SRA fallback); default
  searches \`PATH\`.

- Rproj:

  Initialize an RStudio project (default TRUE; no effect outside
  RStudio).

- force:

  Recreate an existing project database / config (default FALSE).

- ...:

  Passed to \[new_project()\] (e.g. \`custom_seeds_db\`).

## Value

(invisibly) the project path.

## Details

Reads are fetched from ENA when reachable, falling back to NCBI SRA
(\`fasterq-dump\`, which must be on \`PATH\` or given via \`sra_bin\`)
on networks where the ENA fastq mirror is blocked. The full set is
several GB and takes a while; use \`n\` to fetch fewer samples for a
quicker setup. Already-downloaded, valid read files are skipped, so an
interrupted run resumes cleanly.
