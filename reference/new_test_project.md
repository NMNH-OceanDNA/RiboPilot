# Initialize a RiboPilot test project

Sets up a ready-to-run test project from the bundled test set: 10 fish
genome skims (USNM, study SRP313821) plus two synthetic samples, CHIM00
(a clean control) and CHIM01 (fish rDNA joined to a Penicillium 28S).
Downloads each fish sample's paired reads, writes the mapping file, and
initializes the project database and config.

## Usage

``` r
new_test_project(
  path = ".",
  n = Inf,
  full_size = FALSE,
  executor = "docker",
  container = paste0("macguigand/ribopilot:", utils::packageVersion("RiboPilot")),
  min_depth = 500,
  sra_bin = NULL,
  Rproj = TRUE,
  force = FALSE,
  fetch_metadata = TRUE,
  ...
)
```

## Arguments

- path:

  Project directory (created if missing; default current dir).

- n:

  Number of downloaded samples to include (default all 10). The bundled
  synthetic samples CHIM00 and CHIM01 are always included.

- full_size:

  Fetch the complete read sets from ENA/SRA (TRUE) instead of the small
  rDNA-filtered subsets (default FALSE).

- executor:

  Executor / config profile (default "docker").

- container:

  Container image for the pipeline.

- min_depth:

  Minimum number of reads left after trimming, counting both reads of a
  pair (2000000 is 1000000 pairs). A sample below it is not assembled
  and is marked Failed with a note.

- sra_bin:

  Directory holding `fasterq-dump` (for the SRA fallback); default
  searches `PATH`.

- Rproj:

  Initialize an RStudio project (default TRUE; no effect outside
  RStudio).

- force:

  Recreate an existing project database / config (default FALSE).

- fetch_metadata:

  Fetch each sample's NCBI metadata (synthetic samples use their source
  run's) and follow links to GBIF and GEOME (default TRUE). Linking is
  saved as the project setting. Set FALSE when offline.

- ...:

  Passed to
  [`new_project()`](https://nmnh-oceandna.github.io/RiboPilot/reference/new_project.md)
  (e.g. `custom_seeds_db`).

## Value

(invisibly) the project path.

## Details

By default the reads are the rDNA-filtered subsets published with the
RiboPilot GitHub releases (a few MB per sample, seconds to fetch). With
`full_size = TRUE` the complete read sets are fetched instead, from ENA
when reachable and falling back to NCBI SRA (`fasterq-dump`, which must
be on `PATH` or given via `sra_bin`) on networks where the ENA fastq
mirror is blocked; that is several GB and takes a while. Use `n` to
fetch fewer samples. Already-downloaded, valid read files are skipped,
so an interrupted run resumes cleanly.
