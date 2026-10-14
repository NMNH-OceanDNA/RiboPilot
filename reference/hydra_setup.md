# Configure the R session environment for the NMNH Hydra cluster

RStudio Server sessions on Hydra start with a stripped `PATH` that omits
the Univa Grid Engine, Java, and the user's `~/bin` directories, so
Nextflow cannot find `qsub` and every job submission fails. Call this
once per session, before launching the app, to put those directories
back on `PATH`. Has no effect (and warns) when not running on Hydra.

## Usage

``` r
hydra_setup()
```

## Value

Invisibly, `TRUE` if the Hydra environment was applied, otherwise
`FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
RiboPilot::hydra_setup()
RiboPilot::RiboPilot(project = "/scratch/me/my_project")
} # }
```
