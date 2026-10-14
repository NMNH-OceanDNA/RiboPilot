# Open the RiboPilot app

Launches the R-Shiny app for reviewing and curating a RiboPilot project.
Run from (or point `project` at) a project directory containing
`.sqlite`. `RiboPilot()` is the same function.

## Usage

``` r
run_app(
  project = getwd(),
  host = NULL,
  port = NULL,
  launch.browser = NULL,
  headless = NULL
)

RiboPilot(
  project = getwd(),
  host = NULL,
  port = NULL,
  launch.browser = NULL,
  headless = NULL
)
```

## Arguments

- project:

  Project directory containing `.sqlite` (default: current dir).

- host, port, launch.browser:

  Passed to Shiny. To reach the app from another computer, pass a `port`
  and `launch.browser = FALSE`, then open an SSH tunnel to that port;
  the console prints the command.

- headless:

  Offer a cluster job script instead of Run from App. `NULL` (default)
  turns it on when a batch scheduler (sbatch, qsub, or bsub) is on PATH
  or the session is RStudio Server; pass `TRUE` or `FALSE` to override.
