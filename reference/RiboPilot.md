# Open the RiboPilot GUI

Launches the R-Shiny app for reviewing and curating a RiboPilot project.
Run from (or point \`project\` at) a project directory containing
\`.sqlite\`.

## Usage

``` r
RiboPilot(
  project = getwd(),
  host = NULL,
  port = NULL,
  launch.browser = NULL,
  ...
)
```

## Arguments

- project:

  Path to the project directory (default: current dir).

- host, port, launch.browser:

  Passed to Shiny (see \[run_app()\]).

- ...:

  Additional arguments passed to \`run_app()\`.
