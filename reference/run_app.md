# Run the RiboPilot Shiny application

Run the RiboPilot Shiny application

## Usage

``` r
run_app(
  project = getwd(),
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  host = NULL,
  port = NULL,
  launch.browser = NULL,
  ...
)
```

## Arguments

- project:

  Project directory containing \`.sqlite\`.

- onStart, options, enableBookmarking, uiPattern:

  Shiny app options.

- host, port, launch.browser:

  Shiny server options.

- ...:

  Passed to golem options.
