# Screen every unit of one sample

Screen every unit of one sample

## Usage

``` r
screen_sample(con, id, force = FALSE, variants = TRUE)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- force:

  Recompute rows already at the current version.

- variants:

  Also classify on-target units as isomers or paralogs. The pipeline's
  screen step passes FALSE; its triage step runs that check once
  annotation is done.

## Value

(invisibly) the number of units screened, with attribute "written":
data.frame(path, scaffold) of the units whose screen row was rewritten
rather than skipped as current.
