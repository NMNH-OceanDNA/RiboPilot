# Write scored candidate paths for a sample into the `path_scores` table

Write scored candidate paths for a sample into the `path_scores` table

## Usage

``` r
write_path_scores(con, ID, scored)
```

## Arguments

- con:

  Open DBI connection to the project database.

- ID:

  Sample ID.

- scored:

  A scored candidate data.frame from
  [`score_assembly()`](https://nmnh-oceandna.github.io/RiboPilot/reference/score_assembly.md)
  /
  [`score_candidates()`](https://nmnh-oceandna.github.io/RiboPilot/reference/score_candidates.md).

## Value

(invisibly) NULL. Existing rows for `ID` are replaced.
