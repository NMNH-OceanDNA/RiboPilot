# Write scored candidate paths for a sample into the \`path_scores\` table

Write scored candidate paths for a sample into the \`path_scores\` table

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

  A scored candidate data.frame from \[score_assembly()\] /
  \[score_candidates()\].

## Value

(invisibly) NULL. Existing rows for \`ID\` are replaced.
