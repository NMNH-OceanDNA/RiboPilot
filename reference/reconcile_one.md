# Reconcile barrnap + ITSx outputs for one sample; write annotations

Reconcile barrnap + ITSx outputs for one sample; write annotations

## Usage

``` r
reconcile_one(con, id, path = 1L, scaffold = 1L, gff = NULL, positions = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- gff:

  barrnap GFF file (or NULL).

- positions:

  ITSx positions file (or NULL).

## Value

(invisibly) the feature data.frame.
