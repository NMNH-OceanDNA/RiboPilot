# Reconcile barrnap + ITSx outputs for one sample; write annotations

Reconcile barrnap + ITSx outputs for one sample; write annotations

## Usage

``` r
reconcile_one(con, id, gff = NULL, positions = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- gff:

  barrnap GFF file (or NULL).

- positions:

  ITSx positions file (or NULL).

## Value

(invisibly) the feature data.frame.
