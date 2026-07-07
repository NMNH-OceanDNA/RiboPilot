# Compute per-feature sample-vs-reference alignments for one sample

Aligns each rRNA/ITS element that both the sample and its top reference
annotate, writing aligned strings + percent identity to
\`blast_ref_alignment\`.

## Usage

``` r
compute_ref_alignment(con, id)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

## Value

(invisibly) number of element alignments written.
