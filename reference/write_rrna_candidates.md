# Write the rDNA contig search evidence for one sample

Write the rDNA contig search evidence for one sample

## Usage

``` r
write_rrna_candidates(con, id, evidence)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- evidence:

  data.frame from
  [`select_rrna_contigs()`](https://nmnh-oceandna.github.io/RiboPilot/reference/select_rrna_contigs.md).

## Value

(invisibly) rows written.
