# Score and store one unit's screen verdict

Skips a row already stored at the current `SCREEN_VERSION` unless
`force`, so re-running the pipeline is cheap and a scoring change forces
a recompute without a schema migration.

## Usage

``` r
screen_one(con, id, path = 1L, scaffold = 1L, force = FALSE, review = FALSE)
```

## Arguments

- con:

  Open project DB connection.

- id, path, scaffold:

  Assembly unit.

- force:

  Recompute even when the stored row is current.

- review:

  TRUE scores the contig's review copy: hits from `review_blast_hits`,
  verdict to `review_screen`. The variant columns are carried over from
  the original `screen` row, because variant grouping compares every
  contig of a sample and is not rechecked.

## Value

(invisibly) the verdict written, or the stored one when skipped.

## Details

A stored `no_reference` is always recomputed: it records missing
evidence, and the evidence routinely arrives later. A unit screened
before BLAST ran, or before the sample's lineage was resolved, would
otherwise keep that verdict for the life of the project.
