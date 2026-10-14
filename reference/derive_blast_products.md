# Derive a unit's reference panel and contamination screen from stored hits

Reads `blast_hits` rather than a hits file, so thresholds can change
without re-running a remote search. Replaces the unit's rows in
`blast_ref_candidates` and `blast_screen` and updates
`assemblies.blast_*`.

## Usage

``` r
derive_blast_products(
  con,
  id,
  path = 1L,
  scaffold = 1L,
  ts = as.integer(Sys.time()),
  review = FALSE
)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- ts:

  Timestamp to stamp the written rows with.

- review:

  TRUE derives the products of the contig's review copy: it reads
  `review_blast_hits`, replaces the row in `review_blast_screen`, and
  sets or clears the `blast_*` summary in `review_units`. The candidates
  only pick the top hit; `blast_ref_candidates` is never written, and
  the source label still comes from `blast_searched`.

## Value

(invisibly) the ranked candidate data.frame.
