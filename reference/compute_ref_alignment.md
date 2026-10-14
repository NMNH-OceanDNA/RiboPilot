# Compute per-feature sample-vs-reference alignments for one sample

Aligns each rRNA/ITS element that both the sample and its top reference
annotate, writing percent identity and alignment coverage/length (not
the aligned strings themselves) to `blast_ref_alignment`.

## Usage

``` r
compute_ref_alignment(con, id, path = 1L, scaffold = 1L, review = FALSE)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- review:

  TRUE aligns the contig's review copy (its sequence and features) and
  writes `review_ref_alignment`; the reference tables are shared.

## Value

(invisibly) number of element alignments written.
