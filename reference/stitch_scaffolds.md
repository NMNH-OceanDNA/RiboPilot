# Order + orient a path's scaffolds along a reference and stitch them

Order + orient a path's scaffolds along a reference and stitch them

## Usage

``` r
stitch_scaffolds(
  path_fasta,
  reference_fasta,
  max_gap = 2000,
  min_hit = 100,
  min_new = 200
)
```

## Arguments

- path_fasta:

  Fasta of one path's scaffolds.

- reference_fasta:

  Fasta with a single reference rDNA sequence (e.g. the top BLAST hit's
  full 18S-ITS-28S sequence).

- max_gap:

  Cap (bp) on the N-gap inserted between non-abutting scaffolds (default
  2000); larger reference gaps are clamped to this.

- min_hit:

  Minimum BLAST alignment length (bp) to place a scaffold (default 100);
  unplaced scaffolds are appended at the end and flagged.

- min_new:

  Minimum new reference coverage (bp) a scaffold must add to be kept in
  the tiling; scaffolds contained within a larger one are dropped as
  redundant (default 200).

## Value

list(sequence = character stitched sequence, layout =
data.frame(scaffold, ref_start, ref_end, strand, placed, length,
redundant, order), n_gaps, n_unplaced, n_redundant)
