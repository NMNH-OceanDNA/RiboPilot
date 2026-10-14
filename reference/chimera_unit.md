# Run every channel for one unit and write chimera + chimera_windows

Run every channel for one unit and write chimera + chimera_windows

## Usage

``` r
chimera_unit(
  con,
  id,
  path = 1L,
  scaffold = 1L,
  dbarg = NULL,
  taxdump_dir = NULL,
  label_db = NULL,
  threads = 1L,
  workdir = tempfile("chimera_"),
  review = FALSE
)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- dbarg:

  blastn -db argument for the windowed local BLAST channel, or NULL/""
  to fall back to `label_db` in identity-only mode.

- taxdump_dir:

  Local NCBI taxdump directory (or NULL for the bundled one).

- label_db:

  Path to a labels FASTA/BLAST db used when `dbarg` is absent.

- threads:

  blastn thread count.

- workdir:

  Scratch directory for intermediate files.

- review:

  TRUE checks the contig's review copy: its sequence, features, hits,
  and ITSx chimera flag come from the review tables, coverage is cut to
  the copy's window, and the result goes to `review_chimera` and
  `review_chimera_windows`.

## Value

(invisibly) the written `chimera` row, as a list.
