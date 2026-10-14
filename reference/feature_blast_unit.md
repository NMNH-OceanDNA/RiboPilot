# BLAST each core feature of a unit and store the calls

Features under FEATURE_BLAST_MIN_LEN are not searched. With no `dbarg`,
or when blastn fails, searchable features are stored as not_run with a
note.

## Usage

``` r
feature_blast_unit(
  con,
  id,
  path = 1L,
  scaffold = 1L,
  dbarg = NULL,
  taxdump_dir = NULL,
  threads = 1L,
  workdir = tempfile("feature_blast_"),
  review = FALSE,
  blastn = "blastn"
)
```

## Arguments

- review:

  TRUE reads and writes the contig's review copy.

## Value

(invisibly) the stored calls.
