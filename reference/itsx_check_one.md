# Pipeline-side check for one unit from ITSx output already on disk

The annotation step runs stock ITSx next to its own tuned run and hands
the `<out_prefix>.*` files to reconcile, which calls this after the
unit's annotations are written. ITSx exits 0 on fatal errors, so a run
counts only when `<out_prefix>.summary.txt` says it finished and
`<out_prefix>.err` has no error line; otherwise the contig stores
`not_run(ITSx did not finish)`.

## Usage

``` r
itsx_check_one(con, id, path = 1L, scaffold = 1L, out_prefix, review = FALSE)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- out_prefix:

  ITSx output prefix (files `<out_prefix>.positions.txt` etc.).

- review:

  TRUE checks the contig's review copy (its sequence and features) and
  stores the result on the copy instead of the annotate row.

## Value

one-row data.frame as from `itsx_check_units()`, invisibly.
