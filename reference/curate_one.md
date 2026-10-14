# Curate one sample's annotations and write results back to the DB

Reads the reconciled `annotations`, the sample's `curate_opts` ruleset,
any per-feature reference identities (`blast_ref_alignment`), and the
sample vs. hit lineages, runs
[`curate_rrna_core()`](https://nmnh-oceandna.github.io/RiboPilot/reference/curate_rrna_core.md),
then updates each feature's `status`/`flags` and the per-sample
`annotate` summary (`n_flags`, `curate_flags`) in one transaction.
`partial_start`/`partial_stop` are only set (to 1) where still NULL, so
a curator's own choice is kept.

## Usage

``` r
curate_one(con, id, path = 1L, scaffold = 1L, ruleset = NULL, review = FALSE)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- ruleset:

  Optional ruleset override (else read from `curate_opts`).

- review:

  TRUE curates the contig's review copy: it reads and writes the review
  tables and replaces the copy's partial ends on every run (the user's
  choice is kept in `partial_start_user` and `partial_stop_user`).

## Value

(invisibly) the curate_rrna_core() result.
