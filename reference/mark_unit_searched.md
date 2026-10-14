# Record that a unit was searched, and by which search

Record that a unit was searched, and by which search

## Usage

``` r
mark_unit_searched(
  con,
  id,
  path = 1L,
  scaffold = 1L,
  query_md5,
  n_hits = 0L,
  source = "remote",
  db = ""
)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- path, scaffold:

  Assembly unit within the sample.

- query_md5:

  md5 of the query sequence and the search settings, as blast.nf
  computes it.

- n_hits:

  Number of hits the search returned.

- source:

  Which search produced the hits, "local" or "remote".

- db:

  The database that answered: "core_nt", or the local panel's name and
  release.

## Value

(invisibly) TRUE.
