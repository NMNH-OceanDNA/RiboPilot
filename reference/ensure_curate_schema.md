# Ensure the curation columns/tables exist (idempotent)

Thin wrapper over \[migrate_db()\], kept for call sites that only need
the curation schema (e.g. \[curate_one()\]). New projects get these in
\[new_db()\].

## Usage

``` r
ensure_curate_schema(con)
```

## Arguments

- con:

  Open DBI connection.
