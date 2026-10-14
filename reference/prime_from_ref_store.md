# Copy stored references for a project's candidates into the project database

The fetch step selects only candidate accessions absent from
`blast_ref_sequences`, so anything copied here is not fetched again. An
accession the project already holds is left alone.

## Usage

``` r
prime_from_ref_store(con, ref_dir, ids = NULL)
```

## Arguments

- con:

  Open project DB connection.

- ref_dir:

  Reference folder, as resolved by `ref_config()`.

- ids:

  Optional sample IDs to limit to; default is every candidate.

## Value

(invisibly) the number of accessions copied.
