# Re-derive stored BLAST products for every unit that has hits

Replays
[`derive_blast_products()`](https://nmnh-oceandna.github.io/RiboPilot/reference/derive_blast_products.md)
over `blast_hits`, so a changed threshold or a fixed filter reaches
existing projects without re-running a remote search.

## Usage

``` r
rederive_blast_products(con, ids = NULL)
```

## Arguments

- con:

  Open project DB connection.

- ids:

  Optional character vector of sample IDs; default is every sample.

## Value

(invisibly) the number of units re-derived.

## Details

Changing which candidate ranks 1 invalidates `blast_ref_alignment`
([`compute_ref_alignment()`](https://nmnh-oceandna.github.io/RiboPilot/reference/compute_ref_alignment.md)
always aligns against the rank-1 candidate) and any curation derived
from it; this function does not mark those stale.
`inst/nextflow/bin/rp_rederive.R` is a manual tool wired into no `.nf`
process, so after running it the operator must re-run the
reference-fetch and alignment steps by hand.
