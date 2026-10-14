# Rank-aligned lineages (names + taxids) for a set of taxids

One batched index query. Missing ranks are filled so the seven rank
slots always line up positionally, and the parallel taxid string marks a
filled rank with an empty taxid (see `td_rank_align()`).

## Usage

``` r
tk_ranked_lineage(taxids, dir)
```

## Arguments

- taxids:

  Integer vector of NCBI taxids.

- dir:

  Directory holding a pinned taxdump.

## Value

data.frame(taxid, lineage, lineage_taxids).
