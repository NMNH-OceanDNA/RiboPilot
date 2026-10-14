# Is a hit at least as close as `min_rank`?

Is a hit at least as close as `min_rank`?

## Usage

``` r
rank_at_least(d, min_rank = "phylum")
```

## Arguments

- d:

  Distance from
  [`rank_distance()`](https://nmnh-oceandna.github.io/RiboPilot/reference/rank_distance.md).

- min_rank:

  One of `RRNA_LINEAGE_RANKS`; "kingdom" is the most permissive.

## Value

logical.
