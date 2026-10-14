# Taxonomic distance between a hit and the sample, in rank steps

0 = same species, 1 = same genus, 2 = family, 3 = order, 4 = class, 5 =
phylum, 6 = kingdom, 7 = no shared taxon at any rank.

## Usage

``` r
rank_distance(hit_taxids, sample_taxids)
```

## Arguments

- hit_taxids, sample_taxids:

  Rank-aligned taxid vectors from `parse_rank_taxids()`, or ";"-joined
  taxid strings.

## Value

integer in 0..7.

## Details

Distance is deliberately NOT a measure of sequence divergence. A large
distance on its own is the normal state for an under-represented taxon
and must never be read as contamination; it only becomes evidence when
paired with an implausibly high percent identity (see the screen
envelope).
