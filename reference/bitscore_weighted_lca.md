# Bitscore-weighted lowest common ancestor over the top BLAST hits

Defence against a single mislabelled GenBank accession, which a one-hit
screen cannot survive. Considers only hits within `frac` of the top
bitscore (always at least `min_hits` when that many are available) and
returns the deepest rank at which they all agree.

## Usage

``` r
bitscore_weighted_lca(
  taxid_paths,
  bitscores,
  names_paths = NULL,
  frac = 0.1,
  min_hits = 3L
)
```

## Arguments

- taxid_paths:

  Character vector of ";"-joined rank-aligned taxid strings, one per
  hit, or a list of integer vectors.

- bitscores:

  Numeric vector of the same length.

- names_paths:

  Optional matching lineage-name strings, for the returned label.

- frac:

  Bitscore window below the top hit (default 0.10).

- min_hits:

  Always include at least this many top hits when available.

## Value

list(taxid, taxon, rank, n_hits).
