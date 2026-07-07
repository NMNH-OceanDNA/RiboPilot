# Resolve sample Taxon strings to NCBI taxids + lineages (for confirmation)

Prefers a pinned local taxdump (fast, offline, reproducible) when
\`taxdump_dir\` has one; otherwise falls back to live NCBI \`datasets\`.
EVERY row is returned for user confirmation (see \`match_type\`), not
only failures.

## Usage

``` r
resolve_taxonomy(taxa, taxdump_dir = NULL)
```

## Arguments

- taxa:

  Character vector of Taxon strings (duplicates de-duplicated).

- taxdump_dir:

  Directory of a pinned taxdump (see \[ensure_taxdump()\]). If it
  contains a taxdump, resolution is local; otherwise network
  \`datasets\` is used.

## Value

data.frame, one row per unique Taxon: \`match_type\` in
"exact"/"fuzzy_species"/"fuzzy_genus"/"unresolved", plus proposed
\`taxid\`, \`sci_name\`, \`rank\`, \`lineage\`, \`alternatives\`.
