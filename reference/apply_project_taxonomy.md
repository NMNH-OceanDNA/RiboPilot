# Freeze confirmed taxids + lineages from the review CSV into \`samples\`

Reads \`taxonomy_review.csv\`, (re)fetches lineages for any
\`confirmed_taxid\` the user changed, and writes \`taxid\` + \`lineage\`
into the \`samples\` table joined on \`Taxon\`. Rows with no confirmed
taxid are flagged (left NULL).

## Usage

``` r
apply_project_taxonomy(path = ".", taxdump_dir = file.path(path, ".taxonomy"))
```

## Arguments

- path:

  Project directory.

- taxdump_dir:

  Taxdump location (default \`\<path\>/.taxonomy\`).

## Value

(invisibly) the frozen per-Taxon data.frame.
