# Resolve all sample Taxon strings and write a review table for confirmation

Downloads (once) a pinned NCBI taxdump into \`\<path\>/.taxonomy\`,
resolves every unique \`Taxon\` in the project, and writes
\`taxonomy_review.csv\` listing EVERY sample's proposal (not only
failures) for the user to confirm or correct. Edit the
\`confirmed_taxid\` column for any row, then run
\[apply_project_taxonomy()\].

## Usage

``` r
resolve_project_taxonomy(
  path = ".",
  taxdump_dir = file.path(path, ".taxonomy"),
  force_taxdump = FALSE
)
```

## Arguments

- path:

  Project directory.

- taxdump_dir:

  Taxdump location (default \`\<path\>/.taxonomy\`).

- force_taxdump:

  Re-download the taxdump even if present.

## Value

(invisibly) the review data.frame.
