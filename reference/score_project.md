# Score assembled candidate paths for a project and write \`path_scores\`

Score assembled candidate paths for a project and write \`path_scores\`

## Usage

``` r
score_project(
  path = ".",
  out_dir = file.path(path, "out"),
  taxdump_dir = file.path(path, ".taxonomy"),
  blast_db = NULL,
  con = NULL
)
```

## Arguments

- path:

  Project directory.

- out_dir:

  Pipeline output dir (default \`\<path\>/out\`).

- taxdump_dir:

  Pinned taxdump (default \`\<path\>/.taxonomy\`); optional.

- blast_db:

  Optional prebuilt BLAST db prefix for taxonomy of hits. If NULL, one
  is built from each sample's \`seeds_db\` (with species-labeled
  headers).

- con:

  Optional open DBI connection (else opened from \`\<path\>/.sqlite\`).

## Value

(invisibly) NULL.
