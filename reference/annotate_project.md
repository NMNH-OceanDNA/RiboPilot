# Annotate assembled samples for a project

For each sample with a scored/auto-selected path (falls back to path 1),
stitches that path against its top-BLAST-hit reference, runs
\[annotate_seq()\], and writes the \`stitched\`, \`annotations\`, and
\`annotate\` tables.

## Usage

``` r
annotate_project(
  path = ".",
  out_dir = file.path(path, "out"),
  con = NULL,
  cpu = 2
)
```

## Arguments

- path:

  Project directory.

- out_dir:

  Pipeline output dir (default \`\<path\>/out\`) - unused directly but
  kept for symmetry; sequences come from the \`assemblies\` table.

- con:

  Optional open DBI connection.

- cpu:

  CPUs for ITSx.

## Value

(invisibly) NULL.
