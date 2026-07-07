# Score one sample's candidate paths and write path_scores (+ auto_select)

Score one sample's candidate paths and write path_scores (+ auto_select)

## Usage

``` r
score_one(con, id, out_dir, taxdump_dir = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- out_dir:

  Pipeline output dir (assembly fastas + gfa live under
  \`\<out_dir\>/\<id\>/assemble\`).

- taxdump_dir:

  Optional pinned taxdump for hit-lineage resolution.

## Value

(invisibly) the scored data.frame.
