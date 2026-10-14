# Score one sample's candidate paths and write path_scores (+ auto_select)

Score one sample's candidate paths and write path_scores (+ auto_select)

## Usage

``` r
score_one(con, id, dir, taxdump_dir = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- dir:

  Folder holding the sample's `<id>_assembly_*.fasta` files and assembly
  graph (`.gfa`): the score task's own folder.

- taxdump_dir:

  Optional pinned taxdump for hit-lineage resolution.

## Value

(invisibly) the scored data.frame.
