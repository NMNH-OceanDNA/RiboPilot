# Write per-position read depth for one sample's stitched sequence

Write per-position read depth for one sample's stitched sequence

## Usage

``` r
write_coverage(con, id, depth_file)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- depth_file:

  \`samtools depth -a\` output (ref, pos, depth), tab-separated.

## Value

(invisibly) the mean depth.
