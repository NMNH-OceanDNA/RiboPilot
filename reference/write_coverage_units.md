# Split a combined `samtools depth -a` (one path mapped at a time) into per-unit coverage rows

The reference names are the contig fasta headers `path.scaffold`, so
each distinct ref is written as one coverage unit.

## Usage

``` r
write_coverage_units(con, id, depth_file, clips_file = NULL)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- depth_file:

  `samtools depth -a` output (ref, pos, depth), ref = path.scaffold.

- clips_file:

  Optional per-position soft-clip counts (ref, pos, n), tab-separated.

## Value

(invisibly) the number of units written.
