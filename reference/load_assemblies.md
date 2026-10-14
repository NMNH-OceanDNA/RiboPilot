# Load one sample's assembly fastas into the database

Parses `<dir>/<ID>_assembly_*.fasta` (headers
`>ID.path.scaffold topology`), orients each contig to the sense strand,
rotates circular units to start in the IGS, writes per-contig rows into
`assemblies`, and updates the per-sample `assemble` summary.

## Usage

``` r
load_assemblies(con, dir, ID)
```

## Arguments

- con:

  Open DBI connection.

- dir:

  Folder holding the sample's assembly fastas (the score task's own
  folder).

- ID:

  Sample ID.

## Value

(invisibly) the number of contigs loaded.
