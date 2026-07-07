# Stitch one sample's selected path; write the \`stitched\` table + return fasta

Stitch one sample's selected path; write the \`stitched\` table + return
fasta

## Usage

``` r
stitch_one(con, id, fasta_out = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- fasta_out:

  Optional path to also write the stitched sequence as fasta.

## Value

(invisibly) the stitched sequence.
