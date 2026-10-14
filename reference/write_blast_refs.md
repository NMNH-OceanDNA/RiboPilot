# Store every record of a multi-record fasta as a reference

Sequence comes from the fasta; element coordinates come from `ann`
rather than from an NCBI feature table. An accession with no rows in
`ann` is still stored, with no annotations.

## Usage

``` r
write_blast_refs(con, fasta_file, ann = NULL)
```

## Arguments

- con:

  Open project DB connection.

- fasta_file:

  Multi-record fasta, deflines `>ACC.V description`.

- ann:

  Annotation frame carrying seqid, element, type, start, stop, strand,
  or NULL.

## Value

(invisibly) the number of accessions stored.
