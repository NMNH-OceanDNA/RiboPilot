# Fail loudly when a fetched batch stored nothing

Mirrors the old per-accession guard (`WANTED > 0 && STORED == 0`): a
fasta that held at least one input record but stored none means NCBI
returned bytes that then failed to store, which is worth stopping on
even though `refs.fa` itself was non-empty.

## Usage

``` r
stop_if_nothing_stored(fasta_file, n_stored)
```

## Arguments

- fasta_file:

  Path to the fasta passed to
  [`write_blast_refs()`](https://nmnh-oceandna.github.io/RiboPilot/reference/write_blast_refs.md);
  named in the error and re-read here to count input records.

- n_stored:

  Number of records
  [`write_blast_refs()`](https://nmnh-oceandna.github.io/RiboPilot/reference/write_blast_refs.md)
  actually stored.

## Value

(invisibly) `n_stored`, when the check passes.
