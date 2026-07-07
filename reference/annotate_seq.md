# Annotate the nuclear rDNA feature model on a single sequence (all in R)

Annotate the nuclear rDNA feature model on a single sequence (all in R)

## Usage

``` r
annotate_seq(fasta, cpu = 2)
```

## Arguments

- fasta:

  Path to a single-sequence fasta (e.g. a stitched path).

- cpu:

  CPUs for ITSx.

## Value

data.frame(element, type, start, stop, length, strand, tool).
