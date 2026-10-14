# Select nuclear rDNA contigs from a barrnap scan

Given the merged barrnap GFF for one sample's supplied assembly, keep
the contigs that carry a nuclear rDNA anchor. A contig qualifies on an
18S or a 28S hit; a lone 5.8S hit does not, because 5.8S is short enough
that a single hit is as likely to be a degraded pseudogene copy as a
real repeat unit.

## Usage

``` r
select_rrna_contigs(gff_lines, lengths = NULL, max_contigs = 10)
```

## Arguments

- gff_lines:

  Character vector, the merged barrnap GFF.

- lengths:

  Optional named numeric of contig lengths, names = contig ids.

- max_contigs:

  Most contigs to keep, best first (default 10). Fed from the sample's
  `assemble_opts.max_scaffolds`, which is also the gate
  [`load_assemblies()`](https://nmnh-oceandna.github.io/RiboPilot/reference/load_assemblies.md)
  applies afterwards.

## Value

list with `candidates` (contig ids, best first) and `evidence` (one row
per contig scored: contig, length, rrna_features, rrna_bases, rank,
selected, reason)

## Details

Only contigs barrnap hit at all are scored, so the evidence returned is
the handful of rRNA-bearing contigs rather than every contig in a draft
genome.
