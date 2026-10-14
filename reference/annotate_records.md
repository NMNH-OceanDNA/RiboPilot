# Annotate every record in a multi-record fasta

Runs the rRNA anchor caller once over the whole file and ITSx once per
strand group, then reconciles per record. Records longer than
`ITSX_MAX_BP` are kept but get anchors only, because ITSx aborts the
entire run on an oversized target rather than skipping it.

## Usage

``` r
annotate_records(fasta, kingdom = "euk", cpu = 2, taxon = "all")
```

## Arguments

- fasta:

  Path to a multi-record fasta.

- kingdom:

  Anchor-caller kingdom.

- cpu:

  Threads for ITSx.

- taxon:

  ITSx `-t` profile set.

## Value

data.frame(seqid, element, type, start, stop, length, strand, model_cov,
tool).

## Details

barrnap reports coordinates in the input frame; on a genuinely
minus-strand record those anchors run 28S -\> 5.8S -\> 18S in ascending
position, the opposite of the ascending order `reconcile_features()`'s
spacer-clipping assumes. ITSx, run on a reverse-complemented copy of
those records, reports its regions in the frame where the anchors DO
ascend. So minus-strand records are reconciled in that revcomp frame –
barrnap's anchors are flipped into it, ITSx's output is used as reported
– and the finished model is flipped back to the input frame afterward.
