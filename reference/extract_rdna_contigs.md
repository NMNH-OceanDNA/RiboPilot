# Keep the rDNA contigs (BLAST hit vs the label DB) and write assembly fastas

Filters \`contigs_fa\` to sequences long enough and with a BLAST hit
against \`labels_db\`, then writes them as one candidate path in the
standard \`\<id\>\_assembly_1.fasta\` header format (scaffolds ordered
by length).

## Usage

``` r
extract_rdna_contigs(
  contigs_fa,
  labels_db,
  id,
  out_fa,
  min_contig = 300,
  min_pident = 70,
  max_evalue = 1e-10
)
```

## Arguments

- contigs_fa:

  Assembler contigs FASTA.

- labels_db:

  Label DB FASTA (single-locus rRNA records).

- id:

  Sample ID.

- out_fa:

  Output \`\<id\>\_assembly_1.fasta\` path.

- min_contig:

  Minimum contig length to consider.

- min_pident, max_evalue:

  BLAST thresholds to call a contig rDNA.

## Value

Number of contigs kept (0 if none).
