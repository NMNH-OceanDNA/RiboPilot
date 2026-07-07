# Store one fetched NCBI reference (sequence + features) in the project DB

Store one fetched NCBI reference (sequence + features) in the project DB

## Usage

``` r
write_blast_ref(con, accession, fasta_file, ft_file, title = NA_character_)
```

## Arguments

- con:

  Open project DB connection.

- accession:

  NCBI accession.

- fasta_file:

  efetch \`-format fasta\` output.

- ft_file:

  efetch \`-format ft\` output.

- title:

  Optional defline title (for species).

## Value

(invisibly) the reference sequence length.
