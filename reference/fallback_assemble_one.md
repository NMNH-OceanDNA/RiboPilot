# Run an assembly fallback for one sample (Option B2 baiting, or Option A)

Writes `<out_dir>/<id>_assembly_<n>.fasta` (a `_0` sentinel naming the
reason if nothing is recovered), matching the GetOrganelle output
contract so the rest of the pipeline is unchanged. Nothing is assembled
without a readable label database, or, for baiting, a readable seed
database.

## Usage

``` r
fallback_assemble_one(
  id,
  r1,
  r2,
  reference,
  labels_db,
  out_dir = ".",
  method = c("bait", "denovo"),
  threads = 4,
  bowtie2_opts = "",
  spades_opts = ""
)
```

## Arguments

- id:

  Sample ID.

- r1, r2:

  Trimmed paired reads.

- reference:

  Clade reference rRNA FASTA (seeds) for baiting.

- labels_db:

  Label DB FASTA for the rDNA contig filter.

- out_dir:

  Directory to write the assembly fastas into (the process cwd).

- method:

  "bait" (B2, default) or "denovo" (A, assemble all reads).

- threads:

  Threads for mapping/assembly.

- bowtie2_opts, spades_opts:

  Tool flags from the sample's assemble_opts preset.

## Value

(invisibly) number of contigs recovered. When SPAdes was killed for
memory, attribute `retry_status` holds its exit status.
