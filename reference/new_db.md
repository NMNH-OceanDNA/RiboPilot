# Initialize a new RiboPilot project database

Creates the SQLite tracking database. Mirrors MitoPilot's \`new_db()\`
spine (samples / preprocess / pre_opts / assemble / assemble_opts /
blast_opts / assemblies) but with nuclear ribosomal (nrDNA) assembly
options. Annotation, curation, and export tables are added in later
build increments.

## Usage

``` r
new_db(
  db_path = "./.sqlite",
  mapping_fn = NULL,
  mapping_id = "ID",
  mapping_taxon = "Taxon",
  assemble_cpus = 6,
  assemble_memory = 24,
  seeds_db =
    "https://raw.githubusercontent.com/Smithsonian/RiboPilot/main/ref_dbs/getOrganelle/seeds/actinopterygii_rrna_seeds.fasta",
  labels_db =
    "https://raw.githubusercontent.com/Smithsonian/RiboPilot/main/ref_dbs/getOrganelle/labels/actinopterygii_rrna_labels.fasta",
  getOrganelle = paste("-F anonym", "-R 10 -k 21,45,65,85,105,115", "--larger-auto-ws"),
  target = "anonym",
  max_paths = 10,
  max_scaffolds = 10,
  min_assembly_length = 500,
  ...
)
```

## Arguments

- db_path:

  Path to the new database file.

- mapping_fn:

  Path to the mapping CSV (columns ID, Taxon, R1, R2).

- mapping_id:

  Column name to use as the primary key (default "ID").

- mapping_taxon:

  Column name with a taxonomic identifier (default "Taxon").

- assemble_cpus:

  Default number of CPUs for assembly.

- assemble_memory:

  Default memory (GB) for assembly.

- seeds_db:

  Path/URL to the GetOrganelle nr seeds database (built per clade, see
  the reference-DB builder). May be NULL until built.

- labels_db:

  Path/URL to the GetOrganelle nr labels database.

- getOrganelle:

  Default GetOrganelle command-line options.

- target:

  Read-recruitment target passed to GetOrganelle (default "anonym",
  using the custom seed/label DBs).

- max_paths:

  Max assembly paths for a sample to continue (default 10).

- max_scaffolds:

  Max scaffolds within a path to continue (default 10).

- min_assembly_length:

  Minimum contig length (bp) to retain (default 500).

- ...:

  Reserved for forward compatibility.
