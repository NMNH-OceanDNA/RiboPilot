# Load one sample's published assembly fastas into the database

Parses \`out/\<ID\>/assemble/\<ID\>\_assembly\_\*.fasta\` (headers
\`\>ID.path.scaffold topology\`), writes per-scaffold rows into
\`assemblies\`, and updates the per-sample \`assemble\` summary.

## Usage

``` r
load_assemblies(con, out_dir, ID)
```

## Arguments

- con:

  Open DBI connection.

- out_dir:

  The pipeline \`out/\` directory.

- ID:

  Sample ID.

## Value

(invisibly) the number of scaffolds loaded.
