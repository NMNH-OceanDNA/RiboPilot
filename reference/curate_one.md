# Curate one sample's annotations and write results back to the DB

Reads the reconciled \`annotations\`, the sample's \`curate_opts\`
ruleset, any per-feature reference identities (\`blast_ref_alignment\`),
and the sample vs. hit lineages, runs \[curate_rrna_core()\], then
updates each feature's \`status\`/\`flags\` and the per-sample
\`annotate\` summary (\`n_flags\`, \`annotate_notes\`) in one
transaction.

## Usage

``` r
curate_one(con, id, ruleset = NULL)
```

## Arguments

- con:

  Open DBI connection.

- id:

  Sample ID.

- ruleset:

  Optional ruleset override (else read from \`curate_opts\`).

## Value

(invisibly) the curate_rrna_core() result.
