# Curate one reconciled nrDNA feature model against a clade ruleset (pure)

Curate one reconciled nrDNA feature model against a clade ruleset (pure)

## Usage

``` r
curate_rrna_core(
  features,
  params = rrna_params(),
  ref_pid = NULL,
  sample_lineage = NA,
  hit_lineage = NA
)
```

## Arguments

- features:

  data.frame(element, type, start, stop, length, strand, ...) as
  produced by \[reconcile_features()\] / stored in the \`annotations\`
  table.

- params:

  A ruleset from \[rrna_params()\].

- ref_pid:

  Optional named numeric vector (element -\> percent identity to the top
  reference), from the \`blast_ref_alignment\` table. Missing elements
  are skipped (identity check is opportunistic).

- sample_lineage, hit_lineage:

  Optional lineage strings for the contaminant check (sample's assigned
  lineage vs. top BLAST-hit lineage).

## Value

list(features, flags, n_flags, notes) where \`features\` gains
\`status\` ("ok"/"flag") and \`flags\` (comma-separated per-feature
codes), and the sample-level \`flags\`/\`notes\` summarize missing
anchors, ordering, etc.
