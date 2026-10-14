# Curate one reconciled nrDNA feature model against a clade ruleset (pure)

Curate one reconciled nrDNA feature model against a clade ruleset (pure)

## Usage

``` r
curate_rrna_core(
  features,
  params = rrna_params(),
  ref_pid = NULL,
  ref_cov = NULL,
  sample_lineage = NA,
  hit_lineage = NA,
  hit_is_nuclear_rrna = NA,
  sample_taxids = NA_character_,
  hit_taxids = NA_character_
)
```

## Arguments

- features:

  data.frame(element, type, start, stop, length, strand, ...) as
  produced by `reconcile_features()` / stored in the `annotations`
  table.

- params:

  A ruleset from `rrna_params()`.

- ref_pid:

  Optional named numeric vector (element -\> percent identity to the top
  reference), from the `blast_ref_alignment` table. Missing elements are
  skipped (identity check is opportunistic).

- ref_cov:

  Optional named numeric vector (element -\> fraction of the element
  covered by the reference alignment, `blast_ref_alignment.aln_cov`).
  The identity check is skipped when coverage is below
  `min_ref_coverage`, so a partial reference cannot raise a false
  low_identity flag.

- sample_lineage, hit_lineage:

  Optional lineage strings for the contaminant check: the sample's
  assigned lineage vs. the lineage of the BLAST screen's unfiltered top
  hit (`blast_screen.lineage`), not the nuclear-rDNA reference
  candidates. An unknown lineage on either side raises
  `contaminant_unchecked` rather than passing silently.

- hit_is_nuclear_rrna:

  Optional 0/1 (`blast_screen.is_nuclear_rrna`): whether that top hit is
  a nuclear rDNA record at all. 0 raises `not_nuclear_rdna`.

- sample_taxids, hit_taxids:

  Optional rank-aligned taxid strings matching those lineages
  (`samples.lineage_taxids`, `blast_screen.lineage_taxids`). The
  comparison prefers them and falls back to the names when absent.

## Value

list(features, flags, n_flags, notes) where `features` gains `status`
("ok"/"flag"), `flags` (comma-separated per-feature codes), and
`partial_start`/`partial_stop` (1 on the free 5'/3' ends of a truncated
feature, in its own orientation; NA elsewhere), and the sample-level
`flags`/`notes` summarize missing anchors, ordering, etc.
