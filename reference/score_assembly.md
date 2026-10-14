# Score GetOrganelle candidate paths from real assembly output

Score GetOrganelle candidate paths from real assembly output

## Usage

``` r
score_assembly(
  path_fastas,
  gfa_path = NA_character_,
  sample_lineage = NA_character_,
  sample_taxids = NA_character_,
  blast_db = NULL,
  taxdump_dir = NULL,
  with_staxids = FALSE,
  expected_length = c(4000, 12000),
  weights = DEFAULT_SCORE_WEIGHTS
)
```

## Arguments

- path_fastas:

  Per-path fasta files, each named with its path number last
  (`<id>_assembly_<n>.fasta`); `n` becomes the `path` column.

- gfa_path:

  The assembly graph `.gfa` (for depth). May be NA.

- sample_lineage:

  The sample's resolved rank-aligned lineage.

- sample_taxids:

  The sample's rank-aligned lineage taxids. Preferred over the names
  when present: taxids carry no homonyms.

- blast_db:

  BLAST nucleotide db prefix (e.g. the Actinopterygii seeds db, or an
  nt-style db).

- taxdump_dir:

  Pinned taxdump dir (for resolving hit lineages).

- with_staxids:

  TRUE if `blast_db` carries taxids; else taxonomy is taken from the
  seed-DB subject id species.

- expected_length, weights:

  Passed to
  [`score_candidates()`](https://nmnh-oceandna.github.io/RiboPilot/reference/score_candidates.md).

## Value

Scored candidate data.frame (see
[`score_candidates()`](https://nmnh-oceandna.github.io/RiboPilot/reference/score_candidates.md)).
