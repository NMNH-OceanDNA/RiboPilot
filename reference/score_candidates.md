# Score and rank candidate paths, auto-selecting the best

Score and rank candidate paths, auto-selecting the best

## Usage

``` r
score_candidates(
  candidates,
  sample_lineage = NA_character_,
  expected_length = c(4000, 12000),
  weights = DEFAULT_SCORE_WEIGHTS
)
```

## Arguments

- candidates:

  data.frame, one row per candidate path, with columns: \`path\`,
  \`depth\` (mean graph depth), \`blast_pident\` (0-100),
  \`blast_lineage\` ("; " ranked lineage of the top hit), \`anchors\`
  (comma list of distinct anchors found among 18S/5.8S/28S), \`length\`
  (total bp).

- sample_lineage:

  The sample's resolved "; " ranked lineage.

- expected_length:

  Numeric length-2 c(min, max) plausible rDNA repeat unit length for the
  clade (default c(4000, 12000)).

- weights:

  Named weights (see \[DEFAULT_SCORE_WEIGHTS\]).

## Value

\`candidates\` with added sub-score columns (\`s_tax\`, \`s_identity\`,
\`s_anchors\`, \`s_depth\`, \`s_length\`), a \`score\`, a \`rank\` (1 =
best), and an \`auto_select\` logical flagging the top candidate.
