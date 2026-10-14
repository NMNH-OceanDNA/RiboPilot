# Record a triage decision for one assembly unit

Writes the verdict, its justification and its owner, derives
`assemblies.ignore` from the verdict, and reconciles the unit's annotate
row (held on drop, restored on keep) through `reconcile_ignore_unit()`.

## Usage

``` r
set_unit_triage(
  con,
  id,
  path,
  scaffold,
  verdict,
  reason = NA_character_,
  source = c("user", "auto"),
  rule = NA_character_
)
```

## Arguments

- con:

  Open project DB connection.

- id, path, scaffold:

  Assembly unit.

- verdict:

  "keep", "drop", or NA/NULL for undecided.

- reason:

  Short plain-English justification, shown to the curator.

- source:

  "user" or "auto".

- rule:

  Key of the rule that fired, or NA for a user decision.

## Value

(invisibly) the derived ignore flag.
