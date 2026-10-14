# Apply automatic triage to one sample

Skips any unit a human has already decided (`triage_source = 'user'`)
and any unit whose annotation is locked, so re-running is safe. If every
rule is set to "off", this is a complete no-op and writes nothing: see
dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3 for why that special
case exists.

## Usage

``` r
auto_triage_sample(con, id, opts = NULL, overwrite_user = FALSE)
```

## Arguments

- con:

  Open project DB connection.

- id:

  Sample ID.

- opts:

  Named list of per-rule actions; NULL reads the sample's triage_opts
  preset.

- overwrite_user:

  Also re-decide units a curator set by hand.

## Value

(invisibly) the number of units changed.
