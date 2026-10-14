# What would bulk triage do, without doing it

A bulk write across a whole project is the highest-blast-radius action
in the app, so it is never applied without showing the count first –
including how many samples would end up with NOTHING kept, which is the
outcome that silently loses a sample.

## Usage

``` r
triage_preview(con, ids, opts = NULL, overwrite_user = FALSE)
```

## Arguments

- con:

  Open project DB connection.

- ids:

  Sample IDs to consider.

- opts:

  Named list of per-rule actions; NULL reads the sample's triage_opts
  preset.

- overwrite_user:

  Include units a curator decided by hand.

## Value

list with, for the END STATE after Apply: n_units, n_drop, n_undecided,
n_keep, n_samples, emptied, n_user_skipped, n_locked_skipped, by_rule;
and, for what Apply would actually CHANGE: n_changed (units that are not
skipped and whose stored triage/reason/rule differ from the proposed
ones, by the same comparison auto_triage_sample() uses to decide what to
write), broken down by proposed verdict into n_new_drop, n_new_keep and
n_new_flag, which therefore sum to n_changed.

## Details

If every rule is "off" for a sample, that sample is treated as fully
skipped: the preview reports its CURRENT stored state unchanged, with
zero counted in every n_new\_\* counter and in n_changed. See
dev/design/09_TRIAGE_RULES_AND_DEPTH.md section 3.3.
