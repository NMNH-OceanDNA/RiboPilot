# Reference group assignment: design

Date: 2026-09-01

## Problem

The References panel shows one row per clade group with a `Samples` count and
nothing else. Which samples belong to which reference group is invisible in the
app, and the user has no way to influence grouping. The only lever is changing a
group's clade, which splits groups as a side effect.

Users need to steer assignment both before and after the automatic proposer
runs, with a way to say "these samples must stay together" that the proposer
cannot undo.

## Current state

- `samples.ref_id` -> `reference.ref_id` is the only link between samples and
  reference groups. Flat, one column.
- `suggest_reference_groups()` (`R/reference_groups.R:175`) walks each sample's
  lineage upward until a clade has at least `min_records` GenBank records
  (default 1000), capped at `max_records` (default 25000), then folds similar
  groups into shared ancestor clades.
- The proposer runs on demand from the `Propose groups` button
  (`R/app_reference.R:283`) and writes results straight to the database via
  `write_reference_groups()` and `bind_sample_ref()` (`R/reference_build.R:34`,
  `:77`).
- The References table (`R/app_reference.R:16`) shows Status, Clade, Rank,
  Samples, GenBank records, Database, Action.
- `samples.taxonomy_lock` already establishes a lock idiom in this codebase.

## Design

### Data model

Three columns, no new tables.

| Column | Type | Meaning |
| --- | --- | --- |
| `samples.user_group` | TEXT NULL | Label of a user-defined set. All samples sharing a value must land in one reference group. NULL means the sample is free. |
| `samples.ref_source` | TEXT | `auto` or `manual`. Records whether the current `ref_id` came from the proposer or from a hand edit. Drives the flag in the diff modal. |
| `reference.origin` | TEXT | `auto` or `manual`. A group the user created by hand survives a proposer run even when temporarily empty; an `auto` group with no members is dropped. |

Membership itself stays in `samples.ref_id`. Nothing else moves.

User group labels are auto-assigned `A`, `B`, `C`, ... on creation and are
renameable inline. They are labels only; no behavior keys off the text.

A user group of one sample is a no-op under keep-together semantics, so the
`Make user group` control is only enabled for a selection of two or more
samples. A single sample is moved by hand instead; the move is recorded as
`manual` and surfaces in the diff.

### Proposer

`suggest_reference_groups()` gains a pre-step: for each distinct `user_group`,
compute the MRCA of its members' lineages and represent the whole set as a
single pseudo-sample anchored at that MRCA. The existing walk and fold logic
then runs unchanged over a mix of loose samples and user-group pseudo-samples.

Consequences:

- A user group's clade floats. The proposer may walk it up or down to satisfy
  the record window, and may sweep unrelated samples into the same group.
- A user group is never split across two reference groups.
- If a user group's MRCA is very high (for example a sample from Asteroidea
  grouped with one from Scyphozoa), the walk lands on a huge clade. This is not
  refused; it produces a warning line in the diff modal when the resulting
  group exceeds `max_records`.

### UI

A toggle at the top of the References page, beside `Propose groups` and
`Settings`, switches between two views of the same data. Both views write the
same two fields, so either can be ignored entirely.

**Groups view** is today's table with a disclosure caret per row. Expanding a
row reveals its member samples with checkboxes, taxon, and user-group badge,
plus a per-group action bar.

```
REFERENCES        [Groups|Samples]   [Propose groups] [Settings]

 v Ophiuroidea      class    12   8,420   .../ophi.fa   [Build]
     [x] AST_001  Ophiothrix fragilis    User group A
     [x] AST_002  Amphiura filiformis    User group A
     [ ] AST_007  Ophiura albida         hand-set
     ...
     -> [Move selected to v]  [Make user group]  [Dissolve]
 > Scyphozoa        class     5   3,110   (none)        [Set up database]
 > Unassigned                    3     -       -        [Propose groups]
```

All existing group actions are preserved: build from GenBank, use existing
database, skip, change clade.

**Samples view** is a flat, sortable table, one row per sample, with the
reference group as an editable dropdown. Multi-row selection drives the same
`Move to group` / `Make user group` / `Dissolve` actions.

```
 ID       Taxon                  Ref group        User group   Set by
 AST_001  Ophiothrix fragilis    Ophiuroidea  v   A            auto
 AST_007  Ophiura albida         Ophiuroidea  v   -            manual
 OCT_005  Sarcophyton sp.        Unassigned   v   -            -
```

Both views need an `Unassigned` pseudo-group for samples with a NULL `ref_id`.

### Propose groups becomes a preview

`Propose groups` no longer writes directly. It computes the proposed assignment
and opens a diff modal:

```
PROPOSE GROUPS - preview changes

  AST_007  Ophiuroidea  -> Asteroidea   (was hand-set)  [make user group]
  SCY_003  Unassigned   -> Scyphozoa
  OCT_005  Octocorallia -> Anthozoa

  4 groups -> 3 groups
  ! User group B lands on Anthozoa (61,300 records, above max 25,000)

  [Apply]  [Cancel]
```

- One row per sample whose `ref_id` changes.
- Rows whose current assignment is `manual` are flagged, so a hand edit is
  never lost silently.
- Group count delta.
- A warning line per user group whose resulting clade exceeds `max_records`.
- Cancel writes nothing.

On Apply, samples the proposer moved are marked `ref_source = 'auto'`.

### Stale databases

If Apply changes an existing group's clade and that group already has a built
database, the database no longer matches the clade. The row moves to a
`Needs rebuild` status rather than silently keeping the old database. The
existing `Build` action clears it.

## Out of scope

- No schema migration for existing projects. The project is pre-release; new
  columns are added to the DDL only.
- No change to how databases are built, downloaded, or reused.
- No change to the Assemble panel; assignment lives entirely on the References
  page.
