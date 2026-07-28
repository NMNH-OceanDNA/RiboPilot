# 09 - Explicit triage rules, per-rule control, and depth as evidence

Status: design agreed 2026-07-28. Follows Phase 2 (tri-state triage, schema v11).
Supersedes the triage rule set described in
`dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md` section 4.1 while keeping its
governing principles unchanged.

## 0. The one-paragraph summary

Auto-triage today is an anonymous if/else chain inside `triage_rules()`. The
curator sees a single prose sentence in the bulk modal, cannot tell which rule
produced which outcome, cannot turn a rule off, and cannot see what a rule would
cost. Worse, the chain reads "no features" by counting rows in `annotations`,
which returns zero both for a scaffold that was annotated and had nothing and
for a scaffold whose annotation crashed, so a tool failure can drop real data.
This design replaces the chain with an explicit, named, user-controllable rule
registry over the two independent evidence channels, gates every combination
rule on its channels actually being present, and removes depth from the decision
entirely.

## 1. Decisions

1. Rules are named, ordered, and individually controllable: **Off / Flag / Drop**.
2. Drop is offered only on rules that combine both evidence channels.
3. Combination rules are gated on both channels having actually run.
4. Depth never produces a verdict or a reason. It is displayed, not acted on.
5. Rule settings live in a `triage_opts` preset table so the app and the WF1
   `triage` process cannot disagree.

## 2. The evidence channels

Exactly two channels feed triage today:

- **Annotation.** barrnap plus ITSx, reconciled per unit. Present when
  `annotate.n_features IS NOT NULL`, which `reconcile_one()` writes on every
  run including one that found nothing (R/pipeline_steps.R:153-156). `n_features`
  is written only by reconcile, so it survives a re-queue, unlike
  `annotate_switch`.
- **GenBank hit.** The unfiltered best core_nt hit in `blast_screen`, and its
  `is_nuclear_rrna` flag.

The channels are independent: annotation is local and offline, the hit is remote
and reference-dependent. That independence is the entire safety argument for
dropping on their conjunction, and it is why a rule may only drop when both
channels are present.

### 2.1 Why the annotation gate is a correctness fix, not a nicety

`unit_triage_state()` derives the feature count by counting rows in
`annotations` (R/triage.R:72-73). `annotate_tools` carries `errorStrategy
'ignore'` (inst/nextflow/modules/pipeline_r.nf:141) so that one scaffold's
barrnap or ITSx failure cannot sink a sample. A unit whose annotation crashed
therefore has zero rows, reports "no rRNA or ITS features found", and if its
best hit also fails the title regex it is **dropped by rule**. That is the
current `feats == 0 && nuc == 0` branch firing on missing data rather than on
evidence.

Residual caveat, recorded rather than solved here: `ITSx ... || true` followed by
`[ -f unit.itsx.positions.txt ] || : > unit.itsx.positions.txt`
(pipeline_r.nf:172-173) masks an ITSx crash as an empty positions file, so
`n_features = 0` means "barrnap found no rRNA, and ITSx either found no ITS or
failed". barrnap has no `|| true`, so a barrnap failure fails the process and is
caught by the gate. The rRNA half, which is the half the rules depend on, is
sound.

## 3. The rule registry

`triage_rules()` becomes a registry rather than a chain. Each rule is a record:

```
key        stable identifier, stored in assemblies.triage_rule
label      user-facing name, shown in the modal and as the reason prefix
meaning    one-line plain-English description, shown in the modal
channels   which evidence channels it requires
actions    the actions the rule may be set to
default    the action it ships with
test       pure function (row, ctx) -> reason string or NA
```

Evaluation is first match wins over the enabled rules, in the order below.
Nothing firing means `keep`.

| # | key | label | Annotation | GenBank hit | Default | Actions |
|---|---|---|---|---|---|---|
| 1 | `not_annotated` | Annotation has not run | not run | any | Flag | Off / Flag |
| 2 | `no_screen` | GenBank screen has not run | any | absent | Flag | Off / Flag |
| 3 | `no_feat_non_rdna` | No annotations, non-rDNA hit | ran, 0 features | not nuclear rDNA | **Drop** | Off / Flag / Drop |
| 4 | `no_feat_rdna_hit` | No annotations but rDNA hit | ran, 0 features | nuclear rDNA | Flag | Off / Flag / Drop |
| 5 | `feat_non_rdna_hit` | Has annotations but non-rDNA hit | ran, >= 1 feature | not nuclear rDNA | Flag | Off / Flag / Drop |
| 6 | `odd_one_out` | Disagrees with sibling scaffolds | any | present | Flag | Off / Flag |
| - | - | *nothing fires* | ran, >= 1 feature | nuclear rDNA | keep | - |

Rules 3, 4 and 5 are the 2x2 over the two channels; the fourth cell (features
present and an rDNA hit) is the clean state and needs no rule.

Rules 1 and 2 are prerequisites and must precede 3-5, because 3-5 cannot be
evaluated without both channels. A scaffold that was annotated but never
screened is caught by rule 2 rather than being judged against a missing hit.
Rule 6 also needs the screen and so necessarily follows rule 2.

### 3.1 Why only 3-5 may drop

A drop requires that both independent channels have actually reported on the
scaffold. That is the invariant, and it is weaker than concurrence: only rule 3
is a case where both channels independently say the scaffold is not the target
locus, and only rule 3 drops by default. Rules 4 and 5 fire on channel
DISAGREEMENT and are droppable only if a curator opts in, having decided which
channel they trust for their material. Do not describe the scheme as "drops only
when both channels agree"; that is true of the default configuration, not of the
rule set.

Rules 1, 2 and 6 rest on one channel or on its absence, and may never drop:

- `not_annotated` and `no_screen` assert only that evidence is missing.
- `odd_one_out` is pure taxonomic distance, which
  `08_MULTIUNIT_UX_AND_CONTAMINATION.md` section 3.2 forbids as a drop basis by
  name: low identity at large rank distance is the normal state for an
  under-represented taxon, which is the case the tool exists to serve.

A `drop` value hand-written into the DB for rule 1, 2 or 6 is coerced to `flag`
on read, so a bad edit cannot enable an unsupported action.

Note that `is_nuclear_rrna` remains a title regex over a single hit
(`NUCLEAR_RRNA_TITLE_RE`, R/constants.R:54-57) and misses chromosome-level host
records. Rules 3-5 inherit that weakness. What makes rule 3 safe enough to drop
by default is not a better hit test but corroboration: a wrong drop needs the
regex to miss the record **and** barrnap and ITSx to find no rRNA and no ITS on
the same sequence, and those failure modes are uncorrelated. Repairing the regex
is Phase 3 work and will strengthen rules 3-5 without changing their structure.

### 3.2 Reason strings

Each rule owns its reason text, so `triage_reason` gains a one-to-one
relationship with `triage_rule` instead of being assembled ad hoc:

| key | reason |
|---|---|
| `not_annotated` | annotation has not run for this scaffold |
| `no_screen` | the GenBank screen has not run for this scaffold |
| `no_feat_non_rdna` | annotation found no rRNA or ITS feature, and the best GenBank match (%s) is not a nuclear rDNA record |
| `no_feat_rdna_hit` | annotation found no rRNA or ITS feature, though the best GenBank match (%s) is a nuclear rDNA record |
| `feat_non_rdna_hit` | best GenBank match (%s) is not a nuclear rDNA record |
| `odd_one_out` | top GenBank hit (%s) differs from this sample's other scaffolds (%s at %s) |

## 4. Schema v12

```sql
CREATE TABLE triage_opts (
  triage_opts            TEXT NOT NULL,
  rule_not_annotated     TEXT,   -- 'off' | 'flag'
  rule_no_screen         TEXT,   -- 'off' | 'flag'
  rule_no_feat_non_rdna  TEXT,   -- 'off' | 'flag' | 'drop'
  rule_no_feat_rdna_hit  TEXT,   -- 'off' | 'flag' | 'drop'
  rule_feat_non_rdna_hit TEXT,   -- 'off' | 'flag' | 'drop'
  rule_odd_one_out       TEXT,   -- 'off' | 'flag'
  PRIMARY KEY (triage_opts)
);
```

Seeded with a `default` row carrying the Default column of section 3. Plus
`assemble.triage_opts TEXT` defaulting to `'default'`.

Explicit columns rather than a JSON blob, matching `pre_opts`, `assemble_opts`,
`blast_opts`, `annotate_opts` and `curate_opts`: queryable, and a typo'd key
cannot silently disable a rule.

`assemblies.triage_rule TEXT` stores the key that fired, so the modal counts and
the Assemble-table filter work off a key rather than a regex over prose.

`auto_triage_sample()` and `triage_preview()` resolve settings through
`triage_opts_for(con, id)`, so the app's Apply and the WF1 `triage` process read
the same row. `inst/nextflow/bin/rp_triage.R` needs no change; resolution happens
inside `auto_triage_sample()`.

`TRIAGE_DEFAULTS` remains as the fallback when no row exists.
`min_depth_frac` is removed.

## 5. The bulk triage modal

Replaces the single prose paragraph and the depth numeric input.

```
  Rule                             Action              Decided   Matching
 -------------------------------------------------------------------------
  Annotation has not run           [Off|FLAG|  -  ]       1          1
  GenBank screen has not run       [Off|FLAG|  -  ]       0          0
  No annotations, non-rDNA hit     [Off|Flag|DROP ]       3          3
  No annotations but rDNA hit      [Off|FLAG|Drop ]       5          8
  Has annotations but non-rDNA hit [Off|FLAG|Drop ]       2         11
  Disagrees with sibling scaffolds [Off|FLAG|  -  ]       1          4

  Rules are applied top to bottom and the first one that matches decides
  whether to drop or flag the scaffold. "Decided" counts scaffolds this
  rule will drop or flag. "Matching" counts every scaffold that meets the
  rule, even if a higher-priority rule acts on it first.
```

- Each rule's `meaning` string renders under its label or as a tooltip.
- Drop is greyed on rules 1, 2 and 6, with a tooltip giving the reason.
- Changing any control recomputes the preview live.
- `Decided` sums to the unit total. `Matching` overlaps and does not.
- A rule set to Off has `Decided` 0 while `Matching` stays populated, which is
  the "what am I giving up" number.
- Retained from today, below the table: the summary line, the hand-decided and
  annotation-locked skip counts, the "also re-decide scaffolds I set by hand"
  checkbox, and the red "N sample(s) would end up with NO scaffolds kept"
  warning.
- When no rule is set to Drop, the modal says so plainly: Apply will only add and
  clear reasons.

`triage_preview()` gains `by_rule`: a data.frame of
`key, label, meaning, action, n_decided, n_matching`.

## 6. Depth becomes evidence, not a verdict

Depth leaves `triage_rules()` entirely, along with `min_depth_frac`.

Depth never dropped anything today (it was the fifth branch and only ever wrote
an amber reason), but the metric it used was
`median_depth / max(median_depth in sample)`, which is unsound in four separate
ways:

1. **The denominator is arbitrary.** Whatever scaffold happens to be deepest
   becomes the reference. Section 3.5(c) of doc 08 already calls this circular:
   in a fungal bloom, a myxozoan infection or a symbiont-rich coral the deepest
   unit is the contaminant, and every genuine host scaffold then falls below the
   floor.
2. **High depth is frequently collapse, not abundance.** GetOrganelle's graph
   sums depth over segments shared between host and contaminant conserved cores,
   and IGS subrepeats collapse into one inflated segment. The very collapse that
   generates the alternative paths inflates the denominator.
3. **A whole-scaffold median blends incomparable regions.** IGS and ITS are
   repeat-rich and collapse; the 18S and 28S cores are shared across taxa.
4. **Isomer paths split reads.** Two traversals of one molecule each receive
   roughly half the depth. That is path multiplicity, not evidence.

Depth stays visible in the Assemble-modal scaffold table as the existing median
depth column, plus a per-position sparkline rendered from `coverage.depths`,
downsampled server-side to roughly 100 points before it reaches the cell
renderer. The modal explainer gains one line stating why depth is shown and not
acted on.

Ideas considered and deferred rather than dropped, for a later phase that wants
a depth channel back: element-restricted depth over annotated 18S/5.8S/28S spans
only; depth uniformity (CV or MAD over the per-position vector, plus step
detection at element and HSP-window boundaries) which is level-independent and
therefore immune to the collapse problem; a host-concordant denominator per doc
08 section 3.5(c); an absolute read-support gate so a 300 bp scaffold at 10x
cannot produce a ratio at all; and a two-sided reading in which implausibly high
depth with no rRNA features is a collapsed-satellite signal.

## 7. Net behaviour change

With default settings the outcome is identical to today except for two fixes:

- A unit whose annotation crashed can no longer be dropped. It matches rule 1
  instead of the old `feats == 0 && nuc == 0` branch.
- "No features but the hit is rDNA" becomes its own visible rule with its own
  reason, instead of sharing one reason string with the crashed case.

Plus the removal of the depth reason, which flagged but never dropped.

## 8. Tests

- Each rule fires in isolation on a synthetic state frame.
- The 2x2 is exhaustive: all four channel combinations map to the intended rule
  or to `keep`.
- A unit with `annotate.n_features` NULL is never dropped, whatever its hit.
- Setting a rule to Off moves its units to the next matching rule, not to `keep`.
- Setting every combination rule to Off or Flag yields zero drops.
- `sum(n_decided)` equals the unit count; `n_matching >= n_decided` per rule.
- `drop` on rules 1, 2 or 6 in the DB is coerced to `flag` on read.
- Depth never appears in any reason string (regression pin).
- The existing "taxonomic distance alone never drops a scaffold" test still
  passes.
- A user decision remains sticky across an Apply with every rule enabled.

## 9. Out of scope

The identity envelope (doc 08 section 3.2), the `screen` verdict table (3.7),
the repaired `NUCLEAR_RRNA_TITLE_RE`, and the `screen` / `screen_tools` /
`cross_sample_screen` processes (3.9) all remain Phase 3. This design is
deliberately confined to making the existing evidence legible and controllable,
and to removing a channel that was never sound.

## 10. Known follow-ups

Both were found by the final pre-merge review, judged non-blocking, and left
deliberately. Both concern the annotation gate, which is the safety property
this design exists to establish, so they are recorded here rather than lost.

1. **`triage_rules()`'s optional `m` argument can silently mis-attribute.**
   The precomputed rule matrix is accepted without checking that it was built
   from the same frame. An `m` from a different frame of equal dimensions, or a
   row-permuted one, moves a verdict onto the wrong scaffold, including onto a
   crashed-annotation scaffold that must never be dropped. Only an undersized
   `m` errors; an `m` with stripped `dimnames` silently returns `keep` for
   everything. Nothing can trigger this today, since the sole caller builds `m`
   from the same `d` two lines earlier. Fix is one line before the loop:
   `stopifnot(nrow(m) == n, identical(colnames(m), triage_rule_keys()))`.

2. **The "no annotate row at all" test does not pin the gate it appears to.**
   Deleting `!isTRUE(ctx$annotated[i]) ||` from the `no_feat_non_rdna` gate
   leaves the whole suite passing. The test runs only the default action set,
   where `not_annotated` is first in evaluation order and shadows the drop rule
   by first-match. The gate is load-bearing only when `not_annotated` is set to
   `off`, which is a legal setting and is correct in the code. Adding
   `opts = list(not_annotated = "off")` to that test makes it bite.
