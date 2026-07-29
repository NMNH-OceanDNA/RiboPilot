# 12 - Intragenomic rDNA variants: isomers and paralogs

Status: IMPLEMENTED 2026-07-29 on branch `contamination-phase3c`. Plan:
`dev/design/plans/2026-07-29-intragenomic-variants.md`. Five statements in this
design did not survive contact with the code; each is amended in place below
(sections 4.0.1, 4.0.2, 5.1, 8 and 10) rather than quietly implemented
differently.

Design agreed 2026-07-29. Sub-project 3c of the contamination plan in
`dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md`, implementing section 3.8.
Builds on 3b (`dev/design/11_SCREEN_VERDICT_LAYER.md`), which is what makes it
computable.

It takes the 3c slot because the original 3c was retired; the remaining
sub-projects renumber accordingly, so 3d is now the positional and chimeric
channel and 3e is cross-sample screening and index hopping.

**The original sub-project 3c is retired, not deferred.** Its justification was that the
tool-kingdom channel is the only one working fully offline, and users never run
RiboPilot offline. The reasoning is recorded at doc 08 section 3.5, channel (b).
The two salvageable pieces of it are folded into this sub-project, section 8.

## 0. The one-paragraph summary

A sample routinely yields several assembly units, and today the curator is asked
to adjudicate between them with no information about why they differ. Two very
different situations look identical on screen. Some units are the same molecule
assembled two ways, differing only in how many IGS subrepeats the graph
traversed; picking between those is a decision that does not matter. Others are
genuinely different rDNA copies from one individual, which is documented in
Octocorallia, Scyphozoa and several other groups; there the choice determines
what GenBank receives. This sub-project separates the two and says which is
which.

## 1. Why 3b is the prerequisite

The hard part of section 3.8 was never the arithmetic, it was knowing that the
units being compared belong to the same organism. Comparing two scaffolds tells
you they differ; it does not tell you whether that is intragenomic variation or
a contaminant.

3b's `screen_verdict` answers that. Restricting the comparison to units the
screen calls `on_target` means both are host-concordant, so **any difference
between them is intragenomic by construction**. No threshold has to separate
"paralog" from "contaminant", because that separation already happened. That is
why this design needs only one constant where section 3.8 implied several.

## 2. Decisions

1. Only units with `screen_verdict = 'on_target'` are compared. Everything else
   is a contamination question the screen already owns.
2. Comparison is over the annotated ANCHORS only (18S, 5.8S, 28S), never whole
   scaffolds. The spacers are what varies in both cases and so cannot
   discriminate between them.
3. The finding is ORTHOGONAL to `screen_verdict` and gets its own columns. A
   paralogous variant is emphatically `on_target`.
4. Neither class is a defect. No triage rule is added, and nothing is flagged.
5. Export records the finding as submission metadata; it does not gate.

## 3. What is compared, and how

For each sample, take its `on_target` units. For each pair:

1. Build each unit's anchor string by concatenating its annotated `18S`, `5.8S`
   and `28S` spans in element order, from `annotations` and
   `assemblies.sequence`. A unit missing all three anchors is not comparable.
2. Align the two anchor strings with `pairwise_nt(a, b, type = "overlap")` and
   read the identity with `pwalign_pid()`, both already in
   `R/blast_ref_align.R` and both already depending on Biostrings and pwalign.
   No new dependency and no external aligner.
3. Record the identity, and the difference in TOTAL scaffold length.

## 4. Classification

| anchors | total length | `variant_class` | what the curator is told |
|---|---|---|---|
| identity >= `min_isomer_anchor_pid` | differs | `isomer` | these differ only in spacer copy number, pick either |
| identity < `min_isomer_anchor_pid` | any | `paralog` | this sample carries several rDNA copies; which you lock determines what GenBank receives |
| only one on-target unit in the sample | - | `singleton` | nothing |
| anchors absent or unannotated | - | `NA` | not comparable |

`min_isomer_anchor_pid` is the ONLY new constant, and its value is **99.5**
percent, on the 0-100 scale `pwalign_pid()` returns. Across roughly 4 kb of
concatenated 18S+5.8S+28S that tolerates a handful of assembly errors without
admitting a real paralog. It lives in `rrna_params()`'s `screen` sublist so a
clade can retune it, and the reason string always prints the measured identity so
a curator sees the number rather than only the label.

#### 4.0.1 Amendment: the `screen` sublist did not exist, and 3c built it

This design and doc 11 section 5 both sited the screen's tunables in
`rrna_params()`'s `screen` sublist. Neither built it. 3b shipped `SCREEN_WEIGHTS`
and `SCREEN_MIN_CONFIDENCE` as top-level constants in `R/screen.R`, which no
ruleset could reach, so the claim that "a clade can retune them" was false for
the whole of 3b.

3c built the sublist. `SCREEN_PARAMS_DEFAULT` in `R/params_rrna.R` now holds the
three weights, `min_confidence` and `min_isomer_anchor_pid`, and every ruleset
carries it. `SCREEN_WEIGHTS` and `SCREEN_MIN_CONFIDENCE` survive as views on it
so a clade override and the shipped default cannot drift apart. Reaching a
sample's ruleset needed the resolution logic that was inline in `curate_one()`,
so it was lifted verbatim into `resolve_ruleset()` and is now shared.

No value changed. Every ruleset's `screen` sublist is currently identical to the
default, so no stored verdict shifts because of the move.

#### 4.0.2 Amendment: classification is on identity alone

The table above gives `isomer` the condition "identity at or above the threshold
AND total length differs", and gives no row at all for "identity at or above the
threshold AND total length identical". Since the only two classes are `isomer`
and `paralog`, and identity alone already decides between them, the length column
cannot change the answer and the missing row cannot be filled without inventing a
third class.

Implemented on identity alone. Two units with identical anchors and identical
total length are the same molecule assembled twice, for which "pick either" is
exactly the right message. The length difference is recorded in the reason string
as evidence rather than used as a test.

It is a prior, like the confidence weights, and the same caveat applies. It is
however a much safer prior than the identity envelope was, because it compares
two sequences from ONE organism rather than across taxa: the question "are these
two anchors the same sequence" has a far tighter natural answer than "how
similar may two different organisms be".

#### 4.0.0 What 99.5 means in absolute terms, measured

Over a realistic anchor total of 5460 bp (18S 1800 + 5.8S 160 + 28S 3500), with
substitutions only:

| substitutions | identity | class |
|---|---|---|
| 0 | 100.00 | isomer |
| 5 | 99.93 | isomer |
| 20 | 99.74 | isomer |
| 30 | 99.63 | isomer |
| 60 | 99.16 | paralog |
| 200 | 97.18 | paralog |

So 99.5 splits at roughly 27 differences across the three anchors. Anything
tighter than that is read as one molecule assembled twice.

Cost, same measurement: about 0.3 s per pair, and 0 s when the anchors are
identical, because that case short-circuits without aligning. A sample with k
comparable units pays k(k-1)/2 alignments, and k is typically 2 to 5.

### 4.1 Why the IGS is inferred rather than measured

RiboPilot's feature model has no IGS element. `parse_barrnap()` and
`parse_itsx()` between them produce 18S, ITS1, 5.8S, ITS2 and 28S, and nothing
covers the region between one repeat's 28S and the next repeat's 18S.

So the isomer test is "anchors identical, and the unannotated remainder differs
in length", not a direct measurement of subrepeat copy number. That is honest
about what is known and is enough to separate the two cases, but it means the
message says "differ only in spacer copy number" on the strength of an
inference. Annotating the IGS explicitly is a larger piece of work and is
deliberately not attempted here; if it is ever done, this test should be
rewritten to use it.

## 5. Grouping

Pairwise identities are reduced to groups per sample: units whose anchors are
mutually at or above `min_isomer_anchor_pid` form one isomer group; units below
it form separate paralog groups. `variant_group` is a stable label within the
sample, and `variant_identity` is each unit's anchor identity to its group's
representative, which is the longest member. The representative's own
`variant_identity` is NA, since it has nothing to compare against, and so is a
`singleton`'s: in both cases there is no second sequence, and recording 100
would assert a measurement that was never made.

Grouping is transitive-by-threshold and can in principle chain, where A matches
B and B matches C but A does not match C. When that happens the group is still
formed, and the recorded `variant_identity` makes the spread visible rather than
hiding it. Chained groups are expected to be rare at the identities involved and
are not worth a clustering algorithm.

### 5.1 Amendment: `variant_class` is a per-sample property

Section 4's table is written per PAIR. Reduced to one row per unit there is a
case it does not decide: a sample holding an isomer pair A, B plus a distinct
copy C. A and B are mutually isomeric, but the sample carries paralogs.

Implemented so that **every** comparable unit reads `paralog` once the sample has
two or more groups, A and B included. Telling the curator of A "pick either"
would be false while C exists, and the paralog message, "which you lock
determines what GenBank receives", is true of A, B and C alike. `variant_group`
is what tells the curator that A and B are interchangeable, which is the
information the isomer label would have carried.

The consequence is that `variant_class` is constant across a sample's comparable
units and is partly redundant with the group count. That redundancy is deliberate:
the class is what the UI and the export summary read, and it must not say
something false about any single unit.

## 6. Schema v15

Three columns on `screen`, which already carries one row per unit and a
version-gated recompute:

```sql
variant_class    TEXT,   -- isomer | paralog | singleton | NULL
variant_group    TEXT,   -- group label, unique within the sample
variant_identity REAL    -- anchor identity to the group representative
```

`SCREEN_VERSION` moves to `3L`, so existing rows recompute. No new table: the
finding is per unit, the recompute path exists, and a second table would need
its own version gate kept in step with the screen's.

## 7. Surfaces

**App.** The Assemble scaffold table gains a variant column. An isomer group
renders neutrally, because it is telling the curator a decision does NOT matter.
A paralog group renders as information, not as a warning: it is real biology,
and colouring it like a defect is what this sub-project exists to stop.

**Export.** `export_bundle()` records, for each exported scaffold, that it is
one of N intragenomic variants detected in its sample and the anchor identity
between them. It does not gate. This answers doc 08's open question 7 without
imposing a submission policy: the fact is recorded, the choice stays the
curator's.

**Triage.** Unchanged. No rule is added, and neither class contributes to any
verdict.

## 8. Salvage from the retired 3c: DEFERRED, not built

Both items were specified as small parser fixes. Neither is. They were measured
against the code during 3c and deferred, and this section records the
measurements so the next worker does not rediscover them.

**Extending `parse_barrnap()`'s name map retains nothing on its own.** The map at
`R/annotate.R:32` recognises only `18S_rRNA`, `5_8S_rRNA`, `5.8S_rRNA` and
`28S_rRNA`, so a 16S, 23S or 5S call is indeed dropped. But
`reconcile_features()` selects elements by literal name, one at a time
(`R/annotate.R:110-117`: `pick(bar, "18S")`, `pick(bar, "5.8S")`,
`pick(bar, "28S")`), and `reconcile_one()` writes only what `reconcile_features()`
returns (`R/pipeline_steps.R:139-149`). A new element emitted by the parser is
therefore discarded one function later and never reaches the `annotations` table.

Retaining it needs somewhere to put it, which is a schema decision rather than a
parser change, and every candidate has reach: writing non-target calls as ordinary
`annotations` rows puts them into the export feature table, the app's element
strip, `annotate.missing`, and curation's per-element length rules
(`RRNA_LENGTHS_METAZOA`), none of which have a concept of a non-target rRNA. When
this is picked up, decide the store first.

5S remains correct as specified: it is a normal nuclear metazoan gene and must be
treated as eukaryotic, not as a prokaryotic signal.

**The ITSx taxonomic profile never reaches R.** `parse_itsx()` reads only the
coordinate line, and that is all it can read. Both Nextflow processes that run
ITSx declare a single output file: `path("${id}.itsx.positions.txt")` at
`inst/nextflow/modules/annotate.nf:12`, and `path("unit.itsx.positions.txt")` at
`inst/nextflow/modules/pipeline_r.nf:149`. Every other file ITSx writes stays in
the work directory and is discarded with it, so no parser change can reach the
profile.

Harvesting it needs an output declaration in both processes, a new argument
threaded through `rp_reconcile.R` (which today takes exactly `<db> <id> <path>
<scaffold> <gff> <pos>`), and a column to hold the result. The all-in-R path
`run_itsx()` (`R/annotate.R:197`) is cheaper to fix, since it controls its own
output prefix and merely ignores the other files, but fixing only that path would
leave the pipeline and the convenience path disagreeing about what is stored,
which is the divergence the shared-DDL convention exists to prevent.

Neither item feeds a verdict, in this sub-project or in any planned one, which is
why deferring both costs nothing now.

## 9. Testing

- Anchor extraction: a unit with all three anchors, one with two, one with none.
- Two units with identical anchors and different total length are `isomer`.
- Two units with anchors differing beyond the threshold are `paralog`, and both
  remain `screen_verdict = 'on_target'`, which is the orthogonality this design
  turns on.
- A sample with one on-target unit is `singleton`.
- An `off_target` unit is never compared and its variant columns stay NA.
- A unit whose anchors are unannotated is NA, not `singleton`.
- The measured identity appears in the reason string, so a curator sees the
  number.
- `SCREEN_VERSION` at `3L` forces recompute of rows stored at `2L`.
- Export records the variant fact for a paralog and does NOT gate: a paralogous
  unit still exports.
- Grouping: a three-unit isomer set forms one group; a mixed set forms one
  isomer group and one paralog group.

Added during implementation, because the list above did not name them:

- A circularly permuted unit, where 28S precedes 18S by coordinate, produces the
  same anchor string as an unpermuted one. Ordering the concatenation by `start`
  rather than by element makes one molecule read as a paralog of itself.
- An all-minus unit is reverse-complemented before comparison, and a
  MIXED-strand unit is not comparable at all. Concatenating a sense and an
  antisense span produces a string that means nothing, and mixed strands are a
  chimera signal rather than a paralog.
- The group representative is the LONGEST member. A test whose input already
  happens to list the longest unit first passes for the wrong reason.
- Rerunning `screen_sample()` UNFORCED must not stack duplicate variant clauses
  in `screen_reason`. Forced, `screen_one()` rewrites a fresh reason first, so a
  forced rerun cannot exercise the strip and a test that uses `force = TRUE`
  passes while the strip is broken. This one was found by mutation and was a live
  defect: `VARIANT_REASON_SEP` contains `|`, which as a regex is alternation, so
  `sub()` matched a bare space and cut a word out of the verdict's own reason.
- `variant_identity` is asserted in the STORED COLUMN, not only in the reason
  string. Building the clause from the in-memory value and writing NULL to the
  column leaves a reason-only test green.
- `n_variants` in the export summary counts COMPARABLE units, so a fixture must
  contain a screen row with no variant class or the count is indistinguishable
  from the row count.
- The v15 migration must be tested against a screen table built at the V14
  SHAPE. A v8 fixture migrates by creating the table from the current
  `SCREEN_DDL`, which already has the new columns, so `add_col()` is a no-op
  there and the path a real v14 project takes goes untested.

## 10. Out of scope

Annotating the IGS as a first-class element (see 4.1). The identity envelope and
anything depending on it. The positional and chimeric channel, which doc 08
section 3.5(d) requires the envelope for. Cross-sample screening and index
hopping, which need `samples.batch` and do not exist yet. Any change to triage.

Deferred during implementation, with the reasoning in section 8: extending
`parse_barrnap()`'s name map, which needs a store decided first, and harvesting
the ITSx taxonomic profile, which needs Nextflow output plumbing.

Also not built: a version guard on the variant pass. `screen_variants()` recomputes
every time `screen_sample()` runs. The only caller is the per-sample `screen`
pipeline process, once per run, and identical anchors short-circuit to 100 without
aligning, so the common isomer case costs nothing. A guard keyed on "every
on-target unit already has a class" would also be wrong, because a unit with no
annotated anchors legitimately keeps a NULL class forever and would defeat it.
