# 12 - Intragenomic rDNA variants: isomers and paralogs

Status: design agreed 2026-07-29. Sub-project 3c of the contamination plan in
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

`min_isomer_anchor_pid` is the ONLY new constant. It lives in `rrna_params()`'s
`screen` sublist so a clade can retune it, and the reason string always prints
the measured identity so a curator sees the number rather than only the label.

It is a prior, like the confidence weights, and the same caveat applies. It is
however a much safer prior than the identity envelope was, because it compares
two sequences from ONE organism rather than across taxa: the question "are these
two anchors the same sequence" has a far tighter natural answer than "how
similar may two different organisms be".

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

## 8. Salvage from the retired 3c

Two fixes that stand on their own, because the pipeline pays for the computation
and discards the answer:

- `parse_barrnap()`'s name map (R/annotate.R) recognises only `18S_rRNA`,
  `5_8S_rRNA`, `5.8S_rRNA` and `28S_rRNA`, so a 16S, 23S or 5S call is silently
  dropped. Extend the map. 5S is a normal nuclear metazoan gene and must be
  treated as eukaryotic, not as a prokaryotic signal.
- `parse_itsx()` reads only the coordinate line, while ITSx writes its winning
  taxonomic profile to disk. Harvest it.

Neither feeds a verdict in this sub-project. They stop information being thrown
away, and they are what a later channel would need if one is ever built.

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
- `parse_barrnap()` returns 16S, 23S and 5S rather than dropping them, and 5S is
  classified as eukaryotic.
- Grouping: a three-unit isomer set forms one group; a mixed set forms one
  isomer group and one paralog group.

## 10. Out of scope

Annotating the IGS as a first-class element (see 4.1). The identity envelope and
anything depending on it. The positional and chimeric channel, which doc 08
section 3.5(d) requires the envelope for. Cross-sample screening and index
hopping, which need `samples.batch` and do not exist yet. Any change to triage.
