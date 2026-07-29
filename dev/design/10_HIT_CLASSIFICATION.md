# 10 - Three-state hit classification and consensus

Status: design agreed 2026-07-29. Sub-project 3a of the contamination plan in
`dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md`. Implements the classification
half of section 3.4. The identity envelope (section 3.2) moves to sub-project
3b, where the `screen` table consumes it and where element assignment is
available after reconcile.

## 0. The one-paragraph summary

`blast_screen.is_nuclear_rrna` decides whether a scaffold's best GenBank match
is a nuclear rDNA record, and three shipped triage rules read it. It is computed
from a literal-phrase regex over ONE hit's title. In the reference project 26 of
601 distinct hit titles score 0, and every one of them is a chromosome-level or
haplotype-resolved fish assembly, that is, a genuine host rDNA hit landing
inside a whole-genome record. The repair is not a broader regex: a chromosome
title is not evidence against rDNA, it is no evidence either way. This design
replaces the binary with a three-state classification, decides the per-unit flag
by a bitscore-weighted consensus rather than from hit number one, and backfills
both offline so existing projects are repaired without re-running the remote
BLAST.

## 1. What this classifier is and is not

It answers **"does this record describe nuclear rDNA"**. It does NOT answer
"is this the host". A fungal ITS record is correctly `rdna`; a human 18S record
is correctly `rdna`. Detecting that either is a contaminant is the taxonomy
channel's job in 3b, from rank distance and the identity envelope. Anyone
expecting this layer to catch a fungus has misread it, and the code comments
must say so.

## 2. Decisions

1. Titles classify into three states: `rdna`, `other`, `uninformative`.
2. Assembly-level records are `uninformative` and are excluded from the vote,
   not counted against rDNA.
3. The per-unit flag comes from a bitscore-weighted consensus over the top
   hits, not from hit number one.
4. Classification is a pure function of the stored title, so existing projects
   are repaired by an offline backfill, with no new remote BLAST.
5. Only the assembly class and the consensus change behaviour. A title matching
   nothing recognised stays `other`, exactly as today.

## 3. The classification

`classify_hit_title(title)` is pure and vectorised, returning one of `rdna`,
`other`, `uninformative`. Precedence matters and is fixed:

| # | Test | Result | Why it is at this precedence |
|---|---|---|---|
| 1 | organellar (`mitochondri`, `chloroplast`, `plastid`, `apicoplast`, `kinetoplast`) | `other` | A mitochondrial LSU record contains "large subunit ribosomal RNA". The veto must beat the rDNA phrase or every organellar assembly reads as the target locus. This is the single offline-independent signal that can say "this is a mitochondrial assembly, not nuclear rDNA". |
| 2 | `16S` or `23S`, and NOT also `18S`, `28S`, `5.8S` or an ITS phrase | `other` | Guards the broadened patterns below. "16S rRNA gene" must not match a generic "rRNA gene" pattern. A bacterial hit is real negative evidence and must stay negative. The eukaryotic exception exists because a few records name both, for example an operon spanning entry or a title listing several genes. |
| 3 | nuclear rDNA phrasing | `rdna` | See the pattern list below. |
| 4 | assembly-level (`chromosome`, `unplaced`, `scaffold`, `contig`, `whole genome shotgun`, `haplotype`) | `uninformative` | The title describes the container, not the locus. |
| 5 | anything else | `other` | Preserves today's default, so the only titles that change class are the assembly-level ones. |

The nuclear rDNA pattern broadens the current
`NUCLEAR_RRNA_TITLE_RE` (R/constants.R:54-57) to cover the phrasings GenBank
actually uses:

- `internal transcribed spacer`, `ITS1`, `ITS2`, `external transcribed spacer`
- `18S`, `28S`, `5.8S`, `5S` each followed by `ribosomal RNA` or `rRNA`
- `small subunit ribosomal RNA`, `large subunit ribosomal RNA`
- `SSU rRNA`, `LSU rRNA`
- `ribosomal RNA gene`, `rRNA gene` (safe only because rule 2 above has already
  removed 16S and 23S)
- `intergenic spacer` and `IGS` where accompanied by a ribosomal marker

**5S is classified `rdna`.** It is nuclear ribosomal DNA. Whether a given record
is part of the 45S repeat this tool targets, rather than a separate 5S array, is
a locus question that element assignment and the envelope answer in 3b, not
something a title can settle. Recorded here because the choice is arguable.

## 4. Storage and migration (schema v13)

`blast_hits` gains `hit_class TEXT`. `is_nuclear_rrna` keeps its INTEGER type
and its existing 1 and 0 meanings, so nothing that reads it today breaks. It
gains one new state: `NA` for a unit with no informative hits. Callers that
treat it as strictly binary must be checked, and triage is the only one, handled
in section 6. Per-hit `blast_hits.is_nuclear_rrna` is left exactly as it is and
becomes redundant with `hit_class`; it is kept rather than dropped because the
migration contract is additive and pre-v13 rows already carry it.

Because `classify_hit_title()` is a pure function of a column already stored,
the migration can repair history:

- `backfill_hit_class(con)` classifies every `blast_hits` row from its stored
  title and writes `hit_class`.
- It then recomputes `blast_screen.is_nuclear_rrna` by consensus for every unit
  that has `blast_hits` rows, leaving units without hits untouched.

That matters practically: the reference project's remote BLAST took thirteen
hours. The repair must not require repeating it. `blast_genbank` already skips
units already present in `blast_hits`, so a re-run would not recompute them
anyway.

## 5. The consensus

`consensus_hit_class(hit_class, bitscore)` returns `1L`, `0L` or `NA_integer_`.

- Consider hits within 10 percent of the top bitscore, taking at least 3 where
  that many exist. This mirrors `bitscore_weighted_lca()`
  (R/taxonomy_ranks.R:171), which exists for the same reason: one mislabelled
  GenBank accession must not decide a unit.
- Among those, drop `uninformative` hits entirely.
- Sum bitscore for `rdna` against bitscore for `other`. Higher total wins:
  `rdna` gives `1L`, `other` gives `0L`.
- A tie resolves to `1L`. A tie means the evidence is genuinely split, and the
  cost of a false 0 here is a scaffold flagged by rules 4 or 5, while the cost
  of a false 1 is only the absence of that flag.
- No informative hits at all gives `NA_integer_`.

`write_blast_screen()` (R/blast_genbank.R:195-211) stops using
`top$is_nuclear_rrna %||% is_nuclear_rrna_title(top$title)` and calls this
instead over the unit's full hit frame.

The reference candidate panel is unchanged: `parse_blast_candidates()`
(R/blast_genbank.R:117) keeps filtering to `rdna` records only. An assembly-level
record is a poor per-element reference, so `uninformative` hits must not become
reference candidates.

## 6. One consequence for triage

`is_nuclear_rrna` can now be `NA` for a scaffold that was screened perfectly
well and whose hits were all assembly-level. Triage rule `no_screen` fires on
`is.na()`, and its `meaning` string is already correct ("No usable core_nt hit
is recorded, so the scaffold cannot be judged against GenBank"), but its reason
string reads "the GenBank screen has not run for this scaffold", which would be
false. It becomes:

> no GenBank hit says whether this scaffold is nuclear rDNA

covering both "not screened" and "screened, nothing classifiable". No rule
ordering, action set or gate changes.

## 7. The fixture corpus

`tests/testthat/fixtures/hit_titles.csv`, columns `title`, `expected`, `source`,
`note`:

- All 26 real `uninformative` titles from the reference project. These are the
  regression that motivated the work and must be pinned individually.
- A stratified sample across the 575 real positives, covering each observed
  shape: ITS (241 distinct), 28S (237), 18S (69), and the bare
  small/large-subunit phrasing (28).
- Hand-authored rows for the categories the reference project contains none of,
  because it is ten clean fish samples: mitochondrial and chloroplast rRNA,
  prokaryotic 16S and 23S, protein-coding and mRNA records, and the adversarial
  case of an organellar record whose title also contains an rDNA phrase.

`source` records whether a row is real or authored, so a later reader can tell
measured coverage from imagined coverage. This corpus is the beginning of the
labelled fixture set that `08_MULTIUNIT_UX_AND_CONTAMINATION.md` section 5 asks
for as a gate on the envelope constants; it covers classification only, and 3b
must extend it with identity and rank-distance labels before the envelope ships.

## 8. Testing

- Every fixture row asserted against `classify_hit_title()`.
- The organellar veto beats an rDNA phrase in the same title.
- A 16S record does not match the broadened `rRNA gene` pattern.
- A unit whose informative hits are all `rdna` scores 1; all `other` scores 0.
- A unit whose hits are ALL assembly-level scores `NA`, not 0. This is the
  reference project's 26-title case.
- One mislabelled top hit cannot flip a unit whose consensus disagrees, tested
  by giving the outlier the top bitscore.
- A tie resolves to 1.
- `backfill_hit_class()` populates `hit_class` for pre-v13 rows and recomputes
  `blast_screen.is_nuclear_rrna`, leaving units with no hits untouched.
- The migration is idempotent: running it twice changes nothing the second time.
- Reference candidates still exclude `uninformative` records.

## 9. Out of scope

The identity envelope, the `screen` verdict table, the evidence channels, and
the `screen` / `screen_tools` / `cross_sample_screen` processes are 3b onward.
Classifying a record by taxonomy rather than by title is explicitly not this
layer's job, per section 1.
