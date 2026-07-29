# 11 - The screen verdict layer

Status: design agreed 2026-07-29. Sub-project 3b of the contamination plan in
`dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md`, implementing section 3.7 and
the taxonomy terms of channel (a) in section 3.5. Builds on 3a
(`dev/design/10_HIT_CLASSIFICATION.md`).

**The identity envelope (doc 08 section 3.2) is NOT in this sub-project.** It was
deferred after measurement, and the measurement is recorded at doc 08 section
3.2.1: the reference corpus has 1172 hits at class distance, 2 at phylum and
none beyond, so it cannot calibrate a bound that acts at phylum distance and
above; and genuine host 18S reaches 99.6 percent identity at class distance, so
no ceiling can sit below that without flagging real host data. Every term in this
design is therefore free of identity constants.

## 0. The one-paragraph summary

The machine's opinion about a scaffold currently exists only as a substring
inside `annotate.annotate_notes`, a column that `curate_one()` overwrites
wholesale and the modal textarea also writes to. This sub-project gives it a
home: a recomputable, versioned `screen` table holding one row per assembly
unit, populated from the taxonomic agreement of the unit's whole BLAST hit
distribution with the sample's own lineage. It needs no identity threshold,
because "do the hits agree on a taxon, and is that taxon the host's" is a
question about taxonomy alone.

## 1. What this layer decides, and what it cannot

It answers **"does this scaffold's hit distribution agree with the sample's own
lineage"**. A fungal, bacterial or plant contaminant on a metazoan sample fails
that outright, which is the case this catches well.

It cannot separate two metazoans, because a copepod and a fish share Metazoa and
the concordance test is against the sample's rank-aligned path. Doc 08 section
3.5 is explicit that the within-kingdom cases need the depth and cross-sample
channels, which are 3d and 3e. Nothing here should be described as catching a
copepod on a fish.

## 2. Decisions

1. One row per `(ID, path, scaffold)` in a new `screen` table, recomputable and
   versioned, never user-written.
2. Verdicts are `on_target`, `off_target`, `no_reference`, `uncertain`. Chimeric
   and paralogous_variant belong to 3d and 3e.
3. Confidence is computed from the channel's own evidence strength, since with
   one channel there is nothing to renormalise against. The frame is built so
   3c can add a second channel without reshaping the table.
4. The screen becomes a seventh triage rule, shipping as **Flag** with Drop
   available. Doc 08 section 5's validation gate says auto-drop must not ship
   enabled before a labelled fixture set exists, and section 3.2.1 records that
   it does not.
5. Uncultured and environmental records are excluded from every computation
   while remaining visible in the hit list.

## 3. Schema v14

```sql
CREATE TABLE screen (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  screen_verdict    TEXT,
  screen_class      TEXT,
  screen_confidence REAL,
  screen_reason     TEXT,
  screen_version    INTEGER,
  lca_taxid         INTEGER,
  lca_taxon         TEXT,
  lca_rank          TEXT,
  lca_n_hits        INTEGER,
  host_concordant   INTEGER,
  margin            REAL,
  time_stamp        INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);
```

Plus:

- `samples.reference_support INTEGER` - the minimum rank distance across the
  sample's best on-target hits, computed from modal-taxon units only and
  excluding placeholder taxa. When it is 4 or worse the UI must say so
  ("nearest reference is at class distance") rather than staying silent.
- `triage_opts.rule_screen_off_target TEXT` - the new rule's action, seeded
  `'flag'`.
- `SCREEN_VERSION <- 1L`, an internal constant. `screen_one()` rewrites a row
  whose stored `screen_version` is lower, so a scoring change can force a
  recompute without a schema migration.

The `screen` table is created empty by the migration. It is NOT backfilled: the
computation needs `samples.lineage_taxids` and a taxdump, neither of which the
migration can assume, and an empty row is honestly distinguishable from a
computed one. `screen_one()` fills it on the next pipeline run or on demand.

## 4. The channel

`screen_one(con, id, path, scaffold, taxdump_dir = NULL)` is the DB driver.
`screen_channel_taxonomy(hits, sample_taxids)` is the pure core and is what the
tests exercise.

Reads the unit's `blast_hits` rows, then:

1. **Exclude placeholders.** Drop any hit whose lineage fails
   `is_placeholder_taxon()` (R/taxonomy_ranks.R:45). Uncultured and
   environmental records dominate core_nt 18S and resolve to superkingdom only,
   which every rank rule reads as maximum taxonomic distance.
2. **Weighted LCA.** `bitscore_weighted_lca()` over the surviving
   `hit_class = 'rdna'` rows, giving `lca_taxid`, `lca_taxon`, `lca_rank` and
   `lca_n_hits`.
3. **Concordance.** `host_concordant` is 1 when `lca_taxid` appears in the
   sample's rank-aligned taxid path (`parse_rank_taxids(samples.lineage_taxids)`),
   0 when it does not, NA when either side is unresolvable.
4. **Margin.** Best bitscore among concordant hits minus best bitscore among
   discordant hits. NA when either side has no hits.

   Note this uses a SECOND, per-hit notion of concordance, distinct from
   `host_concordant` above. `host_concordant` asks whether the unit's single
   weighted LCA sits on the sample's path; the margin asks the same question of
   each hit individually, so that the two best-supported opposing hits can be
   compared. A hit is concordant here when its own rank-aligned taxid path
   contains the sample's DEEPEST RESOLVED taxid, that is the last non-NA entry
   of `parse_rank_taxids(samples.lineage_taxids)`. Keep the two named apart in
   the code; conflating them silently changes what the margin measures.

## 5. Confidence

Three sub-scores in [0, 1], combined as a weighted mean:

| term | definition | weight |
|---|---|---|
| `support` | `min(1, lca_n_hits / 10)` | 0.30 |
| `resolution` | rank depth of `lca_rank` over `RRNA_LINEAGE_RANKS` (`species` 1.00, `genus` 0.86, `family` 0.71, `order` 0.57, `class` 0.43, `phylum` 0.29, `kingdom` 0.14, `none` 0) | 0.30 |
| `separation` | `min(1, abs(margin) / top_bitscore)`; **1 when every surviving hit falls on one side of the concordance test**, since nothing opposes the winner; 0 when margin is NA for any other reason | 0.40 |

**These weights are priors, not measurements.** That is the same caveat that
deferred the envelope, and it is why the screen rule ships as Flag. They live in
`rrna_params()`'s `screen` sublist so a clade can retune them, and
`screen_reason` always prints each term with its value and its weight, so a
curator can see when a verdict rested on a single term. Any future change to
them must bump `SCREEN_VERSION`.

## 6. Verdicts

Evaluated in this order:

| verdict | condition | `screen_class` |
|---|---|---|
| `no_reference` | nothing resolvable, which includes `host_concordant` being NA; reason code `no_hits`, `no_taxdump`, `unresolved_sample_lineage` or `unresolved_hit_lineages` | NA |
| `uncertain` | resolvable but `screen_confidence < 0.6` | NA |
| `off_target` | `host_concordant = 0` at confidence >= 0.6 | see below |
| `on_target` | `host_concordant = 1` at confidence >= 0.6 | NA |

`no_reference` renders neutral grey, gates nothing, and replaces
`contaminant_unchecked`. It is not a problem state: for an offline project or a
sparse clade it is the normal one, and painting it as needs-attention makes
every sample in the tool's core use case look broken.

`screen_class` for an off-target unit comes from the LCA's kingdom-rank name:
`Bacteria` -> `bacterial`, `Archaea` -> `archaeal`, `Fungi` -> `fungal`,
`Viridiplantae` -> `algal_plant`, `Metazoa` -> `metazoan_other`, anything else
eukaryotic -> `other_euk`. A `Homo sapiens` LCA is subclassed `human`, because a
human handling contaminant is worth naming separately even though it is metazoan.

Gating: a verdict below confidence 0.6 is recorded as `uncertain` with all of
its evidence still stored. The 0.8 threshold governs only whether the UI offers
a drop; it does not change what is stored.

## 7. Pipeline

A `screen` process, per sample, taking `blast_genbank.out`. That is its only
real dependency: without the envelope it needs no annotations and no coverage.
It carries `errorStrategy 'ignore'` so a taxonomy failure cannot sink a sample,
and `cache false` like the other R processes.

`triage` gains a join so a sample is triaged only once both its curation and its
screen exist:

```
    screen(blast_genbank.out)
    triage(curate.out.map { id, upath, uscaf -> tuple(id, "${upath}.${uscaf}") }
        .groupTuple()
        .map { id, _units -> id }
        .join(screen.out))
```

`inst/nextflow/bin/rp_screen.R` mirrors `rp_triage.R`: read the db path and
sample id, call `screen_sample()`, print a one-line summary.

## 8. Triage integration

A seventh rule, `screen_off_target`:

- key `screen_off_target`, label "GenBank taxonomy disagrees with the sample"
- meaning: "The scaffold's BLAST hits agree on a taxon that is not on this
  sample's lineage."
- fires when `screen_verdict = 'off_target'`
- actions Off / Flag / Drop, default **Flag**
- placed after the three combination rules and before `odd_one_out`, which tests
  the same idea more weakly (agreement among siblings rather than against the
  sample's own lineage)

`unit_triage_state()` gains `screen_verdict`, `screen_class` and
`screen_confidence` from a LEFT JOIN on `screen`. A unit with no screen row
leaves them NA and the rule cannot fire, which is correct: absence of a screen is
not evidence.

## 9. Testing

- `screen_channel_taxonomy()` pure tests: a concordant unit, a discordant one, a
  unit whose hits are all placeholders, a unit with one resolvable hit, and a
  unit whose sample lineage is unresolvable.
- Each confidence term independently, plus that the weights sum to 1 and that
  `screen_reason` names every term with its value and weight.
- Every verdict path including each of the four `no_reference` reason codes.
- The 0.6 gate: a resolvable but weakly supported unit is `uncertain`, not
  `off_target`, and its evidence columns are still populated.
- `screen_class` for each subclass, including that a `Homo sapiens` LCA is
  `human` rather than `metazoan_other`.
- `SCREEN_VERSION` forces recompute: a row stored at a lower version is
  rewritten; one at the current version is left alone.
- Triage: the new rule fires only on `off_target`, ships as Flag, sits in the
  documented position, and a unit with no screen row is untouched by it.
- A regression fixture extracted from the reference project: several real units'
  hits, asserting every one comes back `on_target` and host-concordant. That is
  the only real labelled outcome available, and it guards against a change that
  starts calling clean fish scaffolds contaminated.

## 10. Out of scope

The identity envelope and everything that depends on it; the tool-kingdom
channel (3c); the positional, chimeric and depth channels (3d); cross-sample
screening and `paralogous_variant` (3e); and splitting `annotate_notes` into
`curate_summary` plus a durable user field, which doc 08 section 3.7 lists
alongside this work but which touches the annotate UI rather than the screen.
