# 13 - The screen validation fixture set

Status: design agreed 2026-07-30. This is not a contamination channel. It is doc
08 section 5's **validation gate**, turned from a note into something that
exists: a labelled corpus in `tests/` against which the screen's verdicts and
its unfitted confidence weights can be measured rather than asserted.

## 0. The one-paragraph summary

Every number the screen decides with is a prior. Doc 11 section 5 says so of the
three confidence weights, doc 12 section 4 says so of `min_isomer_anchor_pid`,
and doc 08 section 3.2.1 says so of the identity envelope that was deferred
rather than built. Because of that, every rule and verdict added across 3a, 3b
and 3c ships as Flag, and auto-drop stays disabled. This sub-project builds the
labelled corpus that gate names, measures the screen against it, and prints the
numbers. It changes no default: arming auto-drop stays a human decision, taken
after reading the measurement rather than before.

## 1. What this validates, and what it cannot

**It validates** the taxonomy channel's verdict (`screen_channel_taxonomy()`),
the three confidence weights, the 0.6 uncertainty gate, and above all the
FALSE-POSITIVE RATE on clean host data. Doc 08 section 5 is explicit that a false
positive on an under-represented taxon is the failure that makes the tool
useless, so that is the number this exists to produce.

**It cannot validate `min_isomer_anchor_pid`.** Variant comparison runs only on
units the screen calls `on_target` (doc 12 section 2.1), and a spike-in is
`off_target` by construction, so a spiked unit is never compared to anything.
Validating the isomer threshold needs real multi-unit host samples carrying
genuine paralogs, which is the clade data of section 6. The 3c threshold
therefore remains an unvalidated prior after this sub-project, blocked on the
same data it was already blocked on. Recorded here so it is not assumed closed.

**It cannot validate the identity envelope**, which does not exist. See doc 08
section 3.2.1.

## 2. Decisions

1. Ground truth comes from DELIBERATE SPIKE-INS, so it is known by construction
   rather than by anyone's judgement. A fungal rDNA record is a fungal
   contaminant whatever it is queried against; no labelling opinion is involved.
2. The clean-host half of the corpus is REAL: the 19 units and 1833 hit rows of
   the reference project, all of which the screen currently calls `on_target`.
   That is the same 1833 doc 08 section 3.2.1 measures, so the fixture and that
   measurement share one corpus and cannot drift apart.
3. Hit distributions are CAPTURED ONCE by a dev script and committed. Tests are
   offline, deterministic and fast; they never reach the network. This is the
   same split as `backfill_hit_class()` versus the live parse path in 3a. The
   host half is captured from real core_nt results; the spike half from local
   curated databases, for the reasons in section 3.1.1.
4. The fixture stores hit rows in their REAL input shape, so the test exercises
   the same code path a live run does.
5. **No default changes.** The screen rule stays Flag, the export gate is
   untouched, `SCREEN_VERSION` does not move. Nothing in this sub-project alters
   what the pipeline decides.
6. A category the screen is KNOWN to be unable to catch is included and labelled
   as an expected miss, not omitted. Pinning a documented limitation is worth
   more than a corpus that flatters the tool.

## 3. Components

### 3.1 `dev/make_screen_fixture.R` (dev only, never run by tests)

Two jobs, both writing into `tests/testthat/fixtures/`.

**Clean-host capture.** Read `blast_hits` from a real project DB, join the
`screen` verdict and the sample's `lineage_taxids`, and emit every unit whose
stored verdict is `on_target`. No network.

**Spike-in capture.** For each row of a small category table, resolve an
accession with `esearch`/`efetch` (both present in `ribopilot:test`), fetch the
record, blast it LOCALLY against prebuilt curated rRNA databases, and run the
result through the pipeline's own `parse_blast_candidates()`.

#### 3.1.1 Amendment: local databases, not remote core_nt

The spec first said remote `core_nt`, matching the pipeline. That does not work:
NCBI runs BLAST COMPUTE as a shared queue and it stops answering under load.
Measured during implementation: one capture sat **85 minutes at 00:00:00 CPU**
on a query that normally takes minutes, and an independent probe from a separate
container hung too. Their FILE SERVING kept working throughout.

So the generator fetches `SSU_eukaryote_rRNA`, `LSU_eukaryote_rRNA` and
`16S_ribosomal_RNA` once over file serving and searches them locally. Zero NCBI
compute per spike, seconds instead of hours, and reproducible. Each tarball
bundles `taxdb.bt*`, which is what supplies `staxids` and therefore the lineages
the screen scores on, so no separate taxonomy download is needed. Cached outside
the repo at `~/.cache/ribopilot-blastdb`; never committed.

`update_blastdb.pl` is NOT used: the conda-packaged script cannot find curl,
leaves `$curl` uninitialised at its line 599, tries to exec `--user` as a
program, and reports success anyway. The generator uses curl and tar directly.

Two consequences of the curated subject pool, both real:

- **The `human` subclass is unreachable.** `screen_class_for()` returns `human`
  only when the LCA's own label is `Homo sapiens`. `SSU_eukaryote_rRNA` is type
  material, roughly one 18S per taxon, so a human 18S query hits hedgehog,
  opossum, rabbit, platypus and skink at similar identity and the LCA settles at
  Mammalia. `human_on_invert` therefore expects class `metazoan_other`.
  `tests/testthat/test-screen.R:161` covers the `human` branch as a unit; nothing
  covers it end to end.
- The two eukaryote databases date from **2022-12-05**, older than core_nt. The
  database list is recorded in every spike label's `note`.

The clean-host half remains real core_nt, so the fixture retains one corpus that
is genuinely what the pipeline sees.

Accessions are RESOLVED AT GENERATION TIME from a taxon plus a query, not
hardcoded in this document, and the resolved accession is written into the
fixture. Hardcoding an accession that is never verified is how a fixture starts
lying; recording what was actually fetched is how it stays honest.

The script prints what it fetched and what BLAST returned, and refuses to
overwrite an existing fixture without an explicit flag, because regenerating it
silently would erase the provenance of a number the tests assert against.

### 3.2 The fixture files

`tests/testthat/fixtures/screen_hits.csv.gz` - one row per HSP, in the shape
`blast_hits` stores and `screen_channel_taxonomy()` consumes:

```
unit, hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
```

Gzipped because the clean-host half alone is 300 KB of CSV and `read.csv()`
reads `.gz` natively.

`tests/testthat/fixtures/screen_labels.csv` - one row per unit:

```
unit, category, expected_verdict, expected_class, source, sample_taxids,
accession, taxdump_date, note
```

`sample_taxids` is in the fixture because the screen scores a unit's hits
AGAINST its sample's lineage; a hit distribution alone is not a scoreable input.
`taxdump_date` is recorded because lineages come from taxonkit and shift between
releases, so a future failure can be told from a real regression.

### 3.3 `tests/testthat/test-screen-validation.R`

Replays every fixture unit through `screen_channel_taxonomy()` and asserts both
the per-unit expectation and the aggregate floors of section 5.

## 4. Categories

A unit in this fixture is a HIT DISTRIBUTION paired with a HOST LINEAGE. The two
are independent: the distribution comes from BLAST, the host lineage is metadata.
One captured distribution can therefore appear twice against different hosts,
which is how the screen's actual boundary gets measured.

| unit | distribution | host lineage | expected verdict | expected class |
|---|---|---|---|---|
| `host_*` (19 real) | real fish scaffolds | fish | `on_target` | NA |
| `fungal_on_fish` | fungal rDNA | fish | `off_target` | `fungal` |
| `bacterial_on_fish` | bacterial 16S | fish | `no_reference` (BOUNDARY) | NA |
| `algal_on_fish` | plant rDNA | fish | `off_target` | `algal_plant` |
| `arthropod_on_fish` | arthropod rDNA | fish | `off_target` | `metazoan_other` |
| `human_on_fish` | human rDNA | fish | `on_target` (EXPECTED MISS) | NA |
| `human_on_invert` | the SAME human rDNA | octocoral | `off_target` | `metazoan_other` |
| `fungal_weak_on_fish` | top 2 accessions of the fungal capture | fish | `uncertain` | NA |

Fungal and algal_plant are the cross-kingdom cases. Arthropod is a metazoan from
a different phylum, which concordance does catch.

### 4.2 Amendment: bacterial is a boundary, not a catch

Measured: a bacterial 16S spike returns **`no_reference`**, not `off_target`, and
it is right to. Doc 10 classifies a 16S title carrying no eukaryotic marker as
`other` BY DESIGN, precisely so a bacterial record cannot pose as nuclear rDNA.
The taxonomy channel scores only `hit_class = 'rdna'` rows, so all 50 of the
spike's hits are filtered out and it has nothing to vote on.

**The screen's taxonomy channel cannot return `off_target` for a pure bacterial
contaminant.** That case is caught by the triage rules keyed on
`is_nuclear_rrna`, a different channel from the one this fixture validates. Doc
08 section 5 asks for a bacterial spike-in as though the screen would catch it;
it will not, and the fixture asserts the `no_reference` so that a future change
to `off_target` is noticed rather than welcomed unexamined.

### 4.3 Why a derived `uncertain` unit exists

No captured unit exercises the 0.6 gate: every one sits either well above it or
at confidence 0. `fungal_weak_on_fish` is the fungal capture reduced to its top
two accessions, which is a real contaminant that is resolvable but weakly
supported, and it lands at confidence 0.274. It is labelled `source = derived`.
Without it, removing the uncertainty gate entirely breaks no test.

### 4.1 Why human appears twice

`Homo sapiens` is Chordata. So is a ray-finned fish. Concordance is tested at
phylum rank or deeper, so a human contaminant on a FISH sample shares the phylum
slot, is judged host-concordant, and comes back `on_target`. Measured:
`on_target`, `host_concordant = 1`, confidence 0.82. Doc 08 section 5 asks for a
human spike-in as if it were a cross-kingdom case; on a chordate host it is not,
it is the blind spot.

The same distribution against an octocoral host lineage is `off_target`, because
Cnidaria and Chordata share only the shallowest slot. Pairing one capture with two
hosts holds the contaminant constant and varies only the host, so the boundary is
measured rather than argued. It also exercises the `human` subclass end to end,
which needs the LCA's own label to resolve to `Homo sapiens` and so cannot be
reached with authored placeholder lineage strings.

The invertebrate host lineage is a real rank-aligned taxid path resolved from the
taxdump already present in the test project. No new sequence data is needed: only
the host's lineage varies, and that is metadata.

The blind spot is narrower than doc 11 says, and getting it wrong here would
build a fixture that tests the wrong thing.

`host_concordant` asks whether the unit's LCA sits on the sample's path AT PHYLUM
RANK OR DEEPER (doc 11 section 4.4). A metazoan from a DIFFERENT phylum therefore
fails concordance and IS caught: `tests/testthat/test-screen.R:176` pins exactly
that, and doc 11 section 4.5 records the measurement behind it. A metazoan from
the SAME phylum shares the phylum slot, passes concordance, and comes back
`on_target`; `tests/testthat/test-screen.R:248` pins that too.

So the blind spot is same-phylum, not within-kingdom. **Doc 11 section 1's
sentence "a copepod and a fish share Metazoa" names the one metazoan case the
screen handles**, because a copepod is an arthropod and a fish is a chordate.
That sentence predates the phylum-or-deeper rule adopted later in the same
document and was never reconciled with it. Correcting it is a task of this
sub-project.

`metazoan_same_phylum` is included, and expected to be MISSED, so the real
limitation is pinned rather than rediscovered. Its expected verdict is
`on_target`: the screen does not decline, it actively judges the contaminant to
be host-concordant, which is worse than uncertainty and worth recording as such.

## 5. What the test asserts

**Zero false positives on `host`.** A unit in the `host` category that comes
back `off_target` fails the suite outright. This is the assertion the gate is
really about.

**Every catchable spike-in unit is caught, per unit, not in aggregate.**
Each of `fungal_on_fish`, `bacterial_on_fish`, `algal_on_fish`,
`arthropod_on_fish` and `human_on_invert` is asserted separately, so one category degrading to zero
cannot be hidden by the other four. An aggregate recall number is also printed,
but nothing is asserted on it alone: an aggregate floor is exactly what lets a
whole category silently stop working.

**Expected misses still miss.** If `human_on_fish` ever comes back `off_target`,
the test reports it. That is a capability gain and worth noticing,
so it is surfaced rather than tolerated silently.

**The confusion matrix is printed**, so the numbers are visible in the test
output rather than compressed to pass or fail.

An `uncertain` verdict on a spike-in is NOT counted as a catch. The screen
recording `uncertain` means it declined to decide, and counting that as success
would let the 0.6 gate be tuned to look good by refusing to answer.

## 5.1 Measured limits of the corpus itself

A fixture that is not honest about what it cannot catch is worse than no fixture.
Two limits, both found by mutating the scoring core and watching what the corpus
failed to notice.

**No unit exercises placeholder-taxon exclusion.** Doc 11 section 4.2 says
uncultured and environmental records "dominate core_nt 18S", but 0 of the 1833
host hit rows has a placeholder as its deepest lineage field: these are megablast
hits to well-curated fish references. So `is_placeholder_taxon()` is unexercised
by the corpus, and the curated spike databases contain no such records either.

**The 19 host units do not exercise doc 11 section 4.3's congener defect.**
Reintroducing the too-strict rule that doc 11 rejected, `host_concordant` from
`lca$taxid %in% samp_ok` alone, produces ZERO failures against this corpus. None
of these fish units has an LCA sitting off its sample's own path, so the case doc
11 measured at 0.790 does not occur here and a regression to that rule would pass
unnoticed. `tests/testthat/test-screen.R:129` pins it synthetically, so the SUITE
is covered; the CORPUS is not. Closing this needs a host unit whose hits are
dominated by same-genus relatives, which none of the ten samples provides.

**The phylum-or-deeper rule is enforced twice, so no single-point mutation breaks
it.** Both `samp[-1]` and `lca_path[-1]` in `screen_channel_taxonomy()` drop the
kingdom slot, and either alone is sufficient. Removing just one changes no
verdict; removing both makes `human_on_invert` concordant and the fixture catches
it by name. Redundancy rather than a gap, recorded so a future reader does not
delete one as dead code.

**The false-positive rate is measured on actinopterygians only.** Every host unit
is a ray-finned fish from one reference project. Doc 08 section 5 asks for the
rate on under-represented invertebrate clades precisely because that is where a
false positive is most likely and most damaging, and that is the clade half of
the gate, still blocked. A zero here is real but narrow: it says the screen does
not flag clean FISH data.

## 6. Out of scope

**The clade half of doc 08 section 5's gate.** That gate asks for validation
against the asteroidea, scyphozoa and octocorallia invert clades as well. Those
three read sets are not downloaded, their clade reference DBs are not built, and
no invert project has ever been run: `mapping_asteroidea.csv`,
`mapping_scyphozoa.csv` and `mapping_octocorallia.csv` each define a single
sample, and every existing test project is actinopterygian. Satisfying that half
requires read downloads, three reference DBs, three assemblies and three remote
BLAST runs before any labelling can begin, so it is a second deliverable blocked
on data rather than on code. Doc 08 section 5 is amended to record the split.

Arming auto-drop, or enabling the export gate. Both stay disabled. This
sub-project produces the measurement those decisions need; it does not take them.

Broadening what the screen can catch. `human_on_fish` is recorded as a miss,
not fixed. Fixing it needs the depth and cross-sample channels (doc 08
sections 3.5(c) and 3.6), both of which are blocked, and this sub-project must
not grow into them. The chimera detector (doc 21) now exists and its read-clip
and depth-step channels overlap with what would catch this case, but its effect
on `human_on_fish` has not been measured; that needs a pipeline run against
this fixture, not yet done.

`min_isomer_anchor_pid`, per section 1.

Fitting the confidence weights. This measures the shipped priors; it does not
optimise them. Refitting against a corpus this small would overfit to six
categories, and any change to the weights must bump `SCREEN_VERSION`.

## 7. Testing

- Every fixture unit gets its expected verdict, per category.
- The `host` category yields zero `off_target`.
- A spike-in scored `uncertain` is not counted as a catch.
- The fixture files parse, and every unit named in `screen_labels.csv` has at
  least one hit row in `screen_hits.csv.gz`, and vice versa. A label with no
  hits, or hits with no label, is a corrupt fixture and must fail loudly rather
  than be silently skipped.
- `expected_class` is asserted for each off-target category, so a regression in
  `screen_class_for()` is caught here too.
- The test reaches no network: asserted by construction, since it reads only the
  two fixture files.
- `human_on_fish` is missed and `human_on_invert` is caught, from ONE captured
  distribution, so the boundary of the screen's reach is stated as a measurement
  rather than as prose.

## 8. Documentation tasks

- Amend doc 08 section 5 to record that the gate splits into a spike-in half,
  delivered here, and a clade half blocked on data that does not exist locally.
- Correct doc 11 section 1: the screen cannot separate two metazoans OF THE SAME
  PHYLUM. The copepod-and-fish example is a different-phylum case the screen does
  catch, and `tests/testthat/test-screen.R:176` has pinned that since 3b.
- Record in doc 08 section 5 that its requested HUMAN spike-in is not a
  cross-kingdom case on a chordate host, and is caught only when the host is not
  a chordate. See section 4.1.
- Record in doc 12 that `min_isomer_anchor_pid` is still unvalidated after this
  sub-project, and why spike-ins cannot reach it.
