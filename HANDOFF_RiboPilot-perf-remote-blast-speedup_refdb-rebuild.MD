# Handoff: RiboPilot, branch `perf/remote-blast-speedup`

Rewritten 2026-08-17, updated 2026-08-18. Nothing has been pushed. **The tree is NOT clean**: four
files carry uncommitted work from 2026-08-18 (see State).

## State

- Branch `perf/remote-blast-speedup`, **106 commits** ahead of `origin/main`, 110 from the merge base
  `f61912b`, linear, **no merges**.
- HEAD `69fd340` "fix: key the stage markers on the panel definition, not the list alone".
- **Working tree NOT clean. Nothing is committed from 2026-08-18.** Modified: this handoff,
  `dev/refdb/build_refdb.sh`, `dev/refdb/lib.sh`, `tests/testthat/test-refdb-build.R`. Three separate
  pieces of work, all verified, none committed:
  1. the re-review's one **Critical**, fixed (stage_fetch now refuses a list harvested under a
     different panel definition; `rp_marker_present` added to `lib.sh`);
  2. the **`efetch -input` bug**, fixed, plus the stub that let it through;
  3. the **`requested` count semantics**, changed by maintainer ruling (see Contracts).
  **Nothing pushed. Do not push without asking.**
- Suite: **3336 passed / 0 failed / 41 skipped** by default; **3591 passed / 0 failed / 0 skipped**
  with `NOT_CRAN=true`. Nine warnings are pre-existing and expected (tests that deliberately exercise
  warning paths: missing decoy, whitespace decoy dir, corrupt ref store). Do not chase them.

Those counts require real `makeblastdb`, `blastdbcmd`, `seqkit` and `update_blastdb.pl` on PATH. They
live in `/home/dmacguig/miniconda3/envs/ribopilot/bin`; **symlink exactly those four**. Without them
the tool-dependent tests skip cleanly. **Do not put the whole conda env on PATH**: `barrnap` and
`pybarrnap` then become visible and `test-rrna-kingdom.R` fails 4 expectations, pre-existing and
unrelated to this branch.

Run the suite in the foreground, never backgrounded, both ways:

```
Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); r <- testthat::test_dir("tests/testthat", reporter="silent"); d <- as.data.frame(r); cat("failed:", sum(d$failed), "error:", sum(d$error), "skipped:", sum(d$skipped), "passed:", sum(d$passed), "\n")'
```

## THE IMMEDIATE NEXT ACTION

**Commit the 2026-08-18 work, then finish the panel build.** The re-review that this section used to
ask for has been done and is closed; see below.

The scoped re-review of `f3fd2e8..69fd340` ran (six lenses, every finding attacked by two adversarial
verifiers; 34 findings verified, 22 survived, 12 refuted). Verdict was **NOT READY**, on one Critical.
The maintainer set a **stopping rule** at that point: that review was the last, only a Critical gets
fixed, everything else is recorded here as a residual. **The Critical is fixed. There is no further
review round.**

The Critical, for the record: `rp_panel_key` folds the CURRENT query with the accession-list hash, so
recomputing it inside `stage_fetch` proved nothing about which panel definition the list was
*harvested* under. `stage_accessions` is the only stage that knows, and it records that in a marker
**neither `stage_fetch` nor `stage_build` ever read**. Complete a build under query A, edit the query
to B, run `--stage fetch`: the list is unchanged so the resume key matches, every chunk is complete so
nothing is requested, and the stage re-stamps panel A's chunks with key(B). `stage_build`'s
precondition then passes and `blastdb/<NEW_VER>` ships panel A's records under a manifest quoting
query B, exit 0. The full-run path escaped only incidentally, because `stage_accessions` re-harvests
first. The amendment's own new test asserted that laundering as correct.

The three questions the re-review was asked to settle: **Q1 OVERRULED** (`EXPECTED_COUNT` selects no
record and never reaches the manifest, so it buys no protection while arming a destructive re-harvest;
see Residuals), **Q2 confirm-with-caveat** (the caveat *was* the Critical, now closed), **Q3
confirm-with-caveat** (the "requests nothing" claim verified; see Residuals).

Then the build was run for the first time, and found two more defects that no amount of reading had
caught. Both are fixed; both are in Traps. The panel build is the remaining work.

The ledger, `.superpowers/sdd/2026-08-15-refdb-rebuild/progress.md`, is the full record. It is
**gitignored and on disk only**; `git clean -fdx` destroys it. Read it before doing anything.

## What is done

All seven tasks of `dev/design/plans/2026-08-15-refdb-rebuild.md`, both maintainer-directed post-plan
items, the final whole-branch review, its single fix wave, and the narrow amendment.

| | commits | state |
|---|---|---|
| Task 1 orchestrator, preflight, version sync | `d0269ae`, `fc33cc1` | complete, review clean |
| Task 2 accessions stage | `44ef4f0`, `2babc97` | complete, review clean |
| Task 3 taxonomy mapping stage | `5965670` | complete, review clean |
| Task 4 resumable fetch stage | `96c0141`..`8875507` | complete, review clean |
| Task 5 build stage | `aa7c33c`..`a329bd2` | complete, review clean |
| Task 6 decoys and packaging | `a329bd2`..`66249c9` | complete, review clean |
| Task 7 docker wiring | `66249c9`..`ad7c57f` | complete, review clean |
| Post-plan: unexport + docs | `ad7c57f`..`f3fd2e8` | complete, review clean |
| Final review fix wave | `f3fd2e8`..`b3d4921` | complete, re-reviewed |
| Narrow amendment | `69fd340` | re-reviewed 2026-08-18; one Critical, fixed (uncommitted) |

Out of plan: `da6426c` removed `efilter -molecule genomic` from `download_rrna.sh`.

The final whole-branch review's verdict was **ready to merge, with fixes**: 0 Critical, 6 Important,
4 Minor. All six Importants and eight folded minors were fixed in one wave. The wave's own re-review
found one new Important regression, which `69fd340` closes. The scoped re-review of `69fd340` itself
then found **one Critical**, now fixed but uncommitted; see THE IMMEDIATE NEXT ACTION.

**The first real build run, 2026-08-18, found two further defects that six review rounds had not.**
Neither was reachable by reading: `efetch -input` handed the history envelope to a flag that means
"read identifiers from this file", and EDirect's own HTTP/1.0 default made every NCBI call fail. Both
are in Traps. This is the strongest evidence on the branch for its own process note 2, "require the
implementer to hit the real thing".

## The measured findings that reorganised this branch

Everything below was measured against the real historical build tree and the shipped panel, read-only,
before any brief was written. Each corrected a claim in the plan, the spec, or my own text.

**One name-translation bug explains every phantom number in the design docs.** NCBI answers a request
for `28ZU_S2` with a record deflined `>pdb|28ZU|S2`. The old build script never translated it back, so
all 298 PDB records fell out of its reverse reconciliation as "not returned", 11 WGS masters joined
them, and that is the manifest's `unreconciled_accessions: 309`. The same 298 were then dropped a
second time as `pdb_structure_record_nonstandard_defline`. One double-count.

Dropping PDB records *after* the dedup rather than before adds the rest: 188 of them collapse against
each other under the shared `NA` taxid and 110 survive. That is simultaneously the spec's "ledger off
by 188", its "`dropped.tsv` has 110 more rows than `records_dropped`", and its "598,098 then 597,988"
silent-rebuild note. Same bug, three symptoms, three separate places in the documentation.

| | plan/spec claimed | measured |
|---|---|---|
| requested but not returned | 309 | **11** (the WGS masters) |
| returned but not requested | 298 | **0** |

The corrected decomposition balances exactly, and is the build stage's verification target:

```
891382 requested = 597988 kept
                 +     11 wgs_master
                 +    298 pdb_structure
                 + 293085 byte_identical
                 +      0 efetch_unavailable
                 +      0 taxid_unresolved
```

Reference values from the **2026-08-13 dev panel** at `~/ribopilot_refdb_build/blastdb/2026.08.1`
(built by the old scripts, never published): fasta records 597,988; taxid map rows 597,988; distinct
taxids 185,732; `.njs` `number-of-sequences` 597,988, `number-of-letters` 688,449,026,
`number-of-taxids` 185,732.

**The new pipeline's first successful build, 2026-08-18, does not match those numbers**, and should
not be expected to: the source grows ~370 records/day and three days passed. It produced 598,202
records (+214), 688,824,954 letters (+375,928) and 185,708 taxids (-24). It was first built under the
same version string `2026.08.1`, which meant two different panels carried one label. **Resolved by
maintainer ruling: the panel is now versioned `YYYY.MM.DD`, the harvest date** (`a9c8ab4`). The
current panel is **`2026.08.18`**, built and verified; the 2026.08.1 tree under
`dev/refdb/build/blastdb/` was removed as a mislabelled duplicate of it. Nothing was ever published,
so no release tag or image collides.

The scheme change is deliberate, not cosmetic: a hand-bumped counter is a manual step that gets
forgotten, which is exactly how the collision happened. The date is unique per build by construction
and records how old the snapshot is. It stays dotted-numeric because `rp_valid_version` and R's
`numeric_version()` both reject anything else, and `discover_panel()` would then be unable to see the
panel at all. Only two files carry the value (`R/db_dist.R`, `docker/Dockerfile`) and
`build_refdb.sh --check-version` diffs them.

**Two live defects in existing tooling, both confirmed by running the code:**

- `package_release.sh` collected taxonomy files with a `taxdb.*` glob. The main lookup file is
  `taxonomy4blast.sqlite3`, which does not match, so **it had never been included in a release**, while
  the guard meant to catch a missing taxonomy file was satisfied by the two files that did match. Per
  spec section 4, a database without it returns unfiltered `-taxids` rows and exits 0.
- **`run_annotate()` does not exist.** Deleted per `dev/design/15_WF2_DECOMMISSION.md`, yet the
  getting-started vignette still instructed users to call it, so Step 7 hard-errored for every reader.
  The surrounding prose was also backwards: WF1 now annotates every scaffold as it is assembled, and
  the Annotate tab is a review surface. Both fixed in the wave.

## Every ruling I made, and what each costs if wrong

Ordered as made. Each is recorded in the ledger with its evidence.

### Task 4, the fetch stage

- **R0** The plan's "309 want-not-in-got and 298 got-not-in-want" is wrong; measured 11 and 0. Its
  conclusion (no count-equality assert) stands, on the different basis that GenBank can withdraw a
  record inside the fetch window. *Cost: a genuine mass-withdrawal event aborts and needs the override.*
- **R1** `--jobs` defaults to 4, not 1. *Cost: more concurrent NCBI load by default.*
- **R2** The `taxidmap`-concurrent-with-`fetch` overlap of spec section 2 is deferred, not implemented.
  It saves ~4 min of a ~44 min run and costs both the scrub-flush child-reaping hazard and contention
  for the same 10 req/s E-utilities budget. *Cost: the rebuild takes about 4 minutes longer.*
- **R3** One tolerance policy, tight: 3 attempts per chunk, remaining misses recorded as
  `efetch_unavailable`, a zero-record chunk is a hard error, total permanently-missing ceiling 100
  (`RP_MAX_UNAVAILABLE`), set against a measured expectation of 0. *Cost: a mass-withdrawal event needs
  the override.*
- **R4** got-not-in-want is zero-tolerance after `pdb|X|C -> X_C` normalisation. *Cost: a legitimate
  non-PDB source would abort a build.*
- **R5** Artifact names pinned. *Cost: renaming later touches two stages.*
- **R6** A chunk is complete only when both its chunk file and its missing-list exist. *Cost: one chunk
  re-verified after a kill between the two renames.*
- **R7** `cmd_clean` sweeps the new artifacts. *Cost: none.*
- **R8** Every abort path reaps its background children first. *Cost: three lines of trap.*
- **R9** An empty or missing accession list is a hard error. *Cost: none.*
- **R10** The resume key is recomputed from `accessions.txt`, not trusted from the sidecar.
  *Cost: none.*
- The implementer's out-of-brief guard is **kept** (per-chunk missing lists must agree with the chunks
  on disk). *Cost: 2 lines and one test.*
- **`FETCH_ROUNDS=3` means 3 attempts total**, not 1+3. *Cost: a transient outage surviving three
  attempts trips the ceiling and needs a re-run, which resumes from chunks.*
- Reviewer Minor 8 **upgraded to Important** and fixed: Ctrl-C left workers running, and an orphan can
  publish a chunk into a build tree a later run is already using. *Cost: three lines of trap.*

### Task 5, the build stage

- **S0** The historical ledger's brokenness is one PDB double-count; corrected decomposition above.
- **S1** PDB records **are** dropped at the build stage, contradicting one sentence of spec section 3.
  The reconciled reading: section 3 governs the QUERY (never add `NOT pdb[Filter]`), and the records
  are still harvested, fetched, reconciled and recorded as an explicit drop with a reason. Two verified
  reasons they must go: `extract_species()` keeps the first two words of the title, so
  `>pdb|28ZU|S2 Chain S2, 18S ribosomal RNA` displays a species of "Chain S2," in the app; and
  `-parse_seqids` yields a seqid that does not match the taxid map's key, so they would carry no
  staxid. **Audited by the final reviewer and upheld for a stronger reason than I gave**: keeping them
  would make `.njs number-of-taxids` disagree with the restricted map and abort the build outright, and
  since a content change is a new panel version, reproducing 2026.08.1 *requires* the drop.
  *Cost: 298 of 597,988 records absent, named in the shipped ledger, reversible by deleting one filter.*
- **S2** No header row in the shipped `dropped.tsv`. *Cost: none; removes an off-by-one.*
- **S3** `rp_filter_wgs_masters` gains `dropped.accessions.tsv`. *Cost: two lines.*
- **S4** Never land a raw concatenated FASTA. *Cost: chunks are read twice.*
- **S5** The dedup tie-break is explicit `-k1,1 -k2,2 -k3,3`, not left to GNU sort's last-resort
  comparison. *Cost: none; protects 95,997 records of panel content.*
- **S6** An unmapped accession is an abort, not `"NA"`. **This ruling rested on a false premise and the
  implementer corrected it**: `stage_taxidmap` records unresolved accessions but does not remove them
  from `accessions.txt`, so a bare abort would kill a build after the fetch. The shipped behaviour
  skips accessions already in `dropped.taxidmap.tsv` and aborts on any other. *Cost: none; the
  correction was verified.*
- **S7** Chunk count must equal `ceil(final_count / CHUNK_SIZE)`. *Cost: none.*
- **S8** `unpacked_bytes` covers only what `package_release.sh` ships. *Cost: none.*

### Task 6, decoys and packaging

- **T1** The taxonomy file list is explicit, not a glob. *Cost: none; this is the fix for a live bug.*
- **T2** The shell's four-name decoy list is **not** merged into R's `RP_DECOY_DBS`, because that
  constant is composed into a blastn `-db` list and taxdb is a sidecar, not a searchable database. Drift
  is caught by a test instead. *Cost: a test rather than a shared constant.*
- **T3** Version capture is implemented and tested without downloading any database. *Cost: none.*
- **T4** `~/.cache/ribopilot-blastdb` is a read-only dev fixture, never a source. *Cost: none.*
- **T5** Decoy verification requires the four `PANEL_REQUIRED_EXT` extensions, `blastdbcmd -info`, and
  all three taxonomy files. *Cost: none.*
- **T6** Never run `--publish`. *Cost: none.*
- **T7** Task 6 owns the `test-db-release.R` rewrite; Task 7 appends only. *Cost: none.*
- **T8** The on-disk `.tar.gz.md5` is a cache key, never evidence. *Cost: none.*
- **T9** taxdb is fetched last. *Cost: none.*

### Task 7, docker wiring

- **D1** I **rejected** plan Step 1 and spec section 5's instruction to delete the staged database
  tarballs before every build. Three pieces of evidence: the deploy scripts document them as a
  deliberate 82 MB + 77 MB cache and one says outright "Not covered by the rm -f above"; the stale-glob
  hazard does not apply because the Dockerfile ADDs them by exact versioned name while the only glob is
  the package tarball, already handled; and every asset, cached or fresh, is already checksum-verified
  and deleted on mismatch, which is stronger than deleting up front. Replaced with a version-scoped
  sweep. **Independently verified twice, and audited by the final reviewer.** *Cost: a stale
  same-version tarball stays cached, which is the documented intent.*
- **D2** The verification RUN goes in `db-1` after both ADDs, with no pipe. *Cost: none.*
- **D3** Task 1 owns the ARG-equality test; Task 7 does not re-add it. *Cost: none.*
- **D4** Task 7 appends to `test-db-release.R`, never rewrites. *Cost: none.*
- **D5** Verification splits into a cheap host-side behavioural test and one real `--target db-1`
  build. *Cost: none.*
- **D6** Never push, publish, or execute `deploy-dockerhub.sh`. *Cost: none.*

### Post-plan items

- **P1** `NAMESPACE` and `man/` are hand-maintained; `devtools::document()` must not be run.
  *Cost: none.*
- **P2** `README.md` is generated from `README.Rmd`. *Cost: none.*
- **P3** Narrow scope: state the support policy, do **not** delete the conda instructions, do **not**
  change any `executor` default. *Cost: the docs say conda is unsupported while the code still defaults
  to it, until the maintainer decides. See the open question below.*
- **P4** Do not invent instructions. **This ruling's factual premise was stale and the implementer
  overrode it**: the image *is* published as `macguigand/ribopilot`, nine tags, one matching the current
  `DESCRIPTION` version. *Cost: none; corrected before it shipped.*

### Final review

- Finding 6 (zero CI coverage of the refdb and docker tooling) is **recorded as a decision** in spec
  section 7, and the CI job is **not** added. A new workflow cannot be verified without pushing, and
  shipping an unverified workflow would contradict the evidence-before-assertion standard this branch
  held to throughout. *Cost: the gap stays open until the maintainer acts on the recorded
  recommendation.*
- **A narrow amendment was dispatched, deviating from the process's "exactly one fix wave".** The rule
  exists to stop scope spiral, not to force shipping a regression the wave itself introduced, and the
  fix was small. *Cost: one more review cycle, which is the pending next action above.*

## Open question for you, deliberately not resolved

The docs now state that Docker/Singularity is the supported execution path and conda is not supported,
while `new_project()` and `new_test_project()` still default to `executor = "conda"`
(`R/init_project.R:43`, `R/init_test_project.R:104`). The mismatch is stated plainly in four places in
the docs rather than hidden.

Either the defaults move to `"docker"` (a behaviour change affecting every existing workflow, outside
what this branch set out to do) or the docs soften. I did not guess, because this branch already
reversed one large deletion made from a similarly terse directive: `dev/design/19` was going to remove
~4,900 lines of clade-database building, which would have silently assembled every non-fish sample
against hardcoded Actinopterygii seeds, exit 0, with half the path score computed against fish rDNA.

## Residuals, in priority order

Recorded, not fixed, under the maintainer's 2026-08-18 stopping rule (only a Critical gets fixed).
None blocks merge on its own. All are in the ledger with evidence.

1. **`EXPECTED_COUNT` should come out of the panel key** (re-review Q1, OVERRULED). It selects no
   record and never reaches the manifest, so it buys no wrong-panel protection, while making a
   re-centre of that constant invalidate the accessions marker and re-harvest -- overwriting
   `accessions.txt` in place with no backup, which strands a completed multi-hour fetch. Reachable
   with no `--force` and no query edit, by the route the script itself prescribes (`--accept-drift`,
   record the count afterwards). Fix is one line: drop `"$EXPECTED_COUNT"` from `rp_panel_key`'s
   printf. Costs the `EXPECTED_COUNT` case in one test; the query mutant still fails 6.
2. **A `--force` restart can leave "fetch complete" and "nothing dropped" true at once.**
   `: > dropped.fetch.tsv` is not paired with clearing the fetch marker, and the marker gate carries
   no `$STAGE` term, so the skip branch is dead exactly on the path `--force` exists for. Bounded:
   `chunks.missing/` survives and the build stage aborts loudly, except where the same accessions are
   also in `dropped.taxidmap.tsv`. One line.
3. **The getting-started walkthrough dead-ends.** It calls the app "three tabs" against five modes,
   names an "Annotate tab" that is labelled "Annotation Review", and **never mentions the Taxonomy
   step at all** -- while the Assemble tab queries `WHERE taxonomy_lock = 1` and `init_db.R` sets that
   to 0 on every fresh project. A reader who follows the guide spends the multi-hour assembly first,
   then finds an empty tab. Mitigated only because the app opens on Taxonomy and the empty state
   prints the remedy. Pre-existing on `main`, not introduced here.
4. **`cmd_clean` deletes `dropped.accessions.tsv` but leaves the accessions marker green**, so a
   same-version rebuild after `clean` skips the harvest, re-fetches everything, and only then dies
   asking for the deleted file. One word: add `accessions` to the marker-clear loop.
5. **`CHUNK_SIZE` is in neither the panel key nor the resume key.** Changing it with a complete fetch
   present skips the fetch and aborts in the build stage with a message naming neither. Always loud,
   never a wrong panel. Diagnostics only.
6. **Four write-only markers** (`taxidmap`, `build`, `decoys`, `package`) are written and never read.
   Proven inert. A one-line comment saying so would stop the next reader trusting them as gates.
7. The **sourced test drivers disable errexit** (`set +Ee`, `trap - ERR`), so stages are exercised in a
   shell that is not the shell they ship in and the ERR trap is never exercised. Not cheaply fixable.
8. **Zero CI coverage** of the refdb and docker tooling: `.Rbuildignore` carries `^dev$` and
   `^docker$`, and CI leaves `NOT_CRAN` unset, so the branch's largest test file is double-gated out.
   Recorded as a decision in spec section 7; the job was deliberately not added because it cannot be
   verified without pushing.
9. `RiboPilot_test_hydra.R:70` carries a commented, inert `run_annotate(proj)`.
10. **25 minor items** are in the ledger as `minor (deferred)`, each with reasoning.

## Security

On 2026-08-17 a diagnostic printed the maintainer's `NCBI_API_KEY` to a terminal. On 2026-08-18 it
leaked again on its own: EDirect prints the full request URL, api_key included, on every failed call,
and the HTTP/1.0 defect above made every call fail. It is also stored in plaintext at
`RiboPilot_test/run_full_project.R:51`. **Rotate the key.** No file in the repo contains it.

## Maintainer decisions that constrain the work

- **End users must not build reference databases.** Databases ship in the image; the maintainer
  refreshes them alongside each software release.
- **BUT clade building is KEPT.** The maintainer re-answered this after learning what clade databases
  actually are. Read `dev/design/19_REFDB_USER_BUILD_REMOVAL.md`'s new header before touching that file.
- `install_ribopilot_db()` is **kept but unexported**. The roxygen header was removed from `NAMESPACE`,
  which makes roxygen *refuse* to rewrite the file rather than merely discouraging it, so a future
  `document()` run cannot silently re-export it.
- Docker/Singularity is the supported execution path; conda and individual dependency installation are
  explicitly unsupported.
- Pre-release with **no end users**, so no schema migrations are ever written.

## Definitional constants: changing these changes the panel

Live in `dev/refdb/build_refdb.sh`. A change here is a new panel version, never a silent fix. **As of
`69fd340` the stage markers are keyed on `QUERY` and `EXPECTED_COUNT`**, so changing either now forces
a re-harvest instead of silently reusing the previous panel's records.

```
TAXID=33208
QUERY="(txid33208[Organism:exp]) AND (18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR
5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal
RNA[Title] OR internal transcribed spacer[Title]) NOT mitochondrion[Filter] NOT plastid[Filter]"
```

- **Never add `efilter -molecule genomic`.** Measured twice: 892,123 records without it, 741,526 with.
  It removes 100% of RefSeq, whose MolType is `rna`. A test asserts the string is absent.
- Chunk size 5000. Dedup key `(taxid, uppercased sequence)`, byte-exact, survivor tie-break
  lexicographic under `LC_ALL=C`.
- `EXPECTED_COUNT=892123`, tolerance 15000 (~5 weeks of the measured ~370 records/day drift).
- `EXPECTED_COUNT_FLOOR=446000`, overridable in tests via `RP_MIN_ACCESIONS`. **`--accept-drift` cannot
  waive the floor.**
- `MAX_UNAVAILABLE=100`, overridable via `RP_MAX_UNAVAILABLE`. Measured expectation is 0.

## Contracts a later stage must honour

- **Fingerprint.** `accessions.txt` is sorted and `accessions.sha256` is over the sorted list. Do not
  un-sort it.
- **Three counts, deliberately distinct** (was two until 2026-08-18; see below).
  `accessions.esearch_count` is what the source CLAIMED matched, and drives only the collapse and
  drift guards. `accessions.requested_count` is what the history actually returned, and is the ledger
  identity's `requested` and the manifest's `requested_accessions`. `accessions.final_count` is
  post-WGS-filter and is what the chunking divides. There is no file called `accessions.count`.
  **`<Count>` and the number of accessions the history returns are NOT equal** and E-utilities does
  not promise they will be: measured 892,364 reported against 891,825 and 891,876 returned, a ~0.06%
  shortfall uncorrelated with network failures. The historical build never compared them either: it
  took its total from `wc -l` of the list it got (`01_harvest.sh:53`). Divergence is bounded by
  `MAX_UNRETURNED` (`RP_MAX_UNRETURNED`, default 5000) and logged whenever non-zero.
- **Dropped records.** Each stage owns its own file (`dropped.accessions.tsv`, `dropped.taxidmap.tsv`,
  `dropped.fetch.tsv`) and rewrites it wholesale. The build stage merges them into
  `ribopilot_rrna.dropped.tsv` by an explicit precedence: `wgs_master_unavailable`,
  `efetch_unavailable`, `taxid_unresolved`, `pdb_structure_record_nonstandard_defline`,
  `byte_identical_to_kept_record_same_taxid` - how far the record got, earliest failure winning. The
  ledger is checked as a **set**, not a count, so a double-count cannot cancel a miss.
- **Fetch outputs.** `chunks/fa_NNNN` (4-digit, zero-padded, 0-based; nothing else in that directory
  ever), `chunks.missing/fa_NNNN.txt`, `chunks.accessions.sha256`, `chunks.tmp/`. A chunk is complete
  only when both its files exist.
- **Normalisation.** `rp_fasta_accessions` is the one place the `pdb|X|C -> X_C` rule is spelled. Its
  absence is the entire historical ledger defect. Reuse it.
- **Sourcing.** `RP_EXECUTED` gates the scrubber, both traps, the `BUILD_DIR` default and the argument
  loop. A sourced stage calls `rp_require_build_dir` first.
- **Reap children.** Anything the EXIT, ERR or INT/TERM traps can reach must be defined above
  `rp_cleanup`.

## Traps that have already cost time on this branch

- **`efetch` exits 0 with 0 bytes** and `QUERY FAILURE` on total failure. Never branch on its exit
  status.
- **`efetch -input FILE` reads IDENTIFIERS from that file.** It is NOT how you replay an esearch
  history: the `<ENTREZ_DIRECT>` envelope has to arrive on STDIN (`efetch -format acc < xml`, or the
  `esearch | efetch` pipe the historical build used). Hand it the envelope and it has no database and
  dies "Missing -db argument". Two other `-input` uses in the script are correct because they pass
  `-db` and a real list of accessions, which is how the wrong one read as idiomatic for weeks.
- **EDirect 25.5 forces `curl --http1.0`** (`nquire:483`). Against curl 7.81 / OpenSSL 3.0.2 that
  fails EVERY NCBI request (measured 12/12 fail on 1.0, 12/12 succeed on 1.1, API key present so rate
  limiting is excluded). EDirect's internal retries turn it into a slow partial success rather than a
  clean failure. Patched locally to `1.1`; backup at `~/edirect/nquire.bak`. Affects every EDirect
  caller on the host, not just this build.
- **A stub that ignores its arguments is not a stub, it is a hole.** The `efetch` stub for the
  accessions stage `cat`-ed a fixture regardless of how it was called, so the `-input` bug above
  survived implementation, task review, a whole-branch review, a fix wave, and two re-reviews. The
  `xtract` stub two lines above it validates `-input` AND `-pattern` and deliberately fails on a wrong
  one, which is why the equivalent bug in xtract was caught immediately.
- **EDirect wraps xtract output in `<ENTREZ_DIRECT>`**, not `<eSearchResult>`.
- **`$(grep -c . f || echo 0)`** on an existing empty file yields `"0\n0"`; `[` then errors rc 2 and the
  guard is silently skipped. Use `wc -l`.
- **`update_blastdb.pl`'s `.md5` is a cache key, not evidence.** With the marker present and the volumes
  deleted it prints "up to date", installs nothing, exits 0.
- **`makeblastdb` exit status proves nothing.** Parse the generated `.njs` instead.
- **Docker `ADD` does not fail on a corrupt gzip-shaped archive**; it silently plain-copies it. This is
  now demonstrated by a real build, not quoted.
- **Bare errexit** killed the original build whenever zero PDB records matched, a legitimate state.
  `|| true` on every grep used as a filter.
- **Never `pgrep`/`ps`/`/proc` for liveness.**
- **A redirect around a compound command is a deadlock trap under the scrubber.** Bash parks the real
  stderr on fd 10 for `{ ...; } 2>/dev/null`; a trap that exits from inside it leaves
  `rp_scrub_flush`'s `exec 2>&-` closing the wrong descriptor. Use `set +m`, not a redirect.
- **A background subshell inherits the ERR trap under `set -E` but does NOT run the EXIT trap.**
- **`system2(wait = FALSE)` is a silently inert harness for interrupt tests**: a non-interactive shell
  sets SIGINT to `SIG_IGN` for `&` children and bash then refuses the trap.
- **awk compares strnum fields numerically**, so two numeric-looking accessions differing only in
  trailing zeros collapse. Force string comparison.
- **Do not run a full rebuild to test.** 2 h 41 m, ~892k NCBI records.

## Process notes that have paid off

Work ran through `superpowers:subagent-driven-development`: one implementer per task, a review after
each, a scoped re-review after each fix round, and a ledger that survives compaction.

Three habits kept catching real defects and are worth keeping:

1. **The controller measures the claims itself before writing the brief.** Every phantom number on this
   branch was found this way, by reconciling the real historical build tree read-only rather than
   trusting the plan.
2. **Require the implementer to hit the real thing.** A real esearch caught the `<ENTREZ_DIRECT>` shape;
   a real Salvelinus query confirmed the RefSeq loss at 431-vs-233 records; a three-record efetch
   confirmed the `>pdb|28ZU|S2` defline; four real Docker builds proved the corrupt-archive behaviour.
3. **Demand measured mutation results and a second pass.** Task 4's first mutation pass scored **zero
   failures on four of twenty guards**. The final amendment's own author found a surviving mutant in the
   test he had just written, because changing two constants against one build tree let each case pass
   on the other's mismatch.

**Fifteen assertions that could not fail have been found on this branch, and eight defects originated
in plan, spec, or controller text rather than implementer error.** One controller-authored "cheap
minor" was itself a deadlock, caught the same round by a test the controller had just forced the
implementer to strengthen. Two controller instructions contradicted each other. One controller ruling
rested on a false premise about what an earlier stage did; another on a stale fact about a published
image. **Treat the plan, the spec, and the brief as the likeliest sources of bugs**, and require
implementers to check rather than comply - most of these were caught that way.

Documentation accuracy was treated as a correctness concern, not a style one. That is what surfaced
Step 7 of the getting-started guide instructing a function deleted months ago, with surrounding prose
describing the opposite of what the code does. It was invisible to a line-by-line read and only showed
up when someone read the page in order, as a user.
