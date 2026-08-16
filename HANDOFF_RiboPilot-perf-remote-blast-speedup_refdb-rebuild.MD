# Handoff: RiboPilot, branch `perf/remote-blast-speedup`

Written 2026-08-15. Resume point is exact: everything below is committed, the tree is clean,
and nothing has been pushed.

## State

- Branch `perf/remote-blast-speedup`, **86 commits** ahead of `origin/main`, linear, **no merges**.
- HEAD `2babc97` "fix: refuse an empty harvest and stabilise the accession fingerprint".
- Working tree clean. **Nothing pushed. Do not push without asking.**
- Suite: **2909 passed / 0 failed / 20 skipped** by default; **2984 passed / 0 skipped** with
  `NOT_CRAN=true`. Nine warnings are pre-existing and expected: they come from tests that
  deliberately exercise warning paths (missing decoy, whitespace decoy dir, corrupt ref store).
  Do not chase them.

Run the suite exactly this way, in the foreground, never backgrounded:

```
Rscript -e 'suppressMessages(pkgload::load_all(".", quiet=TRUE)); r <- testthat::test_dir("tests/testthat", reporter="silent"); d <- as.data.frame(r); cat("failed:", sum(d$failed), "error:", sum(d$error), "skipped:", sum(d$skipped), "passed:", sum(d$passed), "\n")'
```

and the same prefixed with `NOT_CRAN=true`.

## What this branch delivered before the current work

Five stages, all landed and reviewed. Measured outcomes are recorded in
`dev/design/17_REMOTE_BLAST_PERFORMANCE.md`, which is the spec for all of it.

| stage | result |
|---|---|
| 1 | searching split from interpreting, `qcovs` fix; units with candidates 16 -> 26 retroactively, zero NCBI calls |
| 2 | local barrnap/ITSx annotation of references; usable references 26.7% -> 100% |
| 3 | cross-project reference store; one project satisfied 48 of another's 232 accessions in 0.227 s |
| 4a | local BLAST panel search; a real unit in 1.5 s against a ~2 h remote baseline |
| 4b | `install_ribopilot_db()`, checksummed download, staged extraction, per-project provenance |
| 4c | panel + decoys baked into the Docker image, container-side auto-discovery, zero config |

The panel itself is built and validated: **597,988 records, 688 Mbp, 332 MB unpacked, 76 MB gzipped**,
living at `~/ribopilot_refdb_build/blastdb/2026.08.1/`. Verified end to end in a real image: a
`RP_DB=1` build discovers the panel and returns real hits with no configuration; a `--no-db` build
falls back to remote. 6.33 GB with the database, 5.98 GB without.

## The work in progress: reproducible reference-database rebuild

**Why.** The recipe that built the panel existed only as 508 ad-hoc lines in
`~/ribopilot_refdb_build/*.sh` on one machine. A clean checkout could not rebuild it. That work is
being replaced by committed tooling under `dev/refdb/`, following MitoPilot's split: the recipe is
committed, the multi-hundred-MB artifact is not.

**Authority documents, read in this order:**

1. `dev/design/18_REFDB_REBUILD.md` - the SPEC. Binding. Where the plan and the spec disagree, the spec wins.
2. `dev/design/plans/2026-08-15-refdb-rebuild.md` - the implementation plan, 7 tasks.
3. `.superpowers/sdd/2026-08-15-refdb-rebuild/progress.md` - the ledger. Gitignored, on disk only.
   Every ruling, every review finding, every measured mutation result. **Read it before dispatching
   anything.** Key content is summarised below in case it is lost.

### Done

- **Task 1** (`d0269ae`) orchestrator skeleton, preflight, version sync.
- **Task 1 fixes** (`fc33cc1`) atomic lock and working secret scrubber. Both were Critical.
- **Task 2** (`44ef4f0`) accessions stage.
- **Task 3** (`5965670`) taxonomy mapping stage.
- **Task 2 fixes** (`2babc97`) collapse floor, two named counts, sorted fingerprint, sourcing guard.
- Out of plan: **`da6426c`** removed `efilter -molecule genomic` from `download_rrna.sh`.

### Next, in order

1. **Task 4, the fetch stage.** `dev/refdb/build_refdb.sh:701` is still `stage_fetch() { rp_die ... }`.
   This is the stage where every historical failure lived and it is 2 h 20 m of the 2 h 41 m rebuild.
   Brief: `.superpowers/sdd/2026-08-15-refdb-rebuild/task-4-brief.md` (generate with
   `scripts/task-brief <plan> 4` if absent). Spec section 4, "Fetch".
2. **Task 5**, build stage (`:702`). Spec section 4, "Build / package".
3. **Task 6**, decoys rewrite plus `package_release.sh` patch (`:703`, `:704`).
4. **Task 7**, docker wiring: `rm -f docker/ribopilot_rrna-*.tar.gz docker/decoys-*.tar.gz` in both
   deploy scripts, plus a build-time `blastdbcmd -info` verification in the Dockerfile's `db-1` stage.
5. **Unexport `install_ribopilot_db()`** - `@noRd`, drop from `NAMESPACE` and `man/`. Keep the code
   and its tests. This is the maintainer's explicit decision, see below.
6. **Docs**: state that Docker/Singularity is required and conda is not supported.

## Maintainer decisions that constrain the work

- **End users must not build reference databases.** Databases ship in the image; the maintainer
  refreshes them alongside each software release.
- **BUT clade building is KEPT.** The maintainer re-answered this after learning what clade
  databases actually are. See the next section; this reversal matters.
- `install_ribopilot_db()` is **kept but unexported**, "in case we want it later".
- Docker/Singularity is the supported execution path. Conda and individual dependency installation
  are explicitly unsupported.
- Pipeline is pre-release with **no end users**, so existing artifacts do not need migrating and
  no schema migrations are ever written.

## The thing most likely to be misunderstood

`dev/design/19_REFDB_USER_BUILD_REMOVAL.md` contains an accurate architecture map **and a deletion
list that MUST NOT be executed.** The maintainer reversed that decision. Read this before acting on
that file.

**A clade "database" is not a BLAST database.** It is a pair of FASTA files, GetOrganelle seeds and
labels. They feed ASSEMBLY: GetOrganelle `-s`/`--genes` (`assemble.nf:24`), the bowtie2 bait index
and rDNA contig filter (`assemble_fallback.R:177,188`), and 50% of path-scoring weight
(`scaffold_scoring.R:12`). The shipped panel feeds IDENTIFICATION only (`blast.nf:50,60,109`).
Neither can substitute for the other. GetOrganelle runs in `-F anonym` mode, which requires
user-supplied seeds; there is no built-in-database path in this codebase.

Deleting clade building would have silently assembled every non-fish sample against the hardcoded
Actinopterygii default, exit 0, with half the path score computed against fish rDNA.

## Definitional constants: changing these changes the panel

Live in `dev/refdb/build_refdb.sh`. A change here is a new panel version, never a silent fix.

```
TAXID=33208
QUERY="(txid33208[Organism:exp]) AND (18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR
5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal
RNA[Title] OR internal transcribed spacer[Title]) NOT mitochondrion[Filter] NOT plastid[Filter]"
```

- **Never add `efilter -molecule genomic`.** Measured twice: 892,123 records without it, 741,526
  with. It removes 100% of RefSeq, whose MolType is `rna`. A test asserts the string is absent.
- Chunk size 5000. Dedup key `(taxid, uppercased sequence)`, byte-exact, survivor tie-break
  lexicographic under `LC_ALL=C`.
- `EXPECTED_COUNT=892123`, tolerance 15000 (about 5 weeks of the measured ~370 records/day drift).
- `EXPECTED_COUNT_FLOOR=446000`, overridable in tests via `RP_MIN_ACCESSIONS`. **`--accept-drift`
  cannot waive the floor.** That is the whole point of it.

## Contracts a later stage must honour

- **Fingerprint.** `accessions.txt` is sorted, and `accessions.sha256` is over the sorted list.
  Task 4's chunk cache keys on that hash and refuses to resume on mismatch. Do not un-sort it.
- **Two counts, deliberately distinct.** `accessions.esearch_count` is the pre-filter esearch
  `<Count>`; `accessions.final_count` is post-WGS-filter. They differ by the WGS masters (11 in
  production). Chunk on `final_count`. There is no file called `accessions.count`.
- **Dropped records.** Each stage owns its own file and rewrites it wholesale
  (`dropped.taxidmap.tsv`). The build stage merges them into the shipped `dropped.tsv`. A single
  shared file would let one stage's idempotent rewrite discard another's rows.
- **Sourcing.** `build_refdb.sh` can be sourced safely: `RP_EXECUTED` gates the scrubber, both
  traps, the `BUILD_DIR` default and the whole argument loop. A sourced stage calls
  `rp_require_build_dir` and refuses rather than writing into the developer's real tree. If you add
  a stage, call that guard first.
- **Reap children.** The scrubber waits on filter PIDs at exit; an orphan holding stderr open
  blocks that wait.

## Traps that have already cost time on this branch

- **`efetch` exits 0 with 0 bytes** and `QUERY FAILURE` on total failure. Never branch on its exit
  status. This made the original repair script's entire retry loop inert.
- **EDirect wraps xtract output in `<ENTREZ_DIRECT>`**, not `<eSearchResult>`.
- **`$(grep -c . f || echo 0)`** on an existing empty file yields `"0\n0"`; `[` then errors rc 2 and
  the guard is silently skipped. Use `wc -l`, and test every numeric guard against empty input.
- **`update_blastdb.pl`'s `.md5` is a cache key, not evidence.** With the marker present and the
  volumes deleted it prints "up to date", installs nothing, exits 0. Verify after fetch.
- **`makeblastdb` exit status proves nothing.** A truncated FASTA gives "added 2 sequences" exit 0;
  a taxid map matching nothing prints `Error:` and still exits 0 with every staxid 0. Parse the
  generated `.njs` instead.
- **Docker `ADD` does not fail on a corrupt gzip-shaped archive**; it silently plain-copies it.
- **Bare errexit** killed the original build whenever zero PDB records matched, a legitimate state.
  `|| true` on every grep used as a filter.
- **Never `pgrep`/`ps`/`/proc` for liveness.** The original chain waiter deadlocked because
  `pgrep -f` matched its own argv.
- **Do not run a full rebuild to test.** 2 h 41 m, ~892k NCBI records. Use stubs and fixtures. One
  tiny count-only esearch to sanity-check a query is fine and has caught a real bug.

## Known open items, all deliberate

- Only 2 of 6 version locations are verified by `--check-version`; the other 4 do not exist until
  later tasks build them.
- Fetch concurrency is not yet wired into `main_build`.
- `clean` refuses rather than cleans if `package_release.sh`'s output layout changes.
- A process killed between lock-reclaim steps can leave an unswept `build.pid.stale.<pid>`.
- CI never runs the gated tests: `r-lib/actions/check-r-package@v2` does not set `NOT_CRAN=true`, so
  20 tests including several regression guards never execute there. Pre-existing, not this branch's.
- `blastdb_version` can overstate: a configured-but-unmounted panel records it host-side at launch,
  the run goes fully remote, and nothing corrects it.
- The `db-<version>` GitHub release is **not published**, so a database-bearing image cannot be
  built until it is. `--no-db` builds an honest image without one. Publishing stays manual.

## Process notes that have paid off

Work has run through `superpowers:subagent-driven-development`: one implementer per task, a review
after each, a scoped re-review after each fix round, and a ledger at
`.superpowers/sdd/<plan-basename>/progress.md` that survives compaction.

Two habits are worth keeping because they have repeatedly caught real defects:

1. **Require the implementer to hit the real thing** rather than trust the brief. A real esearch
   caught the `<ENTREZ_DIRECT>` error; a real Salvelinus query confirmed the RefSeq loss at
   431-vs-233 records; a real container proved discovery; real Nextflow caught `--blastdb ''`
   becoming the string `"true"`.
2. **Re-measure mutation claims independently.** One implementer reported a mutation result that did
   not reproduce; the reviewer's re-run showed the tests survived the mutation for an unrelated
   reason. Ask for measured numbers, never predicted ones.

**Forty-plus defects have been caught on this branch, and the large majority originated in plan text
rather than implementer error.** Eight were assertions that could not fail, at least one written
into a brief by the controller. Treat the plan as the likeliest source of bugs.
