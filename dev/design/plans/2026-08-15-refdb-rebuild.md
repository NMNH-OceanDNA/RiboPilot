# Reference Database Rebuild Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the recipe for the rDNA panel and the contaminant databases in the repo, so a clean checkout can rebuild both from their original sources and a rebuild that goes wrong says so instead of shipping a quietly wrong database.

**Architecture:** One orchestrator, `dev/refdb/build_refdb.sh`, runs six stages in process (`accessions, taxidmap, fetch, build, decoys, package`) with a seventh, `publish`, that is never in the default run. It follows MitoPilot's split: the recipe is committed, the artifact is not. The 508 lines of ad-hoc scripts in `~/ribopilot_refdb_build/` are replaced, not ported verbatim; several of their guards are actively broken.

**Tech Stack:** bash (`set -Eeuo pipefail`), NCBI EDirect, `update_blastdb.pl`, BLAST+ 2.17, seqkit, jq, R (>= 4.4.0) for the version-sync test.

**Spec:** `dev/design/18_REFDB_REBUILD.md`. That document is the authority. Where this plan and the spec disagree, the spec wins.

## Global Constraints

- **ASCII only.** No non-ASCII in any script, test, or R file.
- **Minimal comments.** Comment why, not what. A comment explaining a load-bearing guard is wanted; a comment narrating a bugfix is not.
- **No new R package dependencies.** R floor 4.4.0.
- **Never push, never publish, never run `--publish`, never `gh release create`.** Building locally is fine.
- **Never run a full panel rebuild** during implementation. It takes 2 h 41 m and pulls ~600k records from NCBI. Test stages against small fixtures and against the existing build tree at `~/ribopilot_refdb_build/` read-only.
- **`set -Eeuo pipefail` plus an ERR trap in every script, with `|| true` on every grep-used-as-a-filter.** Bare errexit killed the original build whenever zero PDB records matched, which is a legitimate state.
- **`export LC_ALL=C` at the top of every script.** Underscored accessions sort differently under `en_US.UTF-8` and `comm` then silently lies.
- **Never scrub-fail on secrets:** `api_key=[A-Za-z0-9]+` must be scrubbed from anything written to a log. The old `harvest.log` had 35 plaintext hits.
- **`RP_DB_VERSION` in `R/db_dist.R` is the single source of truth** for the version string. Nothing else may define it.

## Definitional Constants, Copy Verbatim

Changing any of these changes what the panel contains. They are a new panel version, never a silent fix.

```
TAXID=33208
RRNA_TERMS='18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal RNA[Title] OR internal transcribed spacer[Title]'
QUERY="(txid${TAXID}[Organism:exp]) AND (${RRNA_TERMS}) NOT mitochondrion[Filter] NOT plastid[Filter]"
```

- **Never `efilter -molecule genomic`.** Measured live: query alone 892,123 records, through the filter 741,526. It removes 100% of RefSeq (148,365 records, MolType=rna).
- Chunk size 5000. Dedup key `(taxid, uppercased sequence)`, byte-exact. Survivor tie-break lexicographic by accession under `LC_ALL=C`.
- `makeblastdb -dbtype nucl -parse_seqids -blastdb_version 5 -taxid_map <map> -title "RiboPilot Metazoa nrDNA <VER>"`.
- Decoys: `16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb`, sharing one definition with `R/blast_local.R`.

---

### Task 1: Orchestrator skeleton, preflight, version sync

**Files:**
- Create: `dev/refdb/build_refdb.sh`
- Create: `dev/refdb/lib.sh`
- Modify: `.gitignore`
- Test: `tests/testthat/test-refdb-build.R`

**Interfaces:**
- Produces: `build_refdb.sh [--version V] [--build-dir D] [--stage NAME] [--jobs N] [--accept-drift] [--force] [--check-version] [status|clean]`; `lib.sh` providing `rp_log`, `rp_die`, `rp_need`, `rp_scrub`, `rp_lock`, `rp_marker_ok`, `rp_marker_set`.

- [ ] **Step 1: Write the failing test**

```r
skip_if_not(file.exists("../../dev/refdb/build_refdb.sh"))

test_that("build_refdb.sh declares the exact Entrez query the panel is defined by", {
  s <- paste(readLines("../../dev/refdb/build_refdb.sh", warn = FALSE), collapse = "\n")
  expect_match(s, "txid33208[Organism:exp]", fixed = TRUE)
  expect_match(s, "NOT mitochondrion[Filter]", fixed = TRUE)
  expect_match(s, "NOT plastid[Filter]", fixed = TRUE)
  # the filter that silently removed every RefSeq record must never reappear
  expect_false(grepl("efilter", s, fixed = TRUE))
})

test_that("every refdb script sets the strict shell options and C locale", {
  for (f in list.files("../../dev/refdb", pattern = "[.]sh$", full.names = TRUE)) {
    s <- paste(readLines(f, warn = FALSE), collapse = "\n")
    expect_match(s, "set -Eeuo pipefail", fixed = TRUE, label = basename(f))
    expect_match(s, "LC_ALL=C", fixed = TRUE, label = basename(f))
  }
})

test_that("the docker ARG default agrees with RP_DB_VERSION", {
  skip_if_not(file.exists("../../docker/Dockerfile"))
  v <- sub('.*"(.*)".*', "\\1",
    grep("^RP_DB_VERSION", readLines("../../R/db_dist.R", warn = FALSE), value = TRUE)[1])
  arg <- grep("^ARG RP_DB_VERSION", readLines("../../docker/Dockerfile", warn = FALSE), value = TRUE)[1]
  expect_equal(sub(".*=", "", arg), v)
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-refdb-build.R")'`
Expected: the file-existence skip fires, or the tests fail. Confirm which; a skip is not a pass.

- [ ] **Step 3: Implement `lib.sh`**

Shared helpers. `rp_need <tool>` aborts naming the tool and how to get it. `rp_scrub` filters `api_key=[A-Za-z0-9]*` to `api_key=REDACTED` and is applied to everything written to the log. `rp_lock` writes a PID lockfile and refuses a second concurrent run; **never use `pgrep` or match on a command line** (the original chain waiter deadlocked because `pgrep -f` matched its own argv). `rp_marker_ok <stage> <input-hash>` returns true only when the stage marker exists and records the same input hash.

- [ ] **Step 4: Implement the orchestrator skeleton**

Argument parsing, stage dispatch, and these preflight checks, all before any stage runs:

1. `rp_need` for: `esearch efetch epost esummary xtract update_blastdb.pl makeblastdb blastdbcmd seqkit jq curl sha256sum tar`.
2. **`command -v curl` must resolve inside the repo shim dir.** Prepend `inst/ref_dbs/getOrganelle/GenBankDownload/bin` to PATH, then assert. The shim's `--http1.1` is the fix for the `curl (56) unexpected eof` failures; conda curl 8.21 fails identically under `--http1.0`. The shim mitigates but does not eliminate it, so retries stay in the fetch stage.
3. **8 GB free** on the build dir filesystem, checked at startup. `sort -S` and `TMPDIR` configurable, never hardcoded `-S 4G`.
4. `NCBI_API_KEY` set, else abort. Without it the rate limit is 3/s instead of 10/s and the fetch takes hours longer.
5. Version: read `RP_DB_VERSION` from `R/db_dist.R`; if `--version` was given and differs, refuse.
6. Refuse to overwrite an existing `blastdb/<VER>` without `--force`. The directory was silently rebuilt in place once already (598,098 then 597,988 records).

Implement `--check-version`, which diffs the version across all six places listed in spec section 5 and exits non-zero on any disagreement. Implement `status` (reads a small state file; never rescans 1.2 GB, never globs `blastdb/*/manifest.json`, which returns false as soon as two version dirs exist) and `clean`.

Version format must be dotted numeric only: `numeric_version("v2026.08.1")` and `"2026.08.1-rc1"` are both `NA`, which makes the panel invisible to `discover_panel()`.

- [ ] **Step 5: gitignore the build tree**

Add `dev/refdb/build/` and `docker/ribopilot_rrna-*.tar.gz`, `docker/decoys-*.tar.gz` if not already present.

- [ ] **Step 6: Run tests, then commit**

```bash
git add dev/refdb/build_refdb.sh dev/refdb/lib.sh .gitignore tests/testthat/test-refdb-build.R
git commit -m "feat: add refdb build orchestrator with preflight and version sync"
```

---

### Task 2: Accessions stage

**Files:**
- Modify: `dev/refdb/build_refdb.sh`
- Test: `tests/testthat/test-refdb-build.R`

- [ ] **Step 1: Implement**

One `esearch` against `nuccore` with the verbatim `QUERY`. Capture the XML once, parse `<Count>`, and replay the same WebEnv/QueryKey into `efetch -format acc`. Write `accessions.txt` and `accessions.sha256`.

Guards, each of which replaces a broken one:

- **Harvested line count must equal `<Count>` exactly.** The old `-lt 500000` floor sits 44% below the true count and passes a stream that dies at 600k. Verified on three clades that `Count == lines == unique`.
- Compare `<Count>` against a committed `EXPECTED_COUNT` with a tolerance. Beyond it, refuse unless `--accept-drift`, which records the new count. The source really moves: 891,382 on 2026-08-13, 892,123 on 2026-08-15, about 370/day.
- **Never `$(grep -c . f || echo 0)`.** On an existing empty file that yields `"0\n0"`, `[` errors with rc 2, the guard is skipped, and an empty list is promoted as success. Use `wc -l`. Test every numeric guard against empty input.
- Filter `^[A-Z]{4}0{6}` (WGS master accessions) up front as a recorded expected drop. Those 11 records are permanently unavailable by design, not transient.

- [ ] **Step 2: Test the guards against empty and truncated input**

Write tests that feed the guard functions an empty file and a truncated file and assert they abort. This is the step that catches the `"0\n0"` class of bug.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add the accession harvest stage"
```

---

### Task 3: Taxonomy mapping stage

**Files:**
- Modify: `dev/refdb/build_refdb.sh`
- Test: `tests/testthat/test-refdb-build.R`

- [ ] **Step 1: Implement**

Stream the three `accession2taxid` dumps from `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/` (`nucl_gb`, `nucl_wgs`, `pdb`), never landing them: about 9.4 GB total, not the 2.5 GB an old comment claims. `curl --fail -sL | zcat | awk`, and **drop the `2>/dev/null` on zcat**: gzip's CRC catches truncation for free, and `--fail` catches the 404-HTML-into-the-pipe case. A 40% truncation silently yields 356,851 taxid-0 records and 5,457 wrongly deleted ones.

Residual accessions newer than the dump go through `epost | esummary | xtract`. Tag each row with its source (`ftp` or `esummary`) for the manifest.

Guards:

- Mapped rows must be at least 99.5% of requested accessions, else abort. Historical coverage was 100%.
- Map validator, **in this order**: every line matches `^[A-Za-z0-9._]+\t[0-9]+$`; the map's accession set is a **subset of the harvest list**; unique accessions equals row count. The subset check is the one that catches the two phantom log lines that the original build injected into the map; "row count equals unique count" passes on that corrupt map and does not catch it.
- Records with no taxid are dropped with reason `taxid_unresolved` in `dropped.tsv`. **Never a shared `"NA"` dedup key**, which collapses different species with identical rRNA and mislabels them as `byte_identical_to_kept_record_same_taxid`.
- Any log helper this stage calls must write to **stderr**. The original bug injected two log lines into the taxid map as fake accessions because a helper wrote to a stdout that was being captured.

- [ ] **Step 2: Test the validator against the historical corrupt map**

Construct a fixture map containing a plausible log line as a fake accession row and assert the subset check rejects it. Confirm the weaker "unique equals row count" check would have passed it.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add the taxonomy mapping stage"
```

---

### Task 4: Fetch stage

**Files:**
- Modify: `dev/refdb/build_refdb.sh`
- Test: `tests/testthat/test-refdb-build.R`

This is the stage that takes 2 h 20 m of the 2 h 41 m, and the one where every historical failure lived.

- [ ] **Step 1: Implement**

`efetch -db nuccore -input <chunk> -format fasta`, 5000 accessions per chunk, `--jobs N` concurrent (default 4; with an API key the limit is 10 req/s).

Guards:

- **Chunk complete means the requested and returned accession SETS match**, after normalising `pdb|XXXX|C` to `XXXX_C`, plus no zero-length sequences, plus a trailing newline. **Never `got >= want` on `>` counts**: a stream cut after the last header passes it.
- Write `.tmp`, verify, atomic rename. Leave nothing else in the chunk dir, and use the strict glob `fa_[0-9][0-9][0-9][0-9]`. The old loose `fa_*` glob eating `.tmp` files is the only reason the original run did not lose 15,000 records.
- **Retry only the missing accessions.** Never refetch or discard a whole chunk. One tolerance policy, not 100% in one script and 95% in another.
- **Never branch on `efetch` exit status.** It exits 0 with 0 bytes and `QUERY FAILURE` on total failure, which made the original repair script's entire `&& break` retry loop inert.
- Chunk cache keyed to the sha256 of the accession list, stored beside `chunks/`. On mismatch, refuse to resume and say the source moved. Positional `split` re-cutting a changed list silently keeps withdrawn records and misses new ones; pairing by name does not catch it.
- **Reverse reconcile**: got-not-in-want must be empty. This is never computed today, and the shrink case shipped 12 withdrawn records in simulation.

Explicitly do **not** add per-stage count-equality asserts: the good build had 309 want-not-in-got and 298 got-not-in-want, so equality would abort a correct build.

- [ ] **Step 2: Test with a stubbed efetch**

Put a fake `efetch` on PATH that returns: a complete chunk; a chunk truncated after the last header; and an exit-0 empty response. Assert the first is accepted and the other two are not. The exit-0-empty case is the one that made the historical retry loop inert, so it is the most important of the three.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add the resumable fetch stage"
```

---

### Task 5: Build stage

**Files:**
- Modify: `dev/refdb/build_refdb.sh`
- Test: `tests/testthat/test-refdb-build.R`

- [ ] **Step 1: Implement**

Concatenate chunks with the strict glob and assert chunk count > 0 (an empty dir currently passes `0 -ne 0`, and a stray non-empty `fa_retry` alone builds and ships a 3-sequence database with a full manifest). Dedup byte-exact on `(taxid, uppercased sequence)`. Run `makeblastdb` with the exact flags above. Write `manifest.json`.

Guards:

- **`makeblastdb` exit status proves nothing.** A truncated FASTA gives "added 2 sequences" and exit 0; a taxid map matching nothing prints `Error:` and still exits 0 with every staxid 0. Instead parse the generated `.njs` and assert: `number-of-sequences` equals the FASTA record count, `number-of-letters` equals the seqkit `sum_len`, `number-of-taxids` equals the distinct taxid count in the restricted map.
- Derive `has_taxids` from a `blastdbcmd -outfmt '%T'` sample. **Never hardcode `true`.**
- **Ledger identity** `requested == kept + sum(dropped by reason)` asserted before the manifest is written, with every count computed from the final artifacts. Today it is off by 188, and `dropped.tsv` has 110 more rows than `records_dropped`.
- `unpacked_bytes` measured over the shipped file list only. It currently reports 1.14 GB for a 331 MiB database.
- Serialise the makeblastdb command into the manifest as a JSON array. The current `${MKCMD[*]}` string is not re-runnable.

- [ ] **Step 2: Test the .njs cross-check**

Build a tiny database from a fixture FASTA, then assert the check passes. Then corrupt the count expectation and assert it fails. A test that only proves the happy path does not cover the failure this guard exists for.

- [ ] **Step 3: Commit**

```bash
git commit -m "feat: add the database build stage with artifact cross-checks"
```

---

### Task 6: Decoys and packaging

**Files:**
- Rewrite: `dev/refdb/fetch_decoys.sh`
- Modify: `dev/refdb/package_release.sh`
- Modify: `dev/refdb/build_refdb.sh` (wire both in as stages)
- Test: `tests/testthat/test-db-release.R`

- [ ] **Step 1: Rewrite `fetch_decoys.sh`**

`fetch_decoys.sh DEST`. Add `--source ncbi` to `update_blastdb.pl`, which skips about 30 s of AWS/GCP metadata probes and stops silent mirror switching. Keep **taxdb last**, with a comment saying why: every decoy tarball bundles a 2024 taxdb that would otherwise overwrite the current one.

- **The on-disk `<db>.tar.gz.md5` is `update_blastdb.pl`'s cache key, not evidence of what landed.** With the marker present and the volumes deleted it prints "up to date", installs nothing, and exits 0. Verify **after** fetch, per database: required extensions present and non-empty, plus `blastdbcmd -info` succeeding, plus all three taxdb files present.
- Missing `taxonomy4blast.sqlite3` means silent wrong answers with exit 0, because `-taxids` then returns unfiltered rows.
- Record real versions in `DECOY_VERSION.txt`: the `.tar.gz.md5` fetched **from source**, HTTP `Last-Modified`, and `blastdb-metadata-1-1.json` `last-updated` and `number-of-sequences`.

- [ ] **Step 2: Patch `package_release.sh`**

Keep its shape (`cd $SRC`, flat tar, allowlist by extension are all correct). Add: an explicit taxdb file list `taxdb.btd taxdb.bti taxonomy4blast.sqlite3`; package `DECOY_VERSION.txt`; assert no `*.tar.gz` is left in `$SRC`; cross-check `$VER` against the `manifest.json` `db_version` and the directory basename.

- [ ] **Step 3: Replace the text-grep tests with a real run**

`tests/testthat/test-db-release.R` currently greps the script text. Replace those with a real invocation against a mock directory of correctly-named empty files, which runs in about a second. Keep `skip_if_not(file.exists(p))` since `dev/` is `.Rbuildignore`d.

- [ ] **Step 4: Commit**

```bash
git commit -m "feat: verify decoys after fetch and package them with real versions"
```

---

### Task 7: Docker wiring

**Files:**
- Modify: `docker/deploy-local.sh`, `docker/deploy-dockerhub.sh`
- Modify: `docker/Dockerfile`
- Test: `tests/testthat/test-db-release.R`

- [ ] **Step 1: Delete stale staged tarballs before every build**

Add `rm -f docker/ribopilot_rrna-*.tar.gz docker/decoys-*.tar.gz` to `deploy-local.sh` and `deploy-dockerhub.sh`. Both already do exactly this for the package tarball. A sibling project shipped a mislabelled image because a stale glob match won by name order, so this is not hypothetical.

- [ ] **Step 2: Verify the database inside the image at build time**

Docker's `ADD` does not fail on a corrupt-but-gzip-shaped archive; it silently falls back to a plain copy. Add a verification `RUN` immediately after the `ADD` layers, in the `db-1` stage:

```dockerfile
RUN blastdbcmd -db /opt/ribopilot/blastdb/${RP_DB_VERSION}/ribopilot_rrna -info > /dev/null && \
    test -s /opt/ribopilot/blastdb/${RP_DB_VERSION}/manifest.json
```

No pipe: a pipe returns 141 on success and would pass for the wrong reason.

- [ ] **Step 3: Pin the ARG default to the R constant**

Add the `expect_equal` test from Task 1 Step 1 asserting `docker/Dockerfile`'s `ARG RP_DB_VERSION` default equals `RP_DB_VERSION` in `R/db_dist.R`.

- [ ] **Step 4: Verify without publishing**

`bash -n` both deploy scripts. Build with `--build-arg RP_DB=0` to confirm the no-database path still works, and confirm the verification `RUN` rejects a deliberately corrupt tarball. Do not push anything.

- [ ] **Step 5: Commit**

```bash
git commit -m "feat: verify the baked database and clear stale staged tarballs"
```

---

## Deferred, with reasons

Carried from spec section 7, all deliberate:

- **Panel build in CI or in `docker build`.** 3.5 h, 1.2 GB, ~600k E-utilities records, and NCBI content changes daily so it is not reproducible anyway. CI gets the cheap half only: `package_release.sh` against a mock directory.
- **Byte-identical reproducibility.** Impossible: +741 records in two days, and `accession2taxid` publishes no dated archive. The contract is to record the fingerprint, compare on re-run, and require `--accept-drift`.
- **Publishing.** `--publish` stays manual. Note the correction: "the private repo is a blocker" was **wrong**, `gh release download` works authenticated. The real reason no release exists is that nobody ran `fetch_decoys.sh` against the panel directory; this plan wires that in so the panel dir is release-ready by construction.
- **Runtime taxdump pinning and setting `BLASTDB`** so the shipped taxdb is actually consulted. Today taxdb is shipped, incomplete, and unreachable under RiboPilot's absolute-`-db` layout, and output is byte-identical with and without it. Fix it or drop it (~300 MB) as a separate deliberate decision.
- **Any change to what the panel contains.** EST/GSS/PDB filters, a RefSeq-preferring dedup tie-break, dropping taxdb: each is a new panel version with its own `EXPECTED_COUNT`.

## Out of scope, but a live bug worth its own branch

`inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh:126,150` still pipes through `efilter -molecule genomic`. That is the same filter measured to remove 100% of RefSeq records, so **user-built clade databases are losing their RefSeq records right now** (`R/build_reference_db.R:57`, `inst/nextflow/bin/rp_ref_db.R:18`). Not part of this plan. It needs its own branch and its own decision about whether existing clade databases should be rebuilt.

## Self-Review

**Spec coverage.** Section 1 (original sources) is Tasks 2, 3, 6. Section 2 (the scripts) is Tasks 1 through 6. Section 3 (definitional constants) is in Global Constraints and pinned by a Task 1 test. Section 4 (guards) is distributed across Tasks 2 through 6, each guard placed in the stage it protects. Section 5 (version sync) is Tasks 1 and 7. Section 6 (cost and resumability) is Tasks 1 and 4. Section 7 (not automating) is the Deferred list.

**Placeholders.** None. Where this plan describes rather than quotes an implementation, it is mechanical bash whose behaviour is fully specified by the guard it must satisfy; every guard states its exact condition and the historical failure it replaces.

**Type consistency.** `lib.sh` helpers are named identically everywhere they appear. `rp_marker_ok <stage> <input-hash>` takes the same two arguments in Tasks 1 and 4. Stage names are the same six strings in the orchestrator, the `--stage` flag, and the marker files.

**Ordering.** Task 1 first; it creates the file every later task edits. Tasks 2 through 5 are strictly sequential, both because they share `build_refdb.sh` and because each consumes the previous stage's output. Task 6 depends on Task 5's build directory layout. Task 7 is independent of 2 through 6 and may run any time after Task 1.
