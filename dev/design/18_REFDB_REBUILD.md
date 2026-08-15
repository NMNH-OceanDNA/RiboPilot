# CONSOLIDATED REFDB REBUILD - SPECIFICATION

## 1. ORIGINAL SOURCES

| Data | Exact original source | Auth | Version capture |
|---|---|---|---|
| Panel accession list | EDirect `esearch -db nuccore -query "$QUERY"` -> eutils.ncbi.nlm.nih.gov. Capture XML once, parse `<Count>`, replay same WebEnv/QueryKey into `efetch -format acc` | `NCBI_API_KEY` env, required, fail fast if unset (3 vs 10 req/s, silent today) | `<Count>`, ISO date, sha256 of sorted list -> manifest |
| Panel FASTA | `efetch -db nuccore -input chunk -format fasta`, 5000 acc/chunk | same | per-chunk requested/returned accession sets; funnel counts -> manifest |
| accession->taxid | `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/{nucl_gb,nucl_wgs,pdb}.accession2taxid.gz`, streamed `curl --fail -sL \| zcat \| awk`, never landed (~9.4 GB, not the 2.5 GB the old comment claims) | none | fetch-time `Last-Modified`, `Content-Length`, source `.md5` sidecar -> manifest. No dated archive upstream: exact refetch impossible by design |
| Residual taxids (newer than dump) | `epost -db nuccore \| esummary \| xtract` | `NCBI_API_KEY` | count + per-row source tag (`ftp` vs `esummary`) -> manifest |
| Decoys 16S/SSU/LSU + taxdb | `update_blastdb.pl --source ncbi --decompress` -> `https://ftp.ncbi.nlm.nih.gov/blast/db/` | none | per-db `.tar.gz.md5` **fetched from source**, HTTP `Last-Modified`, `blastdb-metadata-1-1.json` `last-updated`/`number-of-sequences` -> `DECOY_VERSION.txt` (packaged) + manifest decoy block |
| Runtime taxdump | `R/taxonomy.R:90`, per project at runtime | none | out of scope, see S7 |

Currently taken from derived/cached locations, must stop:
- `~/.cache/ribopilot-blastdb` is a dev fixture (`dev/make_screen_fixture.R:214`), never a source. Its 16S is already stale (27,649 seq vs 27,644 live).
- On-disk `<db>.tar.gz.md5` is update_blastdb.pl's **cache key**, not evidence. With marker present and volumes deleted it prints "up to date", installs nothing, exits 0. The original claim that it records "what landed" was wrong.
- `/home/dmacguig/edirect` and the absolute shim path (`01_harvest.sh:14-15`, `02b:11-12`). Resolve shim from script location; take EDirect from PATH.

## 2. THE SCRIPTS TO WRITE (all under `dev/refdb/`)

- **`build_refdb.sh`** - the orchestrator and the whole panel build, one file, stages run in-process (no chain wrapper, no `pgrep`, no process-name matching anywhere; liveness via PID lockfile).
  `build_refdb.sh [--version V] [--build-dir D] [--stage NAME] [--jobs N] [--accept-drift] [--force] [--check-version] [status|clean]`
  Stages, in order: `accessions, taxidmap, fetch, build, decoys, package[, publish]`. `taxidmap` runs concurrently with `fetch` (depends only on accession list; original build overlapped them 12:03-12:07 inside the 12:03-14:23 fetch window). Default runs all. **Resumable**: `fetch` only, at chunk granularity. **Idempotent**: `accessions`, `taxidmap`, `build`, `package` (re-run from scratch, cheap). `publish` is the only destructive-network stage and is never in the default run.
- **`fetch_decoys.sh`** (rewrite in place) - `fetch_decoys.sh DEST`. Adds `--source ncbi` (skips ~30 s AWS/GCP metadata probes and stops silent mirror switching), keeps taxdb **last** with a comment saying why (every decoy tarball bundles a 2024 taxdb that would otherwise overwrite the current one), post-fetch verification, and real version capture. Idempotent, not resumable.
- **`package_release.sh`** (patch) - unchanged shape (`cd $SRC`, flat tar, allowlist by extension - all correct), plus: explicit taxdb file list `taxdb.btd taxdb.bti taxonomy4blast.sqlite3`, package `DECOY_VERSION.txt`, assert no `*.tar.gz` left in `$SRC`, cross-check `$VER` against `manifest.json` db_version and dir basename.

Deleted: `01_harvest.sh`, `01b_repair.sh`, `02_taxidmap.sh`, `02b_taxid_fallback.sh`, `03_build.sh`, `status.sh` (508 lines). `01b` exists only to undo `01`'s strict check - folded into one tolerance policy. `status` becomes `build_refdb.sh status`, reading a small state file (no 1.2 GB re-scan, no hard-coded 179, no `[ -f blastdb/*/manifest.json ]` glob that returns FALSE as soon as two version dirs exist).

## 3. VERBATIM DEFINITIONAL CONSTANTS

Define **once** in `build_refdb.sh`, echo expanded into manifest at build time:

```
TAXID=33208
RRNA_TERMS='18S ribosomal RNA[Title] OR 28S ribosomal RNA[Title] OR 5.8S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title] OR large subunit ribosomal RNA[Title] OR internal transcribed spacer[Title]'
QUERY="(txid${TAXID}[Organism:exp]) AND (${RRNA_TERMS}) NOT mitochondrion[Filter] NOT plastid[Filter]"
```

Also load-bearing, change = new panel version, never a silent fix:
- **No `efilter -molecule genomic`, ever.** Measured live: query alone 892,123; through the filter 741,526. It kills 100% of RefSeq (148,365, MolType=rna) plus 2,232 non-RefSeq. Do not adapt `inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh` (`:126`, `:150`).
- Chunk size 5000. Dedup key `(taxid, uppercased sequence)`, byte-exact. Survivor tie-break stays **lexicographic by accession under `LC_ALL=C`** - documented as-is, not "improved" to prefer RefSeq, because that changes 95,997 records of panel content.
- makeblastdb flags exactly: `-dbtype nucl -parse_seqids -blastdb_version 5 -taxid_map <map> -title "RiboPilot Metazoa nrDNA <VER>"`. Serialize into manifest as a JSON array (current `${MKCMD[*]}` is not re-runnable).
- Decoy names `16S_ribosomal_RNA SSU_eukaryote_rRNA LSU_eukaryote_rRNA taxdb`, one definition shared with `R/blast_local.R:11`.
- `export LC_ALL=C` at top (underscored accessions sort differently under en_US.UTF-8; `comm` then silently lies).
- Off-target records (1,394 EST/cDNA "similar to ... rRNA", 298 PDB) are **in** the panel by definition. Adding `NOT est[Filter] NOT gss[Filter] NOT pdb[Filter]` is a new panel version, not a bugfix.

## 4. GUARDS

Harvest / accessions
- Harvested line count must equal esearch `<Count>` exactly (verified on 3 clades: Count == lines == unique). Replaces the `-lt 500000` floor, which sits 44% below truth and passes a stream dying at 600k.
- `<Count>` compared to committed `EXPECTED_COUNT` +/- tolerance; beyond it, refuse unless `--accept-drift`, which records the new count. Source really moves: 891,382 (Aug 13) -> 892,123 (Aug 15), ~370/day.
- Never `$(grep -c . f || echo 0)` - on an existing empty file it yields `"0\n0"`, `[` errors rc=2, guard skipped, empty list promoted as success. Use `wc -l`. Test every numeric guard against empty input.
- Filter `^[A-Z]{4}0{6}` (WGS master accessions) up front as a recorded expected drop - the 11 "unavailable" records are permanent by design, not transient.

Fetch
- Chunk complete = requested vs returned **accession set** comparison (after normalizing `pdb|XXXX|C` -> `XXXX_C`), plus no zero-length sequences, plus trailing newline. Never `got >= want` on `>` counts: a stream cut after the last header passes it. Write `.tmp`, verify, atomic rename; leave nothing in the chunk dir (the `fa_*` glob eating `.tmp` files is the only reason the original run did not lose 15,000 records).
- Retry **only the missing accessions**, never refetch or discard a chunk. One tolerance policy, not `100%` in `01` and `95%` in `01b`.
- **Never branch on efetch exit status.** It exits 0 with 0 bytes and `QUERY FAILURE` on total failure, which made `01b`'s entire `&& break` retry loop inert.
- Chunk cache keyed to sha256 of the accession list, stored beside `chunks/`. Mismatch = refuse to resume. Positional `split` re-cutting a changed list silently keeps withdrawn records and silently misses new ones; name-pairing does not catch it.
- Reverse reconcile: got-not-in-want must be empty (never computed today; the shrink case shipped 12 withdrawn records in simulation).
- Assert `command -v curl` resolves inside the repo shim dir before any EDirect call; assert EDirect and `update_blastdb.pl` on PATH. The shim mitigates but does not eliminate curl 56 (it happened with the shim on PATH), so retries stay. PATH prepending is the correct and only lever for EDirect: `nquire:403` does `command -v curl`. Note the original claim that the fix was "guarantee curl at /usr/bin/curl" was wrong - conda curl 8.21 fails identically under `--http1.0`; the shim's `--http1.1` is the fix.
- Scrub `api_key=[A-Za-z0-9]+` from everything written to a log (35 plaintext hits in `harvest.log`). `.gitignore` the build dir.

Taxonomy
- `curl --fail`, drop `2>/dev/null` on zcat, check scan exit status. gzip's own CRC catches truncation and corruption for free; `--fail` catches the 404-HTML-into-the-pipe case. A 40% truncation silently yields 356,851 taxid-0 records and 5,457 wrongly deleted ones.
- Mapped rows >= 99.5% of requested accessions, else abort (historical coverage was 100%).
- Map validator, in this order: every line `^[A-Za-z0-9._]+\t[0-9]+$`; map accession set is a **subset of the harvest list**; unique accessions == row count. The subset check is the one that catches the two phantom log lines - "row count == unique count" passes on the corrupt map and does not.
- No-taxid records are dropped with reason `taxid_unresolved` in `dropped.tsv`. Never a shared `"NA"` dedup key: it collapses different species with identical rRNA, mislabelled as `byte_identical_to_kept_record_same_taxid`.

Build / package
- Strict glob `fa_[0-9][0-9][0-9][0-9]` and assert chunk count > 0. Empty dir currently passes `0 -ne 0`; a stray non-empty `fa_retry` alone builds and ships a 3-sequence database with a full manifest.
- makeblastdb exit status proves nothing: truncated FASTA -> "added 2 sequences", exit 0; a taxid map matching nothing -> prints `Error:`, exit 0, every staxid 0. Parse the generated `.njs`: `number-of-sequences` == FASTA records, `number-of-letters` == seqkit sum_len, `number-of-taxids` == distinct taxids in the restricted map. Derive `has_taxids` from a `blastdbcmd -outfmt '%T'` sample, never hardcode `true`.
- Ledger identity `requested == kept + sum(dropped by reason)` asserted before writing manifest; every count computed from final artifacts. Today it is off by 188 and `dropped.tsv` has 110 more rows than `records_dropped`. `unpacked_bytes` measured over the shipped file list only (currently reports 1.14 GB for a 331 MiB database).
- Decoys verified **after** fetch, per db: required extensions non-empty + `blastdbcmd -info`, plus all three taxdb files present. Missing `taxonomy4blast.sqlite3` = silent wrong answers with exit 0 (`-taxids` returns unfiltered rows).
- `set -Eeuo pipefail` + ERR trap + cleanup trap, **with `|| true` on every grep-as-filter**. Bare errexit kills the real script at `03_build.sh:101` whenever zero PDB records match, which is a legitimate state.
- Explicitly **not** doing: no per-stage count-equality asserts (the good build had 309 want-not-in-got and 298 got-not-in-want; equality aborts a correct build). No `>=15% RefSeq` assert - that figure is the harvest share, panel share is 8.73%, and every clade DB is 0.4-11.9%, so it false-fires on all of them.

## 5. VERSION SYNCHRONISATION

`R/db_dist.R:86 RP_DB_VERSION` is the single source of truth. `build_refdb.sh` reads it; `--version` must match or the run refuses. Same string must appear in: output dir basename `blastdb/<VER>`, `manifest.json` db_version, BLAST `-title`, `ribopilot_rrna-<VER>.tar.gz` / `decoys-<VER>.tar.gz` (+ `.sha256`, format `<hex>  <basename>` - both `sha256sum -c` and `R/db_dist.R:13-24` depend on that shape), release tag `db-<VER>`, `docker/Dockerfile:99 ARG RP_DB_VERSION`. `--check-version` diffs all six and exits non-zero on disagreement.

Format constraint: dotted numeric only. `numeric_version("v2026.08.1")` and `"2026.08.1-rc1"` are both NA, which makes the panel invisible to `discover_panel()` (`R/blast_local.R:87-105`). Ordering is correct without zero-padding.

Refuse to overwrite an existing `blastdb/<VER>` without `--force` (`03_build.sh:26` is bare `mkdir -p`; the directory was silently rebuilt in place once already, 598,098 then 597,988 records).

The Dockerfile ARG default is **not** the desync risk the original analysis claimed - 2 of 3 build cases fail loudly and the third is internally consistent. Fix the real defect instead: add `rm -f docker/ribopilot_rrna-*.tar.gz docker/decoys-*.tar.gz` to `deploy-local.sh:36` and `deploy-dockerhub.sh:70` (both already do this for the package tarball), and add one `expect_equal` in `tests/testthat/test-db-release.R` asserting the ARG default equals `R/db_dist.R:86`.

## 6. COST AND RESUMABILITY

Measured full rebuild: **2 h 41 m**, 99% of it the chunk fetch (179 chunks, ~47 s each = 2 h 20 m). Accession list 19 m, taxid map 4 m (concurrent), fallback 15 s, build 27 s (makeblastdb itself 6.3 s), decoys ~291 MB download. With `--jobs 4` and an API key (10 req/s), fetch should drop to roughly 40 min; that is the one optimization worth having.

Disk: 3.4 GB build tree + up to 2 GB sort spill -> **require 8 GB free**, checked at startup. `sort -S` and `TMPDIR` configurable (hard-coded `-S 4G` is wrong on a shared node). Peak includes ~1.2 GB of pure duplicate (`metazoa_rrna.raw.fasta`, regenerated from chunks every run).

Restart after failure: re-run the identical command. Completed stages skip on marker + input-hash match; `fetch` resumes per chunk and re-requests only missing accessions. `--stage fetch --force` re-runs one stage. If the accession list hash changed, the run refuses to resume and says the source moved - that is the honesty requirement, not a bug. `clean` drops `chunks/`, raw FASTA, and the unshipped `.fasta`/`.taxid.map`/`.dropped.tsv` once the tarballs verify (ship `dropped.tsv.gz` as a release asset; it is the only record of the 293,692 discarded records).

## 7. WHAT I AM NOT AUTOMATING

- **Panel build in CI or in `docker build`.** 3.5 h, 1.2 GB, ~600k E-utilities records, and NCBI content changes daily so it is not reproducible anyway. It would also bust the deliberate layer cache at `docker/Dockerfile:88-90`. CI gets only the cheap half: run `package_release.sh` against a mock directory of correctly-named empty files (~1 s), replacing the text-greps in `tests/testthat/test-db-release.R:5-16,31-32`.
- **Byte-identical reproducibility.** Impossible and I will not pretend otherwise: +741 records in two days, and `accession2taxid` publishes no dated archive. Contract = record the fingerprint, compare on re-run, require `--accept-drift`. Drift *reporting* is a diff of the new accession->taxid map against the previous release's map (dumped in 1.4 s via `blastdbcmd -entry all -outfmt '%a\t%T'`, which reproduces the build-time map md5 exactly). An alarm on the upstream `Last-Modified`/md5 would fire on every rebuild while content moves 0% (0 of 300 sampled taxids changed across two file turnovers).
- **Publishing.** `--publish` stays manual. The original "private repo is a blocker" claim was wrong - `gh release download` works authenticated and `package_release.sh:67` already creates the release. The real reason no release exists is that nobody ran `fetch_decoys.sh` against the panel dir; the consolidated build wires that in so the panel dir is release-ready by construction. Only `install_ribopilot_db()` needs public assets (`R/db_dist.R:103` uses unauthenticated `download.file`); decide at public-release time.
- **Runtime taxdump pinning** (`R/taxonomy.R:90`) and **setting `BLASTDB`** so the shipped taxdb is actually consulted. Both are app-level changes, not rebuild changes. Today taxdb is shipped, incomplete, and unreachable under RiboPilot's absolute-`-db` layout; output is byte-identical with and without it. Fix it or drop it (~300 MB) as a deliberate, separate decision.
- **Any change to what the panel contains**: query filters for EST/GSS/PDB, a RefSeq-preferring dedup tie-break, dropping taxdb. Each is a new panel version with its own `EXPECTED_COUNT`.
- **The clade-DB path.** `download_rrna.sh:126,150` still pipes through `efilter -molecule genomic`, so user clade DBs (`R/build_reference_db.R:57`, `inst/nextflow/bin/rp_ref_db.R:18`) are losing their RefSeq records right now. Real bug, separate branch, out of scope here.
- **Schema migrations** for projects built against an older panel - pre-release, none written by policy.