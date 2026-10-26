# GenBank BLAST check (portal-parity BLAST flags)

Date: 2026-10-26
Branch: feat/genbank-blast-check (off devel at 58b5635)
Status: approved 2026-10-26, not yet built

## 1. Goal

- Predict, before submission, the two BLAST flags GenBank's rRNA/ITS curation
  raises on a record ("Query coverage is less than 75%" and "HSP_problem"),
  so a curator sees them in RiboPilot first and can tell gb-admin up front
  which records were reviewed.
- An additional layer beside the existing checks (screen, chimera windows,
  per-feature BLAST, stock ITSx check). It replaces nothing and gates nothing.
- Modeled on the stock ITSx check (R/itsx_check.R): run the way NCBI runs it,
  store the result per contig, rerun on the recheck, report it on export.

## 2. What GenBank does (email chain, 2026-10-09 to 2026-10-22)

- Source: Andrea Gocke (NCBI), replies to Dan MacGuigan and Matt Girard about
  SUB16547976. Page: https://www.ncbi.nlm.nih.gov/genbank/sequencecheck/
  (16S-centric; no rRNA-ITS thresholds there).
- BLAST databases: `ITS_eukaryote_sequences` (all eukaryotic ITS) and
  `rRNA_blastdb/rRNA_blast` (eukaryotic SSU and LSU, nuclear, mito, chloro).
  Command line or FTP only. ftp.ncbi.nlm.nih.gov/blast/db/ publishes
  `ITS_eukaryote_sequences`, `SSU_eukaryote_rRNA` and `LSU_eukaryote_rRNA`.
  The `rRNA_blast` name is internal; the user's ruling (2026-10-26) is that
  it is the SSU plus LSU eukaryote pair, and the check treats it so.
- Provenance: all three are NCBI RefSeq Targeted Loci "selected GenBank"
  sets: ordinary GenBank accessions that passed NCBI's validation (ITSx for
  the ITS set). The ITS set (77,582 records, median 654 bp, build 2026-10-12)
  is mostly plant and fungal ITS amplicons; 157 records are bony fishes, 584
  stony corals. The SSU and LSU sets (8,784 and 6,575 records) date from
  2022-12-05. NCBI refreshes the ITS set regularly, so a baked copy lags and
  NCBI's own result can drift from RiboPilot's.
- "Default BLAST settings". Two flags per record:
  - low coverage: query coverage of the top hit under 75%.
  - HSP_problem: the top hit has more than one HSP ("potential misassembly").
- Both flags only send a record to a curator. A submitter reply that the
  sequences are correct moves it on to "additional tools for verification".

## 3. Replication (2026-10-26, scratchpad, not in the repo)

- Queries: the review-layer sequence of each contig (what export writes),
  from the Coral_Triangle project database (80 locked, flag-free contigs,
  including the 8 SUB16547976 records NCBI flagged).
- `blastn` 2.17 defaults (megablast), `-evalue 0.05 -max_target_seqs 100`,
  no `-max_hsps`. Top hit = first subject blastn reports. Coverage = `qcovs`.
- Against `ITS_eukaryote_sequences`:
  - HSP_problem: 8 of 8 NCBI-flagged records have 2 or 3 HSPs on the top hit.
  - low coverage: 6 of the 7 NCBI-flagged records are under 75% (70 to 74%).
    07-01-16_S48 scores 76% here; NCBI called it under 75%. 07-01-122_S154,
    which NCBI did not flag for coverage, scores 79%.
  - Across all 80 contigs: 42 under 75% coverage, 76 with more than one HSP.
    The flags are the norm for fishes, not the exception: the ITS database
    has few fish records, so the top hit is a distant fish and ITS1/ITS2 do
    not align, splitting the hit into 18S, 5.8S and 28S HSPs.
- Against `SSU_eukaryote_rRNA` plus `LSU_eukaryote_rRNA` (whole contig):
  every one of the 80 contigs has one HSP and 25 to 65% coverage, because
  each subject record holds one subunit and the query holds two. Applying the
  75% rule to this search would flag all 80; NCBI flagged 7. So NCBI's
  coverage flag does not come from this search as run here.
- Per-region queries (RiboPilot's 18S vs SSU, 28S vs LSU, ITS1-5.8S-ITS2 vs
  ITS set) do not reproduce the pattern either: 18S and 28S give 98 to 99%
  coverage and one HSP for every contig; the ITS region gives one HSP and
  18 to 61% coverage for the flagged eight and under 75% for 70 of 80.
- On the exact submitted FASTA of SUB16547976 (~/Downloads/test_data,
  export label 28SMissingNOITSX_10xTest_Oct26, 10 records of 2837 to 2940 bp),
  the ITS search reproduces NCBI's reply exactly:

  | record | qcovs | HSPs | NCBI low coverage | NCBI HSP_problem |
  |---|---|---|---|---|
  | 07-01-122_S154 | 79 | 2 | no | yes |
  | 07-01-147_S179 | 70 | 2 | yes | yes |
  | 07-01-16_S48 | 69 | 2 | yes | yes |
  | 07-01-61_S93 | 100 | 1 | no | no |
  | 07-01-65_S97 | 70 | 2 | yes | yes |
  | 07-01-77_S109 | 70 | 2 | yes | yes |
  | 07-01-94_S126 | 70 | 2 | yes | yes |
  | 07-02-111_S335 | 100 | 1 | no | no |
  | 07-02-118_S342 | 74 | 2 | yes | yes |
  | 07-02-93_S317 | 71 | 2 | yes | yes |

  The 76% for 07-01-16_S48 in the 80-contig run came from an older, longer
  (5805 bp) version of that contig. The rRNA search gives 62 to 65% and one
  HSP for all ten. In every flagged record the top ITS hit splits into an 18S
  HSP (positions 1 to ~1850) and a 5.8S HSP (~2200 to 2550); ITS1, ITS2 and
  the 28S stub do not align.
- Conclusion: one megablast of the whole export sequence against
  `ITS_eukaryote_sequences` reproduces NCBI's flags exactly on the one
  submission with a known answer.
  The rRNA (SSU plus LSU) search is run and reported beside it because NCBI
  runs it, but its coverage cannot carry the 75% rule.

### 3.1 Validation fixture

- The ten SUB16547976 records are the acceptance test. They go into the
  repo as tests/testthat/fixtures/genbank_check/: the submitted FASTA, the
  two blastn outfmt-6 result files from the 2026-10-26 databases (ITS and
  rRNA), and a small expected-flags table transcribed from NCBI's 2026-10-15
  reply (the two lists above).
- A unit test feeds the canned result files through `genbank_check_calls()`
  and asserts the stored text flags exactly NCBI's 7 and 8 and leaves
  07-01-61_S93 and 07-02-111_S335 clean. No blastn needed.
- An opt-in real-run test, gated by an env var naming a directory that holds
  the three databases (same pattern as the ITSx real-run test), runs
  `genbank_check_unit()` on a scratch project seeded with the ten sequences
  and asserts the same lists; it also records the database build dates it
  used, since a newer ITS set can move a borderline record (07-02-118_S342
  sits at 74%).

## 4. Design

### 4.1 Database shipping

- `ITS_eukaryote_sequences` joins the decoy set that ships beside the panel:
  dev/refdb/fetch_decoys.sh `SEARCHABLE_DBS`, dev/refdb/build_refdb.sh
  `DECOY_DBS`, dev/refdb/package_release.sh, and the Dockerfile verify step.
  It is installed by `install_ribopilot_db()` with the other decoys and found
  by name in the panel's directory (`dirname(panel)`), as the decoys are.
- It is NOT added to `RP_DECOY_DBS` (R/blast_local.R). The whole-contig screen
  and the per-feature BLAST keep searching exactly what they search today.
  `SSU_eukaryote_rRNA` and `LSU_eukaryote_rRNA` already ship as decoys, so
  the rRNA half of the check needs nothing new on disk.
- Needs a new panel/decoy release and an image rebuild, made the normal way:
  a full dev/refdb/build_refdb.sh run under a new version. The panel content
  refreshes from GenBank as a side effect; that is accepted (no re-stamping
  path for a decoys-only release is added). `RP_DB_VERSION` in R/db_dist.R
  moves to the new version.
- tests/testthat/test-db-release.R asserts the shell decoy lists equal
  `RP_DECOY_DBS`; it changes to a shipped list (`RP_DECOY_DBS` plus
  `ITS_eukaryote_sequences`) versus the searched list.
- Without the database on disk (older image, `--no-db` image, panel installed
  before this release) the check stores `not_run(no ITS_eukaryote_sequences
  database)` and nothing else changes.

### 4.2 The check: R/genbank_check.R

- `genbank_check_unit(con, id, path, scaffold, db_dir, workdir, threads,
  review, blastn)`. `db_dir` holds the three NCBI databases (the panel's
  directory). Query: `assemblies.sequence` for the annotate layer,
  `review_units_v.sequence` for the review copy (what export writes).
- Two searches of the whole query, one per NCBI database: "ITS" against
  `ITS_eukaryote_sequences`, "rRNA" against `SSU_eukaryote_rRNA` and
  `LSU_eukaryote_rRNA` together (one `-db` with both prefixes, as
  `local_blast_db_arg()` composes).
- blastn: default task (megablast), `-evalue 0.05 -max_target_seqs 100`,
  outfmt 6 with `saccver stitle pident qcovs qcovhsp evalue bitscore length
  qstart qend sstart send`. No `-max_hsps`. Constants in the file:
  `GENBANK_CHECK_MIN_QCOVS <- 75`, `GENBANK_CHECK_MAX_HSPS <- 1L`.
- A pure function `genbank_check_calls(its_hits, rrna_hits)` turns the two
  parsed hit tables into the stored text, so it is unit-testable on canned
  outfmt-6 files. Per search: top hit = the first subject reported, its
  `qcovs`, its HSP count and the HSP query spans.
  - ITS search: no hits `genbank_no_hit(ITS)`; `qcovs < 75`
    `genbank_low_coverage(ITS <qcovs>%)`; HSPs > 1
    `genbank_multi_hsp(ITS <n> HSPs: <qstart>-<qend>, ...)`.
  - rRNA search: no hits `genbank_no_hit(rRNA)`; HSPs > 1
    `genbank_multi_hsp(rRNA <n> HSPs: ...)`. Its coverage is reported, not
    flagged (section 3): a two-subunit query against one-subunit records
    cannot reach 75%, and NCBI evidently does not flag it. Author's call; a
    one-line change if NCBI says otherwise.
  - Clean: "" (empty), like `itsx_warnings`.
  - Flags joined with "; ", then a trailing note
    ` top hits: ITS <acc> <species> <qcovs>%, rRNA <acc> <species> <qcovs>%`
    whenever either search hit, so the stored text stands on its own.
  - Failures: `not_run(<reason>)` (database missing, blastn exit status, no
    sequence). One missing database fails the whole check; the two ship
    together.
- Stored in one new TEXT column, `genbank_check`, on `annotate` and on
  `review_units` (review view `review_units_v` exposes it like
  `itsx_warnings`: the copy's value, falling back to the original's when the
  copy has none). Schema step 6 (version 6 to 7) in `db_migrations()`, DDL in
  `new_db()`; the backwards-compat parity test covers it. Existing contigs
  read `NULL`, which displays and exports as `not_run(not checked)`. No
  backfill: like the ITSx check at the review-copy upgrade, a locked contig
  keeps `not_run` until it is unlocked and rechecked (the recheck queue takes
  unlocked contigs only). Newly annotated contigs get the check from the
  pipeline.

### 4.3 Where it runs

- Pipeline: inside the existing `feature_blast` process
  (inst/nextflow/bin/rp_feature_blast.R), after `feature_blast_unit()`. The
  process already has the panel path, taxdump and blastn in the container;
  `db_dir` is `dirname()` of the `DBARG` panel prefix (first token), where
  the decoys live. No new process. When
  `DBARG` is empty (`RIBOPILOT_NO_LOCAL_BLAST`, or no panel found) the check
  stores `not_run(no local BLAST)`, the same condition that skips the
  per-feature BLAST; it does not look for the database on its own.
- Recheck: `recheck_unit()` (R/recheck.R) calls it after `feature_blast_unit()`
  with `review = TRUE`, wrapped like the ITSx check so a failure stores
  `not_run(...)` and does not fail the recheck.
- `review_follow_original()` already copies the annotate row's values to an
  unedited copy; the new column follows the same path as `itsx_warnings`.

### 4.4 Where it shows

- Checks column (`unit_checks()`, R/export.R): a new amber entry "GenBank
  BLAST" when `genbank_check` is non-empty and not `not_run(`; never red,
  never read by `has_screen_risk()`. No "not checked" entry (unlike ITSx):
  with most fish contigs flagged, the column must stay readable.
- Scaffold table and Annotation Review table: a `genbank_check` column next
  to `itsx_warnings`, same badge style (`rt_itsx_badge()` reused or
  generalised), hidden by default in the "Check details" group, with a
  help header.
- Export is the primary surface. The point is that the curator sees, at the
  moment of export, which records GenBank is likely to send back and why, so
  the reply to gb-admin can be prepared before the rejection arrives.
  - Export window (`export_check_summary_ui()`, R/app_modals.R): a line
    "GenBank BLAST check: n flagged, n clean, n not checked", amber and bold
    when any are flagged, naming the flagged records with their flags the
    way the ITSx line does (`mp_id_list`, first 10), e.g.
    `07-01-65_S97 (low coverage ITS 70%, 2 HSPs)`. A short display form of
    the flags is used here, not the stored text.
  - `<label>_summary.csv`: `genbank_check` column beside `itsx_check`, the
    full stored text. The HTML export report (inst/export_report.Rmd) adds
    the same column to its record table.
  - README.txt: a "GENBANK BLAST CHECK" section, after the ITSx section,
    with: what NCBI runs (the two databases, default megablast); the two
    flags and the exact wording NCBI uses in its reply ("Query coverage is
    less than 75%", "HSP_problem"); why most fish records trip them (few
    fishes in the ITS set, spacers do not align to distant relatives, so the
    hit splits into 18S, 5.8S and 28S pieces); that NCBI's copy of the ITS
    set is newer than RiboPilot's so NCBI can differ on borderline records;
    and what to do when the submission comes back: reply to gb-admin that
    the listed records were reviewed with RiboPilot's checks and are
    believed correct, after which NCBI runs further verification. The
    flagged seqids are listed in that section so the reply can be pasted
    from it.
  - `not_run(...)` records are counted as "not checked" and the README says
    how to check them again, as for ITSx.
- vignettes/GenBank-Submission.Rmd: a short block after "The ITSx check"
  with the same content as the README section.

### 4.5 Tests

- The SUB16547976 fixture tests (section 3.1) are the acceptance criterion:
  the implementation is done when they pass.
- `genbank_check_calls()` on canned hit tables: no hit in either search,
  clean, ITS low coverage only, ITS multi-HSP only, both, rRNA multi-HSP,
  rRNA low coverage not flagged, HSP span ordering.
- `genbank_check_unit()` with a fake blastn (tests/testthat/helper-fake-tools.R
  pattern): stores text; missing database stores `not_run(...)`; review = TRUE
  writes the copy.
- Schema parity (existing test-backwards-compat.R) and the nf-bin contract
  test for rp_feature_blast.R if its arguments change.
- `unit_checks()` and `export_summary_df()` pick the column up.

## 5. Non-goals

- Not a gate: Lock and Export do not ask about it.
- No per-region (ITSx-split) searches (section 3).
- No attempt to reproduce NCBI's "additional tools for verification", Uchime,
  or vector screening.
- No per-feature variant of this check; the per-feature BLAST already exists.
- No change to the panel, to `RP_DECOY_DBS`, or to the whole-contig screen.

## 6. Decisions (user, 2026-10-26)

1. `ITS_eukaryote_sequences` ships in the decoys tarball of the next panel
   release. No panel release or image rebuild in this work: the dev/refdb
   scripts, `RP_DB_VERSION` and the Dockerfile are updated so the next release
   picks it up, and until then the check stores `not_run(no
   ITS_eukaryote_sequences database)`. Local testing points `db_dir` at a
   directory holding the three databases.
2. Amber warning only, never red, never a gate.
3. No backfill for contigs already locked in existing projects.

## 7. Left to the builder

- Exact column placement and help text wording.
- Whether `rt_itsx_badge()` is reused or a general badge helper is extracted.
- The panel release version string and the image rebuild sequence
  (dev/refdb/package_release.sh, docker/deploy-local.sh).
