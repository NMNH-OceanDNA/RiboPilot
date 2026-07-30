# Screen Validation Fixture Set Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the labelled corpus doc 08 section 5 names as a gate, measure the screen's verdicts and its false-positive rate against it, and print the numbers, without changing any default.

**Architecture:** A dev-only generator captures hit distributions once (the clean-host half straight from a real project DB, the spike-in half through the pipeline's own remote BLAST) and commits two fixture files. An offline, deterministic test replays every fixture unit through `screen_channel_taxonomy()` and asserts per-unit verdicts plus the aggregate floors. A fixture unit is a hit distribution PAIRED WITH a host lineage; the same distribution appears twice against different hosts to measure the screen's boundary.

**Tech Stack:** R 4.4, DBI + RSQLite, testthat 3e, entrez-direct and blastn inside `ribopilot:test`.

**Spec:** `dev/design/13_VALIDATION_FIXTURES.md`

---

## Global Constraints

- ASCII only in all R code, tests, fixtures and docs. No em dashes, no smart quotes.
- Comments minimal, explaining WHY rather than narrating a change.
- `NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only. **This plan adds no exports**, so `NAMESPACE` must end unchanged.
- **No default changes.** `SCREEN_VERSION` stays `4L`, `RIBOPILOT_SCHEMA_VERSION` stays `15L`, `triage_opts.rule_screen_off_target` stays `'flag'`. A diff touching any of those is a defect.
- **No new package code.** Everything ships in `dev/` (excluded by `.Rbuildignore`, verified: 0 `dev/` entries in the built tarball) or `tests/`. `R/` is not modified.
- `inst/nextflow/**` is not modified. No pipeline process is added.
- The test reads ONLY the two fixture files. It must never open a socket, run `blastn`, or read a project DB.
- Verdict strings are exactly `on_target`, `off_target`, `no_reference`, `uncertain`. Class strings are exactly `bacterial`, `archaeal`, `fungal`, `algal_plant`, `other_euk`, `metazoan_other`, `human`.
- An `uncertain` verdict on a spike-in is NOT a catch. Counting it would let the 0.6 gate be tuned to look good by declining to answer.
- `%|NA|%` (`R/utils.R`) is scalar only and errors on a vector. `%||%` tests length, not NA.
- Run tests from the repo root: `Rscript -e 'testthat::test_local()'`. Baseline is **913 passing, 0 failing**.
- Branch is `validation-fixtures`, already created from `main` at `eb452c4`.
- After each task: revert the change, confirm a named test fails, restore, confirm the tree is clean.

---

## File Structure

**Create:**
- `dev/make_screen_fixture.R` - the generator. Dev only, never run by tests, never shipped.
- `tests/testthat/fixtures/screen_hits.csv.gz` - one row per HSP.
- `tests/testthat/fixtures/screen_labels.csv` - one row per unit.
- `tests/testthat/test-screen-validation.R` - the offline replay and the assertions.

**Modify:**
- `dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md` - section 5 gate split, and the human finding.
- `dev/design/11_SCREEN_VERDICT_LAYER.md` - section 1's copepod sentence.
- `dev/design/12_INTRAGENOMIC_VARIANTS.md` - record that 99.5 is still unvalidated.

**Not modified, deliberately:** all of `R/`, all of `inst/`, `NAMESPACE`, `DESCRIPTION`.

---

### Task 1: Capture the clean-host half

No network. Delivers the false-positive corpus, which is the number the gate is really about.

**Files:**
- Create: `dev/make_screen_fixture.R`
- Create: `tests/testthat/fixtures/screen_hits.csv.gz`, `tests/testthat/fixtures/screen_labels.csv`

**Interfaces:**
- Produces: `fixture_capture_host(db_path, taxdump_date = NA_character_)` -> (invisibly) number of units written.
- Produces the two fixture files with the exact column sets below. Later tasks append to them; nothing rewrites them from scratch.

`screen_hits.csv.gz` columns, in this order:
```
unit, hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
```

`screen_labels.csv` columns, in this order:
```
unit, category, expected_verdict, expected_class, source, sample_taxids, accession, taxdump_date, note
```

- [ ] **Step 1: Write the generator's host half**

Create `dev/make_screen_fixture.R`:

```r
#!/usr/bin/env Rscript
# Generator for the screen validation fixture (dev/design/13_VALIDATION_FIXTURES.md).
#
# DEV ONLY. Never run by tests and never shipped: dev/ is in .Rbuildignore, so
# this file is absent from the built package. Tests read the committed fixture
# and never reach the network.
#
# Usage:
#   Rscript dev/make_screen_fixture.R host <project_db>
#   Rscript dev/make_screen_fixture.R spikes <project_db> <taxdump_dir>
#
# Regenerating erases the provenance of numbers the tests assert against, so
# both modes refuse to clobber an existing fixture without --force.

suppressMessages(pkgload::load_all(quiet = TRUE))

FIX_DIR <- file.path("tests", "testthat", "fixtures")
HITS_F <- file.path(FIX_DIR, "screen_hits.csv.gz")
LAB_F <- file.path(FIX_DIR, "screen_labels.csv")

HIT_COLS <- c("unit", "hit_rank", "hsp", "bitscore", "hit_class",
  "lineage_taxids", "lineage")
LAB_COLS <- c("unit", "category", "expected_verdict", "expected_class", "source",
  "sample_taxids", "accession", "taxdump_date", "note")

#' Append rows to the fixture, creating it when absent
fixture_append <- function(hits, labels) {
  dir.create(FIX_DIR, recursive = TRUE, showWarnings = FALSE)
  stopifnot(identical(names(hits), HIT_COLS), identical(names(labels), LAB_COLS))
  if (file.exists(HITS_F)) {
    old_h <- utils::read.csv(gzfile(HITS_F), stringsAsFactors = FALSE)
    old_l <- utils::read.csv(LAB_F, stringsAsFactors = FALSE)
    # a re-captured unit REPLACES its old rows; keeping both would silently
    # double its hits and change every score derived from them
    old_h <- old_h[!old_h$unit %in% hits$unit, , drop = FALSE]
    old_l <- old_l[!old_l$unit %in% labels$unit, , drop = FALSE]
    hits <- rbind(old_h, hits); labels <- rbind(old_l, labels)
  }
  con <- gzfile(HITS_F, "w"); utils::write.csv(hits, con, row.names = FALSE); close(con)
  utils::write.csv(labels, LAB_F, row.names = FALSE)
  cat("fixture:", nrow(labels), "units,", nrow(hits), "hit rows\n")
}

#' Capture every on_target unit of a real project as the clean-host corpus
fixture_capture_host <- function(db_path, taxdump_date = NA_character_) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(con))
  u <- DBI::dbGetQuery(con,
    "SELECT s.ID, s.path, s.scaffold, sm.lineage_taxids
       FROM screen s JOIN samples sm ON sm.ID = s.ID
      WHERE s.screen_verdict = 'on_target'
      ORDER BY s.ID, s.path, s.scaffold")
  if (nrow(u) == 0) stop("no on_target units in ", db_path, call. = FALSE)
  hits <- list(); labs <- list()
  for (i in seq_len(nrow(u))) {
    unit <- sprintf("host_%s_%d.%d", u$ID[i], u$path[i], u$scaffold[i])
    h <- DBI::dbGetQuery(con,
      "SELECT hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
         FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?
        ORDER BY hit_rank, hsp",
      params = list(u$ID[i], u$path[i], u$scaffold[i]))
    if (nrow(h) == 0) next   # a unit with no hits scores no_reference, not on_target
    hits[[length(hits) + 1]] <- cbind(unit = unit, h)
    labs[[length(labs) + 1]] <- data.frame(unit = unit, category = "host",
      expected_verdict = "on_target", expected_class = NA_character_,
      source = "real", sample_taxids = u$lineage_taxids[i],
      accession = NA_character_, taxdump_date = taxdump_date,
      note = "clean host scaffold from the reference project",
      stringsAsFactors = FALSE)
  }
  fixture_append(do.call(rbind, hits)[, HIT_COLS], do.call(rbind, labs)[, LAB_COLS])
  invisible(length(labs))
}

a <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% a
a <- a[a != "--force"]
if (length(a) >= 1 && a[1] == "host") {
  if (file.exists(HITS_F) && !force) stop("fixture exists; pass --force", call. = FALSE)
  if (file.exists(HITS_F)) { unlink(HITS_F); unlink(LAB_F) }
  fixture_capture_host(a[2], taxdump_date = if (length(a) >= 3) a[3] else NA_character_)
}
```

- [ ] **Step 2: Run it against the real project**

```bash
Rscript dev/make_screen_fixture.R host \
  "$HOME/Documents/DJMCompBio/RiboPilot_test/full_project_docker/.sqlite" 2026-07-29
```

Expected: `fixture: 19 units, 1833 hit rows`. If either number differs, stop and find out why before continuing; those two numbers are the corpus doc 08 section 3.2.1 cites.

- [ ] **Step 3: Verify the fixture is well formed and small**

```bash
ls -l tests/testthat/fixtures/screen_hits.csv.gz
Rscript -e 'h <- read.csv(gzfile("tests/testthat/fixtures/screen_hits.csv.gz"));
  l <- read.csv("tests/testthat/fixtures/screen_labels.csv");
  cat("units:", length(unique(h$unit)), "labels:", nrow(l), "\n");
  stopifnot(setequal(h$unit, l$unit), !anyNA(l$sample_taxids))'
```

Expected: units 19, labels 19, no error. Fixture under 100 KB.

- [ ] **Step 4: Commit**

```bash
git add dev/make_screen_fixture.R tests/testthat/fixtures/screen_hits.csv.gz tests/testthat/fixtures/screen_labels.csv
git commit -m "test(screen): capture the clean-host corpus as a fixture"
```

---

### Task 2: The loader, the integrity test, and zero false positives

The whole value of Task 1, asserted. Still no network.

**Files:**
- Create: `tests/testthat/test-screen-validation.R`

**Interfaces:**
- Consumes: `screen_channel_taxonomy(hits, sample_taxids, params)` from `R/screen.R`; the two fixture files from Task 1.
- Produces: `fx_units()` -> named list of `list(hits = data.frame, label = one-row data.frame)`, one entry per fixture unit. Task 5 reuses it.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-screen-validation.R`:

```r
# The validation gate of doc 08 section 5, made real. A labelled corpus the
# screen's verdicts are MEASURED against rather than asserted about.
#
# Offline and deterministic by construction: this file reads two committed
# fixtures and nothing else. It never runs blastn, opens a socket, or touches a
# project DB. The capture that produced the fixtures is dev/make_screen_fixture.R.

fx_hits <- read.csv(gzfile(testthat::test_path("fixtures", "screen_hits.csv.gz")),
  stringsAsFactors = FALSE)
fx_labels <- read.csv(testthat::test_path("fixtures", "screen_labels.csv"),
  stringsAsFactors = FALSE)

fx_units <- function() {
  stats::setNames(lapply(fx_labels$unit, function(u) list(
    hits = fx_hits[fx_hits$unit == u, , drop = FALSE],
    label = fx_labels[fx_labels$unit == u, , drop = FALSE])), fx_labels$unit)
}

# score one fixture unit exactly as the pipeline would
fx_score <- function(x) {
  RiboPilot:::screen_channel_taxonomy(x$hits, x$label$sample_taxids)
}

test_that("the fixture is internally consistent", {
  # a label with no hits, or hits with no label, is a corrupt fixture and must
  # fail loudly rather than be silently skipped
  expect_setequal(unique(fx_hits$unit), fx_labels$unit)
  expect_equal(anyDuplicated(fx_labels$unit), 0L)
  expect_false(anyNA(fx_labels$sample_taxids))
  expect_true(all(nzchar(fx_labels$sample_taxids)))
  expect_true(all(fx_labels$expected_verdict %in% RiboPilot:::SCREEN_VERDICTS))
  ec <- fx_labels$expected_class[!is.na(fx_labels$expected_class)]
  expect_true(all(ec %in% RiboPilot:::SCREEN_CLASSES))
  # every unit must carry at least one hit row, or it scores no_reference for a
  # reason that has nothing to do with the screen
  n <- table(fx_hits$unit)
  expect_true(all(n > 0))
})

test_that("the clean-host corpus yields ZERO false positives", {
  # doc 08 section 5: a false positive on an under-represented taxon is the
  # failure that makes the tool useless. This is the assertion the gate is about.
  u <- fx_units()
  host <- u[vapply(u, function(x) identical(x$label$category, "host"), logical(1))]
  expect_gt(length(host), 0)
  verdicts <- vapply(host, function(x) fx_score(x)$verdict, character(1))
  bad <- names(verdicts)[verdicts == "off_target"]
  expect_equal(bad, character(0))
})

test_that("every host unit gets its expected verdict", {
  u <- fx_units()
  host <- u[vapply(u, function(x) identical(x$label$category, "host"), logical(1))]
  for (nm in names(host)) {
    expect_equal(fx_score(host[[nm]])$verdict,
      host[[nm]]$label$expected_verdict, info = nm)
  }
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e 'testthat::test_local(filter = "screen-validation")'`
Expected: PASS. It should pass immediately, because Task 1 captured only units the live screen already called `on_target`, so this is a tautology at first. Task 5 is what makes it bite; the value now is the integrity assertions and the regression floor.

- [ ] **Step 3: Prove it is not vacuous**

The false-positive test must fail if scoring degrades. Mutate `R/screen.R` so
concordance is tested at the shallowest slot only, which doc 11 section 4.5
records as scoring clean host data `off_target`:

```bash
sed -i 's/^  samp_ok <- samp\[-1\]$/  samp_ok <- samp[1]/' R/screen.R
Rscript -e 'testthat::test_local(filter = "screen-validation")'   # must FAIL
git checkout -- R/screen.R
Rscript -e 'testthat::test_local(filter = "screen-validation")'   # must PASS
```

If the mutation leaves the suite green, the fixture is not exercising the scoring
path and this task is not done.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-screen-validation.R
git commit -m "test(screen): assert zero false positives on the host corpus"
```

---

### Task 3: Resolve and capture the spike-in distributions

The one task that needs the network. Run inside `ribopilot:test`, which has `esearch`, `efetch` and `blastn`.

**Files:**
- Modify: `dev/make_screen_fixture.R`

**Interfaces:**
- Consumes: `fixture_append()`, `HIT_COLS`, `LAB_COLS` from Task 1; `parse_blast_candidates(con, id, path, scaffold, hits_file, taxdump_dir)` from `R/blast_genbank.R` (exported).
- Produces: `fixture_capture_spikes(db_path, taxdump_dir, taxdump_date = NA_character_)` -> (invisibly) number of units written.

The spike table. Accessions are resolved at generation time, never hardcoded:

```r
SPIKES <- data.frame(
  key      = c("fungal",  "bacterial", "algal",       "arthropod", "human"),
  query    = c(
    "Fusarium[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Escherichia coli[Organism] AND 16S ribosomal RNA[Title] AND 1200:1700[SLEN]",
    "Arabidopsis thaliana[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Calanus[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]",
    "Homo sapiens[Organism] AND 18S ribosomal RNA[Title] AND 1500:2500[SLEN]"),
  stringsAsFactors = FALSE)
```

- [ ] **Step 1: Add the spike capture to the generator**

Append to `dev/make_screen_fixture.R`:

```r
#' Resolve one accession for a query, newest first
spike_accession <- function(query) {
  out <- suppressWarnings(system2("sh", c("-c", shQuote(sprintf(
    "esearch -db nuccore -query %s | efetch -format acc | head -1",
    shQuote(query)))), stdout = TRUE, stderr = FALSE))
  acc <- trimws(paste(out, collapse = ""))
  if (!nzchar(acc)) stop("no accession for query: ", query, call. = FALSE)
  acc
}

#' The same remote blastn the pipeline runs (inst/nextflow/modules/blast.nf)
spike_blast <- function(fasta, out_file) {
  fmt <- paste("6 saccver stitle pident qcovs evalue staxids bitscore length",
    "qstart qend sstart send")
  ok <- system2("blastn", c("-remote", "-db", "core_nt", "-query", shQuote(fasta),
    "-outfmt", shQuote(fmt), "-max_target_seqs", "50", "-max_hsps", "5",
    "-task", "megablast"), stdout = out_file)
  if (!identical(ok, 0L) || !file.exists(out_file) || file.size(out_file) == 0) {
    stop("blastn returned nothing for ", fasta, call. = FALSE)
  }
  invisible(out_file)
}

#' Capture each spike-in's hit distribution through the pipeline's own parser
#'
#' A scratch DB is used because parse_blast_candidates() writes into blast_hits
#' and needs the real schema. The scratch sample's lineage is irrelevant here:
#' scoring happens in the test, against whatever host lineage the label carries.
fixture_capture_spikes <- function(db_path, taxdump_dir, taxdump_date = NA_character_) {
  hc <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  fish <- DBI::dbGetQuery(hc,
    "SELECT lineage_taxids t FROM samples WHERE lineage_taxids IS NOT NULL LIMIT 1")$t[1]
  DBI::dbDisconnect(hc)
  if (is.na(fish) || !nzchar(fish)) stop("no host lineage in ", db_path, call. = FALSE)

  tmp <- tempfile("rp_fix_"); dir.create(tmp)
  scratch <- file.path(tmp, ".sqlite")
  sc <- DBI::dbConnect(RSQLite::SQLite(), scratch)
  DBI::dbExecute(sc, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(sc, "CREATE TABLE samples(ID TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(sc, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, sequence TEXT)")
  DBI::dbExecute(sc, "INSERT INTO samples VALUES('SPIKE',?)", params = list(fish))

  hits <- list(); labs <- list()
  for (i in seq_len(nrow(SPIKES))) {
    key <- SPIKES$key[i]
    acc <- spike_accession(SPIKES$query[i])
    fa <- file.path(tmp, paste0(key, ".fa"))
    system2("efetch", c("-db", "nuccore", "-id", acc, "-format", "fasta"), stdout = fa)
    if (!file.exists(fa) || file.size(fa) == 0) stop("efetch failed for ", acc, call. = FALSE)
    cat("spike", key, "accession", acc, "\n")
    hf <- file.path(tmp, paste0(key, ".hits"))
    spike_blast(fa, hf)
    DBI::dbExecute(sc, "DELETE FROM blast_hits")
    parse_blast_candidates(sc, "SPIKE", path = 1L, scaffold = i,
      hits_file = hf, taxdump_dir = taxdump_dir)
    h <- DBI::dbGetQuery(sc,
      "SELECT hit_rank, hsp, bitscore, hit_class, lineage_taxids, lineage
         FROM blast_hits WHERE scaffold = ? ORDER BY hit_rank, hsp",
      params = list(i))
    if (nrow(h) == 0) stop("no parsed hits for spike ", key, call. = FALSE)
    cat("  parsed", nrow(h), "hit rows\n")
    unit <- paste0(key, "_on_fish")
    hits[[length(hits) + 1]] <- cbind(unit = unit, h)
    labs[[length(labs) + 1]] <- data.frame(unit = unit, category = key,
      expected_verdict = NA_character_, expected_class = NA_character_,
      source = "spike", sample_taxids = fish, accession = acc,
      taxdump_date = taxdump_date, note = "", stringsAsFactors = FALSE)
  }
  DBI::dbDisconnect(sc)
  fixture_append(do.call(rbind, hits)[, HIT_COLS], do.call(rbind, labs)[, LAB_COLS])
  invisible(nrow(SPIKES))
}
```

and extend the dispatcher:

```r
if (length(a) >= 1 && a[1] == "spikes") {
  fixture_capture_spikes(a[2], a[3], taxdump_date = if (length(a) >= 4) a[4] else NA_character_)
}
```

- [ ] **Step 2: Run it inside the image**

```bash
P="$HOME/Documents/DJMCompBio/RiboPilot_test/full_project_docker"
docker run --rm -v "$HOME/Documents:$HOME/Documents" -w "$(pwd)" \
  -e NCBI_API_KEY ribopilot:test \
  Rscript dev/make_screen_fixture.R spikes "$P/.sqlite" "$P/.taxonomy" 2026-07-29
```

Expected: one `spike <key> accession <ACC>` line and a non-zero parsed-hit count per spike. Remote BLAST is flaky; on failure re-run rather than lowering the bar, and if a query returns nothing, widen that one query and record what changed.

- [ ] **Step 3: Fill in the expected labels from the spec, by hand**

The generator leaves `expected_verdict` and `expected_class` empty on purpose: a
generator that also writes the expectations proves nothing. Set them from doc 13
section 4, and add the second human unit that reuses the human distribution
against an invertebrate host:

```bash
Rscript -e '
l <- read.csv("tests/testthat/fixtures/screen_labels.csv", stringsAsFactors = FALSE)
exp <- list(fungal = c("off_target","fungal"), bacterial = c("off_target","bacterial"),
            algal = c("off_target","algal_plant"), arthropod = c("off_target","metazoan_other"),
            human = c("on_target", NA))
for (k in names(exp)) {
  i <- which(l$category == k)
  l$expected_verdict[i] <- exp[[k]][1]; l$expected_class[i] <- exp[[k]][2]
}
l$note[l$category == "human"] <- paste("EXPECTED MISS: Homo sapiens is Chordata,",
  "same phylum as the fish host, so concordance passes. See doc 13 section 4.1.")
# the same human distribution against a non-chordate host IS catchable
h <- l[l$unit == "human_on_fish", ]
h$unit <- "human_on_invert"; h$expected_verdict <- "off_target"
h$expected_class <- "human"
h$sample_taxids <- Sys.getenv("RP_INVERT_TAXIDS")
h$note <- "same human distribution, octocoral host lineage; see doc 13 section 4.1"
stopifnot(nzchar(h$sample_taxids))
write.csv(rbind(l, h), "tests/testthat/fixtures/screen_labels.csv", row.names = FALSE)
hh <- read.csv(gzfile("tests/testthat/fixtures/screen_hits.csv.gz"), stringsAsFactors = FALSE)
add <- hh[hh$unit == "human_on_fish", ]; add$unit <- "human_on_invert"
con <- gzfile("tests/testthat/fixtures/screen_hits.csv.gz", "w")
write.csv(rbind(hh, add), con, row.names = FALSE); close(con)
cat("labels now:", nrow(rbind(l, h)), "\n")'
```

`RP_INVERT_TAXIDS` is a real rank-aligned octocoral taxid path. Resolve it from
the taxdump already in the test project rather than inventing one:

```bash
export RP_INVERT_TAXIDS=$(docker run --rm \
  -v "$HOME/Documents:$HOME/Documents" ribopilot:test bash -lc \
  "TAXONKIT_DB='$HOME/Documents/DJMCompBio/RiboPilot_test/full_project_docker/.taxonomy' \
   echo 6134 | taxonkit reformat -I 1 -f '{k};{p};{c};{o};{f};{g};{s}' -t -F | cut -f3")
echo "$RP_INVERT_TAXIDS"
```

Verify it is seven semicolon-separated fields and that field 2 is not the fish
phylum before using it.

- [ ] **Step 4: Commit**

```bash
git add dev/make_screen_fixture.R tests/testthat/fixtures/
git commit -m "test(screen): capture cross-kingdom and metazoan spike-ins"
```

---

### Task 4: Assert every spike-in, per unit

**Files:**
- Modify: `tests/testthat/test-screen-validation.R`

**Interfaces:**
- Consumes: `fx_units()`, `fx_score()` from Task 2.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-screen-validation.R`:

```r
test_that("every spike-in unit gets its expected verdict and class", {
  # per unit, never in aggregate: an aggregate floor is what lets one whole
  # category silently stop working
  u <- fx_units()
  spikes <- u[vapply(u, function(x) identical(x$label$source, "spike"), logical(1))]
  expect_gt(length(spikes), 0)
  for (nm in names(spikes)) {
    r <- fx_score(spikes[[nm]])
    lab <- spikes[[nm]]$label
    expect_equal(r$verdict, lab$expected_verdict, info = nm)
    if (!is.na(lab$expected_class)) {
      expect_equal(r$class, lab$expected_class, info = nm)
    }
  }
})

test_that("an uncertain verdict is never counted as catching a contaminant", {
  # counting uncertain as a catch would let the 0.6 gate be tuned to look good
  # by declining to answer
  u <- fx_units()
  catchable <- u[vapply(u, function(x)
    identical(x$label$source, "spike") &&
      identical(x$label$expected_verdict, "off_target"), logical(1))]
  expect_gt(length(catchable), 0)
  for (nm in names(catchable)) {
    expect_equal(fx_score(catchable[[nm]])$verdict, "off_target", info = nm)
  }
})

test_that("the same human distribution is missed on a fish host and caught on an invertebrate", {
  # holds the contaminant constant and varies only the host, so the screen's
  # boundary is measured rather than argued. Homo sapiens is Chordata; so is a
  # ray-finned fish, and concordance is tested at phylum rank or deeper.
  u <- fx_units()
  expect_true(all(c("human_on_fish", "human_on_invert") %in% names(u)))
  same <- u[["human_on_fish"]]; diff <- u[["human_on_invert"]]
  # identical evidence
  expect_equal(nrow(same$hits), nrow(diff$hits))
  expect_equal(same$hits$bitscore, diff$hits$bitscore)
  # different host, different answer
  expect_equal(fx_score(same)$verdict, "on_target")
  expect_equal(fx_score(diff)$verdict, "off_target")
  expect_equal(fx_score(diff)$class, "human")
})

test_that("the confusion matrix is reported", {
  u <- fx_units()
  got <- vapply(u, function(x) fx_score(x)$verdict, character(1))
  want <- vapply(u, function(x) x$label$expected_verdict, character(1))
  cat("\n--- screen validation ---\n")
  print(table(expected = want, got = got))
  n_host <- sum(vapply(u, function(x) identical(x$label$category, "host"), logical(1)))
  fp <- sum(want == "on_target" & got == "off_target")
  catchable <- want == "off_target"
  cat(sprintf("host units %d, false positives %d, recall %d/%d\n",
    n_host, fp, sum(catchable & got == "off_target"), sum(catchable)))
  expect_equal(fp, 0L)
  expect_equal(sum(catchable & got == "off_target"), sum(catchable))
})
```

- [ ] **Step 2: Run it**

Run: `Rscript -e 'testthat::test_local(filter = "screen-validation")'`
Expected: PASS, and the confusion matrix prints. If a spike-in comes back
`uncertain` rather than `off_target`, do NOT relax the assertion. Record the
measured confidence in doc 13 and raise it with the maintainer: a real
contaminant landing under the 0.6 gate is a finding about the gate.

- [ ] **Step 3: Prove the human pair is not vacuous**

```bash
sed -i 's/^  samp_ok <- samp\[-1\]$/  samp_ok <- samp/' R/screen.R
Rscript -e 'testthat::test_local(filter = "screen-validation")'   # must FAIL
git checkout -- R/screen.R
```

Including the shallowest slot makes the invertebrate host concordant with the
human hits, so `human_on_invert` stops being caught. If that leaves the suite
green, the pair is not testing what it claims.

- [ ] **Step 4: Commit**

```bash
git add tests/testthat/test-screen-validation.R
git commit -m "test(screen): assert every spike-in, and the host-lineage boundary"
```

---

### Task 5: Reconcile the design docs

**Files:**
- Modify: `dev/design/08_MULTIUNIT_UX_AND_CONTAMINATION.md`, `dev/design/11_SCREEN_VERDICT_LAYER.md`, `dev/design/12_INTRAGENOMIC_VARIANTS.md`

- [ ] **Step 1: Amend doc 08 section 5's validation gate**

Record that the gate splits: the spike-in half is delivered by doc 13, the clade
half is blocked on data that does not exist locally (three read sets not
downloaded, three clade reference DBs not built, no invert project ever run).
Record that its requested HUMAN spike-in is not a cross-kingdom case on a
chordate host, citing doc 13 section 4.1.

- [ ] **Step 2: Correct doc 11 section 1**

Line 33 reads "It cannot separate two metazoans, because a copepod and a fish
share Metazoa". Replace with: the screen cannot separate two metazoans OF THE
SAME PHYLUM. A copepod is an arthropod and a fish is a chordate, so that pair is
a different-phylum case the screen DOES catch; `tests/testthat/test-screen.R:176`
has pinned it since 3b, and section 4.5 of the same document records the
measurement. Note the sentence predates the phylum-or-deeper rule adopted later
in the same doc and was never reconciled with it.

- [ ] **Step 3: Note the 3c threshold in doc 12**

Record that `min_isomer_anchor_pid` is still unvalidated after doc 13, because
variant comparison runs only on `on_target` units and a spike-in is `off_target`
by construction, so no spike-in is ever compared. It stays blocked on the clade
data.

- [ ] **Step 4: Commit**

```bash
git add dev/design/
git commit -m "docs(design): the gate splits, and the blind spot is same-phylum"
```

---

### Task 6: Whole-branch verification

- [ ] `Rscript -e 'testthat::test_local()'` -> 0 failures, count recorded against the 913 baseline.
- [ ] `Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'` -> 0 errors, 0 warnings, only the two pre-existing NOTEs.
- [ ] `git diff main --name-only -- R/ inst/ NAMESPACE DESCRIPTION` -> **empty**. This plan changes no package code.
- [ ] `grep -rn "^SCREEN_VERSION <- \|^RIBOPILOT_SCHEMA_VERSION <- " R/` -> still `4L` and `15L`.
- [ ] `git diff main -- R/migrate_db.R` -> empty, which is what keeps `TRIAGE_OPTS_SEED` seeding `rule_screen_off_target` as `'flag'`.
- [ ] `grep -rnP '[^\x00-\x7F]' dev/ tests/ | grep -v Binary` -> no hits.
- [ ] Confirm the test opens no network: `grep -nE "blastn|esearch|efetch|url|http" tests/testthat/test-screen-validation.R` -> no hits.
- [ ] Confirm `dev/` still absent from the built package: `Rscript -e 'devtools::build(path = tempdir(), vignettes = FALSE)'` then `tar tzf` the result and grep for `dev/` -> 0 entries.
- [ ] An independent whole-branch review before the PR, treating every claim in it as a hypothesis to check.
