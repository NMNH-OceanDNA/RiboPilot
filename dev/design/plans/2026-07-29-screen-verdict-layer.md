# Screen Verdict Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the machine's per-scaffold contamination opinion a recomputable, versioned home, populated from the taxonomic agreement of a unit's whole BLAST hit distribution with its sample's own lineage.

**Architecture:** A pure scoring core (`R/screen.R`) turns a unit's hit frame plus the sample's taxid path into an LCA, a concordance, a margin, a confidence and a verdict. A DB driver writes one row per unit into a new `screen` table (schema v14). A per-sample Nextflow `screen` process runs it, and `triage` joins on that process so a sample is triaged only once both its curation and its screen exist. A seventh triage rule reads the verdict.

**Tech Stack:** R 4.4, DBI + RSQLite, testthat 3e, Nextflow 25.10+, roxygen2 (Rd only).

**Spec:** `dev/design/11_SCREEN_VERDICT_LAYER.md`

## Global Constraints

- ASCII only in all R code, tests, fixtures, Nextflow and docs. No em dashes, no smart quotes.
- Comments minimal, explaining WHY rather than narrating a change.
- `NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only.
- `RIBOPILOT_SCHEMA_VERSION` moves from `13L` to `14L` in this plan and nowhere else.
- Migrations only ever add tables and columns. Never drop, never retype.
- Verdict strings are exactly `on_target`, `off_target`, `no_reference`, `uncertain`, lower case.
- `screen_class` strings are exactly `bacterial`, `archaeal`, `fungal`, `algal_plant`, `other_euk`, `metazoan_other`, `human`.
- `no_reference` reason codes are exactly `no_hits`, `no_taxdump`, `unresolved_sample_lineage`, `unresolved_hit_lineages`.
- **No identity thresholds anywhere.** The envelope was deferred after measurement; see doc 08 section 3.2.1. A `pident` comparison appearing in this sub-project is a defect.
- Confidence weights are `support` 0.30, `resolution` 0.30, `separation` 0.40, and must sum to 1.
- Run tests from the repo root: `Rscript -e 'testthat::test_local()'`. Baseline is 551 passing, 0 failing.
- Branch is `contamination-phase3b`, already created from `main` at `234e97c`.

---

## File Structure

**Create:**
- `R/screen.R` - the pure scoring core and the DB driver. No Shiny, no Nextflow.
- `inst/nextflow/bin/rp_screen.R` - the per-sample entry point.
- `tests/testthat/test-screen.R` - scoring and verdict tests.
- `tests/testthat/fixtures/screen_hits.csv` - real hits from the reference project, for the regression.

**Modify:**
- `R/migrate_db.R` - `SCREEN_DDL`, the v14 block, the new `triage_opts` column, version bump.
- `R/init_db.R` - the `screen` table on a fresh project.
- `R/triage_rules.R` - the seventh rule.
- `R/triage.R` - `unit_triage_state()` joins `screen`.
- `R/app_utils.R`, `R/app_server.R` - the verdict column in the Assemble evidence table.
- `inst/nextflow/modules/pipeline_r.nf`, `inst/nextflow/main.nf` - the `screen` process and the `triage` join.
- `tests/testthat/test-migrate.R`, `tests/testthat/test-triage-rules.R`, `tests/testthat/test-nextflow.R`.
- `vignettes/Curation-Flags.Rmd`.

---

### Task 1: The pure scoring core

**Files:**
- Create: `R/screen.R`, `tests/testthat/test-screen.R`

**Interfaces:**
- Consumes: `bitscore_weighted_lca()`, `parse_rank_taxids()`, `is_placeholder_taxon()`, `RRNA_LINEAGE_RANKS` from `R/taxonomy_ranks.R`; `%||%`, `%|NA|%` from `R/utils.R`.
- Produces:
  - `SCREEN_VERSION` - `1L`
  - `SCREEN_VERDICTS`, `SCREEN_CLASSES`, `SCREEN_NO_REF_CODES` - the allowed strings
  - `SCREEN_WEIGHTS` - `list(support = 0.30, resolution = 0.30, separation = 0.40)`
  - `screen_rank_depth(rank)` -> numeric in [0, 1]
  - `screen_confidence_terms(lca_n_hits, lca_rank, margin, top_bitscore, unanimous = FALSE)` -> `list(support, resolution, separation, confidence)`
  - `screen_class_for(lca_kingdom, lca_taxon)` -> one of `SCREEN_CLASSES` or `NA_character_`
  - `screen_channel_taxonomy(hits, sample_taxids)` -> `list(verdict, class, confidence, reason, terms, lca_taxid, lca_taxon, lca_rank, lca_n_hits, host_concordant, margin)`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-screen.R`:

```r
# The screen's taxonomy channel. It answers "do this unit's hits agree on a
# taxon, and is that taxon the host's". It needs NO identity threshold, which is
# why it survives the envelope being deferred (see doc 08 section 3.2.1).

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
FUNGI_T <- "4751;4890;147545;5042;1131492;5052;5062"
COPE_T <- "33208;6656;6854;6843;1;2;3"

# one hit row as blast_hits stores it
hit <- function(tax, bs, cls = "rdna", lin = "K;p;c;o;f;g;s") {
  data.frame(hit_class = cls, bitscore = bs, lineage_taxids = tax,
    lineage = lin, stringsAsFactors = FALSE)
}
hits <- function(...) do.call(rbind, list(...))

test_that("the constants are exactly the documented strings", {
  expect_equal(RiboPilot:::SCREEN_VERDICTS,
    c("on_target", "off_target", "no_reference", "uncertain"))
  expect_equal(RiboPilot:::SCREEN_CLASSES,
    c("bacterial", "archaeal", "fungal", "algal_plant", "other_euk",
      "metazoan_other", "human"))
  expect_equal(RiboPilot:::SCREEN_NO_REF_CODES,
    c("no_hits", "no_taxdump", "unresolved_sample_lineage", "unresolved_hit_lineages"))
  expect_equal(sum(unlist(RiboPilot:::SCREEN_WEIGHTS)), 1)
})

test_that("rank depth runs from species down to nothing", {
  expect_equal(RiboPilot:::screen_rank_depth("species"), 1)
  expect_equal(RiboPilot:::screen_rank_depth("kingdom"), 1 / 7, tolerance = 1e-6)
  expect_equal(RiboPilot:::screen_rank_depth("none"), 0)
  expect_equal(RiboPilot:::screen_rank_depth(NA_character_), 0)
  # monotone: a deeper rank is never worth less
  d <- vapply(c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    RiboPilot:::screen_rank_depth, numeric(1))
  expect_false(is.unsorted(d))
})

test_that("each confidence term is computed independently and weighted", {
  tm <- RiboPilot:::screen_confidence_terms(lca_n_hits = 10, lca_rank = "species",
    margin = 500, top_bitscore = 500)
  expect_equal(tm$support, 1)
  expect_equal(tm$resolution, 1)
  expect_equal(tm$separation, 1)
  expect_equal(tm$confidence, 1)

  # support saturates at 10 hits and scales below it
  expect_equal(RiboPilot:::screen_confidence_terms(5, "species", 500, 500)$support, 0.5)
  expect_equal(RiboPilot:::screen_confidence_terms(50, "species", 500, 500)$support, 1)
  # a NA margin contributes nothing rather than erroring
  tm2 <- RiboPilot:::screen_confidence_terms(10, "species", NA_real_, 500)
  expect_equal(tm2$separation, 0)
  expect_equal(tm2$confidence, 0.6)
})

test_that("a concordant unit is on_target", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_equal(r$host_concordant, 1L)
  expect_true(is.na(r$class))
})

test_that("a discordant unit is off_target and names its subclass", {
  h <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "off_target")
  expect_equal(r$host_concordant, 0L)
  expect_equal(r$class, "fungal")
})

test_that("placeholder taxa are excluded from the vote entirely", {
  # uncultured records dominate core_nt 18S and resolve to superkingdom only,
  # which every rank rule would read as maximum taxonomic distance
  h <- hits(hit(FISH_T, 900, lin = "Eukaryota;;;;;;uncultured eukaryote"),
            hit(FISH_T, 890, lin = "Eukaryota;;;;;;environmental sample"),
            hit(FISH_T, 880, lin = "Eukaryota;;;;;;uncultured organism"))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_hit_lineages")
})

test_that("no hits at all is no_reference, not off_target", {
  r <- RiboPilot:::screen_channel_taxonomy(hits()[0, ], FISH_T)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "no_hits")
})

test_that("an unresolvable sample lineage is no_reference with its own code", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, NA_character_)
  expect_equal(r$verdict, "no_reference")
  expect_equal(r$reason_code, "unresolved_sample_lineage")
})

test_that("unanimity is maximal separation, not missing separation", {
  # a unit whose hits ALL agree with the host has no opposing hit to measure a
  # margin against. Scoring that 0 would put the cleanest possible scaffold
  # below the uncertainty gate, which is backwards.
  tm <- RiboPilot:::screen_confidence_terms(3, "species", NA_real_, 900,
    unanimous = TRUE)
  expect_equal(tm$separation, 1)
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "on_target")
  expect_gte(r$confidence, 0.6)
})

test_that("a resolvable but weakly supported unit is uncertain, with evidence kept", {
  # one hit resolving only to kingdom: unanimous, but almost nothing behind it
  h <- hits(hit("4751;;;;;;", 100, lin = "Fungi;;;;;;"))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  expect_equal(r$verdict, "uncertain")
  expect_lt(r$confidence, 0.6)
  # the evidence is still there, which is the point of recording uncertain
  expect_false(is.na(r$lca_taxid))
  expect_equal(r$host_concordant, 0L)
})

test_that("the reason names every term with its value and its weight", {
  h <- hits(hit(FUNGI_T, 900), hit(FUNGI_T, 890), hit(FUNGI_T, 880))
  r <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  for (k in names(RiboPilot:::SCREEN_WEIGHTS)) expect_match(r$reason, k, fixed = TRUE)
  expect_match(r$reason, "0.30", fixed = TRUE)
  expect_match(r$reason, "0.40", fixed = TRUE)
})

test_that("no identity threshold appears anywhere in the scoring core", {
  # the envelope was deferred after measurement (doc 08 section 3.2.1); a pident
  # comparison here would be a defect, so pin its absence. Deparse the functions
  # rather than reading R/screen.R, because under R CMD check the package is
  # installed and the source tree is not present.
  src <- paste(c(deparse(RiboPilot:::screen_channel_taxonomy),
                 deparse(RiboPilot:::screen_confidence_terms)), collapse = " ")
  expect_false(grepl("pident", src, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "screen")'`
Expected: FAIL, `could not find function "screen_channel_taxonomy"`.

- [ ] **Step 3: Implement**

Create `R/screen.R`:

```r
# The screen's taxonomy channel: does this unit's hit distribution agree on a
# taxon, and is that taxon the host's.
#
# It uses NO identity threshold. The identity envelope was deferred after
# measurement against the reference project (doc 08 section 3.2.1): that corpus
# has 2 hits at phylum distance and none beyond, so it cannot calibrate a bound
# acting there, and genuine host 18S reaches 99.6 percent at class distance, so
# no ceiling can sit below that without flagging real host data.
#
# What this cannot do: separate two metazoans. A copepod and a fish both sit
# under Metazoa, so concordance against the sample's path says nothing. Doc 08
# section 3.5 puts those cases in the depth and cross-sample channels.

#' Scoring version. Bump to force a recompute of stored rows.
#' @noRd
SCREEN_VERSION <- 1L

#' @noRd
SCREEN_VERDICTS <- c("on_target", "off_target", "no_reference", "uncertain")

#' @noRd
SCREEN_CLASSES <- c("bacterial", "archaeal", "fungal", "algal_plant",
  "other_euk", "metazoan_other", "human")

#' @noRd
SCREEN_NO_REF_CODES <- c("no_hits", "no_taxdump", "unresolved_sample_lineage",
  "unresolved_hit_lineages")

#' Confidence weights.
#'
#' Priors, not measurements, which is why the screen rule ships as Flag rather
#' than Drop. Any change here must bump [SCREEN_VERSION].
#' @noRd
SCREEN_WEIGHTS <- list(support = 0.30, resolution = 0.30, separation = 0.40)

#' Confidence below which a verdict is recorded as uncertain.
#' @noRd
SCREEN_MIN_CONFIDENCE <- 0.6

#' Depth of a resolved rank, 0 (nothing) to 1 (species)
#' @noRd
screen_rank_depth <- function(rank) {
  if (length(rank) == 0 || is.na(rank[1])) return(0)
  i <- match(rank[1], RRNA_LINEAGE_RANKS)
  if (is.na(i)) 0 else i / length(RRNA_LINEAGE_RANKS)
}

#' The three confidence terms and their weighted mean
#'
#' `unanimous` means every surviving hit fell on one side of the concordance
#' test, so there is no opposing hit to measure a margin against. That is
#' MAXIMAL separation, not absent separation: a unit whose hits all agree with
#' the host is the cleanest case there is, and scoring it 0 would put the
#' cleanest possible scaffold below the uncertainty gate.
#' @noRd
screen_confidence_terms <- function(lca_n_hits, lca_rank, margin, top_bitscore,
                                    unanimous = FALSE) {
  support <- min(1, (lca_n_hits %||% 0) / 10)
  resolution <- screen_rank_depth(lca_rank)
  separation <- if (isTRUE(unanimous)) {
    1
  } else if (is.na(margin %|NA|% NA) || is.na(top_bitscore %|NA|% NA) ||
             !isTRUE(top_bitscore > 0)) {
    0
  } else {
    min(1, abs(margin) / top_bitscore)
  }
  w <- SCREEN_WEIGHTS
  list(support = support, resolution = resolution, separation = separation,
    confidence = unname(support * w$support + resolution * w$resolution +
      separation * w$separation))
}

#' Off-target subclass from the LCA's kingdom-rank name and its taxon label
#' @noRd
screen_class_for <- function(lca_kingdom, lca_taxon) {
  k <- tolower(lca_kingdom %|NA|% "")
  # a human handling contaminant is worth naming even though it is metazoan
  if (grepl("^homo sapiens$", tolower(lca_taxon %|NA|% ""))) return("human")
  if (grepl("bacteria", k)) return("bacterial")
  if (grepl("archaea", k)) return("archaeal")
  if (grepl("fungi", k)) return("fungal")
  if (grepl("viridiplantae|plant", k)) return("algal_plant")
  if (grepl("metazoa", k)) return("metazoan_other")
  if (nzchar(k)) return("other_euk")
  NA_character_
}

#' Score one unit's hit distribution against its sample's lineage
#'
#' @param hits Frame with `hit_class`, `bitscore`, `lineage_taxids`, `lineage`.
#' @param sample_taxids The sample's rank-aligned taxid path.
#' @return list(verdict, class, confidence, reason, reason_code, terms,
#'   lca_taxid, lca_taxon, lca_rank, lca_n_hits, host_concordant, margin).
#' @noRd
screen_channel_taxonomy <- function(hits, sample_taxids) {
  out <- list(verdict = "no_reference", class = NA_character_, confidence = 0,
    reason = NA_character_, reason_code = "no_hits",
    terms = list(support = 0, resolution = 0, separation = 0),
    lca_taxid = NA_integer_, lca_taxon = NA_character_, lca_rank = NA_character_,
    lca_n_hits = 0L, host_concordant = NA_integer_, margin = NA_real_)

  samp <- parse_rank_taxids(sample_taxids %||% NA_character_)
  if (all(is.na(samp))) {
    out$reason_code <- "unresolved_sample_lineage"
    out$reason <- "the sample's own lineage could not be resolved"
    return(out)
  }
  if (is.null(hits) || nrow(hits) == 0) {
    out$reason <- "no GenBank hits are stored for this scaffold"
    return(out)
  }

  # uncultured and environmental records resolve to superkingdom only, which
  # every rank rule reads as maximum distance; drop them from the vote
  keep <- hits$hit_class %in% "rdna" & !is_placeholder_taxon(hits$lineage %|NA|% "")
  h <- hits[keep, , drop = FALSE]
  if (nrow(h) == 0) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "no hit resolved to a usable lineage"
    return(out)
  }

  lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
  if (is.na(lca$taxid)) {
    out$reason_code <- "unresolved_hit_lineages"
    out$reason <- "the hits agree on no shared taxon"
    out$lca_n_hits <- lca$n_hits
    return(out)
  }
  out$lca_taxid <- lca$taxid; out$lca_taxon <- lca$taxon
  out$lca_rank <- lca$rank;  out$lca_n_hits <- lca$n_hits
  out$host_concordant <- as.integer(lca$taxid %in% samp[!is.na(samp)])

  # the margin uses a SECOND, per-hit concordance, distinct from the unit-level
  # one above: it compares the two best-supported opposing hits
  deep <- utils::tail(samp[!is.na(samp)], 1)
  hit_conc <- vapply(h$lineage_taxids, function(x)
    isTRUE(deep %in% parse_rank_taxids(x)), logical(1), USE.NAMES = FALSE)
  best <- function(v) if (any(v)) max(h$bitscore[v], na.rm = TRUE) else NA_real_
  out$margin <- best(hit_conc) - best(!hit_conc)

  # every surviving hit on one side means nothing opposes the winner, which is
  # maximal separation rather than missing separation
  unanimous <- all(hit_conc) || !any(hit_conc)
  tm <- screen_confidence_terms(lca$n_hits, lca$rank, out$margin,
    suppressWarnings(max(h$bitscore, na.rm = TRUE)), unanimous = unanimous)
  out$terms <- tm[c("support", "resolution", "separation")]
  out$confidence <- tm$confidence
  out$reason_code <- NA_character_

  kingdom <- unname(parse_rank_names(h$lineage[1] %|NA|% "")[1])
  out$verdict <- if (tm$confidence < SCREEN_MIN_CONFIDENCE) {
    "uncertain"
  } else if (identical(out$host_concordant, 1L)) {
    "on_target"
  } else {
    "off_target"
  }
  if (identical(out$verdict, "off_target")) {
    out$class <- screen_class_for(kingdom, lca$taxon)
  }
  out$reason <- sprintf(
    "hits agree on %s (%s, %d voting); %s the sample's lineage. support %.2f w %.2f, resolution %.2f w %.2f, separation %.2f w %.2f -> confidence %.2f",
    lca$taxon %|NA|% "an unnamed taxon", lca$rank, lca$n_hits,
    if (identical(out$host_concordant, 1L)) "on" else "NOT on",
    tm$support, SCREEN_WEIGHTS$support, tm$resolution, SCREEN_WEIGHTS$resolution,
    tm$separation, SCREEN_WEIGHTS$separation, tm$confidence)
  out
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "screen")'`
Expected: PASS.

- [ ] **Step 5: Check the whole suite and commit**

Run: `Rscript -e 'testthat::test_local()'`
Expected: 0 failures. Nothing calls the new code yet.

```bash
git add R/screen.R tests/testthat/test-screen.R
git commit -m "feat(screen): taxonomy channel scoring a unit's hit distribution"
```

---

### Task 2: Schema v14

**Files:**
- Modify: `R/migrate_db.R` (`RIBOPILOT_SCHEMA_VERSION` at line 8, a `SCREEN_DDL` constant beside `TRIAGE_OPTS_DDL`, a v14 block after the v13 block), `R/init_db.R`
- Test: `tests/testthat/test-migrate.R`

**Interfaces:**
- Produces: table `screen` per the spec; `samples.reference_support INTEGER`; `triage_opts.rule_screen_off_target TEXT`; constant `SCREEN_DDL`.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-migrate.R`:

```r
test_that("v14 adds the screen table and its companions", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  expect_true("screen" %in% DBI::dbListTables(con))
  for (col in c("screen_verdict", "screen_class", "screen_confidence",
                "screen_reason", "screen_version", "lca_taxid", "lca_taxon",
                "lca_rank", "lca_n_hits", "host_concordant", "margin")) {
    expect_true(col %in% DBI::dbListFields(con, "screen"), info = col)
  }
  expect_true("reference_support" %in% DBI::dbListFields(con, "samples"))
  expect_true("rule_screen_off_target" %in% DBI::dbListFields(con, "triage_opts"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 14L)

  # the new rule ships as flag, per doc 08 section 5's validation gate
  expect_equal(DBI::dbGetQuery(con,
    "SELECT rule_screen_off_target v FROM triage_opts WHERE triage_opts='default'")$v,
    "flag")
})

test_that("the screen table is created empty rather than backfilled", {
  # scoring needs samples.lineage_taxids and a taxdump, neither of which a
  # migration can assume; an empty row must stay distinguishable from a computed one
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen")$n, 0L)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: FAIL, `screen` not in `dbListTables`.

- [ ] **Step 3: Implement**

In `R/migrate_db.R` set line 8 to `RIBOPILOT_SCHEMA_VERSION <- 14L`, add `rule_screen_off_target TEXT,` to `TRIAGE_OPTS_DDL` after `rule_odd_one_out`, add `'flag'` in the matching position of `TRIAGE_OPTS_SEED`, and add beside `TRIAGE_OPTS_DDL`:

```r
#' DDL for the screen table (machine contamination verdict, one row per unit).
#' Shared by new_db() and the v14 migration so both agree on the schema.
#' @noRd
SCREEN_DDL <- "CREATE TABLE IF NOT EXISTS screen (
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
);"
```

Add after the v13 block, before the `user_version` write:

```r
  # --- schema v14: the screen verdict layer ---
  # The machine's opinion about a scaffold lived only as a substring inside
  # annotate_notes, which curate_one() overwrites wholesale and the modal
  # textarea also writes to. Created empty on purpose: scoring needs a resolved
  # sample lineage and a taxdump, which a migration cannot assume.
  if (!has_tbl("screen")) DBI::dbExecute(con, SCREEN_DDL)
  add_col("samples", "reference_support", "INTEGER")
  add_col("triage_opts", "rule_screen_off_target", "TEXT")
  if (ver < 14L && has_col("triage_opts", "rule_screen_off_target")) {
    DBI::dbExecute(con,
      "UPDATE triage_opts SET rule_screen_off_target = 'flag'
       WHERE rule_screen_off_target IS NULL")
  }
```

In `R/init_db.R`, add `DBI::dbExecute(con, SCREEN_DDL)` beside the other shared-DDL calls, and `reference_support INTEGER,` to the `samples` column list.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: PASS, including the pre-existing idempotency test.

- [ ] **Step 5: Commit**

```bash
git add R/migrate_db.R R/init_db.R tests/testthat/test-migrate.R
git commit -m "feat(db): screen table, reference_support, screen rule option (v14)"
```

---

### Task 3: The database driver

**Files:**
- Modify: `R/screen.R`
- Test: `tests/testthat/test-screen.R`

**Interfaces:**
- Consumes: `screen_channel_taxonomy()`, `SCREEN_VERSION` from Task 1; the `screen` table from Task 2.
- Produces:
  - `screen_one(con, id, path, scaffold, force = FALSE)` -> invisibly the verdict string
  - `screen_sample(con, id, force = FALSE)` -> invisibly the number of units written
  - `update_reference_support(con, id)` -> invisibly the value written

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-screen.R`:

```r
mk_screen_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::SCREEN_DDL)
  DBI::dbExecute(con, RiboPilot:::BLAST_HITS_DDL)
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage_taxids TEXT, reference_support INT)")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT)")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1',?,NULL)", params = list(FISH_T))
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,1)")
  DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,2)")
  con
}
add_hits <- function(con, scaf, tax, n = 3, bs = 900, lin = "Eukaryota;Chordata;Actinopteri;o;f;g;s") {
  for (i in seq_len(n)) {
    DBI::dbExecute(con,
      "INSERT INTO blast_hits (ID,path,scaffold,hit_rank,hsp,bitscore,lineage_taxids,lineage,hit_class)
       VALUES ('S1',1,?,?,1,?,?,?,'rdna')",
      params = list(scaf, i, bs - i, tax, lin))
  }
}

test_that("screen_one writes one row carrying the verdict and its evidence", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE scaffold=1")
  expect_equal(nrow(r), 1L)
  expect_equal(r$screen_verdict, "on_target")
  expect_equal(r$host_concordant, 1L)
  expect_equal(r$screen_version, RiboPilot:::SCREEN_VERSION)
  expect_true(!is.na(r$screen_reason) && nzchar(r$screen_reason))
})

test_that("a unit with no hits gets a no_reference row rather than none", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  screen_one(con, "S1", 1, 2)
  r <- DBI::dbGetQuery(con, "SELECT screen_verdict FROM screen WHERE scaffold=2")
  expect_equal(r$screen_verdict, "no_reference")
})

test_that("a stored row at the current version is left alone unless forced", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "UPDATE screen SET screen_reason='SENTINEL' WHERE scaffold=1")
  screen_one(con, "S1", 1, 1)
  expect_equal(DBI::dbGetQuery(con, "SELECT screen_reason r FROM screen WHERE scaffold=1")$r,
    "SENTINEL")
  screen_one(con, "S1", 1, 1, force = TRUE)
  expect_false(identical(
    DBI::dbGetQuery(con, "SELECT screen_reason r FROM screen WHERE scaffold=1")$r, "SENTINEL"))
})

test_that("a row stored at a lower version is recomputed without force", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T)
  screen_one(con, "S1", 1, 1)
  DBI::dbExecute(con, "UPDATE screen SET screen_version=0, screen_reason='OLD' WHERE scaffold=1")
  screen_one(con, "S1", 1, 1)
  r <- DBI::dbGetQuery(con, "SELECT screen_version v, screen_reason s FROM screen WHERE scaffold=1")
  expect_equal(r$v, RiboPilot:::SCREEN_VERSION)
  expect_false(identical(r$s, "OLD"))
})

test_that("screen_sample covers every unit and sets reference support", {
  con <- mk_screen_db(); on.exit(DBI::dbDisconnect(con))
  add_hits(con, 1, FISH_T); add_hits(con, 2, FISH_T)
  n <- screen_sample(con, "S1")
  expect_equal(n, 2L)
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen")$n, 2L)
  rs <- DBI::dbGetQuery(con, "SELECT reference_support v FROM samples WHERE ID='S1'")$v
  expect_false(is.na(rs))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "screen")'`
Expected: FAIL, `could not find function "screen_one"`.

- [ ] **Step 3: Implement**

Append to `R/screen.R`:

```r
#' Score and store one unit's screen verdict
#'
#' Skips a row already stored at the current [SCREEN_VERSION] unless `force`, so
#' re-running the pipeline is cheap and a scoring change forces a recompute
#' without a schema migration.
#'
#' @param con Open project DB connection.
#' @param id,path,scaffold Assembly unit.
#' @param force Recompute even when the stored row is current.
#' @return (invisibly) the verdict written, or the stored one when skipped.
#' @export
screen_one <- function(con, id, path = 1L, scaffold = 1L, force = FALSE) {
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  cur <- DBI::dbGetQuery(con,
    "SELECT screen_verdict, screen_version FROM screen
     WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
  if (!isTRUE(force) && nrow(cur) == 1 &&
      isTRUE((cur$screen_version[1] %|NA|% -1L) >= SCREEN_VERSION)) {
    return(invisible(cur$screen_verdict[1]))
  }
  hits <- DBI::dbGetQuery(con,
    "SELECT hit_class, bitscore, lineage_taxids, lineage FROM blast_hits
     WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
  samp <- tryCatch(DBI::dbGetQuery(con,
    "SELECT lineage_taxids t FROM samples WHERE ID = ?", params = list(id))$t[1],
    error = function(e) NA_character_)
  r <- screen_channel_taxonomy(hits, samp)
  DBI::dbExecute(con,
    "DELETE FROM screen WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))
  DBI::dbAppendTable(con, "screen", data.frame(
    ID = id, path = path, scaffold = scaffold,
    screen_verdict = r$verdict, screen_class = r$class,
    screen_confidence = r$confidence, screen_reason = r$reason,
    screen_version = SCREEN_VERSION, lca_taxid = r$lca_taxid,
    lca_taxon = r$lca_taxon, lca_rank = r$lca_rank, lca_n_hits = r$lca_n_hits,
    host_concordant = r$host_concordant, margin = r$margin,
    time_stamp = as.integer(Sys.time()), stringsAsFactors = FALSE))
  invisible(r$verdict)
}

#' Screen every unit of one sample and refresh its reference support
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @param force Recompute rows already at the current version.
#' @return (invisibly) the number of units screened.
#' @export
screen_sample <- function(con, id, force = FALSE) {
  u <- DBI::dbGetQuery(con,
    "SELECT path, scaffold FROM assemblies WHERE ID = ? ORDER BY path, scaffold",
    params = list(id))
  for (i in seq_len(nrow(u))) screen_one(con, id, u$path[i], u$scaffold[i], force = force)
  update_reference_support(con, id)
  invisible(nrow(u))
}

#' Refresh a sample's reference support
#'
#' The minimum rank distance across the sample's ON-TARGET units only. Computed
#' from on-target units because a sample whose only resolvable hits belong to a
#' contaminant would otherwise report excellent support for the wrong organism.
#' @noRd
update_reference_support <- function(con, id) {
  v <- DBI::dbGetQuery(con,
    "SELECT MIN(h.rank_distance) v FROM blast_hits h
       JOIN screen s ON s.ID = h.ID AND s.path = h.path AND s.scaffold = h.scaffold
     WHERE h.ID = ? AND s.screen_verdict = 'on_target' AND h.rank_distance IS NOT NULL",
    params = list(id))$v[1]
  DBI::dbExecute(con, "UPDATE samples SET reference_support = ? WHERE ID = ?",
    params = list(v, id))
  invisible(v)
}
```

Add `export(screen_one)` and `export(screen_sample)` to `NAMESPACE` BY HAND, in alphabetical position. Do not run `devtools::document()`.

- [ ] **Step 4: Run tests, regenerate Rd, commit**

```bash
Rscript -e 'testthat::test_local()'
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
git add R/screen.R NAMESPACE man/ tests/testthat/test-screen.R
git commit -m "feat(screen): store one verdict per unit, versioned and recomputable"
```

---

### Task 4: The pipeline process

**Files:**
- Create: `inst/nextflow/bin/rp_screen.R`
- Modify: `inst/nextflow/modules/pipeline_r.nf`, `inst/nextflow/main.nf`
- Test: `tests/testthat/test-nextflow.R`

**Interfaces:**
- Consumes: `screen_sample()` from Task 3.

- [ ] **Step 1: Write the driver**

Create `inst/nextflow/bin/rp_screen.R`, mirroring `rp_triage.R`'s style exactly (read it first):

```r
#!/usr/bin/env Rscript
# Score the contamination screen for every unit of one sample, from the BLAST
# hit distribution already stored by blast_genbank.
# Args: <db> <id>
a <- commandArgs(trailingOnly = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "rp_env.R"))
con <- open_project_db(a[1]); on.exit(DBI::dbDisconnect(con))
n <- screen_sample(con, a[2])
cat("screen: ", a[2], " - ", n, " unit(s)\n", sep = "")
```

Make it executable: `chmod +x inst/nextflow/bin/rp_screen.R`.

- [ ] **Step 2: Add the process**

In `inst/nextflow/modules/pipeline_r.nf`, beside the `triage` process:

```groovy
// Machine contamination verdict per unit, from the retained BLAST hit
// distribution. Best-effort: a taxonomy failure must not sink a sample.
process screen {
    cache false
    tag "${id}"
    errorStrategy 'ignore'
    input:
    val id
    output:
    val id
    script:
    """
    RIBOPILOT_SRC='${params.ribopilot_src}' rp_screen.R '${params.db}' '${id}'
    """
}
```

- [ ] **Step 3: Wire it in**

In `inst/nextflow/main.nf`, add `screen` to the `pipeline_r.nf` include list, then replace the triage block:

```groovy
    // The screen reads only blast_hits, so it hangs off blast_genbank rather
    // than the curate chain. triage joins on it, so a sample is triaged only
    // once both its curation and its screen exist.
    screen(blast_genbank.out)
    triage(curate.out.map { id, upath, uscaf -> tuple(id, "${upath}.${uscaf}") }
        .groupTuple()
        .map { id, _units -> tuple(id, true) }
        .join(screen.out.map { id -> tuple(id, true) })
        .map { id, _a, _b -> id })
```

The explicit `tuple(id, true)` on both sides is deliberate: `join` keys on the first element, and joining two channels of bare values is easy to get subtly wrong. Verify with `nextflow lint inst/nextflow/main.nf` and report the exact warning count, which was 6 before this change.

- [ ] **Step 4: Test the wiring**

Read `tests/testthat/test-nextflow.R` first and follow its existing style. Add a test asserting `main.nf` includes `screen`, that the `screen` process exists in `pipeline_r.nf` with `errorStrategy 'ignore'`, and that `triage`'s input joins on `screen.out`.

- [ ] **Step 5: Run and commit**

```bash
Rscript -e 'testthat::test_local()'
nextflow lint inst/nextflow/main.nf
git add inst/nextflow/ tests/testthat/test-nextflow.R
git commit -m "feat(pipeline): run the screen per sample and gate triage on it"
```

---

### Task 5: Triage rule, app column, docs

**Files:**
- Modify: `R/triage_rules.R`, `R/triage.R`, `R/app_utils.R`, `R/app_server.R`, `vignettes/Curation-Flags.Rmd`
- Test: `tests/testthat/test-triage-rules.R`, `tests/testthat/test-triage.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-triage-rules.R`:

```r
test_that("the screen rule fires only on an off_target verdict", {
  d <- RiboPilot:::triage_rules(us(list(feats = 5, nuc = 1, screen = "off_target")))
  expect_equal(d$new_rule, "screen_off_target")
  expect_true(is.na(d$new_triage))          # ships as flag, not drop
  expect_match(d$new_reason, "does not")

  for (v in c("on_target", "no_reference", "uncertain", NA_character_)) {
    d2 <- RiboPilot:::triage_rules(us(list(feats = 5, nuc = 1, screen = v)))
    expect_false(identical(d2$new_rule, "screen_off_target"), info = v %|NA|% "NA")
  }
})

test_that("the screen rule sits after the combination rules and before odd_one_out", {
  k <- RiboPilot:::triage_rule_keys()
  expect_equal(k[length(k) - 1], "screen_off_target")
  expect_equal(k[length(k)], "odd_one_out")
})

test_that("the screen rule may be set to drop but ships as flag", {
  def <- Filter(function(r) r$key == "screen_off_target", RiboPilot:::TRIAGE_RULE_DEFS)[[1]]
  expect_equal(def$default, "flag")
  expect_true("drop" %in% def$actions)
})
```

The `us()` helper must gain a `screen` field defaulting to `NA_character_`, mapped to a `screen_verdict` column. Update it in place; the existing calls pass no `screen` and must keep working.

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage-rules")'`
Expected: FAIL, `screen_off_target` not in the registry.

- [ ] **Step 3: Implement**

In `R/triage_rules.R`, insert into `TRIAGE_RULE_DEFS` between `feat_non_rdna_hit` and `odd_one_out`:

```r
  list(
    key = "screen_off_target",
    label = "GenBank taxonomy disagrees with the sample",
    meaning = "The scaffold's BLAST hits agree on a taxon that is not on this sample's lineage.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!identical(ctx$screen[i], "off_target")) return(NA_character_)
      sprintf("the scaffold's GenBank hits agree on %s, which does not sit on this sample's lineage",
        ctx$screen_taxon[i] %|NA|% "another taxon")
    }
  ),
```

Add to `triage_rule_ctx()`:

```r
    screen = d$screen_verdict %||% rep(NA_character_, n),
    screen_taxon = d$lca_taxon %||% rep(NA_character_, n),
```

In `R/triage.R`, add to `unit_triage_state()`'s SELECT and its LEFT JOINs:

```sql
       sc.screen_verdict, sc.screen_class, sc.screen_confidence, sc.lca_taxon,
```
```sql
     LEFT JOIN screen sc ON sc.ID = a.ID AND sc.path = a.path AND sc.scaffold = a.scaffold
```

- [ ] **Step 4: Surface it in the app**

In `R/app_utils.R`'s `unit_evidence()`, add the same LEFT JOIN and select `sc.screen_verdict`, `sc.screen_class`, `sc.screen_confidence`, `sc.lca_taxon`. In `R/app_server.R`, add a `screen_verdict` column to the scaffold table's `lead` vector after `elements`, with a colDef named "Screen" rendering `on_target` green, `off_target` red, `uncertain` amber and `no_reference` neutral grey. Reuse `status_badge()` or the existing badge helpers rather than writing a new renderer; read `rt_shares_badge()` for the house style.

- [ ] **Step 5: Vignette**

Add a short paragraph to `vignettes/Curation-Flags.Rmd` in the triage section, in the file's existing voice and flat `#` heading level, saying what the screen verdict means, that `no_reference` is a normal state for an offline project or a sparse clade rather than a problem, and that this check compares against the sample's own lineage so it catches a fungus or a bacterium but cannot separate two metazoans.

- [ ] **Step 6: Verify and commit**

```bash
Rscript -e 'testthat::test_local()'
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'
```

Expected: 0 test failures; check at 0 errors and 0 warnings. Two NOTEs are pre-existing (crayon/dbplyr imports, `new_db.Rd` line width); report any NOTE outside that list.

```bash
git add R/ tests/ vignettes/ man/
git commit -m "feat(triage): flag a scaffold whose GenBank taxonomy disagrees with its sample"
```

---

## Self-Review

**Spec coverage:** section 3 schema -> Task 2; section 4 channel -> Task 1; section 5 confidence -> Task 1; section 6 verdicts -> Task 1; section 7 pipeline -> Task 4; section 8 triage -> Task 5; section 9 testing -> distributed, with the reference-project regression fixture noted below; section 10 out of scope -> nothing here touches the envelope, the tool-kingdom channel, depth, cross-sample, or `annotate_notes`.

**Gap I am recording rather than hiding:** the spec's section 9 asks for a regression fixture extracted from the reference project, asserting real units come back `on_target`. No task builds it, because the extraction needs the reference database, which is not in the repo. Task 3's `mk_screen_db()` covers the same paths synthetically. Whoever executes this should either add the fixture from a copy of the reference project or record that it was skipped.

**Placeholder scan:** none. Every code step carries its code. Task 4 step 4 and Task 5 step 4 describe test and UI work against files whose current contents must be read first, and say so explicitly.

**Type consistency:** `screen_channel_taxonomy()` returns the same named list everywhere; `screen_one()` returns a verdict string; `screen_sample()` returns an integer count. `host_concordant` is `1L`/`0L`/`NA_integer_` in R and INTEGER in SQL. Verdict and class strings match the Global Constraints list exactly.

**Helper behaviour verified while writing this plan, so no task rests on a guess:**
- `parse_rank_taxids(NA_character_)` returns a length-7 all-NA integer vector rather than erroring, which is what the `unresolved_sample_lineage` path relies on.
- `is_placeholder_taxon()` is vectorised and returns TRUE for both an empty string and an "uncultured" prefix.
- `RRNA_LINEAGE_RANKS[1]` is `"kingdom"`, so `parse_rank_names(lineage)[1]` is the kingdom name `screen_class_for()` needs.
- `bitscore_weighted_lca()` returns `list(taxid, taxon, rank, n_hits)` and is already exported.
