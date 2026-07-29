# Intragenomic Variants Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Separate two situations that look identical to a curator today: assembly-graph isomers that differ only in how many IGS subrepeats were traversed, and genuine intragenomic paralogous rDNA copies. The first is a decision that does not matter; the second determines what GenBank receives.

**Architecture:** A pure grouping core turns a set of per-unit anchor strings plus their total scaffold lengths into a class, a group label and an identity per unit. A per-sample DB pass extracts each `on_target` unit's anchor string from `annotations` + `assemblies.sequence`, runs the core, and UPDATEs three new columns on the existing `screen` table (schema v15). It runs inside `screen_sample()` after the per-unit loop, where `update_reference_support()` already sits, because the finding is a property of a sample's whole unit set and cannot be computed one unit at a time.

**Tech Stack:** R 4.4, DBI + RSQLite, Biostrings + pwalign, testthat 3e, roxygen2 (Rd only).

**Spec:** `dev/design/12_INTRAGENOMIC_VARIANTS.md`

---

## Decisions taken before writing this plan

Three of the spec's statements did not survive contact with the code. Each was resolved with the maintainer and each needs a matching amendment to the design docs (Task 9).

**1. `rrna_params()` has no `screen` sublist.** Doc 12 section 4 places `min_isomer_anchor_pid` in it; doc 11 section 5 placed the confidence weights there. Neither was ever built: 3b shipped `SCREEN_WEIGHTS` and `SCREEN_MIN_CONFIDENCE` as top-level constants at `R/screen.R:34` and `R/screen.R:38`. **Resolution: build the sublist** (Task 1) and move the existing constants into it, values unchanged, so both docs become true and a clade can retune the screen as they always claimed.

**2. Doc 12 section 8 is deferred entirely.** Both salvage items are larger than the spec implies:
  - Extending `parse_barrnap()`'s name map retains nothing. `reconcile_features()` selects elements by literal name (`R/annotate.R:110-117`), so a 16S/23S/5S row is discarded one function later. Retaining it needs a store, which is a schema decision touching the app element strip, export markers and curation length rules.
  - The ITSx taxonomic profile is not reachable from R. Both Nextflow processes declare only `path("*.itsx.positions.txt")` as output (`inst/nextflow/modules/annotate.nf:12`, `inst/nextflow/modules/pipeline_r.nf:149`), so every other ITSx file dies with the work directory. Harvesting it needs an output declaration, a new argument through `rp_reconcile.R`, and somewhere to put it.

  **Resolution: neither ships here.** Doc 12 section 8 is rewritten to record both measurements so the next worker does not rediscover them.

**3. `min_isomer_anchor_pid` is `99.5`.** The spec names the constant but never values it. `pwalign_pid()` returns PID1 on a 0-100 scale (`R/blast_ref_align.R:167`), matching `min_anchor_pid`'s existing percent convention. 99.5 across roughly 4 kb of concatenated 18S+5.8S+28S tolerates a handful of assembly errors without admitting a real paralog.

Two further clarifications this plan makes, which the spec left open:

**4. Classification is on identity alone; total length is evidence, not a test.** Doc 12's table gives `isomer` the condition "identity >= threshold AND total length differs", and gives no row at all for "identity >= threshold AND total length identical". Since the only two classes are `isomer` and `paralog`, and identity alone already decides between them, length cannot change the answer. Two units with identical anchors and identical total length are the same molecule assembled twice, for which "pick either" is exactly the right message. The length difference is recorded in the reason string as evidence.

**5. `variant_class` is a per-sample property once groups exist.** The spec's table is written per pair. Reduced to one row per unit: a sample with two or more groups gets `paralog` on **every** comparable unit, including units that are mutually isomeric. Telling the curator of unit A "pick either" would be false when the sample also holds a distinct copy C, and the paralog message ("which you lock determines what GenBank receives") is true of A, B and C alike. `variant_group` is what tells the curator that A and B are interchangeable.

---

## Global Constraints

- ASCII only in all R code, tests, fixtures and docs. No em dashes, no smart quotes.
- Comments minimal, explaining WHY rather than narrating a change.
- `NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd with `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'` only. **This plan adds no exports**, so `NAMESPACE` must end unchanged.
- `RIBOPILOT_SCHEMA_VERSION` moves from `14L` to `15L` in this plan and nowhere else.
- `SCREEN_VERSION` moves from `2L` to `3L` in this plan and nowhere else. Both the new columns and the weights relocating count as a scoring change under doc 11 section 5's rule.
- Migrations only ever add tables and columns. Never drop, never retype.
- `variant_class` strings are exactly `isomer`, `paralog`, `singleton`, lower case, or SQL NULL.
- `%|NA|%` (`R/utils.R`) is scalar only and errors on a vector. `%||%` tests length, not NA. This has caused defects in three separate tasks.
- **Nothing here gates or flags.** No triage rule is added, no `submission_risk` string is added, no export is blocked. A paralog is real biology; colouring it like a defect is what this sub-project exists to stop.
- The finding is orthogonal to `screen_verdict`. A paralogous variant is emphatically `on_target`, and a test must pin that.
- Run tests from the repo root: `Rscript -e 'testthat::test_local()'`. Baseline is **705 passing, 0 failing**.
- Branch is `contamination-phase3c`, already created from `main` at `001fb3f`.
- After each task: revert the fix, confirm its test fails, restore, confirm the tree is clean. A test that passes with its target broken pins nothing. Three tests on the 3b branch looked correct and pinned nothing.

---

## File Structure

**Create:**
- `R/screen_variants.R` - anchor extraction, the pure grouping core, and the per-sample DB pass.
- `tests/testthat/test-screen-variants.R` - all of the above.

**Modify:**
- `R/params_rrna.R` - the `screen` sublist on every ruleset.
- `R/screen.R` - `SCREEN_VERSION` to `3L`; weights read from params; call the variant pass from `screen_sample()`.
- `R/curate_rrna.R` - extract the ruleset resolver so the screen can reuse it.
- `R/migrate_db.R` - `SCREEN_DDL` gains three columns, a v15 block, version bump.
- `R/utils_reactable.R`, `R/app_utils.R`, `R/app_server.R` - the variant column in the Assemble scaffold table.
- `R/export.R` - the variant fact in the per-unit summary.
- `tests/testthat/test-migrate.R`, `test-params.R`, `test-screen.R`, `test-curate.R`, `test-export-units.R`.
- `dev/design/12_INTRAGENOMIC_VARIANTS.md`, `dev/design/11_SCREEN_VERDICT_LAYER.md`, `vignettes/Curation-Flags.Rmd`.

**Not modified, deliberately:** `R/triage_rules.R`, `R/triage.R`, `NAMESPACE`, `inst/nextflow/**`. The screen process already runs per sample and already calls `screen_sample()`, so the pipeline needs no change at all.

---

### Task 1: The `screen` params sublist and a shared ruleset resolver

Makes doc 11's and doc 12's claim true, and gives `min_isomer_anchor_pid` the home the spec names.

**Files:**
- Modify: `R/params_rrna.R`, `R/screen.R`, `R/curate_rrna.R`
- Test: `tests/testthat/test-params.R`, `tests/testthat/test-screen.R`, `tests/testthat/test-curate.R`

**Interfaces:**
- Produces: `rrna_params(ruleset)$screen` = `list(support = 0.30, resolution = 0.30, separation = 0.40, min_confidence = 0.6, min_isomer_anchor_pid = 99.5)`.
- Produces: `resolve_ruleset(con, id, path, scaffold, ruleset = NULL)` -> `list(ruleset, note)`, the logic currently inline at `R/curate_rrna.R:236-256`.
- Produces: `screen_params(con, id)` -> the `screen` sublist for a sample, resolved once per sample rather than per unit.
- `SCREEN_WEIGHTS` and `SCREEN_MIN_CONFIDENCE` remain defined in `R/screen.R` as the default sublist's source of truth, so `screen_channel_taxonomy()` keeps working when called with no params (which is how every pure test calls it).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-params.R`:

```r
test_that("every ruleset carries a screen sublist and its weights still sum to 1", {
  for (rs in RiboPilot:::RRNA_RULESETS) {
    p <- RiboPilot:::rrna_params(rs)
    expect_true(is.list(p$screen), info = rs)
    for (k in c("support", "resolution", "separation", "min_confidence",
                "min_isomer_anchor_pid")) {
      expect_true(!is.null(p$screen[[k]]), info = paste(rs, k))
    }
    expect_equal(sum(unlist(p$screen[c("support", "resolution", "separation")])), 1,
      info = rs)
  }
})

test_that("relocating the screen tunables changed none of their values", {
  # doc 11 section 5: the weights are priors. Building the sublist must not
  # silently retune them, or every stored verdict shifts for a reason nobody chose.
  p <- RiboPilot:::rrna_params("metazoa_generic")$screen
  expect_equal(p$support, 0.30)
  expect_equal(p$resolution, 0.30)
  expect_equal(p$separation, 0.40)
  expect_equal(p$min_confidence, 0.6)
  expect_equal(p$min_isomer_anchor_pid, 99.5)
})

test_that("a clade ruleset can retune the screen without touching the default", {
  base <- RiboPilot:::rrna_params("metazoa_generic")$screen
  for (rs in RiboPilot:::RRNA_RULESETS) {
    expect_equal(RiboPilot:::rrna_params(rs)$screen, base, info = rs)
  }
  # the mechanism exists even though no clade uses it yet: modifyList on the
  # sublist must not drop the untouched keys
  tuned <- utils::modifyList(base, list(min_isomer_anchor_pid = 99.9))
  expect_equal(tuned$support, 0.30)
  expect_equal(tuned$min_isomer_anchor_pid, 99.9)
})
```

Append to `tests/testthat/test-curate.R`:

```r
test_that("the extracted ruleset resolver keeps curate_one's precedence", {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, lineage TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, curate_opts TEXT)")
  DBI::dbExecute(con, "CREATE TABLE curate_opts(curate_opts TEXT, ruleset TEXT)")
  DBI::dbExecute(con, "INSERT INTO curate_opts VALUES('default','metazoa_generic')")
  DBI::dbExecute(con, "INSERT INTO curate_opts VALUES('fish','actinopterygii')")
  DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,1,'default')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Eukaryota; Metazoa; Cnidaria; Scyphozoa')")

  # explicit argument wins over everything
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1, "asteroidea")$ruleset,
    "asteroidea")
  # a generic curate_opts falls through to the lineage
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1)$ruleset, "scyphozoa")
  # a non-generic curate_opts is honoured as-is
  DBI::dbExecute(con, "UPDATE annotate SET curate_opts='fish'")
  expect_equal(RiboPilot:::resolve_ruleset(con, "S1", 1, 1)$ruleset, "actinopterygii")
  # a lineage matching no clade records the fallback rather than hiding it
  DBI::dbExecute(con, "UPDATE annotate SET curate_opts='default'")
  DBI::dbExecute(con, "UPDATE samples SET lineage='Eukaryota; Metazoa; Porifera'")
  r <- RiboPilot:::resolve_ruleset(con, "S1", 1, 1)
  expect_equal(r$ruleset, "metazoa_generic")
  expect_equal(r$note, "no_clade_ruleset")
})
```

Append to `tests/testthat/test-screen.R`:

```r
test_that("the scoring core reads its weights from the params sublist", {
  h <- hits(hit(FISH_T, 900), hit(FISH_T, 890), hit(FISH_T, 880))
  base <- RiboPilot:::screen_channel_taxonomy(h, FISH_T)
  # a ruleset that weighted support alone would score this differently; if the
  # core ignored the argument, both calls would agree
  tuned <- RiboPilot:::screen_channel_taxonomy(h, FISH_T,
    params = list(support = 1, resolution = 0, separation = 0, min_confidence = 0.6))
  expect_false(isTRUE(all.equal(base$confidence, tuned$confidence)))
  expect_equal(tuned$confidence, 0.3)  # lca_n_hits 3 -> support 0.3, weight 1
})

test_that("SCREEN_VERSION is 3 so v2 rows recompute", {
  expect_equal(RiboPilot:::SCREEN_VERSION, 3L)
})
```

- [ ] **Step 2: Run tests to verify they fail**

`Rscript -e 'testthat::test_local(filter = "params|curate|screen")'`
Expected: FAIL on `p$screen` being NULL, on `resolve_ruleset` not existing, and on `SCREEN_VERSION` being `2L`.

- [ ] **Step 3: Implement**

In `R/params_rrna.R`, add to `rrna_params()`'s `base` list, after `check_contaminant`:

```r
    # Screen tunables. Doc 11 section 5 and doc 12 section 4 both site them here
    # so a clade can retune them; until this sub-project they lived as top-level
    # constants in R/screen.R and no clade could reach them.
    screen = SCREEN_PARAMS_DEFAULT
```

and define beside `RRNA_ANCHOR_PID_INVERT`:

```r
#' Screen tunables, shared by every ruleset until a clade needs its own.
#'
#' The three weights are priors, not measurements, which is why the screen rule
#' ships as Flag; `min_isomer_anchor_pid` is a prior too, but a safer one,
#' because it compares two sequences from ONE organism rather than across taxa.
#' Any change to any of them must bump [SCREEN_VERSION].
#' @noRd
SCREEN_PARAMS_DEFAULT <- list(
  support = 0.30, resolution = 0.30, separation = 0.40,
  min_confidence = 0.6,
  min_isomer_anchor_pid = 99.5
)
```

In `R/screen.R`: set `SCREEN_VERSION <- 3L`. Keep `SCREEN_WEIGHTS` and `SCREEN_MIN_CONFIDENCE` but define them from `SCREEN_PARAMS_DEFAULT` so there is one source of truth. Give `screen_confidence_terms()` and `screen_channel_taxonomy()` a `params = SCREEN_PARAMS_DEFAULT` argument and read `params$support` etc. rather than `SCREEN_WEIGHTS$...`. Add:

```r
#' The screen tunables for a sample, resolved through its curation ruleset
#' @noRd
screen_params <- function(con, id) {
  rs <- tryCatch(resolve_ruleset(con, id, 1L, 1L)$ruleset, error = function(e) NULL)
  tryCatch(rrna_params(rs)$screen, error = function(e) SCREEN_PARAMS_DEFAULT) %||%
    SCREEN_PARAMS_DEFAULT
}
```

In `R/curate_rrna.R`, lift lines 236-256 into `resolve_ruleset(con, id, path, scaffold, ruleset = NULL)` returning `list(ruleset =, note =)` and call it from `curate_one()`. The body must be moved verbatim, not rewritten.

Have `screen_one()` pass `screen_params(con, id)` into `screen_channel_taxonomy()`.

- [ ] **Step 4: Verify**

Run the three filters, then the full suite. Then mutate: change `SCREEN_PARAMS_DEFAULT$support` to `0.31` and confirm the sum-to-1 test fails; restore.

- [ ] **Step 5: Commit**

```bash
git add R/params_rrna.R R/screen.R R/curate_rrna.R tests/testthat/test-params.R tests/testthat/test-curate.R tests/testthat/test-screen.R
git commit -m "feat(screen): the params sublist doc 11 promised, and version 3"
```

---

### Task 2: Anchor extraction

**Files:**
- Create: `R/screen_variants.R`, `tests/testthat/test-screen-variants.R`

**Interfaces:**
- Produces: `screen_anchor_string(seq, ann)` -> character(1) or `NA_character_`. Pure: takes the unit's sequence and its annotation rows, returns the concatenated anchors.
- Consumes: `Biostrings::reverseComplement`, `substr`.

Rules, in order:
1. Keep only rows whose `element` is one of `18S`, `5.8S`, `28S`, with usable `start`/`stop`.
2. No anchor rows at all -> `NA_character_`. A unit missing all three anchors is not comparable (doc 12 section 3.1).
3. Mixed strands across the kept anchors -> `NA_character_`. That is a chimera signal, and concatenating antisense and sense spans produces a string that means nothing.
4. All anchors on `-` -> reverse-complement each span, which is that anchor in sense. Do **not** also reverse the concatenation order: the element sort in rule 5 already puts them in transcription order, and reversing again double-reverses. `orient_sequence()` normally leaves a unit on `+`, so this is defensive, but a unit that escaped orientation must not be compared against an oriented one and read as a paralog on the strength of being antisense.
5. Otherwise concatenate in fixed element order `18S`, `5.8S`, `28S`, using only the anchors present. Fixed order, never coordinate order: a circularly permuted assembly would otherwise produce a different string for the same molecule.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-screen-variants.R`:

```r
# Intragenomic variants: are two on-target units the same molecule assembled
# twice, or two genuinely different rDNA copies from one individual. Restricting
# the comparison to on_target units is what makes the difference intragenomic by
# construction (doc 12 section 1), so no threshold has to separate a paralog
# from a contaminant.

ann_row <- function(element, start, stop, strand = "+") {
  data.frame(element = element, start = start, stop = stop, strand = strand,
    stringsAsFactors = FALSE)
}
anns <- function(...) do.call(rbind, list(...))

test_that("the anchor string is 18S + 5.8S + 28S in element order", {
  seq <- paste(c(strrep("A", 10), strrep("C", 5), strrep("G", 8), strrep("T", 4),
                 strrep("A", 12)), collapse = "")
  # laid out 18S(1-10) ITS1(11-15) 5.8S(16-23) ITS2(24-27) 28S(28-39)
  a <- anns(ann_row("18S", 1, 10), ann_row("ITS1", 11, 15),
            ann_row("5.8S", 16, 23), ann_row("ITS2", 24, 27),
            ann_row("28S", 28, 39))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8), strrep("A", 12)))
})

test_that("element order is fixed, not coordinate order", {
  # a circularly permuted assembly puts 28S first; the same molecule must not
  # produce a different anchor string and read as a paralog of itself
  seq <- paste0(strrep("G", 8), strrep("A", 10))
  a <- anns(ann_row("28S", 1, 8), ann_row("18S", 9, 18))
  expect_equal(RiboPilot:::screen_anchor_string(seq, a),
    paste0(strrep("A", 10), strrep("G", 8)))
})

test_that("a unit with two of three anchors is still comparable", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10), ann_row("28S", 11, 18))
  expect_equal(nchar(RiboPilot:::screen_anchor_string(seq, a)), 18L)
})

test_that("a unit with no anchors at all is not comparable", {
  seq <- strrep("A", 30)
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, anns(ann_row("ITS1", 1, 10)))))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, data.frame())))
})

test_that("an all-minus unit is reverse-complemented rather than compared antisense", {
  fwd <- paste0(strrep("A", 10), strrep("G", 8))
  rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(fwd)))
  a <- anns(ann_row("28S", 1, 8, "-"), ann_row("18S", 9, 18, "-"))
  expect_equal(RiboPilot:::screen_anchor_string(rc, a), fwd)
})

test_that("a mixed-strand unit is not comparable", {
  seq <- paste0(strrep("A", 10), strrep("G", 8))
  a <- anns(ann_row("18S", 1, 10, "+"), ann_row("28S", 11, 18, "-"))
  expect_true(is.na(RiboPilot:::screen_anchor_string(seq, a)))
})
```

- [ ] **Step 2: Run test to verify it fails**

`Rscript -e 'testthat::test_local(filter = "screen-variants")'`
Expected: FAIL, `screen_anchor_string` not found.

- [ ] **Step 3: Implement** `screen_anchor_string()` in a new `R/screen_variants.R` with the file header comment explaining that anchors are compared because the spacers are what varies in BOTH cases and so cannot discriminate between them (doc 12 section 2.2).

- [ ] **Step 4: Verify, including mutation**

Break rule 5 by ordering on `start` instead of the fixed vector; confirm the circular-permutation test fails. Restore.

- [ ] **Step 5: Commit**

```bash
git add R/screen_variants.R tests/testthat/test-screen-variants.R
git commit -m "feat(screen): concatenate a unit's anchors into one comparable string"
```

---

### Task 3: The pure grouping core

**Files:**
- Modify: `R/screen_variants.R`, `tests/testthat/test-screen-variants.R`

**Interfaces:**
- Produces: `screen_anchor_pid(a, b)` -> numeric percent, or `NA_real_`. Returns `100` without aligning when the strings are identical, which is the common isomer case and skips an O(n*m) alignment over roughly 4 kb.
- Produces: `screen_variant_groups(units, min_pid)` where `units` is `data.frame(key, anchor, total_length)` -> `data.frame(key, variant_class, variant_group, variant_identity, pid_to_rep, len_delta)`.

Semantics:
- Units with `NA` anchor are dropped from the comparison and come back with all three columns `NA` (doc 12 section 9: "A unit whose anchors are unannotated is NA, not `singleton`").
- Exactly one comparable unit -> `singleton`, one group, `variant_identity` `NA`. There is no second sequence, and recording 100 would assert a measurement that was never made.
- Two or more comparable units: build a graph over pairs with `pid >= min_pid`, take connected components. Grouping is transitive-by-threshold and can chain (A matches B, B matches C, A does not match C); the group still forms, and `variant_identity` makes the spread visible rather than hiding it. Chained groups are expected to be rare and are not worth a clustering algorithm.
- Group representative = the longest member by `total_length`, ties broken by `key` ascending so the label is deterministic.
- `variant_identity` = each member's anchor identity to its group's representative. The representative's own is `NA`.
- `variant_group` = `g1`, `g2`, ... ordered by representative `total_length` descending, then `key` ascending.
- `variant_class` = `isomer` when the sample has exactly one group, `paralog` when it has two or more. See decision 5 above.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-screen-variants.R`:

```r
mk_units <- function(...) {
  d <- do.call(rbind, lapply(list(...), function(x)
    data.frame(key = x[[1]], anchor = x[[2]], total_length = as.integer(x[[3]]),
      stringsAsFactors = FALSE)))
  d
}
A <- strrep("ACGT", 200)
# one substitution in 800 bp is 99.875 percent, above the 99.5 threshold
A1 <- paste0("T", substr(A, 2, nchar(A)))
# a genuinely divergent copy: every 10th base changed, ~90 percent
B <- paste0(vapply(seq_len(200), function(i)
  if (i %% 10 == 0) "AAAA" else "ACGT", character(1)), collapse = "")

test_that("identical anchors score 100 without an alignment", {
  expect_equal(RiboPilot:::screen_anchor_pid(A, A), 100)
})

test_that("two units with near-identical anchors and different lengths are isomers", {
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 6000L), list("u2", A1, 4200L)), 99.5)
  expect_equal(sort(unique(g$variant_class)), "isomer")
  expect_equal(length(unique(g$variant_group)), 1L)
  # the longest member is the representative and has nothing to compare against
  expect_true(is.na(g$variant_identity[g$key == "u1"]))
  expect_gt(g$variant_identity[g$key == "u2"], 99.5)
})

test_that("identical anchors AND identical length are still isomers", {
  # the spec's table has no row for this; two units that agree on everything are
  # the same molecule assembled twice, for which "pick either" is exactly right
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 5000L), list("u2", A, 5000L)), 99.5)
  expect_equal(sort(unique(g$variant_class)), "isomer")
})

test_that("two units differing beyond the threshold are paralogs in separate groups", {
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 6000L), list("u2", B, 6000L)), 99.5)
  expect_equal(sort(unique(g$variant_class)), "paralog")
  expect_equal(length(unique(g$variant_group)), 2L)
})

test_that("a mixed set forms one isomer group and one paralog group", {
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 6000L), list("u2", A1, 4200L), list("u3", B, 5000L)), 99.5)
  # every comparable unit reads paralog once the sample carries two copies:
  # telling u1 "pick either" would be false when u3 exists
  expect_equal(sort(unique(g$variant_class)), "paralog")
  expect_equal(g$variant_group[g$key == "u1"], g$variant_group[g$key == "u2"])
  expect_false(g$variant_group[g$key == "u1"] == g$variant_group[g$key == "u3"])
})

test_that("a three-unit isomer set forms one group", {
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 6000L), list("u2", A1, 5000L), list("u3", A, 4000L)), 99.5)
  expect_equal(length(unique(g$variant_group)), 1L)
  expect_equal(sort(unique(g$variant_class)), "isomer")
})

test_that("a single comparable unit is a singleton with no identity", {
  g <- RiboPilot:::screen_variant_groups(mk_units(list("u1", A, 6000L)), 99.5)
  expect_equal(g$variant_class, "singleton")
  expect_true(is.na(g$variant_identity))
})

test_that("an uncomparable unit is NA, not singleton", {
  g <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 6000L), list("u2", NA_character_, 900L)), 99.5)
  expect_true(is.na(g$variant_class[g$key == "u2"]))
  expect_true(is.na(g$variant_group[g$key == "u2"]))
  # and it does not make the remaining unit look like it had company
  expect_equal(g$variant_class[g$key == "u1"], "singleton")
})

test_that("group labels are deterministic, longest representative first", {
  g1 <- RiboPilot:::screen_variant_groups(
    mk_units(list("u1", A, 4000L), list("u2", B, 9000L)), 99.5)
  g2 <- RiboPilot:::screen_variant_groups(
    mk_units(list("u2", B, 9000L), list("u1", A, 4000L)), 99.5)
  expect_equal(g1$variant_group[g1$key == "u2"], "g1")
  expect_equal(g2$variant_group[g2$key == "u2"], "g1")
})
```

- [ ] **Step 2: Run test to verify it fails**

- [ ] **Step 3: Implement** both functions. `screen_anchor_pid()` wraps `pairwise_nt(a, b, type = "overlap")` in `tryCatch` returning `NA_real_`, then `pwalign_pid()`; a failed alignment must not sink a sample's screen.

- [ ] **Step 4: Verify, including mutation**

Two mutations, both of which must break a test:
  - Change the representative from longest to first; the deterministic-label test must fail.
  - Change `>= min_pid` to `> min_pid`; add a fixture pair at exactly 99.5 if no existing test bites, since the spec says "at or above".

- [ ] **Step 5: Commit**

```bash
git add R/screen_variants.R tests/testthat/test-screen-variants.R
git commit -m "feat(screen): group a sample's units by anchor identity"
```

---

### Task 4: Schema v15

**Files:**
- Modify: `R/migrate_db.R`
- Test: `tests/testthat/test-migrate.R`

`R/init_db.R` needs no edit: it already executes `SCREEN_DDL`, which is the shared constant. That is the whole point of the shared-DDL convention, and a test must pin that a fresh database and a migrated one agree.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-migrate.R`:

```r
test_that("v15 adds the intragenomic variant columns to screen", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  for (col in c("variant_class", "variant_group", "variant_identity")) {
    expect_true(col %in% DBI::dbListFields(con, "screen"), info = col)
  }
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 15L)
})

test_that("a v14 project with an existing screen row gains the columns without losing it", {
  # additive and idempotent: nothing is dropped or retyped, and a stored verdict
  # survives the migration that adds columns beside it
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID, path, scaffold, screen_verdict, screen_version)
     VALUES ('S1',1,1,'on_target',2)")
  migrate_db(con)  # second run must be a no-op
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen")
  expect_equal(nrow(r), 1L)
  expect_equal(r$screen_verdict, "on_target")
  expect_true(is.na(r$variant_class))
})

```

For the fresh-vs-migrated check do **not** add a new test. `test-migrate.R:75` already has "a fresh new_db() matches a migrated one", and its loop at line 93 reads

```r
  for (tbl in c("blast_screen", "blast_ref_candidates", "blast_hits")) {
```

Add `"screen"` to that vector. One line, and it covers the whole table rather than three named columns.

- [ ] **Step 2: Run test to verify it fails**

- [ ] **Step 3: Implement**

`RIBOPILOT_SCHEMA_VERSION <- 15L`. Add to `SCREEN_DDL` before `time_stamp`:

```
  variant_class     TEXT,
  variant_group     TEXT,
  variant_identity  REAL,
```

Add after the v14 block:

```r
  # --- schema v15: intragenomic isomers and paralogs ---
  # Orthogonal to the verdict: a paralogous variant is emphatically on_target,
  # so the finding gets its own columns rather than a new verdict string. On the
  # screen table because the finding is per unit and the recompute path already
  # exists; a second table would need its own version gate kept in step.
  add_col("screen", "variant_class", "TEXT")
  add_col("screen", "variant_group", "TEXT")
  add_col("screen", "variant_identity", "REAL")
```

No backfill and no `UPDATE`: `SCREEN_VERSION` at `3L` makes `screen_one()` rewrite every stored row on the next run, which is what fills them.

- [ ] **Step 4: Verify, including mutation**

Remove one `add_col` line; confirm the column test fails naming that column. Restore.

- [ ] **Step 5: Commit**

```bash
git add R/migrate_db.R tests/testthat/test-migrate.R
git commit -m "feat(db): variant columns on screen (v15)"
```

---

### Task 5: The per-sample DB pass

**Files:**
- Modify: `R/screen_variants.R`, `R/screen.R`
- Test: `tests/testthat/test-screen-variants.R`

**Interfaces:**
- Produces: `screen_variants(con, id)` -> (invisibly) the number of units classified.
- Called from `screen_sample()` after the per-unit loop and before or after `update_reference_support()`. It must be after the loop: `screen_one()` does `DELETE` then `dbAppendTable` (`R/screen.R:300-310`), so a variant column written before the loop would be silently wiped.

Behaviour:
1. Read the sample's units with `screen_verdict = 'on_target'`, joined to `assemblies` for `sequence` and `length`, and their `annotations` rows. Everything else is a contamination question the screen already owns.
2. Reset `variant_class`, `variant_group`, `variant_identity` to `NULL` for every unit of the sample that is **not** `on_target`, so a unit that changed verdict does not keep a stale finding.
3. Build anchor strings, run `screen_variant_groups()` with `screen_params(con, id)$min_isomer_anchor_pid`, and `UPDATE` the three columns per unit. `UPDATE`, never `DELETE` + append, so the pass cannot destroy the verdict columns `screen_one()` owns.
4. Append the measured identity to `screen_reason` so a curator sees the number rather than only the label (doc 12 section 4). Use a fixed separator and strip any existing clause before appending, so the pass is idempotent and a rerun cannot stack duplicates:

```r
VARIANT_REASON_SEP <- " | variant: "
```

  The clause names the class, the group, the identity to the representative and the total-length difference, for example
  `" | variant: isomer, group g1, anchors 99.87 percent identical to the group representative, total length differs by 1800 bp"`.
5. Cost guard: return early without aligning when every `on_target` unit already carries a non-NULL `variant_class` at `screen_version = SCREEN_VERSION` and no unit was rewritten by this call. Pairwise alignment over roughly 4 kb is O(n*m) and a rerun of an unchanged sample should not pay for it.

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-screen-variants.R`. Extend `mk_screen_db()` from `test-screen.R` rather than duplicating it, or add a local `mk_variant_db()` that also creates `annotations` and gives `assemblies` a `sequence` and a `length`.

```r
test_that("two on-target units with identical anchors come back as isomers, both on_target", {
  # the orthogonality this design turns on: the finding is not a verdict
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A, spacer = 200L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict, c("on_target", "on_target"))
  expect_equal(r$variant_class, c("isomer", "isomer"))
  expect_equal(length(unique(r$variant_group)), 1L)
})

test_that("a paralogous unit stays on_target", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = B, spacer = 2000L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict, c("on_target", "on_target"))
  expect_equal(r$variant_class, c("paralog", "paralog"))
})

test_that("an off_target unit is never compared and its variant columns stay NA", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A, spacer = 200L, tax = FUNGI_T)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_equal(r$screen_verdict[2], "off_target")
  expect_true(is.na(r$variant_class[2]))
  # and the surviving unit is a singleton, not half of an isomer pair
  expect_equal(r$variant_class[1], "singleton")
})

test_that("a unit whose anchors are unannotated is NA, not singleton", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A, spacer = 2000L, annotate = FALSE)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT * FROM screen WHERE ID='S1' ORDER BY scaffold")
  expect_true(is.na(r$variant_class[2]))
  expect_equal(r$variant_class[1], "singleton")
})

test_that("the measured identity appears in the reason string", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A1, spacer = 200L)
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT screen_reason FROM screen WHERE scaffold=2")
  expect_true(grepl("variant: isomer", r$screen_reason, fixed = TRUE))
  expect_true(grepl("99\\.[0-9]+ percent", r$screen_reason))
  # the verdict's own evidence survives beside it
  expect_true(grepl("confidence", r$screen_reason, fixed = TRUE))
})

test_that("rerunning does not stack duplicate variant clauses", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A1, spacer = 200L)
  screen_sample(con, "S1")
  screen_sample(con, "S1", force = TRUE)
  r <- DBI::dbGetQuery(con, "SELECT screen_reason FROM screen WHERE scaffold=2")
  expect_equal(lengths(regmatches(r$screen_reason,
    gregexpr("variant:", r$screen_reason, fixed = TRUE))), 1L)
})

test_that("the variant pass touches only its own sample", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A, spacer = 200L)
  add_unit(con, 1, anchor = A, spacer = 2000L, id = "S2")
  add_unit(con, 2, anchor = A, spacer = 200L, id = "S2")
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM screen WHERE ID='S2'")
  expect_equal(r$n, 0L)
})

test_that("SCREEN_VERSION at 3 forces a v2 row to recompute and gain a class", {
  con <- mk_variant_db(); on.exit(DBI::dbDisconnect(con))
  add_unit(con, 1, anchor = A, spacer = 2000L)
  add_unit(con, 2, anchor = A, spacer = 200L)
  DBI::dbExecute(con,
    "INSERT INTO screen (ID,path,scaffold,screen_verdict,screen_version)
     VALUES ('S1',1,1,'on_target',2)")
  screen_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT screen_version, variant_class FROM screen WHERE scaffold=1")
  expect_equal(r$screen_version, 3L)
  expect_equal(r$variant_class, "isomer")
})
```

Note for the implementer: the cross-sample test above pins that `screen_sample("S1")` writes no `S2` rows. That is a real trap and the equivalent test on the 3b branch had symmetric fixtures and passed with the `WHERE` clause removed. Make the two samples asymmetric so the test bites.

- [ ] **Step 2: Run test to verify it fails**

- [ ] **Step 3: Implement** `screen_variants()` and wire it into `screen_sample()`.

- [ ] **Step 4: Verify, including mutation**

Three mutations, each of which must break a named test:
  - Drop the `screen_verdict = 'on_target'` filter -> the off_target test must fail.
  - Move the `screen_variants()` call to before the unit loop -> every variant column comes back NULL.
  - Remove the strip-before-append in the reason clause -> the no-duplicate-clause test must fail.

- [ ] **Step 5: Commit**

```bash
git add R/screen_variants.R R/screen.R tests/testthat/test-screen-variants.R
git commit -m "feat(screen): classify a sample's on-target units as isomers or paralogs"
```

---

### Task 6: The Assemble scaffold table

**Files:**
- Modify: `R/utils_reactable.R`, `R/app_utils.R`, `R/app_server.R`
- Test: `tests/testthat/test-screen-variants.R` (badge rendering rules only; the reactable itself is not unit-testable here)

An isomer group renders **neutrally**, because it is telling the curator a decision does NOT matter. A paralog group renders as **information**, not as a warning: it is real biology, and colouring it like a defect is what this sub-project exists to stop. Neither may reuse the red used by `off_target` in `rt_screen_badge()` (`R/utils_reactable.R:236`).

- [ ] **Step 1** Add `sc.variant_class, sc.variant_group, sc.variant_identity` to the existing `LEFT JOIN screen sc` select in `unit_evidence()` (`R/app_utils.R:113`).

- [ ] **Step 2** Add `rt_variant_badge()` beside `rt_screen_badge()`. NULL/NA renders empty, not a badge: a sample that has never been screened must not sprout a column of placeholders.

- [ ] **Step 3** Add `variant_class` to the `lead` vector and a `colDef` in `R/app_server.R:885` and `R/app_server.R:905`, named "Variant", showing the identity in the cell title so the number is reachable.

- [ ] **Step 4** Add a test asserting the badge JS contains neither the off-target red (`#f8d7da`) nor the word "warning", and returns empty for NA. Assert on the generated JS string, which is what `rt_screen_badge()`'s own coverage does.

- [ ] **Step 5: Commit**

```bash
git add R/utils_reactable.R R/app_utils.R R/app_server.R tests/testthat/test-screen-variants.R
git commit -m "feat(app): show the intragenomic variant class per scaffold"
```

---

### Task 7: Export records the fact and does not gate

**Files:**
- Modify: `R/export.R`
- Test: `tests/testthat/test-export-units.R`

`export_summary_df()` (`R/export.R:247`) gains `variant_class`, `variant_group`, `variant_identity` and `n_variants`, the number of comparable `on_target` units in the unit's sample. It answers doc 08's open question 7 without imposing a submission policy: the fact is recorded, the choice stays the curator's.

**It must not touch `submission_risks()`.** Adding a variant string there would make a paralog warn on export, which is precisely the behaviour this sub-project exists to remove.

- [ ] **Step 1: Write the failing test**

```r
test_that("export records the variant fact for a paralog and does not gate", {
  # a paralogous unit still exports: the fact is recorded, the choice stays the
  # curator's (doc 12 section 7)
  ...set up a project with two on_target paralogous units, both annotate_lock = 1...
  expect_silent(export_bundle(dir, out_dir = file.path(dir, "out"), report = FALSE, zip = FALSE))
  smry <- utils::read.csv(file.path(dir, "out", "export", "export", "export_summary.csv"))
  expect_true(all(c("variant_class", "variant_group", "variant_identity", "n_variants")
    %in% names(smry)))
  expect_equal(sort(smry$variant_class), c("paralog", "paralog"))
  expect_equal(unique(smry$n_variants), 2L)
  expect_true(all(!nzchar(smry$submission_risk)))
})
```

`expect_silent` is the load-bearing assertion: `export_bundle()` `warning()`s on any non-empty `submission_risk` (`R/export.R:446-451`), so a variant string leaking into that column fails here.

- [ ] **Step 2-4** Fail, implement, verify. Mutation: add `variant_class` to `submission_risks()`'s note and confirm `expect_silent` fails. Restore.

- [ ] **Step 5: Commit**

```bash
git add R/export.R tests/testthat/test-export-units.R
git commit -m "feat(export): record the intragenomic variant fact without gating"
```

---

### Task 8: The vignette

**Files:**
- Modify: `vignettes/Curation-Flags.Rmd`

The triage rule table is **unchanged**: no rule is added and neither class contributes to any verdict. Add a short section after the rule 6 paragraph explaining the Variant column: that an isomer means the decision does not matter, that a paralog is real biology documented in Octocorallia, Scyphozoa and several other groups, that neither is a defect, and that a paralog still exports. State plainly that the IGS is inferred, not measured.

- [ ] **Step 1** Write it. ASCII only.
- [ ] **Step 2** `Rscript -e 'testthat::test_local()'` still 0 failures, and confirm no test asserts a rule count that this text contradicts (`test-triage-rules.R` pins seven rules; it must stay seven).
- [ ] **Step 3: Commit**

```bash
git add vignettes/Curation-Flags.Rmd
git commit -m "docs(vignette): the variant column is information, not a flag"
```

---

### Task 9: Amend the design docs to match what was built

**Files:**
- Modify: `dev/design/12_INTRAGENOMIC_VARIANTS.md`, `dev/design/11_SCREEN_VERDICT_LAYER.md`

The 3b branch's most valuable lesson was that several of the worst defects were defects in the spec rather than in the implementations. Record all five divergences from the top of this plan, with the measurement behind each:

- [ ] Doc 12 section 4: `min_isomer_anchor_pid` is `99.5`, and the `screen` sublist now exists because this sub-project built it.
- [ ] Doc 12 section 4: classification is on identity alone; total length is recorded as evidence. Say why the length column cannot change the answer.
- [ ] Doc 12 section 5: `variant_class` is a per-sample property once two groups exist, and why the isomer message would be false otherwise.
- [ ] Doc 12 section 8: rewrite as **deferred, with the reason**. `reconcile_features()` filters by literal element name, so the barrnap map alone retains nothing; the Nextflow processes declare only `*.itsx.positions.txt` as output, so the ITSx profile never reaches R. Quote both file:line references.
- [ ] Doc 12 section 9: add the mixed-strand and circular-permutation cases, which the spec's test list does not name.
- [ ] Doc 11 section 5: the `screen` sublist exists as of 3c. Remove the claim that it already did.
- [ ] Doc 12 section 10: leave the out-of-scope list intact and add the deferred section 8 items to it.

- [ ] **Commit**

```bash
git add dev/design/
git commit -m "docs(design): record what 3c actually built, and what the spec assumed"
```

---

### Task 10: Whole-branch verification

- [ ] `Rscript -e 'testthat::test_local()'` -> 0 failures, count recorded against the 705 baseline.
- [ ] `Rscript -e 'roxygen2::roxygenise(roclets = "rd")'`, then `git status` -> confirm no unexpected Rd churn and that `NAMESPACE` is untouched.
- [ ] `Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'` -> 0 errors, 0 warnings, and only the two pre-existing NOTEs (crayon/dbplyr imports, `new_db.Rd` line width).
- [ ] `grep -rnP '[^\x00-\x7F]' R/ tests/ inst/ vignettes/ dev/design/12_INTRAGENOMIC_VARIANTS.md` -> no hits.
- [ ] Confirm `RIBOPILOT_SCHEMA_VERSION` is `15L` in exactly one place and `SCREEN_VERSION` is `3L` in exactly one place.
- [ ] Confirm no triage rule was added: `test-triage-rules.R` still pins seven.
- [ ] An independent whole-branch review before the PR, treating every claim in it as a hypothesis to check rather than a fact to accept.
