# Triage Rules and Depth Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the anonymous if/else triage chain with a named, ordered, individually controllable rule registry over the two independent evidence channels, and remove depth from the triage decision.

**Architecture:** A pure rule registry (`R/triage_rules.R`) holds the ordered rule definitions and evaluates them into a per-scaffold reason matrix. A thin settings layer (`R/triage_opts.R`) reads each rule's action (`off` / `flag` / `drop`) from a new `triage_opts` preset table so the Shiny app and the WF1 `triage` process cannot disagree. `R/triage.R` keeps only the database layer. The bulk modal renders the registry as a table of three-way controls with live per-rule counts.

**Tech Stack:** R 4.4, DBI + RSQLite, Shiny + shinyWidgets + reactable, testthat 3e, roxygen2 (Rd only), Nextflow 25.10+.

**Spec:** `dev/design/09_TRIAGE_RULES_AND_DEPTH.md`

## Global Constraints

- ASCII only in all R code, including comments and string literals. No em dashes anywhere.
- Comments are minimal. Explain why, never narrate a bugfix.
- `NAMESPACE` is hand-maintained. Never run `devtools::document()`. Regenerate Rd with `roxygen2::roxygenise(roclets = "rd")` only.
- Every new exported function needs a hand-written `NAMESPACE` export line and a roxygen block.
- `RIBOPILOT_SCHEMA_VERSION` moves from `11L` to `12L` in this plan and nowhere else.
- Migrations only ever add tables and columns. Never drop, never retype.
- Run tests with `Rscript -e 'testthat::test_local(filter = "<pattern>")'` from the repo root.
- Rule keys are exactly: `not_annotated`, `no_screen`, `no_feat_non_rdna`, `no_feat_rdna_hit`, `feat_non_rdna_hit`, `odd_one_out`. In that order. They appear in the registry, the `triage_opts` column names (prefixed `rule_`), `assemblies.triage_rule`, and the modal input ids (prefixed `triage_rule_`).
- Actions are exactly `off`, `flag`, `drop`, lower case.
- Branch is `triage-rules-redesign`, already created from `main` at `dbeaa3e`.

---

## File Structure

**Create:**
- `R/triage_rules.R` - the pure rule registry, reason matrix, evaluator, and per-rule counts. No DBI.
- `R/triage_opts.R` - reads and writes rule actions from the `triage_opts` table. DBI only, no rule logic.
- `tests/testthat/test-triage-rules.R` - registry and evaluator tests.
- `tests/testthat/test-triage-opts.R` - settings resolution and coercion tests.

**Modify:**
- `R/triage.R` - `unit_triage_state()` gains `annotated_n`; `triage_rules()` and `TRIAGE_DEFAULTS` move out; `set_unit_triage()` records the rule; `auto_triage_sample()` and `triage_preview()` resolve settings from the DB and drop `min_depth_frac`.
- `R/migrate_db.R` - `TRIAGE_OPTS_DDL`, the v12 block, version bump.
- `R/init_db.R` - `triage_opts` table and seed, `assemble.triage_opts`, `assemblies.triage_rule`.
- `R/app_modals.R` - `triage_modal()` rebuilt around the rules table.
- `R/app_server.R` - per-rule controls, per-rule live count outputs, Apply writes settings then triages; scaffold table gains the depth sparkline column.
- `R/app_utils.R` - `unit_evidence()` selects `coverage.depths` and downsamples it.
- `R/utils_reactable.R` - `rt_depth_sparkline()` cell renderer.
- `R/coverage.R` - `downsample_depths()` helper.
- `tests/testthat/test-triage.R` - updated for the new rules and the removal of depth.
- `tests/testthat/test-migrate.R` - v12 coverage.
- `NAMESPACE`, `man/` - exports and Rd.

---

### Task 1: The pure rule registry

**Files:**
- Create: `R/triage_rules.R`
- Test: `tests/testthat/test-triage-rules.R`

**Interfaces:**
- Consumes: `unit_concordance(d)` from `R/app_utils.R` (takes a frame with `screen_taxids` and `screen_lineage`, returns `data.frame(sib_rank, sib_taxon, odd_one_out)` with the same row count); `%||%` from `R/utils.R`.
- Produces:
  - `TRIAGE_ACTIONS` - `c("off", "flag", "drop")`
  - `TRIAGE_RULE_DEFS` - ordered list of rule definitions, each `list(key, label, meaning, actions, default, test)`
  - `triage_rule_keys()` -> character vector of keys in evaluation order
  - `triage_actions_default()` -> named list key -> action
  - `normalize_triage_actions(x)` -> named list key -> action, illegal values softened
  - `triage_rule_matrix(d)` -> character matrix, `nrow(d)` x 6, colnames = keys, cell is the reason string or `NA`
  - `triage_rules(d, opts = triage_actions_default())` -> `d` plus `new_triage`, `new_reason`, `new_rule`
  - `triage_rule_counts(d, opts = triage_actions_default())` -> `data.frame(key, label, meaning, action, n_decided, n_matching)`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-triage-rules.R`:

```r
# The rule registry is the whole triage decision. Each rule is tested in
# isolation, and the 2x2 over the two evidence channels is tested exhaustively,
# because a gap there is a scaffold that silently keeps or silently drops.

FISH_T <- "33208;7711;186623;7952;2743709;7954;7955"
COPE_T <- "33208;6656;6854;6843;1;2;3"

# One unit row shaped like unit_triage_state() output.
u <- function(...) {
  r <- list(...)
  data.frame(
    path = 1L, scaffold = r$scaffold %||% 1L, length = 5000L, ignore = 0L,
    triage = NA_character_, triage_reason = NA_character_,
    triage_source = NA_character_, triage_rule = NA_character_,
    depth = 100, screen_species = r$sp %||% NA_character_,
    screen_taxids = r$tx %||% NA_character_, rank_distance = NA_integer_,
    is_nuclear_rrna = r$nuc %||% NA_integer_, sample_taxids = FISH_T,
    sci_name = "Danio rerio", n_features = r$feats %||% 5L,
    annotated_n = r$ann %||% 5L, annotate_lock = 0L, reviewed = "no",
    stringsAsFactors = FALSE)
}
us <- function(...) {
  rows <- list(...)
  do.call(rbind, lapply(seq_along(rows), function(i) {
    do.call(u, c(rows[[i]], list(scaffold = i)))
  }))
}

test_that("the rule keys and their order are fixed", {
  expect_equal(RiboPilot:::triage_rule_keys(),
    c("not_annotated", "no_screen", "no_feat_non_rdna",
      "no_feat_rdna_hit", "feat_non_rdna_hit", "odd_one_out"))
})

test_that("annotation that never ran is flagged, never dropped", {
  # ann = NA means annotate.n_features is NULL: reconcile never reached this
  # unit. annotate_tools carries errorStrategy 'ignore', so this is a tool
  # failure, and the old rule dropped it as if it were evidence.
  d <- RiboPilot:::triage_rules(us(list(ann = NA, feats = 0, nuc = 0, sp = "Homo sapiens")))
  expect_equal(d$new_rule, "not_annotated")
  expect_true(is.na(d$new_triage))
  expect_match(d$new_reason, "annotation has not run")
})

test_that("a scaffold with no usable hit is flagged, never dropped", {
  d <- RiboPilot:::triage_rules(us(list(feats = 5, nuc = NA)))
  expect_equal(d$new_rule, "no_screen")
  expect_true(is.na(d$new_triage))
  expect_match(d$new_reason, "GenBank screen has not run")
})

test_that("the 2x2 over the two channels is exhaustive", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T)))
  expect_equal(d$new_rule[1], "no_feat_non_rdna")
  expect_equal(d$new_triage[1], "drop")
  expect_equal(d$new_rule[2], "no_feat_rdna_hit")
  expect_true(is.na(d$new_triage[2]))
  expect_equal(d$new_rule[3], "feat_non_rdna_hit")
  expect_true(is.na(d$new_triage[3]))
  expect_true(is.na(d$new_rule[4]))
  expect_equal(d$new_triage[4], "keep")
})

test_that("taxonomic distance alone never drops a scaffold", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Cyclops sp.", tx = COPE_T)))
  expect_equal(d$new_rule[4], "odd_one_out")
  expect_false(identical(d$new_triage[4], "drop"))
})

test_that("depth never appears in any reason", {
  d <- RiboPilot:::triage_rules(us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T)))
  # %|NA|% is scalar only (R/utils.R:10), so drop the NAs instead of coalescing
  expect_false(any(grepl("depth", stats::na.omit(d$new_reason), ignore.case = TRUE)))
})

test_that("turning a rule off moves its units to the next MATCHING rule", {
  # the odd unit matches both no_feat_non_rdna and odd_one_out, so turning the
  # first off must fall through to the second rather than to keep
  frame <- function() us(
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 0, nuc = 0, sp = "Cyclops sp.", tx = COPE_T))
  d0 <- RiboPilot:::triage_rules(frame())
  expect_equal(d0$new_rule[4], "no_feat_non_rdna")
  expect_equal(d0$new_triage[4], "drop")

  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "off"))
  d1 <- RiboPilot:::triage_rules(frame(), act)
  expect_equal(d1$new_rule[4], "odd_one_out")
  expect_true(is.na(d1$new_triage[4]))

  # when NO other rule matches, turning the only matching rule off IS a keep
  d2 <- RiboPilot:::triage_rules(
    us(list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T)), act)
  expect_true(is.na(d2$new_rule))
  expect_equal(d2$new_triage, "keep")
})

test_that("an action a rule cannot support is softened to flag", {
  act <- RiboPilot:::normalize_triage_actions(
    list(odd_one_out = "drop", no_screen = "drop", not_annotated = "drop"))
  expect_equal(act$odd_one_out, "flag")
  expect_equal(act$no_screen, "flag")
  expect_equal(act$not_annotated, "flag")

  # off is always legal, and drop is honoured on the combination rules
  act2 <- RiboPilot:::normalize_triage_actions(
    list(odd_one_out = "off", no_feat_rdna_hit = "drop"))
  expect_equal(act2$odd_one_out, "off")
  expect_equal(act2$no_feat_rdna_hit, "drop")
})

test_that("setting every combination rule off or flag yields no drops", {
  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "flag"))
  d <- RiboPilot:::triage_rules(us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T)), act)
  expect_false(any(d$new_triage %in% "drop"))
})

test_that("counts partition on Decided and overlap on Matching", {
  d <- us(
    list(feats = 0, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 0, nuc = 1, sp = "Danio rerio", tx = FISH_T),
    list(feats = 5, nuc = 0, sp = "Homo sapiens", tx = FISH_T),
    list(feats = 5, nuc = 1, sp = "Danio rerio", tx = FISH_T))
  cnt <- RiboPilot:::triage_rule_counts(d)
  expect_equal(nrow(cnt), 6L)
  expect_equal(sum(cnt$n_decided), 3L)          # the clean unit fires nothing
  expect_true(all(cnt$n_matching >= cnt$n_decided))

  # a rule set to off decides nothing but still matches
  act <- RiboPilot:::normalize_triage_actions(list(no_feat_non_rdna = "off"))
  cnt2 <- RiboPilot:::triage_rule_counts(d, act)
  off <- cnt2[cnt2$key == "no_feat_non_rdna", ]
  expect_equal(off$n_decided, 0L)
  expect_equal(off$n_matching, 1L)
})

test_that("an empty frame returns empty columns rather than erroring", {
  d <- RiboPilot:::triage_rules(us()[0, ])
  expect_equal(nrow(d), 0L)
  expect_true(all(c("new_triage", "new_reason", "new_rule") %in% names(d)))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage-rules")'`
Expected: FAIL, `could not find function "triage_rule_keys"`.

- [ ] **Step 3: Write the registry**

Create `R/triage_rules.R`:

```r
# Named, ordered triage rules over the two independent evidence channels:
# annotation (barrnap + ITSx, local) and the GenBank screen (remote). Their
# independence is the safety argument for dropping on their conjunction, so a
# rule may only drop when BOTH channels are present. Rules 1-2 assert that a
# channel is missing and therefore must precede the combination rules.

#' Actions a triage rule may be set to.
#' @noRd
TRIAGE_ACTIONS <- c("off", "flag", "drop")

#' Ordered triage rule definitions
#'
#' Each entry is list(key, label, meaning, actions, default, test), where `test`
#' is a pure function (d, i, ctx) returning a reason string or NA.
#' @noRd
TRIAGE_RULE_DEFS <- list(
  list(
    key = "not_annotated",
    label = "Annotation has not run",
    meaning = "Neither annotation tool has reported on this scaffold yet.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (isTRUE(ctx$annotated[i])) return(NA_character_)
      "annotation has not run for this scaffold"
    }
  ),
  list(
    key = "no_screen",
    label = "GenBank screen has not run",
    meaning = "No usable core_nt hit is recorded, so the scaffold cannot be judged against GenBank.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (isTRUE(ctx$screened[i])) return(NA_character_)
      "the GenBank screen has not run for this scaffold"
    }
  ),
  list(
    key = "no_feat_non_rdna",
    label = "No annotations, non-rDNA hit",
    meaning = "Annotation found no rRNA or ITS feature, and the best GenBank match is not a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "drop",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] == 0) || !isTRUE(ctx$nuc[i] == 0)) return(NA_character_)
      sprintf(paste("annotation found no rRNA or ITS feature, and the best GenBank match",
        "(%s) is not a nuclear rDNA record"), ctx$sp[i])
    }
  ),
  list(
    key = "no_feat_rdna_hit",
    label = "No annotations but rDNA hit",
    meaning = "Annotation found no rRNA or ITS feature, though the best GenBank match is a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] == 0) || !isTRUE(ctx$nuc[i] == 1)) return(NA_character_)
      sprintf(paste("annotation found no rRNA or ITS feature, though the best GenBank match",
        "(%s) is a nuclear rDNA record"), ctx$sp[i])
    }
  ),
  list(
    key = "feat_non_rdna_hit",
    label = "Has annotations but non-rDNA hit",
    meaning = "The scaffold has a feature model, but its best GenBank match is not a nuclear rDNA record.",
    actions = TRIAGE_ACTIONS,
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$annotated[i]) || !isTRUE(ctx$screened[i])) return(NA_character_)
      if (!isTRUE(ctx$n_feat[i] > 0) || !isTRUE(ctx$nuc[i] == 0)) return(NA_character_)
      sprintf("best GenBank match (%s) is not a nuclear rDNA record", ctx$sp[i])
    }
  ),
  list(
    key = "odd_one_out",
    label = "Disagrees with sibling scaffolds",
    meaning = "This scaffold's top hit disagrees with the majority of this sample's other scaffolds.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (!isTRUE(ctx$screened[i]) || !isTRUE(ctx$sib$odd_one_out[i])) return(NA_character_)
      sprintf("top GenBank hit (%s) differs from this sample's other scaffolds (%s at %s)",
        ctx$sp[i], ctx$sib$sib_taxon[i] %|NA|% "the majority",
        ctx$sib$sib_rank[i] %|NA|% "an agreed rank")
    }
  )
)

#' Rule keys in evaluation order
#' @noRd
triage_rule_keys <- function() {
  vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "key")
}

#' The shipped action for every rule
#' @noRd
triage_actions_default <- function() {
  stats::setNames(lapply(TRIAGE_RULE_DEFS, `[[`, "default"), triage_rule_keys())
}

#' Coerce a partial or hand-edited action set into a complete, legal one
#'
#' An action a rule cannot support is softened to "flag" rather than honoured or
#' rejected: a hand-edited database must not be able to arm a drop the evidence
#' cannot support, and must not be able to make the app error either.
#' @noRd
normalize_triage_actions <- function(x) {
  out <- triage_actions_default()
  if (is.null(x)) return(out)
  for (def in TRIAGE_RULE_DEFS) {
    v <- x[[def$key]]
    if (is.null(v) || length(v) == 0 || is.na(v[1])) next
    v <- tolower(trimws(as.character(v[1])))
    if (!v %in% TRIAGE_ACTIONS) next
    if (!v %in% def$actions) v <- "flag"
    out[[def$key]] <- v
  }
  out
}

#' Shared per-sample context every rule test reads
#' @noRd
triage_rule_ctx <- function(d) {
  n <- nrow(d)
  nuc <- d$is_nuclear_rrna %||% rep(NA_integer_, n)
  sp <- d$screen_species %||% rep(NA_character_, n)
  list(
    # annotate.n_features is written by reconcile_one() on every run, including
    # one that found nothing, and only by reconcile, so it survives a re-queue.
    annotated = !is.na(d$annotated_n %||% rep(NA_integer_, n)),
    # a screen row with no classifiable hit is not a screen for rule purposes
    screened = !is.na(nuc),
    n_feat = d$n_features %||% rep(0L, n),
    nuc = nuc,
    sp = ifelse(is.na(sp), "unknown", sp),
    sib = unit_concordance(data.frame(
      screen_taxids = d$screen_taxids %||% rep(NA_character_, n),
      screen_lineage = rep(NA_character_, n), stringsAsFactors = FALSE))
  )
}

#' Every rule's reason for every unit, independent of order and of settings
#' @return character matrix nrow(d) x length(TRIAGE_RULE_DEFS), colnames = keys.
#' @noRd
triage_rule_matrix <- function(d) {
  keys <- triage_rule_keys()
  n <- nrow(d)
  m <- matrix(NA_character_, nrow = n, ncol = length(keys),
    dimnames = list(NULL, keys))
  if (n == 0) return(m)
  ctx <- triage_rule_ctx(d)
  for (j in seq_along(TRIAGE_RULE_DEFS)) {
    tst <- TRIAGE_RULE_DEFS[[j]]$test
    for (i in seq_len(n)) m[i, j] <- tst(d, i, ctx)
  }
  m
}

#' Propose a triage verdict for every unit of one sample
#'
#' First enabled rule that matches wins. Nothing matching means keep.
#'
#' @param d Frame from [unit_triage_state()].
#' @param opts Named list key -> action; see [normalize_triage_actions()].
#' @return `d` with `new_triage`, `new_reason` and `new_rule` columns.
#' @noRd
triage_rules <- function(d, opts = triage_actions_default()) {
  n <- nrow(d)
  d$new_triage <- rep(NA_character_, n)
  d$new_reason <- rep(NA_character_, n)
  d$new_rule <- rep(NA_character_, n)
  if (n == 0) return(d)
  act <- normalize_triage_actions(opts)
  m <- triage_rule_matrix(d)
  keys <- colnames(m)
  for (i in seq_len(n)) {
    hit <- NA_integer_
    for (j in seq_along(keys)) {
      if (identical(act[[keys[j]]], "off")) next
      if (is.na(m[i, j])) next
      hit <- j
      break
    }
    if (is.na(hit)) {
      d$new_triage[i] <- "keep"
      next
    }
    d$new_rule[i] <- keys[hit]
    d$new_reason[i] <- m[i, hit]
    if (identical(act[[keys[hit]]], "drop")) d$new_triage[i] <- "drop"
  }
  d
}

#' Per-rule counts for the bulk triage modal
#'
#' `n_decided` is first-match and partitions the units. `n_matching` counts
#' every unit meeting the rule, including ones a rule above acted on first, so
#' a rule set to off still reports what turning it back on would reach.
#' @noRd
triage_rule_counts <- function(d, opts = triage_actions_default()) {
  act <- normalize_triage_actions(opts)
  keys <- triage_rule_keys()
  out <- data.frame(
    key = keys,
    label = vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "label"),
    meaning = vapply(TRIAGE_RULE_DEFS, `[[`, character(1), "meaning"),
    action = vapply(keys, function(k) act[[k]], character(1)),
    n_decided = 0L, n_matching = 0L, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  if (nrow(d) == 0) return(out)
  out$n_matching <- as.integer(colSums(!is.na(triage_rule_matrix(d))))
  fired <- triage_rules(d, act)$new_rule
  out$n_decided <- as.integer(vapply(keys, function(k) sum(fired %in% k), integer(1)))
  out
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "triage-rules")'`
Expected: PASS, 10 tests.

- [ ] **Step 5: Commit**

```bash
git add R/triage_rules.R tests/testthat/test-triage-rules.R
git commit -m "feat(triage): named, ordered rule registry over both evidence channels"
```

---

### Task 2: Schema v12

**Files:**
- Modify: `R/migrate_db.R` (add `TRIAGE_OPTS_DDL` near `BLAST_HITS_DDL` at line 18, add the v12 block after the v11 block at line 310-328, bump `RIBOPILOT_SCHEMA_VERSION` at line 8)
- Modify: `R/init_db.R` (add the table after `curate_opts` at line 550-563, add `assemble.triage_opts` to the `assemble` DDL at line 164, add `assemblies.triage_rule` to the `assemblies` DDL at line 276)
- Test: `tests/testthat/test-migrate.R`

**Interfaces:**
- Consumes: nothing from Task 1 at runtime, but the column names must match `triage_rule_keys()` exactly, prefixed `rule_`.
- Produces: table `triage_opts(triage_opts, rule_not_annotated, rule_no_screen, rule_no_feat_non_rdna, rule_no_feat_rdna_hit, rule_feat_non_rdna_hit, rule_odd_one_out)`; columns `assemble.triage_opts TEXT DEFAULT 'default'` and `assemblies.triage_rule TEXT`; constant `TRIAGE_OPTS_DDL`.

- [ ] **Step 1: Write the failing test**

The existing helper in this file is `mk_v8_db()` (line 7), and it does **not**
create an `assemble` table, so `add_col("assemble", "triage_opts", ...)` would
silently no-op. Add the table to the helper first, since a real v8 project has
one, inserting after the `assemblies` CREATE:

```r
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT PRIMARY KEY, assemble_switch INT,
    assemble_lock INT, assemble_opts TEXT, blast_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1',2,0,'default','default')")
```

Then append to `tests/testthat/test-migrate.R`:

```r
test_that("v12 adds triage_opts, its seed row, and the rule columns", {
  con <- mk_v8_db()
  on.exit(DBI::dbDisconnect(con))
  migrate_db(con)

  expect_true("triage_opts" %in% DBI::dbListTables(con))
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM triage_opts")$n, 1L)

  # every rule key has a column, and the seed carries the shipped defaults
  r <- DBI::dbGetQuery(con, "SELECT * FROM triage_opts WHERE triage_opts = 'default'")
  for (k in RiboPilot:::triage_rule_keys()) {
    expect_true(paste0("rule_", k) %in% names(r))
  }
  expect_equal(r$rule_no_feat_non_rdna, "drop")
  expect_equal(r$rule_odd_one_out, "flag")

  expect_true("triage_opts" %in% DBI::dbListFields(con, "assemble"))
  expect_true("triage_rule" %in% DBI::dbListFields(con, "assemblies"))
  expect_equal(DBI::dbGetQuery(con, "PRAGMA user_version")[[1]], 12L)
})
```

Adding `assemble` to `mk_v8_db()` also exercises the existing v9-v11 migrations
against it, so re-run the whole file, not just the new test.

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: FAIL, `triage_opts` not in `dbListTables`.

- [ ] **Step 3: Implement the migration**

In `R/migrate_db.R`, change line 8:

```r
RIBOPILOT_SCHEMA_VERSION <- 12L
```

Add after the `BLAST_HITS_DDL` block (after line 43):

```r
#' DDL for the triage_opts table (per-rule action presets).
#' Shared by new_db() and the v12 migration so both agree on the schema.
#' @noRd
TRIAGE_OPTS_DDL <- "CREATE TABLE IF NOT EXISTS triage_opts (
  triage_opts            TEXT NOT NULL,
  rule_not_annotated     TEXT,
  rule_no_screen         TEXT,
  rule_no_feat_non_rdna  TEXT,
  rule_no_feat_rdna_hit  TEXT,
  rule_feat_non_rdna_hit TEXT,
  rule_odd_one_out       TEXT,
  PRIMARY KEY (triage_opts)
);"

#' Seed row for triage_opts, matching TRIAGE_RULE_DEFS defaults.
#' @noRd
TRIAGE_OPTS_SEED <- "INSERT INTO triage_opts
  (triage_opts, rule_not_annotated, rule_no_screen, rule_no_feat_non_rdna,
   rule_no_feat_rdna_hit, rule_feat_non_rdna_hit, rule_odd_one_out)
  VALUES ('default', 'flag', 'flag', 'drop', 'flag', 'flag', 'flag');"
```

Add after the v11 block (after line 328, before the `user_version` write):

```r
  # --- schema v12: named, individually controllable triage rules ---
  # The rules were an anonymous if/else chain with no per-rule setting, so the
  # WF1 triage process and the app could not be made to agree on anything the
  # curator changed. triage_opts holds one action per rule; triage_rule records
  # which rule produced a unit's current verdict.
  if (!has_tbl("triage_opts")) DBI::dbExecute(con, TRIAGE_OPTS_DDL)
  if (DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM triage_opts")$n == 0) {
    DBI::dbExecute(con, TRIAGE_OPTS_SEED)
  }
  add_col("assemble", "triage_opts", "TEXT DEFAULT 'default'")
  add_col("assemblies", "triage_rule", "TEXT")
  # unit_evidence() hard-selects coverage.depths for the sparkline, so its
  # absence would break the whole scaffold query rather than one cell
  add_col("coverage", "depths", "TEXT")
  if (ver < 12L && has_col("assemble", "triage_opts")) {
    DBI::dbExecute(con, "UPDATE assemble SET triage_opts = 'default' WHERE triage_opts IS NULL")
  }
```

In `R/init_db.R`, add `triage_opts TEXT DEFAULT 'default',` to the `assemble` CREATE TABLE (line 164 block), add `triage_rule TEXT,` to the `assemblies` CREATE TABLE (line 276 block, next to `triage_source`), and add after the `curate_opts` seed (line 563):

```r
  # triage_opts (per-rule triage actions; see TRIAGE_RULE_DEFS) ----
  DBI::dbExecute(con, TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, TRIAGE_OPTS_SEED)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "migrate")'`
Expected: PASS.

Then check nothing else broke: `Rscript -e 'testthat::test_local()'`
Expected: PASS except `test-triage.R`, which Task 4 updates.

- [ ] **Step 5: Commit**

```bash
git add R/migrate_db.R R/init_db.R tests/testthat/test-migrate.R
git commit -m "feat(db): triage_opts preset table and triage_rule column (v12)"
```

---

### Task 3: Settings resolution

**Files:**
- Create: `R/triage_opts.R`
- Test: `tests/testthat/test-triage-opts.R`

**Interfaces:**
- Consumes: `triage_rule_keys()`, `triage_actions_default()`, `normalize_triage_actions()` from Task 1; `TRIAGE_OPTS_DDL`, `TRIAGE_OPTS_SEED` from Task 2.
- Produces:
  - `triage_opts_row(con, preset = "default")` -> named list key -> action
  - `triage_opts_for(con, id)` -> named list key -> action, following `assemble.triage_opts`
  - `write_triage_opts(con, actions, preset = "default")` -> invisibly `preset`

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-triage-opts.R`:

```r
# Settings must resolve identically for the app and for the WF1 triage process,
# and a hand-edited database must never be able to arm an unsupported drop.

mk_opts_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, triage_opts TEXT)")
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  con
}

test_that("the seeded default resolves to the shipped actions", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::triage_opts_for(con, "S1"),
    RiboPilot:::triage_actions_default())
})

test_that("a sample follows its own preset", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_triage_opts(con,
    list(no_feat_non_rdna = "flag", odd_one_out = "off"), preset = "cautious")
  DBI::dbExecute(con, "UPDATE assemble SET triage_opts = 'cautious' WHERE ID = 'S1'")
  a <- RiboPilot:::triage_opts_for(con, "S1")
  expect_equal(a$no_feat_non_rdna, "flag")
  expect_equal(a$odd_one_out, "off")
  expect_equal(a$no_screen, "flag")
})

test_that("an unsupported drop written into the DB is read back as flag", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "UPDATE triage_opts SET rule_odd_one_out = 'drop' WHERE triage_opts = 'default'")
  expect_equal(RiboPilot:::triage_opts_for(con, "S1")$odd_one_out, "flag")
})

test_that("a missing table, preset or sample falls back to the defaults", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  expect_equal(RiboPilot:::triage_opts_row(con, "nope"),
    RiboPilot:::triage_actions_default())
  expect_equal(RiboPilot:::triage_opts_for(con, "NOT_A_SAMPLE"),
    RiboPilot:::triage_actions_default())

  bare <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(bare), add = TRUE)
  expect_equal(RiboPilot:::triage_opts_for(bare, "S1"),
    RiboPilot:::triage_actions_default())
})

test_that("write then read round-trips every rule", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  want <- RiboPilot:::triage_actions_default()
  want$not_annotated <- "off"
  want$no_feat_rdna_hit <- "drop"
  RiboPilot:::write_triage_opts(con, want)
  expect_equal(RiboPilot:::triage_opts_row(con, "default"), want)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage-opts")'`
Expected: FAIL, `could not find function "triage_opts_for"`.

- [ ] **Step 3: Implement**

Create `R/triage_opts.R`:

```r
# Per-rule triage actions, read from the triage_opts preset table so the Shiny
# app's Apply and the WF1 triage process cannot disagree. All resolution is
# defensive: an old project, a missing preset or a hand-edited value falls back
# to the shipped defaults rather than erroring or arming an unsupported drop.

#' Column name holding a rule's action
#' @noRd
triage_opts_col <- function(key) paste0("rule_", key)

#' Rule actions from one preset
#' @param con Open project DB connection.
#' @param preset Preset name.
#' @return Named list key -> action.
#' @noRd
triage_opts_row <- function(con, preset = "default") {
  r <- tryCatch(
    DBI::dbGetQuery(con, "SELECT * FROM triage_opts WHERE triage_opts = ?",
      params = list(preset %||% "default")),
    error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(triage_actions_default())
  x <- list()
  for (k in triage_rule_keys()) {
    col <- triage_opts_col(k)
    if (col %in% names(r)) x[[k]] <- r[[col]][1]
  }
  normalize_triage_actions(x)
}

#' Rule actions for one sample, following assemble.triage_opts
#' @param con Open project DB connection.
#' @param id Sample ID.
#' @return Named list key -> action.
#' @noRd
triage_opts_for <- function(con, id) {
  p <- tryCatch(
    DBI::dbGetQuery(con, "SELECT triage_opts AS t FROM assemble WHERE ID = ?",
      params = list(id))$t,
    error = function(e) character(0))
  triage_opts_row(con, if (length(p) && !is.na(p[1]) && nzchar(p[1])) p[1] else "default")
}

#' Write rule actions into a preset, creating it if needed
#' @param con Open project DB connection.
#' @param actions Named list key -> action; missing keys take the default.
#' @param preset Preset name.
#' @return (invisibly) the preset name.
#' @noRd
write_triage_opts <- function(con, actions, preset = "default") {
  a <- normalize_triage_actions(actions)
  keys <- triage_rule_keys()
  cols <- vapply(keys, triage_opts_col, character(1))
  # one transaction: a failed insert after a committed delete would silently
  # discard a customised preset and leave the project on the defaults
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM triage_opts WHERE triage_opts = ?", params = list(preset))
    DBI::dbExecute(con,
      sprintf("INSERT INTO triage_opts (triage_opts, %s) VALUES (%s)",
        paste(cols, collapse = ", "),
        paste(rep("?", length(cols) + 1), collapse = ", ")),
      params = c(list(preset), unname(a[keys])))
  })
  invisible(preset)
}
```

The `sprintf()` interpolates only `cols`, which comes from `triage_rule_keys()`,
and the placeholder string. Every caller-supplied value goes through `params`.

Add one test pinning the rollback, using a trigger to induce a failing insert:

```r
test_that("a failed write leaves the stored preset intact", {
  con <- mk_opts_db(); on.exit(DBI::dbDisconnect(con))
  RiboPilot:::write_triage_opts(con, list(odd_one_out = "off"), preset = "keepme")
  DBI::dbExecute(con,
    "CREATE TRIGGER boom BEFORE INSERT ON triage_opts
     BEGIN SELECT RAISE(ABORT, 'boom'); END;")
  expect_error(RiboPilot:::write_triage_opts(con, list(no_screen = "off"), preset = "keepme"))
  expect_equal(RiboPilot:::triage_opts_row(con, "keepme")$odd_one_out, "off")
})
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "triage-opts")'`
Expected: PASS, 5 tests.

- [ ] **Step 5: Commit**

```bash
git add R/triage_opts.R tests/testthat/test-triage-opts.R
git commit -m "feat(triage): resolve per-rule actions from the triage_opts preset"
```

---

### Task 4: Rewire the database layer

**Files:**
- Modify: `R/triage.R` (whole file; `triage_rules()` and `TRIAGE_DEFAULTS` are deleted from here, they now live in `R/triage_rules.R`)
- Modify: `tests/testthat/test-triage.R`
- Modify: `NAMESPACE`, `man/`

**Interfaces:**
- Consumes: `triage_rules()`, `normalize_triage_actions()` from Task 1; `triage_opts_for()` from Task 3.
- Produces:
  - `unit_triage_state(con, id)` gains an `annotated_n` column (integer or NA)
  - `set_unit_triage(con, id, path, scaffold, verdict, reason, source, rule = NA_character_)`
  - `auto_triage_sample(con, id, opts = NULL, overwrite_user = FALSE)`
  - `triage_preview(con, ids, opts = NULL, overwrite_user = FALSE)` gains `by_rule`, loses nothing else

- [ ] **Step 1: Write the failing test**

Rewrite the pure-rule section of `tests/testthat/test-triage.R`. Delete the four `triage_rules()` tests at lines 24-79 (Task 1's file now owns them) and the depth assertions at lines 66-72. Replace the DB helper's `annotate` rows so `n_features` is populated, and add these tests:

```r
test_that("unit_triage_state reports whether annotation actually ran", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  d <- RiboPilot:::unit_triage_state(con, "S1")
  expect_true("annotated_n" %in% names(d))
  expect_false(any(is.na(d$annotated_n)))

  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL WHERE scaffold = 2")
  d2 <- RiboPilot:::unit_triage_state(con, "S1")
  expect_true(is.na(d2$annotated_n[2]))
})

test_that("a unit whose annotation never ran is never dropped", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  # the exact shape that used to drop on a tool failure: no annotation rows and
  # a best hit that fails the nuclear-rDNA title regex
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = NULL WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule, ignore FROM assemblies WHERE scaffold = 1")
  expect_true(is.na(r$triage))
  expect_equal(r$triage_rule, "not_annotated")
  expect_equal(r$ignore, 0L)
})

test_that("the rule that fired is recorded alongside the reason", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0 WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule FROM assemblies WHERE scaffold = 1")
  expect_equal(r$triage, "drop")
  expect_equal(r$triage_rule, "no_feat_non_rdna")
})

test_that("auto-triage honours the sample's triage_opts preset", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0 WHERE scaffold = 1")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens' WHERE scaffold = 1")
  RiboPilot:::write_triage_opts(con, list(no_feat_non_rdna = "flag"))

  auto_triage_sample(con, "S1")
  r <- DBI::dbGetQuery(con, "SELECT triage, triage_rule FROM assemblies WHERE scaffold = 1")
  expect_true(is.na(r$triage))
  expect_equal(r$triage_rule, "no_feat_non_rdna")
})

test_that("preview reports per-rule counts alongside the totals", {
  con <- mk_triage_db(); on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "DELETE FROM annotations")
  DBI::dbExecute(con, "UPDATE annotate SET n_features = 0")
  DBI::dbExecute(con, "UPDATE blast_screen SET is_nuclear_rrna = 0, species = 'Homo sapiens'")

  p <- triage_preview(con, "S1")
  expect_true(is.data.frame(p$by_rule))
  expect_equal(nrow(p$by_rule), 6L)
  expect_equal(sum(p$by_rule$n_decided), p$n_units)
  expect_true(all(p$by_rule$n_matching >= p$by_rule$n_decided))
  expect_equal(p$by_rule$n_decided[p$by_rule$key == "no_feat_non_rdna"], 3L)
})
```

The `mk_triage_db()` helper at line 83 must gain `n_features` on `annotate`, a `triage_rule` column on `assemblies`, and the `triage_opts` table:

```r
mk_triage_db <- function() {
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  DBI::dbExecute(con, "CREATE TABLE assemblies(ID TEXT, path INT, scaffold INT, length INT,
    ignore INT, triage TEXT, triage_reason TEXT, triage_source TEXT, triage_rule TEXT)")
  DBI::dbExecute(con, "CREATE TABLE annotate(ID TEXT, path INT, scaffold INT, annotate_switch INT,
    annotate_lock INT, annotate_opts TEXT, curate_opts TEXT, reviewed TEXT, n_features INT)")
  DBI::dbExecute(con, "CREATE TABLE annotations(ID TEXT, path INT, scaffold INT, element TEXT)")
  DBI::dbExecute(con, "CREATE TABLE coverage(ID TEXT, path INT, scaffold INT, mean_depth REAL, median_depth REAL)")
  DBI::dbExecute(con, "CREATE TABLE blast_screen(ID TEXT, path INT, scaffold INT, species TEXT,
    lineage_taxids TEXT, rank_distance INT, is_nuclear_rrna INT)")
  DBI::dbExecute(con, "CREATE TABLE samples(ID TEXT, sci_name TEXT, lineage_taxids TEXT)")
  DBI::dbExecute(con, "CREATE TABLE assemble(ID TEXT, triage_opts TEXT)")
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_DDL)
  DBI::dbExecute(con, RiboPilot:::TRIAGE_OPTS_SEED)
  DBI::dbExecute(con, "INSERT INTO assemble VALUES('S1','default')")
  DBI::dbExecute(con, "INSERT INTO samples VALUES('S1','Danio rerio',?)", params = list(FISH_T))
  for (s in 1:3) {
    DBI::dbExecute(con, "INSERT INTO assemblies VALUES('S1',1,?,5000,0,NULL,NULL,NULL,NULL)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotate VALUES('S1',1,?,2,0,'default','default','no',1)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO coverage VALUES('S1',1,?,100,100)", params = list(s))
    DBI::dbExecute(con, "INSERT INTO annotations VALUES('S1',1,?,'18S')", params = list(s))
    DBI::dbExecute(con, "INSERT INTO blast_screen VALUES('S1',1,?,'Danio rerio',?,0,1)",
      params = list(s, FISH_T))
  }
  con
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage")'`
Expected: FAIL, `annotated_n` not in names, and `triage_rule` unknown column.

- [ ] **Step 3: Implement**

In `R/triage.R`:

Delete `TRIAGE_DEFAULTS` (line 85-91) and `triage_rules()` (line 93-169) entirely; they live in `R/triage_rules.R` now.

Add `annotated_n` to `unit_triage_state()`'s SELECT, next to the existing annotate subqueries:

```r
       (SELECT an.n_features FROM annotate an
          WHERE an.ID = a.ID AND an.path = a.path AND an.scaffold = a.scaffold) AS annotated_n,
```

Change `set_unit_triage()` to take and store the rule:

```r
set_unit_triage <- function(con, id, path, scaffold, verdict,
                            reason = NA_character_, source = c("user", "auto"),
                            rule = NA_character_) {
  source <- match.arg(source)
  path <- as.integer(path); scaffold <- as.integer(scaffold)
  if (length(verdict) == 0 || is.na(verdict[1]) || !nzchar(verdict[1])) {
    verdict <- NA_character_
  } else if (verdict[1] %in% TRIAGE_VERDICTS) {
    verdict <- verdict[1]
  } else {
    stop("verdict must be one of: ", paste(TRIAGE_VERDICTS, collapse = ", "),
      ", or NA for undecided", call. = FALSE)
  }
  ig <- as.integer(identical(verdict, "drop"))
  DBI::dbExecute(con,
    "UPDATE assemblies SET triage = ?, triage_reason = ?, triage_source = ?,
       triage_rule = ?, ignore = ?
     WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(verdict, reason %||% NA_character_, source,
      rule %||% NA_character_, ig, id, path, scaffold))
  reconcile_ignore_unit(con, id, path, scaffold, ig)
  invisible(ig)
}
```

Change `auto_triage_sample()`:

```r
auto_triage_sample <- function(con, id, opts = NULL, overwrite_user = FALSE) {
  d <- unit_triage_state(con, id)
  if (nrow(d) == 0) {
    return(invisible(0L))
  }
  act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
  d <- triage_rules(d, act)
  changed <- 0L
  for (i in seq_len(nrow(d))) {
    if (!isTRUE(overwrite_user) && identical(d$triage_source[i], "user")) next
    if (isTRUE(d$annotate_lock[i] == 1)) next
    same <- identical(d$triage[i] %|NA|% "", d$new_triage[i] %|NA|% "") &&
      identical(d$triage_reason[i] %|NA|% "", d$new_reason[i] %|NA|% "") &&
      identical(d$triage_rule[i] %|NA|% "", d$new_rule[i] %|NA|% "")
    if (same) next
    set_unit_triage(con, id, d$path[i], d$scaffold[i], d$new_triage[i],
      reason = d$new_reason[i], source = "auto", rule = d$new_rule[i])
    changed <- changed + 1L
  }
  invisible(changed)
}
```

`unit_triage_state()` must also select `a.triage_rule` for that sameness check; add it next to `a.triage_source`.

Change `triage_preview()`:

```r
triage_preview <- function(con, ids, opts = NULL, overwrite_user = FALSE) {
  keys <- triage_rule_keys()
  res <- list(n_units = 0L, n_drop = 0L, n_undecided = 0L, n_keep = 0L,
    n_samples = 0L, emptied = character(0), n_user_skipped = 0L,
    n_locked_skipped = 0L, by_rule = NULL)
  acc <- NULL
  for (id in unique(ids)) {
    d <- unit_triage_state(con, id)
    if (nrow(d) == 0) next
    act <- normalize_triage_actions(opts %||% triage_opts_for(con, id))
    d <- triage_rules(d, act)
    skip_user <- !isTRUE(overwrite_user) & identical_user(d$triage_source)
    skip_lock <- !is.na(d$annotate_lock) & d$annotate_lock == 1
    skipped <- skip_user | skip_lock
    # count only the units Apply would touch, so the per-rule table cannot
    # disagree with the totals rendered beside it
    cnt <- triage_rule_counts(d[!skipped, , drop = FALSE], act)
    if (is.null(acc)) {
      acc <- cnt
    } else {
      acc$n_decided <- acc$n_decided + cnt$n_decided
      acc$n_matching <- acc$n_matching + cnt$n_matching
      # samples may follow different presets; never report one sample's action
      # as if it governed the whole selection
      acc$action[acc$action != cnt$action] <- "mixed"
    }
    res$n_user_skipped <- res$n_user_skipped + sum(skip_user & !skip_lock)
    res$n_locked_skipped <- res$n_locked_skipped + sum(skip_lock)
    eff <- ifelse(skipped, d$triage, d$new_triage)
    res$n_units <- res$n_units + nrow(d)
    res$n_drop <- res$n_drop + sum(eff %in% "drop")
    res$n_keep <- res$n_keep + sum(eff %in% "keep")
    res$n_undecided <- res$n_undecided + sum(is.na(eff))
    res$n_samples <- res$n_samples + 1L
    if (all(eff %in% "drop")) res$emptied <- c(res$emptied, id)
  }
  res$by_rule <- acc %||% triage_rule_counts(
    data.frame(), normalize_triage_actions(opts %||% triage_opts_row(con)))
  res
}
```

Note the roxygen `@return` for `triage_preview()` must gain `by_rule`, and `@param opts` must read "Named list of per-rule actions; NULL reads the sample's triage_opts preset." Same for `auto_triage_sample()`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "triage")'`
Expected: PASS.

Then the whole suite: `Rscript -e 'testthat::test_local()'`
Expected: PASS.

- [ ] **Step 5: Regenerate Rd and commit**

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
git add R/triage.R tests/testthat/test-triage.R man/
git commit -m "feat(triage): gate rules on annotation having run; record the rule that fired"
```

---

### Task 5: The bulk triage modal

**Files:**
- Modify: `R/app_modals.R:81-108` (`triage_modal()`)
- Modify: `R/app_server.R:1059-1101` (the bulk triage block)
- Test: `tests/testthat/test-triage-modal.R` (create)

**Interfaces:**
- Consumes: `TRIAGE_RULE_DEFS`, `triage_rule_keys()`, `triage_actions_default()` from Task 1; `triage_opts_for()`, `write_triage_opts()` from Task 3; `triage_preview()` with `by_rule` from Task 4.
- Produces: `triage_modal(n_samples, actions)` returning a `modalDialog`; server inputs `triage_rule_<key>` carrying `"off"` / `"flag"` / `"drop"`; server outputs `triage_n_<key>` and `triage_preview`.

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-triage-modal.R`:

```r
# The modal is the only place a curator sees what triage will do, so its
# contents are pinned: one control per rule, drop offered only where the
# evidence supports it, and a legend that says what the two counts mean.

test_that("the modal renders one control per rule, in registry order", {
  m <- RiboPilot:::triage_modal(3, RiboPilot:::triage_actions_default())
  html <- as.character(m)
  for (k in RiboPilot:::triage_rule_keys()) {
    expect_match(html, paste0("triage_rule_", k), fixed = TRUE)
  }
  pos <- vapply(RiboPilot:::triage_rule_keys(),
    function(k) regexpr(paste0("triage_rule_", k), html, fixed = TRUE)[1], numeric(1))
  expect_false(is.unsorted(pos))
})

test_that("drop is offered only on the combination rules", {
  # asserts against the RENDERED modal, not against triage_rule_control()'s own
  # inputs: comparing a helper's output to the argument it was given passes
  # whatever the modal actually shows
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  keys <- RiboPilot:::triage_rule_keys()
  pos <- vapply(keys, function(k)
    regexpr(paste0('id="triage_rule_', k, '"'), html, fixed = TRUE)[1], numeric(1))
  expect_false(any(pos < 0))
  ends <- c(pos[-1] - 1, nchar(html))
  for (i in seq_along(keys)) {
    seg <- substr(html, pos[i], ends[i])
    expect_true(grepl('value="off"', seg, fixed = TRUE), info = keys[i])
    expect_true(grepl('value="flag"', seg, fixed = TRUE), info = keys[i])
    expect_equal(grepl('value="drop"', seg, fixed = TRUE),
      "drop" %in% RiboPilot:::TRIAGE_RULE_DEFS[[i]]$actions, info = keys[i])
  }
})

test_that("the legend explains both counts without jargon", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_match(html, "first one that matches decides whether to drop or flag", fixed = TRUE)
  expect_match(html, "even if a higher-priority rule acts on it first", fixed = TRUE)
})

test_that("the modal no longer offers a depth control", {
  html <- as.character(RiboPilot:::triage_modal(1, RiboPilot:::triage_actions_default()))
  expect_false(grepl("triage_depth_frac", html, fixed = TRUE))
  expect_false(grepl("deepest scaffold", html, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "triage-modal")'`
Expected: FAIL, `unused argument (actions)` or `could not find function "triage_rule_control"`.

- [ ] **Step 3: Implement**

In `R/app_modals.R`, replace `triage_modal()`:

```r
#' One rule's three-way action control
#'
#' Off is the disable, so there is no separate tickbox. Drop is offered only on
#' rules that combine both evidence channels; on the others it is absent and the
#' label carries the reason, because a control that looks available and refuses
#' to act is worse than one that is not there.
#' @noRd
triage_rule_control <- function(def, current) {
  choices <- c("Off" = "off", "Flag" = "flag")
  if ("drop" %in% def$actions) choices <- c(choices, "Drop" = "drop")
  shinyWidgets::radioGroupButtons(
    inputId = paste0("triage_rule_", def$key),
    label = NULL, choices = choices, selected = current,
    size = "xs", individual = FALSE)
}

#' Bulk triage modal: re-run automatic triage over the selected samples
#'
#' A bulk write across a project is the highest-blast-radius action in the app,
#' so it always shows what it would do first, per rule, including how many
#' samples would end up with NOTHING kept. Never applied without the counts.
#' @param n_samples Number of samples the action would cover.
#' @param actions Named list key -> action, from [triage_opts_for()].
#' @noRd
triage_modal <- function(n_samples, actions = triage_actions_default()) {
  rows <- lapply(TRIAGE_RULE_DEFS, function(def) {
    tags$tr(
      tags$td(style = "padding:3px 8px 3px 0; vertical-align:top;",
        tags$div(style = "font-weight:600;", def$label),
        tags$div(class = "text-muted", style = "font-size:0.8em;", def$meaning)),
      tags$td(style = "padding:3px 8px; vertical-align:top; white-space:nowrap;",
        triage_rule_control(def, actions[[def$key]])),
      tags$td(style = "padding:3px 8px; vertical-align:top; text-align:right;",
        uiOutput(paste0("triage_n_", def$key), inline = TRUE)))
  })
  modalDialog(
    title = sprintf("Triage scaffolds - %d selected sample(s)", n_samples),
    size = "l", easyClose = TRUE,
    tags$p(class = "text-muted", style = "font-size:0.86em;",
      paste("Proposes keep / drop / undecided per scaffold from the evidence already",
        "in the database. Depth is shown in the scaffold table but never decides:",
        "high depth is often k-mer collapse across the conserved core or IGS",
        "subrepeats, and low depth is often isomer paths splitting reads.")),
    tags$table(style = "width:100%; font-size:0.9em;",
      tags$thead(tags$tr(
        tags$th(style = "text-align:left;", "Rule"),
        tags$th(style = "text-align:left;", "Action"),
        tags$th(style = "text-align:right;", "Decided / Matching"))),
      tags$tbody(rows)),
    tags$p(class = "text-muted", style = "font-size:0.8em; margin-top:6px;",
      paste("Rules are applied top to bottom and the first one that matches decides",
        "whether to drop or flag the scaffold. \"Decided\" counts scaffolds this rule",
        "will drop or flag. \"Matching\" counts every scaffold that meets the rule,",
        "even if a higher-priority rule acts on it first.")),
    div(style = "margin:0.4em 0;",
      shinyWidgets::prettyCheckbox("triage_overwrite_user",
        "Also re-decide scaffolds I set by hand", value = FALSE, status = "danger")),
    tags$hr(style = "margin:8px 0;"),
    uiOutput("triage_preview"),
    footer = tagList(
      actionButton("triage_apply", "Apply", class = "btn-danger btn-sm"),
      modalButton("Cancel"))
  )
}
```

In `R/app_server.R`, replace the bulk triage block (lines 1059-1101):

```r
  # ---- bulk triage over the selected samples ----
  triage_target_ids <- function() {
    ids <- unique(selected_ids())
    if (length(ids)) ids else unique(samples_data()$ID)
  }
  # settings come from the DB so the modal and the WF1 triage process agree; the
  # inputs override them live so the preview reflects what Apply would do.
  triage_actions <- reactive({
    a <- triage_actions_default()
    for (k in triage_rule_keys()) {
      v <- input[[paste0("triage_rule_", k)]]
      if (!is.null(v) && nzchar(v)) a[[k]] <- v
    }
    normalize_triage_actions(a)
  })
  observeEvent(input$triage_btn, {
    ids <- triage_target_ids()
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    cur <- if (length(ids)) triage_opts_for(con, ids[1]) else triage_actions_default()
    showModal(triage_modal(length(ids), cur))
  })
  triage_pv <- reactive({
    ids <- triage_target_ids()
    if (!length(ids)) return(NULL)
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    tryCatch(triage_preview(con, ids, opts = triage_actions(),
      overwrite_user = isTRUE(input$triage_overwrite_user)), error = function(e) NULL)
  })
  # one count output per rule, rendered into the static table so changing a
  # control never rebuilds (and resets) the controls themselves
  for (k in triage_rule_keys()) {
    local({
      key <- k
      output[[paste0("triage_n_", key)]] <- renderUI({
        p <- triage_pv()
        if (is.null(p) || is.null(p$by_rule)) return(span("-"))
        r <- p$by_rule[p$by_rule$key == key, ]
        if (nrow(r) == 0) return(span("-"))
        span(style = "white-space:nowrap;",
          tags$strong(r$n_decided[1]), " / ",
          span(style = "color:#777;", r$n_matching[1]))
      })
    })
  }
  output$triage_preview <- renderUI({
    p <- triage_pv()
    if (is.null(p)) return(em("No samples to triage."))
    no_drop <- !any(triage_actions() %in% "drop")
    tagList(
      tags$p(sprintf("Would drop %d of %d scaffold(s) across %d sample(s); %d kept, %d left undecided.",
        p$n_drop, p$n_units, p$n_samples, p$n_keep, p$n_undecided)),
      if (no_drop) tags$p(style = "color:#555;",
        "No rule is set to Drop, so Apply will only add and clear reasons.")
      else if (p$n_drop == 0) tags$p(style = "color:#555;",
        "A rule is set to Drop, but nothing in this selection will be dropped."),
      if (p$n_user_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d scaffold(s) you decided by hand will be left alone.", p$n_user_skipped)),
      if (p$n_locked_skipped > 0) tags$p(style = "color:#555;",
        sprintf("%d annotation-locked scaffold(s) will be left alone.", p$n_locked_skipped)),
      if (length(p$emptied) > 0) div(style = "color:#842029; font-weight:600; margin-top:6px;",
        sprintf("%d sample(s) would end up with NO scaffolds kept: %s",
          length(p$emptied), paste(utils::head(p$emptied, 12), collapse = ", ")))
      else div(style = "color:#2d6a4f; margin-top:6px;", "No sample would be left with nothing kept."))
  })
  observeEvent(input$triage_apply, {
    ids <- triage_target_ids(); req(length(ids))
    con <- app_con(); on.exit(DBI::dbDisconnect(con))
    act <- triage_actions()
    # persist first: a pipeline re-run must reproduce what the curator just saw
    write_triage_opts(con, act)
    n <- 0L
    for (id in ids) {
      n <- n + tryCatch(auto_triage_sample(con, id, opts = act,
        overwrite_user = isTRUE(input$triage_overwrite_user)), error = function(e) 0L)
    }
    removeModal()
    refresh(refresh() + 1)
    showNotification(sprintf("Triage updated %d scaffold(s) across %d sample(s).", n, length(ids)),
      type = "message")
  })
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "triage")'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/app_modals.R R/app_server.R tests/testthat/test-triage-modal.R
git commit -m "feat(app): per-rule triage controls with live decided/matching counts"
```

---

### Task 6: Depth as display-only evidence

**Files:**
- Modify: `R/coverage.R` (add `downsample_depths()`)
- Modify: `R/app_utils.R:101-120` (`unit_evidence()` selects `c.depths`, adds `depth_track`)
- Modify: `R/utils_reactable.R` (add `rt_depth_sparkline()`)
- Modify: `R/app_server.R:885-915` (scaffold table column defs)
- Test: `tests/testthat/test-coverage-sparkline.R` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces:
  - `downsample_depths(x, n = 100)` -> numeric vector of length `<= n`
  - `depth_track_string(depths_blob, n = 100)` -> comma-separated string, or `NA_character_`
  - `rt_depth_sparkline()` -> `htmlwidgets::JS` cell renderer drawing an inline SVG polyline

- [ ] **Step 1: Write the failing test**

Create `tests/testthat/test-coverage-sparkline.R`:

```r
# coverage.depths is a per-position blob, several kb per scaffold. It is
# downsampled server-side before it reaches a table cell.

test_that("downsampling preserves the ends and bounds the length", {
  x <- as.numeric(1:1000)
  y <- RiboPilot:::downsample_depths(x, n = 100)
  expect_lte(length(y), 100L)
  expect_equal(y[1], 1)
  expect_equal(y[length(y)], 1000)
})

test_that("a short track is returned unchanged", {
  x <- c(3, 4, 5)
  expect_equal(RiboPilot:::downsample_depths(x, n = 100), x)
})

test_that("an empty or unparseable blob yields NA rather than erroring", {
  expect_true(is.na(RiboPilot:::depth_track_string(NA_character_)))
  expect_true(is.na(RiboPilot:::depth_track_string("")))
  expect_true(is.na(RiboPilot:::depth_track_string("not,a,number")))
})

test_that("a comma-separated blob round-trips into a bounded track string", {
  blob <- paste(1:500, collapse = ",")
  s <- RiboPilot:::depth_track_string(blob, n = 50)
  parts <- strsplit(s, ",", fixed = TRUE)[[1]]
  expect_lte(length(parts), 50L)
  expect_equal(as.numeric(parts[1]), 1)
})

test_that("the sparkline renderer draws an SVG and never evaluates its input", {
  js <- as.character(RiboPilot:::rt_depth_sparkline())
  expect_match(js, "svg", fixed = TRUE)
  expect_match(js, "polyline", fixed = TRUE)
  expect_false(grepl("eval(", js, fixed = TRUE))
  expect_false(grepl("innerHTML", js, fixed = TRUE))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_local(filter = "sparkline")'`
Expected: FAIL, `could not find function "downsample_depths"`.

- [ ] **Step 3: Implement**

Append to `R/coverage.R`:

```r
#' Reduce a per-position depth track to a bounded number of points
#'
#' Takes the value at evenly spaced positions rather than a mean, so a dropout
#' stays visible instead of being averaged away.
#' @param x Numeric depths, in position order.
#' @param n Maximum points to keep.
#' @return Numeric vector of length at most `n`, first and last preserved.
#' @noRd
downsample_depths <- function(x, n = 100) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  if (length(x) <= n) return(x)
  x[unique(round(seq(1, length(x), length.out = n)))]
}

#' Downsampled depth track as a comma-separated string for a table cell
#' @param depths_blob The stored `coverage.depths` value.
#' @param n Maximum points to keep.
#' @return Comma-separated string, or NA when there is no usable track.
#' @noRd
depth_track_string <- function(depths_blob, n = 100) {
  if (length(depths_blob) == 0 || is.na(depths_blob[1]) || !nzchar(depths_blob[1])) {
    return(NA_character_)
  }
  v <- suppressWarnings(as.numeric(strsplit(depths_blob[1], "[,;[:space:]]+")[[1]]))
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA_character_)
  paste(round(downsample_depths(v, n), 1), collapse = ",")
}
```

In `R/app_utils.R`, add `c.depths AS depths_raw,` to `unit_evidence()`'s SELECT next to the depth columns, and after `d$unit <- ...` (line 125):

```r
  d$depth_track <- vapply(d$depths_raw %||% rep(NA_character_, nrow(d)),
    depth_track_string, character(1), USE.NAMES = FALSE)
  d$depths_raw <- NULL
```

Append to `R/utils_reactable.R`:

```r
#' Inline depth sparkline from a comma-separated track
#'
#' Shown because a curator deciding keep/drop should see coverage shape, and
#' deliberately not used by any rule: high depth is often k-mer collapse across
#' the conserved core or IGS subrepeats, and low depth is often isomer paths
#' splitting reads.
#' @noRd
rt_depth_sparkline <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var a = String(v).split(',').map(Number).filter(function(x){ return !isNaN(x); });
      if (a.length < 2) return '';
      var w = 90, h = 18, mx = Math.max.apply(null, a);
      if (!(mx > 0)) return '';
      var pts = a.map(function(y, i) {
        var px = (i / (a.length - 1)) * w;
        var py = h - (y / mx) * h;
        return px.toFixed(1) + ',' + py.toFixed(1);
      }).join(' ');
      return '<svg width=\"' + w + '\" height=\"' + h + '\" ' +
        'viewBox=\"0 0 ' + w + ' ' + h + '\" preserveAspectRatio=\"none\">' +
        '<polyline points=\"' + pts + '\" fill=\"none\" stroke=\"#6c757d\" ' +
        'stroke-width=\"1\"/></svg>';
    }")
}
```

`htmlwidgets::JS()` wraps the string exactly once. Do not add a trailing
`|> htmlwidgets::JS()` the way `rt_triage_bttn()` does after its `sprintf()`;
there is no `sprintf()` here and double-wrapping produces a nested `JS` object
that reactable renders as literal text.

In `R/app_server.R`, add `"depth_track"` to the `lead` vector at line 885 immediately after `"depth"`, and add its colDef next to the other scaffold columns:

```r
        depth_track = reactable::colDef(name = "Coverage", width = 100, show = TRUE,
          filterable = FALSE, sortable = FALSE, html = TRUE,
          cell = rt_depth_sparkline()),
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `Rscript -e 'testthat::test_local(filter = "sparkline")'`
Expected: PASS, 5 tests.

Then the whole suite: `Rscript -e 'testthat::test_local()'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/coverage.R R/app_utils.R R/utils_reactable.R R/app_server.R tests/testthat/test-coverage-sparkline.R
git commit -m "feat(app): show per-position coverage at the keep/drop decision"
```

---

### Task 7: Documentation and final verification

**Files:**
- Modify: `vignettes/Curation-Flags.Rmd` (triage section)
- Modify: `dev/design/05_ROADMAP.md`
- Modify: `man/` (regenerate)

- [ ] **Step 1: Update the vignette**

In `vignettes/Curation-Flags.Rmd`, find the triage description and replace it with the rule table from `dev/design/09_TRIAGE_RULES_AND_DEPTH.md` section 3, plus one paragraph stating that depth is displayed and never acted on, and why. Read the surrounding prose first and match its voice; do not paste the design doc verbatim.

- [ ] **Step 2: Note the change in the roadmap**

Add one line to `dev/design/05_ROADMAP.md` recording that triage rules became named and user-controllable, with a pointer to doc 09.

- [ ] **Step 3: Regenerate Rd and run the full suite**

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "rd")'
Rscript -e 'testthat::test_local()'
```

Expected: all tests PASS.

- [ ] **Step 4: Check the package**

```bash
Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'
```

Expected: 0 errors, 0 warnings. Investigate any new NOTE before proceeding.

- [ ] **Step 5: Verify the pipeline path end to end**

The WF1 `triage` process calls `auto_triage_sample(con, id)` with no `opts`, which now resolves from `assemble.triage_opts`. Confirm `inst/nextflow/bin/rp_triage.R` needs no edit, then lint:

```bash
nextflow lint inst/nextflow/main.nf
```

Expected: no errors.

- [ ] **Step 6: Commit**

```bash
git add vignettes/ dev/design/05_ROADMAP.md man/
git commit -m "docs(triage): document the named rules and why depth does not decide"
```

---

## Self-Review

**Spec coverage:**
- Section 1 decisions 1-3 -> Task 1. Decision 4 -> Tasks 4 and 6. Decision 5 -> Tasks 2, 3, 5.
- Section 2 channels and 2.1 annotation gate -> Task 1 (`triage_rule_ctx`) and Task 4 (`annotated_n`).
- Section 3 registry and ordering -> Task 1. Section 3.1 coercion -> Tasks 1 and 3. Section 3.2 reason strings -> Task 1.
- Section 4 schema v12 -> Task 2; `triage_opts_for` -> Task 3; `triage_rule` write -> Task 4.
- Section 5 modal -> Task 5; `by_rule` -> Task 4.
- Section 6 depth -> Task 4 (removal) and Task 6 (display).
- Section 7 behaviour change -> pinned by tests in Tasks 1 and 4.
- Section 8 tests -> distributed across Tasks 1, 3, 4, 5, 6.
- Section 9 out of scope -> nothing in this plan touches the identity envelope, the `screen` table, or `NUCLEAR_RRNA_TITLE_RE`.

**Placeholders:** none. Every code step carries the code.

**Type consistency:** rule keys are identical across `TRIAGE_RULE_DEFS`, `triage_opts` column names (`rule_` prefix), `assemblies.triage_rule` values, and modal input ids (`triage_rule_` prefix). `triage_rules()` keeps its `(d, opts)` signature so existing call sites read the same. `triage_preview()` and `auto_triage_sample()` both take `opts = NULL` meaning "resolve from the DB".

**Known rough edge, deliberately left to the implementer:** Task 5's `for (k in triage_rule_keys()) local({...})` registers outputs at server-construction time, which is correct because the rule set is a compile-time constant. If `TRIAGE_RULE_DEFS` ever becomes dynamic this must move to `observe()`.
