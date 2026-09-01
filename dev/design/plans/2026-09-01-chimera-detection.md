# Chimera Detection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a per-unit chimera detector (per-gene taxonomy, windowed local BLAST, ribotyper, read evidence) with lineage-first storage and a track under the depth plot.

**Architecture:** One new Nextflow process `chimera` after `reconcile`, gated on coverage and BLAST completion. Pure R functions in `R/chimera.R` compute four channels from DB rows plus two tool calls (blastn, ribotyper); a thin driver writes `chimera` and `chimera_windows`. Triage, export risks, the scaffold table, and the annotate modal read those tables.

**Tech Stack:** R (DBI/RSQLite, ggplot2, reactable, testthat), Nextflow DSL2, blastn, samtools, ribovore (bioconda), Docker/mamba.

**Spec:** `dev/design/21_CHIMERA_DETECTION.md`

## Global Constraints

- No non-ASCII characters in R code. No em dashes anywhere.
- Minimal comments. Do not explain fixes in comments.
- Pre-release: new tables are `CREATE TABLE IF NOT EXISTS` DDL constants in `R/db_schema.R`, executed by `new_db()`. No migration code.
- Never `git push`.
- Constants (thresholds) live in `R/chimera.R`: `CHIMERA_WIN = 500L`, `CHIMERA_STEP = 250L`, `CHIMERA_PID_DROP = 90`, `CHIMERA_DEPTH_STEP = 3`, `CHIMERA_CLIP_FRAC = 0.3`, `CHIMERA_JUNCTION_TOL = 300L`, `CHIMERA_TAXONOMY_CAP = 0.7`.
- Labels: `clean` (score < 0.3), `suspect` (0.3 to 0.7 inclusive), `likely` (> 0.7).
- Chimera calls warn. Nothing blocks export.
- Run tests with `Rscript -e 'devtools::test(filter = "chimera")'` from the repo root, or `testthat::test_file()` for one file.
- Commit after each task. Commit messages end with `Claude-Session: https://claude.ai/code/session_01BnCXUHMYRkhpQcaQhv6aPB`. No other attribution.

---

## File map

| File | Responsibility |
|---|---|
| `R/db_schema.R` | `CHIMERA_DDL`, `CHIMERA_WINDOWS_DDL` (new) |
| `R/init_db.R` | execute the two DDLs in `new_db()`; `coverage.clips` column |
| `R/coverage.R` | parse optional `clips.txt`, store `clips` track |
| `R/chimera.R` | all pure channel logic, scoring, reason text, tool wrappers, DB writes (new) |
| `inst/nextflow/bin/rp_chimera.R` | driver (new) |
| `inst/nextflow/bin/rp_coverage.R` | optional 4th arg |
| `inst/nextflow/modules/coverage.nf` | clip pass |
| `inst/nextflow/modules/pipeline_r.nf` | `chimera` process |
| `inst/nextflow/main.nf` | wiring and gate; triage waits on chimera |
| `docker/Dockerfile` | add `ribovore` |
| `R/triage_rules.R` | `chimeric` rule |
| `R/export.R` | `submission_risks()` reads `chimera` |
| `R/curate_rrna.R` | retire `chimeric` curation flag |
| `R/constants.R` | `SCAFFOLD_COL_HELP$chimera` |
| `R/utils_reactable.R` | `rt_chimera_badge()` |
| `R/app_server.R` | scaffold table column; `chimera_plot`, segment tables |
| `R/app_modals.R` | `plot_chimera()`, modal placement |
| `tests/testthat/test-chimera.R` | unit tests (new) |
| `tests/testthat/fixtures/chimera/` | ribotyper short output, synthetic scaffolds (new) |
| `tests/testthat/helper-app-fixture.R` | add the two DDLs |

---

### Task 1: Schema

**Files:**
- Modify: `R/db_schema.R` (append after `SCREEN_DDL`, ~line 102)
- Modify: `R/init_db.R:505-516` (coverage DDL), `R/init_db.R:635` (execute DDLs)
- Modify: `tests/testthat/helper-app-fixture.R:4-10`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces: tables `chimera`, `chimera_windows`, column `coverage.clips`. Constants `CHIMERA_DDL`, `CHIMERA_WINDOWS_DDL`.

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-chimera.R
test_that("new_db creates the chimera tables and coverage.clips", {
  db <- tempfile(fileext = ".sqlite")
  con <- RiboPilot:::new_db(db)
  on.exit(DBI::dbDisconnect(con))
  tabs <- DBI::dbListTables(con)
  expect_true(all(c("chimera", "chimera_windows") %in% tabs))
  expect_true("clips" %in% DBI::dbListFields(con, "coverage"))
  expect_true("segments" %in% DBI::dbListFields(con, "chimera"))
  expect_true("lineage" %in% DBI::dbListFields(con, "chimera_windows"))
})
```

Check how `new_db()` is called elsewhere in tests (grep `new_db(` in `tests/testthat`) and match its signature and return; if it returns a path rather than a connection, open the connection with `RiboPilot:::open_project_db(db)`.

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-chimera.R")'`
Expected: FAIL, tables missing.

- [ ] **Step 3: Add the DDL constants**

Append to `R/db_schema.R`:

```r
#' DDL for chimera: one row per unit, lineage kept as text at every level.
#' @noRd
CHIMERA_DDL <- "CREATE TABLE IF NOT EXISTS chimera (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  score REAL,
  label TEXT,
  junctions TEXT,
  reason TEXT,
  segments TEXT,
  gene_taxa TEXT,
  channels TEXT,
  ribotyper_codes TEXT,
  db_source TEXT,
  time_stamp INTEGER,
  PRIMARY KEY (ID, path, scaffold)
);"

#' DDL for chimera_windows: best local hit per sliding window.
#' @noRd
CHIMERA_WINDOWS_DDL <- "CREATE TABLE IF NOT EXISTS chimera_windows (
  ID TEXT NOT NULL,
  path INTEGER NOT NULL DEFAULT 1,
  scaffold INTEGER NOT NULL DEFAULT 1,
  start INTEGER NOT NULL,
  stop INTEGER,
  accession TEXT,
  species TEXT,
  staxid INTEGER,
  lineage TEXT,
  lineage_taxids TEXT,
  rank_distance INTEGER,
  hit_class TEXT,
  pident REAL,
  bitscore REAL,
  depth_median REAL,
  clip_frac REAL,
  db_source TEXT,
  PRIMARY KEY (ID, path, scaffold, start)
);"
```

In `R/init_db.R` coverage DDL add `clips TEXT,` after `depths TEXT,`. Next to line 635 (`DBI::dbExecute(con, SCREEN_DDL)`) add:

```r
  DBI::dbExecute(con, CHIMERA_DDL)
  DBI::dbExecute(con, CHIMERA_WINDOWS_DDL)
```

In `tests/testthat/helper-app-fixture.R` `fixture_fill_schema()`, add `RiboPilot:::CHIMERA_DDL, RiboPilot:::CHIMERA_WINDOWS_DDL` to the DDL list.

- [ ] **Step 4: Run test to verify it passes**

Run: `Rscript -e 'testthat::test_file("tests/testthat/test-chimera.R")'`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/db_schema.R R/init_db.R tests/testthat/helper-app-fixture.R tests/testthat/test-chimera.R
git commit -m "feat: chimera tables and coverage.clips"
```

---

### Task 2: Soft-clip track in coverage

**Files:**
- Modify: `inst/nextflow/modules/coverage.nf:46-60`
- Modify: `inst/nextflow/bin/rp_coverage.R`
- Modify: `R/coverage.R:59-82` (`write_coverage_units`)
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces: `write_coverage_units(con, id, depth_file, clips_file = NULL)`; `coverage.clips` comma track aligned to `depths` (same length, zeros where no clips).

- [ ] **Step 1: Write the failing test**

```r
test_that("write_coverage_units stores an aligned clip track", {
  db <- tempfile(fileext = ".sqlite")
  con <- RiboPilot:::new_db(db)
  on.exit(DBI::dbDisconnect(con))
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, 'ACGTACGT')")
  dep <- tempfile(); cl <- tempfile()
  writeLines(sprintf("1.1\t%d\t%d", 1:8, c(5, 5, 6, 6, 20, 20, 20, 20)), dep)
  writeLines(c("1.1\t5\t7", "1.1\t6\t1"), cl)
  RiboPilot:::write_coverage_units(con, "S1", dep, cl)
  r <- DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = 'S1'")
  expect_equal(r$depths, "5,5,6,6,20,20,20,20")
  expect_equal(r$clips, "0,0,0,0,7,1,0,0")
})
```

If `assemblies` has NOT NULL columns beyond these, add them to the INSERT (check `R/init_db.R:348`).

- [ ] **Step 2: Run test to verify it fails**

Expected: FAIL, unused argument `cl`.

- [ ] **Step 3: Implement**

In `R/coverage.R`, change the signature to `write_coverage_units <- function(con, id, depth_file, clips_file = NULL)`. After the per-unit depth vectors are built (where `depths` is pasted with commas), add:

```r
  clips <- NULL
  if (!is.null(clips_file) && file.exists(clips_file) && file.size(clips_file) > 0) {
    clips <- utils::read.table(clips_file, sep = "\t", col.names = c("ref", "pos", "n"),
      colClasses = c("character", "integer", "integer"), quote = "", comment.char = "")
  }
  clip_track <- function(ref, len) {
    v <- integer(len)
    if (!is.null(clips)) {
      c1 <- clips[clips$ref == ref & clips$pos >= 1 & clips$pos <= len, , drop = FALSE]
      if (nrow(c1)) v[c1$pos] <- v[c1$pos] + c1$n
    }
    paste(v, collapse = ",")
  }
```

and write `clips = clip_track(ref, length(dep))` alongside `depths` in the row for each unit. Follow whatever loop variable names the existing body uses for the ref key and the depth vector. Include `clips` in the INSERT column list.

In `inst/nextflow/bin/rp_coverage.R`, change usage to `rp_coverage.R <db> <id> <depth_file> [clips_file]` (still `rp_args(a, 3, ...)`) and call `write_coverage_units(con, a[2], a[3], if (length(a) >= 4) a[4] else NULL)`.

In `coverage.nf`, before the loop add `: > clips.txt`, and after `samtools depth -aa aln.bam >> depth.txt` add:

```bash
        # per-position count of reads soft-clipped >= 20 bp at that end
        samtools view aln.bam | awk -v OFS='\t' '{
            cig = $6; pos = $4; ref = $3
            if (match(cig, /^([0-9]+)S/, m) && m[1] + 0 >= 20) print ref, pos
            if (match(cig, /([0-9]+)S$/, m) && m[1] + 0 >= 20) {
                n = 0; s = cig
                while (match(s, /^[0-9]+[MDN=X]/)) { n += substr(s, 1, RLENGTH - 1); s = substr(s, RLENGTH + 1) }
                while (match(s, /[0-9]+[MDN=X]/)) { n += substr(s, RSTART, RLENGTH - 1); s = substr(s, RSTART + RLENGTH) }
                print ref, pos + n - 1
            }
        }' | sort | uniq -c | awk -v OFS='\t' '{print $2, $3, $1}' >> clips.txt
```

`match(..., m)` is gawk. Check the image has gawk (`docker/Dockerfile` base is mambaforge on Debian, which has gawk as `awk`? If `awk --version` is mawk, add `'gawk'` to the mamba install line in Task 5 and call `gawk` explicitly here). Then change the driver call to `rp_coverage.R '${params.db}' '${id}' depth.txt clips.txt`.

- [ ] **Step 4: Run test to verify it passes**

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/coverage.R inst/nextflow/bin/rp_coverage.R inst/nextflow/modules/coverage.nf tests/testthat/test-chimera.R
git commit -m "feat: soft-clip track alongside depth"
```

---

### Task 3: Window tiling and per-gene taxonomy (channels 1 and 2 data prep)

**Files:**
- Create: `R/chimera.R`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces:
  - `chimera_tile(seq_len, win = CHIMERA_WIN, step = CHIMERA_STEP)` returns data.frame(start, stop), 1-based inclusive, last window right-anchored, a scaffold shorter than `win` yields one window.
  - `chimera_windows_fasta(seq, tiles, file)` writes windows named `w<start>_<stop>`.
  - `chimera_gene_lca(hits, feats, min_cov = 0.5)` returns data.frame(element, taxon, rank, taxid, species, pident, n_hits) with `unresolved` rows (taxon NA) for elements with fewer than 3 covering accessions.
  - `chimera_channel_genes(gene_taxa)` returns `list(score, junctions, reason, availability)`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("tiling anchors the last window and handles short scaffolds", {
  t <- RiboPilot:::chimera_tile(1200, win = 500, step = 250)
  expect_equal(t$start, c(1, 251, 501, 701))
  expect_equal(t$stop, c(500, 750, 1000, 1200))
  expect_equal(RiboPilot:::chimera_tile(300, 500, 250), data.frame(start = 1L, stop = 300L))
})

test_that("per-gene LCA uses only HSPs covering half the element", {
  feats <- data.frame(element = c("18S", "28S"), start = c(1L, 3000L), stop = c(1800L, 6000L))
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  hits <- data.frame(
    accession = c("a1", "a2", "a3", "f1", "f2", "f3", "f4"),
    bitscore = c(3000, 2900, 2800, 5000, 4900, 4800, 100),
    qstart = c(1, 5, 10, 3000, 3100, 3050, 1500),
    qend = c(1800, 1790, 1750, 6000, 5900, 5950, 1600),
    lineage_taxids = c(fish, fish, fish, fung, fung, fung, fung),
    lineage = c(rep("Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar", 3),
                rep("Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium sp.", 4)),
    species = c(rep("Salmo salar", 3), rep("Penicillium sp.", 4)),
    pident = c(99, 99, 98, 97, 97, 96, 80),
    stringsAsFactors = FALSE)
  g <- RiboPilot:::chimera_gene_lca(hits, feats)
  expect_equal(g$taxon[g$element == "18S"], "Salmo salar")
  expect_equal(g$taxon[g$element == "28S"], "Penicillium sp.")
  ch <- RiboPilot:::chimera_channel_genes(g)
  expect_equal(ch$score, 0.6)
  expect_equal(ch$junctions, 2400L)
  expect_match(ch$reason, "18S.*Salmo salar.*28S.*Penicillium")
})

test_that("an element with too few covering hits is unresolved and does not vote", {
  feats <- data.frame(element = c("18S", "28S"), start = c(1L, 3000L), stop = c(1800L, 6000L))
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  hits <- data.frame(accession = c("a1", "a2", "a3"), bitscore = c(3, 2, 1),
    qstart = c(1, 5, 10), qend = c(1800, 1790, 1750),
    lineage_taxids = rep(fish, 3), lineage = rep("x;y;z;o;f;g;s", 3),
    species = rep("Salmo salar", 3), pident = c(99, 99, 98), stringsAsFactors = FALSE)
  g <- RiboPilot:::chimera_gene_lca(hits, feats)
  expect_true(is.na(g$taxon[g$element == "28S"]))
  ch <- RiboPilot:::chimera_channel_genes(g)
  expect_equal(ch$score, 0)
  expect_equal(ch$availability, "uninformative")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Expected: FAIL, functions not found.

- [ ] **Step 3: Implement**

Create `R/chimera.R`:

```r
# Chimera detection: per-gene taxonomy, windowed local BLAST, ribotyper,
# read evidence. Design: dev/design/21_CHIMERA_DETECTION.md.

CHIMERA_WIN <- 500L
CHIMERA_STEP <- 250L
CHIMERA_PID_DROP <- 90
CHIMERA_DEPTH_STEP <- 3
CHIMERA_CLIP_FRAC <- 0.3
CHIMERA_JUNCTION_TOL <- 300L
CHIMERA_TAXONOMY_CAP <- 0.7
CHIMERA_LABELS <- c("clean", "suspect", "likely")

#' @noRd
chimera_tile <- function(seq_len, win = CHIMERA_WIN, step = CHIMERA_STEP) {
  seq_len <- as.integer(seq_len)
  if (seq_len <= win) return(data.frame(start = 1L, stop = seq_len))
  starts <- seq.int(1L, seq_len - win + 1L, by = step)
  if (tail(starts, 1) + win - 1L < seq_len) starts <- c(starts, seq_len - win + 1L)
  data.frame(start = as.integer(starts), stop = as.integer(starts + win - 1L))
}

#' @noRd
chimera_windows_fasta <- function(seq, tiles, file) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  subs <- substring(seq, tiles$start, tiles$stop)
  writeLines(as.vector(rbind(paste0(">", ids), subs)), file)
  invisible(ids)
}

#' Bitscore-weighted LCA per annotated element from stored HSPs
#' @noRd
chimera_gene_lca <- function(hits, feats, min_cov = 0.5) {
  out <- lapply(seq_len(nrow(feats)), function(i) {
    s <- feats$start[i]; e <- feats$stop[i]; len <- e - s + 1
    ov <- pmin(hits$qend, e) - pmax(hits$qstart, s) + 1
    h <- hits[!is.na(ov) & ov >= min_cov * len & !is.na(hits$lineage_taxids), , drop = FALSE]
    h <- h[order(-h$bitscore), , drop = FALSE]
    h <- h[!duplicated(h$accession), , drop = FALSE]
    base <- data.frame(element = feats$element[i], taxon = NA_character_, rank = NA_character_,
      taxid = NA_integer_, species = NA_character_, pident = NA_real_, n_hits = nrow(h),
      stringsAsFactors = FALSE)
    if (nrow(h) < 3) return(base)
    lca <- bitscore_weighted_lca(h$lineage_taxids, h$bitscore, h$lineage)
    if (is.null(lca) || is.na(lca$taxid %||% NA)) return(base)
    base$taxon <- lca$taxon; base$rank <- lca$rank; base$taxid <- as.integer(lca$taxid)
    base$species <- h$species[1]; base$pident <- mean(h$pident, na.rm = TRUE)
    base
  })
  do.call(rbind, out)
}

#' @noRd
chimera_unavailable <- function(availability = "unavailable") {
  list(score = 0, junctions = integer(0), reason = NA_character_, availability = availability)
}

```

The class comparison needs each element's lineage taxids. In `chimera_gene_lca()`, add a `lineage_taxids` column to `base` (NA by default) and, once the LCA resolves, set it to `h$lineage_taxids[1]` with ranks deeper than `lca$rank` blanked: split with `strsplit(x, ";")[[1]]`, find `k <- match(lca$rank, RRNA_LINEAGE_RANKS)`, set elements after `k` to `""`, re-join with `;`. Then compare pairs with `rank_distance()` from `R/taxonomy_ranks.R`:

```r
#' Channel 1: two resolved elements whose LCAs sit in different classes
#' @noRd
chimera_channel_genes <- function(gene_taxa, feats = NULL) {
  g <- gene_taxa[!is.na(gene_taxa$taxon), , drop = FALSE]
  if (nrow(g) < 2) return(chimera_unavailable("uninformative"))
  worst <- NULL
  for (i in seq_len(nrow(g) - 1)) for (j in (i + 1):nrow(g)) {
    d <- rank_distance(g$lineage_taxids[i], g$lineage_taxids[j])
    if (!is.na(d) && d >= 5 && (is.null(worst) || d > worst$d)) worst <- list(i = i, j = j, d = d)
  }
  if (is.null(worst)) return(list(score = 0, junctions = integer(0),
    reason = sprintf("elements agree: %s", paste(sprintf("%s %s", g$element, g$taxon), collapse = ", ")),
    availability = "available"))
  a <- g[worst$i, ]; b <- g[worst$j, ]
  junction <- NA_integer_
  if (!is.null(feats)) {
    fa <- feats[feats$element == a$element, ]; fb <- feats[feats$element == b$element, ]
    if (nrow(fa) && nrow(fb)) junction <- as.integer(round((min(fa$stop, fb$stop) + max(fa$start, fb$start)) / 2))
  }
  list(score = 0.6, junctions = junction[!is.na(junction)],
    reason = sprintf("%s resolves to %s (%s, %.1f%%) but %s to %s (%s, %.1f%%)",
      a$element, a$taxon, a$species, a$pident, b$element, b$taxon, b$species, b$pident),
    availability = "available")
}
```

Update the second test to call `chimera_channel_genes(g, feats)` and expect `junctions == 2400L` (midpoint of 1800 and 3000). `rank_distance()` returns 0..7 where 7 means same species; confirm its direction in `R/taxonomy_ranks.R:101-106` and, if it counts shared ranks rather than differing ranks, invert the `d >= 5` test accordingly so that "differ at class or above" is what fires. Write a one-line test pinning that direction with the fish and fungus strings.

- [ ] **Step 4: Run tests to verify they pass**

Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add R/chimera.R tests/testthat/test-chimera.R
git commit -m "feat: chimera tiling and per-gene taxonomy channel"
```

---

### Task 4: Window segments (channel 2)

**Files:**
- Modify: `R/chimera.R`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Consumes: `parse_blast_hits_file()` frame (`qseqid` column named `w<start>_<stop>`), `tk_ranked_lineage()`, `rank_distance()`, `classify_hit_title()`.
- Produces:
  - `chimera_window_hits(hits, tiles, sample_taxids, taxdump_dir, db_source)` returns one row per tile: start, stop, accession, species, staxid, lineage, lineage_taxids, rank_distance, hit_class, pident, bitscore, db_source (NA fields when no hit).
  - `chimera_segments(windows, deciding = "class")` returns data.frame(start, stop, taxon, rank, species, pident, n_windows, db_source, lineage_taxids).
  - `chimera_channel_windows(segments, identity_only = FALSE)` returns channel list.

- [ ] **Step 1: Write the failing tests**

```r
chim_win <- function(starts, lineage, lineage_taxids, pident = 99, species = "sp") {
  data.frame(start = starts, stop = starts + 499L, accession = paste0("a", starts),
    species = species, staxid = 1L, lineage = lineage, lineage_taxids = lineage_taxids,
    rank_distance = NA_integer_, hit_class = "nuclear_rrna", pident = pident,
    bitscore = 500, db_source = "panel", stringsAsFactors = FALSE)
}

test_that("segments collapse adjacent windows sharing a class", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  w <- rbind(
    chim_win(seq(1, 2001, 250), "A;Chordata;Actinopterygii;o;f;g;Salmo salar", fish, species = "Salmo salar"),
    chim_win(seq(2251, 4001, 250), "F;Ascomycota;Eurotiomycetes;o;f;g;Penicillium sp.", fung, 96, "Penicillium sp."))
  s <- RiboPilot:::chimera_segments(w)
  expect_equal(nrow(s), 2)
  expect_equal(s$taxon, c("Actinopterygii", "Eurotiomycetes"))
  expect_equal(s$start, c(1L, 2251L))
  expect_equal(s$stop, c(2500L, 4500L))
  ch <- RiboPilot:::chimera_channel_windows(s)
  expect_equal(ch$score, 0.5)
  expect_equal(ch$junctions, 2375L)
  expect_match(ch$reason, "Actinopterygii")
  expect_match(ch$reason, "Penicillium")
})

test_that("a sliver segment cannot fire the window channel", {
  fish <- "33208;7711;186623;8030;8016;8017;8018"
  fung <- "4751;4890;147545;5042;5043;5044;5045"
  w <- rbind(
    chim_win(seq(1, 5001, 250), "A;Chordata;Actinopterygii;o;f;g;s", fish),
    chim_win(seq(5251, 5501, 250), "F;Ascomycota;Eurotiomycetes;o;f;g;s", fung))
  ch <- RiboPilot:::chimera_channel_windows(RiboPilot:::chimera_segments(w))
  expect_equal(ch$score, 0)
})

test_that("identity-only mode reports junctions but scores zero", {
  w <- chim_win(seq(1, 3001, 250), NA_character_, NA_character_,
    pident = c(rep(99, 6), rep(80, 7)))
  s <- RiboPilot:::chimera_segments(w)
  ch <- RiboPilot:::chimera_channel_windows(s, identity_only = TRUE)
  expect_equal(ch$score, 0)
  expect_length(ch$junctions, 1)
})
```

- [ ] **Step 2: Run tests to verify they fail**

- [ ] **Step 3: Implement**

Append to `R/chimera.R`:

```r
#' Best hit per tile from a parse_blast_hits_file() frame
#' @noRd
chimera_window_hits <- function(hits, tiles, sample_taxids = NA_character_,
                                taxdump_dir = NULL, db_source = NA_character_) {
  ids <- sprintf("w%d_%d", tiles$start, tiles$stop)
  best <- if (nrow(hits)) hits[order(hits$qseqid, -hits$bitscore), , drop = FALSE] else hits
  best <- best[!duplicated(best$qseqid), , drop = FALSE]
  i <- match(ids, best$qseqid)
  lin <- tryCatch(tk_ranked_lineage(best$staxid, taxdump_dir), error = function(e) data.frame())
  li <- if (nrow(lin)) match(best$staxid[i], lin$taxid) else rep(NA_integer_, length(i))
  lt <- if (nrow(lin)) lin$lineage_taxids[li] else rep(NA_character_, length(i))
  rd <- if (!is.na(sample_taxids) && nzchar(sample_taxids)) {
    vapply(lt, function(x) if (is.na(x)) NA_integer_ else as.integer(rank_distance(x, sample_taxids)),
      integer(1), USE.NAMES = FALSE)
  } else rep(NA_integer_, length(i))
  data.frame(start = tiles$start, stop = tiles$stop,
    accession = best$accession[i], species = extract_species(best$title[i]),
    staxid = best$staxid[i], lineage = if (nrow(lin)) lin$lineage[li] else NA_character_,
    lineage_taxids = lt, rank_distance = rd,
    hit_class = classify_hit_title(best$title[i]),
    pident = best$pident[i], bitscore = best$bitscore[i], db_source = db_source,
    stringsAsFactors = FALSE)
}

#' Name at a rank from a ';'-joined 7-field lineage
#' @noRd
chimera_rank_name <- function(lineage, rank = "class") {
  k <- match(rank, RRNA_LINEAGE_RANKS)
  vapply(strsplit(lineage %||% NA_character_, ";", fixed = TRUE), function(p) {
    if (length(p) >= k && nzchar(p[k])) p[k] else NA_character_
  }, character(1))
}

#' Collapse windows into runs sharing a name at the deciding rank
#' @noRd
chimera_segments <- function(windows, deciding = "class") {
  if (!nrow(windows)) return(data.frame())
  key <- chimera_rank_name(windows$lineage, deciding)
  key[is.na(key)] <- ifelse(is.na(windows$pident[is.na(key)]), "nohit",
    ifelse(windows$pident[is.na(key)] >= CHIMERA_PID_DROP, "high_id", "low_id"))
  run <- cumsum(c(TRUE, key[-1] != key[-length(key)]))
  do.call(rbind, lapply(split(seq_along(run), run), function(ix) {
    w <- windows[ix, , drop = FALSE]
    top <- ix[which.max(w$bitscore %||% rep(0, nrow(w)))]
    data.frame(start = w$start[1], stop = w$stop[nrow(w)],
      taxon = key[ix[1]], rank = deciding,
      species = windows$species[top] %||% NA_character_,
      pident = mean(w$pident, na.rm = TRUE), n_windows = nrow(w),
      db_source = w$db_source[1], lineage_taxids = windows$lineage_taxids[top] %||% NA_character_,
      stringsAsFactors = FALSE)
  }))
}

#' Channel 2: adjacent segments of comparable length in different classes
#' @noRd
chimera_channel_windows <- function(segments, identity_only = FALSE) {
  if (!nrow(segments)) return(chimera_unavailable())
  s <- segments[segments$n_windows >= 2, , drop = FALSE]
  if (nrow(s) < 2) return(list(score = 0, junctions = integer(0),
    reason = sprintf("one lineage along the scaffold: %s", segments$taxon[1]), availability = "available"))
  best <- NULL
  for (i in seq_len(nrow(s) - 1)) {
    a <- s[i, ]; b <- s[i + 1, ]
    if (identical(a$taxon, b$taxon)) next
    if (a$taxon %in% c("nohit") || b$taxon %in% c("nohit")) next
    la <- a$stop - a$start + 1; lb <- b$stop - b$start + 1
    if (min(la, lb) < 0.25 * max(la, lb)) next
    if (!identity_only && rank_distance(a$lineage_taxids, b$lineage_taxids) < 5) next
    best <- list(a = a, b = b, junction = as.integer(round((a$stop + b$start) / 2)))
  }
  if (is.null(best)) return(list(score = 0, junctions = integer(0),
    reason = "no class-level switch between comparable segments", availability = "available"))
  desc <- function(x) sprintf("%d-%d %s (%s, %.1f%%)", x$start, x$stop, x$taxon, x$species %||% "", x$pident)
  list(score = if (identity_only) 0 else 0.5, junctions = best$junction,
    reason = paste(desc(best$a), desc(best$b), sep = "; "),
    availability = if (identity_only) "identity_only" else "available")
}
```

The identity-only test builds windows with NA lineage, so the `identity_only` branch must skip `rank_distance()` (it does, by the `!identity_only &&` guard). The `rank_distance` direction check from Task 3 applies here too.

- [ ] **Step 4: Run tests to verify they pass**

- [ ] **Step 5: Commit**

```bash
git add R/chimera.R tests/testthat/test-chimera.R
git commit -m "feat: chimera window segments channel"
```

---

### Task 5: Ribotyper parsing and read evidence (channels 3 and 4)

**Files:**
- Modify: `R/chimera.R`
- Create: `tests/testthat/fixtures/chimera/ribotyper.short.out`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces:
  - `chimera_parse_ribotyper(short_out_file)` returns data.frame(target, pass, codes) where `codes` is a character vector column (list) of failure names without the `*` and the parenthetical.
  - `chimera_channel_ribotyper(rt, spans, itsx_chimeric = FALSE)` where `spans` is data.frame(target, start) mapping target name to scaffold offset. Returns channel list; `codes` joined by `;` in `reason`.
  - `chimera_depth_steps(depths, feats, step = CHIMERA_STEP)` returns integer junction positions.
  - `chimera_clip_pileups(depths, clips, min_depth = 10)` returns integer positions.
  - `chimera_channel_reads(depths, clips, feats)` returns channel list.

- [ ] **Step 1: Write the fixture and failing tests**

`tests/testthat/fixtures/chimera/ribotyper.short.out`:

```
#idx  target  classification         strnd  p/f  unexpected_features
#---  ------  ---------------------  -----  ---  -------------------
1     u18S    SSU.Eukarya            plus   PASS  -
2     u28S    LSU.Eukarya            plus   FAIL  *InconsistentHits:(2);*LowScoreDifference:(0.02<0.10);
```

Tests:

```r
test_that("ribotyper short output is parsed into codes", {
  rt <- RiboPilot:::chimera_parse_ribotyper(test_path("fixtures/chimera/ribotyper.short.out"))
  expect_equal(rt$target, c("u18S", "u28S"))
  expect_equal(rt$pass, c(TRUE, FALSE))
  expect_equal(rt$codes[[2]], c("InconsistentHits", "LowScoreDifference"))
  ch <- RiboPilot:::chimera_channel_ribotyper(rt, data.frame(target = c("u18S", "u28S"), start = c(1L, 3000L)))
  expect_equal(ch$score, 0.8)
  expect_match(ch$reason, "28S.*InconsistentHits")
})

test_that("itsx chimeric alone is a weak ribotyper-channel signal", {
  ch <- RiboPilot:::chimera_channel_ribotyper(NULL, NULL, itsx_chimeric = TRUE)
  expect_equal(ch$score, 0.2)
})

test_that("depth steps ignore feature edges", {
  dep <- c(rep(10, 1000), rep(50, 1000))
  feats <- data.frame(element = "18S", start = 1L, stop = 1000L)
  expect_length(RiboPilot:::chimera_depth_steps(dep, feats), 0)
  expect_equal(RiboPilot:::chimera_depth_steps(dep, data.frame(element = character(0), start = integer(0), stop = integer(0))), 1000L)
})

test_that("clip pileups need depth and fraction", {
  dep <- rep(20, 500); cl <- integer(500); cl[300] <- 8; cl[10] <- 8
  expect_equal(RiboPilot:::chimera_clip_pileups(dep, cl), 300L)
  ch <- RiboPilot:::chimera_channel_reads(dep, cl, data.frame(element = character(0), start = integer(0), stop = integer(0)))
  expect_equal(ch$score, 0.7)
  expect_equal(ch$junctions, 300L)
})
```

- [ ] **Step 2: Run tests to verify they fail**

- [ ] **Step 3: Implement**

Append to `R/chimera.R`:

```r
CHIMERA_RIBOTYPER_SCORES <- c(InconsistentHits = 0.8, DuplicateRegion = 0.8, BothStrands = 0.8,
  MultipleFamilies = 0.6, LowScoreDifference = 0.2)

#' @noRd
chimera_parse_ribotyper <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  ln <- readLines(file, warn = FALSE)
  ln <- ln[!grepl("^#", ln) & nzchar(trimws(ln))]
  if (!length(ln)) return(data.frame(target = character(0), pass = logical(0), codes = I(list())))
  f <- strsplit(trimws(ln), "\\s+")
  target <- vapply(f, `[`, character(1), 2)
  pass <- vapply(f, `[`, character(1), 5) == "PASS"
  feat <- vapply(f, function(x) if (length(x) >= 6) x[6] else "-", character(1))
  codes <- lapply(feat, function(x) {
    if (x == "-") return(character(0))
    p <- strsplit(x, ";", fixed = TRUE)[[1]]
    p <- p[grepl("^\\*", p)]
    sub("^\\*([A-Za-z]+).*$", "\\1", p)
  })
  data.frame(target = target, pass = pass, codes = I(codes), stringsAsFactors = FALSE)
}

#' Channel 3: ribotyper structural failures, ITSx chimeric as a weak add-on
#' @noRd
chimera_channel_ribotyper <- function(rt, spans, itsx_chimeric = FALSE) {
  score <- if (isTRUE(itsx_chimeric)) 0.2 else 0
  parts <- if (isTRUE(itsx_chimeric)) "ITSx chimeric" else character(0)
  if (is.null(rt) || !nrow(rt)) {
    if (score == 0) return(chimera_unavailable())
    return(list(score = score, junctions = integer(0), reason = parts, availability = "available"))
  }
  for (i in seq_len(nrow(rt))) {
    cd <- rt$codes[[i]]; cd <- cd[cd %in% names(CHIMERA_RIBOTYPER_SCORES)]
    if (!length(cd)) next
    score <- max(score, CHIMERA_RIBOTYPER_SCORES[cd])
    parts <- c(parts, sprintf("%s %s", sub("^u", "", rt$target[i]), paste(cd, collapse = ";")))
  }
  list(score = unname(score), junctions = integer(0),
    reason = if (length(parts)) paste(parts, collapse = "; ") else "ribotyper passed",
    availability = "available")
}

#' @noRd
chimera_depth_steps <- function(depths, feats, step = CHIMERA_STEP, ratio = CHIMERA_DEPTH_STEP, edge = 100L) {
  n <- length(depths)
  if (n < 2 * step) return(integer(0))
  b <- seq.int(step, n - step, by = step)
  med <- function(a, z) stats::median(depths[a:z], na.rm = TRUE)
  out <- integer(0)
  for (p in b) {
    l <- med(max(1, p - step + 1), p); r <- med(p + 1, min(n, p + step))
    if (is.na(l) || is.na(r) || min(l, r) <= 0) next
    if (max(l, r) / min(l, r) < ratio) next
    edges <- c(feats$start, feats$stop)
    if (length(edges) && any(abs(edges - p) <= edge)) next
    out <- c(out, as.integer(p))
  }
  out
}

#' @noRd
chimera_clip_pileups <- function(depths, clips, min_depth = 10, frac = CHIMERA_CLIP_FRAC, end = 50L) {
  n <- min(length(depths), length(clips))
  if (n <= 2 * end) return(integer(0))
  i <- (end + 1):(n - end)
  ok <- depths[i] >= min_depth & clips[i] / pmax(depths[i], 1) >= frac
  as.integer(i[which(ok)])
}

#' Channel 4: depth steps and clip pileups
#' @noRd
chimera_channel_reads <- function(depths, clips, feats) {
  if (is.null(depths) || !length(depths)) return(chimera_unavailable())
  steps <- chimera_depth_steps(depths, feats)
  pile <- if (!is.null(clips) && length(clips)) chimera_clip_pileups(depths, clips) else integer(0)
  pile <- chimera_merge_positions(pile, 50L)
  score <- max(if (length(steps)) 0.4 else 0, if (length(pile)) 0.7 else 0)
  parts <- c(if (length(steps)) sprintf("depth step at %s", paste(steps, collapse = ",")),
             if (length(pile)) sprintf("clip pileup at %s", paste(pile, collapse = ",")))
  list(score = score, junctions = sort(unique(c(steps, pile))),
    reason = if (length(parts)) paste(parts, collapse = "; ") else "no read-evidence breakpoints",
    availability = if (is.null(clips) || !length(clips)) "depth_only" else "available")
}

#' Collapse positions closer than `gap` into their first member
#' @noRd
chimera_merge_positions <- function(x, gap) {
  if (!length(x)) return(integer(0))
  x <- sort(unique(x)); keep <- c(TRUE, diff(x) > gap); x[keep]
}
```

- [ ] **Step 4: Run tests to verify they pass**

- [ ] **Step 5: Commit**

```bash
git add R/chimera.R tests/testthat/test-chimera.R tests/testthat/fixtures/chimera/ribotyper.short.out
git commit -m "feat: chimera ribotyper and read-evidence channels"
```

---

### Task 6: Scoring, reason, and the unit runner

**Files:**
- Modify: `R/chimera.R`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces:
  - `chimera_score(channels)` where `channels` is a named list (`genes`, `windows`, `ribotyper`, `reads`) of channel lists. Returns list(score, label, junctions).
  - `chimera_reason(segments, channels, junctions)` returns one string.
  - `chimera_unit(con, id, path, scaffold, dbarg = NULL, taxdump_dir = NULL, blastn = "blastn", ribotyper = "ribotyper", workdir = tempfile())` runs everything and writes both tables. Returns the `chimera` row as a list.

- [ ] **Step 1: Write the failing tests**

```r
chim_ch <- function(score, junctions = integer(0), availability = "available") {
  list(score = score, junctions = junctions, reason = "r", availability = availability)
}

test_that("taxonomy-only evidence caps at suspect", {
  s <- RiboPilot:::chimera_score(list(genes = chim_ch(0.6, 2400L), windows = chim_ch(0.5, 2375L),
    ribotyper = chim_ch(0, availability = "unavailable"), reads = chim_ch(0, availability = "unavailable")))
  expect_equal(s$score, 0.7)
  expect_equal(s$label, "suspect")
})

test_that("corroborated junctions add a bonus and read evidence can reach likely", {
  s <- RiboPilot:::chimera_score(list(genes = chim_ch(0.6, 2400L), windows = chim_ch(0.5, 2375L),
    ribotyper = chim_ch(0), reads = chim_ch(0.7, 2410L)))
  expect_equal(s$score, 0.95)
  expect_equal(s$label, "likely")
  expect_equal(s$junctions, 2375L)
})

test_that("labels follow the thresholds", {
  expect_equal(RiboPilot:::chimera_label(0.29), "clean")
  expect_equal(RiboPilot:::chimera_label(0.3), "suspect")
  expect_equal(RiboPilot:::chimera_label(0.7), "suspect")
  expect_equal(RiboPilot:::chimera_label(0.71), "likely")
})

test_that("chimera_unit writes lineage-bearing rows with mocked tools", {
  db <- tempfile(fileext = ".sqlite")
  con <- RiboPilot:::new_db(db)
  on.exit(DBI::dbDisconnect(con))
  seq <- paste(rep("ACGT", 1000), collapse = "")
  DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('S1', 1, 1, ?)", params = list(seq))
  DBI::dbExecute(con, "INSERT INTO samples (ID, lineage_taxids) VALUES ('S1', '33208;7711;186623;8030;8016;8017;8018')")
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES ('S1',1,1,'18S',1,1800), ('S1',1,1,'28S',2500,4000)")
  DBI::dbExecute(con, "INSERT INTO coverage (ID, path, scaffold, depths, clips) VALUES ('S1',1,1,?,?)",
    params = list(paste(rep(20, 4000), collapse = ","), paste(c(rep(0, 2199), 9, rep(0, 1800)), collapse = ",")))
  fake_blast <- function(query, dbarg, task, out, ...) {
    ids <- sub("^>", "", grep("^>", readLines(query), value = TRUE))
    lines <- vapply(ids, function(w) {
      st <- as.integer(sub("^w(\\d+)_.*", "\\1", w))
      if (st < 2200) sprintf("NC_1\tSalmo salar 18S\t99.0\t100\t0\t8030\t900\t500\t1\t500\t1\t500\t1800\t%s", w)
      else sprintf("KX_1\tPenicillium sp. 28S\t97.0\t100\t0\t5073\t880\t500\t1\t500\t1\t500\t4000\t%s", w)
    }, character(1))
    writeLines(lines, out); 0L
  }
  fake_rt <- function(fasta, outdir, ...) {
    writeLines(c("#idx target class strnd p/f unexpected_features", "1 u18S SSU.Eukarya plus PASS -",
      "2 u28S LSU.Eukarya plus PASS -"), file.path(outdir, "rt.ribotyper.short.out")); 0L
  }
  fake_lin <- function(taxids, dir) data.frame(taxid = c(8030L, 5073L),
    lineage = c("Animalia;Chordata;Actinopterygii;Salmoniformes;Salmonidae;Salmo;Salmo salar",
                "Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;Penicillium sp."),
    lineage_taxids = c("33208;7711;186623;8030;8016;8017;8018", "4751;4890;147545;5042;5043;5044;5045"))
  testthat::local_mocked_bindings(chimera_run_blastn = fake_blast, chimera_run_ribotyper = fake_rt,
    tk_ranked_lineage = fake_lin, .package = "RiboPilot")
  r <- RiboPilot:::chimera_unit(con, "S1", 1L, 1L, dbarg = "panel", taxdump_dir = NULL)
  expect_equal(r$label, "likely")
  w <- DBI::dbGetQuery(con, "SELECT start, lineage FROM chimera_windows WHERE ID='S1' ORDER BY start")
  expect_gt(nrow(w), 10)
  expect_match(w$lineage[1], "Actinopterygii")
  seg <- jsonlite::fromJSON(DBI::dbGetQuery(con, "SELECT segments FROM chimera")$segments)
  expect_equal(seg$taxon, c("Actinopterygii", "Eurotiomycetes"))
  expect_match(DBI::dbGetQuery(con, "SELECT reason FROM chimera")$reason, "Penicillium")
})
```

Check `jsonlite` is in DESCRIPTION Imports (grep). If not, add it. Adjust the INSERTs for any NOT NULL columns in `assemblies`, `samples`, `annotations`.

- [ ] **Step 2: Run tests to verify they fail**

- [ ] **Step 3: Implement**

Append to `R/chimera.R`:

```r
#' @noRd
chimera_label <- function(score) {
  if (is.na(score)) return(NA_character_)
  if (score < 0.3) "clean" else if (score <= 0.7) "suspect" else "likely"
}

#' Combine channels: max, corroboration bonus, taxonomy-only cap
#' @noRd
chimera_score <- function(channels) {
  sc <- vapply(channels, function(x) x$score %||% 0, numeric(1))
  score <- max(sc, 0)
  js <- lapply(channels, function(x) as.integer(x$junctions %||% integer(0)))
  corroborated <- integer(0)
  nm <- names(js)
  for (i in seq_along(js)) for (j in seq_along(js)) if (i < j) for (a in js[[i]]) {
    near <- js[[j]][abs(js[[j]] - a) <= CHIMERA_JUNCTION_TOL]
    if (length(near)) corroborated <- c(corroborated, min(a, near))
  }
  if (length(corroborated)) score <- min(1, score + 0.25)
  strong <- sc[setdiff(nm, c("genes", "windows"))]
  if (!length(strong) || all(strong == 0)) score <- min(score, CHIMERA_TAXONOMY_CAP)
  all_j <- sort(unique(unlist(js)))
  list(score = score, label = chimera_label(score),
    junctions = if (length(corroborated)) sort(unique(corroborated)) else all_j)
}

#' @noRd
chimera_reason <- function(segments, channels, junctions) {
  segtxt <- if (!is.null(segments) && nrow(segments)) paste(sprintf("%d-%d %s (%s, %.1f%%)",
    segments$start, segments$stop, segments$taxon, segments$species %||% "", segments$pident), collapse = "; ")
    else "no window segments"
  jtxt <- if (length(junctions)) {
    who <- vapply(junctions, function(j) paste(names(channels)[vapply(channels, function(x)
      any(abs((x$junctions %||% integer(0)) - j) <= CHIMERA_JUNCTION_TOL), logical(1))], collapse = ","), character(1))
    paste(sprintf("junction near %d (%s)", junctions, who), collapse = "; ")
  } else "no junction"
  extra <- vapply(channels[c("ribotyper", "reads")], function(x) x$reason %||% NA_character_, character(1))
  extra <- extra[!is.na(extra)]
  paste(c(segtxt, jtxt, extra), collapse = "; ")
}

#' @noRd
chimera_run_blastn <- function(query, dbarg, task, out, threads = 1L, blastn = "blastn") {
  system2(blastn, c("-db", shQuote(dbarg), "-query", shQuote(query), "-task", task,
    "-outfmt", shQuote("6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send slen qseqid"),
    "-max_target_seqs", "5", "-max_hsps", "1", "-num_threads", threads), stdout = out, stderr = FALSE)
}

#' @noRd
chimera_run_ribotyper <- function(fasta, outdir, ribotyper = "ribotyper") {
  system2(ribotyper, c("-f", shQuote(fasta), shQuote(outdir)), stdout = FALSE, stderr = FALSE)
}

#' Run every channel for one unit and write chimera + chimera_windows
#' @noRd
chimera_unit <- function(con, id, path = 1L, scaffold = 1L, dbarg = NULL, taxdump_dir = NULL,
                         label_db = NULL, threads = 1L, workdir = tempfile("chimera_")) {
  dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
  seq <- DBI::dbGetQuery(con, "SELECT sequence FROM assemblies WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$sequence[1]
  if (is.na(seq %||% NA) || !nzchar(seq)) return(invisible(NULL))
  feats <- DBI::dbGetQuery(con, "SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
    params = list(id, path, scaffold))
  hits <- tryCatch(DBI::dbGetQuery(con, "SELECT accession, bitscore, qstart, qend, lineage_taxids, lineage, species, pident FROM blast_hits WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  samp <- tryCatch(DBI::dbGetQuery(con, "SELECT lineage_taxids FROM samples WHERE ID = ?", params = list(id))$lineage_taxids[1],
    error = function(e) NA_character_)
  cov <- tryCatch(DBI::dbGetQuery(con, "SELECT depths, clips FROM coverage WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  itsx <- tryCatch(DBI::dbGetQuery(con, "SELECT chimeric FROM annotate WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold))$chimeric[1], error = function(e) NA)
  track <- function(x) if (length(x) && !is.na(x) && nzchar(x)) as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]) else NULL

  # channel 1
  gene_taxa <- if (nrow(hits) && nrow(feats)) chimera_gene_lca(hits, feats) else NULL
  ch_genes <- if (is.null(gene_taxa)) chimera_unavailable() else chimera_channel_genes(gene_taxa, feats)

  # channel 2
  tiles <- chimera_tile(nchar(seq))
  qf <- file.path(workdir, "windows.fa"); chimera_windows_fasta(seq, tiles, qf)
  identity_only <- is.null(dbarg) || !nzchar(dbarg)
  db_use <- if (identity_only) label_db else dbarg
  windows <- NULL; ch_windows <- chimera_unavailable(); segments <- NULL
  if (!is.null(db_use) && nzchar(db_use)) {
    out1 <- file.path(workdir, "mega.txt")
    chimera_run_blastn(qf, db_use, "megablast", out1, threads)
    h <- parse_blast_hits_file(out1)
    missing <- setdiff(sprintf("w%d_%d", tiles$start, tiles$stop), unique(h$qseqid))
    if (length(missing)) {
      keep <- tiles[sprintf("w%d_%d", tiles$start, tiles$stop) %in% missing, , drop = FALSE]
      qf2 <- file.path(workdir, "windows2.fa"); chimera_windows_fasta(seq, keep, qf2)
      out2 <- file.path(workdir, "dc.txt")
      chimera_run_blastn(qf2, db_use, "dc-megablast", out2, threads)
      h <- rbind(h, parse_blast_hits_file(out2))
    }
    windows <- chimera_window_hits(h, tiles, samp, taxdump_dir, db_source = db_use)
    segments <- chimera_segments(windows)
    ch_windows <- chimera_channel_windows(segments, identity_only = identity_only)
  }

  # channel 3
  ch_rt <- chimera_unavailable(); codes <- NA_character_
  rt_feats <- feats[feats$element %in% c("18S", "28S"), , drop = FALSE]
  if (nrow(rt_feats) && nzchar(Sys.which("ribotyper"))) {
    rf <- file.path(workdir, "rt.fa")
    st <- pmax(1L, rt_feats$start - 50L); en <- pmin(nchar(seq), rt_feats$stop + 50L)
    writeLines(as.vector(rbind(paste0(">u", rt_feats$element), substring(seq, st, en))), rf)
    rdir <- file.path(workdir, "rt"); dir.create(rdir, showWarnings = FALSE)
    chimera_run_ribotyper(rf, rdir)
    rt <- chimera_parse_ribotyper(list.files(rdir, "ribotyper\\.short\\.out$", full.names = TRUE)[1])
    ch_rt <- chimera_channel_ribotyper(rt, data.frame(target = paste0("u", rt_feats$element), start = st),
      itsx_chimeric = isTRUE(itsx == 1L))
    if (!is.null(rt)) codes <- paste(unlist(rt$codes), collapse = ";")
  } else if (isTRUE(itsx == 1L)) {
    ch_rt <- chimera_channel_ribotyper(NULL, NULL, itsx_chimeric = TRUE)
  }

  # channel 4
  depths <- if (nrow(cov)) track(cov$depths[1]) else NULL
  clips <- if (nrow(cov) && "clips" %in% names(cov)) track(cov$clips[1]) else NULL
  ch_reads <- chimera_channel_reads(depths, clips, feats)

  channels <- list(genes = ch_genes, windows = ch_windows, ribotyper = ch_rt, reads = ch_reads)
  sc <- chimera_score(channels)
  reason <- chimera_reason(segments, channels, sc$junctions)

  if (!is.null(windows) && nrow(windows) && !is.null(depths)) {
    windows$depth_median <- vapply(seq_len(nrow(windows)), function(i)
      stats::median(depths[windows$start[i]:min(windows$stop[i], length(depths))], na.rm = TRUE), numeric(1))
    windows$clip_frac <- if (!is.null(clips)) vapply(seq_len(nrow(windows)), function(i) {
      ix <- windows$start[i]:min(windows$stop[i], length(clips))
      max(clips[ix] / pmax(depths[ix], 1), na.rm = TRUE) }, numeric(1)) else NA_real_
  } else if (!is.null(windows)) { windows$depth_median <- NA_real_; windows$clip_frac <- NA_real_ }

  ts <- as.integer(Sys.time())
  row <- data.frame(ID = id, path = path, scaffold = scaffold, score = sc$score, label = sc$label,
    junctions = paste(sc$junctions, collapse = ","), reason = reason,
    segments = if (!is.null(segments) && nrow(segments)) jsonlite::toJSON(segments, dataframe = "rows") else NA_character_,
    gene_taxa = if (!is.null(gene_taxa)) jsonlite::toJSON(gene_taxa, dataframe = "rows") else NA_character_,
    channels = jsonlite::toJSON(channels, auto_unbox = TRUE),
    ribotyper_codes = codes, db_source = db_use %||% NA_character_, time_stamp = ts,
    stringsAsFactors = FALSE)
  DBI::dbWithTransaction(con, {
    DBI::dbExecute(con, "DELETE FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbExecute(con, "DELETE FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ?", params = list(id, path, scaffold))
    DBI::dbAppendTable(con, "chimera", row)
    if (!is.null(windows) && nrow(windows)) {
      windows$ID <- id; windows$path <- path; windows$scaffold <- scaffold
      DBI::dbAppendTable(con, "chimera_windows", windows[, c("ID", "path", "scaffold", "start", "stop", "accession",
        "species", "staxid", "lineage", "lineage_taxids", "rank_distance", "hit_class", "pident", "bitscore",
        "depth_median", "clip_frac", "db_source")])
    }
  })
  as.list(row)
}
```

The mocked test skips ribotyper because `Sys.which("ribotyper")` is empty on the test machine; the fake is still registered so the test passes either way. `toJSON` of `segments` must serialize as an array of objects (`dataframe = "rows"`) so `fromJSON` returns a frame with a `taxon` column.

- [ ] **Step 4: Run tests to verify they pass**

- [ ] **Step 5: Commit**

```bash
git add R/chimera.R tests/testthat/test-chimera.R DESCRIPTION
git commit -m "feat: chimera scoring and unit runner"
```

---

### Task 7: Driver, Nextflow process, wiring, Docker

**Files:**
- Create: `inst/nextflow/bin/rp_chimera.R`
- Modify: `inst/nextflow/modules/pipeline_r.nf` (after `screen`, ~line 225)
- Modify: `inst/nextflow/main.nf:160-184`
- Modify: `docker/Dockerfile:11-51`

**Interfaces:**
- Consumes: `chimera_unit()`, `rp_dbarg.R` (prints the DB arg or exits non-zero).
- Produces: process `chimera` emitting `tuple val(id), val(upath), val(uscaf)`.

- [ ] **Step 1: Driver**

`inst/nextflow/bin/rp_chimera.R`:

```r
#!/usr/bin/env Rscript
# Args: <db> <id> <path> <scaffold> [dbarg] [taxdump_dir] [label_db]
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
rp_args(a, 4, "rp_chimera.R <db> <id> <path> <scaffold> [dbarg] [taxdump_dir] [label_db]")
opt <- function(i) if (length(a) >= i && nzchar(a[i])) a[i] else NULL
r <- local({
  con <- open_project_db(a[1])
  on.exit(DBI::dbDisconnect(con))
  chimera_unit(con, a[2], as.integer(a[3]), as.integer(a[4]), dbarg = opt(5),
    taxdump_dir = opt(6), label_db = opt(7), threads = as.integer(Sys.getenv("RP_THREADS", "1")),
    workdir = file.path(getwd(), "chimera_work"))
})
cat("chimera:", a[2], a[3], a[4], "-", r$label %||% "skipped", "\n")
```

Make it executable: `chmod +x inst/nextflow/bin/rp_chimera.R`. Check how `rp_env.R` exposes package functions (it sources or loads `RIBOPILOT_SRC`); `chimera_unit` and `%||%` must be reachable the same way `screen_sample` is in `rp_screen.R`.

- [ ] **Step 2: Process**

Append to `inst/nextflow/modules/pipeline_r.nf`:

```groovy
// Chimera detection per unit: per-gene taxonomy from stored HSPs, windowed
// local BLAST, ribotyper on 18S/28S, depth steps and clip pileups.
process chimera {
    cache false
    tag "${id}.${upath}.${uscaf}"
    errorStrategy 'ignore'
    input:
    tuple val(id), val(upath), val(uscaf)
    output:
    tuple val(id), val(upath), val(uscaf)
    script:
    """
    DBARG=""
    if [ -z "\${RIBOPILOT_NO_LOCAL_BLAST:-}" ]; then
        DBARG=\$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}' 2>/dev/null || echo "")
    fi
    LABELS=\$(sqlite3 -cmd '.timeout 10000' '${params.db}' "
      SELECT COALESCE(ao.labels_db,'') FROM assemble a JOIN assemble_opts ao ON a.assemble_opts = ao.assemble_opts
      WHERE a.ID='${id}' LIMIT 1;" 2>/dev/null || echo "")
    RP_THREADS=${task.cpus} RIBOPILOT_SRC='${params.ribopilot_src}' rp_chimera.R \
      '${params.db}' '${id}' ${upath} ${uscaf} "\$DBARG" '${params.taxdump_dir}' "\$LABELS"
    """
}
```

Verify the `assemble` to `assemble_opts` join columns against `R/init_db.R:208-296` and `pipeline_r.nf:20-27` (which already select `ao.labels_db`); copy that join exactly. The label DB is a FASTA; `chimera_unit` passes it to blastn via `-db`, which needs a BLAST db. Add to `chimera_unit`, before using `label_db`: if `label_db` ends in `.fasta`/`.fa` and no `<label_db>.nin` exists, run `makeblastdb -in <label_db> -dbtype nucl -out <workdir>/labels` and use that prefix. Put this in a helper `chimera_label_prefix(label_db, workdir)` with a test that a temp FASTA yields a prefix whose `.nin` exists (skip the test with `skip_if(!nzchar(Sys.which("makeblastdb")))`).

- [ ] **Step 3: Wiring**

In `inst/nextflow/main.nf`, add `chimera` to the include from `pipeline_r.nf`, and after `screen(blast_genbank.out)`:

```groovy
    // Chimera detection needs annotation (reconcile), and reads depth/clips and
    // blast_hits from the db, so it waits on coverage and blast_genbank the same
    // way reference alignment waits on blast_ref_fetch: a scalar barrier, not a
    // per-unit join.
    chim_done = coverage.out.collect().map { 'ready' }.ifEmpty('ready')
        .combine(blast_genbank.out.collect().map { 'ready' }.ifEmpty('ready'))
        .map { 'ready' }
    chimera(reconcile.out.combine(chim_done).map { id, upath, uscaf, _r -> tuple(id, upath, uscaf) })
```

Then make triage wait on chimera. Replace `curate.out.map { ... }` in the triage input with a channel that requires both curate and chimera per unit:

```groovy
    unit_done = curate.out.map { id, upath, uscaf -> tuple("${id}.${upath}.${uscaf}", id) }
        .join(chimera.out.map { id, upath, uscaf -> tuple("${id}.${upath}.${uscaf}", true) }, remainder: true)
        .map { _k, id, _c -> id }
    triage(unit_done.map { id -> tuple(id, true) }
        .groupTuple()
        .map { id, _units -> tuple(id, true) }
        .join(screen.out.map { id -> tuple(id, true) }, remainder: true)
        .filter { entry -> entry[1] != null }
        .map { id, _a, _b -> id })
```

`remainder: true` keeps a unit whose chimera step was ignored. Check what `coverage.out` and `blast_genbank.out` emit (`coverage.nf`, `blast.nf` output blocks); `collect()` works on any shape.

- [ ] **Step 4: Docker**

In `docker/Dockerfile` first `mamba install` line add `'ribovore=1.0.5'` (check the current bioconda version with `mamba search -c bioconda ribovore` and pin what is available). If Task 2 found `awk` is mawk, add `'gawk'` too. Rebuild is not part of this task; note in the commit that the image needs a rebuild.

- [ ] **Step 5: Syntax check**

Run: `nextflow run inst/nextflow/main.nf -preview --db /dev/null 2>&1 | head` or, if `-preview` is unavailable in the installed version, `nextflow config inst/nextflow/main.nf` to parse the config, and `groovy` parse is covered by `nextflow run ... -stub-run` if stubs exist. At minimum confirm `nextflow run inst/nextflow/main.nf --help` does not print a compile error.

- [ ] **Step 6: Commit**

```bash
git add inst/nextflow/bin/rp_chimera.R inst/nextflow/modules/pipeline_r.nf inst/nextflow/main.nf docker/Dockerfile R/chimera.R tests/testthat/test-chimera.R
git commit -m "feat: chimera process, wiring, ribovore in image"
```

---

### Task 8: Triage rule, export risk, retire the curation flag

**Files:**
- Modify: `R/triage_rules.R:16-102` (rule list), `:144-164` (`triage_rule_ctx`)
- Modify: `R/export.R:190-229`
- Modify: `R/curate_rrna.R:145-148, 335-342`
- Test: `tests/testthat/test-chimera.R`, existing `test-triage*.R`, `test-export*.R`, `test-curate*.R`

**Interfaces:**
- Consumes: `chimera.label`, `chimera.reason`.
- Produces: triage rule key `chimeric`; `risk_flags(chimeric, ann, seq, chimera_label = NA, chimera_reason = NA)`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("triage rule chimeric fires on likely with the reason text", {
  d <- data.frame(ID = "S1", chimera_label = c("likely", "suspect"), chimera_reason = c("1-2500 A; 2600-4000 B", "x"),
    annotated_n = 1L, is_nuclear_rrna = 1L, n_features = 3L, stringsAsFactors = FALSE)
  ctx <- RiboPilot:::triage_rule_ctx(d)
  rule <- Filter(function(r) r$key == "chimeric", RiboPilot:::TRIAGE_RULE_DEFS)[[1]]
  expect_equal(rule$test(d, 1, ctx), "1-2500 A; 2600-4000 B")
  expect_true(is.na(rule$test(d, 2, ctx)))
})

test_that("risk_flags reports suspect and likely chimera labels", {
  r <- RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "suspect", chimera_reason = "why")
  expect_true(r$chimeric)
  expect_match(r$note, "chimeric")
  expect_false(RiboPilot:::risk_flags(0L, NULL, NA_character_, chimera_label = "clean")$chimeric)
})
```

Check what columns `triage_rule_ctx()`'s `d` normally has by reading the triage query in `R/triage.R` or wherever `d` is built (grep `triage_rule_ctx(`), and add `chimera_label` and `chimera_reason` to that query via `LEFT JOIN chimera USING (ID, path, scaffold)`. If `unit_concordance()` inside `triage_rule_ctx()` needs more columns than the test frame provides, add them with NA values in the test.

- [ ] **Step 2: Run tests to verify they fail**

- [ ] **Step 3: Implement**

In `R/triage_rules.R` append to `TRIAGE_RULE_DEFS`:

```r
  list(
    key = "chimeric",
    label = "Likely chimeric assembly",
    meaning = "Two or more independent signals place a junction between different lineages or array copies inside this scaffold.",
    actions = c("off", "flag"),
    default = "flag",
    test = function(d, i, ctx) {
      if (!identical(ctx$chimera[i], "likely")) return(NA_character_)
      ctx$chimera_reason[i] %|NA|% "chimera detector scored this scaffold as likely"
    }
  )
```

In `triage_rule_ctx()` add:

```r
    chimera = d$chimera_label %||% rep(NA_character_, nrow(d)),
    chimera_reason = d$chimera_reason %||% rep(NA_character_, nrow(d)),
```

In `R/export.R`, `submission_risks()` reads the label:

```r
  chm <- tryCatch(DBI::dbGetQuery(con, "SELECT label, reason FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?",
    params = list(id, path, scaffold)), error = function(e) data.frame())
  risk_flags(chim, ann, seq, chimera_label = if (nrow(chm)) chm$label[1] else NA_character_,
    chimera_reason = if (nrow(chm)) chm$reason[1] else NA_character_)
```

and `risk_flags()` gains the two args, with `chim <- isTRUE(chimeric[1] == 1L) || isTRUE(chimera_label %in% c("suspect", "likely"))`, and when `chimera_reason` is not NA append it to `note` as `chimeric(<reason>)` instead of the bare word. Update any existing export test that asserts the exact `note` string.

In `R/curate_rrna.R`, remove the `chimeric` flag emission at lines 145-148 and 335-342 (the ITSx token now feeds channel 3). Fix any curate test asserting that flag.

- [ ] **Step 4: Run tests**

Run: `Rscript -e 'devtools::test()'`
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add R/triage_rules.R R/triage.R R/export.R R/curate_rrna.R tests/testthat
git commit -m "feat: chimeric triage rule and export risk"
```

---

### Task 9: Scaffold table column

**Files:**
- Modify: `R/constants.R:130-153`, `R/utils_reactable.R:255-268` (add after `rt_screen_badge`), `R/app_server.R:1092, 1140-1160`
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces: `rt_chimera_badge()` JS renderer; scaffold table columns `chimera_label`, `chimera_reason`.

- [ ] **Step 1: Write the failing test**

```r
test_that("scaffold column help and badge exist", {
  expect_true(nzchar(RiboPilot:::SCAFFOLD_COL_HELP$chimera))
  js <- RiboPilot:::rt_chimera_badge()
  expect_s3_class(js, "JS_EVAL")
  expect_match(as.character(js), "likely")
})
```

- [ ] **Step 2: Run test to verify it fails**

- [ ] **Step 3: Implement**

`R/constants.R`, add to `SCAFFOLD_COL_HELP`:

```r
  chimera = paste("Chimera check: whether different parts of this scaffold match different",
    "lineages or show read-level breakpoints. clean, suspect, or likely. Hover for the",
    "lineage segments and junction positions. Warns only; never blocks export.")
```

`R/utils_reactable.R`:

```r
rt_chimera_badge <- function() {
  htmlwidgets::JS(
    "function(cellInfo) {
      var v = cellInfo.value;
      if (v === null || v === undefined || v === '' || v === 'NA') return '';
      var bg, fg;
      if (v === 'clean')        { bg = '#d4edda'; fg = '#2d6a4f'; }
      else if (v === 'likely')  { bg = '#f8d7da'; fg = '#842029'; }
      else                      { bg = '#fde8d0'; fg = '#7d4a1e'; }
      var tip = cellInfo.row.chimera_reason || '';
      return '<span title=\"' + String(tip).replace(/\"/g, '&quot;') + '\" style=\"background:' + bg +
        '; color:' + fg + '; padding:1px 6px; border-radius:8px; font-size:0.85em\">' + v + '</span>';
    }")
}
```

`R/app_server.R`: in the scaffold query near line 1092 add `LEFT JOIN chimera ch USING (ID, path, scaffold)` and select `ch.label AS chimera_label, ch.reason AS chimera_reason` (match the alias style used for `screen`). In the colDef list near line 1140 add:

```r
      chimera_label = reactable::colDef(name = "Chimera", width = 90, html = TRUE,
        header = rt_help_header("Chimera", SCAFFOLD_COL_HELP$chimera), cell = rt_chimera_badge()),
      chimera_reason = reactable::colDef(show = FALSE),
```

- [ ] **Step 4: Run tests, then launch the app on the test project and confirm the column renders**

Use the `run` skill for the app launch. Project location is in memory (`ribopilot-progress.md`).

- [ ] **Step 5: Commit**

```bash
git add R/constants.R R/utils_reactable.R R/app_server.R tests/testthat/test-chimera.R
git commit -m "feat: chimera column in scaffold table"
```

---

### Task 10: Chimera track in the annotate modal

**Files:**
- Modify: `R/app_modals.R:317` (placement), append `plot_chimera()` after `plot_coverage()` (~line 680)
- Modify: `R/app_server.R:1567-1569` (add renderers)
- Test: `tests/testthat/test-chimera.R`

**Interfaces:**
- Produces: `plot_chimera(id, path, scaffold)` returning a ggplot; outputs `chimera_plot`, `chimera_segments`, `chimera_windows`.

- [ ] **Step 1: Write the failing test**

```r
test_that("plot_chimera builds from stored windows and segments", {
  db <- tempfile(fileext = ".sqlite")
  con <- RiboPilot:::new_db(db)
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop) VALUES ('S1',1,1,'18S',1,1800)")
  DBI::dbExecute(con, "INSERT INTO chimera (ID, path, scaffold, score, label, junctions, segments) VALUES ('S1',1,1,0.9,'likely','2375',
    '[{\"start\":1,\"stop\":2500,\"taxon\":\"Actinopterygii\",\"rank\":\"class\",\"species\":\"Salmo salar\",\"pident\":99,\"n_windows\":10},
      {\"start\":2251,\"stop\":4500,\"taxon\":\"Eurotiomycetes\",\"rank\":\"class\",\"species\":\"Penicillium sp.\",\"pident\":97,\"n_windows\":9}]')")
  DBI::dbExecute(con, "INSERT INTO chimera_windows (ID, path, scaffold, start, stop, pident, lineage) VALUES ('S1',1,1,1,500,99,'a;b;Actinopterygii'), ('S1',1,1,2251,2750,97,'f;a;Eurotiomycetes')")
  DBI::dbDisconnect(con)
  withr::local_options(list(ribopilot.db = db))
  g <- RiboPilot:::plot_chimera("S1", 1L, 1L, db = db)
  expect_s3_class(g, "ggplot")
  expect_true(any(vapply(g$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))))
})
```

`plot_coverage()` uses `app_query()`, which reads the app's connection. Give `plot_chimera()` a `db = NULL` argument: when supplied, open a connection with `open_project_db()` and query with `DBI::dbGetQuery`; otherwise use `app_query()`. Wrap the two in a local `q <- function(sql, params)` closure.

- [ ] **Step 2: Run test to verify it fails**

- [ ] **Step 3: Implement**

`R/app_modals.R`:

```r
#' Lineage segments, window identity, and junctions under the depth track
#' @noRd
plot_chimera <- function(id, path = 1L, scaffold = 1L, db = NULL) {
  q <- if (is.null(db)) app_query else {
    con <- open_project_db(db); on.exit(DBI::dbDisconnect(con))
    function(sql, params = NULL) DBI::dbGetQuery(con, sql, params = params)
  }
  empty <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", 0, 0, label = msg) + ggplot2::theme_void()
  ch <- q("SELECT label, junctions, segments, channels FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  if (nrow(ch) == 0) return(empty("no chimera check yet (run WF1)"))
  seg <- if (!is.na(ch$segments[1])) jsonlite::fromJSON(ch$segments[1]) else data.frame()
  win <- q("SELECT start, stop, pident, lineage FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start", list(id, path, scaffold))
  feats <- q("SELECT element, start, stop FROM annotations WHERE ID = ? AND path = ? AND scaffold = ?", list(id, path, scaffold))
  junc <- if (nzchar(ch$junctions[1] %||% "")) as.integer(strsplit(ch$junctions[1], ",")[[1]]) else integer(0)
  chan <- if (!is.na(ch$channels[1])) jsonlite::fromJSON(ch$channels[1]) else list()
  jsrc <- vapply(junc, function(j) paste(names(chan)[vapply(chan, function(x)
    any(abs((unlist(x$junctions) %||% integer(0)) - j) <= CHIMERA_JUNCTION_TOL), logical(1))], collapse = ","), character(1))
  g <- ggplot2::ggplot()
  if (nrow(feats)) g <- g + ggplot2::geom_rect(data = feats, ggplot2::aes(xmin = start, xmax = stop, ymin = 0, ymax = 3),
    fill = "grey92", inherit.aes = FALSE)
  if (nrow(seg)) {
    seg$label <- sprintf("%s (%s, %.0f%%)", seg$taxon, seg$species, seg$pident)
    g <- g + ggplot2::geom_rect(data = seg, ggplot2::aes(xmin = start, xmax = stop, ymin = 1.6, ymax = 2.8, fill = label)) +
      ggplot2::geom_text(data = seg[seg$stop - seg$start > 600, ], ggplot2::aes(x = (start + stop) / 2, y = 2.2, label = taxon), size = 3)
  }
  if (nrow(win)) g <- g + ggplot2::geom_rect(data = win, ggplot2::aes(xmin = start, xmax = stop, ymin = 0.2, ymax = 1.2, alpha = pident), fill = "#3a5a8c") +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1), guide = "none")
  if (length(junc)) g <- g + ggplot2::geom_vline(data = data.frame(x = junc, src = jsrc), ggplot2::aes(xintercept = x, colour = src), linetype = "dashed")
  g + ggplot2::scale_y_continuous(breaks = c(0.7, 2.2), labels = c("window identity", "lineage segment"), limits = c(0, 3)) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL, colour = "junction") +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(legend.position = "bottom", panel.grid = ggplot2::element_blank())
}
```

Placement: after line 317's `rp_details("Read coverage", ...)` add:

```r
    rp_details("Chimera check", plotOutput("chimera_plot", height = "170px"),
      reactable::reactableOutput("chimera_segments"),
      rp_details("All windows", reactable::reactableOutput("chimera_windows"), open = FALSE)),
```

Check `rp_details()`'s signature at `R/app_modals.R:274` accepts multiple `...` children and an `open` argument; nest accordingly.

`R/app_server.R` next to `output$coverage_plot`:

```r
  output$chimera_plot <- renderPlot({
    edit_tick(); u <- sel_unit(); req(u); plot_chimera(u$ID, u$path, u$scaffold)
  }, height = 168)
  output$chimera_segments <- reactable::renderReactable({
    u <- sel_unit(); req(u)
    ch <- app_query("SELECT segments FROM chimera WHERE ID = ? AND path = ? AND scaffold = ?", list(u$ID, u$path, u$scaffold))
    req(nrow(ch) == 1, !is.na(ch$segments[1]))
    s <- jsonlite::fromJSON(ch$segments[1])
    reactable::reactable(s[, intersect(c("start", "stop", "taxon", "rank", "species", "pident", "n_windows", "db_source"), names(s))],
      compact = TRUE, pagination = FALSE, defaultColDef = reactable::colDef(minWidth = 70))
  })
  output$chimera_windows <- reactable::renderReactable({
    u <- sel_unit(); req(u)
    w <- app_query("SELECT start, stop, species, lineage, pident, bitscore, rank_distance, depth_median, clip_frac FROM chimera_windows WHERE ID = ? AND path = ? AND scaffold = ? ORDER BY start",
      list(u$ID, u$path, u$scaffold))
    req(nrow(w) > 0)
    reactable::reactable(w, compact = TRUE, pagination = TRUE, defaultPageSize = 15, searchable = TRUE)
  })
```

Junctions as trim suggestions: find where the trim UI reads its default cut positions (grep `trim_plot` and the trim numeric inputs in `R/app_server.R`), and pre-populate a helper text line `sprintf("suggested junctions: %s", junctions)` beneath the trim inputs. Read-only text, no behaviour change.

- [ ] **Step 4: Run tests, then open the app on the test project and check the modal**

- [ ] **Step 5: Commit**

```bash
git add R/app_modals.R R/app_server.R tests/testthat/test-chimera.R
git commit -m "feat: chimera track and tables in annotate modal"
```

---

### Task 11: Synthetic fixtures and end-to-end checks

**Files:**
- Create: `tests/testthat/fixtures/chimera/make_fixtures.R`, generated FASTA files
- Test: `tests/testthat/test-chimera.R`
- Modify: `dev/design/13_VALIDATION_FIXTURES.md:309` (expected outcome)

- [ ] **Step 1: Build fixtures**

`make_fixtures.R` (run once, commit its outputs) takes a real host rDNA scaffold and a fungal rDNA record from the existing test project or `dev/refdb/build/chunks` and writes:

- `host_fungal.fa`: host 1..(end of ITS1) joined to fungal 5.8S..28S.
- `host_variant.fa`: host scaffold with 28S positions 500..1500 (relative to 28S start) mutated at 6% of sites.
- `tandem_18S.fa`: host 18S duplicated in tandem before ITS1.
- `clean.fa`: unmodified host scaffold.

Record each splice position in `fixtures/chimera/expected.tsv` (`file`, `label`, `junction`).

- [ ] **Step 2: Write the test**

```r
test_that("synthetic fixtures get the expected label and a junction within tolerance", {
  skip_if(!nzchar(Sys.which("blastn")) || !nzchar(Sys.which("makeblastdb")), "blast not installed")
  exp <- read.delim(test_path("fixtures/chimera/expected.tsv"))
  panel <- RiboPilot:::rp_local_panel_prefix()
  skip_if(is.null(panel), "no local panel")
  for (i in seq_len(nrow(exp))) {
    db <- tempfile(fileext = ".sqlite"); con <- RiboPilot:::new_db(db)
    seq <- paste(readLines(test_path("fixtures/chimera", exp$file[i]))[-1], collapse = "")
    DBI::dbExecute(con, "INSERT INTO assemblies (ID, path, scaffold, sequence) VALUES ('F',1,1,?)", params = list(seq))
    r <- RiboPilot:::chimera_unit(con, "F", 1L, 1L, dbarg = panel)
    expect_equal(r$label, exp$label[i], info = exp$file[i])
    if (!is.na(exp$junction[i])) expect_true(any(abs(as.integer(strsplit(r$junctions, ",")[[1]]) - exp$junction[i]) <= 300), info = exp$file[i])
    DBI::dbDisconnect(con)
  }
})
```

`rp_local_panel_prefix()` is whatever `rp_dbarg.R` calls to discover the panel (read `inst/nextflow/bin/rp_dbarg.R` and use the same function; if it is only in the script, lift it into `R/blast_local.R` as an exported-internal helper). Without annotations only channel 2 runs, so `host_fungal.fa` is expected `suspect`, not `likely`. Encode that in `expected.tsv`.

- [ ] **Step 3: Nextflow smoke**

Open the test project in the app and run the assembly workflow on one sample (never script `run_assembly`; see memory `run-workflows-from-app.md`). Confirm `chimera` rows exist for every reconciled unit, and that a sample whose coverage failed still gets a row with `reads` marked unavailable in `channels`.

- [ ] **Step 4: Update doc 13 line 309 with the new outcome and commit**

```bash
git add tests/testthat/fixtures/chimera tests/testthat/test-chimera.R dev/design/13_VALIDATION_FIXTURES.md R/blast_local.R
git commit -m "test: chimera synthetic fixtures and smoke"
```

---

## Self-review notes

- Spec section 3.2 dc-megablast fallback: Task 6. Identity-only fallback and label DB makeblastdb: Tasks 4, 7.
- Spec section 3.3 ITSx joins channel 3: Task 5 and 6. Curate flag retired: Task 8.
- Spec section 4 scoring rules: Task 6 with tests for cap, bonus, and thresholds.
- Spec section 5 schema: Task 1. `lineage_taxids` added to `chimera_windows` beyond the spec so segments can be recomputed; harmless.
- Spec section 6 surfaces: Tasks 8, 9, 10. Trim suggestion is text only.
- Spec section 7 pipeline: Tasks 2, 7.
- Spec section 8 tests: Tasks 1 to 6, 11.
- Name consistency: `chimera_unit`, `chimera_score`, `chimera_label`, `chimera_segments`, `chimera_channel_*`, `chimera_run_blastn`, `chimera_run_ribotyper`, `rt_chimera_badge`, `plot_chimera` are used identically across tasks.
- Open risk: `rank_distance()` direction (Task 3 pins it with a test before anything depends on it).
