# Stage 4c: Ship the Panel in the Docker Image

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A Docker user gets fast local BLAST with no configuration at all, because the panel and its decoys are inside the image and the pipeline finds them by itself.

**Architecture:** The panel is baked into the image at `/opt/ribopilot/blastdb/<version>/` with `ENV RIBOPILOT_DB_DIR=/opt/ribopilot/blastdb`. Discovery runs **container-side** in `rp_dbarg.R`, never host-side, so the in-image path never enters `reference_config.blastdb` and therefore never reaches the bind-mount builder. An explicitly configured host panel keeps working exactly as today.

**Tech Stack:** R (>= 4.4.0), Nextflow DSL2, Docker, Singularity, BLAST+ 2.17.0, testthat 3e.

**Spec:** `dev/design/17_REMOTE_BLAST_PERFORMANCE.md`, Stage 4 section.

## Global Constraints

- **ASCII only in R code.** No non-ASCII anywhere in `R/`, `tests/`, `inst/`.
- **Minimal comments.** Comment why, not what. Never narrate a bugfix in a comment.
- **No new package dependencies.**
- **No schema migrations.** Pre-release; guard reads of possibly-missing columns.
- **Never push, never publish, never run a deploy script with a real registry.** Building an image locally is allowed; pushing is not.
- **Discovery must NEVER write into `reference_config.blastdb`.** That column stays NULL for a zero-config project. This is what keeps `reference_mount_dirs()` (`R/reference_nf.R:117-124`) from emitting a bind mount that would shadow the in-image panel. Leave that function completely alone.
- **Do not add a "drop paths under RIBOPILOT_DB_DIR" guard to the mount builder.** `install_ribopilot_db()` installs legitimately to exactly that location on a host, so such a rule breaks the HPC relocation case and a host panel newer than the image.
- **Never use a non-path sentinel (e.g. `"remote"`) in the blastdb column.** `dirname("remote")` is `"."`, which emits `-v .:.` and kills every task with exit 125.

## Measured Facts This Plan Relies On

- `nextflow run x.nf --blastdb ''` sets `params.blastdb` to the **String `"true"`**. Verified on this machine.
- Bind-mounting over an in-image directory **shadows** it on both Docker and Singularity, exit 0, silently.
- `ADD` of a local tarball auto-extracts and does not store the archive: 104 MB vs 204 MB for `COPY` + `RUN tar && rm`. The `rm` saves nothing because the COPY layer is already committed.
- `ARG GIT_VERSION_TAG` at `docker/Dockerfile:9` invalidates every layer below it on each publish. It is dead (2 grep hits) and duplicated by the OCI revision label.
- A partial panel (`.nin` present, `.nsq` truncated) can segfault blastn (exit 139), so `blast_db_exists()` is too weak a completeness test for discovery.
- Plain string sort mis-orders versions: `2026.08.9` sorts above `2026.08.10`.
- `as.character(numeric_version("2026.08.1"))` is `"2026.8.1"`, which is not a real directory name, so discovery must return the **raw basename**.

---

### Task 1: Stop passing empty parameters to Nextflow

**Files:**
- Modify: `R/nf_command.R:49-50`
- Modify: `inst/nextflow/modules/blast.nf:20` and the `if [ -n ... ]` gate
- Test: `tests/testthat/test-nf-command-quoting.R`, `tests/testthat/test-nextflow.R`

This is a standalone bug fix and lands first because everything else depends on `params.blastdb` being genuinely empty when nothing is configured.

**Interfaces:**
- Produces: `nf_command()` omits `--blastdb` and `--ref_dir` entirely when their values are empty.

- [ ] **Step 1: Write the failing tests**

```r
test_that("nf_command omits blastdb and ref_dir when they are empty", {
  args <- c("nextflow", "run", "--projectDir", "/p", "--ref_dir", "", "--blastdb", "")
  # the real assertion is on nf_command output; build a project with neither set
  skip_if_not(exists("mk_cfg_project"))
  expect_true(TRUE)
})
```

Replace that placeholder with a real test against `nf_command()`. Build a project the way `mk_cfg_project()` does in `tests/testthat/test-nf-command-quoting.R:4`, with a config that sets `rawDir` (otherwise `nf_command()` stops), leave `reference_config.blastdb` NULL, and assert:

```r
test_that("an unconfigured project passes no --blastdb flag at all", {
  d <- mk_cfg_project(tempfile(), c("params.rawDir = '/reads'"))
  args <- RiboPilot:::nf_command(d)
  expect_false("--blastdb" %in% args)
  expect_false(any(args == ""))
})

test_that("a configured panel still passes --blastdb with its value", {
  d <- mk_cfg_project(tempfile(), c("params.rawDir = '/reads'"))
  con <- RiboPilot:::open_project_db(file.path(d, ".sqlite"))
  RiboPilot:::set_ref_config(con, blastdb = "/panels/2026.08.1/ribopilot_rrna")
  DBI::dbDisconnect(con)
  args <- RiboPilot:::nf_command(d)
  i <- which(args == "--blastdb")
  expect_length(i, 1)
  expect_equal(args[i + 1], "/panels/2026.08.1/ribopilot_rrna")
})
```

- [ ] **Step 2: Run to verify failure**

Run: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-nf-command-quoting.R")'`
Expected: FAIL, `--blastdb` is present with an empty value.

- [ ] **Step 3: Implement**

In `R/nf_command.R`, replace the two unconditional entries:

```r
    "--taxdump_dir", taxdump,
    if (nzchar(ref_dir)) c("--ref_dir", ref_dir),
    if (nzchar(blastdb)) c("--blastdb", blastdb),
```

In `inst/nextflow/modules/blast.nf`, line 20, the closure can never see a container-side discovery, and the old conditional never worked because `params.blastdb` was the string `"true"`:

```groovy
    cpus 4
```

and drop the outer `if [ -n '${params.blastdb}' ]` gate so `rp_dbarg.R` is always consulted:

```bash
    DBARG=""
    SRC="remote"
    if DBARG=$(RIBOPILOT_SRC='${params.ribopilot_src}' rp_dbarg.R '${params.blastdb}'); then
        SRC="local"
    else
        DBARG=""
    fi
```

Note the "configured local panel is unusable" message is removed here: with discovery, no panel found is the normal unconfigured case, not an error. `rp_dbarg.R` keeps its own stderr warnings for a panel that was configured but is broken.

- [ ] **Step 4: Confirm `params.blastdb` defaults to empty**

`inst/nextflow/main.nf` must declare `params.blastdb = ''` (and `params.ref_dir` similarly) so the pipeline runs when the flag is absent. Verify it does; add it if missing.

- [ ] **Step 5: Run to verify pass, then the whole suite**

- [ ] **Step 6: Commit**

```bash
git add R/nf_command.R inst/nextflow/modules/blast.nf inst/nextflow/main.nf tests/testthat/test-nf-command-quoting.R
git commit -m "fix: omit empty nextflow parameters instead of passing them as true"
```

---

### Task 2: Panel discovery

**Files:**
- Modify: `R/blast_local.R`
- Test: `tests/testthat/test-blast-local.R` (create if absent)

**Interfaces:**
- Consumes: `ribopilot_db_dir()` from `R/db_dist.R`.
- Produces: `discover_panel(db_dir = ribopilot_db_dir())` returning a panel prefix or `NULL`.

- [ ] **Step 1: Write the failing tests**

```r
mk_panel <- function(root, ver, complete = TRUE) {
  d <- file.path(root, ver); dir.create(d, recursive = TRUE)
  ext <- if (complete) c("nin", "nsq", "nhr", "ndb") else c("nin")
  for (e in ext) writeLines("x", file.path(d, paste0("ribopilot_rrna.", e)))
  file.path(d, "ribopilot_rrna")
}

test_that("discover_panel finds a complete panel", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  p <- mk_panel(r, "2026.08.1")
  expect_equal(discover_panel(r), p)
})

test_that("discover_panel orders versions numerically, not as strings", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  mk_panel(r, "2026.08.9")
  newer <- mk_panel(r, "2026.08.10")
  expect_equal(discover_panel(r), newer)
})

test_that("discover_panel returns the raw directory name, not a normalised version", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  p <- mk_panel(r, "2026.08.1")
  got <- discover_panel(r)
  expect_true(file.exists(paste0(got, ".nin")))
  expect_false(grepl("2026[.]8[.]1", got))
})

test_that("discover_panel skips an incomplete panel and falls through", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  good <- mk_panel(r, "2026.08.1")
  mk_panel(r, "2026.08.2", complete = FALSE)
  expect_equal(discover_panel(r), good)
})

test_that("discover_panel returns NULL when there is nothing usable", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  expect_null(discover_panel(r))
  mk_panel(r, "2026.08.1", complete = FALSE)
  expect_null(discover_panel(r))
  expect_null(discover_panel(file.path(r, "absent")))
})

test_that("discover_panel ignores install leftovers", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  good <- mk_panel(r, "2026.08.1")
  mk_panel(r, "2026.08.2.old-123")
  mk_panel(r, ".stage-2026.08.3-99")
  expect_equal(discover_panel(r), good)
})

test_that("discover_panel refuses a path containing whitespace", {
  r <- file.path(tempfile(), "has space"); dir.create(r, recursive = TRUE)
  on.exit(unlink(dirname(r), recursive = TRUE))
  mk_panel(r, "2026.08.1")
  expect_null(discover_panel(r))
})
```

- [ ] **Step 2: Run to verify failure**

Expected: FAIL, "could not find function discover_panel".

- [ ] **Step 3: Implement**

Append to `R/blast_local.R`:

```r
#' Files a panel must have before blastn is pointed at it
#'
#' blast_db_exists() is deliberately looser (it gates a user-supplied prefix).
#' A directory holding only .nin passes that and then segfaults blastn, so
#' discovery, which chooses on the user's behalf, demands the whole set.
#' @noRd
PANEL_REQUIRED_EXT <- c("nin", "nsq", "nhr", "ndb")

#' TRUE when every required volume file is present and non-empty
#' @noRd
panel_complete <- function(prefix) {
  f <- paste0(prefix, ".", PANEL_REQUIRED_EXT)
  all(file.exists(f)) && all(file.info(f)$size > 0)
}

#' Newest usable panel installed under `db_dir`, or NULL
#'
#' Container-side only. The result must never be written to
#' reference_config.blastdb: an in-image path there would make the mount
#' builder bind a host directory over the panel and hide it.
#' @noRd
discover_panel <- function(db_dir = ribopilot_db_dir()) {
  if (length(db_dir) != 1 || is.na(db_dir) || !nzchar(db_dir)) return(NULL)
  db_dir <- path.expand(db_dir)
  if (!dir.exists(db_dir) || grepl("\\s", db_dir)) return(NULL)
  subs <- list.dirs(db_dir, full.names = FALSE, recursive = FALSE)
  subs <- subs[nzchar(subs) & !startsWith(subs, ".")]
  if (!length(subs)) return(NULL)
  v <- suppressWarnings(numeric_version(subs, strict = FALSE))
  keep <- !is.na(v)
  if (!any(keep)) return(NULL)
  subs <- subs[keep]
  subs <- subs[order(v[keep], decreasing = TRUE)]
  for (s in subs) {
    p <- file.path(db_dir, s, "ribopilot_rrna")
    if (panel_complete(p)) return(p)
  }
  NULL
}
```

- [ ] **Step 4: Run to verify pass**

- [ ] **Step 5: Commit**

```bash
git add R/blast_local.R tests/testthat/test-blast-local.R
git commit -m "feat: discover the newest complete installed panel"
```

---

### Task 3: Use discovery when no panel is configured

**Files:**
- Modify: `inst/nextflow/bin/rp_dbarg.R`
- Test: `tests/testthat/test-nextflow.R`

**Interfaces:**
- Consumes: `discover_panel()` from Task 2, `local_blast_db_arg()` from `R/blast_local.R`.
- Produces: `rp_dbarg.R` accepts an empty or absent first argument, meaning "discover".

- [ ] **Step 1: Write the failing test**

A static test on the script's text, in the style already used in `test-nextflow.R`:

```r
test_that("rp_dbarg.R falls back to discovery when given no panel", {
  txt <- readLines(system.file("nextflow/bin/rp_dbarg.R", package = "RiboPilot"), warn = FALSE)
  s <- paste(txt, collapse = "\n")
  expect_match(s, "discover_panel", fixed = TRUE)
  # an empty first argument must not be treated as a usable prefix
  expect_no_match(s, "rp_args(a, 1", fixed = TRUE)
})
```

Also add a behavioural test that runs the script for real, since a text test alone cannot prove it works:

```r
test_that("rp_dbarg.R discovers a panel from RIBOPILOT_DB_DIR with no argument", {
  skip_on_cran()
  skip_if(!nzchar(Sys.which("Rscript")))
  root <- tempfile(); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  d <- file.path(root, "2026.08.1"); dir.create(d)
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", file.path(d, paste0("ribopilot_rrna.", e)))
  script <- system.file("nextflow/bin/rp_dbarg.R", package = "RiboPilot")
  out <- suppressWarnings(system2("Rscript", c(script, "''"),
    env = c(paste0("RIBOPILOT_DB_DIR=", root), paste0("RIBOPILOT_SRC=", getwd())),
    stdout = TRUE, stderr = FALSE))
  expect_match(paste(out, collapse = " "), "2026.08.1", fixed = TRUE)
})
```

- [ ] **Step 2: Run to verify failure**

- [ ] **Step 3: Implement**

Rewrite `inst/nextflow/bin/rp_dbarg.R`:

```r
#!/usr/bin/env Rscript
# Print the composed -db argument for a local search, or exit 1 if unusable.
# Args: [panel_prefix] [decoy_dir]
# An empty or absent prefix means "use whatever panel is installed", which is
# how a containerised run finds the panel baked into the image.
a <- commandArgs(trailingOnly = TRUE)
rp_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(rp_file) != 1) stop("run with Rscript: cannot locate rp_env.R", call. = FALSE)
source(file.path(dirname(rp_file), "rp_env.R"))
local({
  panel <- if (length(a) >= 1 && !is.na(a[1]) && nzchar(trimws(a[1]))) trimws(a[1]) else ""
  if (!nzchar(panel)) {
    found <- RiboPilot:::discover_panel()
    if (is.null(found)) quit(status = 1L)
    panel <- found
  }
  # decoys ship beside the panel, so default to its own directory; the optional
  # second argument exists only to let a site keep them somewhere else
  dd <- if (length(a) >= 2 && nzchar(a[2])) a[2] else dirname(panel)
  arg <- RiboPilot:::local_blast_db_arg(panel, dd)
  if (is.null(arg)) quit(status = 1L)
  cat(arg, "\n", sep = "")
})
```

Note `rp_args()` is no longer called, because zero arguments is now valid.

- [ ] **Step 4: Run to verify pass**

- [ ] **Step 5: Commit**

```bash
git add inst/nextflow/bin/rp_dbarg.R tests/testthat/test-nextflow.R
git commit -m "feat: fall back to the installed panel when none is configured"
```

---

### Task 4: An opt-out, so remote stays reachable

**Files:**
- Modify: `R/blast_opts.R` (or wherever `blast_opts` is defined), `R/app_reference.R:163`
- Modify: `inst/nextflow/bin/rp_dbarg.R`
- Test: `tests/testthat/test-blast-local.R`

Discovery silently overloads the meaning of a blank panel field, which until now meant "search NCBI". A user who genuinely wants the remote search currently has no off-switch: `run_blast = 0` disables BLAST entirely, and `extra_opts = "-remote"` produces a duplicate-flag error that loses the sample.

- [ ] **Step 1: Decide the mechanism**

Add an environment-variable escape hatch honoured by `discover_panel()`, because it works identically for Docker, Singularity, conda, and bare metal, and needs no schema change:

```r
  if (nzchar(Sys.getenv("RIBOPILOT_NO_LOCAL_BLAST", ""))) return(NULL)
```

placed as the first line of `discover_panel()`.

- [ ] **Step 2: Write the failing test**

```r
test_that("RIBOPILOT_NO_LOCAL_BLAST forces discovery to find nothing", {
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  p <- mk_panel(r, "2026.08.1")
  expect_equal(discover_panel(r), p)
  withr::with_envvar(c(RIBOPILOT_NO_LOCAL_BLAST = "1"), {
    expect_null(discover_panel(r))
  })
})

test_that("the opt-out does not affect an explicitly configured panel", {
  # an explicit prefix bypasses discovery entirely
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE))
  p <- mk_panel(r, "2026.08.1")
  withr::with_envvar(c(RIBOPILOT_NO_LOCAL_BLAST = "1"), {
    expect_type(local_blast_db_arg(p, dirname(p)), "character")
  })
})
```

- [ ] **Step 3: Implement, and update the modal copy**

In `R/app_reference.R`, the helper text under the panel input currently says leaving it empty searches NCBI. That is no longer true inside Docker. Change it to say that when empty, RiboPilot uses a panel shipped with the image or installed locally if it finds one, and otherwise searches NCBI.

- [ ] **Step 4: Run to verify pass**

- [ ] **Step 5: Commit**

```bash
git add R/blast_local.R R/app_reference.R tests/testthat/test-blast-local.R
git commit -m "feat: add an opt-out that forces the remote search"
```

---

### Task 5: Bake the panel into the image

**Files:**
- Modify: `docker/Dockerfile`
- Modify: `docker/Dockerfile.dockerignore`
- Modify: `docker/deploy-dockerhub.sh`, `docker/deploy-local.sh`
- Modify: `.gitignore`
- Test: `tests/testthat/test-db-release.R`

**Do not push an image. Building locally is allowed; publishing is not.**

- [ ] **Step 1: Move the cache-busting ARG**

`ARG GIT_VERSION_TAG` at `docker/Dockerfile:9` invalidates every layer below it on each publish, which is why the database layer would otherwise rebuild and re-push every time. It has 2 grep hits repo-wide and the commit already ships as `org.opencontainers.image.revision`. Move it to immediately above the package `COPY` at line 87.

- [ ] **Step 2: Add the database layers**

After the last `mamba install` layer (line 83) and before the package `COPY`:

```dockerfile
# The panel and its decoys are versioned quarterly while the image versions on
# every package bump, so they sit below the package layer and stay cached.
# ADD extracts a local tarball in place and never stores the archive itself.
ARG RP_DB_VERSION=2026.08.1
ADD docker/ribopilot_rrna-${RP_DB_VERSION}.tar.gz /opt/ribopilot/blastdb/${RP_DB_VERSION}/
ADD docker/decoys-${RP_DB_VERSION}.tar.gz         /opt/ribopilot/blastdb/${RP_DB_VERSION}/
ENV RIBOPILOT_DB_DIR=/opt/ribopilot/blastdb
```

Use the exact `${RP_DB_VERSION}` filenames, never a glob: a glob matching zero files is a silent no-op that produces an image with no database and exit 0, and a glob matching two files merges two releases into one directory.

Rewrite the now-false comment at `docker/Dockerfile:92-93` that says no reference databases are baked in.

- [ ] **Step 3: Let the build see the tarballs**

`docker/Dockerfile.dockerignore` is `*` plus a whitelist. Add:

```
!docker/ribopilot_rrna-*.tar.gz
!docker/decoys-*.tar.gz
```

- [ ] **Step 4: Stage the tarballs in both deploy scripts**

In `docker/deploy-dockerhub.sh` before the build, and `docker/deploy-local.sh` likewise:

```bash
RP_DB_VERSION=$(awk -F'"' '/^RP_DB_VERSION/{print $2}' R/db_dist.R)
for A in "ribopilot_rrna-${RP_DB_VERSION}.tar.gz" "decoys-${RP_DB_VERSION}.tar.gz"; do
  if [ ! -f "docker/$A" ]; then
    echo "fetching $A ..."
    curl -fsSL -o "docker/$A" \
      "https://github.com/dmacguigan/RiboPilot/releases/download/db-${RP_DB_VERSION}/$A"
  fi
done
```

and pass `--build-arg "RP_DB_VERSION=${RP_DB_VERSION}"` to `docker build`.

Only fetch when absent: these are 82 MB and 77 MB and must not be re-downloaded per build. Do NOT extend the existing `rm -f docker/RiboPilot_*.tar.gz` cleanup to them.

- [ ] **Step 5: Keep the tarballs out of git**

Add to `.gitignore`:

```
docker/ribopilot_rrna-*.tar.gz
docker/decoys-*.tar.gz
```

Without this the deploy script's dirty-tree check trips and publishes a `-dirty` image with a bogus provenance warning.

- [ ] **Step 6: Extend the image smoke test**

`docker/deploy-dockerhub.sh` already smoke-tests tool presence. Add a check that the panel is really there and readable, with NO pipe (a pipe makes a successful run return 141):

```bash
docker run --rm "$IMAGE" bash -lc \
  'blastdbcmd -db "$RIBOPILOT_DB_DIR/'"${RP_DB_VERSION}"'/ribopilot_rrna" -info > /dev/null'
```

- [ ] **Step 7: Guard test**

Extend `tests/testthat/test-db-release.R` to assert the Dockerfile and the deploy script agree with the installer on asset names and on `RIBOPILOT_DB_DIR`. Guard with `skip_if_not(file.exists(p))` as that file already does for `dev/`.

- [ ] **Step 8: Verify without publishing**

Run `bash -n` on both deploy scripts. Confirm `docker/Dockerfile.dockerignore` permits the tarballs by checking the build context, e.g. `docker build --no-cache -f docker/Dockerfile . 2>&1 | head` fails on a MISSING FILE rather than silently, when the tarball is absent. Do not run a full 6 GB build unless it is quick to do so; if you do build, do not push.

- [ ] **Step 9: Commit**

```bash
git add docker/Dockerfile docker/Dockerfile.dockerignore docker/deploy-dockerhub.sh docker/deploy-local.sh .gitignore tests/testthat/test-db-release.R
git commit -m "feat: ship the rDNA panel inside the docker image"
```

---

### Task 6: Provenance for an in-image panel

**Files:**
- Modify: `R/db_dist.R`
- Modify: `inst/nextflow/bin/rp_dbarg.R` or `inst/nextflow/modules/blast.nf`
- Test: `tests/testthat/test-db-dist.R`

`record_blastdb_provenance()` reads `manifest.json` beside a **host** prefix. With an in-image panel there is no host prefix, so `blastdb_version` stays empty and the "Last run searched against" line never renders. Every zero-config Docker user loses traceability, which is the whole point of Stage 4b Task 5.

- [ ] **Step 1: Write the failing test**

```r
test_that("provenance can be recorded from a discovered panel", {
  con <- mk_prov_db(); on.exit(DBI::dbDisconnect(con))
  r <- tempfile(); dir.create(r); on.exit(unlink(r, recursive = TRUE), add = TRUE)
  d <- file.path(r, "2026.08.1"); dir.create(d)
  for (e in c("nin", "nsq", "nhr", "ndb")) writeLines("x", file.path(d, paste0("ribopilot_rrna.", e)))
  writeLines('{"db_name":"ribopilot_rrna","db_version":"2026.08.1"}', file.path(d, "manifest.json"))
  p <- discover_panel(r)
  expect_equal(record_blastdb_provenance(con, p), "2026.08.1")
  expect_equal(blastdb_provenance(con), "ribopilot_rrna 2026.08.1")
})
```

- [ ] **Step 2: Implement**

`record_blastdb_provenance()` already takes any prefix and reads the manifest beside it, so it needs no change; what is missing is a caller that knows the discovered prefix. Have the pipeline report the prefix it actually used: `rp_dbarg.R` prints the `-db` argument on stdout, so `blast.nf` already has it. Write the resolved panel back to the project database from `rp_searched.R`, which already has the connection and already records the search source, by passing it the first prefix of `DBARG`.

Keep it defensive: a prefix with no manifest records nothing, exactly as today.

- [ ] **Step 3: Run to verify pass, then the whole suite**

- [ ] **Step 4: Commit**

```bash
git add R/db_dist.R inst/nextflow/bin/rp_searched.R inst/nextflow/modules/blast.nf tests/testthat/test-db-dist.R
git commit -m "feat: record provenance for a discovered panel"
```

---

## Deferred, with reasons

- **Dropping `taxdb` from the image.** It is 186 MB and probing showed `staxids` resolves from the v5 `.nto/.not/.ntf` volumes with byte-identical output, so it is not needed for our searches. Kept anyway because it is what the release ships and because `-taxids` filtering would break without it. Revisit if image size becomes a problem on HPC.
- **Singularity bind-parse bug.** `R/reference_nf.R:79` uses `[^\\s,]+`, which captures a trailing quote so `drop_covered()` cannot dedup a config-declared `-B`. Harmless double bind today, one-character fix, independent of this feature.
- **Truncated-but-present panel files.** A truncated `.nsq` passes every check, blastn exits 0 and returns nothing. No cheap detector exists; `blastdbcmd -info` also exits 0. Image-baked panels are immutable so this mainly threatens hand-copied host panels.

## Self-Review

**Spec coverage.** The three approved decisions map to tasks: stage a local tarball (Task 5 Step 4), ship decoys (Task 5 Step 2), auto-discover (Tasks 2 and 3). Task 1 is the prerequisite bug fix that makes an unconfigured project actually pass an empty parameter. Task 4 restores the remote escape hatch that discovery would otherwise remove. Task 6 keeps provenance working for the users this feature is for.

**Placeholders.** Task 1 Step 1 opens with a placeholder test that the step text then explicitly replaces; the real tests follow it. No other step contains one.

**Type consistency.** `discover_panel()` returns a character prefix or `NULL`, and both callers (`rp_dbarg.R`, the Task 6 test) treat `NULL` as "nothing found". `panel_complete()` takes the same prefix shape `blast_db_exists()` does. `PANEL_REQUIRED_EXT` is used only inside `panel_complete()`.

**Ordering.** Task 1 must land first. Task 2 before Task 3 and 4. Task 5 is independent of 1-4. Task 6 needs 2 and 3.
