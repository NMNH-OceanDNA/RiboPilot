# Reference panel: streamlining plan

Branch: `ux/reference-panel-streamline`
Date: 2026-08-12

Review method: five simulated users of differing computational skill (bench
biologist, PhD student, staff bioinformatician, HPC specialist, returning
postdoc), two UX reviewers, and two adversarial verification passes that
checked every claim against the source. Claims that failed verification were
dropped before this plan was written.

Headline: **all five simulated users failed to complete the task unaided.**
Three of them reached a page that *looked* finished while their samples were
still pointed at the built-in ray-finned fish database.

---

## Diagnosis

The page is a per-row state machine driven by a bulk verb bar.

Each row sits at one point of a lifecycle and admits exactly one sensible
verb, but all six verbs live at the top with four different unwritten
selection rules (none / any / exactly one / one-but-silently-uses-the-first).
Before every click the user has to guess which rows are affected, how many are
allowed, and whether the verb is even legal for that row's state. Every
refusal arrives *after* the click.

Meanwhile the row's own state is six unlabelled icons, and the one that means
"waiting for you" draws a running person. Nothing on the page says what the
tab is for, or that doing nothing leaves the samples on a fish database.

### The three findings that matter most

1. **The running-person icon.** A freshly proposed group is
   `reference_switch = 1`, which `REFERENCE_STATE_ICONS` draws as
   `fa fa-person-running` (`R/app_server.R:3-5`). Four of five users read it as
   "the app is already working" and waited. It means the app is waiting for
   *them*. No tooltip, no legend, and Reference mode has no state filter row to
   decode it against.
2. **Skip never says what it costs.** The confirmation reads "their samples
   keep their current database" (`R/app_reference.R:350`). It never says that
   database is a ray-finned fish reference. The bench biologist and the PhD
   student both said they would happily skip and assemble sea stars and
   jellyfish against fish seeds. This is a silent scientific failure, not a
   usability nit.
3. **Suggest is the leftmost, friendliest control and the only irreversible
   one.** `write_reference_groups()` does `DELETE FROM reference` plus
   `UPDATE samples SET ref_id = NULL` (`R/reference_build.R:34-70`), carrying
   state forward only when ref_id **and** clade taxid both match. Any clade the
   user fixed by hand was re-keyed by `unique_split_ref_id()`, so it does not
   match and is silently discarded. Four of five users pressed it a second time
   when unsure, and the returning-user case (new samples added) forces everyone
   through it.

---

## Recommendation

Move the verb into the row, put the state in words, add one sentence of
purpose. The top bar keeps only the two genuinely global controls.

No new package, no new JS helper. Reuses `rt_link`, `rt_icon_link`,
`opts_help`, `conditionalPanel`, `shinyjs`, and a SQL `CASE` for the state
words. The codebase already ships in-cell action controls
(`rt_triage_bttn`, `R/utils_reactable.R:143`), so per-row verbs need no new
machinery.

```
BEFORE
[Reference v] (sync) | Suggest | Edit Clade | Use Existing | Skip | Build | Folder |  3 selected  [CSV] [Work Dirs]
 (i) | ID         | Clade        | Rank   | Samples | GenBank records | Source | Database | Note | Updated
  ^                  ^                                                                      ^
  no legend,         3 doors open one modal (button + Clade cell + ID cell)                 empty link on
  no tooltip                                                                                every healthy row

AFTER
[Reference v] (sync) | Propose groups | Settings |                 3 group(s) selected  [CSV]
Samples assemble against RiboPilot's built-in ray-finned FISH database unless their group has its own.
4 groups: 1 ready, 2 need a database, 1 failed.  6 samples are in no group.

 Status           | Group      | Clade        | Rank   | Samples | GenBank recs | Database | Note      |(log)|                 |
 Needs a database | asteroidea | Asteroidea   | class  | 12      | 4,231        |          |           | [f] | Set up database |
 Building         | cnidaria   | Cnidaria     | phylum | 54      | 91,004       |          |           | [f] | Change          |
 Ready            | echinoidea | Echinoidea   | class  | 9       | 3,110        | /refs/.. |           | [f] | Change          |
 Build failed     | octocoral  | Octocorallia | class  | 7       | 1,204        |          | rate limit| [f] | Try again       |

6 buttons -> 2.
```

Batch work survives: the clicked row folds into the current selection when it
is part of it, so selecting eight rows and clicking any one of them acts on
all eight. No action can be refused after the click.

**Honest runner-up:** leave the bar alone and only fix the table (state words,
column help, log icon, renamed buttons). That is changes 1, 2, 4, 5, 6 without
3, roughly a tenth of the diff and near-zero test churn, and it fixes both
silent scientific failures. It leaves the six-button bar and the arity trap
exactly as they are, which is the original complaint.

---

## The changes

Ordered smallest-diff-first. Each is independently shippable: stop after any
one and the page is better than before.

### 1. State in words, not six mute icons (small)

Kills finding 1. Also deletes the `Source` column, which is never explained
anywhere and flips meaning when Build re-wires from disk.

Files: `R/app_reference.R`, `R/app_server.R`,
`tests/testthat/test-reference-app.R`, `tests/testthat/test-reference-modal.R`

- `reference_table()`: drop `COALESCE(source, '') AS source,` and insert as the
  second selected column:
  ```sql
  CASE COALESCE(reference_switch,1)
    WHEN 2 THEN 'Ready'
    WHEN 3 THEN 'Build failed'
    WHEN 4 THEN 'Building'
    WHEN 5 THEN 'Left on current database'
    ELSE 'Needs a database' END AS status,
  ```
  Doing it in SQL keeps the column order right with no post-hoc reordering.
  Update the empty-frame fallback (lines 19-22) to match.
- `R/app_server.R`: delete `REFERENCE_STATE_ICONS` (lines 1-5). In the
  Reference block (lines 257-270) delete `cols[["source"]]`, replace the state
  override with `cols[[state_col]] <- reactable::colDef(show = FALSE)` (the
  value must stay in the data, the shared `rowClass` JS at 271-274 reads it),
  and add
  `cols[["status"]] <- reactable::colDef(name = "Status", sticky = "left", minWidth = 150)`.
  Leave it filterable: with no state pickers in this mode, the column's own
  filter box becomes the state filter for free.
- Tests: delete the `skipped rows get their own state icon` block
  (test-reference-app.R:17-20); add `expect_equal(d$status[1], "Ready")` to the
  existing `reference_table exposes ref_id as ID` test.
- Check the two sticky-left columns (Status 150 + ID 130) still leave room for
  Clade at 1280px.

Result: the first column of every row reads "Needs a database", "Building",
"Ready", "Build failed", or "Left on current database". The dead hourglass
icon goes with it (nothing ever writes `reference_switch = 0`; verified).

### 2. One line above the table: purpose, and where the project stands (small)

Kills half of finding 2. The only instruction on the page currently lives in
the empty-table string (`R/app_server.R:281`) and vanishes the moment there is
a row to act on.

Files: `R/app_ui.R`, `R/app_reference.R`

- `R/app_ui.R`: after the Annotate filter conditionalPanel (lines 94-96), in
  the same per-mode slot:
  ```r
  shiny::conditionalPanel("input.mode == 'Reference'",
    shiny::div(style = "margin:0.6em 0 0.7em; font-size:0.86em; color:#555;",
      shiny::uiOutput("ref_banner"))),
  ```
  The outer layout is a flex column with the table at `flex:1 1 auto;
  min-height:0`, so this shrinks the table rather than adding a page scrollbar.
- `R/app_reference.R`, inside `reference_handlers()`: add an
  `output$ref_banner <- shiny::renderUI({...})` that reads `deps$refresh()`,
  `reference_table(con)`, and
  `SELECT COUNT(*) FROM samples WHERE ref_id IS NULL`.
  Build the counts with `table(d$status)` so the wording cannot drift from
  change 1.
  - Empty: "Samples assemble against RiboPilot's built-in ray-finned fish
    database unless their group has its own. Press Propose groups to group your
    samples by clade, then give each group a database."
  - With rows: "N groups: a ready, b need a database, c left on the fish
    database, d failed. M sample(s) are in no group (they keep the fish
    database)." Omit the orphan sentence when M is 0.

Result: a permanent one-or-two line summary that states the purpose, names the
fish default out loud, tallies the states, and surfaces samples in no group
(currently invisible on the page and in every later warning).

### 3. Put the verb in the row; cut the bar from six buttons to two (medium)

The main event. Kills the select-then-guess-the-arity dance, `Edit Clade`
silently acting on only the first of several selected rows
(`R/app_reference.R:219`), three doors to one modal, and the two most
clickable-looking cells refusing to select the row. Gives Skip somewhere to
state what it actually costs.

Files: `R/app_ui.R`, `R/app_reference.R`, `R/app_server.R`,
`tests/testthat/test-reference-app.R`, `tests/testthat/test-reference-modal.R`,
`vignettes/Reference-Databases.Rmd`

- **Data.** `reference_table()`: append as the last selected column
  ```sql
  CASE COALESCE(reference_switch,1)
    WHEN 3 THEN 'Try again'
    WHEN 1 THEN 'Set up database'
    WHEN 0 THEN 'Set up database'
    ELSE 'Change' END AS action
  ```
  and add `action = character(0)` to the empty-frame fallback. A colDef with no
  matching data column is dropped at `R/app_server.R:294`, so the column must
  exist in the frame.
- **Column.** `R/app_server.R`, Reference block:
  ```r
  cols[["action"]] <- reactable::colDef(name = "", minWidth = 140,
    sticky = "right", filterable = FALSE, html = TRUE, cell = rt_link("ref_row_action"))
  cols[["ID"]] <- reactable::colDef(name = "Group", minWidth = 130)  # no rt_link: removes the third door
  ```
  Delete the Reference branch of the `open_detail` handler (lines 334-335). Add
  `"action"` to the CSV drop list (line 93). In `output$n_selected` say
  "group(s) selected" in Reference mode.
- **Handler.** Lift the bodies of the `ref_skip_btn` (343-351),
  `ref_existing_btn` (354-362) and `ref_build_btn` (404-437) observers into
  local functions `do_skip(ids)`, `open_existing(ids)`, `do_build(ids)` that
  take ids instead of calling `deps$selected_ids()`. Delete the four bar
  observers (including `ref_edit_btn` at 236). Add:
  ```r
  action_target <- shiny::reactiveVal(character(0))
  shiny::observeEvent(input$ref_row_action, {
    i <- suppressWarnings(as.integer(input$ref_row_action)); d <- deps$table_data()
    if (is.na(i) || i < 1 || i > nrow(d)) return()
    sel <- deps$selected_ids()
    ids <- if (d$ID[i] %in% sel) sel else d$ID[i]   # fold the click into the selection; never refuse
    action_target(ids); shiny::showModal(reference_action_modal(ids))
  })
  ```
  plus three one-line observers on `ref_act_build` / `ref_act_existing` /
  `ref_act_skip` that `removeModal()` then call the matching function.
- **Modal.** New pure function `reference_action_modal(ids)` beside the other
  two modal builders. Three stacked actionButtons, each followed by an
  `opts_help()` line:
  - "Download sequences from GenBank" / "Builds a seed database for this clade
    from GenBank. This is the normal choice."
  - "Use a database I already have" / "Point these samples at seed and label
    FASTA files already on disk, for example one a colleague built."
  - "Leave these samples on the built-in fish database" / "Nothing is built.
    These samples assemble against RiboPilot's ray-finned fish reference, which
    is usually wrong outside fish."
  Footer: `modalButton("Cancel")`.
- **Batch.** `open_existing(ids)` keeps its manifest modal but stores the
  vector; `ref_existing_apply` loops over `existing_target()` and sums the
  counts. Do **not** copy MitoPilot's `open_opts_modal()`: it `req(F)`s when the
  clicked row is outside the selection, which is the silent refusal this change
  exists to kill.
- **Bar.** `R/app_ui.R:58-62` becomes
  ```r
  shiny::conditionalPanel("input.mode == 'Reference'",
    ctrls(bttn("ref_suggest_btn", "Propose groups", "wand-magic-sparkles"),
          bttn("ref_dir_btn", "Settings", "gear"))),
  ```
  Also wrap the Work Dirs actionBttn (82-83) in
  `shiny::conditionalPanel("input.mode != 'Reference'", ...)`: in Reference mode
  `selected_ids()` returns group keys (`R/app_server.R:415`), which it feeds to
  a picker of sample IDs. It is broken there today.
- **Tests.** test-reference-app.R:10-16 asserts only `ref_suggest_btn` and
  `ref_dir_btn`. Expose `do_build`, `do_skip`, `open_existing` in the list
  `reference_handlers()` returns (the same seam `apply_clade` already uses) and
  change the six `session$setInputs(ref_build_btn = 1)` call sites in
  test-reference-modal.R to call `do_build("asteroidea")`.
- **Docs.** Rewrite the button bullets in `vignettes/Reference-Databases.Rmd`
  lines 33-84.

### 4. Make the build log reachable from every row (trivial)

The Note column renders the failure text as its own link label
(`R/app_server.R:266-267`), so on every healthy, running, or cluster-submitted
row it is an empty anchor: the log is unreachable exactly when people want it
(during a run, after an overnight cluster job, on success).

Files: `R/app_reference.R`, `R/app_server.R`

- `reference_table()`: add `'' AS log,` after the fail_reason column (and
  `log = character(0)` to the fallback). `rt_icon_link()` ignores the cell
  value, so the column only has to exist.
- Reference block: change the fail_reason cell from `rt_link("ref_open_log")`
  to `rt_longtext()`, and add
  ```r
  cols[["log"]] <- reactable::colDef(name = "", width = 55, align = "center",
    filterable = FALSE, html = TRUE,
    cell = rt_icon_link("ref_open_log", "fa fa-file-lines", "Build log"))
  ```
  mirroring the Assemble Details/Output columns at 215-218. Add `"log"` to the
  CSV drop list.
- Log modal fallback (line 619): name where a cluster job writes its log
  instead of the dead-end "No log for this group yet."

### 5. Explain the build window in the app's own idiom (trivial)

The tab's most consequential dialog is the only options modal in the app with
no per-field help: an API key, two checkboxes, and a raw Nextflow command, with
nothing saying the defaults are fine. Two users stalled here and one fell back
to Skip, which is the scientifically wrong answer. On a server the Run button
silently disappears and is replaced by a shell script with no explanation.

Files: `R/app_modals.R`, `R/app_reference.R`. Wording only, no behaviour change.

- One muted line above the inputs: "The settings below are fine for most
  projects. Press Run from App." (headless: "...Submit the job script below.")
- `opts_help()` under each control: API key "Optional. A free NCBI key raises
  the download rate limit."; RefSeq "Restricts the download to curated RefSeq
  records, which is usually far fewer and often too few outside model groups.";
  resume "Continue an interrupted run instead of starting over."
- Wrap the Nextflow command preview in
  `tags$details(tags$summary("Advanced: show the Nextflow command"), ...)`.
  `details` is already styled in `inst/app/www/custom.css:69-98`.
- Headless section (356-368): lead with plain language, e.g. "RiboPilot is
  running on a server, so it cannot start the download itself. Submit the job
  script below to the cluster, or save it and hand it to whoever runs jobs for
  you."
- Add `shiny::outputOptions(output, "ref_nf_cmd_text", suspendWhenHidden = FALSE)`
  next to the renderText, or the preview stays blank when the fold is first
  opened.

### 6. Ask before re-proposing over an existing table (trivial)

Kills finding 3.

Files: `R/app_reference.R`, `tests/testthat/test-reference-modal.R`

- Wrap the current `ref_suggest_btn` observer body (156-206) in a local
  `run_suggest()`. The observer becomes a count check: if
  `SELECT COUNT(*) FROM reference` is > 0, show a confirm modal and return;
  otherwise call `run_suggest()` directly. Add an `ref_suggest_ok` observer that
  removes the modal and runs it.
- Body text: "This re-groups every sample from scratch and replaces the N
  group(s) below, including any clade you changed by hand. Databases already
  built stay in the reference folder and are re-wired when you set the group up
  again." Plus: "The first run can take a few minutes and the app will not
  respond while it asks GenBank how many records each candidate clade has."
- Use a `modalDialog`, not a sweet alert: alerts do not queue here and the
  progress alert fires immediately inside `run_suggest()`.
- Tests: the two Suggest tests (test-reference-modal.R:222, :253) run against a
  fixture that pre-inserts a reference row, so each needs
  `session$flushReact(); session$setInputs(ref_suggest_ok = 1)` after the
  `ref_suggest_btn` line.

---

## Rejected, and why

- **A "Set up references" wizard modal.** Wrong container. A wizard models a
  one-shot linear task; this workflow is asynchronous and resumable (a cluster
  build finishes hours later and is wired by the poller; users return the next
  morning or three weeks later with new samples). Every step must stay
  reachable outside it anyway, so you build the row-level surface *and* the
  wizard. Its one good moment, the first run against an empty table, is covered
  by change 2's banner.
- **A JS action cell with several state-dependent buttons** (a `rt_ref_action`
  helper modelled on `rt_triage_bttn`). About 35 lines of new JS to say what one
  `rt_link` cell plus a three-choice modal says with none. The modal is also
  strictly clearer: it can afford full sentences with an `opts_help()` line
  under each choice, which is the only place Skip can honestly state its cost.
  Buttons in a cell force one-word labels, which is how the current bar got into
  trouble.
- **Contextual toolbar (show/hide bar buttons by selection).** Controls that
  appear and vanish stop the user ever building a stable map of the page. It
  degrades badly on a mixed-state multi-row selection, where the toolbar must
  either show the intersection (usually empty, reads as a broken app) or lie.
  It also visibly reflows the CSV / Work Dirs cluster pinned with
  `margin-left:auto`.
- **One primary button plus a `dropdownButton` overflow menu.** Cheapest on the
  list (about 15 lines in app_ui.R, zero server edits, every test still passes)
  but it is a cosmetic reorganisation of a semantic problem. Every hidden
  precondition stays hidden, the arity refusals still fire from behind a menu,
  and it buries Skip, the one action a user most needs to understand. Keep in
  the back pocket if change 3 is rejected outright.
- **Grey out inapplicable bar buttons (`shinyjs::toggleState`).** This is the
  honest fallback if change 3 is rejected, not a competitor: about 25 lines to
  convert every after-the-fact refusal into a before-the-fact constraint. But it
  keeps the verb at the wrong scope, and shinyWidgets' `material-flat`
  `actionBttn` paints its own background, so a disabled button looks identical
  to an enabled one until you add CSS.
- **Hover text on the state icons (`titles=` on `rt_dynamicIcon`).** Best value
  per line in the whole review, but beaten here by change 1: a word needs no
  hover, no legend, and no filter row. Not dead, just out of scope: the other
  four modes still show mute icons and would benefit. See open questions.
- **Linking the vignette from the tab.** There is no route from a running app to
  a vignette. golem serves only `inst/app/www`, and the repo is not public yet.
  Real work, not a one-liner. Change 2's banner carries the sentence people
  needed.

---

## Correctness bugs found, deliberately NOT in this plan

This is a clarity plan. These are real defects the review turned up and should
go on their own branch, in roughly this order:

1. Changing a group's clade after it is built leaves its samples wired to the
   old database while the row shows as unbuilt. Silent provenance error.
2. A group that gains new samples through a re-Suggest keeps its green tick
   while the new members are never wired.
3. One error inside `poll_reference_builds()` fails **every** in-flight build on
   the project: the recovery UPDATE has no `ref_id` predicate
   (`R/app_reference.R:137-140`).
4. The build manifest path is fixed per project, so opening the Build modal can
   rewrite the work list of a job that is still queued.
5. The double-submit guard is a session `reactiveVal`, so a browser reload lets
   the same group be submitted twice.
6. A cancelled or dead cluster job sits at "Building" forever: rows with no pid
   are skipped by the poller.
7. No way to force a rebuild from the app. `start_reference_build()` already
   takes `overwrite=` and no caller passes it, and the `reference_missing_dbs()`
   gate short-circuits in three places before a build starts.

---

## Open questions for the repo owner

1. Should "Propose groups" become **additive** (group only samples that
   currently have no ref_id, leave finished rows alone) instead of a full
   rewrite? That is what the returning-user case needs and it makes change 6's
   confirmation unnecessary, but it is a real change to
   `write_reference_groups()` and its tests, not a UI edit.
2. Colour the Status words (red for failed, green for ready) or leave plain?
   Plain is what is specified above; colour is a one-line JS cell renderer.
3. The other four modes still show the same mute state icons. Add `titles=` to
   `rt_dynamicIcon()` and give all of them hover text, or leave for a separate
   pass?
4. Is the clicked-row-plus-selection fold enough batching for a 400-sample
   project, or keep a single "Set up selected" button in the bar for bulk work?
   Keeping it puts the bar back to three buttons and reintroduces one arity rule.
5. Should the row action modal offer "Put this group back to needing a database"
   for a skipped row? Five-line sibling of `skip_reference_group()`. Note a
   skipped row is not actually locked out today: selecting it and building or
   wiring it already works.
