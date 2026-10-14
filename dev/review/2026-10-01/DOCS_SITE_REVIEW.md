# RiboPilot documentation site review

Date: 2026-10-01. Branch reviewed: `devel` (0.0.0.9019). Nothing in the repo was changed by the review.

## How the review was run

Six independent reviewers, each with a different background:

| Reviewer | Focus |
|---|---|
| Novice (wet-lab grad student) | Read the site start to finish; where a beginner gets stuck |
| Intermediate (bioinformatician, new to rDNA) | Site organization, conflicts between pages, repetition, links, doc claims vs code |
| Expert (rDNA and GenBank specialist) | Scientific accuracy, developer-level detail, per-page length |
| Screenshot auditor | Ran the devel app on a copy of the test project and compared all 21 doc images |
| MitoPilot user | Side-by-side comparison of the two sites; transition friction; what to borrow |
| Web designer | Fresh site build viewed at desktop and phone widths; layout, images, broken links |

Findings below are merged and de-duplicated. "Reviewers" counts how many raised the same point independently. I re-checked the items marked "confirmed" myself.

## Overall verdict

The site is in good shape. It is already shorter and clearer than the MitoPilot site (about 2,900 source lines versus 5,100), cross-links are clean, and almost every number and button name checked against the code was correct. No new pages are needed.

The work falls into five buckets:

1. A handful of quick, clear errors (published to-do list, broken link, one paragraph that contradicts its own table).
2. Screenshots lag the app: 13 of 21 are outdated.
3. Three pages carry most of the bloat: Sample metadata, Quality checks and flags, and Reference databases.
4. The same content is repeated on several pages.
5. Vocabulary: a few words carry two or three meanings.

Estimated effect of the proposed cuts: roughly 35 to 40 percent off the two longest pages, 10 to 30 percent off most others, and about 45 lines of additions in total.

## Urgent, not a docs issue

`RiboPilot_test.R` at the repo root contains a GitHub personal access token, and that file is committed and present on `origin/devel` (confirmed; value not printed anywhere in this review). The screenshot auditor also reported an NCBI API key in the same file; I did not confirm that one. Recommend revoking the token now, then removing it from the file. Removing it from history is a separate decision.

## Priority 1: quick fixes (errors)

| # | Where | Problem | Proposed change | Reviewers |
|---|---|---|---|---|
| 1.1 | Sample-Metadata, "Screenshots to add" (line 585 to end) | Internal to-do list is published and shows in the page contents | Delete the section | 5, confirmed |
| 1.2 | NEWS.md lines 10 to 11 | "Upgrading RiboPilot" link is split across two lines, so it points to `Instal%20lation.html` | Put the URL on one line | 1, confirmed |
| 1.3 | Curation-Flags, "Rule 7 deserves a special note" | Paragraph talks about low identity and taxonomic distance, but rule 7 is "disagrees with the sample's other contigs" | Replace with: "Rule 7 only flags: one contig disagreeing with the others is a hint, not proof." | 3, confirmed |
| 1.4 | Step bar on six pages | Says "4. Annotate"; app tab, page title, and menu say "Annotation Review" | Change the step bar label | 3, confirmed |
| 1.5 | Troubleshooting, step-name table | Lists `find_rrna`; real step names are `find_rrna_scan` and `find_rrna_pick` | Fix the row (or collapse the table, see 3.6) | 2 |
| 1.6 | Troubleshooting, sample run report | Shows 41 minutes, 4 samples, version 9018; rest of site says 16 minutes and 12 samples | Refresh, or trim to the lines being explained | 3 |
| 1.7 | Installation | "Each image takes about 7.7 GB" conflicts with 2.4 GB for Singularity on the Caches page | Say 7.7 GB is Docker only | 1 |
| 1.8 | Get started | "GenBank's Sequence Read Archive" | "NCBI's Sequence Read Archive" | 1 |
| 1.9 | Test-Project-Assemble, update-window caption | Says "with the list of samples to run"; window shows only a count | Reword caption | 1 |
| 1.10 | Citing line 50 | "Annotate Options"; check against the app label | Match the app | 2 |

## Priority 2: screenshots

Checked live against the devel app. 13 outdated, 4 minor drift, 4 match. No orphan or missing image files. The page text is mostly current; the pictures lag.

Three app changes explain most of it: "scaffold" became "contig"; the Metadata column, Metadata button, and Run Reports button were added; and the run window, chimera panel, and Export window were reworked.

| Image | Verdict | Main difference |
|---|---|---|
| taxonomy-locked | Outdated | Run Reports button, Metadata button and column, new ID columns |
| assemble-table | Outdated | Run Reports, Metadata button and column |
| assemble-details-chim01 | Outdated | "Scaffolds" now "Contigs"; chip colors changed |
| assemble-chimera-chim01 | Outdated | Chimera panel reworked; image is also cut off mid-plot |
| failed-sample | Outdated | Reworded message; "BLAST evidence"; new Notes panel |
| run-progress | Outdated | Buttons now "Stop" and "Stop and Close" |
| run-finished | Outdated | New report layout; log path now under `.runs/` |
| troubleshooting/failure-report | Outdated | Same shot as run-finished; picture disagrees with the text block under it |
| annotation-table | Outdated | "Contig #", "Exclude & Lock", new Trim and Restore buttons |
| annotation-details | Outdated | "contig 1" in title, "Exclude & Lock", new Copy sequence button |
| export-window | Outdated | Sections reordered |
| curation/chimera-severe, chimera-mild | Outdated | Chimera panel reworked; plots cut off |
| reference-propose | Minor | Record count changed; crop clips last row and Apply button |
| reference-ready | Minor | Run Reports button missing |
| update-window | Minor | Command text was hand-shortened; app prints more |
| rescue-baited | Minor | One help line reworded |
| export-flagged-prompt | Minor | "GenBank hit" now "BLAST hit"; "scaffold" now "contig" |
| reference-setup, export-assign-group, export-done | Match | export-done shows a different folder name than the other images |

Proposed changes:

- Retake the 13 outdated images in one session from a clean test project at the path the text uses (`~/ribopilot_test`), so paths and the version number are consistent.
- Capture the wide tables at a narrower window (about 1,100 px) or crop them; at 1,600 px the table text is unreadable once shrunk to page width.
- Crop dialog shots to the dialog only (several show a dark strip of the dimmed app).
- Drop two images rather than retake them: `failure-report.png` (reuse `run-finished.png`) and `chimera-severe.png` (same view as `assemble-chimera-chim01.png`). `export-assign-group.png` is a single text box and could also go.
- Add at most one new image: one app screenshot on Get started (it has none, and MitoPilot shows the app in the first screen).
- Fresh comparison captures are in the session scratch folder (`scratchpad/screenshots/NEW_*.png`). They are evidence, not drop-in replacements (long paths, some cropped).

Side note from the auditor: the cached GenBank counts file on this workstation (`~/.local/share/R/RiboPilot/ref_dbs/counts.tsv`, dated Sep 1) makes Propose groups return one group instead of the three the docs describe. Worth clearing before retaking the Reference screenshots.

## Priority 3: cut bloat on three pages

### 3.1 Sample metadata (597 lines, 22 screens, longest page, no images)

All reviewers flagged it. Target: under 250 visible lines, one page (do not split).

- Delete "Screenshots to add" (1.1).
- "Matching a voucher to a GBIF record": replace the three-query walkthrough with two sentences ("links only when exactly one specimen of the right species matches; otherwise paste the gbifID").
- "Where the link is read" and "Rules that keep linking safe": cut or collapse.
- Merge the three "combinations" tables into one.
- Message table: keep only rows where the user must act.
- NCBI API key paragraph: replace with a one-line link to Installation.
- Match the code style of the rest of the site (`RiboPilot::` prefix, full paths).

Two reviewers suggested splitting this page in two. I recommend against it: cutting achieves the same and avoids a new page.

### 3.2 Quality checks and flags (338 lines, 17 screens)

Keep: the three tables (warnings, Checks column, eight rules), chimera levels, and ruleset ranges. Cut or collapse:

- Auto-triage: the "two evidence channels" paragraphs, the rule 2 and rule 6 paragraphs (fold one line into each table cell), the "no-op" and "settings live with the project" sentences. About half the section.
- "Three details worth knowing" and "Depth is shown, not decided": keep first and last sentence.
- "Adding tuning for a new group": replace the pointer to an internal function with "Want ranges for your group? Open an issue."
- Variant section: two sentences per term; rewrite the garbled opening.

### 3.3 Reference databases (124 lines)

- "Setting a group up": cut each of the four bullets to two or three lines.
- "Good to know": repeats earlier points; keep only the 5.8S item.
- Say once: "A new project has no usable database until you build one here."

### 3.4 Smaller trims

| Page | Trim |
|---|---|
| Home (README) | Step 1 is three times longer than the other steps; cut download sizes. Move Installation and Getting started up under the intro. Open scope with "RiboPilot is built for animals." |
| Get started | "Samples, contigs, and records": cut to two sentences. "Open the app" filter paragraph: cut to one or two lines. Drop the MitoPilot history of the `_s` suffix. |
| Test-Project-Assemble | Checks bullet: four short lines plus the existing link |
| Your-Own-Project | Reserved-names list: replace with "RiboPilot stops and names any column it reserves" |
| Troubleshooting | Collapse the step-name table; shorten the `contaminant_unchecked` row to one sentence; "steps that run out of memory are retried up to three times" instead of the named list |
| Caches-and-Offline | Drop the older-versions row and the lock-file row; merge the three container rows |
| Scripting | Move "Build a reference database" to Without-Containers (it needs that route anyway); drop the run-lock and file-path details |
| Clusters | "The app needs the container too": one sentence |
| Installation | NCBI key: keep the `.Renviron` route only. Drop the `NXF_VER` sentence if RiboPilot handles it. |

### 3.5 Repeated content: keep one copy, link from the rest

| Content | Appears on | Keep in |
|---|---|---|
| "No reference database ships" and the fish warning | 6 pages | Reference-Databases, plus one line in Get started |
| Export folder tree | Test-Project-Export, GenBank-Submission | GenBank-Submission |
| Marks table (Reviewed, Exclude & Lock, Partial) | Test-Project-Annotate, GenBank-Submission | GenBank-Submission |
| Checks column colors | Test-Project-Assemble, Curation-Flags | Curation-Flags |
| Keep the tab open, Stop, Resume | Get started, Reference-Databases, Troubleshooting | Get started box and Troubleshooting |
| Hydra "Submit as job" boxes | 7 boxes across the walkthrough, plus Clusters | Clusters; one pointer box in Get started |
| NCBI API key | Installation, Sample-Metadata | Installation |
| Test run time and machine | Get started, Test-Project-Assemble | Get started |
| "NCBI runs ITSx" | 4 places | GenBank-Submission, with one NCBI source link |

### 3.6 Developer detail visible to users

- The 37 internal functions are hidden from the reference index but still get pages and appear in site search. Mark them internal so no pages are built.
- Move the three `*_normalize_*` helpers out of the visible "Sample metadata" group.
- `RiboPilot()` and `run_app()` are the same; list one. "WF1" appears in a function title.

## Priority 4: wording and accuracy

### 4.1 Statements that overpromise

| Where | Problem | Proposed wording |
|---|---|---|
| Curation-Flags and Test-Project-Assemble, green "ok" | "Matched your organism" reads as species confirmed; the check only tests the phylum | "ok: the best matches are in the same phylum as your sample, and nothing is flagged. This does not confirm the species." |
| Curation-Flags, paralog paragraph | Presents a paralog label as certainly real biology; the label only means the genes differ (below 99.5 percent) | Say it can be real variation, a second organism from the same phylum, or a pseudogene; BLAST each before locking |
| GenBank-Submission, "Flanking sequence" | "That is fine" is too reassuring for long unannotated flanks | Advise trimming to gene ends for the Portal unless submitting the flanks another way |
| GenBank-Submission, empty modifiers | "Take the modifier out" risks rejection if the field is required | "Fill in the column, or use one of NCBI's accepted missing-value terms" |
| Curation-Flags, tuning ranges | Presented as tuned per group; code calls them provisional | One sentence: working values, a prompt to look, not a verdict |

The NCBI-side claims (Portal uses ITSx, required modifiers, missing-value terms) were not checked against NCBI documentation by any reviewer. Verify before rewording.

### 4.2 Words with more than one meaning

| Word | Problem | Proposal |
|---|---|---|
| contig, unit, record, scaffold, piece | Five words for nearly one thing; "# Units" in Assemble uses "unit" differently than Get started | "Contig" until export, "record" after; define "unit" once as the app's table word; write `_s<contig>` in doc text |
| reference | The database you build, the closest GenBank sequence used in checks, and the "no reference" Checks value | Qualify each; explain "no reference" as "no match to judge by" |
| group | Taxonomic group, reference group, export group | Always say "reference group" or "export group" |
| off target vs contaminant | Look identical to a user | One sentence on the difference |
| database | Reference database, seed DB, bundled rDNA database, taxonomy copy | One sentence in Test-Project-Reference separating the one you build from the built-in one; mention once that the app calls it the seed database |
| isomer | Organelle vocabulary; confusing for rDNA | Describe as "same copy, assembled two ways" |
| nrDNA vs rDNA | Both used | Use rDNA; define once |
| headless (Clusters) | Used to mean "submits jobs" | "On a cluster, the run window offers Submit as job" |

### 4.3 Where a beginner stalls

- **Annotation Review walkthrough page** describes every button but never says what to do with the test project, and Export then assumes units are locked. Add a three-step "do this now" at the top, plus one sentence on which row to lock when a sample has two paths.
- **Installation**: link each tool name to its install page, say the check commands go in a terminal, and move the "Private preview" token box above the install command that needs it.
- **GenBank-Submission**: open with "Most users: upload `<group>.fasta` at the Submission Portal; ignore the .tbl files", and one sentence on what annotation is still for.
- **Required source fields**: the default header lacks voucher, date, and location. One sentence in Test-Project-Export and one in Your-Own-Project.
- **Test-Project-Reference**: give a time estimate beside the red "takes much longer" warning.
- **Rescue section**: one sentence explaining the baited assembler.
- **Get started**: list all sample states once.
- Walkthrough pages never mention the new Metadata column and button.

### 4.4 Cluster users

- Clusters never links to Caches-and-Offline. Add one line under shared storage.
- Offline steps: clarify when `NXF_VER` matters, and list the two downloads the steps currently omit (taxonomy, and reference sequences fetched during a run).
- `min_depth` (default 2,000,000 reads) is undocumented; a small read set fails every sample with no explanation. Two lines in Your-Own-Project.
- Say whether `executor = "slurm"` in `new_project()` is fine (MitoPilot users will type it).
- Hydra box: mention `Rproj = FALSE`, and where the container cache lands by default.

## Priority 5: site layout and appearance

Based on a fresh local build. The local `docs/` folder is stale (half the articles older than their sources, pages for removed functions); it is not tracked and the live site builds from CI, so this only affects local previews. The live site was not checked.

| # | Problem | Proposed fix | Effort |
|---|---|---|---|
| 5.1 | Table cells have very large padding. Tables cause much of the page length. | One rule in `extra.css` (smaller font and padding). Measured: Curation-Flags 28 percent shorter, Troubleshooting 39 percent, Sample metadata 22 percent. | 2 lines |
| 5.2 | Articles menu is capped in height with an inner scrollbar; 11 of 17 pages are hidden until you scroll inside the menu | Split into real top-level menus (below), or raise the cap | Config |
| 5.3 | Logo overlaps the header rule on every page; step bar wraps to two rows | Smaller logo on article pages (MitoPilot has this rule) | 3 lines |
| 5.4 | Screenshots have no border and captions look like body text | Border and caption rules (MitoPilot has these) | About 10 lines |
| 5.5 | Long code lines scroll sideways, worst on phone | Break long R calls at about 70 characters | Edits |
| 5.6 | "Get started" wraps in the navbar at mid widths; navbar is taller than needed | No-wrap rule; reduce nav height to about 72 px | 2 lines |
| 5.7 | Home page repeat diagram is cut off on phone | Shorten it or use a small figure | Small |
| 5.8 | Tables overflow on phone | Covered by 5.1 plus a phone rule | 2 lines |

Leave dark mode off; it would need about 20 extra rules and no switch is shown today.

### Proposed menu structure

Same pages, no additions; two reviewers proposed nearly the same thing.

- Get started
- Installation (direct link)
- Test project: 1 Taxonomy, 2 Reference, 3 Assemble, 4 Annotation Review, 5 Export
- Guides: Your own project, Reference databases, Quality checks and flags, GenBank submission, Sample metadata, Scripting (advanced)
- Clusters: Running on a cluster; Disk space, caches, and offline use
- Help: Troubleshooting, Citing, Installing without a container (unsupported), Changelog
- Reference

Keep the five walkthrough pages; they mirror the five app steps.

## MitoPilot comparison

**Can a MitoPilot user transition easily?** Mostly yes: same layout, same buttons, same cluster setup. But Get started says only "this will feel familiar", which undersells real differences that will trip habits:

- Two new steps come first, and Assemble stays empty until taxonomy is locked.
- `Taxon` now matters (matched to NCBI; drives the database, checks, and organism name).
- No database ships, fishes included.
- One run does assembly and annotation; Annotation Review has nothing to run.
- The app opens with `RiboPilot("path")`, not `setwd()` then `MitoPilot()`.
- Sample IDs: 18 characters, not 40. An existing mapping file can fail.
- Different words: contig, unit, step, keep or drop.
- GenBank upload is FASTA only.

**Proposal:** replace the "this will feel familiar" paragraph in Get started with a short "Coming from MitoPilot" table (about 14 lines net, no new page). Draft:

| | MitoPilot | RiboPilot |
|---|---|---|
| Steps | Assemble, Annotate, Export | Taxonomy, Reference, Assemble, Annotation Review, Export. A sample shows in Assemble only after its taxonomy is locked. |
| `Taxon` column | For your own use | Matched to NCBI. Picks the reference database and quality checks, and becomes the organism name in GenBank. |
| Reference database | Fish database ships | None ships. Build one per group in the Reference step. |
| Runs | One for Assemble, one for Annotate | One run does both. |
| Open the app | `setwd(project); MitoPilot()` | `RiboPilot("path/to/project")` |
| Your own assemblies | Separate function and project | An `Assembly` column in the mapping file; reads and assemblies can mix. |
| Sample IDs | Up to 40 characters | Up to 18 characters |
| Words | scaffold, assembly, module, ignore | contig, unit, step, keep or drop |
| Extra sequences | Join scaffolds, or ignore paths | Keep or drop contigs; nothing is joined. |
| GenBank upload | FASTA and feature table | FASTA only, through the Submission Portal |
| On Hydra | Call `hydra_setup()` first | The app does it for you. |

The "nothing is joined" and Hydra rows come from the RiboPilot docs, not from testing; confirm before publishing.

**Worth borrowing from MitoPilot** (all small):

| Item | Where | Size |
|---|---|---|
| Workflow overview figure | Home page, beside or replacing the five-step list | One image |
| Lock and state icon table | Get started, replacing prose | 8 lines |
| One-line descriptions under reference index groups | `_pkgdown.yml` | 5 lines |
| Screenshot border, caption style, smaller article logo | `extra.css` | About 13 lines |
| "Run a handful of samples first" tip | Your-Own-Project | 2 lines |
| Editing a shared parameter set affects every sample using it | Rescue section | 3 lines |
| Corrupted project database recovery | Troubleshooting, one row | 1 row |
| Download rows buttons (in the app, undocumented) | Get started | 1 line |

**Deliberately leave out:** FAQ page, separate Hydra and SEDNA pages, the long native-install page, ruleset browser, assembly-methods and difficult-assemblies pages, database inspection, and container R sessions. RiboPilot's single pages cover these better, or they do not apply.

**What RiboPilot already does better (keep):** short walkthrough pages with "What you should see" boxes; Installation with a check per tool and a clear upgrade procedure; the Caches tables; Troubleshooting with real messages; one Clusters page; dedicated GenBank and Citing pages; styles in one file.

## Suggested order of work

1. Revoke the committed token.
2. Priority 1 quick fixes (under an hour).
3. CSS and menu changes (5.1 to 5.4, 5.6): small, and they shorten every page.
4. Cuts and de-duplication (Priority 3), starting with Sample metadata and Quality checks.
5. Wording and accuracy (Priority 4), after checking the NCBI-side claims.
6. Coming from MitoPilot table and the small borrowed items.
7. Retake screenshots last, after wording and app labels settle.

## Decisions for you

- Rename user-facing "isomer", and the triage labels that say "GenBank match"? Both are app labels, so this touches the app, not only the docs.
- Keep Scripting and Without-Containers as visible pages, or tuck both under Help as advanced or unsupported?
- Hydra boxes in the walkthrough: collapse to one pointer, or keep for Smithsonian users?
- Menu restructure as proposed?

## Not covered by this review

- The live website (only sources and a fresh local build).
- NCBI-side claims, external links, and individual function help pages.
- Hydra and "Submit as job" screens; a clean 12-sample first run (screens were captured from a one-sample re-run).
- Dark mode (no switch exists), search box behavior, and real phones.

## Leftovers in the session scratch folder

Comparison captures, the fresh site build, and a 539 MB modified copy of the test project remain in the session scratchpad (`screenshots/`, `appearance/`, `mitopilot_compare/`). The original test project was not touched.
