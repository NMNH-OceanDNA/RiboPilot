# RiboPilot documentation review

Date: 2026-09-23. Docs and code reviewed: `main` at 3362308 (version 0.0.0.9018), read on branch `review/2026-09-full-review`. The pkgdown site was built from a copy of the repo; nothing was published.

This review covers the unpublished pkgdown site: the README (home page), the three vignettes (Get started, Reference databases, Curation flags), the function reference, `docker/README.md`, and the help text inside the app. Ten people worked on it:

- Seven docs reviewers, each reading as a different kind of user:
  - **docs-novice**: a field biologist or museum taxonomist who has never used a terminal.
  - **docs-intermediate**: an ecology grad student who uses R every day and plans 80 gastropods on a Slurm cluster.
  - **docs-hpc**: a research computing engineer installing for a lab (Slurm with Singularity, and a Docker workstation).
  - **docs-domain**: a molecular systematist who has submitted hundreds of rRNA and ITS records to GenBank.
  - **docs-editor**: a senior technical writer and documentation editor.
  - **docs-api**: an R package reference auditor.
  - **docs-appclaims**: a QA analyst checking that the docs describe the app and pipeline that ship.
- Three operators who followed the docs literally on real projects: **op-practice** (the practice project), **op-novel** (two gastropod genome skims from SRA), and **op-edge** (common user mistakes). Every pipeline run was started from the app's Update and Run from App buttons, never scripted.

The novice and intermediate reviews come from an earlier run of this review and were verified then; the other five docs reviews and the three operator runs come from this run. All were checked against the same code.

**Counts.** The reviewers reported 247 findings: 227 from the seven docs reviewers and 20 docs findings from the operators. After verification, 246 are confirmed and 1 is refuted (Appendix A); none are uncertain, disputed, or unverified-low. Merged, the 246 confirmed findings become **99 items: 5 blocker, 43 major, 40 minor, and 11 polish.**

| Section | Items | Blocker | Major | Minor | Polish |
|---|---|---|---|---|---|
| Missing content | 14 | 0 | 13 | 1 | 0 |
| All pages | 2 | 0 | 2 | 0 | 0 |
| README | 11 | 3 | 1 | 6 | 1 |
| Get started | 22 | 2 | 13 | 6 | 1 |
| Reference databases | 3 | 0 | 0 | 3 | 0 |
| Curation flags | 13 | 0 | 8 | 3 | 2 |
| Function reference | 17 | 0 | 5 | 10 | 2 |
| docker README | 4 | 0 | 0 | 2 | 2 |
| In-app help | 7 | 0 | 1 | 4 | 2 |
| Site and repo | 6 | 0 | 0 | 5 | 1 |
| **Total** | **99** | **5** | **43** | **40** | **11** |

How to read the labels:

- **Status.** Every finding in this review was sent to a separate verifier, who re-checked it against the code and often re-ran it. **Confirmed** means the verifier agrees. **Refuted** findings are listed only in Appendix A. No finding ended as uncertain or disputed, and none is **unverified-low** (a low finding not sent to a verifier), so every item below is confirmed.
- **Severity.** **Blocker**: a user following the docs cannot finish, or gets a wrong result without warning. **Major**: serious confusion, wrong or stale information, or a missing essential step. **Minor**: unclear or incomplete but recoverable. **Polish**: wording, style, or formatting. When reviewers or the verifier rated merged findings differently, the item takes the verifier's rating (or the highest confirmed rating) and says so.
- **Items.** Duplicate findings are merged into items DOC-01 to DOC-99. Each item lists every source id; "(part)" means only part of that finding belongs to the item. "Seen by" counts distinct reviewers. "See also" points to the sibling reviews in this folder: [USABILITY_REVIEW.md](USABILITY_REVIEW.md) (ids UX-nn) and [CODE_REVIEW.md](CODE_REVIEW.md) (ids like PKG-H1), which often cover the code side of the same problem.

## The short version

**Verdict: not ready to publish.** None of the seven docs reviewers could finish a real project from the docs alone. The three operators did finish, but only by following the app instead of the guide, and by working around problems the docs never mention. Most fixes are writing. A handful are small code changes that make the documented commands work.

1. **Get started describes an older app** ([DOC-28](#doc-28)). It talks about three tabs and runs the pipeline from the R console. The app has five numbered steps, starts at Taxonomy, and runs everything from its own buttons. Following the guide gives an empty Assemble table, and assemblies made with no taxonomy. Nine of the ten reviewers reported it.
2. **The first commands fail** ([DOC-17](#doc-17), [DOC-18](#doc-18), [DOC-30](#doc-30), [DOC-31](#doc-31), [DOC-94](#doc-94)). `install_github()` cannot find the Bioconductor packages. `new_test_project()` and `new_project()` crash outside RStudio. `RiboPilot("~/my_project")` needs `library()`. The API key placeholder breaks every NCBI call. Code inside the alert boxes renders with curly quotes, so pasting it fails.
3. **While the repo is private, the defaults do not work** ([DOC-19](#doc-19)). Install and the practice reads need `GITHUB_PAT`, which no doc mentions. The built-in fish reference returns 404 even with a token, so the practice project's first run aborts while the app says "Done".
4. **Whole topics are missing** ([DOC-01](#doc-01), [DOC-03](#doc-03), [DOC-04](#doc-04), [DOC-05](#doc-05), [DOC-08](#doc-08), [DOC-09](#doc-09), [DOC-12](#doc-12)): what to install (Docker or Singularity, Java, a Nextflow version), hardware, cluster use, caches, troubleshooting, the mapping-file rules, NEWS, and citation.
5. **What goes to GenBank is not explained, and parts of it are wrong** ([DOC-42](#doc-42), [DOC-44](#doc-44), [DOC-45](#doc-45), [DOC-43](#doc-43), [DOC-39](#doc-39)). Partial versus complete comes from scaffold position, not gene completeness. NCBI's required source modifiers are not listed. Isolate is banned in templates yet added automatically. `[organism=]` uses the raw mapping label. Unannotated flanks are exported, and ITS spacers stop 1 to 4 bp short of the genes.
6. **The Curation guide is out of date, and two scientific statements are wrong** ([DOC-53](#doc-53) to [DOC-60](#doc-60)). It documents 2 of 5 rulesets and 7 of 8 triage rules, misses flags users see, and has no chimera section. The isomer test ignores the ITS, and "cannot be contamination" contradicts the phylum-level screen described on the same page. The "GenBank screen" is, by default, a local panel, so a methods section written from the docs would be wrong.
7. **The function reference is not ready** ([DOC-66](#doc-66), [DOC-67](#doc-67), [DOC-68](#doc-68), [DOC-69](#doc-69), [DOC-72](#doc-72)): 54 topics in one flat list, most of them pipeline internals, raw markdown on every page, and key pages that give wrong advice (for example, the `RiboPilot()` page never says that `launch.browser = FALSE` hides the Run from App button).

There is also plenty to keep: the Reference guide matches the app closely, the writing is plain and friendly, the version claims check out, the site builds cleanly, and the app's evidence views make excellent screenshots. See [What works well](#what-works-well-and-should-be-kept).

**Suggested order of work:**

1. Small code fixes that make the documented commands work: an empty `biocViews:` line ([DOC-17](#doc-17)), the RStudio check ([DOC-18](#doc-18)), `executor = "docker"` as the default ([DOC-16](#doc-16)), the Singularity rawDir bind ([DOC-29](#doc-29)), and gating assembly on the taxonomy lock ([DOC-20](#doc-20)).
2. A release gate: make the repo public, and ship the fish seeds as release assets or inside the image, before the site goes live ([DOC-19](#doc-19)).
3. Rewrite Get started around the app's five steps, using the screenshots in this folder ([DOC-28](#doc-28); [Screenshots for the docs](#screenshots-for-the-docs)).
4. Write the missing pages in the [proposed outline](#proposed-outline-for-the-published-site): Installation, Your own project (mapping file), Preparing a GenBank submission, Running on a cluster, Troubleshooting, Glossary, NEWS, and CITATION.
5. Update the Curation guide and group the function reference.

## Verdict by persona

"Stops" means the first place the reviewer could not go on from the docs alone.

| Persona | Could they finish a real project from the docs? | Where they stop first | Where they stop next | Main items |
|---|---|---|---|---|
| docs-novice: field biologist, never used a terminal | No | Installation. Docker, Java, supported computers, and hardware are never covered, and the only concrete recipe is the unsupported conda one. | `install_github()` (Bioconductor), then `new_test_project()` ("RStudio not running"), then an empty Assemble table, because Get started never mentions Taxonomy. On a 16 GB laptop the run would also stop at the 24 GB memory request. | [DOC-17](#doc-17) [DOC-18](#doc-18) [DOC-28](#doc-28) [DOC-01](#doc-01) [DOC-02](#doc-02) [DOC-03](#doc-03) |
| docs-intermediate: grad student, 80 gastropods on Slurm | No | Get started Steps 5 and 6. `RiboPilot()` needs `library()`, the Assemble tab is empty, and the assemblies already ran with no lineage against fish seeds. | On the cluster: nothing on opening the app, the job script (2 CPUs, 8 GB, 48 h), or caches. At export: raw Taxon in `[organism=]`, isolate refused, `export_project()` errors. No way to add a second batch, and nothing to cite. | [DOC-28](#doc-28) [DOC-30](#doc-30) [DOC-04](#doc-04) [DOC-05](#doc-05) [DOC-45](#doc-45) [DOC-44](#doc-44) [DOC-40](#doc-40) [DOC-11](#doc-11) [DOC-12](#doc-12) |
| docs-hpc: research computing engineer | Partly. A full practice run from the app on Singularity finished cleanly (all 44 tasks, about 5.5 min for 3 samples), but only with steps the docs never mention: a cache folder off home, `NXF_VER`, an SSH tunnel, a hand-edited profile, and local seed files. | Install (Bioconductor). | The practice project (fish seeds 404), assembly-only samples (config parse error on generated Singularity profiles), and a firewalled cluster (had to read the code to learn the network needs). The first Singularity run died when the browser tab closed. | [DOC-17](#doc-17) [DOC-19](#doc-19) [DOC-29](#doc-29) [DOC-06](#doc-06) [DOC-04](#doc-04) [DOC-05](#doc-05) [DOC-36](#doc-36) |
| docs-domain: rDNA systematist and GenBank submitter | Partly: through assembly, triage, and annotation review. | Export. Would have to re-derive partial or complete for every record, fill 1 to 4 bp gaps at the ITS boundaries, decide how to trim ETS and IGS flanks, and rewrite the header for NCBI's required source modifiers. | Could not write an accurate methods paragraph: the "GenBank screen" was a local panel, and for a non-fish group the thresholds applied are recorded nowhere (2 of 5 rulesets documented). | [DOC-42](#doc-42) [DOC-39](#doc-39) [DOC-43](#doc-43) [DOC-44](#doc-44) [DOC-45](#doc-45) [DOC-60](#doc-60) [DOC-54](#doc-54) [DOC-58](#doc-58) [DOC-59](#doc-59) |
| docs-editor: technical writer | Not without reading the code | Install (Bioconductor, and `GITHUB_PAT` while private). | The practice project ("RStudio not running"), then Steps 4 to 6 (Taxonomy missing, empty Assemble table). With the code open, finished through the app. | [DOC-17](#doc-17) [DOC-19](#doc-19) [DOC-18](#doc-18) [DOC-28](#doc-28) |
| docs-api: R reference auditor | No, not from a terminal | The first command (the `Rproj` crash), or install (no `biocViews`). | The `new_project()` help sends Hydra users to `generate_config()`; empty Assemble tab; with `launch.browser = FALSE` the Run button is gone and no page says why; "No export groups assigned"; the useful reference pages are buried among internals. | [DOC-18](#doc-18) [DOC-17](#doc-17) [DOC-67](#doc-67) [DOC-69](#doc-69) [DOC-40](#doc-40) [DOC-66](#doc-66) |
| docs-appclaims: QA analyst | Yes, but only because the app's own empty-table messages rescued them | Step 6. The app opens on Taxonomy, the Assemble table is empty, and the promised endpoint controls are not in Assemble. | Step 9. `export_project()` errors, and after grouping in the app a Problematic unit silently disappears from the bundle. | [DOC-28](#doc-28) [DOC-38](#doc-38) [DOC-40](#doc-40) [DOC-41](#doc-41) |
| op-practice: first-time user, practice project | Finished end to end, with workarounds the docs do not give | `new_test_project()` crashed outside RStudio (worked around with `Rproj = FALSE`). | The first run aborted on the 404 fish reference while the app said "Done" (worked around with local seed files that real users do not have). Learned the Taxonomy and Reference steps from the app. `export_project()` needs groups. Results: 28S fragments exported as "complete sequence", and a clean flounder exported as "chimeric". | [DOC-18](#doc-18) [DOC-19](#doc-19) [DOC-28](#doc-28) [DOC-40](#doc-40) [DOC-42](#doc-42) [DOC-57](#doc-57) |
| op-novel: gastropod grad student, SRA genome skims | Finished, by following the app rather than the guide (only `new_project()` in R) | Creating the project took three tries (the Rproj crash, then R1/R2 in subfolders, which the docs do not cover). | The rescue for a failed sample was found only by exploring Assembly Options; `no_clade_ruleset` on every unit, unexplained; a sample whose rDNA is a sea anemone exported under the snail's name (UX-31 in the usability review). | [DOC-18](#doc-18) [DOC-09](#doc-09) [DOC-37](#doc-37) [DOC-53](#doc-53) [DOC-54](#doc-54) |
| op-edge: support engineer reproducing user mistakes | Not applicable (tested error handling) | No mapping-file rules or troubleshooting page to point users to; many setup mistakes give raw R errors. | Nothing warns that closing or refreshing the tab stops a run, conda is still the default executor, and "Already have an assembly?" says to leave R1 and R2 blank but the columns must exist. | [DOC-09](#doc-09) [DOC-08](#doc-08) [DOC-36](#doc-36) [DOC-16](#doc-16) |

## Top problems to fix before publishing

Ranking rule: severity first, then whether the problem stops users or puts wrong data into GenBank, then how many reviewers hit it. "Fix" says whether the fix is in the docs, the code, or both.

| Rank | Item(s) | Problem | Severity | Seen by | Fix |
|---|---|---|---|---|---|
| 1 | [DOC-28](#doc-28) | Get started describes an older, console-driven app with three tabs | blocker | 9 | docs |
| 2 | [DOC-18](#doc-18) | The documented new_project() and new_test_project() calls crash outside RStudio | blocker | 4 | code + docs |
| 3 | [DOC-17](#doc-17) | The install command fails on a fresh R (no biocViews) | blocker | 5 | code, one line |
| 4 | [DOC-19](#doc-19) | Private repo: GITHUB_PAT never mentioned, and the fish reference 404s even with a token | blocker | 6 | release gate + code |
| 5 | [DOC-29](#doc-29) | "Already have an assembly?" fails on generated Singularity profiles | blocker | 1 | code |
| 6 | [DOC-30](#doc-30), [DOC-31](#doc-31), [DOC-94](#doc-94) | Copy-paste traps: RiboPilot() without the prefix, the API key placeholder, curly quotes in code | major | 5 | docs |
| 7 | [DOC-01](#doc-01), [DOC-02](#doc-02), [DOC-03](#doc-03) | Prerequisites, supported computers, and hardware are not documented | major | 4 | docs |
| 8 | [DOC-16](#doc-16) | conda is the main install recipe and still the default executor | major | 3 | code + docs |
| 9 | [DOC-42](#doc-42), [DOC-44](#doc-44), [DOC-45](#doc-45), [DOC-43](#doc-43), [DOC-39](#doc-39) | What goes to GenBank is not explained, and parts of it are wrong | major | 3 | code + docs |
| 10 | [DOC-40](#doc-40), [DOC-41](#doc-41) | export_project() needs app-made groups; Problematic silently drops units | major | 6 | docs (+ code) |
| 11 | [DOC-53](#doc-53), [DOC-54](#doc-54), [DOC-55](#doc-55), [DOC-56](#doc-56), [DOC-60](#doc-60) | The Curation guide is out of date (flags, rulesets, eight rules, chimeras), and the "GenBank screen" is a local panel | major | 6 | docs |
| 12 | [DOC-58](#doc-58), [DOC-59](#doc-59), [DOC-57](#doc-57) | Wrong scientific guidance on variants and contamination; the chimera false positive | major | 1 | docs + code |
| 13 | [DOC-36](#doc-36) | Closing the tab or run window stops a run; the docs say you can keep working | major | 2 | docs (+ code) |
| 14 | [DOC-04](#doc-04), [DOC-05](#doc-05), [DOC-06](#doc-06), [DOC-69](#doc-69) | No cluster, cache, or network guidance, and the headless rule is described wrongly | major | 4 | docs |
| 15 | [DOC-66](#doc-66), [DOC-67](#doc-67), [DOC-68](#doc-68), [DOC-72](#doc-72) | The function reference: flat list of internals, wrong pages, raw markdown | major | 6 | docs |
| 16 | [DOC-08](#doc-08), [DOC-09](#doc-09) | No troubleshooting page and no mapping-file rules | major | 5 | docs |
| 17 | [DOC-15](#doc-15), [DOC-13](#doc-13) | Docs and app use different names, and there are no screenshots | major | 4 | docs |
| 18 | [DOC-12](#doc-12) | No NEWS, citation, or tool-version list | major | 3 | docs |

## Claim accuracy

The docs-appclaims reviewer checked every claim the README and vignettes make about the app and pipeline; the other reviewers checked claims in their areas. This section lists the statements that are wrong or stale compared with the app and code that ship, then the ones that checked out, so they can be left alone.

### Wrong or stale statements

| # | Claim (as written) | Where | What ships | Item |
|---|---|---|---|---|
| 1 | "samples must be reviewed and locked here before they move on to assembly" | `README.Rmd:44` | The pipeline and Update ignore the taxonomy lock; only the Assemble table filters on it. | [DOC-20](#doc-20) |
| 2 | "comes ready to use for them out of the box" (ray-finned fishes) | `README.Rmd:54` | The default seeds and labels are raw GitHub URLs that return 404 while the repo is private; no token helps. | [DOC-19](#doc-19) |
| 3 | "Checks currently ship for a generic metazoan (animal) setting and a fish-specific setting" | `README.Rmd:56` | Five rulesets ship, picked from each sample's lineage. | [DOC-54](#doc-54) |
| 4 | "The taxonomy snapshot (~72 MB) ... takes about a minute" | `README.Rmd:44` | About 80 MB plus a 400 MB index per project; the app itself says "a few minutes". | [DOC-24](#doc-24) |
| 5 | "RiboPilot refuses to load on an older Nextflow" | `README.Rmd:77` | Only with `library()`; the documented `RiboPilot::` calls check at run or app start. | [DOC-25](#doc-25) |
| 6 | "Run this once per session, before launching the app" (`hydra_setup()`) | `README.Rmd:111`; `man/hydra_setup.Rd` | `run_app()` applies it automatically. | [DOC-21](#doc-21) |
| 7 | The LSU label in the repeat diagram | `README.Rmd:30-34` | Sits under 5.8S, not under 28S. | [DOC-23](#doc-23) |
| 8 | "so setup takes seconds" | `README.md:173-175` only | True for the read subsets, but the first run also pulls a 2.6 GB image, and the paragraph is missing from README.Rmd. | [DOC-22](#doc-22) |
| 9 | "the pin here is the same version, so the conda path and the container path agree" | `environment.yml:1-8` | barrnap 0.9 vs 1.10.5, entrez-direct 24.0 vs 26.0, datasets 18.34.0 vs 18.36.0. | [DOC-26](#doc-26) |
| 10 | "A project moves through four stages ... Assemble ... Annotate ... Review and curate ... Export" | `RiboPilot.Rmd:51-57` | Five app steps: Taxonomy, Reference, Assemble, Annotation Review, Export. | [DOC-28](#doc-28) |
| 11 | "organized into three tabs that follow the workflow: Assemble, Annotate, and Export" | `RiboPilot.Rmd:294-296` | Five numbered steps, opening on Taxonomy; the fourth is "Annotation Review". | [DOC-28](#doc-28) |
| 12 | "the spacers are simply the stretches between those landmarks" | `RiboPilot.Rmd:47-48` | 1 bp and 4 bp are left unannotated between genes and spacers. | [DOC-39](#doc-39) |
| 13 | "In a container, build the same reference data from the app's Reference tab instead" (Step 2) | `RiboPilot.Rmd:169-173` | The Reference tab needs a project (Step 3) and assigned taxonomy first. | [DOC-32](#doc-32) |
| 14 | "If your samples are fish, you can skip to Step 3" | `RiboPilot.Rmd:143-144` | Only ray-finned fishes get the fish seeds and ruleset. | [DOC-33](#doc-33) |
| 15 | `ncbi_api_key = "PASTE_YOUR_KEY_HERE"` | `RiboPilot.Rmd:195` | Passed to NCBI verbatim; NCBI answers HTTP 400 "API key invalid". | [DOC-31](#doc-31) |
| 16 | "You can let it run and come back later" | `RiboPilot.Rmd:258-259` | True only from the console; an app run stops when the tab closes. | [DOC-36](#doc-36) |
| 17 | "If NCBI is unreachable or slow, RiboPilot simply skips that lookup ... not affected" | `RiboPilot.Rmd:264-268` | The default searches the bundled panel first; a missing result adds contaminant_unchecked and a triage flag. | [DOC-60](#doc-60) |
| 18 | "a second attempt at just those samples" | `RiboPilot.Rmd:272-276` | `run_fallback()` also re-assembles good samples that carry one fragment scaffold. | [DOC-37](#doc-37) |
| 19 | R1 and R2 are "the two file names" | `RiboPilot.Rmd:182` | Paths relative to `data_path`; subfolders work and are the usual SRA layout. | [DOC-09](#doc-09) |
| 20 | "Those samples do not need reads: leave their R1 and R2 blank" | `RiboPilot.Rmd:228` | Blank values work, but the R1 and R2 columns must exist, and generated Singularity profiles cannot run such a project at all. | [DOC-09](#doc-09), [DOC-29](#doc-29) |
| 21 | `RiboPilot("~/my_project")` | `RiboPilot.Rmd:291` | Fails: no page calls `library(RiboPilot)`. | [DOC-30](#doc-30) |
| 22 | "switch to an alternative reconstruction or adjust the exact endpoints" (in Assemble) | `RiboPilot.Rmd:303-305` | Assemble has only Keep and Drop; trimming is in Annotation Review. | [DOC-38](#doc-38) |
| 23 | "that sample has to be annotated again. Press Update on the Assemble tab" | `RiboPilot.Rmd:318-322` | A trim re-curates at once, and Keep/Drop needs no re-annotation. | [DOC-38](#doc-38) |
| 24 | "a Curation flags panel ... Each flag links to the feature it concerns" | `RiboPilot.Rmd:332-337` | The panel is "Curation warnings", and only feature-level warnings have select links. | [DOC-15](#doc-15), [DOC-53](#doc-53) |
| 25 | "or export everything from R: `export_project()`" | `RiboPilot.Rmd:342-347` | Stops unless export groups exist (made only in the app); uses the default template; skips Problematic units. | [DOC-40](#doc-40) |
| 26 | "markers/ ... (18S, ITS1, 5.8S, ITS2, 28S, ITS)" | `RiboPilot.Rmd:351-362` | The 5.8S folder is `58S`, and `export_project()` also writes a zip. | [DOC-47](#doc-47) |
| 27 | "Before the bundle is written, RiboPilot runs stock ITSx" | `RiboPilot.Rmd:364-365` | It runs after the FASTA and `.tbl` are written. | [DOC-47](#doc-47) |
| 28 | "The template is the whole defline" | `RiboPilot.Rmd:382-383` | `[isolate=...]` is appended for multi-record samples, and any template containing "isolate" is refused. | [DOC-44](#doc-44) |
| 29 | `[specimen_voucher={Voucher}]`, `[collection_date={collection_date}]` | `RiboPilot.Rmd:382-383` | NCBI spells these specimen-voucher and collection-date. | [DOC-44](#doc-44) |
| 30 | "Locking a sample signs off on it" | `RiboPilot.Rmd:297-300` | Locked units can still be edited, trimmed, and restored. | [DOC-90](#doc-90) |
| 31 | `RiboPilot(project)   # Steps 6-8 ...` | `RiboPilot.Rmd:391` | `project` is never defined, and the Taxonomy and Reference steps are left out. | [DOC-30](#doc-30), [DOC-28](#doc-28) |
| 32 | "The default fish databases are fetched from GitHub at run time" | `docker/README.md:107-109` | The URLs return 404 while the repo is private, even with a token. | [DOC-19](#doc-19) |
| 33 | "without letting it grow past about 25,000 records" | `Reference-Databases.Rmd:41-44` | The first clade with at least 1,000 records is kept even above 25,000. | [DOC-50](#doc-50) |
| 34 | "Each row is one proposed group" | `Reference-Databases.Rmd:47-50` | The tab opens on a Samples view; Propose groups needs Apply; there is a sixth status, "Database missing". | [DOC-51](#doc-51) |
| 35 | "Press Run" | `Reference-Databases.Rmd:58-61` | The button is "Run from App". | [DOC-51](#doc-51) |
| 36 | "four choices" (Set up database) | `Reference-Databases.Rmd:55` | Three when several groups are set up at once. | [DOC-51](#doc-51) |
| 37 | "you can keep working in the app while it runs" | `Reference-Databases.Rmd:61-62` | The build window covers the app, and closing it or the tab stops the build. | [DOC-36](#doc-36) |
| 38 | "If you are using RiboPilot through a remote browser" | `Reference-Databases.Rmd:63-64` | The rule is RStudio Server, a scheduler on PATH, or no browser launch. | [DOC-69](#doc-69) |
| 39 | "from the broadest rank down to the narrowest" | `Reference-Databases.Rmd:73-75` | The ladder is "Shared lineage (deepest first)". | [DOC-51](#doc-51) |
| 40 | "Leave these samples on the built-in fish database. Nothing is built." | `Reference-Databases.Rmd:82-84` | A group already wired to a clade database stays on it. | [DOC-51](#doc-51), [DOC-88](#doc-88) |
| 41 | "Working through the Reference tab is optional" | `Reference-Databases.Rmd:91` | Effectively required outside ray-finned fishes. | [DOC-52](#doc-52) |
| 42 | "its quality checks anchor on the 18S and 28S genes" | `Reference-Databases.Rmd:191` | 5.8S is checked too. | [DOC-93](#doc-93) |
| 43 | The flag table (7 codes), and "those two checks are simply skipped" | `Curation-Flags.Rmd:32-45` | Also truncated, contaminant_unchecked, not_nuclear_rdna, and no_clade_ruleset; a missing lookup raises contaminant_unchecked. | [DOC-53](#doc-53) |
| 44 | "a small count next to each sample in the table" | `Curation-Flags.Rmd:26-28` | One count per unit, and it counts only sample-level flags. | [DOC-53](#doc-53) |
| 45 | "seven named rules ... the other three only flag" | `Curation-Flags.Rmd:69-73` | Eight rules; four only flag. | [DOC-55](#doc-55) |
| 46 | "the **GenBank screen**, the best reference match found remotely"; "searched against all of GenBank" | `Curation-Flags.Rmd:55-57`; `R/constants.R:150` | By default a bundled local panel plus decoy sets; NCBI only when the panel finds nothing. No practice unit was searched remotely. | [DOC-60](#doc-60) |
| 47 | "Picking between isomers is a decision that does not matter either way" | `Curation-Flags.Rmd:127-135` | Isomers are judged on the anchors only, so ITS variants are called isomers. | [DOC-58](#doc-58) |
| 48 | "the difference between them cannot be contamination" | `Curation-Flags.Rmd:120-124` | The screen works at phylum rank; same-phylum contaminants pass. | [DOC-59](#doc-59) |
| 49 | "Low depth can just as easily mean one molecule's reads were split across such paths" | `Curation-Flags.Rmd:161-170` | Reads are mapped per path and multi-mappers are dropped, not split. | [DOC-61](#doc-61) |
| 50 | The rule 7 note about "a low-identity hit at a large taxonomic distance" | `Curation-Flags.Rmd:113-116` | Rule 7 is disagreement with sibling scaffolds. | [DOC-62](#doc-62) |
| 51 | "RiboPilot ships two settings", and fish spacers "kept wide" | `Curation-Flags.Rmd:174-180` | Five rulesets; the fish ITS range is narrower than the generic one. | [DOC-54](#doc-54) |
| 52 | "28S 2000-6000 (the recovered portion varies with sequencing depth)" | `Curation-Flags.Rmd:172-190` | Complete 28S length does not depend on depth. | [DOC-54](#doc-54) |
| 53 | `Rproj`: "no effect outside RStudio" | `man/new_test_project.Rd`, `man/new_project.Rd` | Crashes outside RStudio. | [DOC-18](#doc-18) |
| 54 | "NMNH_Hydra", "NOAA_SEDNA": "a template, not a usable executor"; conda "Ready to use as-is" | `man/new_project.Rd` | Hydra and SEDNA work directly; conda is unsupported. | [DOC-67](#doc-67) |
| 55 | "Annotation and curation are handled by re-running this workflow after you review and lock" | `man/run_assembly.Rd` | Annotation runs in the same run. | [DOC-68](#doc-68) |
| 56 | Headless: "Detected automatically (RStudio Server or a batch scheduler on PATH means yes)" | `man/RiboPilot.Rd:19-23` | Also headless whenever `launch.browser` is not TRUE, which hides Run from App. | [DOC-69](#doc-69) |
| 57 | "the bundled 10-sample test set" | `man/new_test_project.Rd` | 12 samples: 10 real and 2 synthetic. | [DOC-35](#doc-35) |
| 58 | "zip: Also write `<out_dir>/<label>.zip`" | `man/export_bundle.Rd` | Written to `<out_dir>/export/<label>.zip`, and not at all with a relative path. | [DOC-47](#doc-47) |
| 59 | "auto picks docker for local/awsbatch ... or singularity/apptainer/docker" | `man/generate_config.Rd` | The engine is ignored for local, and `queue` for awsbatch. | [DOC-71](#doc-71) |
| 60 | "Annotation, curation, and export tables are added in later build increments" | `man/new_db.Rd` | They already exist. | [DOC-78](#doc-78) |
| 61 | `container = "macguigand/ribopilot:v0.1.0"` | `docker/README.md:94` | Tags are bare versions, for example 0.0.0.9018. | [DOC-83](#doc-83) |
| 62 | "add another `-v` mount to `docker.runOptions`" | `docker/README.md:101-103` | Mounts are added automatically. | [DOC-83](#doc-83) |
| 63 | Selected Path: "Assembly path chosen for annotation" | `R/constants.R:348` | Every kept scaffold of every path is annotated. | [DOC-88](#doc-88) |
| 64 | "Saving re-queues the selected samples: their state becomes Ready to run." (Annotate and Curate options) | `R/app_opts.R:42-44` | The sample becomes In progress while its annotation units are queued. | [DOC-88](#doc-88) |
| 65 | Problematic: "Flag this annotation for another look" | `R/app_server.R:1931` | It also removes the unit from every export. | [DOC-41](#doc-41) |
| 66 | "Min anchor identity (%)" and "Flag contaminants" in Curation options | `R/app_opts.R:302-307` | Not read by the engine. | [DOC-87](#doc-87) |
| 67 | Chimera panel: "A chimera is one scaffold stitched from two organisms ... Overlapping 500 bp windows" | `R/app_modals.R:104-108` | Four channels, and misjoins within one genome are flagged too. | [DOC-56](#doc-56) |
| 68 | README.txt: "itsx_ITS_too_long  misassembly" | `R/export.R:285` | Normal in mammals and birds. | [DOC-48](#doc-48) |
| 69 | README.txt: "NCBI's ITSx evaluation may silently remove sequences" | `R/export.R:263` | NCBI says "You will be alerted to any changes that are made to your sequences". | [DOC-48](#doc-48) |
| 70 | README.txt: "chimeric: the chimera detector or ITSx flagged mixed origins" | `R/export.R:267` | A clean flatfish is flagged by a false positive, and the planted chimera gets the same score. | [DOC-57](#doc-57) |
| 71 | "Workflow run finished." after every run | `R/app_server.R:1017` | Shown even when every task failed. | [DOC-88](#doc-88) |
| 72 | The header Help link | `R/app_server.R:85` | Opens the GitHub repo, not the docs. | [DOC-89](#doc-89) |

### Statements that checked out

These are right as written and can be kept:

- R >= 4.4.0 (matches DESCRIPTION). Nextflow 24.10.x through 25.10.x (matches `R/nextflow_version.R`), with newer versions pinned through `NXF_VER` to 25.10.6.
- `new_project()` picks the image tag from the package version, and there is no `:latest` tag (the tag policy in `docker/README.md`). Docker Hub has 0.0.0.9018 (amd64), and `singularity pull` works.
- "One command runs the whole analysis ... annotates every reconstruction it produces" and "There is no separate annotation command": one run does preprocess through triage.
- Locking gates: an Assemble lock fills Annotation Review, an annotation lock is needed for export, and several rows can be locked together.
- The Reference guide's labels (Propose groups, Settings, the status names, Set up database / Try again / Change, the four setup choices, "records unknown", clicking the Clade cell), the automatic wiring after a build, that locked samples are left alone, and that the app does not call `build_reference_db()`.
- The Curation guide's "Your own decisions stick" (and "Also re-decide scaffolds I set by hand"), "Depth is shown, not decided", the Screen and Variant columns, and that switching every rule off makes Apply a no-op.
- The export bundle: every listed file exists, the ITSx check only warns, and the ITSx Warnings column is in Assemble and Annotation Review. Header templates are checked as you type and saved under a name.
- The GenBank screen really is phylum-level and cannot tell a human from a fish; the guide says so honestly.

## Missing content

The task named install prerequisites, troubleshooting, a glossary, NEWS, and citation. All of them are missing, along with several other topics users needed. The first table shows what exists today; the items that follow give the detail.

| Topic | What the docs have today | Item |
|---|---|---|
| Install prerequisites (Docker or Singularity, Java, Nextflow version) | R and Nextflow only | [DOC-01](#doc-01) |
| Supported operating systems | Nothing | [DOC-02](#doc-02) |
| Hardware, disk, and run time | Nothing | [DOC-03](#doc-03) |
| Running on a cluster | One sentence ("ask your computing support staff") | [DOC-04](#doc-04) |
| Where caches and big files go | "a default location on your computer" | [DOC-05](#doc-05) |
| Network needs and offline use | Nothing | [DOC-06](#doc-06) |
| Upgrading RiboPilot | Nothing | [DOC-07](#doc-07) |
| Troubleshooting | One conda-only box | [DOC-08](#doc-08) |
| Mapping-file rules | Two sentences | [DOC-09](#doc-09) |
| Project folder layout | The export tree only | [DOC-10](#doc-10) |
| Adding samples later | Nothing | [DOC-11](#doc-11) |
| NEWS, citation, and tool versions | Nothing | [DOC-12](#doc-12) |
| Screenshots and an app tour | Nothing (logo and badge only) | [DOC-13](#doc-13) |
| Glossary | Nothing | [DOC-14](#doc-14) |
| The Taxonomy step | One README sentence | [DOC-28](#doc-28) |
| Preparing a GenBank submission (modifiers, completeness, flanks, what to upload) | A few lines in Step 9 | [DOC-44](#doc-44), [DOC-42](#doc-42), [DOC-43](#doc-43), [DOC-46](#doc-46) |
| Chimera detection and submission risk | Nothing | [DOC-56](#doc-56) |
| A private-preview note (`GITHUB_PAT`) until the repo is public | Nothing | [DOC-19](#doc-19) |

<a id="doc-01"></a>**DOC-01. The install list leaves out Docker or Singularity, Java 17+, and a Nextflow version.** Major. Confirmed. ids: D-novice-8, D-intermediate-15, D-hpc-11, D-editor-6. Seen by 4 (docs-editor, docs-hpc, docs-intermediate, docs-novice).
- Where: `README.Rmd:74-78` ("RiboPilot needs three things installed on your computer"); `vignettes/RiboPilot.Rmd:59-119` (Step 0).
- Problem: The list is R, Nextflow, and "the bioinformatics tools, only if you are not running in a container". The supported route needs Docker (usable without sudo) or Singularity/Apptainer, plus Java 17+, and neither is on the list. Step 0 gives no Nextflow version; the linked install page gives 26.x today, which RiboPilot then pins to 25.10.6 by downloading that engine into NXF_HOME (this fails on a compute node with no internet). Nothing says Docker Desktop must be running, that Linux users need the docker group, that `nextflow` must be on the PATH that R sees, or that the first run pulls a 2.6 GB image (7.7 GB on disk).
- Why it matters: A new user cannot tell what to install, and the first run surprises them with a long download or a failed version pin.
- Fix: Add a "Before you start" checklist per operating system: Docker or Singularity/Apptainer, Java 17+, Nextflow with `export NXF_VER=25.10.6` in the shell profile, and R >= 4.4, each with a link and a one-line check (`docker run hello-world`, `java -version`, `nextflow -v`). Say the image is pulled once, on the first run.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-02"></a>**DOC-02. Supported operating systems are never stated.** Major (lowered from blocker by the verifier). Confirmed. ids: D-novice-5. Seen by 1 (docs-novice).
- Where: README Installation; `vignettes/RiboPilot.Rmd` Step 0.
- Problem: RiboPilot calls `bash` and uses `$(id -u)` in the Docker config, and Nextflow does not run natively on Windows, so Windows works only through WSL2 with R inside WSL. The image is amd64 only, so Apple Silicon Macs run it under emulation.
- Why it matters: A Windows user installs R for Windows and gets stuck with no explanation.
- Fix: Add "Supported computers": Linux and macOS (with a note on emulation on Apple Silicon), Windows only through WSL2 with R inside WSL, and clusters through Singularity/Apptainer.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-03"></a>**DOC-03. No hardware, disk, or run-time guidance, and the default resource requests stop small machines.** Major (D-novice-4 lowered from blocker by the verifier). Confirmed. ids: D-novice-4, D-hpc-12, D-intermediate-10, D-editor-6 (part). Seen by 4 (docs-editor, docs-hpc, docs-intermediate, docs-novice).
- Where: README Installation; Get started Step 0 and Step 4; `R/init_db.R:49-50, 192-193`; `inst/nextflow/nextflow.config` (retry rule).
- Problem: Every project, including the practice project, asks for 6 CPUs and 24 GB per assembly and 4 CPUs and 16 GB for read cleaning. A task that runs out of memory is retried with more (48, then 72 GB) up to 3 times and then quietly skipped. On the local executor that the docker profile uses, a task larger than the machine is refused ("Process requirement exceeds available memory"), and that ends the whole run, even under errorStrategy 'ignore' (the verifier re-ran this). No doc gives CPU, RAM, disk, or time figures.
- Why it matters: On a typical 4-core, 16 GB laptop the practice project cannot run as shipped, and nobody can plan an 80-sample project.
- Fix: Add a "Planning a project" table: minimum and recommended machine (for example 8 cores, 32 GB RAM, 20 GB free disk), per-sample needs, and measured times from this review (practice project about 16 min for 12 samples; order-level reference builds about 1.5 min; a Gastropoda-wide build about 8 min for 51,311 records). Say how to lower CPUs and Memory in the Assemble step's Assembly Options. Consider small defaults for the practice project, or `process.resourceLimits` in `config.docker`.
- Status: fixed on `review/2026-09-full-review`. The code half (PIPE-H2) caps oversized requests at the machine's size, and the Installation page now gives the measured run time, disk use, and the resource cap. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-H2

<a id="doc-04"></a>**DOC-04. No "Running on a cluster" page: opening the app on a login node, the head-job script, and shared storage.** Major. Confirmed. ids: D-hpc-6, D-intermediate-9, D-hpc-19, D-hpc-21. Seen by 2 (docs-hpc, docs-intermediate).
- Where: `vignettes/RiboPilot.Rmd:235-243, 290-292` ("see the function help ... or ask your computing support staff"); `README.Rmd:60-67, 109-115`; `man/generate_config.Rd`; `man/RiboPilot.Rd`; `inst/config.slurm:9-10` (comment only).
- Problem:
  - There is no worked `generate_config()` example, and nothing on reaching the app through an SSH tunnel (`RiboPilot(project, port = <port>, launch.browser = FALSE)` plus `ssh -N -L <port>:localhost:<port>`). The tunnel hint prints only in the console, and only when `port` is given.
  - The head-job script asks for 2 CPUs, 8 GB, and 48 h with no partition or account (even if the profile has them), keeps a literal `%%j` in the log name, and uses SGE `#$` lines for PBS without `-V`. The modal says only "Edit the job script if needed", with nothing on what to edit, where logs go (jobs/, logs/, nextflow_logs/), or that there is no live progress.
  - The project, reads, work dir, and the R library must be visible from head and compute nodes (the job runs `main.nf` from the installed package). This is only in template comments. `RIBOPILOT_SQLITE_JOURNAL` is undocumented.
  - The app itself runs the container on the node that serves it: the export ITSx check and the BLAST Options panel summary call `singularity exec` or `docker run` from the R session. If the image is not in the cache the jobs use, they fall back to `docker://<image>` and pull the whole image into `~/.singularity/cache` on the login node (4.9 GB in one pull). The panel summary runs under a 60 s timeout, so that pull is killed and retried each time the window opens.
  - The app has no login and listens on 127.0.0.1, which every user on a shared login node can reach.
- Why it matters: For a Slurm or SGE user this is where the docs run out. A research computing engineer had to read the code to get a working cluster run.
- Fix: New article "Running on a cluster": (1) build a profile with `generate_config()` (worked example), or use `NMNH_Hydra` / `NOAA_SEDNA` as-is; (2) keep the project, reads, and R library on shared storage; (3) open the app with `port` and `launch.browser = FALSE` and tunnel to it; (4) press Update, check partition, account, and walltime in the job script, then Submit as job; (5) where logs go; (6) the app needs singularity or docker on its own node and should share the jobs' image cache; (7) warn that Run from App stops when the browser closes, and that the port is open to other users on the node. Build the job header from the profile's queue and cluster options.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-79](USABILITY_REVIEW.md#ux-79); [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M14, PIPE-L2, APP-L5

<a id="doc-05"></a>**DOC-05. Nothing says where caches and big files go, how big they get, or how to move them.** Major. Confirmed. ids: D-hpc-4, D-novice-23, D-intermediate-27. Seen by 3 (docs-hpc, docs-intermediate, docs-novice).
- Where: README Installation; Get started Step 0; `vignettes/Reference-Databases.Rmd:114` ("RiboPilot uses a default location on your computer"); `man/generate_config.Rd`.
- Problem: These defaults fill a small home quota fast:
  - Singularity image in `<project>/work/singularity` (2.4 GB per project) unless `container_cache` or `NXF_SINGULARITY_CACHEDIR` is set. The built-in Hydra and SEDNA profiles set none.
  - Singularity layer cache `~/.singularity/cache` (grew 4.9 GB during one pull), plus temp space in `SINGULARITY_TMPDIR` (default /tmp).
  - `NXF_HOME` (`~/.nextflow`), which gets engine downloads when NXF_VER pins another version.
  - Reference folder `tools::R_user_dir("RiboPilot", "data")/ref_dbs` (move with `RIBOPILOT_REF_DIR` or `ref_dir`); BLAST panel (`RIBOPILOT_DB_DIR`); saved profiles `~/.config/R/RiboPilot`.
  - Taxonomy about 0.5 GB per project; Docker image 7.66 GB unpacked.
  - Generated profiles and SEDNA also lack the `pullTimeout` that Hydra needed (60 min).
- Why it matters: On clusters with small home quotas the first run or reference build fails on a full disk, and users cannot budget space.
- Fix: Add a "Disk space and caches" table (path, default, env var or argument to move it, typical size). Recommend `container_cache` or `NXF_SINGULARITY_CACHEDIR` on project or scratch space, `SINGULARITY_CACHEDIR` and `SINGULARITY_TMPDIR` off home, and `RIBOPILOT_REF_DIR` in the lab's `.Renviron`. Name `RiboPilot::ribopilot_ref_dir()` in the Reference guide. Add `pullTimeout` to generated profiles.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-L2, REF-L15

<a id="doc-06"></a>**DOC-06. No guidance on which steps need internet, or on compute nodes without it.** Major. Confirmed. ids: D-hpc-5. Seen by 1 (docs-hpc).
- Where: README Installation; Get started Step 0; `docker/README.md:107-109`.
- Problem: The Nextflow head pulls the image, downloads the pinned Nextflow, and stages the default seeds URL. Compute-node tasks call NCBI in `blast_ref_fetch` (efetch), `score` (seed download), remote BLAST when the panel cannot name a hit, and the Reference-tab count and build jobs. The R session downloads the taxonomy from ftp.ncbi.nlm.nih.gov. A working offline recipe (a pre-pulled SIF passed by absolute path as `container`, NXF_VER set and run once online, local seeds) was tested but is not written down.
- Why it matters: On a firewalled cluster an admin has to read the code to know what to open.
- Fix: Add a "Network access" list (what, from which node, which hosts) and the tested offline recipe. Warn that Reference builds and GenBank fetches need NCBI from compute nodes.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-07"></a>**DOC-07. Nothing explains what an upgrade does to existing projects.** Major. Confirmed. ids: D-hpc-7. Seen by 1 (docs-hpc).
- Where: `docker/README.md:58-64`; `README.Rmd:62-65`.
- Problem: `.config` fixes the image tag when the project is made. After an upgrade, old projects run the new `main.nf` and `bin/` scripts from the installed package against the old image's R code, with no check or warning. Every old image (7.66 GB each) must also be kept.
- Why it matters: A lab that upgrades mid-project can get silent mismatches between pipeline scripts and the R code inside the image.
- Fix: Document the rule: finish open projects before upgrading, or change `container` in each project's `.config` to the new tag. Better, warn at run time when the `.config` tag does not match `packageVersion()`.
- Status: fixed on `review/2026-09-full-review` (fourth batch, docs: an Upgrading section). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M15

<a id="doc-08"></a>**DOC-08. No troubleshooting page; the only troubleshooting text is a conda-only BLAST box.** Major. Confirmed. ids: D-intermediate-12, D-novice-45, D-intermediate-33, O-edge-27 (part), D-editor-17 (part). Seen by 4 (docs-editor, docs-intermediate, docs-novice, op-edge).
- Where: `vignettes/RiboPilot.Rmd:318-328` and `401-407` (GNUTLS box, conda only).
- Problem: Not explained anywhere: the five pipeline states (On hold, Ready to run, In progress, Success, Failed) and the State button; that one failed sample does not stop a run; the "Resume previous run?" box; the Work Dirs button; where `.nextflow.log` and job logs live; how to re-queue a failed sample; what "No taxonomy-locked samples", "Needs a database", or "No rawDir set" mean. The GNUTLS box that closes Get started applies only to the unsupported conda route (the image ships BLAST >= 2.12).
- Why it matters: When something fails, users have no path forward except reading code or asking the author. Every failed run in this review also ended with a green "Done".
- Fix: Add a Troubleshooting article with a "When something fails" routine (Notes column, Work Dirs, logs, set State to Ready to run, Update with Resume on) and a "Common errors" table (message, meaning, fix). The op-edge sweep and the usability review's error table are ready-made sources. Move the GNUTLS box there, labelled "conda only".
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-15](USABILITY_REVIEW.md#ux-15), [UX-42](USABILITY_REVIEW.md#ux-42), [UX-54](USABILITY_REVIEW.md#ux-54), [UX-55](USABILITY_REVIEW.md#ux-55)

<a id="doc-09"></a>**DOC-09. The mapping-file rules are not written down.** Major. Confirmed. ids: D-intermediate-6, D-novice-22, O-novel-8, O-edge-26, O-edge-27 (part), D-editor-17 (part). Seen by 5 (docs-editor, docs-intermediate, docs-novice, op-edge, op-novel).
- Where: `vignettes/RiboPilot.Rmd:177-185` (Step 3); `man/new_project.Rd` (mapping_fn).
- Problem: Enforced or applied, but not documented:
  - IDs of 18 characters or fewer, using only letters, digits, `_`, `-`, and `:` (catalog numbers with spaces or dots fail).
  - Paired reads only. Headers are case-sensitive. Comma-separated only (no .xlsx).
  - R1, R2, and Assembly are paths relative to `data_path`, so a subfolder (the usual layout from fasterq-dump) works; Get started says "the two file names".
  - "Already have an assembly?" says to leave R1 and R2 blank, which works, but a file with no R1 and R2 columns at all is refused with a message about paired read files, and a lowercase `assembly` header gives "these have neither".
  - `read.csv` turns spaces and dashes in column names into dots, and user columns named taxid, sci_name, rank, lineage, match_type, export_group, and a few others are overwritten.
  - Leading zeros in IDs are lost; gzipped assemblies cannot be mixed with plain ones.
  - There is no example table. (The optional Assembly column itself is documented under "Already have an assembly?".)
- Why it matters: IDs end up in GenBank records, and most wrong guesses produce raw R errors.
- Fix: Add a "Mapping file" reference section with a 3-row example CSV (one row using a subfolder), the ID rules, path rules, reserved column names, and how names are rewritten. Point to the practice project's `mapping.csv`. In code, treat missing R1 and R2 columns as blank when an Assembly column exists.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-03](USABILITY_REVIEW.md#ux-03), [UX-06](USABILITY_REVIEW.md#ux-06), [UX-55](USABILITY_REVIEW.md#ux-55); [CODE_REVIEW.md](CODE_REVIEW.md) DB-L5, DB-L6, DB-L15

<a id="doc-10"></a>**DOC-10. Where outputs land, and what can be deleted, is not documented.** Major. Confirmed. ids: D-intermediate-11. Seen by 1 (docs-intermediate).
- Where: `vignettes/RiboPilot.Rmd:207-208, 349-362`.
- Problem: The docs show only the export tree. The project also holds `.sqlite` (all results), `.config`, `.taxonomy/`, `out/<ID>/` (fastp reports and assemblies), `work/` (large Nextflow intermediates needed for resume), and the job-script and log folders. Most are hidden dotfiles.
- Why it matters: Users do not know what to back up or what they can delete when disk runs low.
- Fix: Add a "Project folder" section: each item, what it is for, what to back up (`.sqlite`, `.config`, `mapping.csv`), and what is safe to delete (`work/` once results are locked, noting that resume needs it).
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-11"></a>**DOC-11. There is no way to add samples to an existing project, and the docs do not say so.** Major. Confirmed. ids: D-intermediate-13. Seen by 1 (docs-intermediate).
- Where: `man/new_project.Rd` (force); `vignettes/RiboPilot.Rmd` Step 3.
- Problem: MitoPilot has `add_samples()`; RiboPilot has nothing like it, and `force = TRUE` deletes all results.
- Why it matters: Studies sequenced in batches have no documented path, and may lose work with `force = TRUE`.
- Fix: State the limit and a recommended approach (one project per batch), or port `add_samples()`.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-12"></a>**DOC-12. No NEWS, no citation, and no tool-version list for methods sections.** Major. Confirmed. ids: D-intermediate-14, D-domain-21, D-editor-17. Seen by 3 (docs-domain, docs-editor, docs-intermediate).
- Where: repo root (no `NEWS.md`, no `inst/CITATION`, no `CITATION.cff`); README "Behind the scenes"; `vignettes/Curation-Flags.Rmd:137-139`; `docker/Dockerfile`.
- Problem: There is no changelog although Docker Hub already has 16 dev tags (0.0.0.9000 to 0.0.0.9018). There is no "How to cite" and no list of the tools to cite: GetOrganelle, SPAdes, ITSx, barrnap/pybarrnap, BLAST, fastp, bowtie2, Nextflow, and Ribovore/ribotyper (run by the chimera check but never named). The Dockerfile leaves pybarrnap, barrnap, and itsx unpinned, so tool versions are not tied to an image tag. The intragenomic-variation claim ("real, documented biology in Octocorallia, Scyphozoa and other groups") has no references.
- Why it matters: Journal reviewers will ask how to cite RiboPilot and which tool versions produced the records.
- Fix: Add `NEWS.md`, `inst/CITATION` (and `CITATION.cff`, as MitoPilot has), a "Citing RiboPilot" section with tool citations, and a per-tag tool-version table. Pin every tool in the Dockerfile. Cite the paralog literature.
- Status: fixed on `review/2026-09-full-review` (fourth batch, docs: NEWS, CITATION, and a Citing page with a tool-version table). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M12

<a id="doc-13"></a>**DOC-13. No screenshots, and many app features are never described.** Major. Confirmed. ids: D-novice-18, D-appclaims-19, D-editor-17 (part). Seen by 3 (docs-appclaims, docs-editor, docs-novice).
- Where: all vignettes and README (the only images are the logo and a badge).
- Problem: The docs promise a point-and-click app but describe it in words only. Never covered:
  - the Taxonomy step (Assign Taxonomy, Edit Taxonomy, Confidence badges, Lock);
  - the State button and the five states; the Lock, State, and Columns pickers;
  - the Assemble options columns, including the Assembler switch and the Clade ruleset, and that Save re-queues;
  - the Triage window; the Assemble details window columns and panels; the Hit Source column;
  - Reviewed, Problematic, Partial, and Lock & Close;
  - Assign Group, Clear Group, the Mixed assemblers split, the Export option ticks, and Open folder;
  - Work Dirs, the CSV downloads, and the headless Submit as job / Save script route.
- Why it matters: Users learn the app by trial and error, and several of these controls change what gets exported.
- Fix: Add an "App tour" (or fold it into the rewritten walkthrough) with one screenshot per step and a short table per step: toolbar buttons, key columns, and the details window. The "Screenshots for the docs" section below has 15 ready candidates.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-51](USABILITY_REVIEW.md#ux-51), [UX-53](USABILITY_REVIEW.md#ux-53)

<a id="doc-14"></a>**DOC-14. No glossary; sample, path, scaffold, unit, and record are never tied together.** Minor. Confirmed. ids: D-novice-36, D-intermediate-23. Seen by 2 (docs-intermediate, docs-novice).
- Where: all pages.
- Problem: Undefined terms include genome skim, ETS, IGS, SSU, LSU, taxid, container, executor, profile, headless, R1/R2, scaffold, path, unit, seeds and labels database, core_nt, stock ITSx, Portal, defline, source modifiers, feature table, isomer, paralog, and rung. Annotation Review and Export rows are units (ID, path, scaffold); one sample can export several records named `ID_p<path>_s<scaffold>`, each with an automatic `[isolate=]`. Only tooltips explain this.
- Why it matters: Novices stall on vocabulary, and the several-records-per-sample idea surprises users at export.
- Fix: Add a short Glossary article, link terms on first use, and add a "Samples, scaffolds, and records" paragraph to Get started.
- Status: partly fixed on `review/2026-09-full-review` (docs rewrite): Get started ties sample, path, scaffold, unit, and record together. A glossary page is still open.

## Proposed outline for the published site

This outline follows MitoPilot's published site (Setup, walkthrough, Usage, Computing environments, Help, grouped reference), which RiboPilot users coming from MitoPilot will recognize. It puts the app first and gathers console commands and conda material in clearly labelled advanced pages.

**Home page (README.Rmd, about one screen):**

1. One-paragraph pitch and the fixed repeat diagram with its key ([DOC-23](#doc-23)).
2. "How it works": the app's five steps in order (Taxonomy, Reference, Assemble, Annotation Review, Export), with the toolbar screenshot (docs_01).
3. Which organisms: ray-finned fishes built in; other groups through the Reference step; the five rulesets.
4. Install in four lines, linking to Installation.
5. Try the practice project in two lines, linking to Get started.
6. Getting help (docs site and issues), citing RiboPilot, and license. Move the Hydra section to the cluster page.

**Articles:**

- **Setup**
  - *Installation.* Before you start (supported computers, hardware, disk); Docker or Singularity/Apptainer; Java 17+; Nextflow with `NXF_VER`; R and the package; a one-command check that everything works; a short private-preview note to remove at release. ([DOC-01](#doc-01), [DOC-02](#doc-02), [DOC-03](#doc-03), [DOC-17](#doc-17), [DOC-19](#doc-19))
- **Walkthrough**
  - *Get started (rewritten).* The practice project, step by step in the app: create it, open the app, then 1 Taxonomy, 2 Reference, 3 Assemble (Update, Run from App, reading the table and the details window, rescuing a failed sample), 4 Annotation Review, 5 Export. Add a "What you should see" box for CHIM00, CHIM01, and SRR22396893, the "keep the tab open" warning, and one screenshot per step. MitoPilot splits this into Assemble, Annotate, and Export pages, which would also work here if the page grows long. ([DOC-28](#doc-28), [DOC-35](#doc-35), [DOC-36](#doc-36), [DOC-37](#doc-37), [DOC-38](#doc-38), [DOC-13](#doc-13))
- **Usage**
  - *Your own project.* The mapping-file rules with an example CSV, `data_path`, the optional NCBI API key, choosing an executor, "Already have an assembly?", the project folder layout, and the limit on adding samples later. ([DOC-09](#doc-09), [DOC-10](#doc-10), [DOC-11](#doc-11), [DOC-31](#doc-31), [DOC-29](#doc-29))
  - *Reference databases* (the current page, updated). ([DOC-50](#doc-50), [DOC-51](#doc-51), [DOC-52](#doc-52), [DOC-69](#doc-69))
  - *Reviewing results* (the current Curation guide, updated): a warnings table with a "what to do" column, what the rDNA screen searches (local panel first, then NCBI), the five rulesets, the eight triage rules and the Triage window, chimeras and submission risk, variants (isomer versus ITS variant), and depth. ([DOC-53](#doc-53) to [DOC-63](#doc-63))
  - *Preparing a GenBank submission* (new, built from Step 9 and the docs-domain findings): what to upload (Portal FASTA versus table2asn), required source modifiers and the template tokens, which organism name is used, partial versus complete, flanks and trimming, the ITSx check, the submission_risk column, and what Reviewed, Problematic, and Partial do. ([DOC-40](#doc-40) to [DOC-49](#doc-49))
  - *Scripting RiboPilot (advanced).* `run_assembly()`, `run_fallback()`, `build_reference_db()`, `export_project()`, and `export_bundle()`, with a note that the app is the supported route.
- **Computing environments**
  - *Running on a cluster.* Profiles (`generate_config()`, `NMNH_Hydra`, `NOAA_SEDNA`), shared storage, opening the app (RStudio Server or an SSH tunnel), the headless rule, reviewing the job script, and logs. ([DOC-04](#doc-04), [DOC-69](#doc-69), [DOC-21](#doc-21), [DOC-73](#doc-73))
  - *Disk space, caches, and offline use.* The cache table, network needs, and the tested offline recipe. ([DOC-05](#doc-05), [DOC-06](#doc-06))
  - *Upgrading RiboPilot.* ([DOC-07](#doc-07))
  - *Installing tools without a container (unsupported).* All conda material in one place, including the GNUTLS tip. ([DOC-16](#doc-16))
- **Help**
  - *Troubleshooting.* "When something fails" and a common-errors table. ([DOC-08](#doc-08))
  - *Glossary.* ([DOC-14](#doc-14))
  - *FAQ* (optional): for example "Why does my flounder say suspect?" or "Why is my snail on the fish reference?".

**News and citation:** `NEWS.md` (shown under a News menu) and `inst/CITATION` plus `CITATION.cff` ([DOC-12](#doc-12)).

**Function reference groups** (from the docs-api audit; tested on a copy with `check_pkgdown()` reporting "No problems found"):

- Start a project: `new_test_project`, `new_project`, `RiboPilot` (fold `run_app` into this page).
- Reference databases: `build_reference_db`, `ribopilot_ref_dir`.
- Clusters: `generate_config`, `hydra_setup`.
- Export: `export_project`, `export_bundle`.
- Advanced scripting: `run_assembly`, `run_fallback`.
- Pipeline internals, marked `@keywords internal` so they leave the index (exported only because the Nextflow scripts call them): `annotate_records`, `auto_triage_sample`, `chimera_unit`, `compute_ref_alignment`, `curate_one`, `reconcile_one`, `itsx_check_one`, `fallback_assemble_one`, `load_assemblies`, `score_one`, `screen_sample`, `select_rrna_contigs`, `derive_blast_products`, `parse_blast_candidates`, `rederive_blast_products`, `mark_unit_searched`, `prime_from_ref_store`, `push_to_ref_store`, `stop_if_nothing_stored`, `write_blast_ref`, `write_blast_refs`, `write_coverage_units`, `write_rrna_candidates`, `open_project_db`.
- Could be unexported (no caller in `inst/` or `docker/`): `backfill_coverage_median`, `bitscore_weighted_lca`, `curate_rrna_core`, `new_db`, `rank_at_least`, `rank_distance`, `rank_distance_label`, `resolve_taxonomy`, `score_assembly`, `score_candidates`, `screen_one`, `set_unit_triage`, `shared_taxon`, `tk_ranked_lineage`, `triage_preview`, `write_coverage`, `write_path_scores`.

A `_pkgdown.yml` sketch for this outline (article file names are suggestions):

```yaml
navbar:
  structure:
    left: [intro, articles, reference, news]
    right: [search, github]

articles:
- title: Setup
  navbar: Setup
  contents: [Installation]
- title: Walkthrough
  navbar: ~
  contents: [RiboPilot]
- title: Usage
  navbar: Usage
  contents: [Your-Own-Project, Reference-Databases, Curation-Flags, GenBank-Submission, Scripting]
- title: Computing environments
  navbar: Computing environments
  contents: [Clusters, Caches-and-Offline, Upgrading, Without-Containers]
- title: Help
  navbar: Help
  contents: [Troubleshooting, Glossary]

reference:
- title: Start a project
  contents: [new_test_project, new_project, RiboPilot]
- title: Reference databases
  contents: [build_reference_db, ribopilot_ref_dir]
- title: Clusters
  contents: [generate_config, hydra_setup]
- title: Export
  contents: [export_project, export_bundle]
- title: Advanced scripting
  desc: The app's Update and Export buttons do the same work; these are for scripted runs.
  contents: [run_assembly, run_fallback]
```

## Page-by-page issues

Items are grouped by the page they mostly affect; the Where line gives every place a fix must touch. Each page ends with other items that also touch it.

### All pages

<a id="doc-15"></a>**DOC-15. The docs and the app use different names for the same things.** Major. Confirmed. ids: D-editor-10, D-novice-17, D-appclaims-24, O-practice-16. Seen by 4 (docs-appclaims, docs-editor, docs-novice, op-practice).
- Where: all vignettes and README; `R/app_ui.R:40-44`; `R/app_server.R:85, 366-368`; `R/app_modals.R:107, 396`; `R/constants.R:281`.
- Problem:
  - "Annotate tab" (10 times in the docs) vs "Annotation Review" (the app; 0 times in the docs).
  - "tab" (docs) vs "Step:" (toolbar) vs "Taxonomy mode" and "Assemble mode" (empty-table messages).
  - "Curation flags" vs "Curation warnings" and "Warnings".
  - "reconstruction" (Get started) vs path, scaffold, and unit (app and Curation guide).
  - "keep, drop or flag" vs keep / drop / undecided; "Press Run" and "Stop / Close" vs "Run from App" and "Stop / Interrupt"; "set their State to Ready" vs "Ready to run".
  - The Export empty message says "samples" for units; "Neighbouring" is British spelling; the README's "repeat unit" means the rDNA repeat while the app's "unit" means one scaffold.
- Why it matters: Users cannot match instructions to the screen.
- Fix: Adopt the app's toolbar vocabulary (Step names, Warnings, Run from App) in the docs and in app strings, and define the rest in the glossary.
- Status: fixed on `review/2026-09-full-review` (fourth batch, docs: the docs and the app use the same names). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).
- See also: [UX-77](USABILITY_REVIEW.md#ux-77), [UX-50](USABILITY_REVIEW.md#ux-50)

<a id="doc-16"></a>**DOC-16. The unsupported conda route is the main install recipe, and the code still defaults to it.** Major. Confirmed. ids: D-novice-9, D-hpc-10, D-novice-34, D-editor-20, D-hpc-25. Seen by 3 (docs-editor, docs-hpc, docs-novice).
- Where: `README.Rmd:69-83, 104`; `vignettes/RiboPilot.Rmd:65-112, 131, 151-162, 238, 255, 280, 325, 401-407`; `vignettes/Reference-Databases.Rmd:118, 140-160`; `R/init_project.R:52`; `R/init_test_project.R:164`; `R/nextflow_version.R:159-161`; `inst/config.conda:1-4`.
- Problem:
  - The only runnable setup block is `conda env create -f environment.yml`, but `environment.yml` is not in the installed package and lacks R, jq, and sqlite3.
  - Every pipeline command carries a conda comment, and tip boxes explain `conda activate`.
  - `new_project()` and `new_test_project()` still default to `executor = "conda"` with no warning, so the docs explain the workaround four times (README:104; RiboPilot.Rmd:131 and 238; Reference-Databases.Rmd:118).
  - The Nextflow hint says `mamba install -n ribopilot ...`, the GNUTLS alert says `conda install blast` without saying it is conda-only, and `config.conda` still says "HPC/container configs come later".
  - README and Get started repeat the same diagram, install steps, and practice-project text.
  - The README's conda link now redirects to an Anaconda concepts page, not install steps.
- Why it matters: Novices follow the only concrete commands they see, which lead to the unsupported path, and a user who keeps the default gets a run where every sample fails (fastp not found).
- Fix: Default `executor` to "docker" (or make it required) and warn once on conda; then delete the four notes. Change the Nextflow hint to `export NXF_VER=25.10.6`. Move all conda material into one "Installing tools without a container (unsupported)" article. Keep the README short and link to Get started for detail.
- Status: fixed on `review/2026-09-full-review`. The default is docker, conda warns, and all conda material is now on one page marked unsupported (Installing tools without a container), with the known gaps. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- See also: [UX-04](USABILITY_REVIEW.md#ux-04); [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M11

### README (home page)

<a id="doc-17"></a>**DOC-17. The install command fails on a fresh R: the Bioconductor dependencies are not found.** Blocker. Confirmed. ids: D-novice-2, D-intermediate-3, D-hpc-2, D-editor-3, D-api-2. Seen by 5 (docs-api, docs-editor, docs-hpc, docs-intermediate, docs-novice).
- Where: `README.Rmd:87-92`; `vignettes/RiboPilot.Rmd:101-103`; `DESCRIPTION` Imports.
- Problem: Biostrings, DECIPHER, and pwalign are Bioconductor packages in Imports. `remotes` adds the Bioconductor repositories only when DESCRIPTION has a `biocViews` field, and RiboPilot has none (MitoPilot does), so `remotes::install_github("dmacguigan/RiboPilot")` fails on any machine that lacks them. Get started also leaves out `install.packages("remotes")`.
- Why it matters: The very first command fails for new users.
- Fix: Add an empty `biocViews:` line to DESCRIPTION (tested on a copy; the three packages then resolve). Optionally show `BiocManager::install(c("Biostrings", "DECIPHER", "pwalign"))` or `pak::pak("dmacguigan/RiboPilot")`. Add the remotes line to Get started.
- Status: the code fix is on `review/2026-09-full-review` (empty `biocViews:` line); the install instructions need no change. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-H1

<a id="doc-18"></a>**DOC-18. The documented new_test_project() and new_project() calls crash outside RStudio.** Blocker (D-intermediate-2 rated major because RStudio users are not affected). Confirmed. ids: D-novice-3, D-intermediate-2, D-editor-2, D-api-1. Seen by 4 (docs-api, docs-editor, docs-intermediate, docs-novice).
- Where: `README.Rmd:101`; `vignettes/RiboPilot.Rmd:128, 189-197`; `vignettes/Reference-Databases.Rmd:105-111, 166-173`; `man/new_test_project.Rd` and `man/new_project.Rd` (Rproj); `R/init_project.R:98`.
- Problem: `Sys.getenv("RSTUDIO", FALSE)` returns the string "FALSE", so `!isFALSE(...)` is always TRUE and `rstudioapi::initializeProject()` runs. In a terminal, Rscript, R.app, RGui, VS Code, or Positron (with rstudioapi installed) the call stops with "Error: RStudio not running". `new_test_project()` does this after downloading every read, leaving only `mapping.csv`. Inside RStudio it opens a new session without warning. The help says Rproj has "no effect outside RStudio".
- Why it matters: Terminal and HPC users, the main audience for genome skims, fail at step one, and nothing mentions `Rproj = FALSE`.
- Fix: Test `nzchar(Sys.getenv("RSTUDIO"))` or `rstudioapi::isAvailable()`, and fix both help pages. Until the fix ships, add `Rproj = FALSE` to every documented call, and say that `Rproj = TRUE` opens a new RStudio session.
- Status: fixed on `review/2026-09-full-review`: no crash outside RStudio, and the `new_project()` help now matches `new_test_project()` ("no effect outside RStudio"). See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- See also: [UX-01](USABILITY_REVIEW.md#ux-01); [CODE_REVIEW.md](CODE_REVIEW.md) PKG-H2

<a id="doc-19"></a>**DOC-19. While the repo is private, three documented steps need a GitHub token, and the default fish reference fails even with one.** Blocker (D-novice-19 and D-editor-7 rated major, D-api-22 minor; the item takes D-hpc-3's blocker because no token fixes the seeds). Confirmed. ids: D-hpc-3, D-novice-19, D-editor-7, D-api-22, D-intermediate-31 (part), O-novel-7 (part). Seen by 6 (docs-api, docs-editor, docs-hpc, docs-intermediate, docs-novice, op-novel).
- Where: `README.Rmd:16, 24, 54, 87-107, 119-120`; `vignettes/RiboPilot.Rmd:102, 121-144`; `docker/README.md:107-109`; `R/constants.R:24-29`; `R/score_project.R:11-17`; `DESCRIPTION:12-13`; `R/app_server.R:85`.
- Problem:
  - `install_github()` and the practice reads need `GITHUB_PAT`, and no doc says so. Without it, the practice-project error leads with "A slow or interrupted connection is the usual cause" before the token hint.
  - The default fish seeds and labels are raw.githubusercontent.com URLs that return 404. Nextflow stages them without auth and the score task uses `download.file()`, so no token helps. The practice project's first run aborts when Nextflow cannot stage the seeds (op-practice), or scoring fails for every sample (op-edge), and the app still says Done. Every new project, fish or not, also warns that these databases are unreachable.
  - The site URL, every Source link, the GitHub navbar icon, BugReports, the issue links, and the in-app Help link all 404 as well.
- Why it matters: Every pre-release tester stops here, and the first practice run looks as if it succeeded.
- Fix: Treat "repo public and Pages deployed" as a release gate. Ship the fish seeds and labels as checksummed release assets, or inside the package or image, not as a main-branch raw URL. Until then, add a short "Private preview" note (set `GITHUB_PAT`; pass local seeds), and report 404s as "not found or private" rather than a slow connection.
- Status: the fish-reference half is deferred by the author (2026-09-24): most users will build their own reference databases, so shipping the built-in fish reference is not a release gate. The Installation page now has a private-preview note on `GITHUB_PAT` and the unreachable fish reference, to delete at release. The repo-public half (site URL, Source links, BugReports, Help link) still applies.
- See also: [UX-11](USABILITY_REVIEW.md#ux-11); [CODE_REVIEW.md](CODE_REVIEW.md) SCI-M9, REF-L4, REF-L9

<a id="doc-20"></a>**DOC-20. The README says the taxonomy lock gates assembly; the pipeline does not check it.** Major. Confirmed. ids: D-intermediate-4, D-editor-5, D-appclaims-2, D-api-5 (part). Seen by 4 (docs-api, docs-appclaims, docs-editor, docs-intermediate).
- Where: `README.Rmd:44` ("samples must be reviewed and locked here before they move on to assembly"); `inst/nextflow/modules/pipeline_r.nf:40`; `R/nf_command.R:167`; `R/init_db.R:243`.
- Problem: Every sample starts queued (`assemble_switch = 1`). The pipeline's sample query and the Update count select on `assemble_switch` and `assemble_lock` only, so `run_assembly()` or Update assembles samples whose taxonomy is unlocked or unresolved; only the Assemble table hides them. In a fresh project the Assemble table shows 0 rows while Update offers "2 sample(s) to update".
- Why it matters: Samples assemble with no lineage (fish seeds, generic ruleset, contaminant_unchecked), and the user cannot see them in the table.
- Fix: Add `s.taxonomy_lock = 1` to the sample query and `nf_pending_count()` (preferred: it makes the README true), or reword the README to say the Assemble table lists only taxonomy-locked samples.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-09](USABILITY_REVIEW.md#ux-09)

<a id="doc-21"></a>**DOC-21. The Hydra section is out of date: the app runs hydra_setup() itself.** Minor. Confirmed. ids: D-novice-32, D-editor-25, D-hpc-17, D-api-25. Seen by 4 (docs-api, docs-editor, docs-hpc, docs-novice).
- Where: `README.Rmd:109-115`; `man/hydra_setup.Rd:16`; `R/app_run.R:184`; `R/hydra.R:71-76, 101-107`.
- Problem: The README and help say to call `hydra_setup()` "once per session, before launching the app", but `run_app()` calls `ensure_hydra_setup()` when it detects Hydra without qsub. The head-job script runs `source ~/.bashrc`, `module load java`, then bare `nextflow`, so Nextflow must be on PATH through `~/.bashrc` or `~/bin`; the README does not say this. Institution-specific content also sits on the public landing page.
- Why it matters: Hydra users do a step they do not need and miss the one that matters.
- Fix: Say the app applies the setup itself and `hydra_setup()` is only for console runs; say where to install Nextflow on Hydra; move the section to the cluster article.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-22"></a>**DOC-22. README.md has a paragraph that README.Rmd lacks.** Minor. Confirmed. ids: D-hpc-13, D-novice-33, D-editor-19, D-api-27 (part). Seen by 4 (docs-api, docs-editor, docs-hpc, docs-novice).
- Where: `README.md:173-175` vs `README.Rmd:96-107`.
- Problem: The "By default this fetches small rDNA-filtered read subsets ... so setup takes seconds" paragraph (with `full_size`) exists only in the generated README.md, which pkgdown uses for the home page.
- Why it matters: The next knit of README.Rmd silently drops it.
- Fix: Move it into README.Rmd and re-knit. Mention the one-time 2.6 GB image pull there, so "seconds" is not misleading.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L8

<a id="doc-23"></a>**DOC-23. The repeat diagram puts LSU under 5.8S, draws ETS as a gene, and has no key.** Minor. Confirmed. ids: D-domain-17, D-editor-23, D-novice-39. Seen by 3 (docs-domain, docs-editor, docs-novice).
- Where: `README.Rmd:30-34`; `vignettes/RiboPilot.Rmd:39-43`.
- Problem: "LSU" sits at columns 24-26, under 5.8S, not under 28S (columns 34-44), and "SSU" straddles the ETS/18S boundary. ETS is drawn with `=` like a gene. IGS (which holds the NTS and the next unit's 5' ETS), the separate 5S array in animals, SSU, LSU, and "ITS" are never defined. The 81-character line scrolls sideways at both 1400 px and 390 px.
- Why it matters: It is the site's first picture of the biology, and a specialist reader will notice at once.
- Fix: Put LSU under 28S, draw ETS as a spacer, put "(repeated hundreds of times)" on its own line, and add a one-line key: SSU = 18S, LSU = 28S, ITS = ITS1-5.8S-ITS2, IGS = NTS + ETS, and 5S is a separate array in animals.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-24"></a>**DOC-24. The taxonomy snapshot's size and time are understated, and its download can time out.** Minor. Confirmed. ids: D-hpc-15, D-novice-29, D-editor-36. Seen by 3 (docs-editor, docs-hpc, docs-novice).
- Where: `README.Rmd:44` ("~72 MB ... about a minute"); `R/tax_index.R:91` ("a few minutes"); `R/taxonomy.R:84-95`.
- Problem: Each project stores a 76-80 MB tarball plus a 382-400 MB SQLite index, about 0.5 GB. The app says indexing takes "a few minutes"; the README says "about a minute". `ensure_taxdump()` keeps R's 60 s download timeout, so Assign Taxonomy can fail on links slower than about 1.3 MB/s (read from the code; other downloads raise the timeout).
- Why it matters: Users under-budget disk and cannot explain a timeout on a slow link.
- Fix: Say "about 80 MB download and about 0.5 GB per project; indexing takes a few minutes the first time". Raise the timeout the way `download_or_stop()` does.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) REF-L15

<a id="doc-25"></a>**DOC-25. The Nextflow version claim and hint point the wrong way.** Minor. Confirmed. ids: D-novice-31, D-intermediate-32. Seen by 2 (docs-intermediate, docs-novice).
- Where: `README.Rmd:77`; `R/zzz.R:26-35`; `R/nextflow_version.R` (nf_update_hint).
- Problem: "RiboPilot refuses to load on an older Nextflow" is true only for `library(RiboPilot)`; the documented `RiboPilot::` calls check at app start or run time instead. The update hint recommends `mamba install -n ribopilot ...`, the unsupported path, and may even mean going down a version.
- Why it matters: The message sends container users to a conda env they do not have.
- Fix: Reword to "RiboPilot stops with a message when it finds an unsupported Nextflow", and point the hint to the Nextflow install page plus `export NXF_VER=25.10.6`.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L1

<a id="doc-26"></a>**DOC-26. environment.yml says the conda and container paths agree; they do not.** Minor. Confirmed. ids: D-hpc-18. Seen by 1 (docs-hpc).
- Where: `environment.yml:1-8`; `README.Rmd:78`.
- Problem: The image runs barrnap 1.10.5, entrez-direct 26.0, and datasets 18.36.0. The env pins barrnap 0.9 (with a comment that 1.10.x is not what the GFF handling was written against), entrez-direct 24.0, and datasets 18.34.0, and has no R, jq, or sqlite3. The file is also not installed with the package.
- Why it matters: The only concrete install recipe promises a parity it does not have.
- Fix: Pin the Dockerfile to match, or change the header to say the env is for development only and differs from the image.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M12

<a id="doc-27"></a>**DOC-27. The knitted README.md is not ASCII.** Polish. Confirmed. ids: D-hpc-26, D-editor-42. Seen by 2 (docs-editor, docs-hpc).
- Where: `README.md:29, 51, 114, 134` (generated from `README.Rmd`).
- Problem: Pandoc's smart typography turns the straight quotes in README.Rmd into curly quotes (U+2019, U+201C, U+201D) in README.md, which pkgdown and GitHub show as the home page, for example in "your lab's server" and "genome skim".
- Why it matters: The landing page breaks the ASCII-only house style and drifts further from its source.
- Fix: Knit with `md_extensions: -smart` on the `github_document` output.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L8

Also on this page: [DOC-28](#doc-28) (the "What RiboPilot does" order), [DOC-54](#doc-54) (rulesets sentence), [DOC-16](#doc-16) (conda recipe), [DOC-29](#doc-29) ("Already have an assembly?" summary), [DOC-98](#doc-98) (link text).

### Get started (vignettes/RiboPilot.Rmd)

<a id="doc-28"></a>**DOC-28. Get started describes an older, console-driven app with three tabs.** Blocker (D-novice-1, D-intermediate-1, D-editor-1, and D-appclaims-1 blocker; the others major). Confirmed. ids: D-novice-1, D-intermediate-1, D-editor-1, D-appclaims-1, D-novice-10, D-intermediate-5, D-editor-4, D-api-5, O-practice-5, O-practice-6, O-novel-6, O-edge-27, D-hpc-6 (part). Seen by 9 (docs-api, docs-appclaims, docs-editor, docs-hpc, docs-intermediate, docs-novice, op-edge, op-novel, op-practice).
- Where: `vignettes/RiboPilot.Rmd:51-57` ("four stages"), `140-173` (Step 2), `245-310` (Steps 4-6), `294-296` ("three tabs ... Assemble, Annotate, and Export"), `318-328`, `385-394` ("at a glance"); `R/app_ui.R:40-44`; `R/app_server.R:174, 366`.
- Problem:
  - The app has five numbered steps (1 Taxonomy, 2 Reference, 3 Assemble, 4 Annotation Review, 5 Export) and opens on Taxonomy. Get started never mentions Taxonomy, Assign Taxonomy, or locking taxonomy, and gives the Reference tab only as a side note.
  - Steps 4, 5, and 9 run `run_assembly()`, `run_fallback()`, and `export_project()` from the console right after `new_project()`. The pipeline does not check the taxonomy lock, so samples assemble with no lineage: fish seeds, the generic ruleset, and contaminant_unchecked on every unit.
  - The app then opens on Taxonomy, and the Assemble table says "No taxonomy-locked samples. Assign and lock taxonomy in Taxonomy mode first."
  - The app's own run path (Update, then Run from App; the Assembler option for rescue; Assign Group and Export Data) is not described, and the console route skips the app's guards (the missing-reference-database block and the pending-group warning).
- Why it matters: Nine of the ten reviewers reported it, including all three operators. It is the main reason nobody could finish from the docs, and it goes against the goal of an app-first guide.
- Fix: Rewrite Steps 2 to 9 around the app: `new_project()`, open the app, then 1 Taxonomy (Assign Taxonomy, review Confidence, Edit Taxonomy, Lock); 2 Reference (Propose groups, Apply, Set up database, Run from App; ray-finned fishes can skip); 3 Assemble (Update, Run from App, review, Triage, Keep/Drop, Lock); 4 Annotation Review (warnings, Reviewed/Problematic/Partial, trim, Lock); 5 Export (Assign Group, Export Data, templates). One screenshot per step. Move `run_assembly()`, `run_fallback()`, `build_reference_db()`, and `export_project()` into an "Advanced: scripting" article.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-49](USABILITY_REVIEW.md#ux-49), [UX-09](USABILITY_REVIEW.md#ux-09)

<a id="doc-29"></a>**DOC-29. The documented "Already have an assembly?" route fails on every generate_config() Singularity profile.** Blocker. Confirmed. ids: D-hpc-1. Seen by 1 (docs-hpc).
- Where: `vignettes/RiboPilot.Rmd:210-233`; `README.Rmd:66-67`; `R/generate_config.R:38`; `R/init_project.R:177`; `R/nf_command.R:52`; `R/reference_nf.R:159-160`.
- Problem: The generated engine block always binds `${params.rawDir}`. A project with no reads gets `rawDir = "NA"`, and `nf_command()` leaves out `--rawDir`, so Nextflow fails at config parse ("Unknown config attribute singularity.params.rawDir"). The reference build passes `--rawDir NA`, and Singularity dies with "unable to add NA to mount list". The built-in docker, local, and Hydra templates are not affected.
- Why it matters: Cluster users with supplied assemblies cannot run at all.
- Fix: Drop the rawDir bind from `container_engine_block()` (`new_project()` already mounts the read folder through EXTRA_MOUNTS). Until then, tell users to delete `-B ${params.rawDir}` from `.config`.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-H3

<a id="doc-30"></a>**DOC-30. RiboPilot("~/my_project") fails as written, and the summary block uses an undefined variable.** Major. Confirmed. ids: D-novice-12, D-intermediate-20, D-editor-9, D-api-9, D-novice-41. Seen by 4 (docs-api, docs-editor, docs-intermediate, docs-novice).
- Where: `vignettes/RiboPilot.Rmd:291, 391`.
- Problem: Every other call uses `RiboPilot::` and no page calls `library(RiboPilot)`, so this line fails with 'could not find function "RiboPilot"' at the step that opens the app. The "at a glance" block also uses `project` without defining it.
- Why it matters: The single most important command in the guide fails.
- Fix: Use `RiboPilot::RiboPilot("~/my_project")`, and add `project <- "~/my_project"` at the top of the summary block.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-31"></a>**DOC-31. The NCBI API key: the placeholder is sent to NCBI if left in, and the environment variable is ignored.** Major (D-novice-6 lowered from blocker by the verifier; D-editor-22 and O-practice-21 minor). Confirmed. ids: D-novice-6, D-editor-22, O-practice-21. Seen by 3 (docs-editor, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:195` (`ncbi_api_key = "PASTE_YOUR_KEY_HERE"`); `R/init_project.R:181`; the practice project's `.config`.
- Problem: The key is optional, but the example includes a fake value and never says to delete the line. `new_project()` writes it verbatim into `.config`, and it is exported as NCBI_API_KEY. NCBI answers HTTP 400 "API key invalid" to every request that carries it (reference fetch, Reference-tab counts, and builds). The guide also says you paste the key when creating a project, but the Step 1 command has no key argument, and `new_project()` does not read `NCBI_API_KEY` from the environment (the reference build and the app do), so the practice project's `.config` gets an empty key.
- Why it matters: A user without a key silently breaks every NCBI lookup.
- Fix: Comment the line out with "optional: delete this line if you have no key", and check the key format in `new_project()`. Default `ncbi_api_key` to `Sys.getenv("NCBI_API_KEY")`, as the reference code already does.
- Status: docs half fixed on `review/2026-09-full-review` (docs rewrite): no placeholder key in any page, and Installation recommends `NCBI_API_KEY` in `~/.Renviron`. The code default is still open.
- See also: [UX-07](USABILITY_REVIEW.md#ux-07)

<a id="doc-32"></a>**DOC-32. Step 2 sends non-fish users to an unsupported console build before the project exists.** Major (D-novice-7 lowered from blocker by the verifier; O-practice-20 polish). Confirmed. ids: D-novice-7, D-intermediate-16, O-practice-20, D-editor-21 (part). Seen by 4 (docs-editor, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:140-173, 199-205`; `vignettes/Reference-Databases.Rmd:136-173`.
- Problem: The console build uses `conda_env = "ribopilot"` and needs Entrez Direct and biopython on the host; on the supported Docker setup it stops with "Could not locate conda env". It is shown first, and the supported route only as a trailing note (Step 0 does name Step 2 as the one place `conda_env` is used). The note points to the Reference tab, but that needs a project (Step 3) with assigned taxonomy, and Step 3 then uses `db$seeds`, which does not exist on the app route. The gastropod example matches 51,311 records (about twice the app's 25,000 ceiling) and took about 8 minutes.
- Why it matters: Non-fish users either fail on Docker or take an unsupported route with an oversized database.
- Fix: For container users: create the project without reference options, then do Taxonomy and Reference in the app. Move the console build to the scripting article, with a narrower example query and its size and time.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-33"></a>**DOC-33. "If your samples are fish" should say ray-finned fishes.** Polish. Confirmed. ids: D-domain-24. Seen by 1 (docs-domain).
- Where: `vignettes/RiboPilot.Rmd:143-144`.
- Problem: Sharks, rays, lampreys, hagfishes, and lungfishes are fish, but they get the generic ruleset (plus no_clade_ruleset) and non-matching seeds; the ruleset matches only Actinopteri/Actinopterygii.
- Why it matters: Users with other fishes skip a step they need.
- Fix: "If your samples are ray-finned fishes (Actinopterygii), you can skip to Step 3."

<a id="doc-34"></a>**DOC-34. Code fragments need a comma the reader is not told about, and rely on an in-memory object.** Minor. Confirmed. ids: D-novice-37, D-editor-21. Seen by 2 (docs-editor, docs-novice).
- Where: `vignettes/RiboPilot.Rmd:199-205, 217-219, 235-243`; `vignettes/Reference-Databases.Rmd:165-173`.
- Problem: "Add these two lines inside the parentheses" (`custom_seeds_db = db$seeds,` ...) needs a comma after the previous argument, and `db` is gone after an R restart. The executor alert sits after "Already have an assembly?", far from the call it explains.
- Why it matters: Copy and paste fails for novices.
- Fix: Show complete calls with file paths (or note that re-running `build_reference_db()` returns the existing paths), and move the alert next to the call.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-35"></a>**DOC-35. The practice project leaves the reader without a next step, and its samples are unexplained.** Major. Confirmed. ids: D-novice-13, D-appclaims-18, D-domain-22, O-practice-19, D-api-26, D-intermediate-19 (part), D-editor-21 (part). Seen by 7 (docs-api, docs-appclaims, docs-domain, docs-editor, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:121-138`; `README.Rmd:96-107`; `man/new_test_project.Rd:49` ("the bundled 10-sample test set"); `R/init_project.R:199-200`.
- Problem: After Step 1 the guide switches to `~/my_project` and Steps 2-3, which practice users should skip, and never says "now open the app with this path". The project has 12 samples: 10 real SRA subsets plus the synthetic CHIM00 (a copy of an SRR22396758 scaffold) and CHIM01 (that fish rDNA joined to Penicillium). Both are labelled "Upeneus parvus", nothing in the mapping file marks them synthetic, and they export with real-looking GenBank headers (">CHIM00 [organism=Upeneus parvus] ... complete sequence"). No doc says what to expect from them. The help says "the bundled 10-sample test set" while the console prints "Fetching 12 test sample(s)". SRR22396893 has very thin data (about 0.6x) and fails. The command ends with "Please open and review the .config file", a hidden file that needs no edits for Docker.
- Why it matters: The practice project is the first thing users run, and they cannot tell expected results from errors.
- Fix: Add a short practice walkthrough (create, open the app with `~/ribopilot_test`, the five steps) and a "What you should see" box: CHIM01 is caught as off_target with a junction near 2,625 bp, CHIM00 is clean, and SRR22396893 fails and is the rescue example. Say 12 samples (10 downloaded plus 2 synthetic). Mark CHIM00 and CHIM01 as synthetic in the mapping file and leave them out of export by default. Drop or reword the `.config` message.
- Status: docs half fixed on `review/2026-09-full-review` (docs rewrite): Get started explains the 12 samples and ends with next steps, and the new_test_project() help says 12 samples. The "Next: open the app" message at the end of new_project() is still open.
- See also: [UX-11](USABILITY_REVIEW.md#ux-11)

<a id="doc-36"></a>**DOC-36. Closing a run window or the browser tab stops the run, but the docs say you can keep working or come back later.** Major (D-novice-11 rated minor). Confirmed. ids: D-appclaims-6, D-novice-11. Seen by 2 (docs-appclaims, docs-novice).
- Where: `vignettes/Reference-Databases.Rmd:61-62` ("you can keep working in the app while it runs"); `vignettes/RiboPilot.Rmd:258-259` ("You can let it run and come back later"); `R/app_reference.R:980-996`; `R/app_server.R:1029-1032`; `R/app_lifecycle.R:23-30`.
- Problem: The build window is a modal that covers the app, and its Close asks "Stop the build?". Closing or reloading the tab kills both the reference build and an Update run, and the app stops 5 s after the last session closes. "Come back later" is true only for the console route.
- Why it matters: Users lose long runs by closing a laptop lid or a tab, and a dropped SSH tunnel kills a cluster run started from the app. One reviewer's first Singularity run died this way when the driving browser tab closed.
- Fix: Warn in the Assemble and Reference steps: keep the window and tab open and the computer awake; for long runs use the cluster Submit as job route; if it closes, press Update again with Resume on. Better, keep runs alive after the session ends.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-43](USABILITY_REVIEW.md#ux-43); [CODE_REVIEW.md](CODE_REVIEW.md) APP-M1, REF-L17

<a id="doc-37"></a>**DOC-37. Step 5 rescue: the app route is undocumented, and run_fallback() re-assembles good samples.** Major. Confirmed. ids: D-appclaims-8, O-practice-10, D-novice-10 (part), D-intermediate-16 (part). Seen by 4 (docs-appclaims, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:270-283`; `R/run_fallback.R:11-17`; `R/app_opts.R:95-105`.
- Problem: The guide shows only `run_fallback()`. In the app the rescue is: Assembly Options cell, type a new set name, choose "SPAdes, baited to the seed DB", Save, Update. `run_fallback()` picks any unlocked sample with a failed assembly OR any unit missing 18S or 28S, so a sample with one complete unit plus one fragment (SRR33849683 in the practice project) is re-assembled with the bait assembler. Failed rows say only "see Notes", with no hint to try the fallback.
- Why it matters: Users either miss the rescue or re-assemble good samples.
- Fix: Document the app route with a screenshot (docs_10 below). Restrict `fallback_candidates()` to samples with no complete kept unit. Have failed-sample notes say "Try the baited SPAdes assembler in Assembly Options".
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-15](USABILITY_REVIEW.md#ux-15), [UX-16](USABILITY_REVIEW.md#ux-16), [UX-17](USABILITY_REVIEW.md#ux-17); [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-H5

<a id="doc-38"></a>**DOC-38. The Assemble-step instructions describe controls and re-runs that do not exist.** Major. Confirmed. ids: D-appclaims-5. Seen by 1 (docs-appclaims).
- Where: `vignettes/RiboPilot.Rmd:303-305, 318-328`; `R/app_modals.R:121-170, 349`; `R/pipeline_steps.R:479-580`.
- Problem: The guide says to "switch to an alternative reconstruction or adjust the exact endpoints" in the Assemble tab, then press Update so the sample is annotated again. The Assemble details window has only Keep (K) / Drop (D), Copy FASTA, and evidence panels. Trimming is in Annotation Review ("Trim assembly..."), and it shifts coordinates and re-curates at once with no run. Keep/Drop needs no re-annotation, because every unit was annotated in the first run.
- Why it matters: Users look for controls that are not there and start runs they do not need.
- Fix: Say: keep or drop scaffolds in Assemble; trim or nudge boundaries in Annotation Review (applied at once); press Update only after changing an option set or setting State to Ready to run.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-39"></a>**DOC-39. "The spacers are simply the stretches between those landmarks" is not what ships: 1 bp and 4 bp are left unannotated.** Minor (lowered from major by the verifier: 1 to 4 bp at the ends of boundaries that two different tools predict). Confirmed. ids: D-domain-4. Seen by 1 (docs-domain).
- Where: `vignettes/RiboPilot.Rmd:47-48`; `R/annotate.R:3, 97, 106-190`.
- Problem: Reconcile takes the 18S and 28S ends from barrnap and ITS1/ITS2 from ITSx, and only clips overlaps; it never fills gaps. All 10 practice units that have both features leave 1 unannotated base between 18S and ITS1, and 9 of 10 leave 4 between ITS2 and 28S (the code comments even say "gap-free"). Those bases are missing from the ITS marker files and unannotated in the `.tbl`, and no flag reports gaps. Screenshot docs_12 shows one: 18S ends at 2041, ITS1 starts at 2043.
- Why it matters: Submitted ITS records are short by up to 5 bp, and table2asn warns on them.
- Fix: After clipping, extend each spacer to meet its flanking genes (or add a gap flag); then the sentence in the guide is true.
- See also: [UX-84](USABILITY_REVIEW.md#ux-84); [CODE_REVIEW.md](CODE_REVIEW.md) DB-M5

<a id="doc-40"></a>**DOC-40. export_project() does not "export everything".** Major. Confirmed. ids: D-novice-21, D-intermediate-8, D-api-6, D-appclaims-3, O-practice-14, O-practice-26. Seen by 5 (docs-api, docs-appclaims, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:342-347, 393`; `man/export_project.Rd` (title only); `R/export.R:736-752`.
- Problem: With no export groups it stops with "No export groups assigned. Assign samples to a group first.", and groups can only be made in the app. It always uses the "default" header template, skips Problematic units, writes a `<group>.zip` (the app export does not), and stops part-way if a group has no annotation-locked units, after earlier groups were already written.
- Why it matters: The guide's R export fails on every fresh project, and the error does not say where groups come from.
- Fix: Make the app's Export step (Assign Group, then Export Data) the documented path. Describe `export_project()` as "exports each group made in the Export step, with the default template". Skip empty groups, and point the error at the Export step.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-37](USABILITY_REVIEW.md#ux-37)

<a id="doc-41"></a>**DOC-41. Step 9 leaves out behaviors that change what is submitted; Problematic silently drops units.** Major (D-novice-28 and O-practice-17 rated minor). Confirmed. ids: D-intermediate-22, D-appclaims-4, D-editor-16, D-novice-28, O-practice-17. Seen by 5 (docs-appclaims, docs-editor, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:330-383`; `R/app_ui.R:93-94`; `R/app_server.R:1061-1068, 1931`; `R/export.R:545-552`.
- Problem:
  - Units marked Problematic are left out of every export with no message; the tip only says "Flag this annotation for another look". Verified: SRR14433904 (locked, grouped, Problematic) was absent from the bundle, and in op-practice CHIM01 vanished the same way (15 of 16 units exported).
  - Partial marks feature ends partial in the `.tbl` and changes the definition line. Reviewed, Lock & Close (which also marks Reviewed), and the one-path-per-sample lock rule are not described, so a user who only locks gets reviewed = no on every exported record.
  - Groups are split by assembler, and the fact that the Portal takes only the FASTA is only in the bundle README.
- Why it matters: A curator's flag removes a record from the submission without warning.
- Fix: Rewrite Step 9 around Assign Group, Export Data, and the templates. Explain Reviewed, Problematic, and Partial and their export effects. Change the tips to say "Problematic units are left out of export", and show the excluded count in the Export window.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-33](USABILITY_REVIEW.md#ux-33), [UX-34](USABILITY_REVIEW.md#ux-34)

<a id="doc-42"></a>**DOC-42. Partial or complete is decided by scaffold position, not gene completeness, and no doc explains it.** Major (lowered from blocker by the verifier, because the Portal route does its own annotation; the usability review rates the same bug blocker, as UX-32). Confirmed. ids: D-domain-1. Seen by 1 (docs-domain).
- Where: `vignettes/RiboPilot.Rmd` Step 9 (340-383); Partial toggle tooltip `R/app_server.R:1933`; `R/export_opts.R:135-141`; `R/export.R:144-153`.
- Problem: The first feature is marked partial only if it starts at base 1, and the last only if it ends at the last base; the stored model coverage (already used for the truncated flag) is ignored. SRR33849683_p1_s2, a 28S at 45.8 percent of the model, gets "28S ribosomal RNA gene, complete sequence." A trimmed SRR14163830 with 100 percent coverage 18S and 28S gets `<1` and `>7240` and "partial sequence" at both ends. The Partial toggle marks both ends even when only one is incomplete.
- Why it matters: Public records carry the wrong completeness, and users are not told to check it.
- Fix: In code, set `<` and `>` from model coverage at each end (or from the truncated flag), and let Partial choose 5', 3', or both. In Step 9 and the Curation guide, explain the rule, and until the fix ships tell users to check completeness before upload.
- Status: fixed on `review/2026-09-full-review`. The truncated flag sets both `<` and `>` and the "partial sequence" wording, and Preparing a GenBank submission now explains the partial rule. The 5' or 3' choice for Partial is still open. See [Fixed on this branch](REVIEW_SUMMARY.md#fixed-on-this-branch-after-the-review).
- See also: [UX-32](USABILITY_REVIEW.md#ux-32); [CODE_REVIEW.md](CODE_REVIEW.md) DB-H2, DB-L8

<a id="doc-43"></a>**DOC-43. Exports carry unannotated ETS and IGS flanks; trimming and Portal length limits are never mentioned.** Major. Confirmed. ids: D-domain-12. Seen by 1 (docs-domain).
- Where: `vignettes/RiboPilot.Rmd` Step 9; README "GenBank-ready submission package"; `R/export.R:623-627`; `R/constants.R:99-112`.
- Problem: Export writes the whole scaffold. Practice records carry 170-2,150 bp upstream of 18S plus downstream IGS, none of it annotated or named in the definition line. NCBI's Portal help says "Sequences longer than the expected length for a given submission type will be removed." The guides never mention the Trim tool, and trimming triggers the partial problem above. ETS and IGS product strings exist in the code but are never emitted.
- Why it matters: Records may be removed at the Portal, or go out with unexplained flanking sequence.
- Fix: In Step 9, explain the flanks and recommend trimming to the 18S-28S span for Portal rRNA-ITS submissions (after the partial fix). Optionally annotate partial ETS and IGS features for table2asn.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-44"></a>**DOC-44. Required source modifiers are not stated, and isolate is banned in templates but added automatically.** Major (D-appclaims-13 rated minor). Confirmed. ids: D-domain-13, D-appclaims-13, D-intermediate-7. Seen by 3 (docs-appclaims, docs-domain, docs-intermediate).
- Where: `vignettes/RiboPilot.Rmd:378-383` ("The template is the whole defline"); `R/app_modals.R:622-625`; `R/export_opts.R:18, 219-253`.
- Problem:
  - For eukaryotic rRNA and rRNA-ITS, NCBI asks for organism, a unique isolate, strain, clone, or specimen-voucher (or similar), collection-date, and geo_loc_name. The default template has only organism and moltype, and no doc lists the requirements.
  - Any template containing "isolate" is rejected, even for single-record samples, so labs whose unique ID is a field isolate number cannot use it.
  - For samples with more than one record RiboPilot appends `[isolate=<ID>_p<path>_s<scaffold>]`, which labels two scaffolds of one individual as two isolates (INSDC defines /isolate as the individual the sequence came from). The "GenBank rejects" rationale is unsourced.
  - The guide's modifier examples use underscores (`[specimen_voucher=...]`, `[collection_date=...]`), while NCBI's FASTA-modifier list uses hyphens (specimen-voucher, collection-date) and says to write names exactly as shown.
  - The computed tokens (seqid, rrna_definition, completeness, elements, element, product, sci_name, taxid) are not listed, and columns with spaces must be written in dotted form (`{Specimen.Voucher}`).
- Why it matters: Submissions bounce for missing modifiers, or carry a misleading isolate value.
- Fix: List NCBI's required modifiers, in NCBI's spelling, and the full token list in Step 9 and in the template help. Allow a user isolate. Mark variant records from one specimen with a note (for example "rDNA variant 2 of 2"), not with isolate.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-36](USABILITY_REVIEW.md#ux-36), [UX-40](USABILITY_REVIEW.md#ux-40)

<a id="doc-45"></a>**DOC-45. The default [organism=] uses the raw Taxon label, not the name confirmed in the Taxonomy step.** Major. Confirmed. ids: D-domain-5, D-intermediate-7 (part). Seen by 2 (docs-domain, docs-intermediate).
- Where: `R/export_opts.R:18, 118-122`; `README.Rmd:44`; `R/tax_search.R:61, 80`.
- Problem: The default header is `{seqid} [organism={Taxon}] ...`. Corrections made in the Taxonomy step update `sci_name`, not `Taxon`, so typos, synonyms, and informal names still go into `[organism=]`.
- Why it matters: It contradicts the README's "verify or correct it in the app", and wrong organism names reach public records.
- Fix: Default to `{sci_name}` once taxonomy is locked, or explain in Step 9 which column to use and when (for example an undescribed "Genus sp.").
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) DB-M1

<a id="doc-46"></a>**DOC-46. Step 9 does not say which files to upload, and README.txt names a table2asn / BankIt route without steps.** Minor. Confirmed. ids: D-novice-38, D-domain-23. Seen by 2 (docs-domain, docs-novice).
- Where: `vignettes/RiboPilot.Rmd:340-383`; README.txt text in `R/export.R:249-250`.
- Problem: NCBI's Submission Portal takes only `<group>.fasta` for rRNA-ITS; the `.tbl` files are for table2asn. Only the bundle's README.txt says so. README.txt then names an "alternate table2asn / BankIt route" but gives no steps. NCBI's GenBank submission page (updated 2026-06-15) lists only the Submission Portal as the web route and says table2asn output goes by email to gb-sub; whether BankIt still takes rRNA-ITS could not be confirmed.
- Why it matters: Users may upload the wrong files, or follow a route that is not described.
- Fix: Add one paragraph to Step 9: upload `<group>.fasta` at submit.ncbi.nlm.nih.gov; use the `.tbl` files only with table2asn, run it with `-V v`, read the `.val` file, and send the `.sqn` the way NCBI's page says. Drop the BankIt mention unless it is confirmed.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-47"></a>**DOC-47. The export layout block and the export_bundle help have small errors.** Minor. Confirmed. ids: D-appclaims-12, D-intermediate-29, D-novice-40, D-api-15, O-practice-25, O-practice-27, O-practice-33. Seen by 5 (docs-api, docs-appclaims, docs-intermediate, docs-novice, op-practice).
- Where: `vignettes/RiboPilot.Rmd:349-365`; `man/export_bundle.Rd:50`; `R/export.R:518, 611, 643, 659-693, 725-731`; `R/itsx_check.R:424-426`.
- Problem:
  - The 5.8S folder is `markers/58S/` (files `*_58S.fasta`). `<group>_itsx_check.csv` is written only when ITSx ran, and markers, the report, and the ITSx check are optional ticks.
  - `export_project()` writes `<group>.zip`. The help says `<out_dir>/<label>.zip`, but it goes to `<out_dir>/export/<label>.zip`, and with a relative path (including the default ".") zip fails right after "Exported N unit(s)" and the function returns normally with no zip.
  - The ITSx check runs after the FASTA and `.tbl` are written, not "before the bundle is written". Re-export deletes the group folder first.
  - The bundle README.txt's file list omits `_itsx_check.csv` and `itsx_check/`.
  - The summary CSV's variant_class, variant_group, variant_identity, and n_variants columns are not explained in the guide or in README.txt.
- Why it matters: Users cannot match the guide to the files they get.
- Fix: Correct the tree and the help, mark the optional items, mention the zip and the overwrite, normalize the path before `setwd()` and check the zip status, update README.txt's file list, and add a column key for the summary CSV.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-38](USABILITY_REVIEW.md#ux-38), [UX-39](USABILITY_REVIEW.md#ux-39); [CODE_REVIEW.md](CODE_REVIEW.md) DB-L2, DB-L3

<a id="doc-48"></a>**DOC-48. The ITSx notes in README.txt and Step 9 overstate: "ITS too long" is called misassembly, and NCBI is said to remove records "silently".** Minor (D-domain-15 polish, and one of its points was refuted; see the Problem). Confirmed. ids: D-domain-14, D-domain-15. Seen by 1 (docs-domain).
- Where: README.txt text in `R/export.R:237, 260-285`; `vignettes/RiboPilot.Rmd:364-369`.
- Problem:
  - Stock ITSx flags any ITS region over 1,500 bp. The human reference repeat U13369 (ITS1 1,095 bp, ITS2 1,154 bp) trips it, and the fish ruleset allows ITS1 and ITS2 up to 1,200 bp each. README.txt calls the flag "misassembly", and the Get started list of ITSx flags leaves it out.
  - README.txt says "NCBI's ITSx evaluation may silently remove sequences", but NCBI's help says "You will be alerted to any changes that are made to your sequences".
  - Get started says NCBI runs ITSx with "default parameters and no organism group"; NCBI names ITSx but not its parameters, so this needs a source.
  - The worked example accession, PZ692210, is Mallomonas bessudovae, a chrysophyte alga, and its definition line is garbled.
  - (Refuted part of D-domain-15: in ITSx 1.1.3 a 5.8S overlap with SSU or LSU is marked "Chimeric!", so README.txt's itsx_chimeric wording is right.)
- Why it matters: Correct mammal and bird assemblies are labelled misassembled, and users fear a silent removal that NCBI says it reports.
- Fix: Describe itsx_ITS_too_long as "ITS1-5.8S-ITS2 longer than 1,500 bp: normal in mammals, birds, and some other groups; the Portal may still reject it", and add it to the list. Replace "silently" with "may remove (NCBI tells you)". Source the parameter claim or drop it. Use an animal example accession.

<a id="doc-49"></a>**DOC-49. Every N-run in the .tbl claims linkage_evidence "align genus".** Minor. Confirmed. ids: D-domain-16. Seen by 1 (docs-domain).
- Where: `R/export.R:136-141` (the `.tbl` writer); no doc mentions assembly_gap features.
- Problem: "align genus" (AGP 2.1) says the gap was bridged by aligning to a same-genus reference genome. RiboPilot's N-runs come from assembler scaffolding or from user assemblies.
- Why it matters: It is a false provenance statement in public records.
- Fix: Use "paired-ends" where the assembler scaffolded with read pairs, otherwise "unspecified", and document it in Step 9.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) DB-L11

Also on this page: [DOC-17](#doc-17), [DOC-18](#doc-18), and [DOC-19](#doc-19) (install and practice commands), [DOC-23](#doc-23) (diagram), [DOC-15](#doc-15) (terms), [DOC-16](#doc-16) (conda notes), [DOC-94](#doc-94) (curly quotes), [DOC-98](#doc-98) (page length), [DOC-60](#doc-60) (the NCBI lookup, lines 264-268), [DOC-01](#doc-01), [DOC-03](#doc-03), [DOC-08](#doc-08), [DOC-09](#doc-09), [DOC-10](#doc-10), and [DOC-13](#doc-13) (missing content), [DOC-90](#doc-90) (locking), [DOC-93](#doc-93) ("copies"), [DOC-53](#doc-53) (flags panel).

### Reference databases (vignettes/Reference-Databases.Rmd)

<a id="doc-50"></a>**DOC-50. The 25,000-record ceiling is described as a hard limit.** Minor. Confirmed. ids: D-novice-24, D-domain-20, D-editor-24. Seen by 3 (docs-domain, docs-editor, docs-novice).
- Where: `vignettes/Reference-Databases.Rmd:41-44`; `R/reference_groups.R:160-188`.
- Problem: The walk starts at family and stops at the first clade with 1,000 or more records, even if it has more than 25,000 (a warning shows). The maximum matters only when no clade reaches 1,000.
- Why it matters: Users expect a database capped at 25,000 records and get a larger, slower one.
- Fix: "RiboPilot stops at the first clade with at least 1,000 records; 25,000 is a fallback limit, used only when no clade reaches the minimum."
- See also: [UX-12](USABILITY_REVIEW.md#ux-12)

<a id="doc-51"></a>**DOC-51. Reference tab details differ from the app.** Minor. Confirmed. ids: D-novice-25, D-appclaims-17, O-novel-10, D-editor-24 (part). Seen by 4 (docs-appclaims, docs-editor, docs-novice, op-novel).
- Where: `vignettes/Reference-Databases.Rmd:38-81`; `R/app_ui.R:60-76`; `R/app_reference.R:35-46, 132-155, 194, 201, 465-470`; `R/reference_build.R:263-269, 307-345`.
- Problem:
  - The tab opens on a Samples view (Move, Keep together, Release), not on group rows.
  - Propose groups opens a review window that needs Apply.
  - There is a sixth status, "Database missing". The button is "Run from App", not "Run".
  - The guide says the Set up database window offers "four choices"; for several groups at once it shows three (the clade choice appears only for a single group).
  - The clade ladder is "deepest first", not broadest first.
  - A finished build also re-queues the group's unlocked samples for assembly.
  - "Leave these samples on the built-in fish database" keeps an already-wired group on its clade database.
- Why it matters: Users look for buttons and orders that differ from what they see.
- Fix: Describe the Samples/Groups toggle, Apply, all six statuses, the exact button labels, the ladder order, and the re-queue.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) REF-L12

<a id="doc-52"></a>**DOC-52. The Reference tab is called optional, which non-fish users will take literally.** Minor (lowered from major by the verifier). Confirmed. ids: D-novice-20. Seen by 1 (docs-novice).
- Where: `vignettes/Reference-Databases.Rmd:91-94` vs `33-36`.
- Problem: Skipping it assembles non-fish samples against the fish reference, which the same page calls "usually the wrong choice outside fish". The app only warns.
- Why it matters: Non-fish users get worse assemblies without knowing why.
- Fix: "Optional for ray-finned fishes. For any other group, do this step before assembling."
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

Also on this page: [DOC-36](#doc-36) ("you can keep working in the app while it runs"), [DOC-69](#doc-69) ("through a remote browser"), [DOC-05](#doc-05) (the default reference folder), [DOC-32](#doc-32) (the console build), [DOC-18](#doc-18) (the `new_project()` calls), [DOC-91](#doc-91) (the pre-filled API key), [DOC-94](#doc-94) (curly quotes), [DOC-93](#doc-93) (line 191).

### Curation flags (vignettes/Curation-Flags.Rmd)

<a id="doc-53"></a>**DOC-53. The flag table misses flags users see, and misstates what happens without a lookup.** Major (O-novel-9 rated minor). Confirmed. ids: D-novice-14, D-domain-9, D-editor-11, D-appclaims-7, D-intermediate-19, O-novel-9 (part). Seen by 6 (docs-appclaims, docs-domain, docs-editor, docs-intermediate, docs-novice, op-novel).
- Where: `vignettes/Curation-Flags.Rmd:26-45`; `R/curate_rrna.R:55-56, 71, 82, 113-143, 262-266, 334-338`; `R/app_server.R:2046-2073`.
- Problem:
  - Missing codes: truncated (an anchor under 90 percent of its model), contaminant_unchecked, not_nuclear_rdna, and the no_clade_ruleset note (on every gastropod unit in op-novel), plus the chimera and ITSx warnings.
  - "those two checks are simply skipped" is wrong: a missing lookup raises contaminant_unchecked by design.
  - The lookup is the bundled panel with NCBI escalation, not "NCBI"; contaminant means a different phylum, not an "unrelated group".
  - The count shows sample-level flags only: SRR22396758 shows 0 warnings while its 28S carries "length,truncated", and the panel says "No curation warnings" above a table of feature warnings. Only feature-level warnings get a "select" link.
- Why it matters: Users cannot tell which warnings matter or what to do about them.
- Fix: List every emitted code with its threshold and a "what to do" column. Fix the sentences. Count feature flags too, or label the count "sample-level".
- Status: fixed on `review/2026-09-full-review` (fourth batch, docs: the warnings table lists every code and says what happens without a lookup). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).
- See also: [UX-51](USABILITY_REVIEW.md#ux-51)

<a id="doc-54"></a>**DOC-54. Five rulesets ship and are picked from the lineage; the docs describe two.** Major. Confirmed. ids: D-novice-16, D-domain-7, D-editor-12, D-appclaims-10, O-novel-9, D-intermediate-19 (part). Seen by 6 (docs-appclaims, docs-domain, docs-editor, docs-intermediate, docs-novice, op-novel).
- Where: `README.Rmd:56`; `vignettes/Curation-Flags.Rmd:172-194`; `R/params_rrna.R:12-50, 94-139, 157-184`; `R/curate_rrna.R:249-269`; `R/app_opts.R:296-304`.
- Problem:
  - metazoa_generic, actinopterygii, asteroidea, scyphozoa, and octocorallia ship, and `resolve_ruleset()` picks one from each sample's lineage (choosing metazoa_generic in the options is overridden).
  - The fish ITS range (150-1200) is narrower than generic (120-1600), not "kept wide".
  - The per-gene identity floors (95/93/80 generic, 97/95/85 fish, 93/90/78 invertebrate), the 50 percent coverage gate, and the 90 percent truncation cutoff are not documented.
  - "28S 2000-6000 (the recovered portion varies with sequencing depth)" is the wrong reason: complete 28S length does not depend on depth, and the truncated flag covers incomplete genes.
  - The docs never say where the ruleset is set (Curate Options, Clade ruleset).
- Why it matters: Non-fish users cannot tell which thresholds were applied to their samples.
- Fix: Table all five rulesets (length ranges, identity floors, lineage match), explain automatic selection and the override, fix the fish wording and the 28S rationale, and update the README sentence.
- Status: partly fixed on `review/2026-09-full-review` (docs rewrite): the README and the Curation guide name the five rulesets and say they are picked from the lineage. Per-ruleset thresholds are still undocumented.
- See also: [UX-51](USABILITY_REVIEW.md#ux-51)

<a id="doc-55"></a>**DOC-55. The triage section counts seven rules; eight ship.** Major. Confirmed. ids: D-novice-15, D-editor-13, D-domain-19, D-appclaims-11. Seen by 4 (docs-appclaims, docs-domain, docs-editor, docs-novice).
- Where: `vignettes/Curation-Flags.Rmd:69-111`; `R/triage_rules.R:16-114`.
- Problem: Rule 8, "Likely chimeric assembly" (off or flag only), is missing, so four rules only flag, not three. "Each rule can be turned off" never says where (Assemble step, select samples, Triage button).
- Why it matters: The rule table does not match the Triage window.
- Fix: Add rule 8, fix the counts, and name the Triage button and window.
- Status: rule 8 added to the triage table on `review/2026-09-full-review` (docs rewrite).

<a id="doc-56"></a>**DOC-56. Chimera detection is not documented, and the in-app description is incomplete.** Major. Confirmed. ids: D-domain-10, D-editor-13 (part), D-novice-15 (part), D-intermediate-19 (part). Seen by 4 (docs-domain, docs-editor, docs-intermediate, docs-novice).
- Where: no vignette section (only an ITSx mention at `RiboPilot.Rmd:367`); in-app panel `R/app_modals.R:104-108`; `R/chimera.R`; `R/export.R:200-225`.
- Problem: The detector has four channels (per-gene LCA, windowed BLAST, Ribovore ribotyper structure codes, and depth steps plus clip pileups), three labels (clean, suspect, likely), a taxonomy-only cap, and triage rule 8. Both suspect and likely become "chimeric(...)" in the export's submission_risk column, next to gap_in_rDNA and misassembly, which are also undocumented. The panel text covers windows and clips only, and its "two organisms" definition leaves out misjoins between array copies of one genome.
- Why it matters: A label that marks export records "chimeric" is unexplained.
- Fix: Add a "Chimeras and submission risk" section: channels, scores, labels, effects on triage and export, and known false-positive sources (next item). Update the panel text.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-57"></a>**DOC-57. The chimera check scores a clean flatfish the same as the planted chimera, so a good record is exported as chimeric.** Major. Confirmed. ids: D-domain-11. Seen by 1 (docs-domain).
- Where: chimera reason text in the app and in the export's submission_risk; `R/chimera.R`; README.txt text from `R/export.R:267`.
- Problem: SRR33849683 (Ancylopsetta cycloidea) scores 0.6 (suspect) with "18S resolves to Chordata ... but 5.8S to Pleuronectiformes"; Pleuronectiformes is inside Chordata (the pair's rank distance is 5, which fires the test). The window channel on the same scaffold switches to "Mammalia (PREDICTED: Canis, 100.0%)" over 751-2000, a tie in conserved 18S with the neighbouring Paralichthys windows. The synthetic fish-plus-Penicillium chimera CHIM01 also scores 0.6 (suspect); neither reaches "likely", so triage rule 8 never fires on the known chimera. The export then labels the genuine flatfish record chimeric, and README.txt warns that NCBI may remove chimeric sequences.
- Why it matters: It pushes users to drop good records. Screenshots docs_11 and docs_15 show this "suspect" call.
- Fix: Treat an ancestor/descendant LCA pair as agreement; break window ties toward the sample's lineage, or discount windows that sit inside conserved anchors; re-check the label cut-offs against CHIM01; document this false-positive source.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-82](USABILITY_REVIEW.md#ux-82), [UX-83](USABILITY_REVIEW.md#ux-83); [CODE_REVIEW.md](CODE_REVIEW.md) SCI-M3, SCI-M4

<a id="doc-58"></a>**DOC-58. The isomer test ignores the ITS, so real ITS variants are told the pick "does not matter".** Major. Confirmed. ids: D-domain-2. Seen by 1 (docs-domain).
- Where: `vignettes/Curation-Flags.Rmd:127-135, 144-146`; `R/screen_variants.R:1-19, 175-205`.
- Problem: Units whose anchors (18S, 5.8S, 28S) match at 99.5 percent or more are called isomers, and the guide says picking between isomers "does not matter either way". Intragenomic variation sits mostly in ITS1, ITS2, and IGS, so units with identical anchors but different ITS are classic ITS variants, and for an ITS barcode record the choice matters a lot.
- Why it matters: Users may submit an arbitrary ITS variant believing it is equivalent.
- Fix: Compare ITS1 and ITS2 as well. Call units with different spacers "ITS variants" and say the choice changes the submitted ITS. Keep "isomer, pick either" for units whose annotated 18S-28S span is identical.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) SCI-M7

<a id="doc-59"></a>**DOC-59. "Cannot be contamination" contradicts the phylum-level screen described on the same page.** Major. Confirmed. ids: D-domain-3. Seen by 1 (docs-domain).
- Where: `vignettes/Curation-Flags.Rmd:120-124` vs `94-97`; `R/screen.R:10-13, 184, 244`.
- Problem: The guide says on-target variants "cannot be contamination", but the screen tests agreement at phylum rank, and the same page says "a human on a fish sample passes". Same-phylum contaminants (a human, gut-content prey, a parasite) are on-target, and the guide presents them as paralogs that "still export".
- Why it matters: Users skip checking variants that may be contaminants.
- Fix: Say on-target only means the same phylum, and tell users to BLAST each variant's ITS and check it against their taxon before locking.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-86](USABILITY_REVIEW.md#ux-86)

<a id="doc-60"></a>**DOC-60. The "GenBank screen" is by default a local panel, but the docs describe an NCBI search.** Major (D-domain-6 major; the other three rated minor). Confirmed. ids: D-domain-6, D-appclaims-14, D-intermediate-24, D-editor-27. Seen by 4 (docs-appclaims, docs-domain, docs-editor, docs-intermediate).
- Where: `vignettes/Curation-Flags.Rmd:43-45, 55-57, 78`; `vignettes/RiboPilot.Rmd:264-268`; column help `R/constants.R:149-150`; `man/curate_rrna_core.Rd:37`; `man/parse_blast_candidates.Rd:41`; defaults in `R/init_db.R:325-329` and `R/app_opts.R:182-194`.
- Problem:
  - The guide calls it "the GenBank screen, the best reference match found remotely", and the column help says "searched against all of GenBank". By default (remote_blast = 0, remote_fallback = 1) the screen searches a bundled local panel of metazoan nuclear rDNA plus NCBI's 16S, SSU_eukaryote, and LSU_eukaryote decoy sets, and goes to NCBI core_nt only when the panel finds no significant hit. In the practice project 12 units were answered by the local panel and 2 by the seeds; none was searched remotely.
  - The conda route has no panel, so every search goes to NCBI.
  - Get started says a failed lookup leaves results "not affected". In fact units get contaminant_unchecked and the triage flag "GenBank screen has not run".
  - The Hit Source column, the BLAST Options "Search" choice, and `RIBOPILOT_NO_LOCAL_BLAST` are never explained.
- Why it matters: A methods section written from the docs ("screened against GenBank core_nt") would be wrong, and the contaminant and not_nuclear_rdna flags depend on which database answered.
- Fix: Call it the "rDNA screen", and say what it searches by default, when it goes to NCBI, and how to force a remote search. Explain the Hit Source column, show which database answered next to each verdict, fix the column help and the two man pages, and describe what an offline run looks like.
- Status: fixed on `review/2026-09-full-review` (fourth batch, docs: "BLAST screen" and "bundled rDNA database" in the app and docs). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).

<a id="doc-61"></a>**DOC-61. The depth explanation does not match how depth is computed.** Minor. Confirmed. ids: D-domain-18. Seen by 1 (docs-domain).
- Where: `vignettes/Curation-Flags.Rmd:161-170`; triage window `R/app_modals.R:277-280`; column help `R/constants.R:134-136`; `inst/nextflow/modules/coverage.nf:21-54`.
- Problem: Reads are mapped one path at a time (one bowtie2 index per path), and multi-mappers (MAPQ < 10) are dropped, so reads are not "split across paths"; they are dropped wherever sequence repeats inside one path (shared anchors, IGS subrepeats). The Triage window's "isomer paths splitting reads" is stale in the same way. The chimera read channel and the trim tool's suggested junctions read these same dips.
- Why it matters: Users read repeat-driven dips as misjoins.
- Fix: Explain the per-path mapping and the MAPQ filter, and name repeat-driven dips as an expected false signal.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) SCI-M2

<a id="doc-62"></a>**DOC-62. The rule 7 note is about the wrong thing and contradicts the tooltip.** Minor. Confirmed. ids: D-editor-26. Seen by 1 (docs-editor).
- Where: `vignettes/Curation-Flags.Rmd:113-116`; `R/triage_rules.R:89-100`; `R/constants.R:154-155`.
- Problem: The note talks about low-identity hits at a large taxonomic distance, but rule 7 (odd_one_out) is disagreement with sibling scaffolds. The app tooltip calls sibling disagreement "the usual signature of contamination".
- Why it matters: The guide and the app give opposite advice on the same flag.
- Fix: Rewrite the note around sibling disagreement, and align the tooltip with it.

<a id="doc-63"></a>**DOC-63. The page explains what flags mean but not what to do, and the triage prose is written for experts.** Minor. Confirmed. ids: D-novice-27. Seen by 1 (docs-novice).
- Where: `vignettes/Curation-Flags.Rmd:47-116`.
- Problem: Evidence channels, core_nt, feature model, no-op, isomer, and sibling scaffolds are used without definitions, and no flag has a suggested action.
- Why it matters: Novices cannot act on the page.
- Fix: Lead with a "what to do" column, move rule internals to an "In depth" section, and link glossary terms.

<a id="doc-64"></a>**DOC-64. Spaced hyphens are used as dashes.** Polish. Confirmed. ids: D-novice-42, D-editor-34. Seen by 2 (docs-editor, docs-novice).
- Where: `vignettes/Curation-Flags.Rmd:66, 71, 79-83, 96`.
- Problem: " - " stands in for a dash (for example "corroboration - see the table below", "Yes - dropped by default").
- Why it matters: It goes against the house style.
- Fix: Use a comma, colon, semicolon, or two sentences.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-65"></a>**DOC-65. The curation tables use equal column widths.** Polish. Confirmed. ids: D-editor-35. Seen by 1 (docs-editor).
- Where: `vignettes/Curation-Flags.Rmd:32-41, 75-83`.
- Problem: Equal dash counts give each column the same width, so "What it means" is cramped at 1400 px and cut off at 390 px.
- Why it matters: The key column is hard to read, worst on phones.
- Fix: Widen the last column's dashes, or shorten the other cells.

Also on this page: [DOC-15](#doc-15) (Curation flags versus Curation warnings), [DOC-12](#doc-12) (references for the paralog claim), [DOC-87](#doc-87) (the Curation options controls), [DOC-42](#doc-42) (partial versus complete).

### Function reference (man pages, the pkgdown reference index, and console messages)

<a id="doc-66"></a>**DOC-66. The reference index is one flat list of 54 topics, most of them pipeline internals.** Major. Confirmed. ids: D-editor-14, D-api-8, D-api-28, D-api-11, D-novice-35, D-intermediate-30. Seen by 4 (docs-api, docs-editor, docs-intermediate, docs-novice).
- Where: `_pkgdown.yml` (no `reference:` section); `NAMESPACE` (53 exports); `man/*.Rd`; `man/check_container_image.Rd:11` (`R/utils.R:66`); `man/auto_triage_sample.Rd:26`, `man/new_project.Rd:46`, `man/triage_preview.Rd:38` (dev/design links).
- Problem: 54 topics sit under one "All functions" heading. Only 12 exports are for users (new_test_project, new_project, RiboPilot, run_app, build_reference_db, ribopilot_ref_dir, generate_config, hydra_setup, export_project, export_bundle, run_assembly, run_fallback). 24 are exported only because the Nextflow scripts call them by bare name (annotate_records, chimera_unit, write_coverage_units, and so on), and 17 have no caller outside R/ and the tests. Titles use internal codes ("WF1", "Option B2 baiting", "Run every channel for one unit and write chimera + chimera_windows"), and three pages cite unpublished `dev/design` notes. Only hydra_setup, of 54 pages, has an example. check_container_image is not exported but gets a page, because its `@noRd` sits inline in `@param` text, where it shows up as literal text.
- Why it matters: A biologist cannot find the dozen functions they need.
- Fix: Add grouped reference sections (Start a project; Reference databases; Clusters; Export; Advanced scripting). Put `@keywords internal` on the 24 Nextflow helpers (or switch the scripts to `RiboPilot:::` and unexport them), and unexport the 17 with no caller. Tested on a copy, `check_pkgdown()` then reports "No problems found". Move check_container_image's `@noRd` to its own line and drop the dev/design links. Add examples to the user-facing pages: runnable ones for `ribopilot_ref_dir()`, `generate_config(profile_dir = tempdir())`, and `new_project()` in a temp folder with `Rproj = FALSE`, and `\dontrun` for the rest. The outline section below lists which function goes where.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L5

<a id="doc-67"></a>**DOC-67. new_project() help says NMNH_Hydra and NOAA_SEDNA need generate_config() and calls conda ready to use.** Major. Confirmed. ids: D-hpc-9, D-intermediate-18, D-api-3, D-editor-15 (part). Seen by 4 (docs-api, docs-editor, docs-hpc, docs-intermediate).
- Where: `man/new_project.Rd:46-56` (`R/init_project.R:17-27`).
- Problem: The executor text says every scheduler, including "NMNH_Hydra" and "NOAA_SEDNA", "is a template, not a usable executor", but those two hold only tokens `new_project()` fills and work directly (README.Rmd:66 says so), and `generate_config()` cannot even build them. The page calls conda (the default) "Ready to use as-is", does not explain the container tag rule or that `container` can be a SIF path, and cites `dev/design/20_USER_ASSEMBLIES.md`, which is not shipped. Get started sends cluster users to this page.
- Why it matters: A Hydra user would build a generic SGE profile instead of using the tuned built-in one.
- Fix: List NMNH_Hydra and NOAA_SEDNA as ready-made Singularity profiles; say generic schedulers need `generate_config()`; mark conda unsupported; document `container` (the default tag, a local tag, or a SIF path); remove the dev/ reference.

<a id="doc-68"></a>**DOC-68. run_assembly() help says annotation waits for a second run after locking.** Major (D-appclaims-9 rated minor by its verifier). Confirmed. ids: D-api-4, D-intermediate-21, D-appclaims-9, D-editor-15, D-novice-30. Seen by 5 (docs-api, docs-appclaims, docs-editor, docs-intermediate, docs-novice).
- Where: `man/run_assembly.Rd:34-37` (`R/run_assembly.R:40-48`).
- Problem: "Annotation and curation are handled by re-running this workflow after you review and lock the assemblies in the app" contradicts the function's own success message ("Every scaffold is annotated as it is assembled") and Get started. The page also says it "runs fastp + GetOrganelle", leaving out the fallback assemblers and the supplied-assembly search.
- Why it matters: Users re-run for nothing, or think annotation is missing.
- Fix: Say one run preprocesses, assembles, scores, BLASTs, maps coverage, annotates, curates, screens, checks for chimeras, and triages every unit; re-run (or press Update) only after changing a selection, trim, or option.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-69"></a>**DOC-69. The headless (job-script) rule is described wrongly, and launch.browser = FALSE silently removes Run from App.** Major (D-novice-26 and D-intermediate-26 rated minor). Confirmed. ids: D-api-7, D-novice-26, D-intermediate-26. Seen by 3 (docs-api, docs-intermediate, docs-novice).
- Where: `man/RiboPilot.Rd:19-23`; `man/run_app.Rd:25-28`; `vignettes/Reference-Databases.Rmd:63-64` ("If you are using RiboPilot through a remote browser"); `R/app_run.R:16-21, 150-197`; `R/app_modals.R:426-460`.
- Problem: The help says headless mode is "Detected automatically (RStudio Server or a batch scheduler on PATH means yes)", and the Reference guide says "through a remote browser". In the code, `detect_headless()` also returns TRUE whenever `launch.browser` is not TRUE, and `launch.browser` defaults to `interactive()`. So on a plain workstation, `launch.browser = FALSE` (which the SSH-tunnel route needs) or a launch from Rscript gives the headless run window: no Run from App button, only Save script or Submit as job, and the text "No scheduler (sbatch/qsub/bsub) found on PATH ... running from the login node is discouraged". A workstation with sbatch on PATH also gets job-script mode. Neither page mentions the `headless = TRUE/FALSE` override.
- Why it matters: Running from the app is the intended workflow, and a user who opens the app without a browser launch loses the Run button without knowing why.
- Fix: Reword both help pages and the guide: "The app is headless when it does not open a browser (launch.browser = FALSE or a non-interactive session), on RStudio Server, or when a scheduler is on PATH. Headless mode hides Run from App; pass headless = FALSE to run on this machine."
- Status: fixed on `review/2026-09-full-review`: the cluster package changed the headless rule (a scheduler on PATH or RStudio Server; `launch.browser` no longer matters), and Running on a cluster and Reference databases describe it.
- See also: [UX-79](USABILITY_REVIEW.md#ux-79)

<a id="doc-70"></a>**DOC-70. A failed new_project() leaves a half-made project, and the retry silently skips the config.** Major (D-api-14 rated minor). Confirmed. ids: D-intermediate-17, D-api-14. Seen by 2 (docs-api, docs-intermediate).
- Where: `man/new_project.Rd` (executor, force); `R/init_project.R:106-111, 196-203`.
- Problem: With `executor = "slurm"` (no profile), `.sqlite` is written before the config check fails. Running again with a valid profile prints "Database already exists. Use force = TRUE ..." and returns without writing `.config`, so the project cannot run. Neither the page nor the placeholder error mentions `force = TRUE`.
- Why it matters: One typo leaves the user stuck with no clear way out.
- Fix: Check the config before creating the database, and mention `force = TRUE` in the error and the help.
- See also: [UX-02](USABILITY_REVIEW.md#ux-02); [CODE_REVIEW.md](CODE_REVIEW.md) DB-L4

<a id="doc-71"></a>**DOC-71. generate_config() ignores container_engine for "local" and queue for "awsbatch".** Minor (lowered from major by the verifier: a local Singularity profile is a niche setup). Confirmed. ids: D-hpc-8. Seen by 1 (docs-hpc).
- Where: `man/generate_config.Rd` (container_engine, queue); `R/generate_config.R:151-155`; `inst/config.local`; `inst/config.awsbatch`.
- Problem: The local and awsbatch templates have no `<<CONTAINER_ENGINE>>` token, so `generate_config(scheduler = "local", container_engine = "singularity")` silently saves a Docker config. For awsbatch, `queue` is dropped and a warning asks the user to fill BATCH_QUEUE by hand. The help promises both work.
- Why it matters: A Docker-less big server or interactive node cannot get a Singularity local profile.
- Fix: Add the engine token to `config.local` (or stop with an error), map `queue` to BATCH_QUEUE, and say in the help which arguments apply to which scheduler.

<a id="doc-72"></a>**DOC-72. Roxygen markdown is off: raw backticks and [fn()] brackets on every page, and no cross-links.** Minor. Confirmed. ids: D-editor-18, D-api-10, D-intermediate-30 (part), D-hpc-16 (part). Seen by 4 (docs-api, docs-editor, docs-hpc, docs-intermediate).
- Where: `DESCRIPTION` (no `Roxygen: list(markdown = TRUE)`); for example `man/RiboPilot.Rd:19`, `man/new_project.Rd:54`, `man/generate_config.Rd:62`.
- Problem: 44 of 54 Rd files show literal backticks and 26 show literal `[fn()]`; no Rd has a `\link`. Turning markdown on (tested on a copy) creates 16 links from exported pages to unexported topics (detect_headless, ribopilot_config_dir, rrna_params, risk_flags, and others).
- Why it matters: The manual looks broken and cannot be navigated.
- Fix: Reword those 16 references as plain code, add `Roxygen: list(markdown = TRUE)`, and re-document with the rd roclet only (keeping the hand-maintained NAMESPACE).
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-M1

<a id="doc-73"></a>**DOC-73. generate_config() help has no example and does not say where profiles are saved or how to share one.** Minor. Confirmed. ids: D-hpc-16, D-api-17. Seen by 2 (docs-api, docs-hpc).
- Where: `man/generate_config.Rd`; `R/generate_config.R:10-12, 29-36, 290-296`.
- Problem: The save folder is given as `ribopilot_config_dir()`, which is unexported, so the path (`~/.config/R/RiboPilot`) is never shown. There is no advice on sharing a lab profile (`profile_dir` or `R_USER_CONFIG_DIR`, and passing `profile_dir` to `new_project()`). `container_run_options` is said to apply to Singularity/Apptainer only, but it also applies to Docker. Nothing says Hydra and SEDNA need no profile, or that awsbatch still needs hand edits. The saved profile keeps the template header "Turn it into a reusable profile with generate_config(...)", and the page opens "Mirrors MitoPilot".
- Why it matters: Get started sends cluster users here, and they leave without a working recipe.
- Fix: Add a slurm example followed by `new_project(executor = <name>)`, state the default folder, add the sharing and Hydra/SEDNA/awsbatch notes, fix the argument text, and rewrite the header when saving.
- Status: docs half fixed on `review/2026-09-full-review` (docs rewrite): Running on a cluster has a worked example, the save folder, and sharing with `profile_dir`. The help page example is still open.

<a id="doc-74"></a>**DOC-74. new_project() help omits ID rules, the return value, and the data loss from force, and typos pass silently.** Minor. Confirmed. ids: D-api-13, D-novice-30 (part). Seen by 2 (docs-api, docs-novice).
- Where: `man/new_project.Rd:33-87`; `R/init_db.R:63, 79-90`.
- Problem: There is no `\value` (it returns NULL). mapping_fn omits the ID rules and that R1/R2 are paths inside `data_path`. `force` does not say all results are deleted. `...` passes to `new_db()`, whose own `...` swallows typos (`assembel_cpus = 99` was accepted). `custom_*_db` does not say NULL means the fish default. The page cites `dev/design/20_USER_ASSEMBLIES.md` and opens "Mirrors MitoPilot".
- Why it matters: The main setup function hides its rules and its most destructive option.
- Fix: Add `\value`, the ID rules, a data-loss warning on `force`, and the default database. Drop the dev/ and MitoPilot references. Make `new_db()` error on unknown arguments.

<a id="doc-75"></a>**DOC-75. The RiboPilot() and run_app() pages lack a return value and explain defaults poorly.** Minor. Confirmed. ids: D-api-12. Seen by 1 (docs-api).
- Where: `man/run_app.Rd`; `man/RiboPilot.Rd:19-21`; `R/app_run.R:175-206`.
- Problem: Neither page has `\value` (both return a Shiny app object that starts when printed). `launch.browser` defaults to `interactive()`; `headless` is TRUE whenever `launch.browser` is not TRUE; the tunnel hint prints only when `port` is set. run_app's description repeats its title, and two pages cover one entry point.
- Why it matters: Remote users cannot work out from the help how to launch the app.
- Fix: Fold run_app into the RiboPilot page (`@rdname`) or hide it. Add `\value`, the defaults, the headless rule, and a `\dontrun` example with `port` for remote use.

<a id="doc-76"></a>**DOC-76. build_reference_db(): count_only reports success when tools are missing, and a retry printed the API key.** Minor. Confirmed. ids: D-api-16, D-intermediate-28. Seen by 2 (docs-api, docs-intermediate).
- Where: `man/build_reference_db.Rd:29-54`; `R/build_reference_db.R:69-76`; `inst/ref_dbs/getOrganelle/GenBankDownload/download_rrna.sh:117`.
- Problem: Without esearch on PATH, `count_only` prints "esearch: command not found" and still returns NA counts with no error. On an HTTP 429 retry, edirect echoed the full request, including the NCBI API key, to the console and log, and 4 of 51,311 records were lost while the build reported success. The page does not say the key falls back to `NCBI_API_KEY`, or plainly that this route is unsupported and the Reference tab is the supported way.
- Why it matters: Silent failures, and a key in logs that labs may share.
- Fix: Stop on a nonzero `count_only` status with a conda_env hint; report parsed versus expected counts; document the env var; call the route unsupported and point to the Reference tab; add a `\dontrun` example.
- See also: [UX-47](USABILITY_REVIEW.md#ux-47)

<a id="doc-77"></a>**DOC-77. Pipeline jargon in user-facing titles and messages.** Minor. Confirmed. ids: D-api-19. Seen by 1 (docs-api).
- Where: `man/run_assembly.Rd:5`; `man/run_fallback.Rd:5, 23, 41`; `R/run_assembly.R:70-72`.
- Problem: "WF1", "Option B2 baiting / Option A de-novo", "derived assemble_opts preset", and "self-contained DB I/O" are design-note terms. A missing `.config` gives 'file.exists(file.path(path, ".config")) is not TRUE'.
- Why it matters: Users cannot connect the help pages to the app's labels.
- Fix: Use plain names ("baited SPAdes", "de novo SPAdes"), and replace `stopifnot()` with a message naming the missing file.

<a id="doc-78"></a>**DOC-78. new_db() is exported with a stale description and misleading defaults.** Minor. Confirmed. ids: D-api-20. Seen by 1 (docs-api).
- Where: `man/new_db.Rd:61-85`; `R/init_db.R:102-103, 527-614`.
- Problem: It says annotation, curation, and export tables "are added in later build increments", but they already exist. "May be NULL until built" really means NULL is the fish URL. `...` is ignored. Users never call it directly.
- Why it matters: Noise in the manual, and the text is wrong.
- Fix: Unexport it or mark it internal, and fix the text.

<a id="doc-79"></a>**DOC-79. Console messages point to .config, conda, or maintainer scripts instead of the app.** Minor. Confirmed. ids: D-api-21, D-editor-28, D-intermediate-31, O-novel-7, D-novice-13 (part). Seen by 5 (docs-api, docs-editor, docs-intermediate, docs-novice, op-novel).
- Where: `R/init_project.R:199-220`; `R/utils.R:79-112`; `R/run_assembly.R:82`; `R/nextflow_version.R:160`.
- Problem: `new_project()` ends with "Please open and review the .config file ...", which no doc explains and which needs no edits for Docker; the next step (open the app) is never named. The unreachable-database warning first says "Fix it from the app's reference database step", then gives an Asteroidea `build_reference_db()` example (the conda route), even for a fish or snail project, so a snail user is told their fish database is broken. The container warning tells users to run `docker/deploy-dockerhub.sh`, which is not installed. `run_assembly()` says "WF1 complete ...", and the Nextflow hint suggests `mamba install -n ribopilot`.
- Why it matters: The first messages a user sees send them the wrong way.
- Fix: End with 'Next: RiboPilot::RiboPilot("<path>")'. Skip the fish-database check when the project will get its databases from the Reference step, or point it there with a neutral example. Point the container warning to the `container` argument. Drop WF1.
- See also: [UX-56](USABILITY_REVIEW.md#ux-56); [CODE_REVIEW.md](CODE_REVIEW.md) REF-L18

<a id="doc-80"></a>**DOC-80. The ribopilot_ref_dir() page lacks a return value, an example, and the real default path.** Minor. Confirmed. ids: D-api-18. Seen by 1 (docs-api).
- Where: `man/ribopilot_ref_dir.Rd:10-12`; `vignettes/Reference-Databases.Rmd:114`.
- Problem: There is no `\value` or example, and the default (`tools::R_user_dir("RiboPilot", "data")/ref_dbs`, in the home directory) is never stated. There is no HPC note on `RIBOPILOT_REF_DIR`.
- Why it matters: Cluster users fill their home quota without knowing why.
- Fix: Add `\value`, a runnable example, the Linux and macOS paths, and an HPC note; name the function in the Reference guide.

<a id="doc-81"></a>**DOC-81. roxygen2 warnings from glued comment blocks.** Polish. Confirmed. ids: D-editor-40, D-api-27 (part). Seen by 2 (docs-api, docs-editor).
- Where: `R/app_modals.R:723`; `R/app_run.R:68`; `R/generate_config.R:55`; `R/screen.R:96`.
- Problem: "@noRd must not be followed by any text": each `@noRd` sits before the title text, so text attaches to the wrong object.
- Why it matters: Build noise, and help text can land on the wrong function.
- Fix: Move `@noRd` to the end of each block, or separate the blocks.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L4

<a id="doc-82"></a>**DOC-82. Small help-page and metadata issues.** Polish. Confirmed. ids: D-api-27. Seen by 1 (docs-api).
- Where: `man/new_project.Rd:13`; `man/chimera_unit.Rd:38`; `DESCRIPTION`; `NAMESPACE` header; `R/zzz.R:36`.
- Problem: Usage renders `min_depth = 2e+06`. chimera_unit's page says it returns "(invisibly)", but it returns visibly. DESCRIPTION sets LazyData with no `data/` folder. There is no package-level help topic, so `?RiboPilot` opens the app launcher's page. The NAMESPACE header says not to run `document()`, but roxygen reproduces it exactly. The startup message prints "{RiboPilot}" in bold white, which is invisible on light terminals. (The other parts of this finding are in the README-sync, roxygen-warning, and DESCRIPTION items.)
- Why it matters: Small polish before the manual is public.
- Fix: Write the default as `2000000L`, fix the chimera_unit value text, drop LazyData, add a `"_PACKAGE"` doc, update the NAMESPACE header, and use a startup colour that shows on light and dark terminals.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PKG-L2, PKG-L5

Also on this page: [DOC-18](#doc-18) (the Rproj text), [DOC-21](#doc-21) (`hydra_setup`), [DOC-35](#doc-35) (`new_test_project` "10-sample"), [DOC-40](#doc-40) (`export_project` is a bare title), [DOC-47](#doc-47) (`export_bundle` zip path), [DOC-19](#doc-19) (`new_db` default URLs 404 while private).

### docker/README.md and the config templates

<a id="doc-83"></a>**DOC-83. docker/README.md examples use tags that do not exist, trigger a false warning, and give stale mount advice.** Minor. Confirmed. ids: D-hpc-14, D-api-24, D-editor-38, D-intermediate-35, D-novice-44, D-appclaims-25. Seen by 6 (docs-api, docs-appclaims, docs-editor, docs-hpc, docs-intermediate, docs-novice).
- Where: `docker/README.md:4-6, 24-25, 52-53, 88-103`; `R/utils.R:67-84` (`check_container_image`).
- Problem:
  - The `new_project()` example uses `container = "macguigand/ribopilot:v0.1.0"`; published tags are bare versions (0.0.0.9000 to 0.0.0.9018, no `latest`, no `v` tags), and the same file explains that policy.
  - The local-tag example (`ribopilot:0.0.0.9004`) names another version, and `new_project()` looks a namespace-less tag up on Docker Hub as `library/<name>`, so it warns falsely that the image "is not published on Docker Hub" and points to a maintainer script.
  - The advice to add `-v` mounts for outside reference databases is out of date: `new_project()` and `nf_command()` now add read, assembly, and seed/label mounts automatically.
  - The tool list leaves out pybarrnap, the default anchor tool.
- Why it matters: Copying the example asks for a tag that does not exist, and users edit mounts they do not need.
- Fix: Use `macguigand/ribopilot:<packageVersion>` (for example 0.0.0.9018) in examples, say mounts are added automatically, and add pybarrnap to the tool list. Skip the Docker Hub check for a namespace-less tag, or run `docker image inspect` first.

<a id="doc-84"></a>**DOC-84. The generated .config has the mount list pasted into a comment.** Polish. Confirmed. ids: D-novice-43, D-intermediate-34, D-api-24 (part). Seen by 3 (docs-api, docs-intermediate, docs-novice).
- Where: `inst/config.docker:11` (a comment that contains the `<<EXTRA_MOUNTS>>` token).
- Problem: The token in the comment is filled too, so line 11 of every generated `.config` has the `-v ...` mount list in the middle of a sentence, and `new_project()` then asks users to review that file.
- Why it matters: The one file users are told to read looks broken.
- Fix: Reword the comment so it does not contain the literal token.

<a id="doc-85"></a>**DOC-85. Every Singularity task log starts with duplicate-bind and RIBOPILOT_SRC noise.** Minor. Confirmed. ids: D-hpc-20. Seen by 1 (docs-hpc).
- Where: generated Singularity profiles; task `.command.err`; `inst/nextflow/bin/rp_env.R`.
- Problem: The profile binds `${params.rawDir}` and `new_project()` adds the same folder again, so every task warns "destination is already in the mount point list". `rp_env.R` also prints a RIBOPILOT_SRC note on every task, even for normal installed packages (92 and 35 such lines in one practice run).
- Why it matters: Users looking in `.command.err` for a real error see noise first.
- Fix: The rawDir fix in the "Already have an assembly?" item removes the duplicate bind. Print the RIBOPILOT_SRC note only for a source checkout.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) PIPE-M15

<a id="doc-86"></a>**DOC-86. The Docker Hub page is blank, and config.local has a garbled sentence.** Polish. Confirmed. ids: D-hpc-22, D-hpc-23. Seen by 1 (docs-hpc).
- Where: Docker Hub repository page (linked from `README.Rmd:63`); `inst/config.local:5-6`.
- Problem: The Docker Hub repository has no overview (only the short description "Container for RiboPilot"), so users who follow the README link land on a blank page. `config.local` says "the project and directory is bind-mounted".
- Why it matters: Small, but the Hub page is where many users land first.
- Fix: Add a short Hub overview (bare version tags, no latest, amd64 only, and a link to the docs), and delete "and".

Also on this page: [DOC-07](#doc-07) (upgrades and the tag policy), [DOC-06](#doc-06) (network needs), [DOC-29](#doc-29) (the rawDir bind in generated Singularity profiles).

### In-app help (tooltips, window text, and the export bundle's README.txt)

<a id="doc-87"></a>**DOC-87. The Curation options "Min anchor identity" and "Flag contaminants" controls do nothing.** Major (D-intermediate-25 rated minor). Confirmed. ids: D-domain-8, D-intermediate-25. Seen by 2 (docs-domain, docs-intermediate).
- Where: Curation options window, `R/app_opts.R:297-309`; `R/params_rrna.R:84-110`; `R/curate_rrna.R:285-303`.
- Problem: The window shows and saves both controls (default 85 percent, and a checkbox), but the engine reads only the ruleset; `rrna_params()` always supplies per-gene floors and `check_contaminant = TRUE`. A comment in `params_rrna.R` admits the column "is not read by the engine".
- Why it matters: A user who tightens identity or turns off the contaminant check sees no effect and may believe the setting applied.
- Fix: Wire both through, or remove them and show the real per-gene floors in the window's help text.
- See also: [CODE_REVIEW.md](CODE_REVIEW.md) APP-M13

<a id="doc-88"></a>**DOC-88. Several tooltips and option-window texts describe behavior the app no longer has.** Minor. Confirmed. ids: D-appclaims-15, D-appclaims-21. Seen by 1 (docs-appclaims).
- Where: `R/constants.R:123-125, 160-162, 347-348, 362-372`; `R/app_reference.R:205-207, 318-336`; `R/app_modals.R:629`; `R/app_server.R:1017, 1196-1210`; `R/app_opts.R:42-44, 589-592`; `R/app_ui_helpers.R:159`; `R/tax_search.R:103-104`.
- Problem:
  - Selected Path says "Assembly path chosen for annotation", but every kept scaffold of every path is annotated.
  - "# Kept" and the Keep column help speak of annotation, but dropping only leaves a scaffold out of Annotation Review and export.
  - The Warnings tip says "sample-level" on unit rows and lists only some flags.
  - The Reference skip button, "Leave these samples on the built-in fish database", is wrong for a group already wired to a clade database. The ladder help says "offers to split", but the split is automatic.
  - The Export markers help omits ITS, and the "Written to ... as:" line omits `_itsx_check.csv` and `itsx_check/`.
  - The Annotate and Curate options windows say "Saving re-queues the selected samples: their state becomes Ready to run", and the toast says the same; for these option types the sample becomes In progress while its annotation units are queued.
  - The run window title adds samples and units together; "Workflow run finished." shows even when every task failed; the Taxonomy search help says "then apply" next to an "Assign" button; the Taxonomy Lock tip talks about options and updates that Taxonomy does not have.
- Why it matters: Tooltips are where users look first, and these contradict the app.
- Fix: Update each string to current behavior, and rename the skip action "Keep these samples on their current database".
- See also: [UX-42](USABILITY_REVIEW.md#ux-42); [CODE_REVIEW.md](CODE_REVIEW.md) APP-M2

<a id="doc-89"></a>**DOC-89. Help opens the GitHub repo, the Samples view says "groups", and the Triage window uses jargon.** Minor. Confirmed. ids: D-editor-29, D-appclaims-25 (part). Seen by 2 (docs-appclaims, docs-editor).
- Where: `R/app_server.R:85, 91`; `R/app_modals.R:276-279`.
- Problem: The header Help link opens https://github.com/dmacguigan/RiboPilot (private for now, and code-facing), not the docs site. The Reference Samples view says "Showing 12 of 12 groups". The Triage window explains depth with "k-mer collapse" and "isomer paths".
- Why it matters: The one help entry point in the app is a dead end.
- Fix: Point Help at https://dmacguigan.github.io/RiboPilot/ (ideally at the article for the current step), say "samples" in the Samples view, and simplify the Triage lead text.
- See also: [UX-52](USABILITY_REVIEW.md#ux-52), [UX-78](USABILITY_REVIEW.md#ux-78); [CODE_REVIEW.md](CODE_REVIEW.md) APP-L8

<a id="doc-90"></a>**DOC-90. Locking a unit in Annotation Review does not freeze it, and the one-path rule depends on which button you press.** Minor. Confirmed. ids: D-appclaims-16. Seen by 1 (docs-appclaims).
- Where: `vignettes/RiboPilot.Rmd:297-300`; `R/app_server.R:1925, 2097-2279, 2335-2353`; `R/app_modals.R:200-222`.
- Problem: On a locked unit, Edit boundaries, Trim assembly, Restore assembly, and Notes still work; only the Mark buttons are disabled, so a signed-off or already exported record can change. The toolbar Lock offers "Lock anyway" for a second path in a sample, while Lock & Close refuses. Neither is documented.
- Why it matters: "Locking a sample signs off on it" (the guide) is not true.
- Fix: Disable edit, trim, and restore on locked units; use one rule for both lock buttons; document the one-path-per-sample rule.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).
- See also: [UX-20](USABILITY_REVIEW.md#ux-20), [UX-25](USABILITY_REVIEW.md#ux-25); [CODE_REVIEW.md](CODE_REVIEW.md) APP-M7, APP-M15

<a id="doc-91"></a>**DOC-91. The run windows show the NCBI API key in plain text.** Minor. Confirmed. ids: D-appclaims-22. Seen by 1 (docs-appclaims).
- Where: Reference build window and Update run window; `R/nf_command.R:58`; `vignettes/Reference-Databases.Rmd` (the build options).
- Problem: The build window pre-fills the key from the environment and prints it in plain text, both in the key field and in the "Nextflow command" box; the Update command adds `--ncbi_api_key` the same way. The Reference guide does not say the field is pre-filled.
- Why it matters: A screenshot shared when asking for help leaks the key; two operator screenshots were kept out of this review's docs set for that reason.
- Fix: Mask the key in the displayed command and field, and say in the Reference guide that it comes from `NCBI_API_KEY`.
- See also: [UX-47](USABILITY_REVIEW.md#ux-47)

<a id="doc-92"></a>**DOC-92. Jargon and developer-facing text in app messages.** Polish. Confirmed. ids: D-appclaims-23, D-editor-43. Seen by 2 (docs-appclaims, docs-editor).
- Where: `R/tax_search.R:103-124`; `R/constants.R:336`; `R/triage_rules.R:31`; `R/app_modals.R:276-280, 465`; `R/run_assembly.R:71, 79`; `R/utils_reactable.R:394-405, 545-555`; `R/help_utils.R:32-33`.
- Problem: "taxdump", "core_nt", "k-mer collapse ... IGS subrepeats", "WF1" and "DB I/O", raw codes as badges (on_target, off_target, no_reference, fuzzy_species, fuzzy_genus, unresolved), and an unexplained "Resume previous run?". When a tool's bundled help is missing, the "?" popover says "Regenerate with tools/capture_tool_help.sh against the RiboPilot image", a developer script that is not installed (the bundled help itself matches the image's tool versions).
- Why it matters: Biologist users meet terms no doc defines.
- Fix: Use plain words ("local copy of the NCBI taxonomy", "GenBank (NCBI core_nt)"), show codes as readable labels, explain Resume in one line, word the tool-help fallback for users ("see the tool's website"), and put the rest in the glossary.
- See also: [UX-76](USABILITY_REVIEW.md#ux-76), [UX-80](USABILITY_REVIEW.md#ux-80)

<a id="doc-93"></a>**DOC-93. Small naming and wording mismatches in app and guide text.** Polish. Confirmed. ids: D-domain-25. Seen by 1 (docs-domain).
- Where: `vignettes/Reference-Databases.Rmd:191`; `vignettes/RiboPilot.Rmd:35, 360`; `R/app_modals.R:629-630`; `R/constants.R:153-155`.
- Problem: (a) "its quality checks anchor on the 18S and 28S genes", but 5.8S is also checked (missing_5.8S, length, identity floor 93). (b) The markers folder is `58S`. (c) The export help calls ITS1 and ITS2 "rRNA elements" and omits the combined ITS file. (d) The odd_one_out tooltip calls it "the usual signature of contamination" while the guide says it is not proof. (e) The Portal rewrites definition lines in small/large subunit wording (PZ692210 shows this), so RiboPilot's 18S/28S wording will not survive the Portal route. (f) Get started says every cell "copies" its ribosomal RNA genes; the README already says "carries".
- Why it matters: Small slips undercut trust in the rest.
- Fix: Fix each phrase, and add one sentence on the Portal's SSU/LSU definition style.
- Status: one phrase fixed in Reference databases on `review/2026-09-full-review` (docs rewrite): the checks look at the 18S, 5.8S, and 28S genes. The rest is still open.

Also on this page: [DOC-41](#doc-41) (Problematic tip), [DOC-42](#doc-42) (Partial tip), [DOC-60](#doc-60) (the screen_species column help), [DOC-56](#doc-56) (chimera panel text), [DOC-57](#doc-57) (chimera reason text), [DOC-62](#doc-62) (odd_one_out tip), [DOC-46](#doc-46), [DOC-47](#doc-47), and [DOC-48](#doc-48) (README.txt), [DOC-15](#doc-15) (names).

### Site configuration, metadata, and repo files

<a id="doc-94"></a>**DOC-94. Code in the alert boxes renders with curly quotes and breaks copy and paste, and the alert CSS is copied into each vignette.** Minor (lowered from major by the verifier: the runnable code blocks keep straight quotes; only inline code in the notes is affected). Confirmed. ids: D-editor-8, D-editor-41. Seen by 1 (docs-editor).
- Where: `vignettes/RiboPilot.Rmd:14-19, 65-83, 105-112, 235-243, 401-407`; `vignettes/Reference-Databases.Rmd:14-16, 116-122, 154-160`; `vignettes/Curation-Flags.Rmd:14-17`.
- Problem: Code inside the raw-HTML `div.alert` blocks goes through pandoc's smart quotes, and 9 code spans on the built site contain curly quotes. Pasting `conda_env = "ribopilot"` gives an R syntax error, and the pasted conda line turns `>` into a shell redirect that writes a file named `=2.12` plus a quote mark. Each vignette also carries its own `<style>` block for the alerts: Curation-Flags defines alert classes but has no alert boxes, RiboPilot.Rmd defines `.alert-tip` and `.alert-danger` without using them, and the custom `.alert` rule overrides Bootstrap's own class with hard-coded light colors.
- Why it matters: Copy and paste of documented code fails in confusing ways.
- Fix: Use pandoc fenced divs (`::: {.alert .alert-info}`) with backtick code, which fixes both problems, or turn off smart quotes for the vignettes. Put any custom alert CSS once in `pkgdown/extra.css`.
- Status: fixed on `review/2026-09-full-review` (docs rewrite).

<a id="doc-95"></a>**DOC-95. check_pkgdown fails on the DESCRIPTION URL, the favicons are stale, and the theme upper-cases every heading.** Minor. Confirmed. ids: D-editor-32, D-api-23, D-intermediate-30 (part). Seen by 3 (docs-api, docs-editor, docs-intermediate).
- Where: `DESCRIPTION:12`; `_pkgdown.yml`; `pkgdown/favicon`.
- Problem: The site report says "In DESCRIPTION, URL is missing package url (https://dmacguigan.github.io/RiboPilot)" and "Package logo is newer than favicons". Lux upper-cases all headings, including long step titles; the Articles menu is not upper-cased while the other navbar items are; the home sidebar has no citation or license links.
- Why it matters: The build is not clean, and long upper-case headings are hard to read.
- Fix: `URL: https://github.com/dmacguigan/RiboPilot, https://dmacguigan.github.io/RiboPilot/`; rerun `build_favicons()`; override the heading `text-transform` or pick another theme; add citation and license to the sidebar.
- Status: fixed on `review/2026-09-full-review` (fourth batch, reproducibility: site URL in DESCRIPTION). See [Fourth batch](REVIEW_SUMMARY.md#fourth-batch-working-at-scale-reproducibility-and-docs).

<a id="doc-96"></a>**DOC-96. The deploy workflow only builds when a commit message contains "BUILD SITE".** Minor. Confirmed. ids: D-editor-30. Seen by 1 (docs-editor).
- Where: `.github/workflows/pkgdown.yaml:18`.
- Problem: `if: contains(github.event.head_commit.message, 'BUILD SITE')`: release, workflow_dispatch, and pull_request events have no head_commit, so they always skip. The rule is not written down. GitHub Pages from a private repo needs a paid plan and would be public.
- Why it matters: The site may silently never deploy.
- Fix: Document the trigger, or gate on workflow_dispatch plus release.

<a id="doc-97"></a>**DOC-97. A stale Hydra driver script and internal notes sit at the repo root.** Minor. Confirmed. ids: D-editor-31. Seen by 1 (docs-editor).
- Where: `RiboPilot_test_hydra.R:37-67`; repo root (`BUG_*.MD`, `HANDOFF_*.MD`, `update-version.log`).
- Problem: The driver's seed path no longer ships (`system.file()` returns ""), it references `run_annotate()` (which does not exist), makes a manual `<<EXTRA_MOUNTS>>` edit that is now automatic, and runs `run_assembly()` instead of the app. The internal notes are tracked at the root and would show on the GitHub landing page.
- Why it matters: Visitors see internal files, and the driver teaches the wrong route.
- Fix: Update the driver to stop at `new_project()` and open the app; move the internal notes under `dev/`.

<a id="doc-98"></a>**DOC-98. Article group names, titles, and link text do not match, and Get started is too long.** Minor. Confirmed. ids: D-editor-39, D-editor-33. Seen by 1 (docs-editor).
- Where: `_pkgdown.yml` (articles); `README.Rmd:54, 56`; `vignettes/RiboPilot.Rmd:166, 337, 396-399`.
- Problem: The two guides sit under a navbar group called "Usage", which says nothing. Their titles ("Using RiboPilot with your group of organisms", "Understanding the automatic quality checks") do not match the link text other pages use ("reference database guide", "curation guide"). Get started is the longest page (407 source lines, Steps 0 to 9 plus a summary) and mixes install, practice, reference building, setup, console runs, review, and export.
- Why it matters: Readers cannot find the linked page in the navbar, and the main walkthrough is hard to follow.
- Fix: Rename the group (for example "Guides"), use short titles that match the cross-links, and split Get started as in the proposed outline.

<a id="doc-99"></a>**DOC-99. The DESCRIPTION Title and Description undersell the package and never mention containers.** Polish. Confirmed. ids: D-editor-37, D-hpc-24, D-api-27 (part). Seen by 3 (docs-api, docs-editor, docs-hpc).
- Where: `DESCRIPTION:2, 7-10`.
- Problem: The Title is sentence case, and the Description ("An R package for ... uses Nextflow for execution and includes an R-Shiny app") leaves out taxonomy checks, GenBank export, and the supported way to run it, in a Docker or Singularity container.
- Why it matters: This text is the one-line pitch on the site and in package listings.
- Fix: Put the Title in Title Case, and add one clause on taxonomy checks, GenBank export, and running in a Docker or Singularity container.

Also on this page: [DOC-66](#doc-66) and [DOC-72](#doc-72) (reference index and Rd markdown), [DOC-12](#doc-12) (NEWS and citation), [DOC-19](#doc-19) (site and repo links 404 while private).

## Screenshots for the docs

These 15 operator screenshots were copied into [screenshots/](screenshots/) as `docs_*.png` (all 1600 px wide and under 250 KB each, so none needed downscaling). They show the app from a clean practice run and a gastropod project, and cover every step of the proposed walkthrough.

| File | Taken by | What it shows | Where it belongs | Notes |
|---|---|---|---|---|
| [docs_01](screenshots/docs_01_taxonomy_assigned_and_locked.png) | op-practice (screen 05) | Taxonomy step after Assign Taxonomy: all 12 practice samples resolved as "exact" and locked. The toolbar shows the five numbered steps. | Get started, Step 1 Taxonomy; also the home page ("How it works"). | |
| [docs_02](screenshots/docs_02_taxonomy_confidence_badges.png) | op-edge (screen 11) | Assign Taxonomy on common Taxon mistakes: exact, fuzzy_species, ambiguous, and unresolved (a typo, a common name, a blank, and a numeric taxid). | Your own project (the Taxon column) and Troubleshooting; the Confidence explanation in Step 1. | Also shows "Morus" resolved to the plant genus (UX-08), a good "check this" example. |
| [docs_03](screenshots/docs_03_reference_propose_groups.png) | op-novel (screen 07) | The Propose groups review window: two gastropods proposed into Littorinimorpha and Neogastropoda, with Apply. | Reference databases, "The Reference tab"; Get started, Step 2. | |
| [docs_04](screenshots/docs_04_reference_setup_choices.png) | op-practice (screen 19) | The Set up database window for three fish groups: download from GenBank, use an existing database, or leave on the fish database. | Reference databases, "Setting a group up". | |
| [docs_05](screenshots/docs_05_reference_groups_ready.png) | op-novel (screen 14) | The Groups view after a build: both groups Ready, with record counts, the Samples/Groups toggle, Propose groups, and Settings. | Reference databases (statuses and the Samples/Groups toggle). | The Database column shows a review folder path; retake with a neutral path. |
| [docs_06](screenshots/docs_06_run_window_progress.png) | op-practice (screen 10) | The Update run window: the Resume box, the Nextflow command, live progress, and Stop / Interrupt. | Get started, Step 3 (run from the app); Troubleshooting (Resume). | Shows long local paths and a "Nextflow 26.04.6 is available" line. No API key is visible. |
| [docs_07](screenshots/docs_07_assemble_table_after_run.png) | op-practice (screen 28) | The Assemble table after a run: kept scaffolds, ITSx warnings, scores, Notes, and a failed sample with a red triangle. | Get started, Step 3 (reading the table); Troubleshooting (a failed sample). | |
| [docs_08](screenshots/docs_08_assemble_genbank_evidence.png) | op-practice (screen 30_1) | The Assemble details window for CHIM01: best match anywhere in core_nt (Penicillium) and the nuclear rDNA hit table. | Get started, Step 3 (the details window); Reviewing results (the GenBank screen). | |
| [docs_09](screenshots/docs_09_assemble_chimera_check.png) | op-practice (screen 30_3) | The chimera check for the planted chimera CHIM01: lineage segments along the scaffold and the suspected junction at 2,625 bp. | Reviewing results, the new "Chimeras and submission risk" section ([DOC-56](#doc-56)). | A true positive, which makes it the right example. |
| [docs_10](screenshots/docs_10_rescue_baited_assembler.png) | op-practice (screen 35) | Assembly Options for a failed sample: a new set name and "SPAdes, baited to the seed DB". | Get started, rescuing a failed sample in the app ([DOC-37](#doc-37)). | The seed path shows a local repo path; retake. |
| [docs_11](screenshots/docs_11_annotation_review_table.png) | op-practice (screen 46) | The Annotation Review table: one row per unit, with Hit Source, Warnings, ITSx Warnings, Screen, Chimera, and the Reviewed, Problematic, and Partial buttons. | Get started, Step 4; Reviewing results (the Warnings count). | SRR33849683 shows "suspect", the false positive in [DOC-57](#doc-57); retake after that fix. |
| [docs_12](screenshots/docs_12_annotation_review_detail.png) | op-practice (screen 49a) | Annotation Review details for a clean sample: all five features, coverage with feature bands, Trim and Restore, the Mark buttons, the boundary editor, and Lock & Close. | Get started, Step 4. | Also shows the 1 bp gap (18S ends at 2041, ITS1 starts at 2043) from [DOC-39](#doc-39); retake after that fix. |
| [docs_13](screenshots/docs_13_curation_warnings_panel.png) | op-practice (screen 48) | The Curation warnings panel for CHIM01: sample-level warnings, plus per-feature warnings with select links. | Reviewing results, "What the flags mean" ([DOC-53](#doc-53)). | The Reference alignment panel above it gives the "run Update again" advice that UX-28 calls wrong; retake after that fix. |
| [docs_14](screenshots/docs_14_export_window.png) | op-practice (screen 57) | The Export data window: export group, header template, per-marker files, HTML report, and the final ITSx check. | Get started, Step 5; Preparing a GenBank submission (templates). | A local output path shows at the bottom. |
| [docs_15](screenshots/docs_15_export_split_by_assembler.png) | op-practice (screen 55) | The Mixed assemblers window after Assign Group, offering to split by assembler, with Screen and Chimera for each unit. | Get started, Step 5 (export groups). | Also shows the "suspect" flounder rows. |

Other good candidates that were not copied (to stay within 15): op-practice screen 59 (Export complete, with the ITSx and submission-risk summary), 49_SRR14433899_scroll900 (reference alignment), 62_report_align_ITS1 (HTML report alignments), and 09_run_modal; op-novel screen 37 (a contaminated sample flagged off_target) and 42 (a fragmented rescued unit); op-edge screen 04 (a config error shown in the run window, for Troubleshooting) and 26 (the existing-database picker). The docs-editor reviewer also has screenshots of the current site at 1400 px and 390 px in `docs-editor/shots/`, useful as a "before" set.

Two cautions for the final set:

- **Do not publish op-novel's run or build windows** (for example screens 13 and 16). They show a real NCBI API key in plain text, in both the Nextflow command and the key field. They were deliberately not copied. The app showing the key at all is [DOC-91](#doc-91) here and UX-47 in the usability review.
- **Retake the published images from a fresh practice project at a neutral path** (for example `~/ribopilot_test`) after the fixes land. The current images show review folder paths, the chimera false positive, and the 1 bp ITS gap.

## What works well and should be kept

- **The Reference guide matches the app closely**: the labels, the status names, the four setup choices, and "records unknown". op-novel built two order-level databases from the Reference tab in about 1.5 minutes with no errors, the high point of that run.
- **The writing is plain and friendly.** The README intro and the "you do not need to be a programmer" framing work, and the Curation guide is honest about the limits of its checks (the phylum-level screen, "Depth is shown, not decided", "Your own decisions stick").
- **The in-app help is mostly good.** The column tooltips (`R/constants.R:119-159`) and the option-window help are clear and specific, and the Mixed assemblers dialog explains itself.
- **The export bundle is well described**: every listed file exists, the ITSx check only warns, and the bundle README explains how the NCBI Portal handles rRNA-ITS records.
- **Version and tag claims check out**: R >= 4.4.0, Nextflow 24.10.x through 25.10.x with the `NXF_VER` pin, and image tags that equal the package version with no `latest`.
- **The build is healthy.** pkgdown builds with no rendering warnings, `man/` and `NAMESPACE` match roxygen output, `codoc`, `undoc`, `checkDocFiles`, and the Rd cross-reference checks are clean, all internal site links resolve, and all external links work except those that 404 while the repo is private.
- **The app's evidence views** (GenBank evidence, coverage with feature bands, the chimera plot, reference alignments, and the HTML report) are clear and make excellent documentation screenshots.
- The `GITHUB_PAT` hint in the `new_test_project()` error is helpful, once it comes first.

## Appendix A: refuted findings

- **Re-annotating a unit is undocumented; the "run Update again" hint is wrong** (D-appclaims-20). Refuted: the finding said re-annotation needs the assembly to be unlocked first, which takes the unit out of Annotation Review. In fact the Annotation Review State button can set a unit to Ready, and Update then re-annotates it without unlocking the assembly. The verifier agreed that the hint is vague and the loop is undocumented; the usability review covers the hint as UX-28.
- Part of D-domain-15 (that README.txt describes ITSx "chimeric" wrongly) was also refuted; the rest of that finding stands and is in [DOC-48](#doc-48).

## Appendix B: where each source finding went

| Source id | Item |
|---|---|
| D-api-1 | [DOC-18](#doc-18) |
| D-api-2 | [DOC-17](#doc-17) |
| D-api-3 | [DOC-67](#doc-67) |
| D-api-4 | [DOC-68](#doc-68) |
| D-api-5 | [DOC-28](#doc-28), [DOC-20](#doc-20) (part) |
| D-api-6 | [DOC-40](#doc-40) |
| D-api-7 | [DOC-69](#doc-69) |
| D-api-8 | [DOC-66](#doc-66) |
| D-api-9 | [DOC-30](#doc-30) |
| D-api-10 | [DOC-72](#doc-72) |
| D-api-11 | [DOC-66](#doc-66) |
| D-api-12 | [DOC-75](#doc-75) |
| D-api-13 | [DOC-74](#doc-74) |
| D-api-14 | [DOC-70](#doc-70) |
| D-api-15 | [DOC-47](#doc-47) |
| D-api-16 | [DOC-76](#doc-76) |
| D-api-17 | [DOC-73](#doc-73) |
| D-api-18 | [DOC-80](#doc-80) |
| D-api-19 | [DOC-77](#doc-77) |
| D-api-20 | [DOC-78](#doc-78) |
| D-api-21 | [DOC-79](#doc-79) |
| D-api-22 | [DOC-19](#doc-19) |
| D-api-23 | [DOC-95](#doc-95) |
| D-api-24 | [DOC-83](#doc-83), [DOC-84](#doc-84) (part) |
| D-api-25 | [DOC-21](#doc-21) |
| D-api-26 | [DOC-35](#doc-35) |
| D-api-27 | [DOC-82](#doc-82), [DOC-22](#doc-22) (part), [DOC-81](#doc-81) (part), [DOC-99](#doc-99) (part) |
| D-api-28 | [DOC-66](#doc-66) |
| D-appclaims-1 | [DOC-28](#doc-28) |
| D-appclaims-2 | [DOC-20](#doc-20) |
| D-appclaims-3 | [DOC-40](#doc-40) |
| D-appclaims-4 | [DOC-41](#doc-41) |
| D-appclaims-5 | [DOC-38](#doc-38) |
| D-appclaims-6 | [DOC-36](#doc-36) |
| D-appclaims-7 | [DOC-53](#doc-53) |
| D-appclaims-8 | [DOC-37](#doc-37) |
| D-appclaims-9 | [DOC-68](#doc-68) |
| D-appclaims-10 | [DOC-54](#doc-54) |
| D-appclaims-11 | [DOC-55](#doc-55) |
| D-appclaims-12 | [DOC-47](#doc-47) |
| D-appclaims-13 | [DOC-44](#doc-44) |
| D-appclaims-14 | [DOC-60](#doc-60) |
| D-appclaims-15 | [DOC-88](#doc-88) |
| D-appclaims-16 | [DOC-90](#doc-90) |
| D-appclaims-17 | [DOC-51](#doc-51) |
| D-appclaims-18 | [DOC-35](#doc-35) |
| D-appclaims-19 | [DOC-13](#doc-13) |
| D-appclaims-20 | refuted (Appendix A) |
| D-appclaims-21 | [DOC-88](#doc-88) |
| D-appclaims-22 | [DOC-91](#doc-91) |
| D-appclaims-23 | [DOC-92](#doc-92) |
| D-appclaims-24 | [DOC-15](#doc-15) |
| D-appclaims-25 | [DOC-83](#doc-83), [DOC-89](#doc-89) (part) |
| D-domain-1 | [DOC-42](#doc-42) |
| D-domain-2 | [DOC-58](#doc-58) |
| D-domain-3 | [DOC-59](#doc-59) |
| D-domain-4 | [DOC-39](#doc-39) |
| D-domain-5 | [DOC-45](#doc-45) |
| D-domain-6 | [DOC-60](#doc-60) |
| D-domain-7 | [DOC-54](#doc-54) |
| D-domain-8 | [DOC-87](#doc-87) |
| D-domain-9 | [DOC-53](#doc-53) |
| D-domain-10 | [DOC-56](#doc-56) |
| D-domain-11 | [DOC-57](#doc-57) |
| D-domain-12 | [DOC-43](#doc-43) |
| D-domain-13 | [DOC-44](#doc-44) |
| D-domain-14 | [DOC-48](#doc-48) |
| D-domain-15 | [DOC-48](#doc-48) |
| D-domain-16 | [DOC-49](#doc-49) |
| D-domain-17 | [DOC-23](#doc-23) |
| D-domain-18 | [DOC-61](#doc-61) |
| D-domain-19 | [DOC-55](#doc-55) |
| D-domain-20 | [DOC-50](#doc-50) |
| D-domain-21 | [DOC-12](#doc-12) |
| D-domain-22 | [DOC-35](#doc-35) |
| D-domain-23 | [DOC-46](#doc-46) |
| D-domain-24 | [DOC-33](#doc-33) |
| D-domain-25 | [DOC-93](#doc-93) |
| D-editor-1 | [DOC-28](#doc-28) |
| D-editor-2 | [DOC-18](#doc-18) |
| D-editor-3 | [DOC-17](#doc-17) |
| D-editor-4 | [DOC-28](#doc-28) |
| D-editor-5 | [DOC-20](#doc-20) |
| D-editor-6 | [DOC-01](#doc-01), [DOC-03](#doc-03) (part) |
| D-editor-7 | [DOC-19](#doc-19) |
| D-editor-8 | [DOC-94](#doc-94) |
| D-editor-9 | [DOC-30](#doc-30) |
| D-editor-10 | [DOC-15](#doc-15) |
| D-editor-11 | [DOC-53](#doc-53) |
| D-editor-12 | [DOC-54](#doc-54) |
| D-editor-13 | [DOC-55](#doc-55), [DOC-56](#doc-56) (part) |
| D-editor-14 | [DOC-66](#doc-66) |
| D-editor-15 | [DOC-68](#doc-68), [DOC-67](#doc-67) (part) |
| D-editor-16 | [DOC-41](#doc-41) |
| D-editor-17 | [DOC-12](#doc-12), [DOC-08](#doc-08) (part), [DOC-09](#doc-09) (part), [DOC-13](#doc-13) (part) |
| D-editor-18 | [DOC-72](#doc-72) |
| D-editor-19 | [DOC-22](#doc-22) |
| D-editor-20 | [DOC-16](#doc-16) |
| D-editor-21 | [DOC-34](#doc-34), [DOC-32](#doc-32) (part), [DOC-35](#doc-35) (part) |
| D-editor-22 | [DOC-31](#doc-31) |
| D-editor-23 | [DOC-23](#doc-23) |
| D-editor-24 | [DOC-50](#doc-50), [DOC-51](#doc-51) (part) |
| D-editor-25 | [DOC-21](#doc-21) |
| D-editor-26 | [DOC-62](#doc-62) |
| D-editor-27 | [DOC-60](#doc-60) |
| D-editor-28 | [DOC-79](#doc-79) |
| D-editor-29 | [DOC-89](#doc-89) |
| D-editor-30 | [DOC-96](#doc-96) |
| D-editor-31 | [DOC-97](#doc-97) |
| D-editor-32 | [DOC-95](#doc-95) |
| D-editor-33 | [DOC-98](#doc-98) |
| D-editor-34 | [DOC-64](#doc-64) |
| D-editor-35 | [DOC-65](#doc-65) |
| D-editor-36 | [DOC-24](#doc-24) |
| D-editor-37 | [DOC-99](#doc-99) |
| D-editor-38 | [DOC-83](#doc-83) |
| D-editor-39 | [DOC-98](#doc-98) |
| D-editor-40 | [DOC-81](#doc-81) |
| D-editor-41 | [DOC-94](#doc-94) |
| D-editor-42 | [DOC-27](#doc-27) |
| D-editor-43 | [DOC-92](#doc-92) |
| D-hpc-1 | [DOC-29](#doc-29) |
| D-hpc-2 | [DOC-17](#doc-17) |
| D-hpc-3 | [DOC-19](#doc-19) |
| D-hpc-4 | [DOC-05](#doc-05) |
| D-hpc-5 | [DOC-06](#doc-06) |
| D-hpc-6 | [DOC-04](#doc-04), [DOC-28](#doc-28) (part) |
| D-hpc-7 | [DOC-07](#doc-07) |
| D-hpc-8 | [DOC-71](#doc-71) |
| D-hpc-9 | [DOC-67](#doc-67) |
| D-hpc-10 | [DOC-16](#doc-16) |
| D-hpc-11 | [DOC-01](#doc-01) |
| D-hpc-12 | [DOC-03](#doc-03) |
| D-hpc-13 | [DOC-22](#doc-22) |
| D-hpc-14 | [DOC-83](#doc-83) |
| D-hpc-15 | [DOC-24](#doc-24) |
| D-hpc-16 | [DOC-73](#doc-73), [DOC-72](#doc-72) (part) |
| D-hpc-17 | [DOC-21](#doc-21) |
| D-hpc-18 | [DOC-26](#doc-26) |
| D-hpc-19 | [DOC-04](#doc-04) |
| D-hpc-20 | [DOC-85](#doc-85) |
| D-hpc-21 | [DOC-04](#doc-04) |
| D-hpc-22 | [DOC-86](#doc-86) |
| D-hpc-23 | [DOC-86](#doc-86) |
| D-hpc-24 | [DOC-99](#doc-99) |
| D-hpc-25 | [DOC-16](#doc-16) |
| D-hpc-26 | [DOC-27](#doc-27) |
| D-intermediate-1 | [DOC-28](#doc-28) |
| D-intermediate-2 | [DOC-18](#doc-18) |
| D-intermediate-3 | [DOC-17](#doc-17) |
| D-intermediate-4 | [DOC-20](#doc-20) |
| D-intermediate-5 | [DOC-28](#doc-28) |
| D-intermediate-6 | [DOC-09](#doc-09) |
| D-intermediate-7 | [DOC-44](#doc-44), [DOC-45](#doc-45) (part) |
| D-intermediate-8 | [DOC-40](#doc-40) |
| D-intermediate-9 | [DOC-04](#doc-04) |
| D-intermediate-10 | [DOC-03](#doc-03) |
| D-intermediate-11 | [DOC-10](#doc-10) |
| D-intermediate-12 | [DOC-08](#doc-08) |
| D-intermediate-13 | [DOC-11](#doc-11) |
| D-intermediate-14 | [DOC-12](#doc-12) |
| D-intermediate-15 | [DOC-01](#doc-01) |
| D-intermediate-16 | [DOC-32](#doc-32), [DOC-37](#doc-37) (part) |
| D-intermediate-17 | [DOC-70](#doc-70) |
| D-intermediate-18 | [DOC-67](#doc-67) |
| D-intermediate-19 | [DOC-53](#doc-53), [DOC-35](#doc-35) (part), [DOC-54](#doc-54) (part), [DOC-56](#doc-56) (part) |
| D-intermediate-20 | [DOC-30](#doc-30) |
| D-intermediate-21 | [DOC-68](#doc-68) |
| D-intermediate-22 | [DOC-41](#doc-41) |
| D-intermediate-23 | [DOC-14](#doc-14) |
| D-intermediate-24 | [DOC-60](#doc-60) |
| D-intermediate-25 | [DOC-87](#doc-87) |
| D-intermediate-26 | [DOC-69](#doc-69) |
| D-intermediate-27 | [DOC-05](#doc-05) |
| D-intermediate-28 | [DOC-76](#doc-76) |
| D-intermediate-29 | [DOC-47](#doc-47) |
| D-intermediate-30 | [DOC-66](#doc-66), [DOC-72](#doc-72) (part), [DOC-95](#doc-95) (part) |
| D-intermediate-31 | [DOC-79](#doc-79), [DOC-19](#doc-19) (part) |
| D-intermediate-32 | [DOC-25](#doc-25) |
| D-intermediate-33 | [DOC-08](#doc-08) |
| D-intermediate-34 | [DOC-84](#doc-84) |
| D-intermediate-35 | [DOC-83](#doc-83) |
| D-novice-1 | [DOC-28](#doc-28) |
| D-novice-2 | [DOC-17](#doc-17) |
| D-novice-3 | [DOC-18](#doc-18) |
| D-novice-4 | [DOC-03](#doc-03) |
| D-novice-5 | [DOC-02](#doc-02) |
| D-novice-6 | [DOC-31](#doc-31) |
| D-novice-7 | [DOC-32](#doc-32) |
| D-novice-8 | [DOC-01](#doc-01) |
| D-novice-9 | [DOC-16](#doc-16) |
| D-novice-10 | [DOC-28](#doc-28), [DOC-37](#doc-37) (part) |
| D-novice-11 | [DOC-36](#doc-36) |
| D-novice-12 | [DOC-30](#doc-30) |
| D-novice-13 | [DOC-35](#doc-35), [DOC-79](#doc-79) (part) |
| D-novice-14 | [DOC-53](#doc-53) |
| D-novice-15 | [DOC-55](#doc-55), [DOC-56](#doc-56) (part) |
| D-novice-16 | [DOC-54](#doc-54) |
| D-novice-17 | [DOC-15](#doc-15) |
| D-novice-18 | [DOC-13](#doc-13) |
| D-novice-19 | [DOC-19](#doc-19) |
| D-novice-20 | [DOC-52](#doc-52) |
| D-novice-21 | [DOC-40](#doc-40) |
| D-novice-22 | [DOC-09](#doc-09) |
| D-novice-23 | [DOC-05](#doc-05) |
| D-novice-24 | [DOC-50](#doc-50) |
| D-novice-25 | [DOC-51](#doc-51) |
| D-novice-26 | [DOC-69](#doc-69) |
| D-novice-27 | [DOC-63](#doc-63) |
| D-novice-28 | [DOC-41](#doc-41) |
| D-novice-29 | [DOC-24](#doc-24) |
| D-novice-30 | [DOC-68](#doc-68), [DOC-74](#doc-74) (part) |
| D-novice-31 | [DOC-25](#doc-25) |
| D-novice-32 | [DOC-21](#doc-21) |
| D-novice-33 | [DOC-22](#doc-22) |
| D-novice-34 | [DOC-16](#doc-16) |
| D-novice-35 | [DOC-66](#doc-66) |
| D-novice-36 | [DOC-14](#doc-14) |
| D-novice-37 | [DOC-34](#doc-34) |
| D-novice-38 | [DOC-46](#doc-46) |
| D-novice-39 | [DOC-23](#doc-23) |
| D-novice-40 | [DOC-47](#doc-47) |
| D-novice-41 | [DOC-30](#doc-30) |
| D-novice-42 | [DOC-64](#doc-64) |
| D-novice-43 | [DOC-84](#doc-84) |
| D-novice-44 | [DOC-83](#doc-83) |
| D-novice-45 | [DOC-08](#doc-08) |
| O-edge-26 | [DOC-09](#doc-09) |
| O-edge-27 | [DOC-28](#doc-28), [DOC-08](#doc-08) (part), [DOC-09](#doc-09) (part) |
| O-novel-6 | [DOC-28](#doc-28) |
| O-novel-7 | [DOC-79](#doc-79), [DOC-19](#doc-19) (part) |
| O-novel-8 | [DOC-09](#doc-09) |
| O-novel-9 | [DOC-54](#doc-54), [DOC-53](#doc-53) (part) |
| O-novel-10 | [DOC-51](#doc-51) |
| O-practice-5 | [DOC-28](#doc-28) |
| O-practice-6 | [DOC-28](#doc-28) |
| O-practice-10 | [DOC-37](#doc-37) |
| O-practice-14 | [DOC-40](#doc-40) |
| O-practice-16 | [DOC-15](#doc-15) |
| O-practice-17 | [DOC-41](#doc-41) |
| O-practice-19 | [DOC-35](#doc-35) |
| O-practice-20 | [DOC-32](#doc-32) |
| O-practice-21 | [DOC-31](#doc-31) |
| O-practice-25 | [DOC-47](#doc-47) |
| O-practice-26 | [DOC-40](#doc-40) |
| O-practice-27 | [DOC-47](#doc-47) |
| O-practice-33 | [DOC-47](#doc-47) |

## Appendix C: sources and evidence

- Reviewer notes: `/media/dmacguig/Mockra/DJMCompBio/temp/RiboPilot_review/reports/docs-*.md` (novice, intermediate, hpc, domain, editor, api, appclaims) and `reports/op-*.md` (practice, novel, edge). Verified findings from the earlier novice and intermediate run: `RiboPilot_review/tmp/docs_prior.json`. The finding ids in this review (for example D-hpc-3) are the ids from the verified runs; the reviewers' Markdown notes use their own labels (for example B1 or M3).
- Evidence per reviewer (logs, scripts, repro folders, and screenshots): `RiboPilot_review/<reviewer>/`. The built site is in `docs-editor/site/`, with page screenshots in `docs-editor/shots/`. The operators' final project snapshots are in `op-practice/snapshot_final/`, `op-novel/snapshot_final/`, and `op-edge/projects/`.
- The review folder is temporary and will be deleted after the review. The screenshots used here were copied into [screenshots/](screenshots/).
- Sibling reviews in this folder: [CODE_REVIEW.md](CODE_REVIEW.md) and [USABILITY_REVIEW.md](USABILITY_REVIEW.md).
