# 08 - Multi-scaffold UX and contamination screening

Design review, 2026-07-24, branch `pipeline-review-2026-07`.

Scope: (A) UX/UI for multiple paths and scaffolds per sample, especially when
scaffolds disagree on their best BLAST hit; (B) detecting contamination and
off-target assemblies, and keeping only assemblies close to the target taxon.

Method: six parallel code surveys (assemble UX, annotate/curate UX, BLAST +
Nextflow data flow, export, MitoPilot comparison, docs-vs-code), three design
tracks, two adversarial critics (feasibility against the real code, and
biology). Everything below that carries a `file:line` was re-verified against
the working tree before being written down. Proposals that the critics killed
are recorded in section 6 so they are not re-proposed later.

---

## 0. The one-paragraph summary

The per-scaffold pivot (b0f1b12 .. 28fab91) moved the data model to
`(ID, path, scaffold)` but left the app, the export path, and the contamination
check at sample grain. The result is that RiboPilot **already computes almost
everything needed to answer "which of these scaffolds is my organism"** and then
renders none of it: `blast_screen` (the unfiltered top hit over all of core_nt,
with staxid, lineage, and `is_nuclear_rrna`) appears in **zero** app files;
`blast_ref_candidates` holds up to 50 ranked hits per unit and the app shows
one; `coverage.mean_depth` is per unit and is never compared across a sample's
scaffolds; `assemblies.depth` and `assemblies.gc` are declared and never
written. Meanwhile the only contamination test that ships,
`lineages_conflict()` (R/curate_rrna.R:188-197), is a zero-shared-token test
that **structurally cannot fire on any metazoan contaminant in a metazoan
sample**, because both lineages share the token `metazoa`. So the tool is blind
to the single commonest museum case (a co-extracted copepod, an epibiont, gut
contents) while flagging every offline project as `contaminant_unchecked`.

The fix is not more BLAST. It is (1) close a short list of real data-loss bugs,
(2) render the evidence that already exists at the moment the keep/drop decision
is made, and (3) replace the token test with a rank-aware verdict whose
governing principle is inverted from the obvious one: **never flag on low
identity, flag on implausibly high identity at large taxonomic distance.**

---

## 1. Bugs and dead code, in priority order

These are defects, not features. Several destroy user work. Fix before building
anything on top.

### 1.1 A single misclick deletes locked, reviewed curation with no confirmation

`R/app_server.R:868-882` flips `assemblies.ignore` and immediately calls
`reconcile_ignore_unit()`, which `DELETE`s the unit's rows from `annotate`,
`annotations` and `coverage` (R/pipeline_steps.R:472-476). There is no
confirmation dialog and no check of `annotate_lock`, `reviewed`, `problematic`,
`partial`, `export_group` or `annotate_notes`. The trigger is a one-word text
link in a 130px column (`rt_link("scaf_ignore")`, R/app_server.R:863). The only
feedback is a green success toast. The Annotate-mode lock guard
(R/app_server.R:473-493) cannot see this path, because the deletion happens from
the Assemble modal.

Fix: split into a doer plus a guarded observer. Hard-refuse when
`annotate_lock = 1`; confirm when `reviewed = 'yes'` or `export_group` is set.
Better, make dropping non-destructive (section 4.1), after which this reduces to
the lock check alone.

This also matters biologically, not just as a data-loss bug. "Drop it, look
again, put it back" is the *normal* motion when comparing rDNA scaffolds,
because whether a divergent scaffold is a contaminant or a real intragenomic
array copy is genuinely ambiguous and is resolved by toggling. A destructive
toggle makes the correct workflow impossible.

### 1.2 Re-assembly silently destroys every keep/drop decision

`load_assemblies()` does `DELETE FROM assemblies WHERE ID = ?`
(R/pipeline_steps.R:275) and then re-appends rows built with a hardcoded
`ignore = 0L` (R/pipeline_steps.R:240). Any sample re-queued to
`assemble_switch IN (1,4)` loses its ignore flags, and `opts_apply` re-queues
automatically whenever a curator changes `assemble_opts` (R/app_opts.R:349-350).

This is the highest-priority item in the whole review, because **any** triage
design built on `assemblies` evaporates on the first re-run without it. Fix:
preserve-and-remerge keyed on `(path, scaffold)` inside `load_assemblies()`.

### 1.3 Export ignores the per-unit selection and emits unstable SeqIDs

The Export tab renders one row per `(ID, path, scaffold)` and Assign/Clear Group
correctly use `selected_units()` (R/app_server.R:529, :538). Export does not:
`export_ids()` (R/app_server.R:759-769) collapses to bare sample IDs, and
`export_bundle`'s filter is a sample-level membership test,
`units[units$ID %in% ids, ]` (R/export.R:357). **Selecting scaffold 2 of a
four-scaffold sample exports all four.**

Worse, `n_per_id <- table(units$ID)` (R/export.R:362) is computed *after* both
the group and ids filters, so the same physical sequence exports as
`USNM123458` in one bundle and `USNM123458_p1_s1` in another. Once an accession
is bound to a SeqID, silently renaming records between runs is unrecoverable.

Fix: add a `units` argument to `export_bundle()`; compute the SeqID denominator
from an unfiltered query over the *correct* exportable set
(`annotate_lock = 1 AND COALESCE(asm.ignore,0) = 0 AND COALESCE(an.problematic,'') != 'yes'`,
matching R/export.R:353-355); add the same two predicates to the Export tab
query (R/app_server.R:123-131) so a visible, selectable row cannot be silently
absent from the bundle. Grey excluded rows rather than hiding them, so the
curator can see what was left out and why. Record the SeqID scheme in the
bundle: `_p1_s2` encodes assembly-graph coordinates that mean nothing to GenBank
and will otherwise be read as biological replicates.

### 1.4 "Set as selected path" is a no-op with a harmful side effect

`R/app_server.R:883-896` writes `path_scores.auto_select`. WF2's unit query
never reads `path_scores` (inst/nextflow/modules/pipeline_r.nf:60-72 filters on
locks plus `COALESCE(asm.ignore,0) = 0` only), and `load_assemblies()` seeds an
annotate unit for every non-ignored scaffold of every kept path
(R/pipeline_steps.R:297-305). So the modal's only primary-styled button changes
no pipeline behaviour, and the success toast asserts behaviour that does not
happen.

Its one real effect is a bug:
`UPDATE annotate SET annotate_switch = 1 WHERE ID = ?` (R/app_server.R:892)
re-queues **every** unit of the sample, including annotate-locked ones.

Fix now: delete line 892. It is an unambiguous bug and needs no replacement.
Defer the replacement action ("keep only this path") until drops are
non-destructive.

### 1.5 The Assemble table's BLAST columns are permanently blank

`R/app_server.R:117` selects `a.blast_species, a.blast_accession, a.blast_pident,
a.blast_qcovs` from the per-sample `assemble` table. Verified: **no statement
anywhere writes those columns.** The only `UPDATE assemble SET` sites are
R/pipeline_steps.R:222, :290 and R/run_fallback.R:58, none of which touch
`blast_*`. The per-unit truth lives in `assemblies`
(written at R/pipeline_steps.R:68 from the local seed DB, then possibly
overwritten at R/blast_genbank.R:97 from the remote search). So the entire
"BLAST" entry of the Columns picker (R/app_server.R:253) toggles four blank
columns, at the exact moment the reviewer is asked to judge a sample.

Two adjacent defects in the same query: `a.scaffolds` is
`max(table(df$path))` (R/pipeline_steps.R:280), the largest *single path's*
scaffold count rather than the total unit count, and `a.length` is a
`;`-joined `unique()` list (R/pipeline_steps.R:281) rendered in a 90px
`wrap = FALSE` column. Neither is updated when a scaffold is ignored.

Fix: stop selecting the dead columns and add `units` / `kept` counts from
`assemblies`. **Do not drop them from the DDL** - R/migrate_db.R:41-42 documents
"New columns are only ever added, never dropped or retyped", and dropping them
in `init_db.R` only would make fresh projects diverge from migrated ones. Fix
the `assemble.scaffolds` writer at source (R/pipeline_steps.R:280) as well, or
the CSV download (R/app_server.R:82-90) keeps emitting the wrong number.

### 1.6 Reopening the same unit's Annotate modal wipes `annotate_notes`

`sel_unit` is a `reactiveVal` and `row_unit()` (R/app_server.R:15-20) rebuilds a
structurally identical list, so a *second* open of the same unit assigns an
`identical()` value, the reactiveVal does not invalidate, and the reset observer
(R/app_server.R:906-918) never reaches its `updateTextAreaInput()` at :916. But
`annotate_modal()` renders a fresh `textAreaInput` with no `value` argument
(R/app_modals.R:147-149), so the client reports `""`, and the unguarded write
observer (R/app_server.R:1236-1241) persists it. Since `curate_one()` writes the
machine flag summary into the same column (R/curate_rrna.R:305), a second modal
open destroys either the curator's note or the flag summary.

Fix: seed the textarea from the DB in `annotate_modal()`, plus a guard in the
write observer that never blanks a non-empty stored note. The deeper collision
(one column serving as both machine summary and human notepad) needs a schema
split, section 3.7.

### 1.7 `n_flags` omits every feature-level flag

`curate_rrna_core()` maintains two disjoint namespaces. `length`
(R/curate_rrna.R:88), `low_identity` (:103) and `truncated` (:115) call
`add_feat_flag()` only and never touch `sample_flags`, while `n_flags` is
`length(sample_flags)` (:159). So a unit whose 18S is short **and** below its
identity floor **and** truncated shows `n_flags = 0`, and
`output$curate_notes` (R/app_server.R:1053-1062) renders a green
"No curation flags for this sample" directly above `flag_tbl`
(R/app_server.R:1063-1080) listing all three.

**Blocked on a prerequisite.** Folding feature flags into the badge without
first fixing `missing_<element>` makes the badge worse, not better.
`curate_rrna_core()` emits `missing_<element>` for every one of
18S/ITS1/5.8S/ITS2/28S absent from the unit, evaluated **per scaffold**
(R/curate_rrna.R:50-51). Post-pivot, a legitimately 28S-only fragment - an
ordinary product of a tandem-repeat assembly, not a defect - already collects
four `missing_*` flags. Adding feature flags on top puts the typical *clean*
unit at 4-7 flags, at which point the badge stops discriminating and the curator
learns to ignore it.

Fix, in this order: (a) evaluate repeat completeness **per sample across kept
units** ("does this sample have all five elements between its scaffolds?") and
emit at most one per-sample incompleteness flag, demoting per-scaffold
`missing_*` to informational and out of `n_flags`; (b) then fold feature flags
in; (c) suppress `low_identity` and `truncated` on features whose `edited` flag
is 1, since `save_edits`/`trim`/`restore` re-run `curate_one()` against
`blast_ref_alignment` values computed from pre-edit coordinates.

### 1.8 Path/Scaf columns vanish exactly when the project reaches its goal state

`multi_unit` (R/app_server.R:168-174) is
`nrow(unique(d[, c("ID","path","scaffold")])) > length(unique(d$ID))`. Once the
curator prunes each sample to one unit, total units equals distinct IDs, and
Path/Scaf **disappear** - so the record of *which* isomer was locked is lost
from the table. It is also computed inside `isolate(samples_data())` within
`renderReactable` (R/app_server.R:156) while `refresh()` pushes data only
(:280-282), so it cannot change mid-session: un-ignoring a scaffold produces two
visually identical rows for one sample with different locks and flags.

Fix: MitoPilot's rule (hide only when `all(path == 1 & scaffold == 1)`) via
mounted columns plus a CSS class driven from an `output$col_css` observer, which
is what MitoPilot/R/app_annotate.R:269-293 does. Do **not** use a plain
non-isolated reactive: Shiny reactives do not compare values, so reading one
inside `renderReactable` re-renders on every `refresh()` and defeats the
`isolate()` that preserves client-side filter/sort/page. Add a sticky `unit`
column (`"p1.s3"`) after ID so a row's identity is always one glance away.

Which GetOrganelle isomer path was locked is part of the analytical record for a
submitted rDNA sequence, and a reviewer needs it to reproduce the result.

### 1.9 Smaller confirmed defects

| Defect | Location | Note |
|---|---|---|
| `curate_opts.min_anchor_pid` and `check_contaminant` written by the UI, read by nothing | R/app_opts.R:340-342 writes; R/curate_rrna.R:241-244 selects only `co.ruleset` | `check_contaminant` is hardcoded `TRUE` in base params (R/params_rrna.R:77), so the checkbox cannot turn the check off |
| The `min_anchor_pid` UI is a single scalar defaulting to 85 | R/init_db.R:548, R/app_opts.R:133 | The engine uses a per-element named vector via `anchor_pid_floor()`: metazoa 18S 95 / 5.8S 93 / 28S 80; fish 97/95/85; sparse inverts 93/90/78 (R/params_rrna.R:60-63). A scalar 85 is simultaneously far too loose for 18S and too tight for 28S |
| Explicit `metazoa_generic` treated as "unset" and silently re-derived from lineage | R/curate_rrna.R:247-249 | A user cannot force the generic ruleset. Distinguish NULL/NA (auto) from an explicit choice |
| BLAST hits ranked by `pident * qcovs` | R/blast_genbank.R:71-76 | Systematically favours **short** high-identity hits: a 400 bp conserved-18S-core HSP against an unrelated organism scores 99*95, a full 12 kb repeat from a genuine relative scores 95*95. Rank by bitscore, which scales with alignment length |
| `assemblies.blast_*` is mixed-provenance with no record of which | R/pipeline_steps.R:68 (seed DB) then R/blast_genbank.R:97 (core_nt) | `blast_pident` silently means two different things. Add `blast_source TEXT` |
| `blast.nf` exits 1 for the whole per-sample task on any unit's failure, under a shared 3-retry budget | inst/nextflow/modules/blast.nf:18-19, :56 | A 30-unit sample re-submits already-completed searches to NCBI. Cheap fix: skip units that already have hit rows |
| ~~`assemble_opts.max_scaffolds` stored and editable, never enforced~~ FIXED | R/pipeline_steps.R:292-320 | Both `max_paths` and `max_scaffolds` are now GATES, not caps: a sample over either is failed with a note, its scaffolds parked (`ignore = 1`), and the sequences kept so raising the limit and re-running recovers it without re-assembling. `max_scaffolds` counts scaffolds within one path (worst path), not the sample total. NA or <= 0 means no limit. This doc argued against a flat count (degraded museum extracts legitimately fragment one repeat across many short contigs) and proposed capping by cumulative length or anchor content; the maintainer chose the flat count anyway, reasoning that a user can raise the limit and re-run. Recorded so the objection is not re-litigated blind |
| `assemble_opts.min_assembly_length` stored and editable, never enforced | R/init_db.R:43, R/app_opts.R:119, R/app_opts.R:389 write it; nothing reads it | The remaining dead option of the pair. Either join it to the gate above or drop it from the schema and the modal |
| `score_candidates()` coerces missing sub-scores to 0 | R/scaffold_scoring.R:135-146 | A path whose BLAST never ran scores identically to one that BLASTed at 0% identity. Renormalize over available sub-scores |
| `write_blast_screen()` extracts the ~500 MB taxdump per unit with no cleanup | R/blast_genbank.R:112-119, cf. R/taxonomy.R:379-380 | Not a cumulative leak (blast.nf:58 spawns a fresh R process per unit and R clears its tempdir), but a 30-unit sample does 30 serial 500 MB untars. Cache per process |
| `compute_feature_alignment()` reads `stitched` and queries `annotations` with no path/scaffold predicate | R/blast_ref_align.R:11 | Under the per-scaffold model this mixes features from different units of one sample into one alignment, producing wrong per-element identities and therefore wrong `low_identity` flags. Zero callers, absent from NAMESPACE - delete |
| `main_fallback.nf` broken post-pivot | inst/nextflow/main_fallback.nf:43-44 | `stitch.out` is `(id, fasta)`, joined with `(id, reads)` gives a 3-tuple fed to `coverage`'s 2-input. `run_fallback()` is silently broken |
| `rt_highlight_row()` returns `#FFFFFF` for unselected rows | R/utils_reactable.R:5-12 | Silently cancels the `striped = TRUE` at R/app_server.R:264, so there is not even zebra banding to separate one sample's units from the next |
| `selected_ids()` has no `unique()` | R/app_server.R:367-375 | Selecting three scaffolds of one sample reports "3 selected" |
| `opts_apply` writes are hard-wired to `WHERE ID = ?` | R/app_opts.R:347-351 | Changing one unit's curation ruleset rewrites and re-queues every unit of the sample |
| Docs describe a pipeline stage that no longer runs | R/run_assembly.R:78, R/run_fallback.R:70, inst/nextflow/modules/blast.nf:1 | All still say "assemble -> score -> stitch -> BLAST"; main.nf:4 says "Assemblies are NOT stitched". dev/design/07:152-166 says the Level-1 stitch is built and just needs integrating - the opposite of what happened |
| `app_con()` is a bare `dbConnect` with no `busy_timeout` and no `journal_mode` | R/app_utils.R:71-74, cf. `open_project_db()` R/pipeline_steps.R:185-203 | Every new per-unit surface multiplies concurrent writers against the same SQLite file Nextflow tasks write to, on the network filesystems these projects live on |
| No indexes anywhere beyond implicit PKs | - | Every proposed per-unit join (`assemblies` + `blast_screen` + `coverage` + `annotate`) re-runs on each `samples_data()` invalidation |

Four surfaces every proposal below builds on have **zero test coverage**:
`score_candidates`, `taxonomic_match`, `parse_gfa_depth`,
`parse_blast_candidates`. `dev/test_scoring.R:44` already contains a fungal
contaminant fixture - the only labelled contamination case anywhere in the repo.
Promote it into `tests/testthat/` as the first screen regression test.

---

## 2. Quick wins on data the database already holds

Highest value per line. All use existing columns; none needs a schema change or
a pipeline change.

### 2.1 Render `blast_screen` where the decision is made

`blast_screen` holds, per `(ID, path, scaffold)`, the top hit over **all** of
core_nt with species, staxid, resolved lineage and `is_nuclear_rrna` (written at
R/blast_genbank.R:80, :126). It is the only place in the schema that can say
"this scaffold is a symbiont" or "this is not nuclear rDNA at all".

Verified: `grep -rn "blast_screen|is_nuclear_rrna|blast_lineage" R/app_*.R
R/utils_reactable.R` returns **nothing**. The evidence reaches the user only as
a substring of `annotate.annotate_notes`, inside a panel that ships collapsed
(R/app_modals.R:143), two stages after the keep/drop decision, and only if WF2
has already run. The column the modal *does* show
(`assemblies.blast_species`, R/app_server.R:848) is nuclear-rDNA-title-filtered
and falls back to the clade seed DB, which is host-clade by construction - so a
contaminant scaffold and the target scaffold render as two rows of identical
plain black text.

Put `blast_screen` next to `blast_ref_candidates` in the Assemble modal. Two
adjacent blocks - "best match anywhere in core_nt" versus "best match among
records whose title says nuclear rDNA" - are the whole contamination call, and
the database already holds both. When they name different organisms, that *is*
the signal.

**Do not render this as a green "host" badge.** The obvious implementation
drives the badge from `lineages_conflict()`, which cannot fire within a kingdom
(section 3.1). A green "host" on the co-extracted copepod scaffold converts an
absence of evidence into a positive assertion, and leaves the curator worse off
than showing nothing. Label states by what they actually assert:
`shares: Notropis (genus)` / `shares: Eukaryota` / `no reference`. Reserve red
for positive evidence only.

### 2.2 Show the ranked hit list that is already stored

`blast_ref_candidates` holds `rank, accession, species, title, pident, qcovs,
evalue` for up to `max_target_seqs` (default 50, R/init_db.R:267) hits **per
unit** (R/blast_genbank.R:87-93). The app reads it nowhere: its only consumers
are `compute_feature_msa()` (R/blast_ref_align.R:47-49) and the rank-1 pick in
`compute_ref_alignment()` (:110-113), so the user sees it only as MSA row
labels. A reviewer asking "is this really *Danio*, or is *Danio* just the best
of a bad list?" - the question that matters most for a taxon with sparse
references - has no way to look.

Add HSP length and `qcovs` as visible columns. A 99% identity over a 300 bp
conserved 18S core carries almost no taxonomic information, and without the
length beside it a biologist reads it as a strong match.

Handle the empty case explicitly. `parse_blast_candidates()` returns early when
the nuclear-rDNA title regex matched nothing (R/blast_genbank.R:84-86), which is
common in sparse invertebrate clades and is **itself a finding**. Render
"no nuclear rDNA references matched this scaffold", not a bare empty table that
reads as a pipeline failure.

### 2.3 Make the scaffold table usable and evidence-bearing

`output$assemble_scaffolds` (R/app_server.R:846-865) is the only per-scaffold
view and the only per-scaffold control in the app: `path, scaffold, length,
topology`, one BLAST hit, and a text-link toggle, at `defaultPageSize = 10` with
no `filterable`, no `searchable`, no `sortable`, no selection. With
`max_paths = 10` (R/init_db.R:41) and `max_scaffolds = 10` per path, a routine
rDNA sample is 20-40 rows, i.e. 2-4 pages of blind paging. (Both are gates now,
so an over-limit sample is failed rather than trimmed; the row count for a
sample that passes is unchanged.)

Add, all from existing tables: `coverage.mean_depth` and a within-sample
relative depth, `assemblies.blast_qcovs` (written at R/blast_genbank.R:97, never
displayed), the annotated element list from `annotations`, and the path's
`rank`/`score` from `path_scores`. Turn on `filterable`, `searchable`,
`sortable`, page size 20.

**Hard dependency, not advisory.** The `scaf_ignore` handler re-queries with
`ORDER BY path, scaffold` (R/app_server.R:870) and maps the clicked row *index*
back to a unit, and the function on the other end `DELETE`s curation. Changing
this table's ordering without first refactoring both the renderer and the
handler onto one shared reactive silently toggles the **wrong** scaffold. Do
both in the same change. Better, key clicks on a `uid` column
(`paste(ID, path, scaffold, sep = "|")`) rather than a row index - IDs are
already constrained to `^[a-zA-Z0-9_:-]+$` (R/init_db.R:66), so `|` cannot
collide.

Label the relative-depth column **"Depth (rel, within sample)"** and say so in
the tooltip. `coverage.nf` builds one bowtie2 index over *all* of a sample's
non-ignored scaffolds at once (inst/nextflow/modules/coverage.nf:22-32), and
GetOrganelle paths are near-identical rDNA isomers, so bowtie2 splits multireads
across them arbitrarily. Cross-path depth ratios are a multi-mapping artifact,
not an abundance comparison. Within a path they are fine.

### 2.4 Sibling concordance: the odd scaffold out

This is the highest signal-per-line contamination check available, and it fires
in precisely the structural blind spot the current test has. It needs no
threshold, no new column, and no network. MitoPilot already does it
(`.score_concordance`, MitoPilot/R/assembly_path_scoring.R:221-236).

The idea: "nine scaffolds say *Notropis*, one says a copepod" is catchable
without any absolute taxonomic threshold, by asking each unit whether it is the
odd one out among its siblings.

But the naive genus-level implementation is wrong for exactly the target taxa,
and both critics independently killed it. Nine scaffolds of one octocoral with
no congeneric GenBank record top-hit **nine different genera**, because each
scaffold's best hit is determined by which region it covers and which records
happen to exist for that region. Genus-level concordance therefore flags every
unit on precisely the samples the tool exists for.

Correct form: **concord at the deepest rank at which a strict majority exists,
and report that rank.** This preserves the copepod catch (majority resolves at
phylum: 9 Chordata vs 1 Arthropoda) while going silent on sparse clades (all
nine agree at class, no flag). Requirements:

- `n >= 3` units with a resolved taxon, and a **strict** majority (> 50%), not a
  tie-break. A 2-2 split has no majority.
- Rank distance from taxid-path intersection, not name tokens (section 3.2).
- Strip `PREDICTED:`, `UNVERIFIED:`, `MAG:`, `TSA:` and `Uncultured` prefixes
  before deriving any taxon token. `extract_species()` keeps the first two
  defline words verbatim (R/blast_genbank.R:138-145), so RefSeq deflines produce
  a pseudo-genus `predicted:` - and if enough units hit RefSeq records, the
  *modal* genus becomes `predicted:` and every correctly-named unit is flagged
  odd. That is the exact inverse of the intended behaviour.
- Exclude taxids under `unclassified` and `environmental samples` entirely
  (section 3.3).

---

## 3. Contamination and off-target detection: one coherent scheme

### 3.1 Why the current check cannot work

```r
# R/curate_rrna.R:188-197
LINEAGE_STOPWORDS <- c("cellular organisms", "eukaryota", "opisthokonta")
lineages_conflict <- function(a, b) {
  ...
  length(intersect(ta, tb)) == 0
}
```

A conflict requires **zero** shared tokens after dropping three stopwords and
tokens of `nchar <= 3`. A fish sample and a copepod contaminant both carry
`metazoa` and `bilateria`, so the intersection is non-empty and the copepod is
reported clean. The function's own comment concedes it only catches
wrong-kingdom hits. Every metazoan-on-metazoan contaminant - gut contents, a
parasitic copepod, an epibiont bryozoan, a human handling contaminant - is
invisible. That is the common museum case, not the exotic one.

At the same time, `curate_one()` raises `contaminant_unchecked` whenever
**either** lineage is unresolvable (R/curate_rrna.R:127-131), so an offline HPC
run - the exact deployment the run modal is built for - flags every single unit
and the signal is worthless. `vignettes/Curation-Flags.Rmd:43-45` tells the user
these checks are "simply skipped" when NCBI is unavailable, which is the
opposite of what the code does.

### 3.2 Governing principle: flag implausibly HIGH identity, never low

The tool exists because target taxa lack close references. So "93% to a distant
echinoderm" must read as **normal**, not suspect. The documented alternative
(dev/design/07_SCAFFOLD_SELECTION.md:87-89 asks for rank-depth scoring that
would flag a phylum-only match) is catastrophic in the opposite direction: it
would flag every asteroid, scyphozoan and octocoral whose nearest GenBank 18S is
in a different family, i.e. exactly the taxa in the invert test set.

Neither is right, because both read identity as a scalar. **Identity and rank
distance must be read jointly, and the informative direction is upward.**

- A **true host** sequence from an under-represented taxon shows LOW identity at
  LARGE rank distance. That is normal and must never be flagged.
- A **contaminant** shows HIGH identity at LARGE rank distance, because the
  contaminant's own species is well represented in GenBank. A fish scaffold
  matching a bird 18S at 99.8% is not a divergent fish 18S. It is the bird's.

This inversion also survives the objection that would otherwise sink it:
accelerated-rate lineages (Myxozoa, Strepsiptera, Neodermata, Acoela, Rotifera,
many Platyhelminthes and Nematoda) push 18S identity **down**, sometimes below
inter-phylum values. An upper bound is never tripped by low identity.

Low identity at large rank distance gets its own **non-punitive** verdict:
`no_reference`. It renders neutral grey and gates nothing.

**Critical granularity caveat.** The envelope is applied to BLAST `pident`,
which is computed over an HSP, and megablast HSPs preferentially land in the
**conserved core** of 18S and 28S - that is where the exact 28-mers are. So the
motivating example is true for whole-gene alignment (fish vs tetrapod 18S is
~92-94% over 1800 bp) and **false** for a 300 bp core HSP, which routinely runs
98-99.5% between any two vertebrates with no contamination present.
Consequences:

- Gate the envelope on HSP length as a **fraction of the element** (>= 60% of
  the element's expected length from `rrna_params` lengths), not an absolute bp
  count. Short core-only HSPs get **no** envelope verdict at all.
- Having done that, the envelope's discriminating power is concentrated on
  **28S** - whose mosaic of near-universal core and hypervariable D-expansion
  segments makes whole-gene identity fall much faster - and on the margin and
  LCA terms, not on 18S.
- 5.8S gets no envelope at all: 157 bp and near-invariant across metazoan
  classes, so 99% is one or two differences. Confirm-only.
- ITS: a categorical rule, but weaker than the obvious one. Require the HSP to
  be **contained within** the ITS with >= 150 bp of ITS-internal alignment, and
  trigger at `d >= 5`, not `d >= 3`. Many metazoan ITS2 are only 200-350 bp so a
  "60% local coverage" test is trivially met by an HSP anchored on the conserved
  3'-18S or 5'-28S flank spilling inward; anthozoan ITS is unusually conserved
  across genera in some octocoral groups; and slow-evolving groups do align at
  family and sometimes order.

Every envelope number must be recorded **together with the taxa it was derived
from and the alignment regime it applies to** (full-gene vs HSP). Constants
recorded without their granularity get reused at the wrong granularity, which is
exactly how the "18S does not stay at 99.8% across classes" example ends up
true in principle and false in the pipeline.

All numbers live in `rrna_params()` as a per-clade `screen` sublist, so the
sparse-GenBank rulesets (R/params_rrna.R:102-104) can widen them. The per-clade
widening must relax the ITS rule too, not only the envelope tolerance.

#### 3.2.1 Measured against the reference project, 2026-07-29

Before building the envelope, its constants were checked against the only real
identity data available: 1833 BLAST hits from the ten-sample fish project. Two
findings, both of which must be respected by any future attempt.

**The corpus cannot calibrate the envelope.** The envelope acts on high identity
at LARGE rank distance, and the corpus is almost entirely close-range: 12 hits at
species, 48 genus, 388 family, 211 order, 1172 class, **2 at phylum, and none
beyond**. There is no data in the band where the envelope would fire.

**The near side rules out a naive ceiling.** Mapping 1110 HSPs onto annotated
element spans at >= 60% element coverage, maximum identity for GENUINE host hits:

| element | family (d2) | order (d3) | class (d4) |
|---|---|---|---|
| 18S | 100.0 | 99.7 | 99.6 |
| 28S | 99.8 | 99.5 | 98.0 |
| 5.8S | 99.4 | 100.0 | 99.4 |
| ITS1 / ITS2 | 83.9 | 87.3 | 83.0 |

Genuine host 18S reaches 99.6% at class distance, so any envelope firing on high
identity below phylum distance would flag real host data, which is the failure
section 3.2 says would make the tool useless. Section 3.2 predicts the
discriminating power concentrates on 28S; measured, even 28S sits at 98.0% at
class distance. Only ITS shows usable slope, and only 16 HSPs of it were
observed.

**Constraints for any future envelope:** it must not fire below phylum distance
for 18S, 28S or 5.8S; ITS is the only element with measured slope; and the
constants must come from data in the d5-d7 band, which does not yet exist. This
is why sub-project 3b was built WITHOUT the envelope.

### 3.3 Taxonomy primitives (prerequisite for everything else)

Every taxonomic comparison in the codebase is built on a lineage string that
silently loses rank alignment. `tk_taxid_info()` reformats to
`{K};{p};{c};{o};{f};{g};{s}` and then collapses empty slots with
`gsub("; *;", ";", ...)` (R/taxonomy.R:169), so a taxon missing an order rank
has its family sitting in the order slot. `taxonomic_match()`
(R/scaffold_scoring.R:88-100) compares those positionally as an unbroken prefix,
so it under-counts shared ranks - and this hits hardest on the under-described
taxa RiboPilot targets, i.e. precisely backwards.

Needed:

1. **`{K}` is superkingdom in taxonkit, not kingdom.** The stored string for a
   fish begins `Eukaryota; Chordata; ...` and **Metazoa never appears in it**.
   A copepod, a fungus and a diatom all reduce to "shares: Eukaryota", which
   collapses the metazoan-contaminant case back into the cross-kingdom bucket.
   Use `{k}`, or emit both, and verify with `taxonkit reformat --help` before
   writing tests against it. Until then do not label the superkingdom slot
   "kingdom" anywhere in the UI.
2. **`rank_distance()` by exact set intersection over the taxid PATH**
   (`taxonkit lineage -t`), returning 0=species .. 7=no shared taxon. Taxid
   intersection is immune to synonyms, homonyms and name-string formatting,
   which token matching is not. This replaces both `taxonomic_match()` and
   `lineages_conflict()`.
3. **`-F` placeholders must map to NA.** `taxonkit reformat -F` fills missing
   ranks with `unclassified <parent>` strings; two unrelated taxa can otherwise
   "match" on an identical placeholder and `rank_distance` reports a shared rank
   that does not exist. Probe for the flag (`taxonkit reformat --help`) rather
   than parsing a version.
4. **`bitscore_weighted_lca()`** over hits within 10% of the top bitscore
   (minimum 3 when available). This is the cheap defence against a single
   mislabelled GenBank accession, which is the failure mode a one-hit
   `blast_screen` cannot survive.
5. **Exclude uncultured and environmental-sample taxids from every verdict,
   concordance and reference-support computation, while still displaying them.**
   They dominate core_nt 18S, `extract_species()` yields "Uncultured eukaryote",
   and their lineages resolve to superkingdom only - which every rank rule reads
   as maximum taxonomic distance. Mark them at the parser.
6. Store `lineage_taxids` alongside `lineage` wherever a lineage is frozen
   (`samples`, `blast_screen`).

### 3.4 Retain the hit distribution (zero extra remote round trips)

`inst/nextflow/modules/blast.nf:51-52` asks NCBI for
`-outfmt "6 saccver stitle pident qcovs evalue staxids" -max_hsps 1`: no query
coordinates, no bitscore, one HSP per subject. `parse_blast_candidates()` then
persists exactly one row of up to 50 hits into `blast_screen`
(R/blast_genbank.R:80) and discards everything the title regex rejects (:83-86).
The pipeline pays the full cost of a remote core_nt search and throws away every
signal a contamination call needs.

Change the search, not the number of searches:

```
-outfmt "6 saccver stitle pident qcovs evalue staxids bitscore length qstart qend sstart send"
-max_hsps 5     # was 1
```

~150 lines per unit instead of 50, zero extra round trips. This buys, for free:

- **bitscore**, the statistic that actually orders alignments (section 1.9).
- **the hit distribution**: is the top hit an outlier or a consensus?
- **the margin** between the best on-target and the best off-target hit.
- **positional structure**, which is exactly what chimera detection needs:
  host 18S fused to a fungal 28S produces two disjoint HSP sets with different
  subjects. Windowed *remote* BLAST is explicitly rejected - it would multiply
  the pipeline's slowest, most failure-prone step for information the HSP
  coordinates already contain.

Add a bounded `dc-megablast` fallback when megablast returns nothing or best
`qcovs < 25`. This is biologically necessary, not optional: megablast word size
28 requires a 28 bp exact match, which a divergent fungal or bacterial
contaminant frequently does not provide, so "no hits" currently reads as
`contaminant_unchecked` rather than as a signal.

Store in a new `blast_hits` table. **Keep `blast_screen` and
`blast_ref_candidates` as real tables** derived from it in the same transaction,
not as SQL views: `parse_blast_candidates()` `DELETE`s from and
`dbAppendTable()`s into both (R/blast_genbank.R:45-48, :93, :126-132), views are
not writable, and `migrate_db`'s `has_tbl()` uses `dbListTables` which returns
views too - so the v7/v8 guards (R/migrate_db.R:228-249) would see a view and
skip the rebuild.

Recompute `is_nuclear_rrna` from the **consensus of the top hits**, or better
from **HSP overlap with an annotated 18S/5.8S/28S/ITS region**, rather than from
hit #1's title. `NUCLEAR_RRNA_TITLE_RE` (R/constants.R:54-57) is a literal-phrase
regex that misses "18S rRNA gene", "SSU rRNA", and any WGS or chromosome record
containing the array - and host chromosome-level assemblies are now among the
commonest best hits, so a genuine host scaffold currently scores
`is_nuclear_rrna = 0`. Broaden the regex as a fallback, but **keep the
`ORGANELLAR_TITLE_RE` exclusion** (R/constants.R:61): it is the single
offline-independent signal that can say "this is a mitochondrial assembly, not
nuclear rDNA".

### 3.5 Evidence channels

Each returns a signed score in [-1, 1] (negative = on-target) plus an
availability flag. Confidence is renormalized over whichever channels are
available - the MitoPilot pattern RiboPilot's scorer lacks.

**(a) BLAST taxonomy** (strongest). Rank distance, the identity envelope, the
on-target/off-target bitscore margin, and the bitscore-weighted LCA.

**(b) Tool kingdom** (only channel that works fully offline). The annotation
tools already emit a taxonomic signal per scaffold and the pipeline discards it.

- `parse_barrnap()`'s name map (R/annotate.R:32) recognises only
  18S/5_8S/5.8S/28S and drops everything else, so 16S/23S are discarded
  regardless of kingdom. **That is the real blocker**, not the `-k` flag:
  inst/nextflow/modules/pipeline_r.nf:137-141 runs `pybarrnap -k "$KINGDOM"`
  from `annotate_opts.rrna_kingdom`, which already accepts bac/arc/euk/all
  (R/init_db.R:523-525) and is already exposed in the options UI. Extend the
  parse map behind a flag; run a separate `-k all` screen pass so WF2 keeps its
  euk feature model.
- **Two traps.** Metazoan mt-SSU and mt-LSU are bacterial-derived, so a host
  mitochondrial contig - which a GetOrganelle nr run can recruit if the seed DB
  contains any mt sequence - reads as 16S/23S and would be labelled bacterial.
  Run barrnap's mito model explicitly and require any prokaryotic call to be
  corroborated by a prokaryotic taxid path. And 5S rRNA is a normal nuclear
  metazoan gene, sometimes linked to the 45S array: add it to the parse map as
  a **eukaryotic** feature, not a prokaryotic signal.
- ITSx is run with `-t all` (inst/nextflow/modules/pipeline_r.nf:143-144),
  scoring against ~18 taxon HMM profile sets, and RiboPilot parses only the
  coordinate line and the literal string `Chimeric!` (R/annotate.R:58-88). The
  winning profile is written to disk and thrown away. Harvest it.
- **But do not promote SSU/LSU profile disagreement to a chimera detector.**
  ITSx was built for fungal ITS barcoding and its Metazoa profile set is far
  weaker than Fungi, so a divergent asteroid, scyphozoan or octocoral routinely
  gets its SSU flank assigned to one group and its LSU flank to another with no
  chimerism whatsoever. Demote to weak, cross-kingdom-only corroboration.
- **Hard constraint on this channel.** For a metazoan contaminant in a metazoan
  sample, ITSx and barrnap return Metazoa/eukaryote for **both**. If the channel
  then returns a negative (on-target) score it actively votes to *clear* the
  contaminant. So: **tool_kingdom may only ever push toward off_target, never
  toward on_target.** Same kingdom returns score 0 with
  `availability = "uninformative"`, distinct from `"unavailable"`.

**Status 2026-07-29: channel (b) is RETIRED, not deferred.** Its stated
justification here is that it is the only channel working fully offline, and the
maintainer has confirmed users never run RiboPilot offline. What remains is
covered elsewhere or is redundant: a bacterial or fungal contaminant is already
caught by channel (a)'s taxonomy, which this section itself calls "heavily
correlated" with (b); the prokaryotic call would in any case need corroboration
"by a prokaryotic taxid path", which is channel (a) again; and a mitochondrial
assembly is caught offline-independently by `ORGANELLAR_TITLE_RE` in
`dev/design/10_HIT_CLASSIFICATION.md`. The only case genuinely uncovered is a
scaffold with zero usable GenBank hits, which the screen already reports
honestly as `no_reference` rather than guessing.

Two pieces are worth salvaging on their own merits, because the pipeline pays
for the computation and discards the answer: `parse_barrnap()`'s name map
(R/annotate.R) recognises only 18S, 5.8S and 28S, so a 16S or 23S call is
silently dropped; and ITSx writes its winning taxonomic profile to disk while
`parse_itsx()` reads only the coordinate line.

**(c) Depth**, as a within-sample ratio, never deciding alone.

The intuition is right but the usual arithmetic is wrong by more than an order
of magnitude. rDNA read fraction scales as
`(copies x repeat length) / genome size`, **not** as DNA mass fraction. A fungus
at 1% of extracted DNA has ~100 rDNA copies in a ~30 Mb genome versus a fish's
~200 copies in a ~1 Gb genome, so it contributes roughly **15-20%** of the
host's rDNA depth, not 1%. Bacterial per-genome-equivalent rDNA fraction is
several times a vertebrate's. Fungal growth on museum tissue and bacterial load
are the two commonest contaminants in this material, so the "real and large"
separation does not exist for them. State this explicitly wherever the
thresholds are recorded, or someone will later tighten them to match the wrong
model.

- Populate `assemblies.depth` (currently declared at R/init_db.R:283 and written
  by nothing) from `scaffold_depth(header, parse_gfa_depth(gfa))`. The segment
  ids survive: `format_assembly.sh:20-21` keeps the
  `scaffold_N--<edge ids>` header as the third token and `load_assemblies`' `hdr`
  is the whole line (R/pipeline_steps.R:232-236).
- **`rdna_depth_ref` must be the deepest unit whose taxonomy is
  host-concordant**, not merely the deepest anchor-bearing unit. Otherwise the
  definition is circular exactly when it matters: a heavy fungal bloom, a
  myxozoan-infected fish, a symbiont-rich coral. In those cases every genuine
  host scaffold falls below the floor and receives "strong off-target evidence",
  inverting the call on the whole sample. When no unit is host-concordant, mark
  the channel unavailable rather than picking a reference.
- Within a path use read depth; across paths use graph depth, and **record that
  graph depth is also biased**: the 18S and 28S conserved cores of host and
  contaminant share exact k-mers at GetOrganelle's k values, so those segments
  carry summed depth. The very collapse that generates the alternative paths
  also inflates a contaminant's graph depth. Neither basis is unbiased; record
  which was used in `depth_basis`.
- Skip the channel entirely when the sample has fewer than three units - a
  2-unit ratio has no meaningful reference - or when `coverage` did not run
  (coverage.nf carries `errorStrategy 'ignore'`).
- Thresholds as evidence weight, not verdict: `>= 0.33` no penalty (fragments of
  one array legitimately differ 2-3x from conserved-core read collapse and
  contig-edge effects); 0.05-0.33 weak; `< 0.05` strong; `< 0.01` very strong.
  High ratio is **never** on its own evidence of being on-target.

**(d) Positional / chimeric.** From the HSP windows (section 3.4), not from new
searches. But an off-target window must **also exceed the identity envelope for
its own rank distance** before it counts. Without that, the rule fires from
**reference availability asymmetry**, not chimerism: GenBank's 18S coverage
across metazoan phyla is far denser than its 28S coverage, so a scyphozoan's 18S
window best-hits a scyphozoan at d=1 while its 28S window best-hits whatever
distant taxon has a 28S record matching the conserved core - a phylum-level
subject change across windows on a perfectly clean, single-organism scaffold.
Require comparable HSP lengths on both sides so a 200 bp sliver cannot outvote a
2 kb on-target block. Never call a chimera on a *gap*.

**(e) Cross-sample.** See 3.6.

**(f) GC**, weight 0.05, tie-break only, never able to move a verdict alone.
Include a test asserting that GC alone can never produce a non-uncertain
verdict. Why it is weak: bacterial 16S sits at 50-56% GC essentially
independent of genome GC and overlaps the metazoan rDNA range directly, so the
whole-genome blobplot intuition does not transfer to a conserved rRNA locus.
And whole-scaffold rDNA GC is dominated by **28S expansion-segment content**
(D-segments reach 75-80% in vertebrates, near 50% in many invertebrates), so GC
tracks which clade a scaffold is from and how much 28S expansion it happens to
contain, far more strongly than whose it is. Once annotation exists, compute GC
on the **anchor regions only**.

**Calibration note.** Channels (a) and (b) are heavily correlated - both are
essentially "is it a fungus" - so the effective weight on the cross-kingdom axis
is large and the **within-kingdom cases are decided almost entirely by depth and
cross-sample, the two weakest channels**. That is the design's known soft spot.
`screen_reason` must list every channel with its score **and** its weight, so a
user can see when a verdict rested on one channel.

### 3.6 Cross-sample contamination and index hopping

Nothing anywhere compares a sample against the other samples in its own project.
ExAmp index hopping moves ~0.1-2% of reads between multiplexed libraries, rDNA
is the highest-copy locus in a skim, and museum projects are pooled at high
multiplex - so hopped rDNA is both the most abundant hopped sequence and the
most likely to assemble. A hopped scaffold looks perfectly clean to every
per-sample check: real metazoan taxonomy, high identity, Metazoa ITSx profile.

Mechanically cheap: build one local BLAST db from every non-ignored scaffold in
the project (headers `ID.path.scaffold`), blast each scaffold against it. A few
hundred sub-15 kb sequences, seconds, entirely local.

But the naive taxonomic guard fires systematically on the normal project shape.
"Not congeneric (`rank_distance >= 2`)" still lets an 18S identity rule apply to
same-family and same-order pairs, and metazoan 18S is routinely > 99.5%
identical within a family and frequently across families within an order. **A
museum project is almost always one clade on one plate, so this is the common
case.** And the congeneric fallback is logically inverted: it is invoked exactly
when the samples ARE congeneric, which includes **conspecific** samples -
multiple specimens of one species on one plate, the commonest plate composition
in a systematics project - whose ITS is genuinely ~99.9-100% identical.

Corrected rule:

- Require the **ITS-based** test whenever `rank_distance <= 4` (same class), not
  only for congeners. 18S carries no discriminating information at that range.
- Require `>= 150 bp` of **ITS-internal** alignment (many metazoan ITS2 are
  shorter than 300 bp).
- Replace fixed cutoffs with "observed identity exceeds the expected
  within-clade identity for this pair's rank distance", reusing the envelope
  machinery rather than inventing a second scalar.
- Depth asymmetry `<= 0.02` (tightened), required in **both** branches. Hopping
  donates a small fraction, so the contaminated copy is always the shallow one.
- **Ship `samples.batch` in the same change.** Without run/lane/plate metadata
  the screen compares specimens that were never on the same flowcell, and every
  hit it produces there is meaningless by construction.
- Evidence only. Never a verdict on its own.

**The free win nobody proposed.** Restrict that same all-vs-all database to
pairs **within one sample** and it answers, quantitatively and at zero marginal
cost, the exact question dev/design/07_SCAFFOLD_SELECTION.md:10-17 names as the
core one: same region redundant / same region divergent / different region /
different organism. See 3.8 and 4.4.

### 3.7 The verdict, and who owns the per-unit opinion

There is no per-unit contamination verdict anywhere in the schema. The machine's
judgement exists only as a substring inside `annotate.annotate_notes`
(R/curate_rrna.R:126, :130, :136 -> :305), a column `curate_one()` overwrites
wholesale on every run and the modal textarea also writes to
(R/app_server.R:1239). The two clobber each other.

**Two layers, cleanly separated** (the design tracks proposed two competing
models; this is the reconciliation):

- **`screen`** (machine): one row per `(ID, path, scaffold)` holding
  `screen_verdict`, `screen_class`, `screen_confidence`, a human-readable
  `screen_reason`, and the evidence record. Recomputable, versioned
  (`screen_version`), never user-written.
- **`assemblies.triage`** (action): `keep` / `drop` / NULL(undecided), plus
  `triage_reason` and `triage_source` (`auto` | `user`). `triage_source = 'user'`
  is sticky and no re-run may overwrite it. `assemblies.ignore` stays the single
  downstream gate and becomes derived, never hand-written.

Verdicts: `on_target`, `off_target` (sub-classed bacterial / archaeal / fungal /
algal_plant / other_euk / metazoan_other / human / cross_sample),
`chimeric`, `paralogous_variant` (see 3.8), `not_rdna`, `no_reference`,
`uncertain`.

Gating: emit a contamination verdict only at confidence `>= 0.6`; only at
`>= 0.8` does the UI **offer** to drop. Below 0.6 the verdict is `uncertain`
with the evidence still stored. The point is to direct attention, not to decide.

**Auto-drop only on positive evidence, never on absence of a close reference:**

1. barrnap `-k all` calls 16S+23S with no 18S **and** the hit lineage is
   Bacteria/Archaea by taxid intersection (and not the host's mitochondrion).
2. Cross-superkingdom taxid intersection confirmed **and** identity exceeds the
   envelope for that rank distance.
3. No rRNA anchor from any kingdom, no ITSx region, best hits are
   genomic/organellar, and the scaffold is short.

Everything driven by rank distance alone becomes **undecided/amber with a reason
string**, never a drop. The bulk-action label reads "Flag for review", not
"Drop".

Also in this change:

- Split `annotate_notes` into `curate_summary` (machine, overwritten) and
  `annotate_notes` (user, durable). Backfill by copying the existing value into
  **both** rather than guessing which it was.
- Deprecate `contaminant` / `contaminant_unchecked` / `not_nuclear_rdna` rather
  than deleting them: `tests/testthat/test-curate.R:64-75` asserts all three by
  name and `curate_rrna_core()` is exported. Keep them emitting for one release
  behind a params switch, keep the arguments accepted-and-ignored, and update
  the test in the same commit.
- `no_reference` reason codes replace `contaminant_unchecked`: `no_hits` /
  `no_taxdump` / `unresolved_sample_lineage` / `unresolved_hit_lineages` /
  `only_distant_low_identity`. Neutral grey. Gates nothing. Tell the user in the
  UI that the tool_kingdom channel still runs in all of those states, so an
  offline project still catches bacteria and fungi and only loses the ability to
  separate two metazoans.
- Sample-level `reference_support` (minimum rank distance across the sample's
  best on-target hits). When `>= 4` (class or worse), relax sample-wide and
  **say so visibly**: "nearest reference is at class distance; identity checks
  relaxed". Compute it from **modal-taxon units only**, recompute iteratively
  after any unit is marked off_target, and exclude uncultured/environmental
  taxids - otherwise a sample whose only hits are uncultured clones scores as
  having reference support at superkingdom and never relaxes.
- A loud sample-level warning when **the deepest unit's taxon does not match
  `samples.lineage`**. One query. This is the only defence against the
  high-abundance-symbiont case (Symbiodinium in a coral, a myxozoan infection, a
  fungal bloom) where the depth channel and the taxonomy channel fail together
  and every other mechanism computes its reference from the contaminant.

### 3.8 The third outcome: intragenomic rDNA variation

Every proposal in the first pass framed multiple scaffolds as either "graph
isomer of one molecule" or "contaminant". **The third case is common in the
stated target clades and is neither.** Metazoan rDNA arrays are frequently not
concertedly evolved: Octocorallia, Scyphozoa and several other groups carry ITS
copies differing by several percent within one individual. GetOrganelle emits
those as distinct paths or scaffolds. They are real host biology, not artifacts
to reconcile and not contaminants to drop, and **which copy the curator locks
determines what GenBank receives.**

Add `paralogous_variant` as a first-class verdict, driven by within-sample
pairwise identity between host-concordant units over their overlapping interval
(free, from 3.6's database restricted to one sample). Without it,
sibling-concordance and any pseudogene test both mislabel it. `paralog_suspect`
in particular must compute its anchor-divergence floor from **the within-sample
spread among host-concordant units**, not from a fixed per-element value, or
real intragenomic array variation is reported as pseudogenes.

Related, and also missing: the IGS is nowhere named as the **cause** of the
isomer paths. GetOrganelle emits multiple paths for rDNA chiefly because the
intergenic spacer carries its own tandem subrepeats the graph cannot resolve, so
alternative traversals differ only in subrepeat copy number and are the same
molecule. Two scaffolds differing solely in IGS length currently render as
identical anchor tracks with different lengths and no explanation. The curator
should be told "these paths differ only in IGS subrepeat copy number, pick
either" rather than being asked to adjudicate.

### 3.9 Where it plugs into the pipeline

```
                                  main.nf (WF1)
  preprocess -> assemble -> score -> blast_genbank -> blast_ref_fetch
                             |
                             +---> coverage
                             |
                             +---> screen_tools        (NEW, local only)
                                        |
     blast_genbank.out.join(coverage.out).join(screen_tools.out, remainder: true)
                                        |
                                     screen            (NEW, local only)
                                        |
                            cross_sample_screen        (NEW, .collect() fan-in)
```

- `screen_tools`: `pybarrnap -k all` (plus the mito model) and
  `ITSx --detailed_results T` per unit. **Local only** - this is the arm that
  still works on an offline HPC node.
- `screen`: reads `blast_hits`, `coverage`, `assemblies`, `samples.lineage` and
  the tool columns; writes the verdict. Its dependency on coverage is soft
  (coverage.nf carries `errorStrategy 'ignore'`), so use
  `.join(..., remainder: true)` and have `screen_one()` mark absent channels
  unavailable rather than dropping the sample silently.
- `cross_sample_screen`: `errorStrategy 'ignore'` plus standalone re-runnability
  from R, so a `.collect()` fan-in never becomes the critical path. Nothing
  downstream consumes it in the same run.
- WF2 refinement folds into `curate_one()`, not a new process:
  `refine_screen_verdict()` re-runs the envelope against **actual annotated
  element boundaries** and `blast_ref_alignment`'s per-element pid/aln_cov,
  which are more precise than HSP-derived windows. This is where "the 18S is
  on-target but the 28S is fungal" becomes definitive.
- Do **not** wire the screen into `main_fallback.nf` until it is repaired
  (section 1.9).

**Explicitly not recommended: blobtools.** This design already produces both
blobplot axes (per-unit depth, per-unit bitscore-weighted LCA) at finer
resolution, without a diamond/nt dependency. Build the blob *view* in the app
from DB columns; skip the tool. Remove or annotate the commented blobtools line
at environment.yml:29.

---

## 4. Multi-scaffold UX: the target design

### 4.1 Tri-state, non-destructive triage

Replace the binary, reason-free, destructive `assemblies.ignore` toggle:

```
              +-----------+   click D           +-----------+
  (auto)      | KEEP      |-------------------->| DROP      |
   ---------->| ignore=0  |<--------------------| ignore=1  |
              | green     |   click K           | grey      |
              +-----------+                     +-----------+
                    ^                                 |
                    |         +--------------+        |
                    +---------| UNDECIDED    |<-------+
                              | triage NULL  |
                              | amber        |
                              +--------------+

  annotate / annotations / coverage rows SURVIVE every transition.
  triage_source='user' is sticky: auto re-runs skip those rows.
```

`reconcile_ignore_unit()` stops `DELETE`ing and instead holds the annotate row
(`annotate_switch = 0`). The row is already invisible to the app, to WF2 and to
export - all three filter `COALESCE(asm.ignore,0) = 0` - so deleting it buys
nothing and costs the curator's work. `reconcile_ignore_unit()` has exactly one
caller (R/app_server.R:877) and zero test coverage, so this refactor is low risk.

Narrow the orphan cleanup at R/pipeline_steps.R:307-313 so a WF1 re-run deletes
only units whose scaffold genuinely vanished. This is only meaningful **after**
section 1.2 - otherwise the ignore flags are already gone by the time it runs.

### 4.2 Assemble table: sample rollup with an agreement strip

Keep Assemble at sample grain (correct for scanning 100 samples). Make every
column earn its width.

```
 Assemble  [State][Lock][Triage v][Update]     Show: (o) Needs review ( ) Clean ( ) All
 +--+--+---------------+--------------+---------+---------------+--------------------+-----------------+---+
 |  |  | ID            | Taxon        | Units   | Agreement     | Best hit (rank)    | Repeat covered  |   |
 +--+--+---------------+--------------+---------+---------------+--------------------+-----------------+---+
 |[]|OK| USNM 123456   | Etheostoma.. | 4 / 4   | * * * *       | Etheostoma (genus) |18S I1 5.8 I2 28S| > |
 |[]|OK| USNM 123457   | Etheostoma.. | 3 / 12  | * * * x x x x | Etheostoma (genus) |18S I1 5.8 I2 28S| > |
 |[]|! | USNM 123458   | Notropis lu..| 2 / 3   | * ? x         | Notropis (genus)   |18S I1 5.8 -- -- | > |
 |[]|! | USNM 123459   | Hypentelium..| 0 / 6   | x x x x x x   | Saprolegnia (none) |-- -- --- -- --  | > |
 |[]|  | USNM 123460   | Octocorallia | 5 / 5   | . . . . .     | no close reference |18S I1 5.8 I2 28S| > |
 +--+--+---------------+--------------+---------+---------------+--------------------+-----------------+---+

  Agreement strip, one dot per unit in path.scaffold order:
    *  green   on_target
    ?  amber   undecided / odd-one-out          <- needs attention
    .  GREY    no_reference / offline            <- NORMAL, not a problem
    x  red     off_target / not_rdna
    o  hollow  dropped
  hover any dot -> "2.1  Notropis lutrellus  5.9 kb  shares: Notropis (genus)"

  row 1: nothing to do, Lock it.
  row 3: one undecided + one off-target, expand.
  row 4: everything off-target, the assembly failed; reason on hover.
  row 5: sparse clade with no close reference. Grey, NOT amber. This is the
         tool's core use case and must not arrive looking broken.
```

Notes that matter:

- **`no_reference` must render neutral grey, visually distinct from amber.** For
  an offline project or a sparse clade this is the normal state, and if it reads
  as "needs attention" every sample in the tool's core use case arrives looking
  broken. Same for the anchor strip: a hollow ITS1 pill on a 28S-only fragment
  means "this scaffold is a different region", not "something is missing".
- **Do not use "modal species" as the consensus column.** For the target taxa
  each scaffold's single top hit is driven by which region it covers and which
  records exist for that region, so the modal species across six scaffolds of
  one octocoral is arbitrary and a "(5/6)" suffix reads as a contamination
  signal on a clean sample. Show the **deepest-majority taxon plus its rank**
  ("Notropis (genus, 6/6)" / "Cnidaria (phylum, 5/6)"), and a species-level
  string only when a genus-or-better majority exists.
- **Keep a length column.** A 12 kb complete repeat versus a 900 bp
  conserved-core fragment is the primary sort key for deciding what to look at.
  Use `max_len` (longest kept unit) with the full list in a tooltip, replacing
  the `;`-joined string.
- **Label "Repeat covered" as the union across kept units**, so a curator does
  not read 18S+ITS1+5.8S+ITS2+28S as one complete molecule when it is five
  fragments.
- `units` reads "kept / total", red when kept is 0.

### 4.3 Where the keep/drop action lives

The obvious design is expandable reactable sub-rows with keep/drop buttons in
them. **Verified against the installed reactable 0.4.5, that does not work as
assumed:**

- The `details` callback is **eager**, not lazy. It fires for every row at
  widget build time, whether or not anything is expanded. So caching the
  per-unit frame in one query per refresh and subsetting in R is **mandatory**,
  not an optimization.
- Sub-table HTML is baked into the column definitions at render time, not into
  `data`. So the data-only `updateReactable` at R/app_server.R:280-282 does
  **not** refresh sub-rows: after a triage click they show stale keep/drop state
  until a full re-render, which clears the client-side filter/sort/page that the
  `isolate()` at R/app_server.R:156 exists to protect.

**Recommendation: sub-rows are read-only evidence; the action lives in the
modal**, where a re-render is free. This resolves the conflict cleanly and is
also better ergonomics: the sub-row answers "why is this sample amber" in place,
and the modal is where you act.

```
 |[]|! | USNM 123458   | Notropis lu..| 2 / 3   | * ? x     | Notropis (genus) | v |
 +--+--+---------------------------------------------------------------------------+
 |     |Path.Scaf| Length | rRNA genes       |Depth| Top GenBank hit    | Shares   |
 |     |   1.1   |  6,120 |18S I1 5.8 I2 --  |#### | Notropis lutrellus | genus    |
 |     |   1.2   |  1,442 |--- -- --- -- 28S |##   | Notropis lutrellus | genus    |  amber, undecided
 |     |   2.1   |    890 |--- -- --- -- --  |#    | Saprolegnia paras. | none     |  red, auto-dropped
 |     |         |        |                  |     | why: 99.1% to an oomycete, implausible at        |
 |     |         |        |                  |     |      cross-kingdom distance; depth 0.008x        |
 +--+--+---------------------------------------------------------------------------+
```

### 4.4 The Assemble modal becomes the evidence surface

**Status 2026-07-30: the implementation has moved past this sketch.** The modal
(R/app_modals.R:42, R/app_server.R:875) now has a scaffold table where clicking
a row examines it, highlighted via reactable meta rather than a row index, and
the two evidence panels below the table (GenBank screen plus nuclear-rDNA
candidates, and read coverage) follow the examined row. The tick boxes are a
separate multi-select used by the Copy FASTA button, which copies the ticked
scaffolds headed by sample, taxon and unit. The old per-row "hits" text link is
gone. `set_selected` is still present, contrary to the paragraph after the
sketch. Read the sketch below as intent, not as a description of the app.

```
+-- Assemble - USNM 123458 -------------------------------------------------[x]+
| Notropis lutrellus - taxid 123456 - species                                  |
| Auto-triage kept 2 of 3 scaffolds. 1 undecided. Repeat coverage: 18S I1 5.8  |
|                                                                              |
| Scaffolds                             [filter] [search]      [ ] select all  |
| +--------------------------------------------------------------------------+ |
| |Keep |P.S | Length| rRNA genes       |Depth| Top hit          |Shares |Why | |
| |[K][]|1.1 |  6,120|18S I1 5.8 I2 --  |#### |Notropis lutrellus|genus  |    | |
| |[][ ]|1.2 |  1,442|--- -- --- -- 28S |##   |Notropis lutrellus|genus  |    | |
| |[][D]|2.1 |    890|--- -- --- -- --  |#    |Saprolegnia paras.|none   |... | |
| +--------------------------------------------------------------------------+ |
|                                                                              |
| v Read coverage (selected scaffold)      [plot_coverage(), already exists]   |
|                                                                              |
| v GenBank hits for scaffold 2.1                                              |
|   SCREEN - best match anywhere in core_nt (NOT restricted to nuclear rDNA):  |
|     Saprolegnia parasitica 18S-ITS-28S  MH123456  99.1% id  87% cov  3.1 kb  |
|     nuclear rDNA: YES   lineage: Eukaryota; Stramenopiles; Oomycota; ...     |
|   PANEL - best matches among nuclear-rDNA-titled records:                    |
|     1 MH123456 Saprolegnia parasitica  99.1  87  3.1 kb                      |
|     2 MH999888 Saprolegnia diclina     97.4  85  3.0 kb                      |
|   ^ these two blocks naming different organisms IS the contamination call    |
|                                                                              |
| v Compare selected scaffolds                                                 |
|   path 1              0        2k       4k       6k                          |
|   1.1  [====18S====][I1][5.8][I2][=====28S=====]   120x  Notropis  genus     |
|   1.2                               [==28S==]       48x  Notropis  genus     |
|                                     ^ overlaps 1.1's 28S over 1,380 bp       |
|                                       at 99.9% identity -> same region,      |
|                                       redundant, drop either                 |
|   path 2                                                                     |
|   2.1  [~~~~~ no rRNA genes ~~~~~]                    9x  Saprolegnia  none  |
|                                                                              |
|           [Keep selected] [Drop selected] [Keep ONLY selected] [Close]       |
+------------------------------------------------------------------------------+
```

Delete `output$sel_path_ui` and `input$set_selected` (section 1.4). Make
`sel_path_ui`, `assemble_scores` and the scaffold table all read `refresh()`
through one shared reactive - today none of them do, so they display stale state
after a write.

**The pairwise-identity number in the compare panel is the load-bearing part.**
A track diagram alone cannot distinguish "two scaffolds cover the same 28S
region" at 0.1% divergence (assembly noise, pick either) from 3% divergence
(real paralogous array copies, and your submitted sequence is one of several).
One pairwise identity over the overlapping interval settles it, and 3.6's
database provides it free. Report identity, overlap length and strand.

**Do not port MitoPilot's MSA conflict navigator**
(MitoPilot/R/app_assemble_coverage_details.R:409-733: DECIPHER MSA, conflict-block
RLE, per-column depth/error/GC projection, Path 0 consensus builder, ~1200
lines). That subsystem exists because a mitogenome is one molecule whose isomers
must be reconciled into a single record. Tandem-repeat rDNA scaffolds are
alternatives to **choose between**, and the curator's question is regional
overlap, not nucleotide conflict. A track diagram plus one identity number
answers the real question at roughly 5% of the cost.

### 4.5 Bulk triage with an honest preview

```
+-- Triage scaffolds - 37 selected samples ---------------------------------[x]+
|                                                                              |
|  (o) Re-run automatic triage with the current settings                       |
|  ( ) Drop scaffolds whose barrnap kingdom is bacterial / archaeal            |
|  ( ) Drop scaffolds with no rRNA genes found and no ITS region               |
|  ( ) FLAG FOR REVIEW: best hit is a different  [ phylum   v ]                |
|  ( ) Mark every undecided scaffold as Keep                                   |
|                                                                              |
|  [ ] also change scaffolds I decided by hand      (12 would change)          |
|  [ ] also change scaffolds whose annotation is locked (3 would change)       |
|                                                                              |
|  Would drop 214 of 396 scaffolds across 31 samples:                          |
|     188 have a closer on-target hit available  (foreign organism)            |
|      26 have NO closer reference anywhere      -> marked needs-review, kept  |
|  4 samples would end up with NO scaffolds kept: USNM 123459, ...             |
|                                                                              |
|                                         [ Apply ]  [ Cancel ]                |
+------------------------------------------------------------------------------+
```

The preview is the safety design and none of it should be traded away. Two
things the first-pass version got wrong:

- The rank predicate must **flag**, not drop. "Drop scaffolds whose best GenBank
  hit is a different phylum" hands a biologist a one-click way to apply the
  rank-distance-as-contamination error across a whole project. A project of
  asteroids run through "different phylum" loses most of its real scaffolds in
  one transaction.
- The preview must **split the would-drop count** into "has a closer on-target
  hit available" (real evidence of a foreign organism, droppable) versus "no
  closer reference exists anywhere in the hit list" (`no_reference`, the normal
  state, marked needs-review only).

Add "drop scaffolds whose barrnap kingdom is bacterial/archaeal" as a rule: it
is positive evidence, works offline, and is the rule most projects actually
need. Keep the `annotate_lock` skip default-on. Apply inside one transaction.

### 4.6 Annotate and Export tables

- Sticky `unit` column ("p1.s3") after ID; Path/Scaf on the MitoPilot rule via
  reactive CSS (section 1.8).
- `defaultSorted = list(ID = "asc", path = "asc", scaffold = "asc")` so a
  sample's units are adjacent, plus a `blk` parity column and a group stripe so
  each sample reads as one band. Currently `defaultSorted` is
  `time_stamp desc` (R/app_server.R:247), interleaving units from different
  samples, and `rt_highlight_row()` cancels `striped` (section 1.9), so there is
  no separator at all. Keep Last Updated sortable.
- **Make flags filterable.** `n_flags` is a bare integer with
  `filterable = FALSE` (R/app_server.R:185-186) and the flag *names* live in
  `annotate_notes`, which the Annotate query does not select
  (R/app_server.R:135-146). So "show me the contaminants" is currently
  impossible. Select the notes, render as a filterable `rt_longtext()`. This is
  the single most-requested triage in a 200-sample project, and it is also what
  rescues the flag column in the interim while `missing_*` semantics are still
  per-scaffold (section 1.7).
- Carry `shares` (deepest shared taxon + rank) and a `not_rdna` badge into
  Annotate and Export, so the verdict formed at Assemble time stays visible
  while curating and at the moment of submission.
- Export tab: add the SeqID column read-only, so what GenBank will see is
  visible **before** the bundle is written. Carry `triage_reason` /
  `screen_reason` into `export_summary_df()` (R/export.R:245-262): why a
  scaffold was left out is part of the submission record for a systematist.

### 4.7 Export gating

Add `off_target` and `not_rdna` at confidence `>= 0.8` as hard blocks, after
`effective_verdict()` has consulted the user override. Emit a hard error naming
the units, the verdict, the confidence **and the single channel that contributed
most**, so the user can judge it in one read.

**`chimeric` must be a loud warning, not a block.** Every chimera detector here
has a false-positive profile driven by the target taxa: ITSx profile
disagreement fires on divergent invertebrates, HSP-window rank disparity fires
from 18S-vs-28S reference availability asymmetry, and neither can see the
commonest real rDNA chimera (host variant crossed with host variant). A hard
block converts each of those into lost submissions, and buys little, because
GenBank runs its own contamination screening on submitted records anyway.

**`no_reference` must never block.** For the tool's core use case it is the
expected state.

Move the one-locked-path-per-sample invariant out of the two Shiny observers
(R/app_server.R:485-491, :1293-1301) and into `export_bundle()`, so a
programmatic export cannot emit alternative resolutions of one rDNA repeat as
independent GenBank records.

---

## 5. Larger and optional

- **Per-position minor-allele fraction.** A submitted rDNA scaffold is a
  consensus over 100-10,000 non-identical array copies, and that is the first
  question a systematist gets asked about a genome-skim rDNA submission.
  `coverage.depths` already stores the per-position track; a per-position minor
  allele fraction over the mapped reads would flag positions where the array is
  genuinely heterogeneous and let the curator emit IUPAC ambiguity codes rather
  than an arbitrary majority base. Submission quality, not contamination.
- **Museum DNA damage.** Cytosine deamination in old and formalin-exposed
  material inflates apparent divergence, systematically pushing museum samples
  toward `low_identity` flags and `no_reference`, and away from tripping the
  upward-flagging envelope. The per-element floors in R/params_rrna.R:60-63 are
  calibrated for fresh tissue. At minimum document this next to the floors;
  better, estimate damage per sample from C-to-T vs G-to-A asymmetry at read
  termini and widen the floors accordingly.
- **Local tiled BLAST against the clade seed DB.** `score_one()` already builds
  it per sample (R/pipeline_steps.R:21-22), so tiling costs only CPU. Useful to
  **corroborate** an HSP-window conflict and localize a junction so the user can
  trim there (`trim_assembly()` and the trim UI already exist). It must never
  trigger alone: a legitimately divergent host region (a 28S expansion segment,
  an ITS in an under-sampled clade) drops seed support and mimics a junction.
- **A second, taxonomy-free junction detector** for the host-variant chimera,
  from within-sample pairwise identity plus a per-position depth discontinuity.
  Taxonomic windowing is blind to it by construction (both windows name the same
  species) and it is the case that most damages a submitted reference sequence.
- **Doc and dead-code reconciliation.** Rewrite
  dev/design/07_SCAFFOLD_SELECTION.md for the per-scaffold model (it currently
  says the Level-1 stitch is built and needs integrating, the opposite of what
  b0f1b12 did) and record the decision to drop stitching, which no document
  currently records. Deprecate rather than delete the exported pre-pivot
  functions (`annotate_project`, `score_project`, `stitch_scaffolds`;
  dev/prototype_stitch_viz.R still calls one). Delete
  `compute_feature_alignment()` outright - zero callers, absent from NAMESPACE,
  and actively wrong under the per-scaffold model.

  **Status 2026-07-30: the doc half is DONE, the dead-code half is PART DONE and
  the rest is blocked on open question 2.**

  Done: `07_SCAFFOLD_SELECTION.md` now carries a status section recording that
  stitching was dropped by `b0f1b12`, why, and what the per-scaffold model
  replaced it with. That section is the record this bullet says no document
  holds. Its three stale in-body claims (the Level-1 heading, "Remaining:
  integrate stitch", and P1) are marked in place rather than rewritten, so the
  pre-pivot reasoning stays readable.

  Done: `compute_feature_alignment()` is deleted. Its successor
  `compute_feature_msa()` in the same file already does the job with full `ID`,
  `path` and `scaffold` predicates.

  **Correction: "deprecate the exported pre-pivot functions" understates what is
  involved for `stitch_scaffolds`.** It is not merely called by
  `dev/prototype_stitch_viz.R`. `R/pipeline_steps.R:394` and
  `R/annotate_project.R:89` both call it in package code, and it is exported
  alongside `stitch_one`. Neither caller is reachable from the main pipeline -
  the first only through `main_fallback.nf`, the second only if a user invokes
  the legacy `annotate_project()` driver by hand - but "unreachable from
  `main.nf`" is not the same as "dead", and deprecating a live-in-fallback
  export is exactly the piecemeal move open question 2 warns against. Left in
  place deliberately, pending answers to open questions 1 and 2.

  **New consequence, recorded rather than acted on.** Deleting
  `compute_feature_alignment()` removed the last reader of the `stitched` table.
  `R/init_db.R:447` still creates it on every fresh project and no migration
  drops it, so it is now created, never written on any reachable path, and never
  read. That makes dropping it part of open question 2's bundle, not a separate
  cleanup.

  **Closed 2026-07-31.** Open question 2 was answered and the whole bundle,
  including the exports this bullet proposed only deprecating, was deleted.
  Schema v18 drops the `stitched` table in both `new_db()` and the migration.
- **Validation gate.** The envelope numbers and channel weights are unfitted
  priors, and a false positive on an under-represented taxon is the failure that
  makes the tool useless. **Auto-drop rules and the export gate should ship
  disabled by default until a labelled fixture set exists in `tests/`**,
  validated against the asteroidea / scyphozoa / octocorallia invert clades plus
  a deliberate fungal and human spike-in. Make that a gate, not a note.

  **Status 2026-07-30: the gate SPLITS into two deliverables.**

  The spike-in half is delivered by `dev/design/13_VALIDATION_FIXTURES.md`: a
  labelled corpus in `tests/testthat/fixtures/`, built from the reference
  project's 19 clean-host units plus deliberate cross-kingdom spike-ins whose
  ground truth is known by construction. It measures the false-positive rate,
  which is the number this bullet is really about. It changes no default.

  The clade half is BLOCKED ON DATA, not on code. The three invert read sets are
  not downloaded, their clade reference DBs are not built, and no invert project
  has ever been run: `mapping_asteroidea.csv`, `mapping_scyphozoa.csv` and
  `mapping_octocorallia.csv` each define a single sample, and every test project
  to date is actinopterygian. Satisfying it needs read downloads, three reference
  DBs, three assemblies and three remote BLAST runs before any labelling starts.

  **The human spike-in this bullet asks for is not a cross-kingdom case.** `Homo
  sapiens` is Chordata and so is a ray-finned fish, so on a fish host a human
  contaminant shares the phylum slot, passes concordance, and is MISSED. It is
  caught only on a non-chordate host, which is another reason the clade half
  matters. Doc 13 section 4.1 has the measurement; doc 11 section 1 is corrected
  to match.

---

## 6. Considered and rejected

| Idea | Why not |
|---|---|
| blobtools | This design produces both blobplot axes at finer resolution with no diamond/nt dependency. A blobplot's power comes from a wide GC spread across whole-genome contigs; a conserved rRNA locus has almost none (bacterial 16S sits at 50-56% GC, overlapping metazoan rDNA directly). Build the *view*, skip the tool |
| Windowed **remote** BLAST for chimera detection | Multiplies the pipeline's slowest, most failure-prone step for information the HSP coordinates already contain, once `-max_hsps 1` is lifted |
| Porting MitoPilot's MSA conflict navigator | ~1200 lines answering a question rDNA does not ask. Isomers are alternatives to choose between, not variants to merge. A track diagram plus a pairwise identity number costs ~60 lines |
| Dropping `assemble.blast_*` and `assembly_blast` from the DDL | Violates the migration contract at R/migrate_db.R:41-42; fresh projects would diverge from migrated ones. Stop selecting them instead |
| SQL views for `blast_screen` / `blast_ref_candidates` | Both are written by `DELETE` + `dbAppendTable` (R/blast_genbank.R:45-48, :93, :126-132); views are not writable, and `migrate_db`'s `has_tbl()` uses `dbListTables`, which returns views, so the v7/v8 guards would skip their rebuild |
| Genus-level sibling concordance | Flags every unit on sparse clades, where nine scaffolds of one octocoral top-hit nine different genera. Use deepest-rank-with-a-strict-majority instead |
| Auto-dropping on rank distance alone | The exact failure mode that makes the tool useless for the taxa it exists to serve. Flag, never drop |
| A green "host" badge from `lineages_conflict()` | Converts absence of evidence into a positive assertion, on a test that cannot fire within a kingdom |
| `hit_taxa = COUNT(DISTINCT blast_species)` as a free contamination hint | `extract_species()` takes the first two defline words, so "Danio rerio", "PREDICTED: Danio" and "Danio sp." are three distinct values, and a clean sample's 18S-bearing and 28S-bearing scaffolds legitimately top-hit different congeners. A false-positive generator |
| ITSx SSU-vs-LSU profile disagreement as a primary chimera detector | ITSx was built for fungal ITS barcoding; its Metazoa profile set is far weaker, so divergent invertebrates get discordant flank assignments with no chimerism. Weak cross-kingdom corroboration only |
| Reactable sub-rows as the triage control surface | The `details` callback is eager and sub-row HTML is baked at render time, so triage clicks show stale state until a full re-render that clears filter/sort/page. Read-only evidence; act in the modal |

---

## 7. Open questions

1. **Is the fallback workflow still supported?** `main_fallback.nf:43-44` is
   broken post-pivot and `assemble_fallback.R` writes `assemblies` through a
   path none of this traverses. Repair or formally mark `run_fallback()`
   unsupported - it should not stay silently broken by deferral.

   **ANSWERED 2026-07-30: supported, and it stops being a workflow.**
   `dev/design/14_ASSEMBLER_SWAP.md`. The fallback becomes a per-sample
   assembler inside WF1, selected by a new `assemble_opts.assembler` column
   (schema v16), and `main_fallback.nf` is deleted. `run_fallback()` keeps its
   signature and becomes a wrapper that repoints samples and runs the ordinary
   assembly workflow.

   Two corrections to the question as posed. The coverage arity bug it names
   was real but minor; the substantial defect was that `main_fallback.nf` never
   ran annotation, screening or triage at all, because it was a copy of the
   pre-pivot tail. And `assemble_fallback.R` does NOT write `assemblies`
   through a stray path - it writes no table at all, only FASTAs, which the
   shared `score` process loads exactly as it loads GetOrganelle's. Its headers
   were already the correct per-scaffold shape.

   The general lesson, recorded because it is the reason the bug survived: the
   duplication was the defect. Repairing the copy would have preserved the
   condition that produced it.
2. **Does `path_scores` survive?** Post-pivot it gates nothing: WF2 never reads
   it, `auto_select` drives only a dead button and two display columns. Either
   retire it with `score_candidates()`, or repurpose it as a *scaffold* scorer.
   Same question for the `stitched` table, `stitch_one`, `rp_stitch.R` and the
   stitch process - decide together, not piecemeal.

   **NARROWED 2026-07-30, not answered.** After `14_ASSEMBLER_SWAP.md` no
   shipped workflow calls the stitch machinery: `main_fallback.nf` was its last
   caller and is gone, and `query_fallback_samples` went with it. So the bundle
   is now unreferenced rather than half-referenced, which is what made it
   undecidable before. It is still not deleted: `stitch_scaffolds` and
   `stitch_one` are exported public API, so removing them is a breaking change
   with its own question.

   **ANSWERED 2026-07-31 for the stitch half: deleted.** `R/stitch_scaffolds.R`,
   `stitch_one()`, `stitch_selected()`, `inst/nextflow/bin/rp_stitch.R`, the
   `stitch` process, `annotate_project()` and `dev/prototype_stitch_viz.R` are
   gone, and schema v18 drops the `stitched` table from a fresh database and a
   migrated one alike, so the two cannot diverge. The three exports were removed
   as a breaking change on the grounds that none was reachable, none was tested,
   and `annotate_project()` had been broken since v7 anyway (it writes
   `annotations` without `path`/`scaffold`). See
   `plans/2026-07-31-stitch-removal.md`.

   `path_scores` is the exception and was NOT swept in; it is the live remainder
   of this question. It
   is still written by the live `score` process on every WF1 run
   (`R/pipeline_steps.R:37`) and read by the app (`R/app_server.R:121-122`,
   `:830`, `:842-857`). Only its stitch-side reader died. Deciding it is a
   question about the app's Selected Path / Best Score columns, not about
   stitching.
3. **How aggressive should auto-triage be out of the box?** The recommendation
   here is conservative: auto-drop only on positive cross-kingdom or
   prokaryotic evidence, everything else amber. Is that the right default for a
   200-sample fish project, where the extra amber rows cost real clicks?
4. **Is `_p<path>_s<scaffold>` the right SeqID suffix for GenBank?** It encodes
   assembly-graph coordinates that mean nothing to a downstream reader and will
   be interpreted as biological replicates. And what is the migration for
   projects that already exported bare IDs?
5. **Should the envelope numbers be shipped at all before validation?** They are
   the most consequential constants in the design and are a defensible prior,
   not a measurement. Alternative: ship the machinery with all rank-based rules
   disabled, and let the invert test clades plus a spike-in set the defaults.
6. **Is there run/lane/plate metadata available?** `samples.batch` makes the
   cross-sample screen meaningful and is meaningless without it.
7. **Should `paralogous_variant` block export, warn, or be purely
   informational?** For a clade with documented non-concerted rDNA evolution the
   curator arguably *should* report which copy was submitted, but that is a
   submission-policy call.
