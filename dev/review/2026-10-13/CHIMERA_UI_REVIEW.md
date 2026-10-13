# Chimera check UI/UX review

Date: 2026-10-13. Branch: devel (cfad2a3). Read-only review, no code changed.

Four reviewers, each with a different lens:

1. Senior data visualization designer (colorblind accessibility, measured CIE dE)
2. First-time user (second-year bio grad student, knows BLAST and Geneious)
3. Principal product designer / Shiny architect (genome-browser tools)
4. Senior R/Shiny engineer and QA lead (trimmed-contig behaviour)

Key claims below were spot-checked against the code by the coordinator.

## What the panel is today

`chimera_panel()` (R/app_modals.R:149) renders, in both the Assemble details modal and the Annotation Review modal, collapsed by default:

- a pill ("Chimera call: severe") with the full reason hidden in a hover title
- a help paragraph and a "windows searched against: ..." line listing internal DB names
- a static 3-facet ggplot (`plot_chimera()`, R/app_modals.R:1041): gene boxes, a lineage segment bar, a % identity dot track, dashed junction lines with labels
- three tables: segments (rows tinted as the legend), possible junctions, all windows (collapsed)

## Consensus diagnosis

1. **The plot answers "how good is each BLAST hit" when the user asks "which part is foreign and where does it break".** The tall facet is % identity; the foreign-vs-own story sits in a thin pastel bar with no legend. In CHIM01 the Penicillium half is the flattest, highest line, so it reads as the "good" half.
2. **Colors fail.** Two of the five lineage colors are indistinguishable from each other (`#BBCC33` vs `#AAAA00`, dE 9.6) and one is indistinguishable from the "own lineage" grey under deuteranopia (`#44BB99` vs `#B3B3B3`, dE 4.8). Overlapping segments (by design, 200 bp windows at step 100) stack alpha 0.6 fills into a darker patch that looks like a fourth category, right where the junction sits.
3. **The verdict is hover-only.** The reason string is the one sentence users most need; it is a `title` attribute on the pill and on the table badge. Nothing inline says "positions 1-2500 are your fish, 2501-4371 are a mold, cut at 2625".
4. **It is a second, static genome browser.** Genes, coverage, and clipped-read bars already live in the zoomable seqviewer canvas directly above, and the seqviewer already uses the chimera check's clip cutoff. The chimera plot has no depth track, so the "reads" channel evidence can only be cross-read between two panels on different x scales.
5. **The action is far from the evidence.** Trim lives in a different window (Annotation Review > Trim assembly) and the only chimera guidance there is a bare text line "suggested junctions: 2625", not drawn on the brush plot.

## Trimmed contigs: verdict is mostly "fine", with one surprise

The engineer's core finding: **the app never shifts chimera results on trim.** `review_trim()` calls `review_invalidate()` (R/review.R:109), which deletes the review-layer chimera and window rows and sets `check_switch = 1`. Recheck then re-windows the trimmed sequence from position 1 (R/chimera.R:462). `chimera_shift()` and `chimera_restore()` are only called by the schema-4 upgrade (R/review_upgrade.R:147-153). App Restore is `review_drop()`, which deletes the copy so views fall back to originals.

Consequences:

- After trim, before Update: no stale or mis-shifted data can reach the panel. `chim_unit()` (R/app_server.R:2685) and `plot_chimera()` (R/app_modals.R:1046) gate on the same `review_unit_state()`. Tests cover it (test-chimera.R:626-666).
- After Update: x axis, windows, segments, seqviewer, and trim plot all use trimmed coordinates. No mismatch found.
- Trim-modal hint "suggested junctions: N": N is always in the coordinates the user sees (verified for untrimmed, trimmed+rechecked). Hidden when the check is not done, which is correct but silent.
- Annotation Review table badge and Checks column read the review layer and show "not run" after a trim. Triage reads originals by design.

Gaps:

| # | Kind | Where | Finding |
|---|---|---|---|
| T1 | UX gap | R/app_server.R:2704-2708, 2785-2789 | Panel goes almost blank after a trim. Only the plot says "Not run yet. Press Update."; pill, segments, junctions, windows all `req()` out. Curation warnings panel prints a "Checks not run" line in the same state; chimera should too. |
| T2 | UX gap | R/app_server.R:3156-3163, R/app_modals.R:988-1003 | Trim modal: junction is a number, not a line on the brush plot. After a trim and before Update the hint vanishes with no message. |
| T3 | UX gap, minor | R/app_modals.R:149-170 | Panel never says which layer it shows (Assemble = original, Annotation Review = review copy). Docs and trim help say so; the panel does not. |
| T4 | Latent bug | R/chimera.R:589-599 | `chimera_shift()` does not shift `channels`, so after a shift over 300 bp the junction table loses its Checks and Why columns, and `reason` is rebuilt without ribotyper/reads text. Only reachable via the upgrade fallback for old projects without a snapshot. |
| T5 | Latent bug | R/chimera.R:588, 599 | `chimera_shift()` downgrades label to "clean" when all junctions fall outside the kept window even if two foreign-class segments remain inside it. Dead in production (upgrade path never drops anything). |
| T6 | Doc nit | R/chimera.R:570-574 | Comment says "after its sequence was trimmed or restored"; trims never call it. Say "upgrade shift-back". |

## Visual findings (ranked)

| # | Sev | Finding |
|---|---|---|
| V1 | high | `#44BB99` (5th foreign colour) vs `#B3B3B3` (own lineage) at alpha 0.6: dE 4.8 deutan, 7.2 protan. Fifth foreign taxon reads as "own". |
| V2 | high | `#BBCC33` vs `#AAAA00`: dE 9.6 for everyone. 3rd and 4th foreign taxa merge. |
| V3 | high | Overlapping segments stack alpha into a darker patch (CHIM01, 2500-2750) that looks like a separate category; junction line lands inside the muddy band. |
| V4 | high | % identity is the dominant visual; foreign-vs-own is the question. Semantics inverted. |
| V5 | med | No in-plot legend. Segments table doubles as legend but is rarely on screen with the bar (300 px plot plus help text). "lineage" strip label does not say lineage of what. |
| V6 | med | Gene colors (ITS1 `#E69F00`, 28S `#D55E00` darkened at alpha 0.28) sit directly above yellow `#EEDD88` at 0.6 (dE 19). Same hue family. |
| V7 | med | Junction labels are two lines, ~150 px wide, all at one y. Two junctions closer than ~700 bp overlap. |
| V8 | low | Red vs amber junction line: dE 17 deutan, but "chimera"/"possible" text and the pill carry severity too. Fine. |
| V9 | low | No-hit dots at y = 68 sit on the axis under the 70 tick, half in the margin; dotted range lines clip at the panel edge. Looks like a render glitch. |

## First-time user findings (ranked by what would help most)

1. No plain verdict. Wants: "Positions 1-2500 match ray-finned fish (your sample). Positions 2501-4371 match Penicillium (a mold). Junction near 2625 bp, supported by gene hits, windows, and a read-depth step."
2. Cannot tell from the plot which side is foreign. Nothing names a taxon on the plot; nothing marks "your lineage".
3. Jargon with no in-panel definition: ribotyper, clipping, windows, unresolved, nohit, `high_id`, `low_id` (R/chimera.R:190), rank dist, bitscore, clip frac, "channels" (genes, windows, reads).
4. "minor" with "No possible junction was found" and an all-grey plot reads as contradictory. The explanation (genes agree only at phylum) is hover-only.
5. No path from evidence to action. Wants "Trim here" or at least a pointer to Trim assembly / Keep-Drop / Mark reviewed with the junction position carried across.
6. "windows searched against: ribopilot_rrna, 16S_ribosomal_RNA, ..." teaches nothing. Hide or say "bundled rDNA database".
7. Reason string for nohit segments prints "(, NA%)" (`chimera_reason()`, R/chimera.R:402-403: `species %||% ""` and `%.1f%%` of NA). Looks broken.
8. Docs say 200 bp windows; the CHIM01 screenshot (vignettes/figures/get-started/assemble-chimera-chim01.png) says 500 bp. Screenshot is stale.
9. Panel is buried: 3rd collapsed section in Assemble, 6th in Annotation Review. Nothing on the table row points at it.

## Redundancy map (architect)

| Fact | Shown where | Verdict |
|---|---|---|
| Label | Annotate table badge, Assemble scaffold table, Checks column, curation warnings, panel pill | Fine; badges are the summary layer |
| Reason | Hover on badge (x2) and pill | Bad; junction table "why" is the readable form, surface that |
| Junctions | Plot vlines + labels, junction table, reason, trim hint text | 4 copies; seqviewer (the only zoomable view) has none |
| Lineage segments | Plot bar, segments table tint, windows table tint, reason | Bar + table are one unit; bar belongs in the seqviewer |
| Per-gene origin | Feature BLAST table vs chimera "genes" channel | Near duplicate; cross-reference, do not merge |
| Gene boxes | Chimera plot facet, seqviewer, feature table, badges | Plot facet is pure duplicate; drop |
| Depth and clips | Seqviewer, trim plot, sparkline | Chimera plot has none, yet "reads" channel is about depth steps and clip pileups |
| db_source line | Panel | Low value; move to tooltip |

## Recommended target design

```
Annotation Review - SRR... [path 1 / contig 1]
+---------------------------------------------------------------------------+
| [severe] Chimera. BLAST windows switch lineage at 2625 bp (Actinopteri -> |
| Eurotiomycetes). Genes' hits disagree: 28S = Penicillium. 61% of reads    |
| clipped at 2625.         [Trim at 2625...] [Mark reviewed]                |
+---------------------------------------------------------------------------+
Feature model table                       [Trim assembly...] [Restore]
v Sequence and read coverage
  [goto][Fit][Whole][+][-][Copy] [x]Coverage [x]Clipped [x]Origin [x]Nucleotides
  ruler    500    1000    1500    2000   |2625    3000    3500    4000
  depth    ..######################|______________________
  clipped  ____________________________|______________________
  origin   [no hit][ Actinopteri ][no hit][Actin]|[ Eurotiomycetes       ]
  lanes            [18S >>>>>>>>>>>>>>>>>]|[28S >>>>>>>>>>>>>>>>>>>>>>]
                                  junction 2625 (genes, windows, reads)
> Edit active feature boundaries
> Reference alignment
> Feature BLAST
> Curation warnings
> Chimera evidence
    segments table (row expands to its windows); possible junctions table
> Notes

Trim modal:
Trim to: (o) Annotated  ( ) Annotated + flank  ( ) Keep 1..2624  ( ) Keep 2625..4371  ( ) Custom
[coverage plot: grey trimmed flanks, dashed from/to, red dotted junction line labelled 2625]
```

Assemble modal: same card with "[Drop contig]" in place of trim, same seqviewer track, same evidence section. Both modals keep the same evidence; only the action buttons differ. Verdict card is open by default when label != clean; the evidence section stays collapsed.

## Phased path (smallest valuable step first)

### Phase 1, small, no new UI surfaces (do first regardless of later phases)

- Palette: cut foreign colors to three with pairwise dE >= 26 under deutan/protan and >= 27 from grey: `#EEDD88`, `#77AADD`, `#AAAA00`. Drop `#BBCC33`, `#44BB99`, `#99DDFF`. (V1, V2)
- Real ggplot legend at the top of the plot, keyed by segment taxon ("Actinopteri (your lineage)", "Eurotiomycetes", "no hit", "unresolved"); bump plot height to ~330 px. In-bar taxon names were rejected: narrow and overlapping segments make them unreadable. (V5, user 2)
- Clip overlapping segments for drawing only (`stop[i] <- start[i+1] - 1`); tables keep true coordinates. (V3)
- Label only foreign segments, and only when wider than ~10% of the contig after de-overlap; own-lineage and no-hit segments never get text. Optional: a `#` column in the segments table with the digit drawn at the bar midpoint. (user 2)
- Draw fills at alpha 1 via the existing `tint_colors()` rather than alpha 0.6, which also removes the stacking. (V3)
- Foreign dots bigger with a dark ring, own-lineage dots small and grey, no-hit hollow; optionally shade the identity panel background per segment at alpha 0.12. (V4)
- Junction label one line ("junction 2625"), stagger vjust by index. (V7)
- Show `reason` inline under the pill (or better, the junction table "why" sentences), drop the "hover for the full reason" line. (user 1, architect)
- Fix "(, NA%)" for nohit segments in `chimera_reason()`. (user 7)
- After trim, before Update: pill renders a grey "not run" with the same text the curation panel uses; junctions section prints a line instead of nothing. (T1)
- Replace the db_source line with "searched against the bundled rDNA database", details in a tooltip. (user 6)
- Re-shoot `assemble-chimera-chim01.png` (says 500 bp). (user 8)

### Phase 2, small: trim modal

- `plot_trim()` draws dashed junction vlines, labelled. (T2)
- `trim_modal()` gains "Keep 1..J-1" / "Keep J..n" radio presets per junction; `trim_preset()` handles them. Delete `trim_chimera_hint`.
- When the check is not done, say "junction hint unavailable until Update" instead of hiding.

### Phase 3, small-medium: verdict card and actions

- Replace `chimera_result_text()` with a card at the top of the modal body (next to `failure_banner`), only when label != clean: pill, the junction "why" sentences, buttons. Annotation Review: "Trim at J..." (opens trim modal with the preset selected), "Mark reviewed". Assemble: "Drop contig".
- Checks column text gains the position ("severe chimera at 2625").
- One help line: "Annotation Review shows the review copy; Assemble shows the original." (T3)
- In-panel glossary for ribotyper, windows, unresolved, no hit; hide or rename `high_id`/`low_id`; drop rank dist / bitscore / clip frac from the default windows columns. (user 3)

### Phase 4, medium: origin track in the seqviewer, retire the static plot

- `seqview_payload()` gains `segments` [{pos1, pos2, color, taxon, species, pident}] and `junctions` [{pos, severe, checks}]; `seqview_send()` pulls them with the queries `plot_chimera()` uses today.
- seqviewer.js: one more lane (~16 px) drawn with the same column loop as `drawTrack`, dashed junction line across all tracks with a ruler label, tooltip reports segment taxon at the base, "Origin" checkbox next to Coverage / Clipped. Roughly 70 lines of JS, 20 of R.
- Open the seqviewer section in the Assemble modal when label != clean (`examine_unit` observer already manages collapse).
- Retire `plot_chimera()`, `output$chimera_plot`, and the `fixed_gutter()` call. Fold "All windows" into the segments table as reactable `details`. `chimera_panel()` shrinks to segments + junctions. Add a % identity track only if curators miss the dots.

Reused: `chimera_seg_colors()`, `CHIMERA_LINEAGE_COLORS`, `chimera_junction_table()`, `output$chimera_segments`, `output$chimera_junctions`, `plot_trim()`, `plot_coverage()`, `trim_preset()`, `seqview_*`, `drawTrack`. Removed: `plot_chimera()`, `output$chimera_plot`, hover-only reason, `trim_chimera_hint`, db_source line, All windows subsection.

### Separately: latent bugs in `chimera_shift()` (T4, T5, T6)

Only reachable through the schema-4 upgrade fallback. Either fix (shift `channels` junctions too, keep extra reasons, do not auto-downgrade to clean) or leave and fix the doc comment. Low priority.

## Open decisions for the user

1. Phase 4 removes the static plot entirely. Keep it as an optional "static view" for screenshots in docs/papers, or drop it?
2. Verdict card at the top of the modal (always visible) vs inside the collapsed section (status quo location, just inline instead of hover).
3. Trim presets per junction: radio options in the trim modal, or a "Trim at J" button on the verdict card that pre-fills Custom range?
4. Palette: three foreign colors is enough if segments are labeled; confirm that four or more foreign taxa on one contig is rare enough to accept color recycling.

## Addendum: merge into the seqviewer instead of polishing the ggplot (user direction)

Decision direction after discussion: skip the static-plot polish and move the chimera evidence into the canvas sequence viewer (inst/app/www/seqviewer.js, 394 lines). Reasons:

- Hover gives the taxon at any base and zoom makes narrow segments wide, so no in-bar labels are needed.
- Depth, clips, lineage, genes, and junctions share one x axis; today the "reads" evidence and the "windows" evidence sit in two panels with different scales.
- A standalone JS chimera plot would duplicate the ruler, pan, zoom, tooltip, and resize code the viewer already has.
- `plot_chimera()` has no caller outside `output$chimera_plot`; the export report does not use it. Its four tests move to `seqview_payload()` tests.

Design:

- Origin track: flat band under the clipped-read track, per-pixel fill from the segment color (own = grey, no hit = white, foreign = 3-color palette), de-overlapped in R before sending.
- Junctions: dashed vertical line through all tracks, red (severe) or amber, labelled in the ruler; tooltip near it lists the checks.
- Optional identity track: windows expanded to per-base max identity, drawn with the existing `drawTrack` (vmax 100); no-hit = gap.
- "Origin" checkbox beside Coverage / Clipped / Nucleotides. Tooltip gains "origin: <taxon> (<best match>, <pident>)".
- R: `seqview_payload()` gains `segments`, `junctions`, optional `ident`; `seqview_send()` reuses `plot_chimera()`'s queries, honoring the review flag and the check state (sends nothing when not run). ~30 lines. JS ~70 lines.
- Retire `plot_chimera()`, the `fixed_gutter()` call, `output$chimera_plot`. `chimera_panel()` shrinks to pill, inline verdict, segments table with windows as expandable rows, junctions table.

Revised order:

1. Origin and junction tracks in the seqviewer; retire the static plot. Palette fix and draw-only de-overlap land here.
2. Verdict inline under the pill; "not run" text after a trim; fold windows into the segments table; hide raw DB names.
3. Trim modal: junction lines on the brush plot and keep-left / keep-right presets. The trim plot stays ggplot (brush).

Caveat: the Assemble modal opens every section collapsed, so a non-clean call needs the seqviewer section auto-opened or a verdict card above the sections.

### Labels and tables after the merge

Origin segment labels follow the gene-lane rule already in seqviewer.js (line ~208): draw the taxon name only when the visible part of the segment is wider than `measureText(label).width + 8`; otherwise nothing, hover carries it. Zoom makes narrow segments labelled. Tiered label when room allows: "Actinopteri" first, then "Actinopteri, <best match> <pident>%". Own-lineage segments are labelled "(your lineage)". Segments are de-overlapped in R before sending, so labels never collide.

Tables keep their job; only the plot goes:

- Segments: main table, gains a `#` column, windows fold in as expandable `details` rows, row tint stays as the track legend.
- Possible junctions: kept; its "why" sentences also become the inline verdict text under the pill.
- All windows subsection: dropped (folded into segments).

Optional: clicking a segment row zooms the viewer to it (`fit()`/`goto()` via one custom message), mirroring feature-row selection.
