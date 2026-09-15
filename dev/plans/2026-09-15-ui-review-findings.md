# RiboPilot UI/UX review, 2026-09-15

Four independent passes: MitoPilot gap sweep (Sonnet), novice-user wording
(Haiku), UX designer layout (Sonnet), senior Shiny dev design clarity (Opus).
Consolidated and de-duplicated. Line numbers as of commit daf8618.

## MitoPilot gap sweep

No generic UI/UX change in MitoPilot since 2026-08-01 is missing from
RiboPilot. Everything else in that window is map-to-ref, seqview, synteny,
circularize, or outlier review, which have no RiboPilot equivalent.

## High: behaviour and safety

1. app_ui.R:49-89: state_btn, lock_btn, run_btn emitted once per mode inside
   conditionalPanel. All copies stay in the DOM with one id; each keeps its own
   click counter, so the first click in a second mode can resend an unchanged
   value and be dropped. Verified: Shiny binds all copies and only warns.
   Fix: per-mode ids (tax_lock_btn ...) like filter_pickers(), branch in server.
2. app_server.R:944, app_reference.R:978: "Close" on the run and reference
   build modals calls kill_nf(). Label says close, action kills the pipeline.
   Fix: confirm when a process is alive, or close without killing.
3. app_server.R:1940-2135: boundary nudges write to the DB immediately. Close,
   Lock & Close, and sample nav skip the dirty check, so edits persist with
   stale flags and no Reset. Fix: route all exits through edit_dirty(), or
   buffer nudges until Save.
4. app_server.R:678: bulk Lock is a per-row toggle. Mixed selection locks some,
   unlocks others, one toast, and unlocked rows vanish from the next tab.
   Fix: one direction per click, confirm unlocks, say what leaves the next tab
   (MitoPilot write_lock pattern).
5. app_server.R:231 vs 402: Export table ID cell is still an open_detail link
   but open_detail ignores Export mode. Dead link. Fix: plain text in Export.
6. app_server.R:88: status line counts visible_mask() only; reactable column
   filters and search are not counted. "Showing 50 of 50" over three rows.
   Fix: read reactable pageInfo, or drop the visible count.

## Medium: state, flows, consistency

7. app_server.R:212: column defs (multi_unit, ID width, ref_details) frozen at
   render; updateReactable pushes new rows. Path/Scaffold stay hidden after a
   sample becomes multi-unit until mode switch.
8. app_server.R:950, 1028, app_modals.R:593: export counts are selection
   narrowed but labelled as "locked & ready"; modal count never updates when
   the group picker changes.
9. app_server.R:703: toolbar Reviewed only sets yes; Problem/Partial toggle;
   modal Mark/Clear toggles. Three behaviours for one flag. Use mp_flag_next()
   for all three toolbar buttons.
10. app_server.R:1497, 1095: "Drop a reviewed scaffold?" and export overwrite
    cancels close the parent modal too. Reopen parent on cancel (trim_cancel
    already does this).
11. app_opts.R:428: options modal silently drops locked samples from a mixed
    selection and reports the smaller count. Name the skipped samples.
12. app_guards.R:1031, utils_reactable.R:242: help says "Unlock in the Lock
    column" but the lock cell is a static glyph. Name the Lock button, or make
    the cell a toggle.
13. app_reference.R:160: Move radio defaults to first clade; pressing Move
    untouched reassigns silently. selected = character(0), refuse until picked.
14. app_reference.R:458: proposal Apply overwrites hand-set memberships with a
    plain primary button. Danger tone, count in label, or exclude hand-set.
15. app_server.R:717: empty group name closes the dialog and does nothing.
16. app_server.R:2098: notes autosave per keystroke with no feedback; clearing
    text does not clear the stored note. Say so under the box.
17. app_server.R:737: CSS row hiding leaves rows in the reactable model; pages
    can render empty. Filter server-side and updateReactable.
18. app_modals.R:349: three exit styles across assemble (easyClose), annotate
    (Close button), trim (neither); none has the title X. One rule.
19. app_server.R:691: Annotate state modal reads state from ids[1] without
    path/scaffold; preselect from an arbitrary unit.
20. app_server.R:2082: Restore assembly discards trim without confirm.
21. app_notify.R vs app_alerts.R: rp_alert/rp_toast and mp_alert/mp_toast
    both live and mix in one flow. Keep one pair.
22. app_modals.R:16-19, 104-105, 679-681, utils_reactable.R (ten cell
    renderers): three or four hand-typed hex palettes for success/warning/
    danger that drift from the --mp-* tokens (#842029 vs #c0392b). One shared
    palette, or mp-pill classes.
23. app_modals.R (20 inline font-size declarations), app_ui_helpers.R:216:
    four ad hoc help-text stylings beside opts_help() / p.mp-help-text.
24. app_modals.R:562, 442-455, 500-513: footers mix default-size and btn-sm
    buttons in one row. Bake one size into mp_footer().
25. app_reference.R:300-334: reference action modal stacks four full-width
    body buttons with no clear primary, unlike every other modal.
26. app_modals.R:317 vs 362-403: assemble panels all closed by default,
    annotate panels mostly open. Pick one policy.
27. app_server.R:184-204: Annotate table has ~20 columns and no column-group
    picker; Assemble has one.
28. app_ui.R:108-110: CSV download buttons lack mp-toolbar-btn, so a different
    height from neighbours. app_ui_helpers.R:83 mp_csv_download_row unused.

## Low

29. app_ui.R:22 ctrls() unused. app_ui.R inline flex styles instead of classes.
30. app_modals.R:456 / app_reference.R:514 duplicated run-modal title blocks.
31. app_opts.R:95,133,147 same border box style three times.
32. app_modals.R:314 nudge buttons have no aria-label; app_server.R:475 nav
    arrows bypass mp_toolbar_button (no aria-label).
33. utils_reactable.R:254 rt_link sets id=index+1, repeated across columns.
34. app_server.R:1367 unit_candidates and chimera outputs lack refresh().
35. app_server.R:1521 triage fallback branch unreachable.
36. app_modals.R:32-60 seed-hit caveat only visible inside a collapsed panel.
37. app_reference.R:170 rung labels do not show whether a count is inside the
    min/max window set elsewhere.

## Wording (novice pass, keep the ones that matter)

- "Problem" button vs "Problematic" column and pill (app_ui.R:92).
- "Annotate" tab vs "Annotation Review" step and modal title (app_ui.R:42, 87).
- WF1 (app_modals.R:122), ITSx (constants.R:70,135), taxid (app_modals.R:64),
  triage (constants.R:122), unit/path/scaffold (constants.R:120, 127) used
  without a first-use definition. One short gloss each, or a glossary popover.
- "No unlocked samples are queued" (app_server.R:851) says what, not how.
- Chimera blurb (app_modals.R:121) is one 38-word sentence. Split.
- "Opts." column headers (constants.R:274-280): spell out Options.
- Lock button tooltip wording differs per mode for the same action
  (app_ui.R:83-95).

## Status (applied 2026-09-15 on fix/ui-review-2026-09)

Done: 1-16, 19-24, 26-28, 30-35, 37, and the wording items (Problem ->
Problematic, WF1 -> "the pipeline (Assemble > Update)", ITSx/triage/unit
glosses, queued-samples message says how, chimera blurb split, Opts. ->
Options, Reviewed tooltip says mark or clear).

Notes on how:
- 1: tax_/asmb_/annot_ prefixed toolbar ids; one observer over the list.
- 4: any locked row in the selection means the click unlocks all of it, else
  locks all; unlock confirmed; toast says locked/unlocked.
- 6: custom.js reports the rows reactable shows (after column filters, search,
  and the lock/state CSS) as input$tbl_visible_n; the server mask is the
  fallback before the client reports.
- 21: rp_alert/rp_toast kept as names, now thin wrappers over mp_alert/mp_toast.
- 22: table badges use the mp-pill classes; remaining hex in
  utils_reactable.R are component colours (keep/drop buttons, depth bar,
  match-type dots), not the success/warning/danger tones.
- 35: fallback kept; the scaffold-uid tests open Triage with no selection.

Skipped:
- 17 (server-side row filtering): CSS hiding is what keeps the selection and
  column filters alive across lock/state writes; 6 covers the visible symptom.
- 18 (one dismissal rule): annotate Close is now a guarded server button;
  assemble stays easyClose since it has nothing to lose. Title X not added.
- 25 (reference action modal footer): it already has one primary; left as a
  choice list.
- 29 (inline flex styles in app_ui.R): cosmetic; ctrls() removed.
- 36 (seed caveat in a collapsed panel): the panel is where the hit is read.
