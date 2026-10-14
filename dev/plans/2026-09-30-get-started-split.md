# Get started split (MitoPilot format) Implementation Plan

**Goal:** Split the single Get started page into a short Get started page plus
five test-project step pages, with a pill bar on top of each step page and
a Next button at the bottom, matching the MitoPilot docs site.

**Approach:** Text moves as is, by section. Button styles live once in
`pkgdown/extra.css` (MitoPilot pastes them into every page; we do not). The
five new pages get their own group in the Articles menu.

**Design:** approved in chat on 2026-09-30 (five pages, Next button only, no
Previous button, no restyle of existing callouts).

## Global constraints

- ASCII only in new text. No em dashes.
- Do not reword moved text, except the seam edits listed below.
- Figures stay in `vignettes/figures/get-started/`; paths do not change.
- `vignettes/RiboPilot.Rmd` keeps its file name (navbar "Get started" and
  README links depend on it) and keeps the `# Samples, scaffolds, and records`
  heading (Curation-Flags links to that anchor).
- Branch `docs/site-cleanup`. No push.

## Files

| File | Action |
|---|---|
| `pkgdown/extra.css` | create: `.step-nav`, `.next-step` styles |
| `_pkgdown.yml` | replace "Walkthrough" article group |
| `vignettes/RiboPilot.Rmd` | trim to lines 1-101, add hub pills and "Where to go next" |
| `vignettes/Test-Project-Taxonomy.Rmd` | create from old lines 105-133 |
| `vignettes/Test-Project-Reference.Rmd` | create from old lines 137-168 |
| `vignettes/Test-Project-Assemble.Rmd` | create from old lines 172-293, `## ` promoted to `# ` |
| `vignettes/Test-Project-Annotate.Rmd` | create from old lines 297-354 |
| `vignettes/Test-Project-Export.Rmd` | create from old lines 358-405 |
| `vignettes/Scripting.Rmd`, `Troubleshooting.Rmd`, `GenBank-Submission.Rmd` | re-aim one link each |

Line numbers refer to `git show HEAD:vignettes/RiboPilot.Rmd`.

## Task 1: styles and menu (cheap model)

- [x] Create `pkgdown/extra.css`:

```css
/* Walkthrough navigation: pill bar on top, Next button at the bottom. */
.step-nav {
  margin: 0 0 2.5em;
  padding-bottom: 1.2em;
  border-bottom: 1px solid #dde3e6;
}
.step-nav.no-rule {
  border-bottom: none;
  padding-bottom: 0;
}
/* pandoc wraps the row contents in a <p>, so that is the flex container. */
.step-nav p {
  margin: 0;
  display: flex;
  flex-wrap: wrap;
  align-items: center;
  gap: 8px;
}
.step-nav .step-label {
  font-weight: 600;
  color: #6b7b85;
  margin-right: 4px;
}
.step-nav a,
.step-nav .current {
  display: inline-block;
  padding: 6px 16px;
  border-radius: 6px;
  font-weight: 600;
  font-size: 0.95rem;
  text-decoration: none !important;
  border: 1px solid #c8dde5;
}
.step-nav a {
  background-color: #e4f0f4;
  color: #1f6f8b !important;
}
.step-nav a:hover,
.step-nav a:focus {
  background-color: #d3e7ee;
  color: #17566b !important;
}
.step-nav .current {
  background-color: #1f6f8b;
  color: #ffffff;
  border-color: #1f6f8b;
}
.next-step {
  text-align: right;
  margin: 2.5em 0 1em;
}
.next-step a {
  display: inline-block;
  background-color: #1f6f8b;
  color: #ffffff !important;
  font-weight: 600;
  font-size: 1.05rem;
  padding: 10px 22px;
  border-radius: 6px;
  text-decoration: none !important;
}
.next-step a:hover,
.next-step a:focus {
  background-color: #17566b;
}
```

- [x] In `_pkgdown.yml`, replace

```yaml
- title: "Walkthrough"
  navbar: ~
  contents:
  - 'RiboPilot'
```

with

```yaml
- title: "Test project walkthrough"
  navbar: Test project
  contents:
  - 'Test-Project-Taxonomy'
  - 'Test-Project-Reference'
  - 'Test-Project-Assemble'
  - 'Test-Project-Annotate'
  - 'Test-Project-Export'
```

## Task 2: split the page (cheap model)

Every step page has this shape:

```
---
title: "Test project: <Title>"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Test project: <Title>}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<pill bar>

<moved lines, verbatim>

<p class="next-step"><a href="<next file>.html">Next: <next label> &rarr;</a></p>
```

Pill bar, with the current page as `<span class="current">` and the rest as links:

```html
<div class="step-nav">
<span class="step-label">Test project:</span>
<a href="Test-Project-Taxonomy.html">1. Taxonomy</a>
<a href="Test-Project-Reference.html">2. Reference</a>
<a href="Test-Project-Assemble.html">3. Assemble</a>
<a href="Test-Project-Annotate.html">4. Annotate</a>
<a href="Test-Project-Export.html">5. Export</a>
</div>
```

| File | Title | Old lines | Next href | Next label |
|---|---|---|---|---|
| Test-Project-Taxonomy.Rmd | Taxonomy | 105-133 | Test-Project-Reference.html | Reference |
| Test-Project-Reference.Rmd | Reference | 137-168 | Test-Project-Assemble.html | Assemble |
| Test-Project-Assemble.Rmd | Assemble | 172-293 | Test-Project-Annotate.html | Annotation Review |
| Test-Project-Annotate.Rmd | Annotation Review | 297-354 | Test-Project-Export.html | Export |
| Test-Project-Export.Rmd | Export | 358-405 | Your-Own-Project.html | Your own project |

Seam edits (the only wording changes):

- Taxonomy page: between the pill bar and the moved text, add
  `This walkthrough picks up where [Get started](RiboPilot.html) leaves off, so start there if you have not yet created the test project and opened the app.`
- Assemble page: the four `## ` headings become `# ` (other RiboPilot pages use `# ` for sections).
- `RiboPilot.Rmd`: title and index entry become `Get started`. Keep lines
  9-101. In old lines 47-49, `the Reference step below builds` becomes
  `the [Reference step](Test-Project-Reference.html) builds`. Delete old lines
  102-412 and append:

```
# Test project walkthrough

The app has five steps. Each page below walks through one of them with the
test project. Do them in order:

<div class="step-nav no-rule">
<a href="Test-Project-Taxonomy.html">1. Taxonomy</a>
<a href="Test-Project-Reference.html">2. Reference</a>
<a href="Test-Project-Assemble.html">3. Assemble</a>
<a href="Test-Project-Annotate.html">4. Annotate</a>
<a href="Test-Project-Export.html">5. Export</a>
</div>

# Where to go next

After the walkthrough, [Your own project](Your-Own-Project.html) shows how to
write a mapping file and set up a project for your samples. If you work on a
computing cluster, read [Running on a cluster](Clusters.html); if something
goes wrong, see [Troubleshooting](Troubleshooting.html).
```

## Task 3: links, build, check (main session)

- [x] `Scripting.Rmd:30-31` and `Troubleshooting.Rmd:100-101`: "the rescue
  section of [Get started](RiboPilot.html)" becomes "the rescue section of the
  [Assemble walkthrough](Test-Project-Assemble.html#rescue-a-sample-that-failed)".
- [x] `GenBank-Submission.Rmd:13-14`: "[Get started](RiboPilot.html) shows the
  Export step itself." becomes "The [Export walkthrough](Test-Project-Export.html)
  shows the Export step itself."
- [x] Moved-text check: old lines 105-133, 137-168, 172-293, 297-354, 358-405
  appear verbatim in the new pages (Assemble after `## ` to `# `).
- [x] `Rscript -e 'pkgdown::check_pkgdown(); pkgdown::build_articles(lazy = FALSE)'`
  finishes clean.
- [x] Every relative `href` and `src` in the six built pages resolves to a
  file in `docs/articles/`, and every `#anchor` exists in its target page.
- [x] Screenshot Get started, one middle step page, and Export; eyeball pill
  bar and Next button.

## Review focus

- Pill bar markup depends on pandoc wrapping the row in `<p>`; confirm in built HTML.
- Pages with no section headings (Taxonomy, Reference, Annotate, Export) still render with a sane sidebar.
- `RiboPilot` is no longer listed under `articles:`; confirm pkgdown does not flag it as missing.
