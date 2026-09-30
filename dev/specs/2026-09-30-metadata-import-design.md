# Sample metadata import from GEOME, GBIF, and NCBI

Date: 2026-09-30
Status: approved (design agreed in chat; full port, three stages)

## Goal

Port MitoPilot's specimen metadata feature into RiboPilot using the same
method: identifiers in the mapping file, per-source fetch chains over public
APIs, every fetched field stored as is, GenBank-ready export tokens, a viewer
in the app, link-following between databases, and conflict detection against
the mapping file. Where MitoPilot code can be copied as is, copy it; where
RiboPilot's structure differs, adapt the wiring, not the method.

Source of truth for behaviour: MitoPilot v1.5.8 at
`/home/dmacguig/Documents/GitHub/MitoPilot` (files listed under "Files to
port"). Its guide `vignettes/Specimen-Metadata.Rmd` describes the user-facing
behaviour this spec expects, except where "Deviations" says otherwise.

## Scope by stage

| Stage | Delivers |
|---|---|
| 1 | Identifier columns, normalizers, fetch chains, storage tables, schema step, console functions (`fetch_geome`, `fetch_gbif`, `fetch_biosample`, `remove_metadata`), setup-time fetch in `new_project()`, `add_samples()`, `update_sample_metadata()`, the viewer modal with Fetch, Refresh all, and Remove, export field selection ("Set Export Metadata"), export columns and token chips, tests with recorded fixtures. |
| 2 | Link-following (`link_sources`, `meta_links`), the Compare tab, `set_metadata_columns()`, the conflict warning before export. |
| 3 | The "Metadata" table-column picker (`meta_view_fields`), the guide page, docs updates, NEWS. |

Each stage ends with the full test suite passing and is committed on the
branch `feat/metadata-import` before the next starts.

## Data model

### `samples` columns

Three TEXT columns, added in `new_db()` and by the schema step: `GEOME_BCID`,
`GBIF_ID`, `BioSample`. A mapping file may carry them under these names or under
custom names given by `mapping_geome`, `mapping_gbif`, `mapping_biosample`
(default `"GEOME_BCID"`, `"GBIF_ID"`, `"BioSample"`). Values are stored
normalized; an invalid value is stored raw so it stays visible and its fetch
fails with a clear message. These three names are not added to
`SAMPLE_BUILTIN_COLS`: the mapping file may supply them, and `read_mapping()`
must accept them.

### Tables (verbatim from MitoPilot `meta_db.R` and `meta_view.R`)

```sql
CREATE TABLE IF NOT EXISTS meta_records (
    ID TEXT NOT NULL, source TEXT NOT NULL, level TEXT NOT NULL, depth INTEGER NOT NULL,
    ref TEXT, field TEXT NOT NULL, value TEXT, PRIMARY KEY (ID, source, depth, field));
CREATE TABLE IF NOT EXISTS meta_status (
    ID TEXT NOT NULL, source TEXT NOT NULL, ref TEXT, status TEXT NOT NULL, message TEXT,
    fetched_at INTEGER, PRIMARY KEY (ID, source));
CREATE TABLE IF NOT EXISTS meta_export_fields (key TEXT NOT NULL, PRIMARY KEY (key));
CREATE TABLE IF NOT EXISTS meta_csv_map (concept TEXT NOT NULL, column TEXT, PRIMARY KEY (concept));
CREATE TABLE IF NOT EXISTS meta_options (key TEXT NOT NULL, value TEXT, PRIMARY KEY (key));
CREATE TABLE IF NOT EXISTS meta_links (
    ID TEXT NOT NULL, source TEXT NOT NULL, ref TEXT, via TEXT, note TEXT, PRIMARY KEY (ID, source));
CREATE TABLE IF NOT EXISTS meta_view_fields (key TEXT NOT NULL, shown INTEGER NOT NULL, PRIMARY KEY (key));
```

All seven tables and the three columns are created in stage 1 (one schema
step, even though stages 2 and 3 fill some of them), so the schema changes
once.

### Schema step

- DDL constants live in `R/db_schema.R` (`META_RECORDS_DDL`, ... one per
  table), following the existing `CREATE TABLE IF NOT EXISTS` idiom.
- `new_db()` executes them next to `ensure_export_opts(con)` and adds the three
  `samples` columns before the `created` history row.
- `db_migrations()` gets its first step (version 1 to 2):
  `list(note = "sample metadata (GEOME/GBIF/NCBI) tables and columns", fn = ...)`
  which creates the tables and adds the columns when missing.
- `tests/testthat/fixtures/schema_v1.sql` is untouched. The parity test in
  `test-backwards-compat.R` must pass. `fixture_fill_schema()` in
  `helper-app-fixture.R` learns the new DDL constants and the new columns.
- `.meta_ensure_tables()` is kept as an internal guard (idempotent) so console
  functions work on a project that was just upgraded, but no code relies on lazy
  creation. The MitoPilot legacy `geome_*` migration inside it is dropped.

### Storage semantics (unchanged from MitoPilot)

- One row per (sample, source, depth, field) in `meta_records`; scalar fields
  only; numbers formatted with 15 digits, never scientific.
- `meta_status` holds one row per (sample, source): `ok` or `failed`, message,
  `fetched_at`.
- A failed refetch of the same identifier keeps the previous records.
- Changing or clearing an identifier drops that source's records and status for
  the sample at once.
- Nothing is merged across sources.

## Sources and fetch chains

`META_SOURCES` registry (`meta_db.R`) with `col`, `label`, `id_label`, `arg`,
`normalize`, `invalid`, `chain`. Ported verbatim:

- GEOME: `records/<bcid>?includeParent=true` up the parent chain (loop and
  depth guard 10), then `projects/<pid>/expeditions/<code>`, then
  `projects?includePublic=true` filtered client-side.
- GBIF: `occurrence/<key>`, `dataset/<datasetKey>`, `organization/<orgKey>`;
  NMNH EZIDs resolved through `occurrence/search`.
- NCBI: optional SRA (`esearch` + `esummary`), then BioSample `efetch` XML, then
  BioProject `efetch` XML.

HTTP: `httr2` with user agent `RiboPilot (https://github.com/dmacguigan/RiboPilot)`,
30 s timeout, 3 tries on 429/500/502/503/504, status-to-message mapping as in
MitoPilot. Three low-level functions `.gbif_get`, `.geome_get`, `.ncbi_get` are
the only places that touch the network, so tests mock them.

In-memory cache per batch (environment) so shared ancestors, datasets,
organizations, and bioprojects are fetched once.

## Console API

Exported, added to the hand-maintained `NAMESPACE`:

```r
fetch_geome(path = ".", ids = NULL, bcids = NULL, link_sources = NULL)
fetch_gbif(path = ".", ids = NULL, gbifs = NULL, link_sources = NULL)
fetch_biosample(path = ".", ids = NULL, biosamples = NULL, link_sources = NULL)
remove_metadata(path = ".", sources = names(META_SOURCES), ids = NULL)
set_metadata_columns(path = ".", ...)            # stage 2
geome_normalize_bcid(x); gbif_normalize_id(x); ncbi_normalize_id(x)
```

Setup entry points gain the MitoPilot arguments, with the same defaults:
`new_db()`, `new_project()`, `new_test_project()` (pass-through), `add_samples()`,
and `update_sample_metadata()` get `mapping_geome`, `fetch_geome`, `mapping_gbif`,
`fetch_gbif`, `mapping_biosample`, `fetch_biosample`, `link_sources`.
`update_sample_metadata()` keeps RiboPilot's existing behaviour (Taxon change
re-queues taxonomy) and adds `.meta_take_cols` and `.meta_sync_changed`:
refetch only changed IDs, drop data for cleared IDs, skip unchanged.

Mapping checks: `read_mapping()` accepts the three columns; `check_mapping()`
warns about invalid identifiers ("these samples will show a failed GEOME
fetch"); `new_project()` does the reachability pre-flight only when a column
holds at least one valid ID and fetching is on.

## Export

- `meta_export_fields` holds the ticked keys (`<source>:combo:<name>` or
  `<source>:raw:<Level>:<field>`); `.meta_key_col()` turns a key into a column
  name such as `geome_lat_lon` or `gbif_Occurrence_occurrenceID`.
- Combos ported verbatim from `geome_export.R`, `gbif_export.R`,
  `ncbi_export.R`: `lat_lon`, `collection_date`, `geo_loc_name`,
  `specimen_voucher`, `collected_by`, `identified_by`, `sex`, `dev_stage`,
  `tissue_type`, plus the plain NCBI `biosample` and `bioproject`. Nearest level
  wins within a source; no cross-source fallback.
- `.meta_join(dat, con)` adds the ticked columns to the per-unit header data in
  `export_header_data()` (`R/export_opts.R`) and to the export summary rows in
  `R/export.R`, so they work as `{token}`s in header templates and appear in
  `<group>_summary.csv`. A ticked column that collides with a computed token
  name is shadowed by the computed one, as today.
- "Set Export Metadata" toolbar button in the Export step opens the fields modal
  (combos as checkboxes with example and sample count, raw fields in a
  searchable table).
- The Export window's "Available columns" block is replaced by the token chip
  groups from `app_export_tokens.R`: Basics, Your mapping columns, GEOME, GBIF,
  NCBI, and RiboPilot columns. A combo chip inserts `[label={col}]`, a plain
  chip inserts `{col}`, chips with empty values in the chosen group turn
  orange, and "Choose fields" per source opens the fields modal and returns.
- Empty values leave empty modifiers; RiboPilot's existing `empty_modifiers()`
  warning covers them.

## App

RiboPilot's app is one table with modes. Wiring:

- A **Metadata** column (`specimen_col_def`, logos per linked source, orange
  badge on failure or conflict) in the Taxonomy, Assemble, Annotation Review,
  and Export views, placed after Taxon. Data comes from
  `.specimen_status_join()` on the table's data frame.
- Clicking it opens the viewer modal (`specimen_viewer_server`, ported): sample
  list, per-source tabs with ID input, Fetch, status line, level cards with
  external links, Refresh all, Remove panel; stage 2 adds the Compare tab and
  the "Follow links between databases" checkbox.
- The **Metadata** button in the filter row (stage 3) opens the column picker;
  chosen columns are added to the four views through `meta_view_join()` with
  `meta_view_col_defs()`.
- Refresh across views uses a `reactiveVal` bump in `app_server.R` instead of
  gargoyle.
- Connections: the ported functions take `con`; the app passes `app_con()` and
  disconnects on exit, or uses `app_query()`/`app_write()` where a single
  statement suffices.
- Progress and errors as in MitoPilot: `withProgress` around fetches,
  `showNotification` for failures, `rp_alert` for hard errors.
- Assets: `inst/app/www/specimen/{geome_g,gbif_leaf,ncbi_helix}.png` and the
  `mp-spec-*`, `mp-meta-*`, `mp-token-*` CSS appended to
  `inst/app/www/custom.css`; `inst/extdata/countries.csv`.

## Stage 2: links and conflicts

Ported verbatim from `meta_link.R` and `specimen_reconcile.R`:

- Linking off by default; `link_sources` at setup or the project option
  `meta_options.link_sources`; the checkbox in the viewer is shared by all tabs.
- Fills empty ID columns only, never replaces a user's ID; provenance in
  `meta_links`; a differing ID becomes a note.
- Rules: NCBI to GEOME (bcid), NCBI to GBIF (EZID, then voucher search), GEOME
  to NCBI (child with BioSample), GEOME to GBIF (EZID, then voucher search),
  GBIF to NCBI (associatedSequences), GBIF to GEOME (ARK in occurrence fields).
- `specimen_conflicts()`: nine concepts, statuses `single`, `agree`,
  `conflict`, `not checked`; comparison rules as in MitoPilot (0.01 degrees,
  ymd parts, ISO2 country lookup, voucher normalization, first two words of
  taxon). Mapping-file column detection with `set_metadata_columns()`
  overrides in `meta_csv_map`.
- Compare tab and "Mapfile columns..." panel in the viewer; before an export,
  conflicts on concepts the templates use trigger a confirm dialog "Sample
  metadata disagrees" with "Export anyway".

## Deviations from MitoPilot

1. NCBI key: `project_ncbi_key(project)` (project `.config`, then
   `NCBI_API_KEY`), falling back to `ENTREZ_KEY`. Rate gap 0.11 s with a key,
   0.34 s without, as in MitoPilot.
2. No gargoyle; a `reactiveVal` bump.
3. No user-assembly project variant (RiboPilot has one project type).
4. No legacy `geome_*` table migration.
5. Tables and columns come from `new_db()` and a schema step, not lazy creation.
6. The Metadata column also appears in the Taxonomy view.
7. Guide figures are not copied from MitoPilot; the guide ships with a list of
   screenshots to take in RiboPilot.

## Tests

Port MitoPilot's test files and fixtures (recorded JSON and XML under
`tests/testthat/fixtures/{gbif,geome,ncbi}/`, helpers `helper-gbif.R`,
`helper-geome.R`, `helper-ncbi.R`) and adapt names. Network is mocked with
`local_mocked_bindings(.gbif_get = ..., .package = "RiboPilot")`. Live tests
are guarded by `skip_on_cran()` and `skip_if_offline()`. New tests for the
schema step (parity), the mapping-file columns, `update_sample_metadata()`
refetch semantics, export tokens through `export_header_data()`, and the viewer
via `shiny::testServer(RiboPilot:::app_server, ...)`.

## Docs (stage 3)

- New page `vignettes/Sample-Metadata.Rmd`, "Sample metadata: GEOME, GBIF, and
  NCBI", following MitoPilot's outline with RiboPilot names; no figures yet.
- `_pkgdown.yml`: article under Usage; reference entries for the new functions.
- `Your-Own-Project.Rmd` mapping-file section mentions the three ID columns.
- `GenBank-Submission.Rmd` "Required source modifiers" points at the fetched
  tokens.
- NEWS entry.

## Dependencies

`DESCRIPTION` Imports gain `httr2` and `xml2`. `curl` is needed by
`skip_if_offline()` in tests (Suggests).

## Non-goals

Merging values across sources, automatic modifier insertion into templates,
private GEOME records, caching to disk, and any change to how RiboPilot
resolves taxonomy.

## Files to port (MitoPilot `R/`)

`gbif_api.R`, `gbif_db.R`, `gbif_export.R`, `geome_api.R`, `geome_db.R`,
`geome_export.R`, `ncbi_api.R`, `ncbi_db.R`, `ncbi_export.R`, `meta_db.R`,
`meta_export.R`, `meta_link.R`, `meta_view.R`, `app_meta_view.R`,
`app_specimen.R`, `specimen_reconcile.R`, `app_export_tokens.R`, and the
metadata parts of `update_sample_metadata.R`, `init_checks.R`, `init_db.R`,
`init_project.R`, `add_samples.R`, `app_export.R`, `app_export_utils.R`,
`export.R`, `app_ui.R`. Internal function names stay as in MitoPilot so the two
packages can be diffed later.
