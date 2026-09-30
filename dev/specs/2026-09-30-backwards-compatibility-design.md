# Backwards compatibility for project databases

Date: 2026-09-30
Status: approved

## Goal

From the first beta release (0.1.0) on, old projects can be brought up to the
installed RiboPilot with one call, `backwards_compatibility()`. Every project
database records which schema version it is on and a history of how it got
there. Databases made before 0.1.0 are not supported.

## Versions and history

- `DB_SCHEMA_VERSION` (integer) counts table-layout changes. It starts at 1 and
  only goes up when the schema changes; package versions move independently.
- One table, `schema_history`, one row per event:

  | column            | type    | meaning                                   |
  |-------------------|---------|-------------------------------------------|
  | event             | TEXT    | `created`, `upgraded`, or `config`        |
  | schema_version    | INTEGER | schema version after the event            |
  | ribopilot_version | TEXT    | installed package version at the time     |
  | time_stamp        | INTEGER | unix seconds                              |
  | note              | TEXT    | one-line description                      |

- The current schema version is `schema_version` on the newest row. No table
  means an unversioned (pre-0.1.0) database.
- `new_db()` creates the table and writes a `created` row.

## Migration steps

- `R/db_migrate.R` holds `DB_MIGRATIONS`, a list where element `k - 1` takes a
  version `k - 1` database to version `k`. Each element is
  `list(note = "...", fn = function(con) ...)`.
- `DB_SCHEMA_VERSION` is `1 + length(DB_MIGRATIONS)`. At 0.1.0 the list is empty.
- Rule for contributors: any table or column change edits `new_db()` and adds
  a step. Released steps are never edited.

## `backwards_compatibility(path = ".", dry_run = FALSE)`

1. Read the project version. Stop with a clear message when the database is
   unversioned (recreate the project) or newer than the installed package
   (install a newer RiboPilot).
2. `dry_run = TRUE` prints pending steps and the container change and writes
   nothing.
3. When steps are pending, copy `.sqlite` to `.sqlite.bak.<timestamp>`.
4. Each step runs in its own transaction, followed by an `upgraded` history
   row in the same transaction. A failing step rolls back; completed steps and
   their history stay, and the error names the backup.
5. If `.config` names the official `macguigand/ribopilot:<tag>` image and the
   tag is not the installed version, back up `.config` to
   `.config.bak.<timestamp>`, rewrite only the container line, and write a
   `config` history row. Custom images are left alone.
6. Nothing to do: report "up to date" and write nothing.

## App start

`run_app()` stops before launching when the database is unversioned, newer
than the package, or behind it, naming the fix. The container mismatch alert
stays and now points to `backwards_compatibility()`. Pipeline steps inside the
container are not checked; runs only start from the app.

## Testing

- New project: one `created` row at `DB_SCHEMA_VERSION`.
- Injected fake step: version rises, `upgraded` row, backup file present.
- Failing fake step: database unchanged, no history row for it.
- Unversioned and too-new databases stop with the right message.
- Dry run writes nothing.
- Container line rewritten; custom image untouched.
- Parity: a v1 database, rebuilt from a saved v1 schema fixture, run through
  every step matches a fresh `new_db()`
  table and column set (trivial at v1).

## Out of scope

Pre-0.1.0 databases; the release itself (version bump, NEWS, PR).
