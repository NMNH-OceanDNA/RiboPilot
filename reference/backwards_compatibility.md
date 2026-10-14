# Update an old project to the installed RiboPilot

Upgrades the project database to the current schema one step at a time
and points the project's `.config` at the matching container. Each
change is logged in the database's `schema_history` table. The database
is copied to `.sqlite.bak.<timestamp>` before any step runs, and
`.config` to `.config.bak.<timestamp>` before it is edited. Only the
container line of `.config` changes, and only for the official
`macguigand/ribopilot` image. Projects made before RiboPilot 0.1.0 are
not supported.

## Usage

``` r
backwards_compatibility(path = ".", dry_run = FALSE)
```

## Arguments

- path:

  Project directory (default: working directory).

- dry_run:

  Report what would change without changing anything.

## Value

Invisibly, TRUE when something changed (or would change).
