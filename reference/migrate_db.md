# Bring a project database up to the current schema (idempotent, additive)

Adds any tables/columns introduced after a project was created, without
touching existing data. Safe to call repeatedly. New columns are only
ever added, never dropped or retyped.

## Usage

``` r
migrate_db(con)
```

## Arguments

- con:

  Open DBI connection to a project database.

## Value

(invisibly) the schema version after migration.
