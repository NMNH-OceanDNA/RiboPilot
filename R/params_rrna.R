# Per-clade curation rulesets for the nuclear rDNA repeat. Mirrors MitoPilot's
# params_<clade> design: a shared generic default plus thin clade overrides.
# A ruleset supplies expected length ranges per element (min, max) and the
# minimum reference identity for a conserved anchor before it is flagged.
#
# Adding a clade = add a switch() branch in rrna_params(), its key to
# RRNA_RULESETS, a lineage token in ruleset_for_lineage(), and (optionally) a
# curate_<clade>_rrna() wrapper in curate_rrna.R.

#' Generic broad-metazoan length tolerances (bp), deliberately wide.
#' @noRd
RRNA_LENGTHS_METAZOA <- list(
  "18S"  = c(1600L, 2100L),
  "5.8S" = c(140L, 180L),
  "28S"  = c(2000L, 6000L), # recovered portion varies a lot from skims
  "ITS1" = c(120L, 1600L),
  "ITS2" = c(120L, 1600L)
)

#' Invertebrate (echinoderm + cnidarian) anchor tolerances.
#'
#' 18S is tightened from measured complete records in the target clades
#' (Asteroidea / Scyphozoa / Octocorallia label DBs): complete 18S tops out
#' around 1810-1850 bp, so the generic upper bound of 2100 is loose. 5.8S is
#' near-universal. 28S is left wide on purpose: complete 28S is long and
#' variable and GenBank records for these clades are dominated by partial
#' D1-D3 barcode fragments, so a tight range would false-flag correct data
#' (distinguishing "correctly partial" from "wrong" needs the partial concept
#' that is deferred to a later pass, see below). ITS bounds are provisional -
#' the spacers are hypervariable and clean per-region length measurements for
#' these clades are scarce.
#' @noRd
RRNA_LENGTHS_METAZOA_INVERT <- utils::modifyList(RRNA_LENGTHS_METAZOA, list(
  "18S"  = c(1700L, 1900L),
  "5.8S" = c(150L, 170L),
  "ITS1" = c(120L, 1200L),
  "ITS2" = c(120L, 1200L)
))

#' Per-element reference-identity floors (%). The anchors evolve at very
#' different rates: 18S and 5.8S are near-invariant within a phylum, while 28S is
#' a mosaic of conserved core and hypervariable expansion segments whose overall
#' identity between families in one class routinely drops below 85%. A single
#' scalar therefore both false-flags a correct 28S and lets a divergent 18S (a
#' contamination signature) pass. Values are a defensible starting point, not a
#' published standard, and are meant to be tuned per clade.
#' @noRd
RRNA_ANCHOR_PID_METAZOA <- c("18S" = 95, "5.8S" = 93, "28S" = 80)
RRNA_ANCHOR_PID_FISH    <- c("18S" = 97, "5.8S" = 95, "28S" = 85)
RRNA_ANCHOR_PID_INVERT  <- c("18S" = 93, "5.8S" = 90, "28S" = 78)

#' Screen tunables, shared by every ruleset until a clade needs its own.
#'
#' The three weights are priors, not measurements, which is why the screen rule
#' ships as Flag rather than Drop. `min_isomer_anchor_pid` is a prior too, but a
#' safer one: it compares two sequences from ONE organism rather than across
#' taxa, and "are these two anchors the same sequence" has a far tighter natural
#' answer than "how similar may two different organisms be". Any change to any
#' value here must bump `SCREEN_VERSION`.
#'
#' Defined here rather than in R/screen.R because a package with no Collate
#' field sources alphabetically, so this file loads first.
#' @noRd
SCREEN_PARAMS_DEFAULT <- list(
  support = 0.30, resolution = 0.30, separation = 0.40,
  min_confidence = 0.6,
  min_isomer_anchor_pid = 99.5
)

#' Curation ruleset for a clade
#'
#' Length ranges are for a COMPLETE, correctly assembled element; the flag/ok
#' model has no partial concept yet, so a genuinely partial anchor is flagged
#' short. `min_anchor_pid` is a NAMED PER-ELEMENT vector (see
#' [RRNA_ANCHOR_PID_METAZOA]); the `curate_opts.min_anchor_pid` column and its
#' options-modal control are still a single scalar and are not read by the
#' engine, which is tracked separately.
#'
#' @param ruleset Clade key; one of [RRNA_RULESETS]. Unknown keys error rather
#'   than silently falling back, so a typo cannot mis-curate a sample.
#' @return list(name, lengths = named list el -> c(min,max), min_anchor_pid,
#'   min_ref_coverage, min_anchor_coverage, check_contaminant, screen).
#' @noRd
rrna_params <- function(ruleset = "metazoa_generic") {
  ruleset <- ruleset %||% "metazoa_generic"
  if (!nzchar(ruleset)) ruleset <- "metazoa_generic"
  base <- list(
    name = "metazoa_generic",
    lengths = RRNA_LENGTHS_METAZOA,
    # Named per-element identity floors; an anchor is flagged low_identity only
    # when the reference alignment also covers enough of it (min_ref_coverage),
    # so a partial reference cannot trigger a false flag.
    min_anchor_pid = RRNA_ANCHOR_PID_METAZOA,
    min_ref_coverage = 0.5,   # fraction of the sample element the ref must cover
    min_anchor_coverage = 90, # % of the rRNA model an anchor must cover (else "truncated")
    check_contaminant = TRUE,
    # Screen tunables. Docs 11 and 12 both site them here so a clade can retune
    # them; until 3c they were top-level constants no ruleset could reach.
    screen = SCREEN_PARAMS_DEFAULT
  )
  invert <- function(name) {
    utils::modifyList(base, list(name = name, lengths = RRNA_LENGTHS_METAZOA_INVERT,
      min_anchor_pid = RRNA_ANCHOR_PID_INVERT))
  }
  switch(ruleset,
    "metazoa_generic" = base,
    "actinopterygii" = utils::modifyList(base, list(
      name = "actinopterygii",
      # Ray-finned fishes: tighter cores; ITS stays wide (varies within genera).
      lengths = utils::modifyList(base$lengths, list(
        "18S"  = c(1750L, 1950L),
        "5.8S" = c(140L, 175L),
        "28S"  = c(2500L, 5000L),
        "ITS1" = c(150L, 1200L),
        "ITS2" = c(150L, 1200L)
      )),
      min_anchor_pid = RRNA_ANCHOR_PID_FISH
    )),
    # Sparse-GenBank invertebrate clades: the nearest reference is often in a
    # different family, and 28S diverges fast, so lower per-element floors avoid
    # drowning the user in false low_identity calls. Contamination is caught by
    # the GenBank screen (blast_screen), not this threshold.
    "asteroidea"   = invert("asteroidea"),
    "scyphozoa"    = invert("scyphozoa"),
    "octocorallia" = invert("octocorallia"),
    stop("Unknown ruleset: ", ruleset, ". Known: ",
      paste(RRNA_RULESETS, collapse = ", "), call. = FALSE)
  )
}

#' The reference-identity floor for an element (handles a per-element named
#' vector or a single scalar for back-compatibility)
#' @noRd
anchor_pid_floor <- function(min_anchor_pid, element) {
  if (is.null(min_anchor_pid)) return(NA_real_)
  if (!is.null(names(min_anchor_pid)) && element %in% names(min_anchor_pid)) {
    return(unname(min_anchor_pid[[element]]))
  }
  if (length(min_anchor_pid) == 1 && is.null(names(min_anchor_pid))) {
    return(unname(min_anchor_pid[1]))
  }
  NA_real_
}

#' Available curation rulesets (for app selectors / validation).
#' @noRd
RRNA_RULESETS <- c("metazoa_generic", "actinopterygii",
                   "asteroidea", "scyphozoa", "octocorallia")

#' Map a resolved lineage string to a clade ruleset
#'
#' Greps the sample's frozen lineage for known clade tokens so a clade sample is
#' curated against its own ruleset instead of silently falling back to
#' `metazoa_generic`. Returns `metazoa_generic` when no clade token matches (the
#' caller records a `no_clade_ruleset` note so the fallback is visible). Octocoral
#' order names from the 2022 reclassification (Scleralcyonacea / Malacalcyonacea)
#' and the older Alcyonacea/Pennatulacea are all matched so the result is robust
#' to which classification NCBI's lineage uses.
#' @param lineage A "; "- or ","-joined lineage string.
#' @return a key in [RRNA_RULESETS].
#' @noRd
ruleset_for_lineage <- function(lineage) {
  if (length(lineage) == 0 || is.na(lineage[1]) || !nzchar(trimws(lineage[1]))) {
    return("metazoa_generic")
  }
  l <- tolower(lineage[1])
  if (grepl("octocorallia|alcyonacea|scleralcyonacea|malacalcyonacea|pennatul|helioporacea", l)) {
    return("octocorallia")
  }
  if (grepl("scyphozoa", l)) return("scyphozoa")
  if (grepl("asteroidea", l)) return("asteroidea")
  if (grepl("actinopteri|actinopterygii", l)) return("actinopterygii")
  "metazoa_generic"
}
