# Per-clade curation rulesets for the nuclear rDNA repeat. Mirrors MitoPilot's
# params_<clade> design: a shared generic default plus thin clade overrides.
# A ruleset supplies expected length ranges per element (min, max) and the
# minimum reference identity for a conserved anchor before it is flagged.
#
# Adding a clade = add a branch here + (optionally) a curate_<clade>_rrna()
# wrapper in curate_rrna.R.

#' Generic broad-metazoan length tolerances (bp), deliberately wide.
#' @noRd
RRNA_LENGTHS_METAZOA <- list(
  "18S"  = c(1600L, 2100L),
  "5.8S" = c(140L, 180L),
  "28S"  = c(2000L, 6000L), # recovered portion varies a lot from skims
  "ITS1" = c(120L, 1600L),
  "ITS2" = c(120L, 1600L)
)

#' Curation ruleset for a clade
#'
#' @param ruleset Clade key: "metazoa_generic" (default) or "actinopterygii".
#' @return list(name, lengths = named list el -> c(min,max), min_anchor_pid,
#'   check_contaminant).
#' @noRd
rrna_params <- function(ruleset = "metazoa_generic") {
  ruleset <- ruleset %||% "metazoa_generic"
  base <- list(
    name = "metazoa_generic",
    lengths = RRNA_LENGTHS_METAZOA,
    min_anchor_pid = 85,      # ref identity (%) below which a conserved anchor is flagged
    check_contaminant = TRUE
  )
  if (identical(ruleset, "actinopterygii")) {
    base$name <- "actinopterygii"
    # Ray-finned fishes: tighter cores; ITS stays wide (varies within genera).
    base$lengths <- utils::modifyList(base$lengths, list(
      "18S"  = c(1750L, 1950L),
      "5.8S" = c(140L, 175L),
      "28S"  = c(2500L, 5000L),
      "ITS1" = c(150L, 1200L),
      "ITS2" = c(150L, 1200L)
    ))
    base$min_anchor_pid <- 88
  }
  base
}

#' Available curation rulesets (for app selectors / validation).
#' @noRd
RRNA_RULESETS <- c("metazoa_generic", "actinopterygii")
