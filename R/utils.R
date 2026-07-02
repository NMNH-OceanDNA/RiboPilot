#' If x is `NULL` (or length 0), return y, otherwise return x
#' @noRd
"%||%" <- function(x, y) {
  if (length(x) == 0L) y else x
}

#' If x is `NA`, return y, otherwise return x
#' @noRd
"%|NA|%" <- function(x, y) {
  if (is.na(x)) y else x
}

#' Inverted %in%
#' @noRd
`%nin%` <- Negate(`%in%`)
