.meta_chr <- function(x) {
  if (!is.numeric(x)) return(as.character(x))
  vapply(x, function(v) if (is.na(v)) NA_character_ else format(v, scientific = FALSE, digits = 15),
         character(1), USE.NAMES = FALSE)
}

.meta_flatten <- function(x, level, depth, ref) {
  keep <- vapply(x, function(v) {
    length(v) == 1L && !is.list(v) && !is.na(v) && nzchar(as.character(v))
  }, logical(1))
  x <- x[keep]
  if (!length(x)) return(NULL)
  data.frame(level = level, depth = as.integer(depth), ref = ref,
             field = names(x), value = vapply(x, .meta_chr, ""),
             row.names = NULL)
}
