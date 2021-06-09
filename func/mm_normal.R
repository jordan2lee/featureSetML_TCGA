mm_normal <- function(x) {
  #' Calculate min-max normalization of input values
  return ((x - min(x)) / (max(x) - min(x)))
}
