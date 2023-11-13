

#' vector replicating
#' @param vec a vector
#' @param reps vector of integers indicating how many times vec[i] is replicated
#' @export
vec_rep <- function(vec, reps) {
  n <- length(vec)
  stopifnot("vec and reps must be the same length"= n==length(vec),
            "reps must be an integer"= is.integer(reps <- as.integer(reps)))
  m <- sum(reps)
  k <- 0L
  out <- vector(mode(vec), m)
  for (i in seq_along(vec)) {
    k <- k + 1L
    k_ <- k + reps[i] - 1L
    out[k:k_] <- vec[i]
    k <- k_
  }
  out
}