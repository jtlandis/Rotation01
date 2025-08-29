#' vector replicating
#' @param vec a vector
#' @param reps vector of integers indicating how many times vec[i] is replicated
#' @export
vec_rep <- function(vec, reps) {
  n <- length(vec)
  stopifnot(
    "vec and reps must be the same length" = n == length(vec),
    "reps must be an integer" = is.integer(reps <- as.integer(reps))
  )
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

comb2 <- function(vec) {
  n <- length(vec)
  seq <- (n - 1):1
  reps_1 <- rep(vec, c(seq, 0))
  reps_2 <- vector(mode(reps_1), length(reps_1))
  ind <- cumsum(seq)
  for (i in n:2) {
    reps_2[ind] <- vec[i]
    ind <- ind[seq_len(i - 2L)] - 1L
  }
  matrix(c(reps_1, reps_2), nrow = 2, byrow = TRUE)
}
