
box::use(
  cli[console_width]
)

base_color <- cli::col_red
ciclic_col <- list(cli::col_cyan, cli::col_grey, cli::col_yellow)


#' for the purpose of this module
#' we do not care about the units when
#' calculating time differences whenever 
#' we calculate a difference.
time <- function(){
  .Internal(Sys.time())
}

time_env <- new.env(parent = emptyenv())
time_env$.start <- NULL
time_env$.stamps <- NULL
time_env$.labels <- NULL
time_env$.end <- NULL
time_env$.i <- NULL
time_env$.last_stamp <- NULL
time_env$.diffs <- NULL

mod_env <- new.env(parent = emptyenv())
mod_env$width <- console_width()

value_width <- 4L
unit_width <- 4L

u <- c("ns", "\u00B5s", "ms", "s", "min", "hr", "day")
t <- c(1e-9, 1e-6, 1e-3, 1, 60, 3600, 86400)

unit <- function(vec) {
  m <- min(vec)
  i <- sum(m > t)
  out <- u[i]
  attr(out, "value") <- vec/t[i]
  out
}



con <- cli::cli_output_connection()

#' @export
time_init <- function(stamp_labels = character()) {
  time_env$.i <- 0L
  time_env$.labels <- stamp_labels
  n <- length(stamp_labels)
  .stamps <- numeric(n)
  time_env$.stamps <- .stamps
  mod_env$width <- console_width()
  time_env$.diffs <- numeric(n + 2L)
  time_env$.last_stamp <- time()
  invisible(time_env)
}

#' @export
time_start <- function() {
  stamp <- time()
  i <- 1L
  time_env$.diffs[i] <- stamp - time_env$.last_stamp
  time_env$.last_stamp <- stamp
  time_env$.i <- i
  invisible(time_env)
}

#' @export
time_stamp <- function(label = NA_character_) {
  stamp <- time()
  i <- time_env$.i + 1L
  # time_env$.labels[i] <- label
  time_env$.diffs[i] <- stamp - time_env$.last_stamp
  time_env$.last_stamp <- stamp
  time_env$.i <- i 
  invisible(time_env)
}

#' @export
time_end <- function() {
  stamp <- time()
  # browser()
  i <- time_env$.i + 1L
  time_env$.diffs[i] <- stamp - time_env$.last_stamp
  time_env$.i <- 0L
  time_env$.last_stamp <- stamp
  diffs <- time_env$.diffs
  overall_time <- sum(diffs[-1L], na.rm = T)
  prop <- diffs[-1L]/overall_time
  bar_space <- mod_env$width - 20L - value_width - unit_width
  bar_est <- ceiling(prop * bar_space)
  remain <- sum(bar_est) %% bar_space
  i <- which.max(bar_est)
  bar_est[i] <- bar_est[i] - remain
  bars <- character(length = length(bar_est))
  # bars[1L] <- base_color(strrep("|", bar_est[1L]))
  for (i in 1:length(bars)) {
    col <- ciclic_col[[(i %% 3L) + 1L]]
    bars[i] <- col(strrep("|", bar_est[i]))
  }
  overall_time <- unit(overall_time)
  overhead <- unit(diffs[1])
  out <- sprintf("\r[%s] %4.1f %-3s (oh: \033[31m%4.1f %-3s\033[39m)", 
          paste(bars, collapse = ""),
          attr(overall_time, "value"),
          overall_time,
          attr(overhead, "value"),
          overhead
          )
  
  cat(out, file = con)
  invisible(out)
  
}
#time_init ---> time_start ---> [time_stamp]... ---> time_end
#                  ^____________________________________/

#' @export
time_finalize <- function() {
  
}
