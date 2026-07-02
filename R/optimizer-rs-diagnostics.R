#' Record one elapsed RS timer checkpoint
#'
#' @noRd
.gl_record_rs_timer_step <- function(
    timer,
    timer_start,
    label,
    elapsed_sec = NULL,
    difftime_fn = difftime,
    sys_time_fn = Sys.time) {
  if (is.null(elapsed_sec)) {
    elapsed_sec <- as.numeric(difftime_fn(sys_time_fn(), timer_start, units = "secs"))
  }

  timer <- c(timer, elapsed_sec)
  names(timer)[length(timer)] <- label
  timer
}

#' Plot RS optimisation progress diagnostics
#'
#' @noRd
.gl_plot_rs_progress <- function(
    log_lik_history,
    par_history,
    par_count,
    true_val,
    par_fn = graphics::par,
    plot_fn = graphics::plot,
    abline_fn = graphics::abline) {
  plot_count <- 3 + par_count
  sides <- round(sqrt(plot_count))

  par_fn(mfrow = c(sides + 1, sides))

  plot_fn(log_lik_history[, 3], type = "l", main = "LogLik - Overall")
  plot_fn(log_lik_history[, 1], type = "l", main = "LogLik - Margin")
  plot_fn(log_lik_history[, 2], type = "l", main = "LogLik - Copula")

  for (i in 1:length(colnames(par_history))) {
    if (!all(is.na(true_val))) {
      plot_fn(
        par_history[, i],
        type = "l",
        main = colnames(par_history)[i],
        xlab = "Iteration",
        ylab = "Parameter estimate",
        ylim = range(c(par_history[, i], true_val[i]))
      )
      abline_fn(h = true_val[i], col = "red")
    } else {
      plot_fn(
        par_history[, i],
        type = "l",
        main = colnames(par_history)[i],
        xlab = "Iteration",
        ylab = "Parameter estimate"
      )
    }
  }

  invisible(TRUE)
}
