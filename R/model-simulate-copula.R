.gl_simulate_copula_matrix <- function(object, diag_data, nsim, fit_data = NULL, time_levels = NULL) {
  if (is.null(fit_data)) {
    fit_data <- .copula_v2_fit_data(object)
  }

  n <- nrow(fit_data)

  out <- matrix(NA_real_, nrow = n, ncol = nsim)

  if (is.null(time_levels)) {
    time_levels <- .gl_simulation_time_levels(fit_data$time)
  }

  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))

  fit_data$.time_idx <- unname(time_lookup[as.character(fit_data$time)])

  fit_data$.row_id <- seq_len(n)

  fit_data_ordered <- fit_data[order(fit_data$subject, fit_data$.time_idx), , drop = FALSE]
  edge_groups <- .gl_simulation_edge_groups(fit_data, fit_data_ordered)
  family_args <- attr(diag_data$params, "family_args", exact = TRUE) %||% list()
  batch_size <- .gl_simulation_batch_size(n, nsim)

  for (batch_start in seq.int(1L, nsim, by = batch_size)) {
    batch_columns <- batch_start:min(nsim, batch_start + batch_size - 1L)
    batch_n <- length(batch_columns)
    target_ordered <- matrix(stats::runif(n * batch_n), nrow = n, ncol = batch_n)
    u <- matrix(NA_real_, nrow = n, ncol = batch_n)
    u[fit_data_ordered$.row_id, ] <- target_ordered

    for (edges in edge_groups) {
      current_rows <- edges$current
      left_rows <- edges$left
      theta <- fit_data$theta_fit[left_rows]
      dependent <- is.finite(theta)
      if (!any(dependent)) next

      current_rows <- current_rows[dependent]
      left_rows <- left_rows[dependent]
      theta <- theta[dependent]
      zeta <- fit_data$zeta_fit[left_rows]
      zeta[!is.finite(zeta)] <- 0

      inverted <- .sim_invert_hfunc1(
        u1 = as.vector(u[left_rows, , drop = FALSE]),
        target = as.vector(u[current_rows, , drop = FALSE]),
        family = object$copula_dist,
        par = rep(theta, times = batch_n),
        par2 = rep(zeta, times = batch_n)
      )
      u[current_rows, ] <- matrix(inverted, nrow = length(current_rows), ncol = batch_n)
    }

    margin_args <- c(
      .gl_simulation_repeat_args(family_args, n, batch_n),
      .gl_simulation_repeat_args(diag_data$params, n, batch_n)
    )
    out[, batch_columns] <- matrix(
      .gl_simulation_quantiles(
        family = diag_data$family,
        p = as.vector(u),
        args = margin_args
      ),
      nrow = n,
      ncol = batch_n
    )
  }

  out
}

.gl_simulation_edge_groups <- function(fit_data, fit_data_ordered) {
  row_ids <- fit_data_ordered$.row_id
  previous_rows <- c(NA_integer_, head(row_ids, -1L))
  same_subject <- c(
    FALSE,
    as.character(head(fit_data_ordered$subject, -1L)) ==
      as.character(tail(fit_data_ordered$subject, -1L))
  )
  adjacent <- same_subject &
    fit_data$.time_idx[row_ids] == fit_data$.time_idx[previous_rows] + 1L
  adjacent[is.na(adjacent)] <- FALSE

  current_rows <- row_ids[adjacent]
  left_rows <- previous_rows[adjacent]
  time_index <- fit_data$.time_idx[current_rows]
  edge_frame <- data.frame(current = current_rows, left = left_rows)
  unname(split(edge_frame, time_index))
}

.gl_simulation_batch_size <- function(n, nsim, max_cells = 250000L) {
  if (n < 1L || nsim < 1L) return(1L)
  max(1L, min(as.integer(nsim), as.integer(max_cells %/% n)))
}

.gl_simulation_repeat_args <- function(args, n, batch_n) {
  lapply(args, function(value) {
    rep(rep(value, length.out = n), times = batch_n)
  })
}
