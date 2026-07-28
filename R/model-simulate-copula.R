.gl_simulate_copula_matrix <- function(object, diag_data, nsim, fit_data = NULL, time_levels = NULL) {
  if (is.null(fit_data)) {
    fit_data <- .copula_v2_fit_data(object)
  }

  n <- nrow(fit_data)

  qfun <- get(paste0("q", diag_data$family), envir = asNamespace("gamlss.dist"), inherits = FALSE)

  out <- matrix(NA_real_, nrow = n, ncol = nsim)

  if (is.null(time_levels)) {
    time_levels <- .gl_simulation_time_levels(fit_data$time)
  }

  time_lookup <- stats::setNames(seq_along(time_levels), as.character(time_levels))

  fit_data$.time_idx <- unname(time_lookup[as.character(fit_data$time)])

  fit_data$.row_id <- seq_len(n)

  fit_data_ordered <- fit_data[order(fit_data$subject, fit_data$.time_idx), , drop = FALSE]

  split_rows <- split(fit_data_ordered, fit_data_ordered$subject)

  for (j in seq_len(nsim)) {
    u <- rep(NA_real_, n)

    for (subject_rows in split_rows) {
      row_ids <- subject_rows$.row_id

      if (length(row_ids) == 0L) next

      u[row_ids[[1L]]] <- stats::runif(1L)

      if (length(row_ids) > 1L) {
        for (k in 2:length(row_ids)) {
          left_row <- row_ids[[k - 1L]]

          current_row <- row_ids[[k]]

          if (fit_data$.time_idx[[current_row]] != fit_data$.time_idx[[left_row]] + 1L) {
            u[current_row] <- stats::runif(1L)

            next
          }

          theta <- fit_data$theta_fit[[left_row]]

          zeta <- fit_data$zeta_fit[[left_row]]

          target <- stats::runif(1L)

          if (!is.finite(theta)) {
            u[current_row] <- target
          } else {
            u[current_row] <- .sim_invert_hfunc1(
              u1 = u[left_row],
              target = target,
              family = object$copula_dist,
              par = theta,
              par2 = if (is.finite(zeta)) zeta else 0
            )
          }
        }
      }
    }

    args <- c(list(p = u), diag_data$params)

    args <- args[names(args) %in% formalArgs(qfun)]

    out[, j] <- do.call(qfun, args)
  }

  out
}
