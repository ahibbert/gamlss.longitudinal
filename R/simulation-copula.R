.sim_copula_uniform_matrix <- function(n, n_time, family, copula_params, long_data, subject_var, time_var) {
  u <- matrix(stats::runif(n * n_time), nrow = n, ncol = n_time)
  if (n_time == 1L) {
    return(u)
  }

  edge_data <- .sim_edge_data(long_data, n_time, subject_var, time_var)
  theta <- .sim_resolve_copula_theta(copula_params, edge_data, n, n_time - 1L, family)
  zeta <- .sim_resolve_copula_zeta(copula_params, edge_data, n, n_time - 1L, family)

  theta_mat <- matrix(theta, nrow = n, ncol = n_time - 1L, byrow = TRUE)
  zeta_mat <- matrix(zeta, nrow = n, ncol = n_time - 1L, byrow = TRUE)

  for (time_index in 2:n_time) {
    target <- stats::runif(n)
    edge_index <- time_index - 1L
    for (subject_index in seq_len(n)) {
      u[subject_index, time_index] <- .sim_invert_hfunc1(
        u1 = u[subject_index, time_index - 1L],
        target = target[subject_index],
        family = family,
        par = theta_mat[subject_index, edge_index],
        par2 = zeta_mat[subject_index, edge_index]
      )
    }
  }

  attr(u, "edge_truth") <- list(theta = theta_mat, zeta = zeta_mat)
  u
}

.sim_edge_data <- function(long_data, n_time, subject_var, time_var) {
  edge_data <- long_data[long_data$.sim_time_index < n_time, , drop = FALSE]
  edge_data$.sim_edge_index <- edge_data$.sim_time_index
  edge_data$time_left <- edge_data[[time_var]]
  edge_data$time_right <- long_data[[time_var]][long_data$.sim_time_index > 1L]
  edge_data
}

.sim_resolve_copula_theta <- function(copula_params, edge_data, n, n_edge, family) {
  if (!is.null(copula_params$tau)) {
    tau <- .sim_eval_edge_param(copula_params$tau, edge_data, n, n_edge, "copula_params$tau")
    return(.copula_tau_to_par(family, tau))
  }
  theta_spec <- copula_params$theta
  if (is.null(theta_spec)) {
    theta_spec <- copula_params$par
  }
  if (is.null(theta_spec)) {
    stop("copula_params must include theta, par, or tau.", call. = FALSE)
  }
  .sim_eval_edge_param(theta_spec, edge_data, n, n_edge, "copula_params$theta")
}

.sim_resolve_copula_zeta <- function(copula_params, edge_data, n, n_edge, family) {
  zeta_spec <- copula_params$zeta
  if (is.null(zeta_spec)) {
    zeta_spec <- copula_params$par2
  }
  if (is.null(zeta_spec)) {
    zeta_spec <- if (identical(family, "t")) 4 else 0
  }
  .sim_eval_edge_param(zeta_spec, edge_data, n, n_edge, "copula_params$zeta")
}

.sim_invert_hfunc1 <- function(u1, target, family, par, par2) {
  eps <- sqrt(.Machine$double.eps)
  target <- min(max(target, eps), 1 - eps)
  objective <- function(u2) {
    .copula_hfunc1(u1, u2, family = family, par = par, par2 = par2) - target
  }
  lower_value <- objective(eps)
  upper_value <- objective(1 - eps)
  if (!is.finite(lower_value) || !is.finite(upper_value)) {
    stop("Non-finite copula h-function encountered during simulation.", call. = FALSE)
  }
  if (lower_value >= 0) {
    return(eps)
  }
  if (upper_value <= 0) {
    return(1 - eps)
  }
  stats::uniroot(objective, interval = c(eps, 1 - eps), tol = 1e-10)$root
}
