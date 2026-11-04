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
    u[, time_index] <- .sim_invert_hfunc1(
      u1 = u[, time_index - 1L],
      target = target,
      family = family,
      par = theta_mat[, edge_index],
      par2 = zeta_mat[, edge_index]
    )
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

.sim_invert_hfunc1 <- function(u1, target, family, par, par2 = 0) {
  eps <- sqrt(.Machine$double.eps)
  family <- .copula_family_code(family)
  values <- .copula_recycle(
    .copula_clamp01(u1),
    pmin(pmax(as.numeric(target), eps), 1 - eps),
    as.numeric(par),
    as.numeric(par2)
  )
  u1 <- values[[1L]]
  target <- values[[2L]]
  par <- values[[3L]]
  par2 <- values[[4L]]

  if (any(!is.finite(c(u1, target, par)))) {
    stop("Non-finite copula input encountered during simulation.", call. = FALSE)
  }

  if (.copula_backend() != "native") {
    .copula_require_vinecopula("the delegated VineCopula inverse h-function backend")
    out <- VineCopula::BiCopHinv1(
      u1,
      target,
      family = .copula_family_number(family),
      par = par,
      par2 = par2
    )
    return(pmin(pmax(out, eps), 1 - eps))
  }

  out <- switch(family,
    N = .sim_invert_gaussian_hfunc1(u1, target, par),
    C = .sim_invert_clayton_hfunc1(u1, target, par),
    F = .sim_invert_frank_hfunc1(u1, target, par),
    t = .sim_invert_t_hfunc1(u1, target, par, par2),
    G = .sim_invert_archimedean_hfunc1(u1, target, family, par, par2),
    J = .sim_invert_archimedean_hfunc1(u1, target, family, par, par2)
  )

  if (any(!is.finite(out))) {
    stop("Non-finite copula inverse h-function encountered during simulation.", call. = FALSE)
  }

  pmin(pmax(out, eps), 1 - eps)
}

.sim_invert_gaussian_hfunc1 <- function(u1, target, par) {
  rho <- .copula_gaussian_rho(par)
  stats::pnorm(
    rho * stats::qnorm(u1) + sqrt(1 - rho^2) * stats::qnorm(target)
  )
}

.sim_invert_t_hfunc1 <- function(u1, target, par, par2) {
  if (any(!is.finite(par2))) {
    stop("Non-finite t copula degrees of freedom encountered during simulation.", call. = FALSE)
  }

  rho <- .copula_gaussian_rho(par)
  df <- .copula_t_df(par2)
  x <- stats::qt(u1, df = df)
  scale <- sqrt((df + x^2) * (1 - rho^2) / (df + 1))
  stats::pt(
    rho * x + scale * stats::qt(target, df = df + 1),
    df = df
  )
}

.sim_invert_clayton_hfunc1 <- function(u1, target, par) {
  theta <- .copula_clayton_theta(par)
  out <- target
  dependent <- theta > 1e-10

  if (any(dependent)) {
    theta_dep <- theta[dependent]
    u1_dep <- u1[dependent]
    target_dep <- target[dependent]
    u1_power <- exp(-theta_dep * log(u1_dep))
    inside <- 1 + u1_power * expm1(
      -theta_dep / (theta_dep + 1) * log(target_dep)
    )
    out[dependent] <- exp(-log(inside) / theta_dep)
  }

  out
}

.sim_invert_frank_hfunc1 <- function(u1, target, par) {
  theta <- .copula_frank_par(par)
  out <- target
  dependent <- abs(theta) > 1e-8

  if (any(dependent)) {
    theta_dep <- theta[dependent]
    target_dep <- target[dependent]
    exp_u1 <- exp(-theta_dep * u1[dependent])
    numerator <- target_dep * expm1(-theta_dep)
    denominator <- exp_u1 - target_dep * expm1(-theta_dep * u1[dependent])
    out[dependent] <- -log1p(numerator / denominator) / theta_dep
  }

  out
}

.sim_invert_archimedean_hfunc1 <- function(u1, target, family, par, par2, tol = 1e-10) {
  eps <- sqrt(.Machine$double.eps)
  out <- target
  dependent <- par > 1 + 1e-10

  if (!any(dependent)) return(out)

  u1_dep <- u1[dependent]
  target_dep <- target[dependent]
  par_dep <- par[dependent]
  par2_dep <- par2[dependent]
  lower <- rep(eps, length(target_dep))
  upper <- rep(1 - eps, length(target_dep))

  lower_value <- .copula_hfunc1(
    u1_dep, lower, family = family, par = par_dep, par2 = par2_dep
  )
  upper_value <- .copula_hfunc1(
    u1_dep, upper, family = family, par = par_dep, par2 = par2_dep
  )

  if (any(!is.finite(c(lower_value, upper_value)))) {
    stop("Non-finite copula h-function encountered during simulation.", call. = FALSE)
  }

  at_lower <- lower_value >= target_dep
  at_upper <- upper_value <= target_dep
  active <- !(at_lower | at_upper)
  result <- rep(NA_real_, length(target_dep))
  result[at_lower] <- eps
  result[at_upper] <- 1 - eps

  max_iterations <- ceiling(log2((1 - 2 * eps) / tol))
  for (iteration in seq_len(max_iterations)) {
    if (!any(active)) break
    middle <- (lower[active] + upper[active]) / 2
    value <- .copula_hfunc1(
      u1_dep[active],
      middle,
      family = family,
      par = par_dep[active],
      par2 = par2_dep[active]
    )
    if (any(!is.finite(value))) {
      stop("Non-finite copula h-function encountered during simulation.", call. = FALSE)
    }
    move_lower <- value < target_dep[active]
    active_index <- which(active)
    lower[active_index[move_lower]] <- middle[move_lower]
    upper[active_index[!move_lower]] <- middle[!move_lower]
  }

  result[active] <- (lower[active] + upper[active]) / 2
  out[dependent] <- result
  out
}
