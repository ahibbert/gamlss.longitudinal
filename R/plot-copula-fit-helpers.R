#' @keywords internal
#' @noRd
.plot_copula_selection_spec <- function(copula) {
  if (inherits(copula, "copula_selection")) {
    copula <- as.data.frame(copula)[1L, , drop = FALSE]
  }
  if (is.character(copula) && length(copula) == 1L) {
    family <- .copula_family_code(copula)
    par <- switch(family,
      N = 0,
      F = 0,
      t = 0,
      C = 1e-8,
      G = 1 + 1e-8,
      J = 1 + 1e-8,
      0
    )
    par2 <- if (identical(family, "t")) 4 else 0
    return(list(family = family, par = par, par2 = par2, tau = 0))
  }
  if (is.data.frame(copula) && nrow(copula) >= 1L) {
    family <- .copula_family_code(as.character(copula$family[1L]))
    par <- if ("par" %in% names(copula)) as.numeric(copula$par[1L]) else NA_real_
    par2 <- if ("par2" %in% names(copula)) as.numeric(copula$par2[1L]) else 0
    tau <- if ("tau" %in% names(copula)) as.numeric(copula$tau[1L]) else NA_real_
    if (!is.finite(par) && is.finite(tau)) {
      par <- .copula_tau_to_par(family, tau)
    }
    return(list(family = family, par = par, par2 = par2, tau = tau))
  }
  stop("'copula' must be a copula_selection result, one-row selection data frame, or family code.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.plot_copula_candidate_families <- function(copula) {
  if (!is.character(copula) || length(copula) < 1L || any(is.na(copula))) {
    stop("'copula_dist' must be a copula_selection result, one-row selection data frame, or family code.", call. = FALSE)
  }
  if (length(copula) == 1L && tolower(copula) %in% c("best", "auto")) {
    return(c("N", "C", "F", "G", "J", "t"))
  }
  vapply(copula, .copula_family_code, character(1), USE.NAMES = FALSE)
}

#' @keywords internal
#' @noRd
.plot_copula_resolve_spec <- function(copula, pair_data = NULL, min_pairs = 3L) {
  if (is.character(copula)) {
    if (is.null(pair_data)) {
      return(list(spec = .plot_copula_selection_spec(copula[1L]), selection = NULL))
    }
    pair_data <- pair_data[is.finite(pair_data$u1) & is.finite(pair_data$u2), , drop = FALSE]
    if (nrow(pair_data) < min_pairs) {
      stop("Need at least ", min_pairs, " finite pseudo-observation pairs to fit the copula overlay.", call. = FALSE)
    }
    selection <- select_copula(
      u1 = pair_data$u1,
      u2 = pair_data$u2,
      families = .plot_copula_candidate_families(copula),
      min_pairs = min_pairs
    )
    return(list(spec = .plot_copula_selection_spec(selection), selection = selection))
  }
  list(spec = .plot_copula_selection_spec(copula), selection = if (inherits(copula, "copula_selection")) copula else NULL)
}

#' @keywords internal
#' @noRd
.plot_copula_density_for_spec <- function(pair_data, spec, grid_n, max_pairs_overlay) {
  pair_data$theta_pair <- if ("theta_pair" %in% names(pair_data)) pair_data$theta_pair else rep(spec$par, nrow(pair_data))
  pair_data$zeta_pair <- if ("zeta_pair" %in% names(pair_data)) pair_data$zeta_pair else rep(spec$par2, nrow(pair_data))
  .copula_v2_average_density_grid(
    family_num = spec$family,
    pair_data = pair_data,
    grid_n = grid_n,
    max_pairs_overlay = max_pairs_overlay
  )
}

#' @keywords internal
#' @noRd
.plot_copula_transform_grid <- function(density_grid, transform) {
  if (identical(transform, "normal")) {
    z1 <- stats::qnorm(.copula_v2_clamp01(density_grid$u1))
    z2 <- stats::qnorm(.copula_v2_clamp01(density_grid$u2))
    density_grid$density <- density_grid$density * stats::dnorm(z1) * stats::dnorm(z2)
    density_grid$u1 <- z1
    density_grid$u2 <- z2
  }
  density_grid
}

.plot_copula_pair_time_group <- function(pair_data) {
  if ("time_pair" %in% names(pair_data)) {
    return(factor(pair_data$time_pair, levels = unique(pair_data$time_pair)))
  }
  if ("copula_time" %in% names(pair_data)) {
    return(factor(pair_data$copula_time, levels = unique(pair_data$copula_time)))
  }
  stop("'by_time = TRUE' requires pair data with adjacent time-pair labels.", call. = FALSE)
}
