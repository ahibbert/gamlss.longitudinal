# Internal BCPE/log marginal observed curvature. Row arithmetic is ported from
# the independently checked mu-curvature.R; no design, penalty or copula terms.
.hessian_bcpe_mu_plain <- function(x) {
  typeof(x) %in% c("integer", "double") &&
    is.null(attr(x, "class", exact = TRUE)) && is.null(dim(x)) &&
    length(x) > 0L && all(is.finite(x))
}

.hessian_bcpe_mu_log_abs_expm1 <- function(x) {
  # Caller supplies one finite, nonzero x. No exp(x) overflow for positive x.
  if (x > 0) x + log(-expm1(-x)) else log(-expm1(x))
}

.hessian_bcpe_mu_curvature <- function(y, mu, sigma, nu, tau, eta_radius = NULL) {
  v <- list(y = y, mu = mu, sigma = sigma, nu = nu, tau = tau)
  if (!all(vapply(v, .hessian_bcpe_mu_plain, logical(1L))))
    stop("Finite plain numeric vectors required")
  n <- max(vapply(v, length, integer(1L)))
  if (any(!vapply(v, length, integer(1L)) %in% c(1L, n)))
    stop("Only scalar broadcasting is supported")
  v <- lapply(v, function(x) if (length(x) == 1L) rep(x, n) else x)
  if (any(v$y <= 0 | v$mu <= 0 | v$sigma <= 0 | v$tau <= 0))
    stop("Positive y, mu, sigma and tau required")
  if (!is.null(eta_radius)) {
    if (!.hessian_bcpe_mu_plain(eta_radius) || !length(eta_radius) %in% c(1L, n) ||
        any(eta_radius < 0)) stop("Finite nonnegative eta radius required")
    if (length(eta_radius) == 1L) eta_radius <- rep(eta_radius, n)
  }
  value <- rep(NA_real_, n); available <- rep(FALSE, n)
  reason <- rep("uncomputed", n); log_ratio <- rep(NA_real_, n)
  zero_in_eta_interval <- rep(NA, n)
  for (i in seq_len(n)) {
    yi <- v$y[i]; mi <- v$mu[i]; si <- v$sigma[i]
    ni <- v$nu[i]; ti <- v$tau[i]
    # Near one, avoid cancellation of two logs. Else avoid ratio overflow.
    lr <- if (yi >= mi / 2 && yi / 2 <= mi)
      log1p((yi - mi) / mi) else log(yi) - log(mi)
    log_ratio[i] <- lr
    if (!is.null(eta_radius)) zero_in_eta_interval[i] <- abs(lr) <= eta_radius[i]
    lc <- if (ti == 2) 0 else
      0.5 * (-(2 / ti) * log(2) + lgamma(1 / ti) - lgamma(3 / ti))
    cs <- exp(lc)
    if (!is.finite(lr) || !is.finite(lc) || !is.finite(cs) || cs <= 0) {
      reason[i] <- "unrepresentable_transform_or_scale"; next
    }
    if (lr == 0 && yi != mi) {
      reason[i] <- "log_ratio_rounded_to_zero"; next
    }
    if (yi == mi || (ni == 0 && ti == 2)) {
      if (ti > 2) {
        value[i] <- 0; available[i] <- TRUE; reason[i] <- "smooth_zero"; next
      }
      if (ti < 2) {
        reason[i] <- "no_finite_second_derivative_at_zero"; next
      }
      magnitude <- exp(-2 * log(si))
      if (!is.finite(magnitude) || magnitude == 0) {
        reason[i] <- "curvature_overflow_or_underflow"; next
      }
      value[i] <- -magnitude; available[i] <- TRUE; reason[i] <- "normal_shape_identity"; next
    }
    log_r <- ni * lr
    if (!is.finite(log_r) || (ni != 0 && log_r == 0)) {
      reason[i] <- "unrepresentable_power_transform"; next
    }
    log_abs_z <- if (ni == 0) log(abs(lr)) - log(si) else
      .hessian_bcpe_mu_log_abs_expm1(log_r) - log(abs(ni)) - log(si)
    abs_z <- exp(log_abs_z)
    if (!is.finite(log_abs_z) || !is.finite(abs_z) || abs_z == 0) {
      reason[i] <- "z_overflow_or_underflow"; next
    }
    if (ni == 0 && ti == 1) {
      value[i] <- 0; available[i] <- TRUE; reason[i] <- "linear_log_density_off_zero"; next
    }
    log_tau_r <- log(ti) + log_r
    if (!is.finite(log_tau_r) || log_tau_r == 0) {
      # A true smooth zero is possible, but a rounded cancellation cannot certify it.
      reason[i] <- "unresolved_log_bracket_cancellation"; next
    }
    log_abs_value <- log(ti) - log(2) - ti * lc +
      (ti - 2) * log_abs_z + log_r - 2 * log(si) +
      .hessian_bcpe_mu_log_abs_expm1(log_tau_r)
    magnitude <- exp(log_abs_value)
    if (!is.finite(log_abs_value) || !is.finite(magnitude) || magnitude == 0) {
      reason[i] <- "curvature_overflow_or_underflow"; next
    }
    value[i] <- -sign(log_tau_r) * magnitude
    available[i] <- TRUE; reason[i] <- "available"
  }
  list(curvature = value, available = available, reason = reason,
       log_ratio = log_ratio, zero_in_eta_interval = zero_in_eta_interval)
}

# Return an eta-scale override only for the exact family/link it describes.
.hessian_bcpe_curvature_eta_rows <- function(eta_inv, margin_dist, response) {
  if (!identical(unname(margin_dist$family[1L]), "BCPE") ||
      !identical(margin_dist$mu.link, "log")) return(eta_inv)
  # calc_eta() retains n-by-1 products for fitted smooth contributions. Adapt
  # that internal representation only; the curvature kernel remains vector-only.
  for (pn in c("mu", "sigma", "nu", "tau")) {
    x <- eta_inv[[pn]]
    shape <- dim(x)
    if (is.null(shape)) next
    if (!typeof(x) %in% c("integer", "double") ||
        !is.null(attr(x, "class", exact = TRUE)) || length(shape) != 2L ||
        shape[2L] != 1L || !shape[1L] %in% c(1L, length(response))) {
      stop("BCPE curvature internal parameter '", pn,
        "' must be a plain numeric scalar or response-sized column matrix.", call. = FALSE)
    }
    eta_inv[[pn]] <- as.numeric(x)
  }
  eta_inv
}

.hessian_bcpe_log_mu_eta_curvature <- function(eta_inv, margin_dist, response) {
  if (!identical(unname(margin_dist$family[1L]), "BCPE") ||
      !identical(margin_dist$mu.link, "log")) return(NULL)
  if (!typeof(response) %in% c("integer", "double") || !is.null(dim(response)) ||
      !is.null(attr(response, "class", exact = TRUE)))
    stop("BCPE curvature requires a numeric response vector.", call. = FALSE)
  n <- length(response)
  observed <- is.finite(response)
  value <- numeric(n)
  if (!any(observed)) return(value)
  fail_rows <- function(rows, reasons) {
    stop(structure(list(
      message = paste0("BCPE/log marginal mu curvature unavailable at observed row(s) ",
        paste(utils::head(rows, 10L), collapse = ", "), "."),
      call = NULL, rows = rows, reasons = reasons,
      context = "BCPE/log marginal eta-mu curvature"),
      class = c("gamlss_longitudinal_score_evaluation_error", "error", "condition")))
  }
  args <- list(y = response[observed])
  for (pn in c("mu", "sigma", "nu", "tau")) {
    x <- eta_inv[[pn]]
    if (!typeof(x) %in% c("integer", "double") || !is.null(dim(x)) ||
        !is.null(attr(x, "class", exact = TRUE)) || !length(x) %in% c(1L, n))
      stop("BCPE curvature parameter '", pn, "' must be a plain scalar or row vector.",
        call. = FALSE)
    if (length(x) == 1L) x <- rep(x, n)
    args[[pn]] <- x[observed]
  }
  invalid <- !is.finite(args$mu) | !is.finite(args$sigma) |
    !is.finite(args$nu) | !is.finite(args$tau) |
    args$y <= 0 | args$mu <= 0 | args$sigma <= 0 | args$tau <= 0
  invalid[is.na(invalid)] <- TRUE
  if (any(invalid)) fail_rows(which(observed)[invalid],
    rep("invalid_observed_BCPE_domain", sum(invalid)))
  rows <- do.call(.hessian_bcpe_mu_curvature, args)
  bad <- !rows$available | !is.finite(rows$curvature)
  if (any(bad)) fail_rows(which(observed)[bad], rows$reason[bad])
  value[observed] <- rows$curvature
  value
}

