# Fixed synthetic objects, not optimizer fits or saved-model replays.
curvature_model_fixture <- function(family = "N", tail = FALSE) {
  n <- 32L
  angle <- (seq_len(n) - 0.5) * 2 * pi / n
  basis <- vapply(1:4, function(j) sqrt(2) * cos(j * angle), numeric(n))
  x <- basis[, 1L]
  rho <- 0.4
  e <- cbind(basis[, 2L],
    rho * basis[, 2L] + sqrt(1 - rho^2) * basis[, 3L],
    rho^2 * basis[, 2L] + rho * sqrt(1 - rho^2) * basis[, 3L] +
      sqrt(1 - rho^2) * basis[, 4L])
  y <- 0.3 + 0.7 * x + exp(-0.2) * e
  if (tail) y[1L, 2L] <- 0.3 + 0.7 * x[1L] - 5.5 * exp(-0.2)
  data <- data.frame(response = as.vector(t(y)), x = rep(x, each = 3L),
    time = rep(1:3, n), subject = rep(seq_len(n), each = 3L))
  theta <- if (family %in% c("N", "t")) rho else 2
  # Write transforms independently from the package's link helpers.
  theta_beta <- switch(family, N = atanh(theta), t = atanh(theta),
    C = log(theta), F = theta, G = log(theta - 1), J = log(theta - 1))
  beta <- c(mu.intercept = 0.3, mu.x = 0.7, sigma.intercept = -0.2,
    theta.intercept = theta_beta)
  matrices <- list(mu = cbind(intercept = 1, x = data$x),
    sigma = matrix(1, nrow(data), 1L, dimnames = list(NULL, "intercept")),
    theta = matrix(1, 2L * n, 1L, dimnames = list(NULL, "intercept")))
  if (family == "t") {
    beta <- c(beta, zeta.intercept = log(3))
    matrices$zeta <- matrices$theta
  }
  empty <- setNames(rep(list(list()), length(matrices)), names(matrices))
  structure(list(par = beta, model_matrix = list(x = matrices, s = empty),
    margin_dist = gamlss.dist::NO(), copula_dist = family,
    response = data$response, response_margin = data$time,
    response_subject = data$subject, dataset = data, dataset_original = data,
    par_s = empty, lambda_s = empty, df_s = empty, weights = empty,
    response_var = "response", time_var = "time", subject_var = "subject",
    var_map = list(response = "response", time = "time", subject = "subject"),
    missingness = list(objective = "ordinary", mode = "fail", criteria_status = "available"),
    likelihood_contract = list(objective = "ordinary", criteria_status = "available"),
    include_dlcopdpar = TRUE, optimization = "joint", optim_method = "RS",
    convergence = list(converged = TRUE, method = "synthetic_stationary_fixture",
      stop_reason = "not_an_optimizer_run"),
    synthetic_fixed_object = TRUE), class = "gamlss.longitudinal")
}

# Independent defining joint likelihood. No package likelihood, copula backend,
# Hessian, pair cache, link helper or model-matrix helper is called here.
curvature_reference_clusters <- function(beta, object) {
  y <- matrix(object$response, ncol = 3L, byrow = TRUE)
  x <- object$dataset$x[seq(1L, length(object$response), 3L)]
  mu <- beta[["mu.intercept"]] + beta[["mu.x"]] * x
  if (!is.null(object$synthetic_smooth_covariates)) {
    mu <- mu + as.numeric(object$synthetic_smooth_covariates %*%
      beta[colnames(object$synthetic_smooth_covariates)])
  }
  sigma <- exp(beta[["sigma.intercept"]])
  z <- (y - mu) / sigma
  out <- rowSums(stats::dnorm(z, log = TRUE) - log(sigma))
  family <- object$copula_dist
  p <- beta[["theta.intercept"]]
  theta <- switch(family, N = tanh(p), t = tanh(p), C = exp(p), F = p,
    G = min(1 + exp(p), 17), J = 1 + exp(p))
  nu <- if (family == "t") 2 + exp(beta[["zeta.intercept"]]) else 5
  # Explicit native likelihood conventions for boundary fixtures. Ordinary
  # fixtures remain in the common native/VineCopula interior.
  backend <- object$reference_backend
  if (is.null(backend)) backend <- "native"
  if (backend == "native") {
    if (family %in% c("N", "t")) theta <- min(max(theta, -0.999999), 0.999999)
    if (family == "t") nu <- max(nu, 2.000001)
  }
  for (j in 1:2) {
    zu <- z[, j]; zv <- z[, j + 1L]
    u <- stats::pnorm(zu); v <- stats::pnorm(zv)
    active_u <- u > 1e-12 & u < 1 - 1e-12
    active_v <- v > 1e-12 & v < 1 - 1e-12
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
    # Use defining normal tail identities rather than round-tripping a CDF
    # near one. Otherwise the finite-difference reference itself loses curvature.
    lu <- log(u); lv <- log(v)
    su <- log1p(-u); sv <- log1p(-v)
    lu[active_u] <- stats::pnorm(zu[active_u],log.p=TRUE)
    lv[active_v] <- stats::pnorm(zv[active_v],log.p=TRUE)
    su[active_u] <- stats::pnorm(zu[active_u],lower.tail=FALSE,log.p=TRUE)
    sv[active_v] <- stats::pnorm(zv[active_v],lower.tail=FALSE,log.p=TRUE)
    if (family == "N" || (family == "t" && backend == "native" && nu >= 1e6)) {
      a <- stats::qnorm(u); b <- stats::qnorm(v)
      a[active_u] <- zu[active_u]; b[active_v] <- zv[active_v]
      lc <- -log(1 - theta^2) / 2 +
        (2 * theta * a * b - theta^2 * (a^2 + b^2)) / (2 * (1 - theta^2))
    } else if (family == "t") {
      a <- stats::qt(u, nu); b <- stats::qt(v, nu)
      upper_u <- active_u & zu > 0; upper_v <- active_v & zv > 0
      a[upper_u] <- stats::qt(exp(su[upper_u]),nu,lower.tail=FALSE)
      b[upper_v] <- stats::qt(exp(sv[upper_v]),nu,lower.tail=FALSE)
      lc <- lgamma((nu + 2) / 2) + lgamma(nu / 2) - 2 * lgamma((nu + 1) / 2) -
        log(1 - theta^2) / 2 - (nu + 2) / 2 *
        log1p((a^2 - 2 * theta * a * b + b^2) / (nu * (1 - theta^2))) +
        (nu + 1) / 2 * (log1p(a^2 / nu) + log1p(b^2 / nu))
    } else if (family == "C") {
      s <- log1p(expm1(-theta * lu) + expm1(-theta * lv))
      lc <- log1p(theta) - (theta + 1) * (lu + lv) - (2 + 1 / theta) * s
    } else if (family == "F") {
      den <- expm1(-theta) + expm1(-theta * u) * expm1(-theta * v)
      lc <- log(-theta * expm1(-theta)) - theta * (u + v) - 2 * log(abs(den))
    } else if (family == "G") {
      a <- -lu; b <- -lv
      s <- (a^theta + b^theta)^(1 / theta)
      lc <- -s + (theta - 1) * (log(a) + log(b)) +
        (1 - 2 * theta) * log(s) + log(s + theta - 1) - lu - lv
    } else {
      a <- exp(theta * su); b <- exp(theta * sv)
      s <- a + b - a * b
      lc <- (theta - 1) * (su + sv) +
        (1 / theta - 2) * log(s) + log(s + theta - 1)
    }
    out <- out + lc
  }
  out
}

curvature_reference_hessian <- function(object) {
  numDeriv::hessian(function(beta) sum(curvature_reference_clusters(beta, object)),
    object$par, method = "Richardson", method.args = list(d = 1e-3, r = 4))
}

curvature_reference_covariance <- function(object) {
  hessian <- curvature_reference_hessian(object)
  bread <- solve(-hessian)
  dimnames(bread) <- list(names(object$par), names(object$par))
  scores <- numDeriv::jacobian(function(beta) curvature_reference_clusters(beta, object),
    object$par, method = "Richardson", method.args = list(d = 1e-3, r = 4))
  n <- length(object$response); g <- nrow(scores); p <- length(object$par)
  sandwich <- bread %*% crossprod(scores) %*% bread * (g / (g - 1)) * ((n - 1) / (n - p))
  list(hessian = hessian, model = bread, sandwich = sandwich, scores = scores)
}
boundary_followup_compare <- function(object, key) {
  actual <- calc_analytical_hessian(object, progress=FALSE)
  reference <- curvature_reference_hessian(object)
  error <- sqrt(sum((actual-reference)^2))/sqrt(sum(reference^2))
  if (exists("validation_details", inherits=TRUE))
    validation_details$hessians[[key]] <- list(relative_error=error,
      actual=actual,reference=reference)
  expect_true(error < 2e-4,info=key)
}
