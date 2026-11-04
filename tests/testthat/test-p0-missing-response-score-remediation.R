missing_score_margin_spec <- function(family, n = 8L) {
  specs <- list(
    NO = list(response = seq(-1, 1, length.out = n), mu = 0.1, sigma = 0.9),
    GA = list(response = seq(0.5, 3, length.out = n), mu = 2, sigma = 0.5),
    GG = list(response = seq(0.5, 3, length.out = n), mu = 2, sigma = 0.4, nu = 1),
    BCPE = list(response = seq(0.5, 3, length.out = n), mu = 2, sigma = 0.25, nu = 0.5, tau = 2),
    LOGNO = list(response = exp(seq(-0.5, 0.8, length.out = n)), mu = 0.2, sigma = 0.5),
    PO = list(response = rep(0:3, length.out = n), mu = 2),
    NBI = list(response = rep(0:3, length.out = n), mu = 2, sigma = 0.5),
    DEL = list(response = rep(0:3, length.out = n), mu = 3, sigma = 0.7, nu = 0.4),
    ZINBI = list(response = rep(0:3, length.out = n), mu = 2, sigma = 0.5, nu = 0.2),
    ZIBNB = list(response = rep(0:3, length.out = n), mu = 2, sigma = 0.4, nu = 1, tau = 0.15),
    ZIPIG = list(response = rep(0:3, length.out = n), mu = 2, sigma = 0.7, nu = 0.2),
    ZISICHEL = list(response = rep(0:3, length.out = n), mu = 2, sigma = 0.5, nu = -0.5, tau = 0.15)
  )
  specs[[family]]
}

missing_score_family <- function(family) {
  get(family, envir = asNamespace("gamlss.dist"))()
}

missing_score_direct_cdf <- function(family, spec, response) {
  pfun <- get(paste0("p", family), envir = asNamespace("gamlss.dist"))
  formals <- formalArgs(pfun)
  vapply(seq_along(response), function(i) {
    if (!is.finite(response[i])) return(NA_real_)
    args <- c(list(q = response[i]), spec[intersect(names(spec), formals)])
    args$response <- NULL
    as.numeric(do.call(pfun, args))[1L]
  }, numeric(1))
}

test_that("missing-safe CDF evaluation covers paper margins and missingness patterns", {
  skip_if_not_installed("gamlss.dist")
  families <- c(
    "NO", "GA", "GG", "BCPE", "LOGNO",
    "PO", "NBI", "DEL", "ZINBI", "ZIBNB", "ZIPIG", "ZISICHEL"
  )
  patterns <- list(
    complete = integer(),
    leading = c(1L, 2L),
    terminal = c(7L, 8L),
    intermittent = c(2L, 6L),
    segmented = c(4L, 5L)
  )

  for (family in families) {
    margin <- missing_score_family(family)
    spec <- missing_score_margin_spec(family)
    response0 <- spec$response
    eta_inv <- spec[setdiff(names(spec), "response")]
    eta_inv <- lapply(eta_inv, rep, length(response0))

    for (pattern in names(patterns)) {
      response <- response0
      response[patterns[[pattern]]] <- NA_real_
      actual <- gamlss.longitudinal:::calc_F_x(eta_inv, NULL, margin, response)
      expected <- missing_score_direct_cdf(family, spec, response)
      expect_equal(actual[is.finite(response)], expected[is.finite(response)], tolerance = 1e-10,
                   info = paste(family, pattern))
      expect_true(all(is.na(actual[!is.finite(response)])), info = paste(family, pattern))
    }

    if (identical(margin$type, "Discrete")) {
      lower_response <- response0 - 1
      lower <- gamlss.longitudinal:::calc_F_x(eta_inv, NULL, margin, lower_response)
      expect_identical(lower[lower_response < 0], rep(0, sum(lower_response < 0)), info = family)
    }
  }
})

test_that("observed invalid CDF values are classed errors and missing rows are not evaluated", {
  skip_if_not_installed("gamlss.dist")
  margin <- gamlss.dist::NO()
  bad_cdf <- function(q, mu) ifelse(q == 2, Inf, stats::pnorm(q - mu))

  expect_error(
    gamlss.longitudinal:::.gl_evaluate_margin_cdf(
      input = list(q = c(1, 2), mu = c(0, 0)), margin_dist = margin,
      response = c(1, 2), FUN = bad_cdf
    ),
    class = "gamlss_longitudinal_margin_cdf_error"
  )
  expect_equal(
    gamlss.longitudinal:::.gl_evaluate_margin_cdf(
      input = list(q = c(1, NA), mu = c(0, 0)), margin_dist = margin,
      response = c(1, NA), FUN = bad_cdf
    ),
    c(stats::pnorm(1), NA_real_)
  )
})

test_that("missing-safe CDF evaluation retries observed rows for scalar-only families", {
  response <- c(-0.5, NA_real_, 0.25)
  scalar_only <- function(q, mu) {
    if (length(q) != 1L || length(mu) != 1L) {
      stop("scalar inputs required")
    }
    stats::pnorm(q, mean = mu)
  }
  value <- gamlss.longitudinal:::.gl_evaluate_margin_cdf(
    input = list(q = response, mu = c(0.1, 99, -0.2)),
    margin_dist = gamlss.dist::NO(),
    response = response,
    FUN = scalar_only,
    FUN_args = formalArgs(scalar_only)
  )

  expect_equal(
    value[c(1L, 3L)],
    c(stats::pnorm(-0.5, 0.1), stats::pnorm(0.25, -0.2)),
    tolerance = 1e-12
  )
  expect_true(is.na(value[[2L]]))
})

missing_score_discrete_oracle <- function(family, copula, missing = integer()) {
  spec <- missing_score_margin_spec(family, n = 12L)
  response <- spec$response
  response[missing] <- NA_real_
  n <- length(response)
  time <- rep(seq_len(3L), times = 4L)
  subject <- rep(seq_len(4L), each = 3L)
  margin_par <- setdiff(names(spec), "response")
  mm <- setNames(lapply(margin_par, function(x) {
    matrix(1, nrow = n, dimnames = list(NULL, "intercept"))
  }), margin_par)
  mm$theta <- matrix(1, nrow = n, dimnames = list(NULL, "intercept"))
  eta_inv <- lapply(spec[margin_par], rep, n)
  eta_inv$theta <- rep(if (identical(copula, "C")) 1.2 else 0.3, n)
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)
  calc_lik <- gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv, mm = mm, margin_dist = missing_score_family(family),
    copula_dist = copula, response = response, response_margin = time,
    response_subject = subject, pair_cache = pair_cache
  )
  analytical <- gamlss.longitudinal:::.calc_discrete_rectangle_scores(
    eta_inv, mm, missing_score_family(family), copula, response, time, subject,
    pair_cache = pair_cache, calc_lik = calc_lik, method = "analytical", h = 1e-4
  )
  finite <- gamlss.longitudinal:::.calc_discrete_rectangle_scores(
    eta_inv, mm, missing_score_family(family), copula, response, time, subject,
    pair_cache = pair_cache, calc_lik = calc_lik, method = "finite", h = 1e-4
  )
  list(analytical = analytical, finite = finite, response = response)
}

test_that("discrete missing-data analytical scores agree with numerical oracles", {
  skip_if_not_installed("gamlss.dist")
  fixtures <- list(
    zisichel_intermittent = list("ZISICHEL", "G", c(2L, 8L)),
    zisichel_terminal = list("ZISICHEL", "G", c(3L, 6L, 9L, 12L)),
    nbi_complete = list("NBI", "C", integer()),
    del_missing = list("DEL", "C", c(5L, 11L))
  )
  for (fixture in names(fixtures)) {
    args <- fixtures[[fixture]]
    oracle <- missing_score_discrete_oracle(args[[1]], args[[2]], args[[3]])
    for (parameter in names(oracle$finite)) {
      expect_equal(
        oracle$analytical[[parameter]], oracle$finite[[parameter]],
        tolerance = 1e-2, info = paste(fixture, parameter)
      )
    }
    for (parameter in intersect(names(oracle$analytical), c("mu", "sigma", "nu", "tau"))) {
      expect_identical(
        oracle$analytical[[parameter]][!is.finite(oracle$response)],
        rep(0, sum(!is.finite(oracle$response))),
        info = paste(fixture, parameter, "missing rows")
      )
    }
  }
})

missing_score_continuous_oracle <- function(
    family, copula, missing, response = NULL, direction = NULL) {
  spec <- missing_score_margin_spec(family, n = 12L)
  if (is.null(response)) response <- rep(spec$response, length.out = 12L)
  response[missing] <- NA_real_
  n <- length(response)
  time <- rep(seq_len(3L), times = 4L)
  subject <- rep(seq_len(4L), each = 3L)
  margin <- missing_score_family(family)
  copula_link <- gamlss.longitudinal:::get_copula_dist(copula)$copula_link
  margin_par <- setdiff(names(spec), "response")
  model_par <- c(margin_par, "theta", if (identical(copula, "t")) "zeta")
  mm_x <- setNames(lapply(model_par, function(parameter) {
    out <- matrix(1, nrow = n, ncol = 1L)
    colnames(out) <- "(Intercept)"
    out
  }), model_par)
  mm <- list(x = mm_x, s = setNames(lapply(model_par, function(x) list()), model_par))
  par_s <- setNames(lapply(model_par, function(x) NA_real_), model_par)
  beta <- numeric(length(model_par))
  names(beta) <- paste0(model_par, ".(Intercept)")
  for (parameter in margin_par) {
    beta[paste0(parameter, ".(Intercept)")] <-
      margin[[paste0(parameter, ".linkfun")]](spec[[parameter]])
  }
  theta_natural <- switch(copula, N = 0.3, C = 1.2, F = 2, G = 1.5, J = 1.5, t = 0.3)
  beta["theta.(Intercept)"] <- copula_link$theta.linkfun(theta_natural)
  if (identical(copula, "t")) beta["zeta.(Intercept)"] <- copula_link$zeta.linkfun(5)
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)
  evaluate <- function(coef, gradient = FALSE) {
    eta <- gamlss.longitudinal:::calc_eta(coef, mm, margin, copula_link, par_s)
    lik <- gamlss.longitudinal:::calc_likelihood_minimal(
      eta$eta_inv, mm = mm$x, margin_dist = margin, copula_dist = copula,
      calc_d2 = FALSE, response = response, response_margin = time,
      response_subject = subject, pair_cache = pair_cache
    )
    if (!gradient) return(as.numeric(lik$log_lik["joint"]))
    gamlss.longitudinal:::.cg_analytical_gradient(
      coef, mm, eta, lik, margin, copula, TRUE, response, time, subject
    )
  }
  analytical <- evaluate(beta, gradient = TRUE)
  finite <- vapply(seq_along(beta), function(i) {
    h <- 1e-6 * max(1, abs(beta[i]))
    plus <- minus <- beta
    plus[i] <- plus[i] + h
    minus[i] <- minus[i] - h
    (evaluate(plus) - evaluate(minus)) / (2 * h)
  }, numeric(1))
  normalized_error <- max(abs(analytical - finite)) / max(1, max(abs(finite)))
  if (is.null(direction)) direction <- seq_along(beta) / sum(seq_along(beta))
  h <- 1e-5
  numerical_direction <- (evaluate(beta + h * direction) - evaluate(beta - h * direction)) / (2 * h)
  analytical_direction <- sum(analytical * direction)
  direction_error <- abs(analytical_direction - numerical_direction) /
    max(1, abs(numerical_direction))
  list(
    analytical = analytical, finite = finite,
    normalized_error = normalized_error, direction_error = direction_error
  )
}

test_that("continuous missing-data coefficient scores agree across paper routes", {
  skip_if_not_installed("gamlss.dist")
  fixtures <- c(
    setNames(lapply(c("N", "C", "F", "G", "J", "t"), function(cp) {
      list("NO", cp, c(2L, 8L))
    }), paste0("NO_", c("N", "C", "F", "G", "J", "t"))),
    list(
      BCPE_t_dropout = list("BCPE", "t", c(3L, 6L, 9L, 12L)),
      BCPE_t_intermittent = list("BCPE", "t", c(2L, 8L)),
      GA_G = list("GA", "G", c(2L, 8L)),
      GG_C = list("GG", "C", c(2L, 8L)),
      LOGNO_N = list("LOGNO", "N", c(2L, 8L))
    )
  )
  for (fixture in names(fixtures)) {
    args <- fixtures[[fixture]]
    oracle <- missing_score_continuous_oracle(args[[1]], args[[2]], args[[3]])
    expect_true(oracle$normalized_error <= 1e-3, info = fixture)
    expect_true(oracle$direction_error <= 1e-3, info = fixture)
  }
})

missing_score_hessian_fixture <- function(family, copula, missing) {
  spec <- missing_score_margin_spec(family, n = 12L)
  response <- spec$response
  response[missing] <- NA_real_
  n <- length(response)
  time <- rep(seq_len(3L), times = 4L)
  subject <- rep(seq_len(4L), each = 3L)
  margin <- missing_score_family(family)
  copula_link <- gamlss.longitudinal:::get_copula_dist(copula)$copula_link
  margin_par <- setdiff(names(spec), "response")
  model_par <- c(margin_par, "theta", if (identical(copula, "t")) "zeta")
  mm_x <- setNames(lapply(model_par, function(parameter) {
    out <- matrix(1, nrow = n, ncol = 1L)
    colnames(out) <- "(Intercept)"
    out
  }), model_par)
  mm <- list(x = mm_x, s = setNames(lapply(model_par, function(x) list()), model_par))
  par_s <- setNames(lapply(model_par, function(x) NA_real_), model_par)
  beta <- numeric(length(model_par))
  names(beta) <- paste0(model_par, ".(Intercept)")
  for (parameter in margin_par) {
    beta[paste0(parameter, ".(Intercept)")] <-
      margin[[paste0(parameter, ".linkfun")]](spec[[parameter]])
  }
  theta_natural <- if (identical(copula, "G")) 1.5 else 0.3
  beta["theta.(Intercept)"] <- copula_link$theta.linkfun(theta_natural)
  if (identical(copula, "t")) beta["zeta.(Intercept)"] <- copula_link$zeta.linkfun(5)
  object <- list(
    response = response, response_margin = time, response_subject = subject,
    margin_dist = margin, copula_dist = copula, model_matrix = mm,
    par = beta, par_s = par_s
  )
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(response, time, subject)
  evaluate <- function(coef) {
    eta <- gamlss.longitudinal:::calc_eta(coef, mm, margin, copula_link, par_s)
    as.numeric(gamlss.longitudinal:::calc_likelihood_minimal(
      eta$eta_inv, mm = mm$x, margin_dist = margin, copula_dist = copula,
      calc_d2 = FALSE, response = response, response_margin = time,
      response_subject = subject, pair_cache = pair_cache,
      calc_margin_deriv = FALSE
    )$log_lik["joint"])
  }
  list(object = object, evaluate = evaluate, beta = beta)
}

test_that("analytical Hessians match numerical missing-data controls", {
  skip_if_not_installed("gamlss.dist")
  skip_if_not_installed("numDeriv")
  fixtures <- list(
    continuous = missing_score_hessian_fixture("NO", "N", c(2L, 8L)),
    discrete = missing_score_hessian_fixture("ZISICHEL", "G", c(2L, 8L))
  )
  for (fixture in names(fixtures)) {
    x <- fixtures[[fixture]]
    analytical <- gamlss.longitudinal:::calc_analytical_hessian(
      x$object, progress = FALSE, h = 1e-4
    )
    numerical <- numDeriv::hessian(x$evaluate, x$beta, method.args = list(eps = 1e-4))
    dimnames(numerical) <- dimnames(analytical)
    disagreement <- sqrt(sum((analytical - numerical)^2)) /
      max(sqrt(sum(numerical^2)), .Machine$double.eps)
    expect_true(disagreement <= 0.05, info = paste(fixture, disagreement))
  }
})

missing_score_tiny_fit <- function(family, copula, response, missing = integer()) {
  n <- length(response)
  dat <- data.frame(
    id = rep(seq_len(n / 3L), each = 3L),
    time = rep(seq_len(3L), times = n / 3L),
    y = response
  )
  dat$y[missing] <- NA_real_
  suppressWarnings(gamlss_longitudinal(
    data = dat,
    margin_dist = missing_score_family(family),
    copula_dist = copula,
    time_var = "time",
    subject_var = "id",
    mu.formula = "y ~ 1",
    sigma.formula = "~ 1",
    nu.formula = "~ 1",
    tau.formula = "~ 1",
    theta.formula = "~ 1",
    zeta.formula = "~ 1",
    method = "RS",
    optimizer_control = gamlss_longitudinal_control(
      max_outer_iter = 2L,
      rs = list(max_inner_iter = 2L, warm_start_joint = FALSE)
    ),
    compute_vcov = FALSE,
    include_dlcopdpar = TRUE,
    verbose = 0
  ))
}

test_that("tiny end-to-end missing and complete fits retain valid metadata", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("gamlss.dist")
  fits <- list(
    zisichel_g = missing_score_tiny_fit(
      "ZISICHEL", "G", rep(c(0, 1, 2, 3, 1, 0), 3L), c(2L, 8L)
    ),
    bcpe_t_complete = missing_score_tiny_fit(
      "BCPE", "t", exp(seq(-0.4, 0.8, length.out = 18L))
    ),
    bcpe_t_terminal_dropout = missing_score_tiny_fit(
      "BCPE", "t", exp(seq(-0.4, 0.8, length.out = 18L)),
      c(3L, 6L, 9L)
    ),
    bcpe_t_intermittent = missing_score_tiny_fit(
      "BCPE", "t", exp(seq(-0.4, 0.8, length.out = 18L)), c(2L, 8L)
    ),
    bcpe_t_delayed_entry = missing_score_tiny_fit(
      "BCPE", "t", exp(seq(-0.4, 0.8, length.out = 18L)),
      c(1L, 4L, 7L)
    ),
    nbi_c = missing_score_tiny_fit(
      "NBI", "C", rep(c(0, 1, 2, 3, 1, 0), 3L)
    )
  )
  for (fixture in names(fits)) {
    fit <- fits[[fixture]]
    expect_true(inherits(fit, "gamlss.longitudinal"), info = fixture)
    expect_true(is.finite(as.numeric(fit$model_selection["LogLik", "joint"])), info = fixture)
    expect_true(is.logical(fit$convergence$converged) && length(fit$convergence$converged) == 1L,
                info = fixture)
    expect_identical(fit$vcov_meta$stationarity_status, "not_checked", info = fixture)
    eta <- gamlss.longitudinal:::calc_eta(
      fit$par, fit$model_matrix, fit$margin_dist,
      gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link,
      fit$par_s
    )
    likelihood <- gamlss.longitudinal:::calc_likelihood_minimal(
      eta$eta_inv, mm = fit$model_matrix$x, margin_dist = fit$margin_dist,
      copula_dist = fit$copula_dist, response = fit$response,
      response_margin = fit$response_margin, response_subject = fit$response_subject
    )
    expected_type <- if (identical(fit$margin_dist$type, "Discrete")) {
      "discrete_rectangle"
    } else {
      "continuous_density"
    }
    expect_identical(likelihood$likelihood_type, expected_type, info = fixture)

    evaluate <- function(beta) {
      eta_beta <- gamlss.longitudinal:::calc_eta(
        beta, fit$model_matrix, fit$margin_dist,
        gamlss.longitudinal:::get_copula_dist(fit$copula_dist)$copula_link,
        fit$par_s
      )
      as.numeric(gamlss.longitudinal:::calc_likelihood_minimal(
        eta_beta$eta_inv, mm = fit$model_matrix$x,
        margin_dist = fit$margin_dist, copula_dist = fit$copula_dist,
        response = fit$response, response_margin = fit$response_margin,
        response_subject = fit$response_subject, calc_margin_deriv = FALSE
      )$log_lik["joint"])
    }
    base_value <- evaluate(fit$par)
    direct_gradient <- vapply(seq_along(fit$par), function(index) {
      step <- 1e-6 * max(1, abs(fit$par[[index]]))
      plus <- minus <- fit$par
      plus[[index]] <- plus[[index]] + step
      minus[[index]] <- minus[[index]] - step
      plus_value <- evaluate(plus)
      minus_value <- evaluate(minus)
      if (is.finite(plus_value) && is.finite(minus_value)) {
        (plus_value - minus_value) / (2 * step)
      } else if (is.finite(plus_value) && is.finite(base_value)) {
        (plus_value - base_value) / step
      } else if (is.finite(minus_value) && is.finite(base_value)) {
        (base_value - minus_value) / step
      } else {
        NA_real_
      }
    }, numeric(1L))
    expect_true(all(is.finite(direct_gradient)), info = paste(fixture, "direct gradient"))
  }
})

test_that("stationarity failures remain visible and nonblocking", {
  H <- -diag(c(2, 3))
  dimnames(H) <- list(c("mu.x", "sigma.x"), c("mu.x", "sigma.x"))
  expect_warning(
    failed <- gamlss.longitudinal:::.gl_validate_hessian_inference(
      H, parameter_names = colnames(H),
      gradient = list(gradient = c(1, 1), steps = c(1, 1), scaled_max = 1)
    ),
    class = "gamlss_longitudinal_stationarity_warning"
  )
  expect_identical(failed$hessian_diagnostics$status, "available_with_warnings")
  expect_identical(failed$hessian_diagnostics$stationarity_status, "failed")
  expect_equal(failed$hessian_diagnostics$stationarity_scaled_gradient_max, 1)
  expect_equal(failed$hessian_diagnostics$stationarity_tolerance, 0.25)
  expect_identical(
    failed$hessian_diagnostics$stationarity_warning_code,
    "fitted_gradient_too_large"
  )
  expect_true(all(is.finite(failed$vcov)))

  expect_warning(
    unchecked <- gamlss.longitudinal:::.gl_validate_hessian_inference(
      H, parameter_names = colnames(H), gradient = NULL
    ),
    class = "gamlss_longitudinal_stationarity_warning"
  )
  expect_identical(unchecked$hessian_diagnostics$stationarity_status, "not_checked")
  expect_identical(
    unchecked$hessian_diagnostics$stationarity_warning_code,
    "fitted_gradient_not_checked"
  )
})
