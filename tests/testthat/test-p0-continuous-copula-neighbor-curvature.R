test_that("native df argument derivatives do not cross the Gaussian approximation join", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  limit <- .copula_t_gaussian_limit_df
  for (u in c(1e-8, 0.2, 0.7)) {
    out <- .copula_parameter_argument_mixed(u, 0.37, "t", 0.6,
      limit + c(0.01, 1, 5), parameter = "par2")
    # The native density is exactly Gaussian throughout this branch, hence
    # its argument derivatives are constant in df even next to the join.
    expect_equal(out$u1, rep(0, 3), tolerance = 0)
    expect_equal(out$u2, rep(0, 3), tolerance = 0)
  }
  expect_error(.copula_parameter_argument_mixed(0.2, 0.37, "t", 0.6,
    limit, parameter = "par2"), "nondifferentiable")
  x <- limit - c(0.01, 1, 500, .Machine$double.eps * 2^floor(log2(limit)))
  domain <- .copula_t_df_derivative_domain(x)
  seen <- numeric()
  polynomial <- function(p) {seen <<- c(seen, p); (p / limit)^2}
  expect_true(all(x <= domain$upper))
  out <- .copula_bounded_log_derivatives(polynomial, x, rep(1000, length(x)),
    domain$lower, domain$upper)
  expect_true(all(seen < limit))
  expect_true(all(seen >= 2.000001))
  expect_equal(out$first * limit, 2 * x / limit, tolerance = 1e-6)
  expect_equal(out$second * limit^2, rep(2, length(x)), tolerance = 1e-6)
  seen <- numeric()
  first <- .copula_bounded_first_difference(polynomial, x, rep(10, length(x)),
    domain$lower, domain$upper)
  expect_true(all(seen < limit))
  expect_equal(first * limit, 2 * x / limit, tolerance = 1e-6)
})

test_that("argument-parameter curvature stays at dependence boundaries", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  settings <- list(N = c(-0.999999, -0.8, 0.6, 0.999999),
    C = c(0, 1e-6, 0.5, 2), F = c(-5, -0.001, 0, 0.001, 5),
    G = c(1, 1 + 1e-6, 1.5, 5), J = c(1, 1 + 1e-6, 1.5, 5),
    t = c(-0.999999, -0.8, 0.6, 0.999999))
  for (backend in c("native", "vinecopula")) {
  options(gamlss.longitudinal.copula_backend = backend)
  for (family in names(settings)) for (theta in settings[[family]]) {
    if (backend == "vinecopula" &&
        ((family %in% c("C", "F") && theta == 0) || (family == "J" && theta == 1))) {
      expect_error(.copula_parameter_argument_mixed(0.2, 0.37, family, theta, 5),
        "parameter|Parameter")
      next
    }
    for (u in c(1e-8, 0.2, 0.7, 1 - 1e-6)) {
      ref <- copula_curvature_reference(u, 0.37, family, theta, 5)
      actual <- .copula_parameter_argument_mixed(u, 0.37, family, theta, 5)
      for (name in c("u", "v")) {
        expected <- ref[[paste0("density_", name, "theta")]]
        value <- actual[[if (name == "u") "u1" else "u2"]]
        expect_equal(value / max(1, abs(expected)),
          expected / max(1, abs(expected)), tolerance = 1e-4,
          info = paste(backend, family, theta, u, name))
      }
    }
  }
  }
})

test_that("df-argument curvature uses the actual df near the lower boundary", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (nu in c(2.000001, 2.00001, 5, 40)) for (u in c(1e-8, 0.2, 0.7)) {
      score_density <- function(df, arg) {
        ref <- copula_curvature_reference(u, 0.37, "t", 0.6, df)
        ref[["density"]] * ref[[paste0("log_", arg)]]
      }
      h <- 1e-4 * nu
      actual <- .copula_parameter_argument_mixed(u, 0.37, "t", 0.6, nu,
        parameter = "par2")
      for (arg in c("u", "v")) {
        # Richardson-extrapolated differences of the independent expression;
        # the reference's mathematical df domain extends below backend clamp.
        coarse <- (score_density(nu + h, arg) - score_density(nu - h, arg)) / (2 * h)
        fine <- (score_density(nu + h / 2, arg) - score_density(nu - h / 2, arg)) / h
        ref <- (4 * fine - coarse) / 3
        value <- actual[[if (arg == "u") "u1" else "u2"]]
        expect_equal(value / max(1, abs(ref)), ref / max(1, abs(ref)),
          tolerance = 2e-4, info = paste(backend, nu, u, arg))
      }
    }
  }
})

test_that("backend clamp and excluded pair semantics remain explicit", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (family in c("N", "C", "F", "G", "J", "t")) {
      theta <- if (family %in% c("N", "t")) 0.6 else 2
      for (u in c(1e-15, 1 - 1e-15)) {
        clamped <- min(max(u, 1e-12), 1 - 1e-12)
        ref <- copula_curvature_reference(clamped, 0.37, family, theta, 5)[["density_uv"]]
        value <- .copula_density_mixed_arguments(u, 0.37, family, theta, 5)
        expect_equal(value / max(1, abs(ref)), ref / max(1, abs(ref)),
          tolerance = 2e-6, info = paste(backend, family, u))
      }
      excluded <- .copula_hessian_pair_derivatives(c(0.2, NA), c(0.37, NA),
        family, c(theta, NA), c(5, NA), c(TRUE, FALSE), family == "t")
      expect_identical(excluded$c_val[2L], 1)
      expect_identical(excluded$d2cdu1u2[2L], 0)
      expect_identical(excluded$d2cdthu1[2L], 0)
    }
  }
})

test_that("native argument curvature matches symbolic reference in tails", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  settings <- list(N = 0.6, C = 2, F = 5, G = 1.5, J = 1.5, t = 0.6)
  for (family in names(settings)) for (u in c(1e-8, 0.2, 1 - 1e-6)) {
    theta <- settings[[family]]
    ref <- copula_curvature_reference(u, 0.37, family, theta, 5)
    actual <- .copula_hessian_pair_derivatives(u, 0.37, family, theta, 5,
      TRUE, family == "t")$d2cdu1_2
    expected <- ref[["density_uu"]]
    expect_equal(actual / max(1, abs(expected)),
      expected / max(1, abs(expected)), tolerance = 1e-4,
      info = paste(family, theta, u))
  }
})
