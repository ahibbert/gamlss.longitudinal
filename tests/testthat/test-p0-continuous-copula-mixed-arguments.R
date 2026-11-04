test_that("all six mixed density derivatives match defining symbolic densities", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  parameters <- list(N = c(-0.8, 0, 0.6, 0.95), C = c(1e-6, 0.5, 2, 10),
    F = c(-20, -5, -0.001, 0.001, 5, 20),
    G = c(1, 1.000001, 1.5, 5), J = c(1, 1.000001, 1.5, 5),
    t = c(-0.8, 0, 0.6, 0.95))
  pairs <- rbind(cbind(c(1e-10, 1e-6, 0.2, 0.7, 1 - 1e-8), 0.37),
    c(1e-8, 2e-8), c(1 - 1e-8, 1 - 2e-8))
  for (backend in c("native", "vinecopula")) {
  options(gamlss.longitudinal.copula_backend = backend)
  for (family in names(parameters)) for (theta in parameters[[family]]) {
    if (backend == "vinecopula" && family == "J" && theta == 1) {
      expect_error(.copula_density_mixed_arguments(0.2, 0.37, family, theta, 5),
        "Joe.*interval")
      next
    }
    for (nu in if (family == "t") c(2.000001, 5, 40) else 5) {
      for (i in seq_len(nrow(pairs))) {
        u <- pairs[i, 1L]
        v <- pairs[i, 2L]
        ref <- copula_curvature_reference(u, v, family, theta, nu)
        actual <- .copula_density_mixed_arguments(u, v, family, theta, nu)
        label <- paste(backend, family, theta, nu, u, v)
        expect_true(is.finite(ref[["density_uv"]]), info = label)
        expect_true(is.finite(actual), info = label)
        scale <- max(.Machine$double.xmin, abs(ref[["density_uv"]]))
        expect_equal(actual / scale,
          ref[["density_uv"]] / scale,
          tolerance = 2e-6, info = label)
      }
    }
  }
  }
})

test_that("symbolic references agree with convergent transformed differences", {
  settings <- list(N = 0.6, C = 2, F = 5, G = 1.5, J = 1.5, t = 0.6)
  for (family in names(settings)) for (u in c(1e-8, 0.2, 1 - 1e-6)) {
    v <- 0.37
    theta <- settings[[family]]
    nu <- 5
    # Numerical density derivatives in logit coordinates avoid a displaced
    # tail centre. The mixed Jacobian has no cross-coordinate second term.
    x <- stats::qlogis(u)
    y <- stats::qlogis(v)
    f <- function(a, b) copula_curvature_reference(
      stats::plogis(a), stats::plogis(b), family, theta, nu)[["density"]]
    numeric_mixed <- function(h) {
      if (family %in% c("F", "J")) {
        # For smooth densities at u near zero, a logit u perturbation gives
        # a density change below double precision. Complex-step in u keeps
        # that information; central differences in v then recover c_uv.
        # These defining density formulas do not use symbolic derivatives.
        density_complex <- function(a, b) {
          if (family == "F") {
            eu <- exp(-theta * a)
            ev <- exp(-theta * b)
            return(-theta * expm1(-theta) * eu * ev /
              (expm1(-theta) + (exp(-theta * a) - 1) * (exp(-theta * b) - 1))^2)
          }
          aa <- (1 - a)^theta
          bb <- (1 - b)^theta
          s <- aa + bb - aa * bb
          (1 - a)^(theta - 1) * (1 - b)^(theta - 1) *
            s^(1 / theta - 2) * (s + theta - 1)
        }
        return((Im(density_complex(u + 1e-20i, v + h)) -
          Im(density_complex(u + 1e-20i, v - h))) / (2 * h * 1e-20))
      }
      (f(x + h, y + h) - f(x + h, y - h) -
        f(x - h, y + h) + f(x - h, y - h)) /
          (4 * h^2 * u * (1 - u) * v * (1 - v))
    }
    ref <- copula_curvature_reference(u, v, family, theta, nu)[["density_uv"]]
    coarse <- numeric_mixed(0.01)
    fine <- numeric_mixed(0.003)
    scale <- max(1, abs(ref))
    expect_true(abs(fine - ref) / scale < 0.002, info = paste(family, u))
    expect_true(abs(fine - ref) / scale <=
      abs(coarse - ref) / scale + 2e-5, info = paste(family, u))
  }
})

test_that("mixed argument derivatives preserve recycling, symmetry and masking", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N", "C", "F", "G", "J", "t")) {
    theta <- if (family %in% c("N", "t")) 0.6 else 2
    u <- c(0.2, 1e-8, 0.7)
    v <- c(0.37, 0.37, 0.8)
    vector <- .copula_density_mixed_arguments(u, v, family, theta, 5)
    scalar <- vapply(seq_along(u), function(i)
      .copula_density_mixed_arguments(u[i], v[i], family, theta, 5), numeric(1))
    expect_equal(vector, scalar)
    expect_equal(vector, .copula_density_mixed_arguments(v, u, family, theta, 5))
    out <- .copula_hessian_pair_derivatives(u, v, family, theta, 5,
      c(TRUE, FALSE, TRUE), family == "t")
    expect_equal(out$d2cdu1u2, replace(vector, 2L, 0))
  }
})

test_that("VineCopula and native mixed derivatives agree on their common domain", {
  expect_true(requireNamespace("VineCopula", quietly = TRUE))
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N", "C", "F", "G", "J", "t")) {
    theta <- if (family %in% c("N", "t")) 0.6 else 2
    u <- c(0.2, 1e-8, 1 - 1e-8)
    v <- c(0.37, 0.37, 0.37)
    options(gamlss.longitudinal.copula_backend = "native")
    native <- .copula_density_mixed_arguments(u, v, family, theta, 5)
    options(gamlss.longitudinal.copula_backend = "vinecopula")
    vine <- .copula_density_mixed_arguments(u, v, family, theta, 5)
    scale <- pmax(1, abs(native))
    expect_equal(vine / scale, native / scale, tolerance = 1e-6, info = family)
  }
})
