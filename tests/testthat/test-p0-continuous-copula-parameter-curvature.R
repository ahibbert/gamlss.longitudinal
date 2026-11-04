test_that("parameter-only Hessian curvature matches independent densities", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  settings <- list(N = c(-0.8, 0, 0.6, 0.95), C = c(0.01, 0.5, 2),
    F = c(-5, -0.01, 0.01, 5), G = c(1, 1.000001, 1.5, 5),
    J = c(1, 1.000001, 1.5, 5), t = c(-0.8, 0, 0.6, 0.95))
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (family in names(settings)) for (theta in settings[[family]]) {
      if (backend == "vinecopula" && family == "J" && theta == 1) next
      for (u in c(1e-8, 0.2, 1 - 1e-6)) {
        ref <- copula_curvature_reference(u, 0.37, family, theta, 5)
        score <- ref[["density_theta"]] / ref[["density"]]
        curvature <- ref[["density_thetatheta"]] / ref[["density"]] - score^2
        actual <- .copula_hessian_pair_derivatives(u, 0.37, family, theta, 5,
          TRUE, family == "t")
        expect_equal(actual$dldth / max(1, abs(score)), score / max(1, abs(score)),
          tolerance = 1e-4, info = paste(backend, family, theta, u, "score"))
        expect_equal(actual$d2ldth2_pair / max(1, abs(curvature)),
          curvature / max(1, abs(curvature)), tolerance = 1e-4,
          info = paste(backend, family, theta, u, "curvature"))
      }
    }
  }
})

test_that("Student-t df curvature has an independent reference at its floor", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (nu in c(2.000001, 2.00001, 5, 40)) for (u in c(1e-8, 0.2, 0.7)) {
      ref <- function(df) {
        r <- copula_curvature_reference(u, 0.37, "t", 0.6, df)
        c(logdensity = log(r[["density"]]), theta = r[["density_theta"]] / r[["density"]])
      }
      h <- 0.003 * nu
      differentiate <- function(step) {
        a <- ref(nu - step); b <- ref(nu); c <- ref(nu + step)
        c(score = (c[[1L]] - a[[1L]]) / (2 * step),
          curvature = (c[[1L]] - 2 * b[[1L]] + a[[1L]]) / step^2,
          cross = (c[[2L]] - a[[2L]]) / (2 * step))
      }
      coarse <- differentiate(h); fine <- differentiate(h / 2)
      expected <- (4 * fine - coarse) / 3
      actual <- .copula_hessian_pair_derivatives(u, 0.37, "t", 0.6, nu, TRUE, TRUE)
      value <- c(score = actual$zeta$dldz_pair,
        curvature = actual$zeta$d2ldz2_pair, cross = actual$zeta$d2ldthdz_pair)
      scale <- pmax(1, abs(expected))
      expect_equal(value / scale, expected / scale, tolerance = 2e-4,
        info = paste(backend, nu, u))
    }
  }
})

test_that("unrepresentable density and the native t approximation join are explicit", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    expect_error(.copula_hessian_pair_derivatives(1e-12, 0.9, "N", 0.999999,
      0, TRUE, FALSE), "underflowed")
  }
  options(gamlss.longitudinal.copula_backend = "native")
  expect_error(.copula_hessian_pair_derivatives(0.2, 0.37, "t", 0.6, 1e6,
    TRUE, TRUE), "nondifferentiable")
  out <- .copula_hessian_pair_derivatives(0.2, 0.37, "t", 0.6, 1e6 + 1000,
    TRUE, TRUE)
  expect_equal(out$zeta$dldz_pair, 0)
  expect_equal(out$zeta$d2ldz2_pair, 0)
  expect_equal(out$zeta$d2ldthdz_pair, 0)
  ref <- copula_curvature_reference(0.2, 0.37, "N", 0.6)
  expect_equal(out$dldth, ref[["density_theta"]] / ref[["density"]], tolerance = 1e-8)
})
