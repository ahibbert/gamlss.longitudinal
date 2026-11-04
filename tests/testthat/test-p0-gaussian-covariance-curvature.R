# The narrow patch's original implementation/tests remain in Git history.
# All-family routing and cache checks now live in test-p0-continuous-curvature-cache.R.
test_that("the Gaussian compatibility path retains independent tail curvature", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (u in c(1e-8, 4.46758e-6, 0.2)) {
      ref <- copula_curvature_reference(u, 0.0068, "N", 0.749)[["density_uv"]]
      actual <- .gl_covariance_gaussian_pair_derivatives(
        u, 0.0068, "N", 0.749, 0, TRUE, FALSE)$d2cdu1u2
      expect_equal(actual / max(1, abs(ref)), ref / max(1, abs(ref)),
        tolerance = 1e-6)
    }
  }
})
