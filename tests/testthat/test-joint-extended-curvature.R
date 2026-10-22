# Fit-free joint-coordinate checks. Source the separately pinned independent
# helper-continuous-model-reference.R before this file. No optimizer or model IO.
joint_extended_fixture <- function(family = "N", tail = FALSE) {
  object <- curvature_model_fixture(family, tail)
  environment <- new.env(parent = baseenv())
  environment$s <- mgcv::s
  smooth <- .gl_build_smooth_model_matrices(~ s(x, bs = "cr", k = 4),
                                            object$dataset, environment)
  term <- names(smooth)[[1L]]
  B <- smooth[[term]]
  object$model_matrix$s$mu <- smooth
  object$par_s$mu[[term]] <- stats::setNames(seq_len(ncol(B)) * 0.015,
                                           paste("mu", term, seq_len(ncol(B)), sep = "."))
  object$lambda_s$mu[[term]] <- 0.7
  object
}

test_that("joint augmented smooth curvature matches six defining copula likelihoods", {
  stopifnot(requireNamespace("numDeriv", quietly = TRUE),
            requireNamespace("VineCopula", quietly = TRUE))
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N", "C", "F", "G", "J", "t")) for (tail in c(FALSE, TRUE)) {
    object <- joint_extended_fixture(family, tail)
    original <- object
    joint <- .gl_joint_design(object, inference_control())
    reference <- object
    reference$par <- joint$q
    B <- object$model_matrix$s$mu[[1L]]
    reference$synthetic_smooth_covariates <- B[seq(1L, nrow(B), by = 3L), , drop = FALSE]
    colnames(reference$synthetic_smooth_covariates) <- names(object$par_s$mu[[1L]])
    expected <- curvature_reference_hessian(reference)
    dimnames(expected) <- list(names(joint$q), names(joint$q))
    smooth_ids <- names(object$par_s$mu[[1L]])
    expect_equal(as.numeric(joint$penalty[smooth_ids, smooth_ids]),
      as.numeric(0.7 * attr(B, "penalty")), tolerance = 1e-12)
    for (backend in c("native", "vinecopula")) {
      options(gamlss.longitudinal.copula_backend = backend)
      actual <- calc_analytical_hessian(joint$object, progress = FALSE)
      error <- sqrt(sum((actual - expected)^2)) / sqrt(sum(expected^2))
      expect_true(is.finite(error) && error < 2e-4,
                  info = paste(family, backend, "tail", tail, "relative error", error))
      evaluate <- .gl_joint_likelihood_evaluator(joint$object)
      expect_equal(evaluate(joint$q)$loglik,
        sum(curvature_reference_clusters(joint$q, reference)), tolerance = 2e-8)
      expect_identical(object, original)
    }
  }
})

joint_nbi_reference <- function(q, object, B) {
  # NBI size=1/sigma, mean=mu; exact Clayton CDF rectangles. This reference
  # calls no package likelihood, link, pair-cache, derivative or rectangle helper.
  mu <- exp(q[["mu.intercept"]] + q[["mu.x"]] * object$dataset$x +
              as.numeric(B %*% q[colnames(B)]))
  sigma <- exp(q[["sigma.intercept"]])
  theta <- exp(q[["theta.intercept"]])
  y <- object$response
  upper <- stats::pnbinom(y, size = 1 / sigma, mu = mu)
  lower <- stats::pnbinom(y - 1, size = 1 / sigma, mu = mu)
  mass <- stats::dnbinom(y, size = 1 / sigma, mu = mu)
  cdf <- function(u, v) {
    out <- numeric(length(u)); positive <- u > 0 & v > 0
    out[positive] <- (u[positive]^(-theta) + v[positive]^(-theta) - 1)^(-1 / theta)
    out
  }
  left <- which(object$response_margin < 3L); right <- left + 1L
  rectangle <- cdf(upper[left], upper[right]) - cdf(lower[left], upper[right]) -
    cdf(upper[left], lower[right]) + cdf(lower[left], lower[right])
  stopifnot(all(is.finite(rectangle)), all(rectangle > 0), all(mass > 0))
  sum(log(mass)) + sum(log(rectangle) - log(mass[left]) - log(mass[right]))
}

test_that("joint NBI smooth curvature matches independent Clayton rectangle likelihood", {
  object <- joint_extended_fixture("C")
  object$margin_dist <- gamlss.dist::NBI()
  object$par[c("mu.intercept", "mu.x", "sigma.intercept", "theta.intercept")] <-
    c(log(3), 0.15, log(0.5), log(1.2))
  object$response <- rep(c(0, 1, 3, 2, 8, 5, 1, 4, 2, 6, 3, 0), length.out = length(object$response))
  object$dataset$response <- object$response
  original <- object
  joint <- .gl_joint_design(object, inference_control())
  B <- object$model_matrix$s$mu[[1L]]
  colnames(B) <- names(object$par_s$mu[[1L]])
  objective <- function(q) joint_nbi_reference(q, object, B)
  # Differentiate at zero offsets so eps is an absolute link-coordinate step.
  # Fractional d=0.001 at small smooth coefficients suffered cancellation;
  # these three steps were independently checked before changing this test.
  # Every reference must agree: no best-step selection or tolerance relaxation.
  absolute_reference <- function(eps) {
    zero <- stats::setNames(numeric(length(joint$q)), names(joint$q))
    evaluated <- 0L
    largest_step <- numeric(length(joint$q))
    shifted <- function(delta) {
      evaluated <<- evaluated + 1L
      stopifnot(evaluated <= 250L, is.numeric(delta),
        length(delta) == length(joint$q), all(is.finite(delta)),
        all(abs(delta) <= eps + 1e-12))
      largest_step <<- pmax(largest_step, abs(delta))
      objective(joint$q + delta)
    }
    H <- numDeriv::hessian(shifted, zero, method = "Richardson",
      method.args = list(eps = eps, d = 1e-3, r = 4, v = 2))
    stopifnot(all(is.finite(H)), sum(H^2) > 0,
      max(abs(largest_step - eps)) < 1e-12)
    dimnames(H) <- list(names(joint$q), names(joint$q))
    H
  }
  steps <- c(0.005, 0.01, 0.02)
  references <- lapply(steps, absolute_reference)
  relative <- function(a, b) sqrt(sum((a - b)^2)) / sqrt(sum(b^2))
  pairs <- utils::combn(seq_along(references), 2L)
  for (i in seq_len(ncol(pairs))) {
    a <- pairs[1L, i]; b <- pairs[2L, i]
    error <- relative(references[[a]], references[[b]])
    expect_true(is.finite(error) && error < 2e-3,
      info = paste("NBI Clayton reference stability", steps[a], steps[b], error))
  }
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    evaluate <- .gl_joint_likelihood_evaluator(joint$object)
    expect_equal(evaluate(joint$q)$loglik, objective(joint$q), tolerance = 1e-9)
    actual <- calc_analytical_hessian(joint$object, progress = FALSE)
    for (i in seq_along(references)) {
      error <- relative(actual, references[[i]])
      expect_true(is.finite(error) && error < 2e-3,
        info = paste(backend, "NBI Clayton absolute step", steps[i], "relative error", error))
    }
    expect_identical(object, original)
  }
})
