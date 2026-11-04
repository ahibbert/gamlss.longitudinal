# Fixed synthetic designs only: no optimizer, saved fit, simulation study or RDS read.
joint_gaussian_reference_fixture <- function() {
  m <- 12L
  x <- seq(-1, 1, length.out = m)
  z <- (-1)^seq_len(m) + 0.7 * x
  dat <- data.frame(x = rep(x, each = 2L), z = rep(z, each = 2L))
  environment <- new.env(parent = baseenv())
  environment$s <- mgcv::s
  term <- 's(x, bs = "cr", k = 4)'
  smooth <- .gl_build_smooth_model_matrices(stats::as.formula(paste("~", term)), dat, environment)
  term <- names(smooth)[1L]
  B <- smooth[[term]]
  subject_design <- cbind(intercept = 1, z = z, B[seq(1L, 2L * m, by = 2L), , drop = FALSE])
  Q <- qr.Q(qr(subject_design), complete = TRUE)
  rank <- ncol(subject_design)
  residual <- Q[, rank + 1:2, drop = FALSE] * sqrt(m)
  response <- as.vector(t(residual))
  fixed <- c(mu.intercept = 0, mu.z = 0, sigma.intercept = 0, theta.intercept = 0)
  b <- stats::setNames(rep(0, ncol(B)), paste("mu", term, seq_len(ncol(B)), sep = "."))
  empty <- list(mu = list(), sigma = list(), theta = list())
  par_s <- empty
  par_s$mu[[term]] <- b
  lambda <- empty
  lambda$mu[[term]] <- 0.7
  matrices <- list(
    x = list(mu = cbind(intercept = 1, z = dat$z),
      sigma = matrix(1, 2L * m, 1L, dimnames = list(NULL, "intercept")),
      theta = matrix(1, m, 1L, dimnames = list(NULL, "intercept"))),
    s = list(mu = smooth, sigma = list(), theta = list()))
  object <- structure(list(par = fixed, par_s = par_s, lambda_s = lambda,
    model_matrix = matrices, response = response,
    response_subject = rep(seq_len(m), each = 2L), response_margin = rep(1:2, m),
    margin_dist = gamlss.dist::NO(), copula_dist = "N",
    convergence = list(converged = TRUE, stop_reason = "synthetic_exact_stationary_point"),
    likelihood_contract = list(objective = "ordinary"), missingness = list(objective = "ordinary")),
    class = "gamlss.longitudinal")
  # Independent Gaussian likelihood algebra: every subject has two independent
  # N(mu,1) observations, both design rows D. Residual columns are orthogonal to
  # D and each other, with squared norm m: all nuisance cross-blocks vanish.
  q <- c(fixed, b)
  I <- P <- matrix(0, length(q), length(q), dimnames = list(names(q), names(q)))
  location <- c("mu.intercept", "mu.z", names(b))
  I[location, location] <- 2 * crossprod(subject_design)
  I["sigma.intercept", "sigma.intercept"] <- 4 * m
  I["theta.intercept", "theta.intercept"] <- m
  P[names(b), names(b)] <- 0.7 * attr(B, "penalty")
  list(object = object, information = I, penalty = P, covariance = solve(I + P),
       term = term, data = dat, q = q)
}

test_that("joint inverse information matches independent Gaussian block algebra", {
  x <- seq(-2, 2)
  B <- x^2 - mean(x^2) + 0.5 * x
  D <- cbind(intercept = 1, x = x, smooth = B)
  sigma2 <- 2
  lambda <- 3
  I <- crossprod(D) / sigma2
  P <- diag(c(0, 0, lambda))
  dimnames(P) <- dimnames(I)
  y <- c(-1, 2, 4, 1, 3)
  q <- as.numeric(solve(I + P, crossprod(D, y) / sigma2))
  names(q) <- colnames(D)
  raw_score <- as.numeric(crossprod(D, y - D %*% q) / sigma2)
  names(raw_score) <- names(q)
  # Closed-form block inverse, independent of the implementation's solver.
  F_inv <- diag(c(sigma2 / length(x), sigma2 / sum(x^2)))
  cross <- matrix(c(sum(B), sum(x * B)) / sigma2, ncol = 1L)
  variance_b <- 1 / (sum(B^2) / sigma2 + lambda - as.numeric(t(cross) %*% F_inv %*% cross))
  V <- rbind(cbind(F_inv + F_inv %*% cross %*% t(cross) %*% F_inv * variance_b,
                  -F_inv %*% cross * variance_b),
             cbind(-t(cross) %*% F_inv * variance_b, variance_b))
  dimnames(V) <- dimnames(I)
  result <- .gl_joint_solve(-I, P, q, raw_score, inference_control())
  expect_equal(result$vcov, V, tolerance = 1e-12)
  expect_gt(abs(result$vcov["x", "smooth"]), 0.001)
  expect_gt(max(abs(raw_score)), 0.1)
  expect_lt(max(abs(result$hessian_diagnostics$penalized_gradient)), 1e-12)
  expect_identical(result$hessian_diagnostics$stationarity_status, "passed")
  V_frequentist <- V %*% I %*% V
  expect_equal(V - V_frequentist, V %*% P %*% V, tolerance = 1e-12)
  expect_gt(max(abs(V - V_frequentist)), 0.001)
  contrast <- c(1, 0.2, 0.7)
  expect_equal(as.numeric(t(contrast) %*% result$vcov %*% contrast),
               as.numeric(t(contrast) %*% V %*% contrast), tolerance = 1e-12)
})

test_that("public joint covariance matches a full Gaussian likelihood reference", {
  if (!requireNamespace("mgcv", quietly = TRUE) ||
      !requireNamespace("VineCopula", quietly = TRUE)) {
    stop("Joint inference verification requires both mgcv and VineCopula; backend coverage cannot be skipped.")
  }
  fixture <- joint_gaussian_reference_fixture()
  original <- fixture$object
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    for (method in c("analytical_only", "numderiv")) {
      result <- vcov.gamlss.longitudinal(fixture$object, coefficients = "joint",
        method = method, details = TRUE, progress = FALSE)
      expect_equal(result$likelihood_hessian, -fixture$information, tolerance = 2e-4)
      expect_equal(result$vcov$overall, fixture$covariance, tolerance = 2e-4)
      expect_equal(result$penalty, fixture$penalty, tolerance = 1e-12)
      expect_identical(result$coefficient_map$name, names(fixture$q))
      expect_identical(result$inference_contract$contract_id, "joint_penalized_inverse")
      expect_identical(result$covariance_type, "penalized_inverse")
      expect_identical(result$hessian_diagnostics$stationarity_status, "passed")
      expect_equal(result$vcov$fixed_smooth,
        fixture$covariance[names(original$par), names(original$par_s$mu[[fixture$term]]), drop = FALSE],
        tolerance = 2e-4)
      expect_identical(fixture$object, original)
    }
  }
})

test_that("joint coordinate names determine covariance despite fixed-vector permutation", {
  fixture <- joint_gaussian_reference_fixture()
  permuted <- fixture$object
  permuted$par <- rev(permuted$par)
  a <- .gl_joint_design(fixture$object, inference_control())
  b <- .gl_joint_design(permuted, inference_control())
  expect_identical(b$map$name[seq_along(permuted$par)], names(permuted$par))
  expect_equal(b$penalty[names(a$q), names(a$q)], a$penalty)
  out <- vcov.gamlss.longitudinal(permuted, coefficients = "joint", method = "numderiv",
                                 details = TRUE, progress = FALSE)
  expect_equal(out$vcov$overall[names(a$q), names(a$q)], fixture$covariance, tolerance = 2e-4)
})

test_that("joint scope excludes unsupported objectives and unverifiable smooth metadata", {
  fixture <- joint_gaussian_reference_fixture()
  fit <- fixture$object
  request <- function(x, ...) vcov.gamlss.longitudinal(x, coefficients = "joint", progress = FALSE, ...)
  expect_error(request(fit, method = "sandwich"), class = "gamlss_longitudinal_inference_unavailable")
  expect_error(request(fit, par = list(par = fit$par, par_s = fit$par_s)), "fitted coefficients only")
  segmented <- fit
  segmented$likelihood_contract$objective <- "segmented"
  expect_error(request(segmented), "segmented or intermittent")
  invalid <- fit
  invalid$convergence$converged <- FALSE
  expect_error(request(invalid), class = "gamlss.longitudinal_nonconvergence_error")
  invalid <- fit
  attr(invalid$model_matrix$s$mu[[fixture$term]], "gl_smooth_spec") <- NULL
  expect_error(request(invalid), "retained single-component")
  invalid <- fit
  attr(invalid$model_matrix$s$mu[[fixture$term]], "gl_smooth_penalty_count") <- 2L
  expect_error(request(invalid), "retained single-component")
  invalid <- fit
  invalid$lambda_s$mu[[fixture$term]] <- -1
  expect_error(request(invalid), "penalty/lambda")
  invalid <- fit
  names(invalid$par_s$mu[[fixture$term]])[1] <- "mu.unmapped.1"
  expect_error(request(invalid), "Unmatched smooth basis names")
})

test_that("joint SPD and rank checks cannot be bypassed by positive inverse diagonals", {
  q <- c(a = 0, b = 0, c = 0)
  A <- matrix(110 / 3, 3, 3) - diag(100, 3)
  dimnames(A) <- list(names(q), names(q))
  P <- A * 0
  expect_true(all(diag(solve(A)) > 0))
  expect_error(.gl_joint_solve(-A, P, q, q, inference_control()), "positive definite")
  q <- c(a = 0, b = 0)
  A <- matrix(c(1, 1 - 1e-12, 1 - 1e-12, 1), 2, 2, dimnames = list(names(q), names(q)))
  expect_error(.gl_joint_solve(-A, A * 0, q, q, inference_control(rank_tol = 1e-8)), "rank deficient")
})

test_that("stored fitted smooth predicts the same coordinates on reordered and new rows", {
  fixture <- joint_gaussian_reference_fixture()
  B <- fixture$object$model_matrix$s$mu[[fixture$term]]
  predict_basis <- function(data) .gl_joint_smooth_prediction_matrix(fixture$object, "mu", fixture$term, data)
  expect_equal(as.numeric(predict_basis(fixture$data)), as.numeric(B), tolerance = 1e-10)
  subset <- c(9, 1, 22, 4)
  expect_equal(as.numeric(predict_basis(fixture$data[subset, , drop = FALSE])),
               as.numeric(B[subset, , drop = FALSE]), tolerance = 1e-10)
  expect_equal(as.numeric(predict_basis(fixture$data[1, , drop = FALSE])), as.numeric(B[1, ]), tolerance = 1e-10)
  shifted <- data.frame(x = c(-0.85, -0.2, 0.3, 0.9))
  together <- predict_basis(shifted)
  individually <- do.call(rbind, lapply(seq_len(nrow(shifted)), function(i) predict_basis(shifted[i, , drop = FALSE])))
  expect_equal(as.numeric(together), as.numeric(individually), tolerance = 1e-10)
})

test_that("joint requests cannot reuse fixed caches or be blocked by stale fixed failures", {
  fixture <- joint_gaussian_reference_fixture()
  fit <- fixture$object
  for (legacy_cache in list(NULL, diag(length(fit$par)), 1, "legacy")) {
    fit$vcov <- legacy_cache
    expect_false(.can_use_cached_vcov(fit))
    expect_false(.can_use_cached_vcov(fit, extra_args = list(coefficients = "joint")))
  }
  fit$vcov <- list(coefficients = "fixed", vcov = list(overall = diag(length(fit$par))))
  expect_false(.can_use_cached_vcov(fit, extra_args = list(coefficients = "joint")))
  fit$vcov$coefficients <- "joint"
  expect_false(.can_use_cached_vcov(fit))
  fit$vcov_meta <- list(inference_status = "unavailable", method = "analytical",
    hessian_diagnostics = list(status = "unavailable", failure_codes = "covariance_nonfinite"))
  out <- .resolve_vcov(fit, extra_args = list(coefficients = "joint", method = "numderiv", progress = FALSE))
  expect_identical(out$coefficients, "joint")
  expect_equal(out$vcov$overall, fixture$covariance, tolerance = 2e-4)
})

test_that("default fixed behavior and no-smooth joint covariance agree", {
  fixture <- joint_gaussian_reference_fixture()
  fit <- fixture$object
  fit$model_matrix$s$mu <- list()
  fit$par_s$mu <- list()
  fit$lambda_s$mu <- list()
  original <- fit
  default <- vcov.gamlss.longitudinal(fit, method = "analytical_only", progress = FALSE)
  explicit <- vcov.gamlss.longitudinal(fit, coefficients = "fixed", method = "analytical_only", progress = FALSE)
  joint <- vcov.gamlss.longitudinal(fit, coefficients = "joint", method = "analytical_only", progress = FALSE)
  expect_identical(default, explicit)
  expect_equal(as.numeric(joint), as.numeric(default), tolerance = 1e-10)
  expect_identical(rownames(default), names(fit$par))
  expect_identical(fit, original)
})
