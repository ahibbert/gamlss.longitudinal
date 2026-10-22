# Synthetic coefficient/basis algebra only. No fit, optimizer, RDS or simulations.
# The covariance below is a declared synthetic SPD matrix, not a model Hessian.
contrast_test_bundle <- function(q, map, V, penalty = V * 0, lambda = list()) {
  list(coefficients = "joint", joint_inference_version = 1L,
    covariance_type = "penalized_inverse", coefficient_estimates = q,
    coefficient_map = map, vcov = list(overall = V), penalty = penalty,
    smoothing_parameters = lambda, conditioning = "fitted lambda and constrained basis",
    hessian_diagnostics = list(status = "available", stationarity_status = "passed",
      failure_codes = character(), blocking_failure_codes = character(), warning_codes = character()),
    inference_contract = list(contract_id = "joint_penalized_inverse", validity_status = "available",
      observed_failures = character(), omitted_uncertainty = "smoothing parameter uncertainty"))
}

contrast_basis_fixture <- function() {
  if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv is required; retained-basis tests cannot skip.")
  dat <- data.frame(s1 = seq(0.03, 0.82, length.out = 24L))
  env <- new.env(parent = baseenv())
  env$s <- mgcv::s
  formulas <- list(mu = ~ s(s1, bs = "cr", k = 4) + s(s1, bs = "ps", k = 5),
    sigma = ~ s(s1, bs = "cr", k = 4), theta = ~ s(s1, bs = "cr", k = 4))
  smooth <- par_s <- lambda <- list(mu = list(), sigma = list(), theta = list())
  X <- list()
  for (p in names(smooth)) {
    d <- if (p == "theta") dat[seq(1L, 23L, by = 2L), , drop = FALSE] else dat
    X[[p]] <- matrix(1, nrow(d), 1L, dimnames = list(NULL, "intercept"))
    smooth[[p]] <- .gl_build_smooth_model_matrices(formulas[[p]], d, env)
    for (term in names(smooth[[p]])) {
      k <- ncol(smooth[[p]][[term]])
      par_s[[p]][[term]] <- stats::setNames(seq_len(k) / 10, paste(p, term, seq_len(k), sep = "."))
      lambda[[p]][[term]] <- 0.7
    }
  }
  object <- structure(list(par = c(mu.intercept = 2, sigma.intercept = -0.3, theta.intercept = 0.2),
    par_s = par_s, lambda_s = lambda, model_matrix = list(x = X, s = smooth),
    convergence = list(converged = TRUE, stop_reason = "synthetic_contrast_mapping_only"),
    likelihood_contract = list(objective = "ordinary"), missingness = list(objective = "ordinary")),
    class = "gamlss.longitudinal")
  design <- .gl_joint_design(object, inference_control())
  k <- length(design$q)
  v <- seq_len(k) / (10 * k)
  V <- diag(seq_len(k) / k + 0.5) + tcrossprod(v)
  dimnames(V) <- list(names(design$q), names(design$q))
  list(object = object, design = design,
    joint = contrast_test_bundle(design$q, design$map, V, design$penalty, lambda))
}

test_that("named contrasts retain fixed-smooth cross covariance and align columns", {
  q <- c(mu.intercept = 2, mu.s1.1 = 3)
  map <- data.frame(name = names(q), parameter = "mu", kind = c("fixed", "smooth"),
    term = c("intercept", "s1"), basis_index = c(NA_integer_, 1L))
  V <- matrix(c(4, -1, -1, 9), 2, dimnames = list(names(q), names(q)))
  joint <- contrast_test_bundle(q, map, V)
  L <- matrix(c(1, 0.5), 1L, dimnames = list("reference", names(q)))
  out <- .gl_joint_linear_contrast(joint, L)
  expect_equal(out$table$estimate, 3.5)
  expect_equal(out$table$std_error^2, 4 + 0.25 * 9 + 2 * 0.5 * (-1))
  expect_gt(abs(out$table$std_error^2 - (4 + 0.25 * 9)), 0.5)
  expect_identical(.gl_joint_linear_contrast(joint, L[, 2:1, drop = FALSE]), out)
  expect_equal(out$table$conf_high, 3.5 + stats::qnorm(0.975) * sqrt(5.25))
  expect_identical(out$inference_contract$covariance_type, "penalized_inverse")
  expect_true(out$table$uncertainty_available)
  expect_true(out$table$clean_uncertainty_eligible)
  joint$hessian_diagnostics$status <- joint$inference_contract$validity_status <- "available_with_warnings"
  joint$hessian_diagnostics$warning_codes <- "score_warning"
  warned <- .gl_joint_linear_contrast(joint, L)
  expect_equal(warned$table$std_error, out$table$std_error)
  expect_true(warned$table$uncertainty_available)
  expect_false(warned$table$clean_uncertainty_eligible)
})

test_that("contrasts reject incomplete weights, invalid metadata and indefinite covariance", {
  f <- contrast_basis_fixture()
  q <- f$joint$coefficient_estimates
  L <- matrix(0, 1L, length(q), dimnames = list("mu", names(q)))
  L[1, 1] <- 1
  bad <- L; colnames(bad)[2] <- colnames(bad)[1]
  expect_error(.gl_joint_linear_contrast(f$joint, bad), "unique named")
  expect_error(.gl_joint_linear_contrast(f$joint, L[, -1, drop = FALSE]), "unique named")
  bad <- L; colnames(bad)[1] <- "unknown"
  expect_error(.gl_joint_linear_contrast(f$joint, bad), "unique named")
  bad <- L; dimnames(bad) <- NULL
  expect_error(.gl_joint_linear_contrast(f$joint, bad), "unique named")
  bad <- L; bad[1] <- Inf
  expect_error(.gl_joint_linear_contrast(f$joint, bad), "unique named")
  expect_error(.gl_joint_linear_contrast(f$joint, L * 0), "unique named")
  expect_error(.gl_joint_linear_contrast(f$joint, L, level = 1), "confidence level")
  bad <- f$joint; bad$hessian_diagnostics <- NULL
  expect_error(.gl_joint_linear_contrast(bad, L), "inconsistent validity")
  bad <- f$joint; bad$vcov$overall[1, 1] <- -1
  expect_error(.gl_joint_linear_contrast(bad, L), "positive definite")
  bad <- f$joint; bad$coefficient_map$name <- rev(bad$coefficient_map$name)
  expect_error(.gl_joint_linear_contrast(bad, L), "map/estimates")
})

test_that("BCPE sample reference equals independent mean fitted contributions", {
  f <- contrast_basis_fixture()
  original <- f$object
  out <- .gl_bcpe_sample_intercept_contrasts(f$object, f$joint)
  manual <- matrix(0, 3L, length(f$design$q), dimnames = dimnames(out$L))
  values <- numeric(3L)
  for (i in seq_along(out$reference$parameters)) {
    p <- out$reference$parameters[i]
    manual[i, paste0(p, ".intercept")] <- 1
    eta <- rep(f$object$par[[paste0(p, ".intercept")]], nrow(f$object$model_matrix$x[[p]]))
    for (term in names(f$object$par_s[[p]])) {
      b <- f$object$par_s[[p]][[term]]
      B <- f$object$model_matrix$s[[p]][[term]]
      eta <- eta + as.numeric(B %*% b)
      manual[i, names(b)] <- apply(B, 2L, mean)
    }
    values[i] <- mean(eta)
  }
  expect_equal(out$L, manual, tolerance = 1e-12)
  expect_equal(out$table$estimate, values, tolerance = 1e-12)
  expect_equal(out$covariance, manual %*% f$joint$vcov$overall %*% t(manual), tolerance = 1e-12)
  expect_identical(unname(out$reference$row_counts), c(24L, 24L, 12L))
  expect_length(out$reference$terms$mu, 2L)
  expect_identical(f$object, original)
})

test_that("missingness reference uses all retained smooths on the fixed population grid", {
  f <- contrast_basis_fixture()
  out <- .gl_missingness_population_intercept_contrasts(f$object, f$joint)
  grid <- data.frame(s1 = seq(0, 1, length.out = 10001L))
  manual <- matrix(0, 3L, length(f$design$q), dimnames = dimnames(out$L))
  for (i in seq_along(out$reference$parameters)) {
    p <- out$reference$parameters[i]
    manual[i, paste0(p, ".intercept")] <- 1
    for (term in names(f$object$par_s[[p]])) {
      # Independent call to the retained mgcv specification, not our mapper.
      spec <- attr(f$object$model_matrix$s[[p]][[term]], "gl_smooth_spec", exact = TRUE)
      B <- mgcv::PredictMat(spec, grid)
      manual[i, names(f$object$par_s[[p]][[term]])] <- apply(B, 2L, mean)
    }
  }
  expect_equal(out$L, manual, tolerance = 1e-12)
  expect_equal(out$table$estimate, as.numeric(manual %*% f$design$q), tolerance = 1e-12)
  expect_equal(out$covariance, manual %*% f$joint$vcov$overall %*% t(manual), tolerance = 1e-12)
  expect_identical(out$reference$grid$points, 10001L)
  expect_true(out$reference$grid$endpoints)
  expect_identical(out$reference$grid$weights, "equal")
  expect_length(out$reference$terms$mu, 2L)
  expect_gt(max(abs(out$L - .gl_bcpe_sample_intercept_contrasts(f$object, f$joint)$L)), 1e-5)
  changed <- f$object
  changed$response <- seq_len(24L) + 100
  changed$unrelated_mask <- rep(c(TRUE, FALSE), 12L)
  # Basis weights do not depend on response values; bundle origin remains the caller's responsibility.
  expect_identical(.gl_missingness_population_intercept_contrasts(changed, f$joint)$L, out$L)
})

test_that("paper helpers refuse incompatible objectives and unsupported reference bases", {
  f <- contrast_basis_fixture()
  call <- function(object, joint = f$joint) .gl_missingness_population_intercept_contrasts(object, joint)
  bad <- f$object; bad$likelihood_contract$objective <- "segmented"
  expect_error(call(bad), "Segmented/intermittent")
  bad <- f$object; bad$missingness$n_subjects_with_gaps <- 1L
  expect_error(call(bad), "Segmented/intermittent")
  bad <- f$object; bad$convergence$converged <- FALSE
  expect_error(call(bad), class = "gamlss.longitudinal_nonconvergence_error")
  bad <- f$joint; bad$coefficient_estimates[1] <- 500
  expect_error(call(f$object, bad), "do not match")
  term <- names(f$object$model_matrix$s$mu)[1]
  bad <- f$object; attr(bad$model_matrix$s$mu[[term]], "gl_smooth_spec") <- NULL
  expect_error(call(bad), "retained single-component")
  bad <- f$object
  attr(bad$model_matrix$s$mu[[term]], "gl_smooth_spec")$by <- "group"
  expect_error(call(bad), "without by variables")
  bad <- f$object
  attr(bad$model_matrix$s$mu[[term]], "gl_smooth_spec")$term <- "other_covariate"
  expect_error(call(bad), "univariate s1")
  comparator <- .gl_joint_unvalidated_gamlss2_contrast()
  expect_false(comparator$uncertainty_available)
  expect_false(comparator$clean_uncertainty_eligible)
  expect_true(is.na(comparator$std_error))
})
