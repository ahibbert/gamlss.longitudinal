# Fit-free tests. Closed-form oracles are independent of the production kernel.
test_that("BCPE log-mu curvature retains exact signed and regularity oracles", {
  normal <- .hessian_bcpe_mu_curvature(c(.25, 1, 2, 4), 1, c(1, 2, 1, 2), 0, 2)
  expect_true(all(normal$available))
  expect_equal(normal$curvature, c(-1, -.25, -1, -.25), tolerance = 1e-14)
  signed <- .hessian_bcpe_mu_curvature(c(.25, 2), 1, 1, 1, 2)
  expect_true(all(signed$available))
  expect_equal(signed$curvature, c(.125, -6), tolerance = 1e-14)
  zero <- .hessian_bcpe_mu_curvature(1, 1, 2, c(0, 1, -.5, 2), c(3, 1, 1.5, 2))
  expect_identical(zero$available, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(zero$curvature[c(1, 4)], c(0, -.25), tolerance = 1e-14)
  expect_identical(zero$reason[2:3], rep("no_finite_second_derivative_at_zero", 2))
  linear <- .hessian_bcpe_mu_curvature(c(.5, 2), 1, 1, 0, 1)
  expect_true(all(linear$available))
  expect_identical(linear$curvature, c(0, 0))
  # tau*r=1 here: do not mislabel unresolved floating cancellation as a cusp.
  cancellation <- .hessian_bcpe_mu_curvature(.5, 1, 1, 1, 2)
  expect_false(cancellation$available)
  expect_identical(cancellation$reason, "unresolved_log_bracket_cancellation")
})

test_that("BCPE eta override uses the observed subset and refuses unavailable rows", {
  family <- list(family = c("BCPE", "Box-Cox power exponential"), mu.link = "log")
  eta <- list(mu = c(1, NA_real_, 1), sigma = c(1, -1, 2),
    nu = c(0, Inf, 0), tau = c(2, NA_real_, 2))
  # Row 3 is observed even when it has no complete copula pair.
  expect_equal(.hessian_bcpe_log_mu_eta_curvature(eta, family, c(.5, NA_real_, 3)),
    c(-1, 0, -.25), tolerance = 1e-14)
  expect_identical(.hessian_bcpe_log_mu_eta_curvature(list(), family, c(NA_real_, NA_real_)), c(0, 0))
  expect_error(.hessian_bcpe_log_mu_eta_curvature(modifyList(eta, list(mu = 1:2)),
    family, c(1, NA_real_, 2)), "plain scalar or row vector")
  invalid <- tryCatch(.hessian_bcpe_log_mu_eta_curvature(eta, family, c(.5, 2, 3)), error = identity)
  expect_s3_class(invalid, "gamlss_longitudinal_score_evaluation_error")
  expect_identical(invalid$rows, 2L)
  expect_identical(invalid$reasons, "invalid_observed_BCPE_domain")
  eta <- list(mu = 1, sigma = 1, nu = 0, tau = c(2, 1.5, 2))
  unavailable <- tryCatch(.hessian_bcpe_log_mu_eta_curvature(eta, family, c(.5, 1, 2)), error = identity)
  expect_s3_class(unavailable, "gamlss_longitudinal_score_evaluation_error")
  expect_identical(unavailable$rows, 2L)
  expect_identical(unavailable$reasons, "no_finite_second_derivative_at_zero")
  family$mu.link <- "identity"
  expect_null(.hessian_bcpe_log_mu_eta_curvature(list(), family, -1))
  family$mu.link <- "log"; family$family <- "BCPEo"
  expect_null(.hessian_bcpe_log_mu_eta_curvature(list(), family, -1))
  family$family <- "NO"
  expect_null(.hessian_bcpe_log_mu_eta_curvature(list(), family, -1))
})

test_that("eta override replaces the marginal chain rule in every fixed and smooth mu block", {
  X <- cbind(intercept = 1, x = c(-1, .5, 2), smooth = c(.4, -.7, 1.2))
  colnames(X)[3] <- "s(s1, k = 10).1"
  qnames <- paste0("mu.", colnames(X))
  object <- list(par = setNames(numeric(3), qnames), weights = list(mu = c(100, 200, 300)))
  mu <- c(1, 2, 3); B <- c(.125, -6, -1)
  c1 <- c(.4, -.1, .2); c2 <- c(-.2, .3, .1); cross <- .7
  cop <- list(row_id1 = 1L, row_id2 = 2L, pair_ok = TRUE,
    cop_d1l_margin = list(mu = c1), cop_d2l_margin = list(mu = list(mu = c2)),
    cross_pair_contribs = list(mu = list(mu = cross)))
  args <- list(margin_d1l = list(mu = rep(NA_real_, 3)),
    margin_d2l = list(mu = list(mu = rep(NA_real_, 3))), copula_hess = cop,
    eta_dr = list(mu = mu), eta_d2 = list(mu = mu), mm = list(x = list(mu = X)),
    pair_cache = list(), margin_mu_eta_d2 = B)
  actual <- do.call(.assemble_covariate_hessian, c(list(object = object), args))
  expected <- crossprod(X, X * (B + c2 * mu^2 + c1 * mu)) +
    cross * mu[1] * mu[2] * (tcrossprod(X[1, ], X[2, ]) + tcrossprod(X[2, ], X[1, ]))
  dimnames(expected) <- list(qnames, qnames)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_identical(dimnames(actual), list(qnames, qnames))
  expect_equal(actual, t(actual), tolerance = 0)
  object$weights$mu <- c(.001, 99, 0)
  expect_identical(do.call(.assemble_covariate_hessian, c(list(object = object), args)), actual)
  # The observed unpaired third row contributes to the marginal block.
  args$margin_mu_eta_d2[3] <- 0
  without_third <- do.call(.assemble_covariate_hessian, c(list(object = object), args))
  expected_third <- B[3] * tcrossprod(X[3, ])
  dimnames(expected_third) <- list(qnames, qnames)
  expect_equal(actual - without_third, expected_third, tolerance = 1e-12)
  args$margin_mu_eta_d2 <- c(1, 2)
  expect_error(do.call(.assemble_covariate_hessian, c(list(object = object), args)), "equally sized")
  args$margin_mu_eta_d2 <- c(1, NA_real_, 2)
  args$mm$x$mu[2, ] <- 0
  expect_error(do.call(.assemble_covariate_hessian, c(list(object = object), args)), "finite")
  overflow <- list(pa = "mu", pb = "mu", xa = rep(1e200, 2), xb = rep(1e200, 2),
    dra = c(1, 1), drb = c(1, 1), margin_d1l = list(), margin_d2l = list(),
    copula_hess = list(cross_pair_contribs = list()), eta_dr = list(mu = c(1, 1)),
    eta_d2 = list(mu = c(1, 1)), pair_ok = logical(), id1_ok = integer(),
    id2_ok = integer(), margin_mu_eta_d2 = c(0, 0))
  # Finite inputs can overflow before multiplication by zero: never drop NaN.
  expect_error(do.call(.hessian_assembly_margin_margin_block, overflow), "row products.*representable")
  overflow$xa <- c(1e308, 1e308); overflow$xb <- c(1, 1)
  overflow$margin_mu_eta_d2 <- c(1, 1)
  # Each contribution is finite, but their aggregate is not representable.
  expect_error(do.call(.hessian_assembly_margin_margin_block, overflow), "row sum.*representable")
})

test_that("mixed and default natural assembly paths retain their chain rule", {
  args <- list(pa = c(mu.x = "mu"), pb = c(sigma.x = "sigma"),
    xa = c(1, 2), xb = c(3, 4), dra = c(2, 3), drb = c(4, 5),
    margin_d1l = list(mu = c(9, 8)),
    margin_d2l = list(mu = list(mu = c(-2, -3), sigma = c(.5, -.2)),
      sigma = list(mu = c(.5, -.2), sigma = c(-1, -1))),
    copula_hess = list(cross_pair_contribs = list(),
      cop_d1l_margin = list(mu = c(1, 2)),
      cop_d2l_margin = list(mu = list(mu = c(.1, .2), sigma = c(.3, .4)))),
    eta_dr = list(mu = c(2, 3), sigma = c(4, 5)), eta_d2 = list(mu = c(2, 3)),
    pair_ok = logical(), id1_ok = integer(), id2_ok = integer(), margin_mu_eta_d2 = c(-100, -200))
  expect_equal(do.call(.hessian_assembly_margin_margin_block, args),
    sum(c(1, 2) * c(3, 4) * (c(.5, -.2) + c(.3, .4)) * c(2, 3) * c(4, 5)), tolerance = 1e-12)
  args$pb <- c(mu.z = "mu"); args$drb <- args$dra; args$margin_mu_eta_d2 <- NULL
  expected <- sum(args$xa * args$xb * ((c(-2, -3) + c(.1, .2)) * c(2, 3)^2 +
    (c(9, 8) + c(1, 2)) * c(2, 3)))
  expect_equal(do.call(.hessian_assembly_margin_margin_block, args), expected, tolerance = 1e-12)
})

test_that("mu diagonal bypass omits only its two density evaluations", {
  fd <- .calc_margin_d2l_fd
  fake <- new.env(parent = environment(fd)); contexts <- character()
  fake$.hessian_margin_density_setup <- function(eta_inv, mm, margin_dist, response)
    list(margin_pars = c("mu", "sigma"), dfun = NULL, base_args = eta_inv)
  fake$.hessian_margin_log_density_eval <- function(base_args, dfun, margin_dist, response, context) {
    contexts <<- c(contexts, context)
    # Known natural Hessian: mu-mu=-2, sigma-sigma=-6, mu-sigma=2.
    with(base_args, -mu^2 - 3 * sigma^2 + 2 * mu * sigma)
  }
  environment(fd) <- fake
  eta <- list(mu = c(1, 2), sigma = c(2, 3))
  family <- list(mu.link = "identity", sigma.link = "identity")
  full <- fd(eta, list(), family, c(1, 2), h = .25)
  all_contexts <- contexts; contexts <- character()
  reduced <- fd(eta, list(), family, c(1, 2), h = .25, skip_mu_diagonal = TRUE)
  expect_length(all_contexts, 9L); expect_length(contexts, 7L)
  expect_identical(setdiff(all_contexts, contexts), c(
    "Positive marginal log-density perturbation for 'mu'",
    "Negative marginal log-density perturbation for 'mu'"))
  expect_null(reduced$mu$mu)
  expect_named(reduced$mu, c("mu", "sigma"))
  expect_equal(full$mu$mu, c(-2, -2), tolerance = 0)
  expect_equal(reduced$sigma$sigma, c(-6, -6), tolerance = 0)
  expect_equal(reduced$mu$sigma, c(2, 2), tolerance = 0)
  expect_identical(reduced$mu$sigma, full$mu$sigma)
  expect_identical(reduced$sigma$mu, full$sigma$mu)
  expect_error(fd(eta, list(), family, c(1, 2), skip_mu_diagonal = NA), "nonmissing logical")
})

test_that("continuous Hessian entry threads the override and diagonal bypass together", {
  run <- calc_analytical_hessian
  fake <- new.env(parent = environment(run)); seen <- NULL
  fake$get_copula_dist <- function(...) list(copula_link = "identity")
  fake$calc_eta <- function(...) list(eta = list(mu = c(0, 0)),
    eta_inv = list(mu = c(1, 1), sigma = 1, nu = 0, tau = 2), eta_dr = list(mu = c(1, 1)))
  fake$.calc_eta_d2_linkinv <- function(...) list(mu = c(1, 1))
  fake$.warn_gg_near_zero_nu <- function(...) invisible(NULL)
  fake$build_copula_pair_cache <- function(...) list()
  fake$calc_likelihood_minimal <- function(...) list(likelihood_type = "continuous", margin_p = c(.2, .8), margin_deriv = list())
  fake$.calc_dFdpar <- fake$.calc_d2Fdpar2 <- function(...) list()
  fake$.calc_d2Fdpar_cross <- function(...) list()
  fake$.calc_copula_hessian_contributions <- function(...) list()
  fake$.get_margin_d1l <- function(...) list(mu = c(111, 222))
  fake$.calc_margin_d2l_fd <- function(..., skip_mu_diagonal = FALSE) {
    seen <<- skip_mu_diagonal
    list(mu = list(mu = if (skip_mu_diagonal) NULL else c(-2, -2)))
  }
  fake$.assemble_covariate_hessian <- function(..., margin_mu_eta_d2 = NULL)
    list(override = margin_mu_eta_d2, skip = seen)
  environment(run) <- fake
  object <- list(response = c(.5, 2), margin_dist = list(family = "BCPE", mu.link = "log"), copula_dist = "N")
  expect_identical(run(object, FALSE), list(override = c(-1, -1), skip = TRUE))
  # Smooth products in the real calc_eta path retain n-by-1 dimensions.
  # Exercise the analytical entry, not just the boundary adapter in isolation.
  columns <- list(eta = list(mu = matrix(0, 2, 1)),
    eta_inv = list(mu = matrix(1, 2, 1), sigma = matrix(1, 2, 1),
      nu = 0, tau = matrix(2, 1, 1)), eta_dr = list(mu = c(1, 1)))
  unchanged <- columns
  fake$calc_eta <- function(...) columns
  expect_identical(run(object, FALSE), list(override = c(-1, -1), skip = TRUE))
  expect_identical(columns, unchanged)
  object$margin_dist$mu.link <- "identity"
  expect_identical(run(object, FALSE), list(override = NULL, skip = FALSE))
})

test_that("only BCPE log analytical caches acquire the marginal curvature marker", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  family <- list(family = "BCPE", type = "Continuous", mu.link = "log")
  fit <- structure(list(copula_dist = "t", margin_dist = family,
    vcov = list(vcov = list(overall = diag(2)), method = "analytical",
      hessian_diagnostics = list(status = "available")),
    vcov_meta = list(cache_version = .gl_vcov_cache_version(), method = "analytical", numderiv = FALSE)),
    class = "gamlss.longitudinal")
  marker <- .gl_vcov_continuous_curvature_marker(fit)
  expect_identical(marker$bcpe_log_mu_marginal_curvature_version, 1L)
  expect_identical(marker$version, 2L)
  expect_identical(.gl_vcov_cache_version(), 3L)
  fit$vcov$continuous_curvature <- marker
  expect_true(.can_use_cached_vcov(fit, method = "analytical"))
  fit$vcov$continuous_curvature$bcpe_log_mu_marginal_curvature_version <- NULL
  expect_false(.can_use_cached_vcov(fit, method = "analytical"))
  fit$vcov$continuous_curvature$bcpe_log_mu_marginal_curvature_version <- 0L
  expect_false(.can_use_cached_vcov(fit, method = "analytical"))
  fit$margin_dist$mu.link <- "identity"
  fit$vcov$continuous_curvature <- .gl_vcov_continuous_curvature_marker(fit)
  expect_null(fit$vcov$continuous_curvature$bcpe_log_mu_marginal_curvature_version)
  expect_true(.can_use_cached_vcov(fit, method = "analytical"))
  fit$margin_dist$family <- "NO"; fit$margin_dist$mu.link <- "log"
  expect_null(.gl_vcov_continuous_curvature_marker(fit)$bcpe_log_mu_marginal_curvature_version)
  build <- .gl_vcov_build_result
  fake <- new.env(parent = environment(build))
  fake$.gl_compute_smooth_vcov <- function(...) list(smooth_vcov = list(), smooth_se = list())
  fake$.gl_enrich_vcov_contract <- function(out, object) out
  environment(build) <- fake
  fit$margin_dist <- family
  result <- function(method, diagnostics = list(status = "available"))
    build(fit, list(), numeric(), diag(2), c(1, 1), method, method, diagnostics)
  expect_identical(result("analytical")$continuous_curvature, marker)
  expect_null(result("numderiv")$continuous_curvature)
  expect_null(result("sandwich_cluster", list(bread_method = "numderiv"))$continuous_curvature)
  expect_identical(result("sandwich_cluster", list(bread_method = "analytical"))$continuous_curvature, marker)
})
