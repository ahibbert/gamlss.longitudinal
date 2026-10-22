# Permanent fit-free regression candidate. Requires the unchanged independent
# curvature_model_fixture and joint_nbi_reference function declarations.
cg_nbi_chain_fixture <- function() {
  object <- curvature_model_fixture("C")
  n_subject <- length(unique(object$response_subject))
  angle <- (seq_len(n_subject) - 0.5) * 2 * pi / n_subject
  object$dataset$z <- rep(sqrt(2) * cos(2 * angle), each = 3L)
  object$dataset$response <- rep(c(0, 1, 3, 2, 8, 5, 1, 4, 2, 6, 3, 0), length.out = nrow(object$dataset))
  object$response <- object$dataset$response
  object$dataset_original <- object$dataset
  object$margin_dist <- gamlss.dist::NBI()
  object$par <- c(mu.intercept = log(3), mu.x = 0.15, mu.z = -0.12,
                  sigma.intercept = log(0.5), theta.intercept = log(1.2))
  object$model_matrix$x$mu <- cbind(intercept = 1, x = object$dataset$x, z = object$dataset$z)
  object
}

test_that("named same-margin chain terms preserve all coefficient blocks", {
  x <- c(-1, 0.5, 2); z <- c(0.4, -0.7, 1.2)
  X <- cbind(intercept = 1, x = x, z = z)
  mu <- as.numeric(exp(X %*% c(0.2, 0.3, -0.1)))
  scalar <- function(a) stats::setNames("mu", paste0("mu.", a))
  args <- list(dra = mu, drb = mu, margin_d1l = list(mu = rep(-1, 3)),
    margin_d2l = list(mu = list(mu = numeric(3))),
    copula_hess = list(cross_pair_contribs = list()), eta_dr = list(mu = mu),
    eta_d2 = list(mu = mu), pair_ok = logical(), id1_ok = integer(), id2_ok = integer())
  actual <- matrix(0, 3, 3, dimnames = list(colnames(X), colnames(X)))
  for (a in colnames(X)) for (b in colnames(X)) actual[a, b] <- do.call(
    .hessian_assembly_margin_margin_block,
    c(list(pa = scalar(a), pb = scalar(b), xa = X[, a], xb = X[, b]), args))
  # Poisson y=0: l(beta)=-sum(exp(X beta)), with independently known Hessian.
  expected <- -crossprod(X, X * mu)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_true(abs(expected["x", "z"]) > 0.1)
  expect_equal(actual, t(actual), tolerance = 0)
  args$dra <- args$drb <- rep(1, 3)
  args$eta_dr$mu <- rep(1, 3); args$eta_d2$mu <- numeric(3)
  args$margin_d2l$mu$mu <- rep(-1, 3)
  identity <- do.call(.hessian_assembly_margin_margin_block,
    c(list(pa = scalar("x"), pb = scalar("z"), xa = x, xb = z), args))
  expect_equal(identity, -sum(x * z), tolerance = 1e-12)
  args$margin_d2l$mu$sigma <- c(0.5, -0.2, 0.4)
  args$margin_d2l$sigma <- list(mu = args$margin_d2l$mu$sigma, sigma = numeric(3))
  args$eta_dr$sigma <- rep(2, 3); args$eta_d2$mu <- mu; args$eta_d2$sigma <- rep(2, 3)
  args$drb <- rep(2, 3)
  mixed <- do.call(.hessian_assembly_margin_margin_block,
    c(list(pa = scalar("x"), pb = c("sigma.z" = "sigma"), xa = x, xb = z), args))
  expect_equal(mixed, sum(x * z * c(0.5, -0.2, 0.4) * 2), tolerance = 1e-12)
  expect_equal(.gl_limit_cg_step(c(3, 4), trust_radius = 2, max_delta = 9), c(1.2, 1.6), tolerance = 1e-12)
  expect_equal(.gl_limit_cg_step(c(3, 4), trust_radius = 9, max_delta = 1), c(0.75, 1), tolerance = 1e-12)
  expect_equal(.gl_limit_cg_step(c(3, 4), trust_radius = 2, max_delta = 1), c(0.75, 1), tolerance = 1e-12)
})

test_that("CG nonlinear-margin Hessian and one bounded search match the defining likelihood", {
  stopifnot(requireNamespace("numDeriv", quietly = TRUE), requireNamespace("VineCopula", quietly = TRUE))
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  object <- cg_nbi_chain_fixture(); untouched <- object
  B <- matrix(object$dataset$z, ncol = 1L, dimnames = list(NULL, "mu.z"))
  objective <- function(beta) joint_nbi_reference(beta, object, B)
  beta0 <- object$par
  reference_start <- objective(beta0)
  zero <- stats::setNames(numeric(length(beta0)), names(beta0))
  reference_hessian <- function(eps) {
    calls <- 0L
    offset <- function(delta) {
      calls <<- calls + 1L
      if (calls > 250L || !is.numeric(delta) || length(delta) != 5L ||
          any(!is.finite(delta)) || any(abs(delta) > eps + 1e-15))
        stop("Independent CG reference exceeded its call or perturbation budget")
      objective(beta0 + delta)
    }
    numDeriv::hessian(offset, zero, method = "Richardson",
      method.args = list(eps = eps, d = 1e-3, r = 4, v = 2))
  }
  Href <- reference_hessian(0.01)
  Hcheck <- reference_hessian(0.02)
  dimnames(Href) <- dimnames(Hcheck) <- list(names(beta0), names(beta0))
  relative <- function(a, b) sqrt(sum((a - b)^2)) / max(1, sqrt(sum(b^2)))
  expect_lt(relative(Href, Hcheck), 2e-4)
  h <- 1e-5 * pmax(1, abs(beta0))
  gref <- vapply(seq_along(beta0), function(i) {
    plus <- beta0; plus[i] <- plus[i] + h[i]
    (objective(plus) - objective(beta0)) / h[i]
  }, numeric(1))
  newton_reference <- -as.numeric(solve(Href, gref))
  independent_limit <- function(delta) {
    norm <- sqrt(sum(delta^2))
    if (norm > 0.5) delta <- delta * 0.5 / norm
    if (max(abs(delta)) > 0.5) delta <- delta * 0.5 / max(abs(delta))
    delta
  }
  penalty <- matrix(0, length(beta0), length(beta0), dimnames = dimnames(Href))
  hessian_count <- 0L
  for (backend in c("native", "vinecopula")) {
    options(gamlss.longitudinal.copula_backend = backend)
    runtime <- .gl_build_cg_runtime_helpers(par_cov = object$par, par_s = object$par_s,
      mm = object$model_matrix, margin_dist = object$margin_dist,
      copula_link = get_copula_dist("C")$copula_link, copula_dist = "C", dataset = object$dataset,
      pair_cache = build_copula_pair_cache(object$response, object$response_margin, object$response_subject),
      margin_eval_cache = .build_margin_eval_cache(object$margin_dist, calc_d2 = FALSE),
      cg_gradient_method = "forward", cg_hessian_method = "analytical", verbose = 0)
    model <- runtime$build_model(object$model_matrix, object$par, object$par_s)
    expect_identical(model$beta, beta0)
    eval_start <- runtime$evaluate(model$beta, model$mm)
    expect_equal(eval_start$loglik, objective(beta0), tolerance = 1e-9)
    gradient <- runtime$gradient(model$beta, eval_start$loglik, model$mm)
    expect_equal(as.numeric(gradient), gref, tolerance = 2e-6)
    tmp <- .gl_build_cg_hessian_object(object$dataset, object$margin_dist, "C", model$mm, model$beta)
    # Same observed-Hessian dispatch used by the runtime, with finite fallback
    # deliberately refused so this cannot pass via another curvature algorithm.
    actual <- .gl_cg_observed_hessian(tmp, model$beta, model$mm,
      analytical_fn = function(obj) {
        hessian_count <<- hessian_count + 1L
        if (hessian_count > 2L) stop("CG regression Hessian budget exceeded")
        calc_analytical_hessian(obj, progress = FALSE)
      }, finite_fn = function(...) stop("Analytical CG regression must not fall back"),
      hessian_method = "analytical", verbose = 0)
    error <- sqrt(sum((actual - Href)^2)) / sqrt(sum(Href^2))
    expect_identical(dimnames(actual), list(names(beta0), names(beta0)))
    expect_lt(error, 2e-3)
    candidates <- .gl_build_cg_candidate_steps(gradient, actual, trust_radius = 0.5)
    expect_true(length(candidates) >= 2L)
    newton <- candidates[[2L]]
    expect_lt(relative(newton, newton_reference), 2e-3)
    expect_equal(.gl_limit_cg_step(newton, 0.5, 0.5), independent_limit(newton), tolerance = 1e-12)
    # One search episode per backend with two independently specified candidate
    # directions: Newton first, ascent fallback. No optimizer iteration loop.
    actual_steps <- list(newton, candidates[[1L]])
    reference_steps <- list(newton_reference, as.numeric(0.5 * gref / sqrt(sum(gref^2))))
    search <- function(steps, evaluator, gradient_value) {
      calls <- 0L
      counted <- function(beta, mm) {
        calls <<- calls + 1L
        if (calls > 60L) stop("CG regression line search budget exceeded")
        evaluator(beta, mm)
      }
      out <- .gl_run_cg_line_search(candidate_steps = steps, beta_vec = beta0, mm_cg = model$mm,
        penalty_current = penalty, obj_start = objective(beta0), trust_radius = 0.5,
        max_delta = 0.5, armijo_c1 = 1e-4, line_search = "best", max_line_search_evals = 60L,
        use_backtracking = TRUE, backtracking_max_halves = 8L, eval_fn = counted,
        objective_fn = function(beta, loglik, penalty) loglik,
        gradient = gradient_value, verbose = 0)
      expect_identical(out$line_eval_count, calls)
      out
    }
    chosen <- search(actual_steps, runtime$evaluate, gradient)
    # Independent oracle: no package candidate evaluator, clipping, objective,
    # Armijo or best-selection function is used here. Check selection on the
    # EXACT supplied actual directions, then on independent Newton directions.
    oracle <- function(steps, gradient_value) {
      best <- NULL; calls <- 0L
      for (direction in steps) for (halves in 0:8) {
        calls <- calls + 1L
        if (calls > 60L) stop("Independent CG search budget exceeded")
        delta <- independent_limit(direction / 2^halves)
        beta <- beta0 + delta
        ll <- objective(beta)
        improvement <- ll - reference_start
        slope <- sum(gradient_value * delta)
        threshold <- if (is.finite(slope) && slope > 0) 1e-4 * slope else 0
        if (is.finite(improvement) && improvement > max(1e-8, threshold) &&
            (is.null(best) || improvement > best$improvement))
          best <- list(beta = beta, improvement = improvement, loglik = ll)
      }
      list(best = best, line_eval_count = calls)
    }
    exact_oracle <- oracle(actual_steps, gradient)
    expected <- oracle(reference_steps, gref)
    expect_identical(chosen$line_eval_count, exact_oracle$line_eval_count)
    expect_false(is.null(exact_oracle$best))
    if (!is.null(chosen$best) && !is.null(exact_oracle$best)) {
      expect_equal(as.numeric(chosen$best$beta), as.numeric(exact_oracle$best$beta), tolerance = 1e-12)
      expect_equal(chosen$best$improvement, exact_oracle$best$improvement, tolerance = 1e-9)
    }
    expect_false(is.null(chosen$best)); expect_false(is.null(expected$best))
    if (!is.null(chosen$best) && !is.null(expected$best)) {
      expect_gt(chosen$best$improvement, 0)
      expect_lte(chosen$line_eval_count, 60L)
      expect_equal(chosen$best$eval$loglik, objective(chosen$best$beta), tolerance = 1e-9)
      expect_lt(relative(chosen$best$beta, expected$best$beta), 2e-3)
      expect_lt(abs(chosen$best$improvement - expected$best$improvement) /
        max(1, abs(expected$best$improvement)), 2e-3)
    }
    expect_identical(object, untouched)
    if (exists("cg_validation_output", inherits = TRUE)) {
      out <- get("cg_validation_output", inherits = TRUE)
      stopifnot(dir.exists(out))
      utils::write.csv(actual, file.path(out, paste0("cg-", backend, "-hessian.csv")))
      utils::write.csv(Href, file.path(out, paste0("cg-", backend, "-reference.csv")))
      utils::write.csv(data.frame(coordinate = names(beta0), forward_gradient = gradient,
        reference_gradient = gref, newton = newton, reference_newton = newton_reference),
        file.path(out, paste0("cg-", backend, "-gradient-newton.csv")), row.names = FALSE)
      saveRDS(list(backend = backend, hessian_relative_error = error,
        reference_step_relative_error = relative(Href, Hcheck), hessian_condition = kappa(Href),
        chosen = chosen, exact_oracle = exact_oracle, expected = expected, hessian_calls_so_far = hessian_count),
        file.path(out, paste0("cg-", backend, "-search.rds")), version = 3L)
    }
  }
  expect_identical(hessian_count, 2L)
})
