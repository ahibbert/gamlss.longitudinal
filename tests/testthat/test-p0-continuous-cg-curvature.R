curvature_cg_runtime <- function(object) {
  .gl_build_cg_runtime_helpers(par_cov = object$par, par_s = object$par_s,
    mm = object$model_matrix, margin_dist = object$margin_dist,
    copula_link = get_copula_dist(object$copula_dist)$copula_link,
    copula_dist = object$copula_dist, dataset = object$dataset,
    pair_cache = build_copula_pair_cache(object$response, object$response_margin, object$response_subject),
    margin_eval_cache = .build_margin_eval_cache(object$margin_dist, calc_d2 = FALSE),
    cg_gradient_method = "forward", cg_hessian_method = "analytical", verbose = 0)
}

test_that("CG curvature, candidate steps and line search agree with an independent likelihood", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N", "C", "F", "G", "J", "t")) {
    object <- curvature_model_fixture(family)
    object$par <- object$par + rep(c(0.05, -0.03, 0.02, 0.04, -0.02), length.out = length(object$par))
    original <- object
    ref_hessian <- curvature_reference_hessian(object)
    ref_gradient <- numDeriv::grad(function(beta)
      sum(curvature_reference_clusters(beta, object)), object$par)
    p <- length(object$par)
    penalty <- matrix(0, p, p, dimnames = list(names(object$par), names(object$par)))
    for (backend in c("native", "vinecopula")) {
      options(gamlss.longitudinal.copula_backend = backend)
      runtime <- curvature_cg_runtime(object)
      model <- runtime$build_model(object$model_matrix, object$par, object$par_s)
      eval <- runtime$evaluate(model$beta, model$mm)
      expect_equal(eval$loglik, sum(curvature_reference_clusters(model$beta, object)), tolerance = 1e-8)
      gradient <- runtime$gradient(model$beta, eval$loglik, model$mm)
      # Preserve CG's configured first-order forward algorithm. Its expected
      # truncation is h*H_ii/2; compare against that prediction and the same
      # forward stencil of the independent defining likelihood.
      h <- 1e-5 * pmax(1, abs(model$beta))
      forward_reference <- vapply(seq_along(model$beta), function(i) {
        plus <- model$beta; plus[i] <- plus[i] + h[i]
        (sum(curvature_reference_clusters(plus, object)) -
          sum(curvature_reference_clusters(model$beta, object))) / h[i]
      }, numeric(1))
      expect_equal(as.numeric(gradient), forward_reference, tolerance = 2e-6)
      predicted <- ref_gradient + h * diag(ref_hessian) / 2
      expect_lt(sqrt(sum((gradient - predicted)^2)) / max(1, sqrt(sum(ref_gradient^2))), 1e-6)
      tmp <- .gl_build_cg_hessian_object(object$dataset, object$margin_dist,
        family, model$mm, model$beta)
      actual <- runtime$observed_hessian(tmp, model$beta, model$mm)
      expect_lt(sqrt(sum((actual - ref_hessian)^2)) / sqrt(sum(ref_hessian^2)), 2e-4)
      candidates <- .gl_build_cg_candidate_steps(ref_gradient, actual, trust_radius = 0.5)
      reference_candidates <- .gl_build_cg_candidate_steps(ref_gradient, ref_hessian, trust_radius = 0.5)
      # Independently check the first Newton candidate before the algorithm's
      # symmetric alternatives and ridge ladder are added.
      newton <- -as.numeric(solve(ref_hessian, ref_gradient))
      expect_equal(.gl_limit_cg_step(candidates[[2L]], 0.5, 0.5),
        .gl_limit_cg_step(newton, 0.5, 0.5), tolerance = 2e-4)
      search <- function(steps, evaluator) .gl_run_cg_line_search(
        candidate_steps = steps, beta_vec = model$beta, mm_cg = model$mm,
        penalty_current = penalty, obj_start = eval$loglik,
        trust_radius = 0.5, max_delta = 0.5, armijo_c1 = 1e-4,
        line_search = "best", max_line_search_evals = 60L,
        use_backtracking = TRUE, backtracking_max_halves = 50L,
        eval_fn = evaluator, objective_fn = runtime$objective,
        gradient = ref_gradient, verbose = 0)
      chosen <- search(candidates, runtime$evaluate)
      expected <- search(reference_candidates, function(beta, mm)
        list(loglik = sum(curvature_reference_clusters(beta, object))))
      expect_false(is.null(chosen$best))
      expect_false(is.null(expected$best))
      expect_gt(chosen$best$improvement, 0)
      expect_lte(chosen$line_eval_count, 60L)
      expect_equal(as.numeric(chosen$best$beta), as.numeric(expected$best$beta), tolerance = 2e-4)
      expect_equal(chosen$best$improvement, expected$best$improvement, tolerance = 2e-4)
      expect_identical(object, original)
      if (exists("validation_details", inherits = TRUE)) {
        validation_details$cg[[paste(backend, family)]] <- list(family = family,
          backend = backend, hessian_relative_error = sqrt(sum((actual - ref_hessian)^2)) / sqrt(sum(ref_hessian^2)),
          forward_gradient_relative_error = sqrt(sum((gradient - ref_gradient)^2)) / max(1, sqrt(sum(ref_gradient^2))),
          chosen_beta = chosen$best$beta, reference_beta = expected$best$beta,
          improvement = chosen$best$improvement, reference_improvement = expected$best$improvement,
          evaluations = chosen$line_eval_count)
      }
    }
  }
})

test_that("CG smooth EDF uses independently verified augmented likelihood curvature", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  object <- curvature_model_fixture()
  angle <- (seq_len(32L) - 0.5) * 2 * pi / 32L
  subject_basis <- cbind(mu.s1.1 = sqrt(2) * cos(5 * angle), mu.s1.2 = sqrt(2) * cos(6 * angle))
  B <- subject_basis[rep(seq_len(32L), each = 3L), , drop = FALSE]
  attr(B, "penalty") <- diag(2)
  object$model_matrix$s$mu <- list(s1 = B)
  object$par_s$mu <- list(s1 = c(mu.s1.1 = 0.05, mu.s1.2 = -0.03))
  runtime <- curvature_cg_runtime(object)
  model <- runtime$build_model(object$model_matrix, object$par, object$par_s)
  reference_object <- object
  reference_object$par <- model$beta
  reference_object$synthetic_smooth_covariates <- subject_basis
  ref <- curvature_reference_hessian(reference_object)
  dimnames(ref) <- list(names(model$beta), names(model$beta))
  tmp <- .gl_build_cg_hessian_object(object$dataset, object$margin_dist,
    "N", model$mm, model$beta)
  H <- runtime$observed_hessian(tmp, model$beta, model$mm)
  expect_lt(sqrt(sum((H - ref)^2)) / sqrt(sum(ref^2)), 2e-4)
  lambda <- object$lambda_s; lambda$mu <- list(s1 = 2.5)
  P <- runtime$build_penalty(names(model$beta), lambda)
  idx <- names(object$par_s$mu$s1)
  information <- -ref[idx, idx]
  ridge <- max(1e-8, 1e-8 * max(1, max(abs(diag(information)))))
  expected <- sum(diag(solve(information + P[idx, idx] + diag(ridge, 2), information)))
  actual <- runtime$smooth_edf_list(H, P, names(model$beta))$mu$s1
  expect_equal(actual, expected, tolerance = 2e-5)
  expect_true(actual > 0 && actual < 2)
})
test_that("CG reloads a partial standalone Hessian dependency set", {
  ready <- FALSE
  sourced <- character()
  expect_invisible(.gl_ensure_cg_hessian_available(package_root = "pkgroot",
    exists_fn = function(name, mode) ready || name == "calc_analytical_hessian",
    source_fn = function(file, local) {
      sourced <<- c(sourced, basename(file))
      if (basename(file) == "hessian-analytical.R") ready <<- TRUE
    }, file_exists_fn = function(paths) rep(TRUE, length(paths))))
  expect_true(all(c("copula-backend-parameter-stencil.R",
    "hessian-copula-pair-eval.R", "hessian-assembly-margin-copula-blocks.R") %in% sourced))
  expect_error(.gl_ensure_cg_hessian_available(
    exists_fn = function(name, mode) name == "calc_analytical_hessian",
    source_fn = function(...) stop("Must not source an incomplete file set"),
    file_exists_fn = function(paths) basename(paths) != "copula-backend-parameter-stencil.R"),
    "complete analytical Hessian helpers")
})

test_that("standalone CG source loading produces the independently checked Hessian", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  loaded <- new.env(parent = environment(.gl_ensure_cg_hessian_available))
  loaded$calc_analytical_hessian <- function(...) stop("Stale partial Hessian")
  .gl_ensure_cg_hessian_available(package_root = test_path("..", ".."),
    exists_fn = function(name, mode) exists(name, envir = loaded, mode = mode, inherits = FALSE),
    source_fn = function(file, local) sys.source(file, envir = loaded))
  object <- curvature_model_fixture("t")
  # This separately sourced function is outside the namespace tracer. Count
  # its one evaluation explicitly when running under the bounded test driver.
  count <- getOption("curvature.validation.counters")
  if (!is.null(count)) {
    count$hessian_calls <- count$hessian_calls + 1L
    if (count$hessian_calls > 80L) stop("Synthetic Hessian budget exceeded")
  }
  actual <- loaded$calc_analytical_hessian(object, progress = FALSE)
  reference <- curvature_reference_hessian(object)
  expect_lt(sqrt(sum((actual-reference)^2)) / sqrt(sum(reference^2)), 1e-4)
  expect_true(exists(".copula_t_df_derivative_domain", envir = loaded, inherits = FALSE))
})
