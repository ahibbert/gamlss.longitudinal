#' Decide whether the one-shot RS BFGS refinement is due
#'
#' @noRd
.gl_should_run_rs_bfgs_refinement <- function(
    completed_outer,
    enabled,
    refine_after,
    already_refined,
    refinements_completed = as.integer(isTRUE(already_refined)),
    max_refinements = 1L,
    lambda_update_scheduled = FALSE,
    objective_change = Inf,
    outer_tol = 0,
    prior_refinement_incomplete = FALSE) {
  if (!isTRUE(enabled) || isTRUE(already_refined) ||
      refinements_completed >= max_refinements || completed_outer < refine_after) {
    return(FALSE)
  }
  if (refinements_completed == 0L) return(TRUE)
  isTRUE(lambda_update_scheduled) &&
    (isTRUE(prior_refinement_incomplete) ||
      !is.finite(objective_change) || abs(objective_change) > outer_tol)
}

#' Rescue weak large-magnitude fixed-effect coordinates after capped BFGS
#'
#' @noRd
.gl_run_rs_coordinate_rescue <- function(
    beta,
    fixed_names,
    value_fn,
    threshold = 1.5,
    radius_floor = 4,
    max_coordinates = 10L,
    tolerance = 1e-3,
    check_elapsed_budget = function(stage) NULL,
    optimize_fn = stats::optimize) {
  fixed_names <- intersect(fixed_names, names(beta))
  selected <- fixed_names[abs(beta[fixed_names]) >= threshold]
  if (length(selected)) {
    selected <- selected[order(abs(beta[selected]), decreasing = TRUE)]
    selected <- utils::head(selected, as.integer(max_coordinates))
  }
  current_value <- value_fn(beta)
  initial_value <- current_value
  trace <- vector("list", length(selected))
  for (index in seq_along(selected)) {
    check_elapsed_budget("RS BFGS coordinate rescue")
    coefficient <- selected[[index]]
    start <- beta[[coefficient]]
    start_value <- current_value
    radius <- max(radius_floor, abs(start))
    result <- tryCatch(
      optimize_fn(
        function(value) {
          candidate <- beta
          candidate[[coefficient]] <- value
          value_fn(candidate)
        },
        interval = start + c(-radius, radius),
        maximum = TRUE,
        tol = tolerance
      ),
      error = function(error) {
        if (inherits(error, "gamlss.longitudinal_time_limit_error")) stop(error)
        error
      }
    )
    accepted <- !inherits(result, "error") && is.finite(result$objective) &&
      result$objective >= current_value
    if (accepted) {
      beta[[coefficient]] <- result$maximum
      current_value <- result$objective
    }
    trace[[index]] <- data.frame(
      coefficient = coefficient,
      start = start,
      end = beta[[coefficient]],
      start_value = start_value,
      end_value = current_value,
      improvement = current_value - start_value,
      accepted = accepted,
      stringsAsFactors = FALSE
    )
  }
  list(
    beta = beta,
    value = current_value,
    trace = if (length(trace)) do.call(rbind, trace) else data.frame(),
    used = length(selected) > 0L,
    improved = is.finite(current_value) && is.finite(initial_value) &&
      current_value > initial_value
  )
}

#' Refine an RS state jointly with a penalized analytical-gradient BFGS step
#'
#' The returned state is always followed by ordinary RS cycles. BFGS therefore
#' accelerates correlated blocks but cannot itself satisfy the RS convergence
#' contract or bypass scheduled smoothing-parameter updates.
#'
#' @noRd
.gl_run_rs_bfgs_refinement <- function(
    mm,
    par_cov,
    par_s,
    lambda_s,
    margin_dist,
    copula_link,
    copula_dist,
    dataset,
    pair_cache,
    margin_eval_cache,
    include_dlcopdpar,
    method = "BFGS",
    max_iter,
    reltol,
    coordinate_rescue = FALSE,
    coordinate_threshold = 1.5,
    coordinate_radius = 4,
    coordinate_max = 10L,
    coordinate_tol = 1e-3,
    repolish_iter = max_iter,
    check_elapsed_budget,
    verbose = 0L,
    build_model_fn = .gl_build_cg_model,
    build_penalty_fn = .gl_build_cg_penalty,
    evaluate_fn = .gl_evaluate_cg_beta,
    gradient_fn = .cg_analytical_gradient,
    optim_fn = stats::optim) {
  started <- proc.time()[["elapsed"]]
  augmented <- build_model_fn(mm, par_cov, par_s)
  beta_start <- augmented$beta
  penalty <- build_penalty_fn(
    beta_names = names(beta_start),
    lambda_current = lambda_s,
    par_s = par_s,
    mm = mm
  )

  # stats::optim() commonly requests fn and gr at exactly the same coefficient
  # vector.  A joint discrete likelihood evaluation is expensive and already
  # contains every quantity needed by both closures, so retain the most recent
  # state instead of calculating it twice.  Exact coefficient equality is
  # intentional: this changes no numerical approximation or optimizer
  # semantics, and invalid evaluations are cached as well.
  last_beta <- NULL
  last_evaluation <- NULL
  has_last_evaluation <- FALSE
  evaluation_calls <- 0L
  evaluation_cache_hits <- 0L
  evaluation_cache_misses <- 0L

  evaluate <- function(beta) {
    check_elapsed_budget("RS BFGS refinement")
    evaluation_calls <<- evaluation_calls + 1L
    beta_numeric <- as.numeric(beta)
    cache_hit <- isTRUE(has_last_evaluation) &&
      length(beta_numeric) == length(last_beta) &&
      isTRUE(all(beta_numeric == last_beta))
    if (cache_hit) {
      evaluation_cache_hits <<- evaluation_cache_hits + 1L
      return(last_evaluation)
    }
    evaluation_cache_misses <<- evaluation_cache_misses + 1L
    evaluation <- suppressWarnings(evaluate_fn(
      beta_vec = beta,
      mm_cg = augmented$mm,
      par_cov_template = par_cov,
      par_s_template = par_s,
      margin_dist = margin_dist,
      copula_link = copula_link,
      copula_dist = copula_dist,
      dataset = dataset,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache
    ))
    last_beta <<- beta_numeric
    last_evaluation <<- evaluation
    has_last_evaluation <<- TRUE
    evaluation
  }
  penalized_value <- function(beta, evaluation) {
    as.numeric(
      evaluation$loglik - 0.5 * sum(beta * as.numeric(penalty %*% beta))
    )
  }

  initial_evaluation <- evaluate(beta_start)
  if (is.null(initial_evaluation)) {
    return(list(
      par_cov = par_cov,
      par_s = par_s,
      diagnostics = list(
        enabled = TRUE,
        used = TRUE,
        accepted = FALSE,
        reason = "invalid_start",
        evaluation_calls = evaluation_calls,
        evaluation_cache_hits = evaluation_cache_hits,
        evaluation_cache_misses = evaluation_cache_misses,
        elapsed_sec = proc.time()[["elapsed"]] - started
      )
    ))
  }
  initial_value <- penalized_value(beta_start, initial_evaluation)

  objective <- function(beta) {
    evaluation <- evaluate(beta)
    if (is.null(evaluation)) return(1e100)
    -penalized_value(beta, evaluation)
  }
  gradient <- function(beta) {
    evaluation <- evaluate(beta)
    if (is.null(evaluation)) {
      stop("BFGS requested a gradient at an invalid likelihood state.", call. = FALSE)
    }
    score <- gradient_fn(
      beta,
      augmented$mm,
      evaluation$eta_out,
      evaluation$calc_lik,
      margin_dist,
      copula_dist,
      isTRUE(include_dlcopdpar),
      dataset$response,
      dataset$time,
      dataset$subject
    )
    -stats::setNames(
      as.numeric(score - penalty %*% beta),
      names(beta)
    )
  }

  run_optim <- function(start, iterations) {
      # Dense BFGS otherwise starts with an identity inverse Hessian on the raw
      # coefficient scale. Sparse interaction coefficients can legitimately be
      # several times larger than dense main effects, making progress along
      # those weak directions extremely slow. Optimizing beta / parscale keeps
      # the initial metric commensurate without changing the objective.
      parameter_scale <- pmin(pmax(1, abs(as.numeric(start))), 25)
      optim_control <- if (identical(method, "L-BFGS-B")) {
        list(
          maxit = as.integer(iterations),
          factr = max(1, reltol / .Machine$double.eps),
          pgtol = 0,
          lmm = 20L,
          parscale = parameter_scale,
          trace = if (verbose > 2L) 1L else 0L,
          REPORT = 10L
        )
      } else {
        list(
          maxit = as.integer(iterations),
          reltol = reltol,
          parscale = parameter_scale,
          trace = if (verbose > 2L) 1L else 0L,
          REPORT = 10L
        )
      }
      result <- optim_fn(
      par = start,
      fn = objective,
      gr = gradient,
      method = method,
      control = optim_control
    )
      attr(result, "parameter_scale") <- parameter_scale
      result
  }
  optimized <- tryCatch(
    run_optim(beta_start, max_iter),
    error = function(error) {
      if (inherits(error, "gamlss.longitudinal_time_limit_error")) {
        stop(error)
      }
      error
    }
  )
  if (inherits(optimized, "error")) {
    return(list(
      par_cov = par_cov,
      par_s = par_s,
      diagnostics = list(
        enabled = TRUE,
        used = TRUE,
        accepted = FALSE,
        reason = "optim_error",
        error = conditionMessage(optimized),
        evaluation_calls = evaluation_calls,
        evaluation_cache_hits = evaluation_cache_hits,
        evaluation_cache_misses = evaluation_cache_misses,
        elapsed_sec = proc.time()[["elapsed"]] - started
      )
    ))
  }

  final_evaluation <- evaluate(optimized$par)
  final_value <- if (is.null(final_evaluation)) -Inf else {
    penalized_value(optimized$par, final_evaluation)
  }
  rescue <- list(used = FALSE, improved = FALSE, trace = data.frame())
  repolish <- NULL
  if (isTRUE(coordinate_rescue) &&
      (optimized$convergence %||% 0L) != 0L &&
      is.finite(final_value) && final_value >= initial_value) {
    rescue <- .gl_run_rs_coordinate_rescue(
      beta = optimized$par,
      fixed_names = names(par_cov),
      value_fn = function(beta) {
        evaluation <- evaluate(beta)
        if (is.null(evaluation)) return(-1e100)
        penalized_value(beta, evaluation)
      },
      threshold = coordinate_threshold,
      radius_floor = coordinate_radius,
      max_coordinates = coordinate_max,
      tolerance = coordinate_tol,
      check_elapsed_budget = check_elapsed_budget
    )
    if (isTRUE(rescue$improved)) {
      rescue_evaluation <- evaluate(rescue$beta)
      rescue_value <- if (is.null(rescue_evaluation)) -Inf else {
        penalized_value(rescue$beta, rescue_evaluation)
      }
      if (is.finite(rescue_value) && rescue_value >= final_value) {
        optimized$par <- rescue$beta
        final_evaluation <- rescue_evaluation
        final_value <- rescue_value
      }
      repolish <- tryCatch(
        run_optim(rescue$beta, repolish_iter),
        error = function(error) {
          if (inherits(error, "gamlss.longitudinal_time_limit_error")) stop(error)
          error
        }
      )
      if (!inherits(repolish, "error")) {
        repolish_evaluation <- evaluate(repolish$par)
        repolish_value <- if (is.null(repolish_evaluation)) -Inf else {
          penalized_value(repolish$par, repolish_evaluation)
        }
        if (is.finite(repolish_value) && repolish_value >= final_value) {
          optimized <- repolish
          final_evaluation <- repolish_evaluation
          final_value <- repolish_value
        }
      }
    }
  }
  accepted <- is.finite(final_value) && final_value >= initial_value
  if (accepted) {
    unpacked <- .gl_unpack_cg_beta(
      optimized$par,
      par_cov_template = par_cov,
      par_s_template = par_s
    )
    par_cov <- unpacked$par_cov
    par_s <- unpacked$par_s
  }

  list(
    par_cov = par_cov,
    par_s = par_s,
    diagnostics = list(
      enabled = TRUE,
      used = TRUE,
      accepted = accepted,
      reason = if (accepted) "accepted" else "no_penalized_improvement",
      method = method,
      initial_penalized_logLik = initial_value,
      final_penalized_logLik = final_value,
      improvement = final_value - initial_value,
      convergence = optimized$convergence,
      message = optimized$message %||% NULL,
      counts = optimized$counts,
      evaluation_calls = evaluation_calls,
      evaluation_cache_hits = evaluation_cache_hits,
      evaluation_cache_misses = evaluation_cache_misses,
      parscale_min = min(attr(optimized, "parameter_scale")),
      parscale_max = max(attr(optimized, "parameter_scale")),
      coordinate_rescue = rescue,
      repolish = if (inherits(repolish, "error")) {
        list(used = TRUE, error = conditionMessage(repolish))
      } else if (is.null(repolish)) {
        list(used = FALSE)
      } else {
        list(
          used = TRUE,
          convergence = repolish$convergence,
          counts = repolish$counts
        )
      },
      elapsed_sec = proc.time()[["elapsed"]] - started
    )
  )
}
