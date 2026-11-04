# Synthetic objective/control fixtures: no retained data or distribution fits.
test_that("RS accepts penalized ascent and rejects a single downhill try", {
  mk <- function(b) list(par_cov = numeric(), par_s = list(mu = list(s = c(b = b))),
    calc_lik_out_end = list(log_lik = c(joint = -(b - 2)^2 / 2), valid = TRUE))
  raw <- function(x) unname(x$calc_lik_out_end$log_lik["joint"])
  Q <- function(x) raw(x) - unname(x$par_s$mu$s["b"])^2 / 2
  attempt <- function(objective, backtracking = TRUE, proposed = mk(1)) {
    .gl_rs_accept_backfitting_step(proposed, mk(2), 1, backtracking, 4,
      function(alpha) mk(2 - alpha), 1, 1, 1, "mu", objective_fn = objective)
  }
  old_target <- attempt(raw)
  corrected <- attempt(Q)
  expect_true(old_target$rejected)
  expect_equal(old_target$trace_row$accepted_change, 0)
  expect_false(corrected$rejected)
  expect_equal(corrected$accepted_results$par_s$mu$s, c(b = 1))
  expect_equal(corrected$trace_row$accepted_change, -0.5)
  expect_equal(corrected$change_objective, 1)
  expect_true(attempt(raw, FALSE)$rejected)
  invalid <- mk(1); invalid$calc_lik_out_end$valid <- FALSE
  expect_true(attempt(Q, FALSE, invalid)$rejected)
  expect_equal(Q(mk(1.5)) - Q(mk(2)), 0.75)
})

test_that("RS rebases baseline, trials and rejected EDF to selected lambda", {
  X <- matrix(1, 1, 1, dimnames = list(NULL, "mu.s.1"))
  cache <- list(mu = list(X = X, fixed_names = character(),
    smooth_penalty_meta = list(s = list(idx = 1L, S_base = matrix(1), B = X))))
  mk <- function(b) list(par_cov = numeric(), par_s = list(mu = list(s = c(mu.s.1 = b))),
    calc_lik_out_end = list(log_lik = c(joint = -(b - 2)^2 / 2), valid = TRUE),
    df_s = list(mu = list(s = 0.25)))
  trial_lambdas <- list()
  out <- .gl_run_rs_backfitting_acceptance(
    lambda_s = list(mu = list(s = 1)), par_name = "mu",
    par_s = mk(1)$par_s, par_cov = numeric(), beta_start = c(mu.s.1 = 1), K = 2,
    margin_dist = NULL, copula_dist = NULL, dataset = data.frame(), mm = list(),
    copula_link = NULL, df_s = list(mu = list(s = 0.5)), step_size = 1,
    rs_update_lambda = TRUE, inner_run_counter = 1L, outer_only_run_counter = 2L,
    outer_run_counter = 1L, verbose = 0, backfitting_fn = function(lambda_s, ...) {
      trial_lambdas[[length(trial_lambdas) + 1L]] <<- lambda_s
      mk(3) # deliberately downhill under the selected lambda
    }, calc_lik_out = mk(1)$calc_lik_out_end, rs_block_trace = list(), run_counter = 1L,
    timer_start = Sys.time(), use_backtracking = TRUE, backtracking_max_halves = 2L,
    plot_results = FALSE, log_lik_history = matrix(numeric(), 0, 3),
    par_history = matrix(numeric(), 0, 0), true_val = NULL, penalty_cache = cache,
    prepare_proposal_fn = function(...) list(lambda_s = list(mu = list(s = 3)),
      proposed_results = mk(3),
      lambda_search = list(status = "selected", searches = list(list(status = "selected")))))
  expect_equal(out$lambda_s$mu$s, 3)
  expect_equal(out$par_s$mu$s, c(mu.s.1 = 1))
  expect_equal(out$df_s$mu$s, 0.25)
  expect_equal(out$coefficient_step_status, "rejected")
  row <- out$rs_block_trace[[1]]
  expect_equal(row$start_objective, -2) # not the old-lambda value -1
  expect_true(row$lambda_changed)
  expect_true(row$rejected)
  expect_true(all(vapply(trial_lambdas, function(x) identical(x$mu$s, 3), logical(1))))
})

test_that("RS lambda search failures never install an invalid candidate", {
  lambda <- list(mu = list(s = 1))
  run <- function(code, candidate = 3, value = 1) .gl_rs_update_smoothing_parameters(
    lambda, "mu", list(), numeric(), numeric(), 2, NULL, NULL, data.frame(), list(),
    NULL, list(), 1, TRUE, 1L, 0, function(...) list(),
    optim_fn = function(...) list(par = candidate, convergence = code,
                                  value = value, message = "fixture"), details = TRUE)
  for (x in list(run(1L), run(0L, Inf), run(0L, 1e7), run(0L, 3, NA_real_))) {
    expect_equal(x$lambda_s, lambda)
    expect_equal(x$status, "failed")
    expect_equal(x$searches[[1]]$message, "fixture")
  }
  expect_equal(run(0L)$lambda_s$mu$s, 3)
})

.rs_conditional_runner_fixture <- function(max_outer = 40L, adaptive = FALSE,
    stop_early = TRUE, reject = FALSE, inner_limit = FALSE, search_failure = FALSE,
    reject_once = FALSE) {
  B <- matrix(1, 1, 1); attr(B, "penalty") <- matrix(1)
  mm <- list(x = list(mu = matrix(1, 1, 1, dimnames = list(NULL, "a")),
                          theta = matrix(numeric(), 1, 0)),
             s = list(mu = list(), theta = list(s = B)))
  initial_s <- list(mu = list(), theta = list(s = c(theta.s.1 = 0)))
  cache <- .gl_build_rs_design_cache(mm, initial_s)
  ell <- function(a, b) unname(-((a + b - 3)^2 + (a - 1)^2) / 2)
  search_calls <- 0L
  fixed_starts <- list()
  final_refreshes <- 0L
  acceptance_trace <- list()
  out <- .gl_run_rs_optimizer(
    first_outer_run = TRUE, outer_log_lik_change = 0, outer_stop_crit = 1e-12,
    outer_only_run_counter = 1L, max_outer_iter = max_outer,
    check_elapsed_budget = function(stage) NULL, verbose = 0, mm = mm,
    rs_calc_eta = function(...) stop("synthetic runner must not evaluate a distribution"),
    par_cov = c(mu.a = 0), par_s = initial_s, margin_dist = NULL, copula_dist = NULL,
    dataset = data.frame(), pair_cache = NULL, margin_eval_cache = NULL,
    outer_start_log_lik = 0, log_lik_history = matrix(numeric(), 0, 3),
    par_history = matrix(numeric(), 0, 1), copula_link = NULL,
    discrete_score_method = "analytical", include_dlcopdpar = TRUE,
    check_dlcopdpar_gradient = FALSE, rs_design_cache = cache,
    rs_smooth_trust_radius = Inf, lambda_penalty_K = 2,
    lambda_s = list(mu = list(), theta = list(s = 1)), df_s = list(),
    step_size = 1, rs_update_lambda = adaptive, outer_run_counter = 1L,
    rs_block_trace = list(), use_backtracking = TRUE, backtracking_max_halves = 4L,
    plot_results = FALSE, true_val = NULL, inner_stop_crit = 1e-12,
    max_inner_iter = 10L, outer_negative_streak = 0L, step_adjustment = 1,
    max_steps = 1L, start_step_size = 1, max_negative_outer_streak = 1L,
    weights_final = list(), stop_on_convergence = stop_early,
    parameter_iterations_fn = function(par_cov, par_s, lambda_s, rs_update_lambda,
      outer_only_run_counter, rs_block_trace, outer_run_counter, ...) {
      a <- unname(par_cov["mu.a"]); b <- unname(par_s$theta$s[1])
      start <- ell(a, b); searches <- list()
      rejected_b <- FALSE
      rejected_joint <- FALSE
      if (rs_update_lambda) {
        if (outer_only_run_counter > 1L) {
          search_calls <<- search_calls + 1L
          if (!search_failure) lambda_s$theta$s <- 3
          searches <- list(list(status = if (search_failure) "failed" else "selected"))
        }
      } else {
        fixed_starts[[length(fixed_starts) + 1L]] <<- c(a = a, b = b)
        rejected_b <- reject_once && length(fixed_starts) == 1L
        if (!reject) {
          # The injected analytic sweep must obey the real acceptance contract,
          # including when continued sweeps reach floating-point objective ties.
          candidate_a <- (4 - b) / 2
          candidate_b <- if (rejected_b) b else
            (3 - candidate_a) / (1 + lambda_s$theta$s)
          proposal <- function(alpha) {
            next_cov <- par_cov; next_s <- par_s
            next_cov["mu.a"] <- a + alpha * (candidate_a - a)
            next_s$theta$s[1] <- b + alpha * (candidate_b - b)
            list(par_cov = next_cov, par_s = next_s,
              calc_lik_out_end = list(log_lik = c(joint = ell(
                next_cov["mu.a"], next_s$theta$s[1])), valid = TRUE))
          }
          Q <- function(result) .gl_rs_objective(result$calc_lik_out_end,
            result$par_s, lambda_s, cache)
          accepted <- .gl_rs_accept_backfitting_step(proposal(1), proposal(0),
            1, TRUE, 4L, proposal, outer_only_run_counter, 1L, outer_run_counter,
            "joint_fixture", objective_fn = Q)
          row <- accepted$trace_row
          row$lambda <- lambda_s$theta$s
          row$start_a <- a; row$start_b <- b
          row$accepted_a <- unname(accepted$accepted_results$par_cov["mu.a"])
          row$accepted_b <- unname(accepted$accepted_results$par_s$theta$s[1])
          row$accepted_likelihood_valid <- isTRUE(
            accepted$accepted_results$calc_lik_out_end$valid)
          acceptance_trace[[length(acceptance_trace) + 1L]] <<- row
          rs_block_trace[[length(rs_block_trace) + 1L]] <- row
          rejected_joint <- isTRUE(accepted$rejected)
          a <- unname(accepted$accepted_results$par_cov["mu.a"])
          b <- unname(accepted$accepted_results$par_s$theta$s[1])
        }
      }
      par_cov["mu.a"] <- a; par_s$theta$s[1] <- b
      list(par_cov = par_cov, par_s = par_s, lambda_s = lambda_s, df_s = list(),
        calc_lik_out_end = list(log_lik = c(joint = ell(a, b)), valid = TRUE),
        rs_block_trace = rs_block_trace, outer_run_counter = outer_run_counter + 1L,
        weights_final = list(), first_outer_run = FALSE, outer_start_log_lik = start,
        log_lik_history = matrix(numeric(), 0, 3), par_history = matrix(numeric(), 0, 1),
        block_status = list(mu = if (reject || rejected_joint) "rejected" else if (inner_limit)
          "max_iterations" else "objective_plateau",
          theta = if (rejected_b || rejected_joint) "rejected" else "objective_plateau"),
        complete_sweep = TRUE, lambda_searches = searches,
        lambda_search_failed = search_failure && length(searches) > 0L)
    },
    final_weights_fn = function(par_cov, par_s, ...) {
      final_refreshes <<- final_refreshes + 1L
      list(weights = list(mu = 1, theta = 1),
        likelihood = list(log_lik = c(joint = ell(par_cov[1], par_s$theta$s[1])), valid = TRUE))
    })
  list(out = out, searches = search_calls, starts = fixed_starts, refreshes = final_refreshes,
    acceptance_trace = acceptance_trace)
}

test_that("RS settles all coupled blocks at frozen lambda within one total budget", {
  fit <- .rs_conditional_runner_fixture()
  out <- fit$out
  expect_true(out$rs_protocol$objective_plateau)
  expect_equal(unname(out$par_cov), 5 / 3, tolerance = 1e-6)
  expect_equal(unname(out$par_s$theta$s), 2 / 3, tolerance = 1e-6)
  a <- unname(out$par_cov); b <- unname(out$par_s$theta$s)
  expect_equal(c(4 - 2*a - b, 3 - a - 2*b), c(0, 0), tolerance = 2e-6)
  expect_equal(solve(matrix(c(2, 1, 1, 2), 2)), matrix(c(2, -1, -1, 2), 2) / 3)
  expect_equal(fit$searches, 0L)
  expect_equal(fit$refreshes, 1L)
  expect_true(is.na(out$rs_protocol$lambda_converged))
  expect_equal(out$rs_protocol$stationarity_status, "not_checked")
  adaptive <- .rs_conditional_runner_fixture(adaptive = TRUE)
  expect_equal(adaptive$searches, 1L)
  expect_equal(adaptive$starts[[1]], c(a = 0, b = 0))
  expect_equal(adaptive$out$lambda_s$theta$s, 3)
  expect_equal(unname(adaptive$out$par_cov), 13 / 7, tolerance = 1e-6)
  expect_equal(unname(adaptive$out$par_s$theta$s), 2 / 7, tolerance = 1e-6)
  expect_true(adaptive$out$rs_protocol$objective_plateau)
  short <- .rs_conditional_runner_fixture(max_outer = 2L, adaptive = TRUE)
  expect_false(short$out$rs_protocol$objective_plateau)
  expect_equal(short$out$rs_protocol$stop_reason, "max_iterations")
  expect_length(short$starts, 0L)
  continued <- .rs_conditional_runner_fixture(max_outer = 35L, stop_early = FALSE)
  trace <- do.call(rbind, continued$acceptance_trace)
  protocol <- continued$out$rs_protocol
  n <- nrow(trace)
  expect_lte(n, 35L)
  expect_length(protocol$sweeps, n)
  expect_length(continued$starts, n)
  expect_equal(continued$out$outer_only_run_counter, n + 1L)
  expect_true(all(trace$accepted_likelihood_valid))
  expect_true(all(is.finite(trace$start_logLik) & is.finite(trace$accepted_logLik)))
  expect_true(all(is.finite(trace$start_objective) & is.finite(trace$accepted_objective)))
  expect_true(all(trace$accepted_objective >= trace$start_objective))
  expect_true(all(trace$lambda == 1))
  expect_true(all(vapply(protocol$sweeps,
    function(x) identical(x$phase, "fixed_lambda"), logical(1))))
  # FALSE suppresses successful convergence stopping, while safeguards remain.
  successful <- which(vapply(protocol$sweeps,
    function(x) isTRUE(x$objective_plateau), logical(1)))
  expect_gt(length(successful), 0L)
  first_plateau <- successful[[1]]
  expect_equal(first_plateau, length(out$rs_protocol$sweeps))
  expect_identical(continued$acceptance_trace[seq_len(first_plateau)], fit$acceptance_trace)
  expect_gt(n, first_plateau)
  rejected <- which(trace$rejected)
  if (length(rejected)) {
    # One genuine rejected joint attempt gives a zero-change rejected plateau.
    expect_equal(rejected, n)
    expect_equal(protocol$stop_reason, "max_stall")
    expect_equal(protocol$detail_reason, "terminal_block_rejected")
    expect_false(protocol$objective_plateau)
    expect_identical(unlist(protocol$block_status, use.names = FALSE),
      c("rejected", "rejected"))
    last <- trace[n, , drop = FALSE]
    expect_true(last$rejected)
    expect_true(last$proposed_likelihood_valid)
    expect_identical(last$accepted_step_size, 0)
    expect_identical(last$accepted_objective_change, 0)
    expect_identical(last$accepted_objective, last$start_objective)
    expect_equal(last$backtracking_attempts, 4L)
    expect_equal(last$max_backtracking_attempts, 4L)
    expect_identical(c(last$accepted_a, last$accepted_b), c(last$start_a, last$start_b))
    expect_identical(c(a = unname(continued$out$par_cov),
      b = unname(continued$out$par_s$theta$s)), continued$starts[[n]])
  } else {
    expect_equal(n, 35L)
    expect_equal(continued$out$outer_only_run_counter, 36L)
    expect_false(any(trace$rejected))
    expect_true(protocol$objective_plateau)
    expect_equal(protocol$stop_reason, "converged")
    expect_equal(protocol$detail_reason, "conditional_objective_plateau")
  }
  expect_equal(unname(continued$out$par_cov), 5 / 3, tolerance = 1e-6)
  expect_equal(unname(continued$out$par_s$theta$s), 2 / 3, tolerance = 1e-6)
  recovered <- .rs_conditional_runner_fixture(reject_once = TRUE)$out$rs_protocol
  expect_equal(recovered$sweeps[[1]]$block_status$theta, "rejected")
  expect_true(recovered$objective_plateau)
})

test_that("RS rejection, inner limits and failed search cannot become zero-change success", {
  stalled <- .rs_conditional_runner_fixture(reject = TRUE)$out$rs_protocol
  expect_false(stalled$objective_plateau)
  expect_equal(stalled$stop_reason, "max_stall")
  limited <- .rs_conditional_runner_fixture(inner_limit = TRUE)$out$rs_protocol
  expect_false(limited$objective_plateau)
  expect_equal(limited$stop_reason, "max_iterations")
  failed <- .rs_conditional_runner_fixture(adaptive = TRUE, search_failure = TRUE)$out
  expect_equal(failed$lambda_s$theta$s, 1)
  expect_equal(failed$rs_protocol$detail_reason, "lambda_search_failed")
  expect_false(failed$rs_protocol$objective_plateau)
  empty <- .rs_conditional_runner_fixture(max_outer = 0L, adaptive = TRUE)
  expect_equal(empty$out$rs_protocol$stop_reason, "max_iterations")
  expect_equal(empty$out$rs_protocol$settling_status, "not_started")
  expect_false(empty$out$rs_protocol$objective_plateau)
  expect_false(empty$out$rs_protocol$final_edf_refreshed)
  expect_length(empty$out$rs_protocol$sweeps, 0L)
  expect_equal(empty$refreshes, 0L)
})

test_that("RS raw decline is diagnostic when the comparable conditional objective rises", {
  x <- .gl_update_rs_outer_iteration_state(
    list(log_lik = c(joint = -0.5)), 0, 1L, 0L, 1, 1L, 1, 1L, 1e-6, 0,
    objective_change = 1, objective_kind = "conditional_penalized",
    comparable = TRUE, report_convergence = FALSE)
  expect_equal(unname(x$outer_log_lik_change), -0.5)
  expect_equal(x$outer_negative_streak, 0L)
})

test_that("RS objective plateau does not certify independent stationarity", {
  b <- 2 - 1e-9
  mk <- function(b) list(par_s = list(mu = list(s = b)),
    calc_lik_out_end = list(log_lik = c(joint = -(b - 2)^2 / 2), valid = TRUE))
  Q <- function(x) x$calc_lik_out_end$log_lik[["joint"]] - x$par_s$mu$s^2 / 2
  tiny <- .gl_rs_accept_backfitting_step(mk(b), mk(2), 1e-9, FALSE, 0L,
    function(alpha) stop("single try"), 1L, 1L, 1L, "mu", objective_fn = Q)
  expect_false(tiny$rejected)
  expect_gt(tiny$change_objective, 0)
  expect_lt(tiny$change_objective, 1e-6)
  protocol <- list(objective = Q(mk(b)), objective_change = tiny$change_objective, objective_plateau = TRUE,
    stop_reason = "converged", target = "joint_loglik_minus_half_quadratic_penalty")
  info <- .gl_build_convergence_info("RS", 0, 1e-6, 2L, 10L,
    NA_character_, NA_real_, NA_real_, -Inf, NA_integer_, NA_real_, 10,
    "forward", "analytical", "analytical", objective = mk(b)$calc_lik_out_end$log_lik[["joint"]],
    rs_protocol = protocol)
  expect_true(info$converged)
  expect_equal(info$objective, Q(mk(b)))
  expect_equal(info$logLik, mk(b)$calc_lik_out_end$log_lik[["joint"]])
  expect_true(is.na(info$gradient_criterion_met))
  H <- matrix(-1, 1, 1, dimnames = list("b", "b"))
  P <- matrix(1, 1, 1, dimnames = list("b", "b"))
  expect_warning(inference <- .gl_joint_solve(H, P, c(b = b), c(b = 2 - b),
    .gl_normalize_inference_control(inference_control("standard"))))
  expect_equal(inference$hessian_diagnostics$stationarity_status, "failed")
  expect_equal(inference$hessian_diagnostics$stationarity_scaled_gradient_max,
    (2 * b - 2) / sqrt(2))
})

test_that("RS final weights use one common endpoint without coefficient updates", {
  seen <- list()
  eta <- list(eta = list(mu = c(2, 2), theta = c(3, 3)),
    eta_inv = list(mu = c(2, 2), theta = c(3, 3)),
    eta_dr = list(mu = c(1, 1), theta = c(1, 1)))
  out <- .gl_refresh_rs_final_weights(c(mu.a = 2), list(),
    list(x = list(mu = matrix(1, 2), theta = matrix(1, 2))),
    function(par_cov_current, par_s_current) eta, NULL, NULL,
    data.frame(response = 1:2, time = 1:2, subject = 1L), NULL, NULL, NULL,
    "analytical", TRUE, function(stage) NULL,
    likelihood_fn = function(...) list(log_lik = c(joint = -1), valid = TRUE),
    score_state_fn = function(par_name, eta, ...) {
      seen[[par_name]] <<- eta
      list(d1 = if (par_name == "mu") c(2, 2) else c(3, 3))
    })
  expect_equal(out$weights, list(mu = c(4, 4), theta = c(9, 9)))
  expect_identical(seen$mu, eta$eta)
  expect_identical(seen$theta, eta$eta)
  expect_equal(out$likelihood$log_lik, c(joint = -1))
})

test_that("invalid raw or conditional objectives override claimed RS convergence reasons", {
  make_info <- function(raw = 0, Q = -2, change = 0) .gl_build_convergence_info(
    "RS", 0, 1e-6, 2L, 10L, NA_character_, NA_real_, NA_real_, -Inf,
    NA_integer_, NA_real_, 10, "forward", "analytical", "analytical", objective = raw,
    rs_protocol = list(objective = Q, objective_change = change, objective_plateau = TRUE,
      stop_reason = "converged", target = "joint_loglik_minus_half_quadratic_penalty"))
  invalid_raw <- make_info(raw = NA_real_)
  expect_false(invalid_raw$converged)
  expect_equal(invalid_raw$stop_reason, "invalid_likelihood")
  invalid_Q <- make_info(Q = NA_real_)
  expect_false(invalid_Q$converged)
  expect_equal(invalid_Q$stop_reason, "numerical_failure")
  inconsistent <- make_info(change = 1)
  expect_false(inconsistent$converged)
  expect_equal(inconsistent$stop_reason, "numerical_failure")
})

test_that("composed RS carry path freezes lambda and preserves rejection and total budget", {
  # Real runner -> parameter loop -> inner iteration -> backfitting -> acceptance.
  # Only scalar likelihood, raw score, coefficient proposal and lambda candidate
  # are analytic injections. This is not a GAMLSS distribution fit.
  run_composed <- function(max_outer, reject_selected = FALSE) {
    selected <- 0L
    in_search <- FALSE
    proposals <- list()
    budget <- character()
    B <- matrix(1, 1, 1); attr(B, "penalty") <- matrix(1)
    mm <- list(x = list(mu = matrix(numeric(), 1, 0)), s = list(mu = list(s = B)))
    par_s <- list(mu = list(s = c(mu.s.1 = 1)))
    cache <- .gl_build_rs_design_cache(mm, par_s)
    eta_fn <- function(par_cov_current, par_s_current) {
      b <- unname(par_s_current$mu$s)
      list(eta = list(mu = b), eta_inv = list(mu = b), eta_dr = list(mu = 1))
    }
    scalar_likelihood <- function(eta_inv, ...) {
      b <- unname(eta_inv$mu)
      ll <- -(b - 2)^2 / 2
      list(log_lik = c(marginal = ll, copula = 0, joint = ll), valid = TRUE,
        margin_deriv = list(dldm = 2 - b), margin_d = 1, margin_p = 0.5,
        copula_d = 1, copula_p = 0.5, Fx_1_2 = matrix(0.5, 1, 2),
        order_copula = 1L)
    }
    scalar_score <- function(par_name, eta, timer, ...) {
      list(d1 = 2 - eta$mu, timer = timer)
    }
    likelihood_state <- function(...) .gl_evaluate_rs_iteration_likelihood_state(
      ..., likelihood_context_fn = function(current_calc_lik_out = NULL, ...) {
        # Always use the injected analytic likelihood, including later inner
        # iterations; do not enter the family's cached derivative evaluator.
        .gl_rs_likelihood_context(..., current_calc_lik_out = NULL,
          likelihood_fn = scalar_likelihood)
      })
    proposal_runner <- function(design_info, w_k_vec, ...) {
      function(par_s, par_cov, beta_start, lambda_s, step_size, K, df_s, ...) {
        b <- unname(beta_start)
        lambda <- lambda_s$mu$s
        proposals[[length(proposals) + 1L]] <<- list(
          start = b, lambda = lambda, search = in_search, step = step_size)
        target <- if (reject_selected && lambda == 3) b + 1 else 2 / (1 + lambda)
        next_b <- b + step_size * (target - b)
        next_s <- par_s
        next_s$mu$s <- stats::setNames(next_b, names(par_s$mu$s))
        likelihood <- scalar_likelihood(eta_fn(par_cov, next_s)$eta_inv)
        meta <- design_info$smooth_penalty_meta$s
        df_s$mu$s <- .gl_rs_smooth_term_edf(meta$B, w_k_vec, meta$S_base, lambda)
        list(par_cov = par_cov, par_s = next_s, calc_lik_out_end = likelihood,
          df_s = df_s, GAIC_lambda_k = -2 * likelihood$log_lik[["joint"]] + K * df_s$mu$s)
      }
    }
    known_lambda <- function(par, fn, ...) {
      selected <<- selected + 1L
      in_search <<- TRUE
      on.exit(in_search <<- FALSE)
      value <- fn(3)
      list(par = 3, convergence = 0L, value = value, message = "analytic candidate")
    }
    acceptance <- function(...) .gl_run_rs_backfitting_acceptance(
      ..., prepare_proposal_fn = function(...) .gl_rs_prepare_backfitting_proposal(
        ..., lambda_update_fn = function(..., details = FALSE) {
          .gl_rs_update_smoothing_parameters(..., optim_fn = known_lambda, details = details)
        }))
    backfitting <- function(...) .gl_run_rs_inner_backfitting_step(
      ..., runner_fn = proposal_runner, acceptance_fn = acceptance)
    out <- .gl_run_rs_optimizer(
      first_outer_run = TRUE, outer_log_lik_change = 0, outer_stop_crit = 1e-10,
      outer_only_run_counter = 1L, max_outer_iter = max_outer,
      check_elapsed_budget = function(stage) budget <<- c(budget, stage),
      verbose = 0, mm = mm, rs_calc_eta = eta_fn,
      par_cov = numeric(), par_s = par_s, margin_dist = NULL, copula_dist = NULL,
      dataset = data.frame(response = 2, time = 1, subject = 1),
      pair_cache = NULL, margin_eval_cache = NULL,
      outer_start_log_lik = 0, log_lik_history = matrix(numeric(), 0, 3),
      par_history = matrix(numeric(), 0, 0), copula_link = NULL,
      discrete_score_method = "analytical", include_dlcopdpar = TRUE,
      check_dlcopdpar_gradient = FALSE, rs_design_cache = cache,
      rs_smooth_trust_radius = Inf, lambda_penalty_K = 2,
      lambda_s = list(mu = list(s = 1)), df_s = list(mu = list(s = 0)),
      step_size = 1, rs_update_lambda = TRUE, outer_run_counter = 1L,
      rs_block_trace = list(), use_backtracking = TRUE, backtracking_max_halves = 2L,
      plot_results = FALSE, true_val = NULL, inner_stop_crit = 1e-10,
      max_inner_iter = 4L, outer_negative_streak = 0L, step_adjustment = 1,
      max_steps = 1L, start_step_size = 1, max_negative_outer_streak = 1L,
      weights_final = list(), likelihood_state_fn = likelihood_state,
      score_state_fn = scalar_score, backfitting_step_fn = backfitting,
      final_weights_fn = function(...) .gl_refresh_rs_final_weights(
        ..., likelihood_fn = scalar_likelihood, score_state_fn = scalar_score))
    list(out = out, searches = selected, proposals = proposals, budget = budget)
  }

  accepted <- run_composed(6L)
  x <- accepted$out
  rows <- do.call(rbind, x$rs_block_trace)
  expect_equal(x$par_s$mu$s, c(mu.s.1 = 0.5)) # exact optimum 2/(1+3)
  expect_equal(x$lambda_s$mu$s, 3)
  expect_true(x$rs_protocol$objective_plateau)
  expect_equal(x$rs_protocol$selection_stop_reason, "raw_likelihood_plateau")
  expect_equal(x$rs_protocol$settling_status, "objective_plateau")
  expect_equal(x$outer_only_run_counter, 5L) # 3 selection sweeps + 1 frozen sweep
  expect_equal(accepted$searches, 2L)
  expect_equal(x$rs_protocol$lambda_epoch, 1L)
  selected_row <- rows[rows$outer_iteration == 2L & rows$inner_iteration == 1L, ]
  expect_equal(selected_row$start_objective, -2) # rebase at lambda3, never old -1
  expect_equal(selected_row$accepted_objective, -1.5)
  expect_equal(selected_row$accepted_change, -0.625) # legitimate raw decline
  expect_equal(selected_row$accepted_objective_change, 0.5)
  carried <- rows[rows$outer_iteration == 2L & rows$inner_iteration == 2L, ]
  expect_equal(carried$start_objective, -1.5)
  expect_equal(carried$accepted_objective_change, 0)
  expect_true(any(vapply(accepted$proposals, function(p) {
    !p$search && p$lambda == 3 && p$start == 0.5
  }, logical(1))))
  frozen <- rows[rows$phase == "fixed_lambda", ]
  expect_equal(nrow(frozen), 1L)
  expect_true(all(frozen$lambda_search_status == "not_searched"))
  expect_true(all(vapply(frozen$lambda_used, function(z) identical(z$mu$s, 3), logical(1))))
  expect_equal(x$weights_final$mu, 2.25)
  expect_equal(unname(x$df_s$mu$s), 3 / 7) # final common-state weight/(weight+lambda)
  expect_equal(sum(accepted$budget == "RS outer iteration"), 4L)
  expect_equal(sum(accepted$budget == "RS inner iteration"), nrow(rows))

  exhausted <- run_composed(3L)
  expect_equal(exhausted$out$outer_only_run_counter, 4L)
  expect_equal(exhausted$out$rs_protocol$stop_reason, "max_iterations")
  expect_equal(exhausted$out$rs_protocol$settling_status, "not_started")
  expect_false(exhausted$out$rs_protocol$objective_plateau)
  expect_equal(sum(exhausted$budget == "RS outer iteration"), 3L)

  rejected <- run_composed(6L, reject_selected = TRUE)
  y <- rejected$out
  rejection_rows <- do.call(rbind, y$rs_block_trace)
  expect_equal(y$par_s$mu$s, c(mu.s.1 = 1))
  expect_equal(y$lambda_s$mu$s, 3)
  expect_equal(rejected$searches, 1L) # no new search in the frozen end phase
  expect_equal(y$rs_protocol$stop_reason, "max_stall")
  expect_equal(y$rs_protocol$detail_reason, "terminal_block_rejected")
  expect_equal(y$rs_protocol$block_status$mu, "rejected")
  expect_false(y$rs_protocol$objective_plateau)
  terminal <- tail(rejection_rows, 1)
  expect_equal(terminal$phase, "fixed_lambda")
  expect_true(terminal$rejected)
  expect_equal(terminal$accepted_change, 0)
  expect_equal(terminal$accepted_objective_change, 0)
  expect_equal(terminal$start_objective, -2)
  expect_equal(y$weights_final$mu, 1)
  expect_equal(unname(y$df_s$mu$s), 0.25) # selected lambda remains in EDF after rejection
})
