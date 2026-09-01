make_jss001_likelihood_inputs <- function(response, discrete = FALSE) {
  n <- length(response)
  stopifnot(n %% 2L == 0L)
  n_subject <- n %/% 2L
  time <- rep(1:2, times = n_subject)
  subject <- rep(seq_len(n_subject), each = 2L)
  if (isTRUE(discrete)) {
    list(
      response = response,
      time = time,
      subject = subject,
      margin_dist = gamlss.dist::PO(),
      mm = list(
        mu = matrix(1, n, 1),
        theta = matrix(1, n_subject, 1)
      ),
      eta_inv = list(mu = rep(2, n), theta = rep(1, n_subject))
    )
  } else {
    list(
      response = response,
      time = time,
      subject = subject,
      margin_dist = gamlss.dist::NO(),
      mm = list(
        mu = matrix(1, n, 1),
        sigma = matrix(1, n, 1),
        theta = matrix(1, n_subject, 1)
      ),
      eta_inv = list(
        mu = rep(0, n),
        sigma = rep(1, n),
        theta = rep(0.25, n_subject)
      )
    )
  }
}

evaluate_jss001_inputs <- function(x, eta_inv = x$eta_inv, pair_cache = NULL,
                                   copula_dist = "C") {
  gamlss.longitudinal:::calc_likelihood_minimal(
    eta_inv = eta_inv,
    mm = x$mm,
    margin_dist = x$margin_dist,
    copula_dist = copula_dist,
    response = x$response,
    response_margin = x$time,
    response_subject = x$subject,
    pair_cache = pair_cache
  )
}

test_that("likelihood inclusion counts depend only on response missingness", {
  skip_if_not_installed("gamlss.dist")
  x <- make_jss001_likelihood_inputs(c(0, 1, NA_real_, 2))
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
    x$response, x$time, x$subject
  )

  valid <- evaluate_jss001_inputs(x, pair_cache = pair_cache)
  invalid_eta <- x$eta_inv
  invalid_eta$sigma[2] <- 0
  invalid <- suppressWarnings(evaluate_jss001_inputs(
    x, eta_inv = invalid_eta, pair_cache = pair_cache
  ))

  expect_equal(
    unname(valid$contribution_counts[c("marginal_included", "pair_included")]),
    c(3, 1)
  )
  expect_equal(
    invalid$contribution_counts[c("marginal_included", "pair_included")],
    valid$contribution_counts[c("marginal_included", "pair_included")]
  )
  expect_identical(valid$margin_included, c(TRUE, TRUE, FALSE, TRUE))
  expect_identical(valid$pair_included, c(TRUE, FALSE))
  expect_true(valid$valid)
  expect_false(invalid$valid)
  expect_identical(unname(invalid$log_lik["joint"]), -Inf)
  expect_true("invalid_margin_contribution" %in% invalid$failure$codes)
})

test_that("continuous marginal and copula terms use stable log evaluation", {
  skip_if_not_installed("gamlss.dist")

  single <- make_jss001_likelihood_inputs(c(40, NA_real_))
  single$eta_inv$theta <- 0
  margin_tail <- evaluate_jss001_inputs(single)
  expect_equal(margin_tail$margin_d[1], 0)
  expect_true(is.finite(margin_tail$margin_log_d[1]))
  expect_true(margin_tail$valid)
  expect_true(is.finite(margin_tail$log_lik["joint"]))

  opposite_tails <- make_jss001_likelihood_inputs(c(-7, 7))
  opposite_tails$eta_inv$theta <- 0.999
  copula_tail <- evaluate_jss001_inputs(opposite_tails, copula_dist = "N")
  expect_equal(copula_tail$copula_d[1], 0)
  expect_true(is.finite(copula_tail$copula_log_d[1]))
  expect_true(copula_tail$valid)
  expect_true(is.finite(copula_tail$log_lik["joint"]))
})

test_that("invalid included copula inputs fail instead of becoming independence", {
  skip_if_not_installed("gamlss.dist")
  x <- make_jss001_likelihood_inputs(c(0, 1))
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
    x$response, x$time, x$subject
  )
  invalid_eta <- x$eta_inv
  invalid_eta$theta[] <- Inf

  out <- evaluate_jss001_inputs(x, eta_inv = invalid_eta, pair_cache = pair_cache)

  expect_identical(out$pair_included, TRUE)
  expect_false(out$pair_input_valid[1])
  expect_true(is.na(out$copula_d[1]))
  expect_false(out$valid)
  expect_identical(unname(out$log_lik["joint"]), -Inf)
  expect_true("invalid_pair_input" %in% out$failure$codes)

  expect_error(
    evaluate_jss001_inputs(
      within(x, response <- c(0, NA_real_)),
      eta_inv = invalid_eta,
      pair_cache = pair_cache
    ),
    "different response missingness pattern"
  )
})

test_that("invalid discrete contributions fail with fixed inclusion counts", {
  skip_if_not_installed("gamlss.dist")
  x <- make_jss001_likelihood_inputs(c(1, 2), discrete = TRUE)
  valid <- evaluate_jss001_inputs(x)
  invalid_eta <- x$eta_inv
  invalid_eta$mu[1] <- 0
  invalid <- suppressWarnings(evaluate_jss001_inputs(x, eta_inv = invalid_eta))

  expect_identical(valid$likelihood_type, "discrete_rectangle")
  expect_true(valid$valid)
  expect_equal(
    unname(valid$contribution_counts[c("marginal_included", "pair_included")]),
    c(2, 1)
  )
  expect_equal(
    invalid$contribution_counts[c("marginal_included", "pair_included")],
    valid$contribution_counts[c("marginal_included", "pair_included")]
  )
  expect_false(invalid$valid)
  expect_identical(unname(invalid$log_lik["joint"]), -Inf)
  expect_true("invalid_margin_contribution" %in% invalid$failure$codes)
})

test_that("CG and RS reject likelihood-invalid proposals consistently", {
  skip_if_not_installed("gamlss.dist")
  dataset <- data.frame(response = c(0, 1), time = 1:2, subject = c("a", "a"))
  mm_cg <- list(x = list(mu = cbind(`(Intercept)` = c(1, 1))), s = list(mu = list()))
  par_cov <- c(`mu.(Intercept)` = 0)

  cg <- gamlss.longitudinal:::.gl_evaluate_cg_beta(
    beta_vec = par_cov,
    mm_cg = mm_cg,
    par_cov_template = par_cov,
    par_s_template = list(mu = list()),
    margin_dist = gamlss.dist::NO(),
    copula_link = "identity",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = list(),
    margin_eval_cache = list(),
    calc_eta_fn = function(...) {
      list(eta_inv = list(mu = c(0, 0), sigma = c(1, 1), theta = 0))
    },
    likelihood_fn = function(...) {
      list(valid = FALSE, log_lik = c(marginal = -1, copula = -Inf, joint = -Inf))
    }
  )
  expect_null(cg)

  make_result <- function(loglik, valid, id) {
    list(
      par_cov = c(value = id),
      par_s = list(mu = list()),
      calc_lik_out_end = list(
        valid = valid,
        log_lik = c(marginal = loglik, copula = 0, joint = loglik)
      ),
      GAIC_lambda_k = NA_real_,
      df_s = list(mu = list())
    )
  }
  current <- make_result(-10, TRUE, 1)
  proposed <- make_result(-Inf, FALSE, 2)

  rs_no_backtrack <- gamlss.longitudinal:::.gl_rs_accept_backfitting_step(
    proposed_results = proposed,
    current_results = current,
    nominal_step_size = 1,
    use_backtracking = FALSE,
    backtracking_max_halves = 0,
    proposal_fn = function(...) stop("not called"),
    outer_iteration = 1,
    inner_iteration = 1,
    global_inner_iteration = 1,
    parameter = "mu"
  )
  expect_true(rs_no_backtrack$rejected)
  expect_equal(rs_no_backtrack$accepted_results$par_cov, current$par_cov)
  expect_false(rs_no_backtrack$trace_row$proposed_likelihood_valid)

  rs_backtrack <- gamlss.longitudinal:::.gl_rs_accept_backfitting_step(
    proposed_results = proposed,
    current_results = current,
    nominal_step_size = 1,
    use_backtracking = TRUE,
    backtracking_max_halves = 2,
    proposal_fn = function(step_size) make_result(-9, TRUE, step_size),
    outer_iteration = 1,
    inner_iteration = 1,
    global_inner_iteration = 1,
    parameter = "mu"
  )
  expect_false(rs_backtrack$rejected)
  expect_true(rs_backtrack$backtracking_applied)
  expect_equal(rs_backtrack$accepted_step_size, 0.5)
})
