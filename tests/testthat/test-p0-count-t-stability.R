test_that("t-copula startup screens complex margins on the response scale", {
  captured_margin <- NULL
  testthat::local_mocked_bindings(
    .starting_margin_parameter_values = function(...) {
      list(
        margin_par = c(mu = log(4), sigma = log(0.7), nu = stats::qlogis(0.2)),
        margin_par_already_eta = TRUE
      )
    },
    .select_t_copula_zeta_start = function(..., margin_par) {
      captured_margin <<- margin_par
      4
    },
    .package = "gamlss.longitudinal"
  )

  dat <- data.frame(
    response = c(0, 1, 2, 3, 1, 4),
    time = rep(1:3, 2),
    subject = rep(1:2, each = 3)
  )
  start <- gamlss.longitudinal:::get_starting_values(
    copula_dist = "t",
    margin_dist = gamlss.dist::ZIPIG(),
    dataset = dat,
    eta_transform = TRUE
  )

  expect_equal(unname(captured_margin), c(4, 0.7, 0.2), tolerance = 1e-12)
  expect_equal(unname(start[c("mu", "sigma", "nu")]),
               c(log(4), log(0.7), stats::qlogis(0.2)))
})

test_that("copula likelihood broadcasts scalar parameters and skips invalid pairs", {
  evaluated <- integer()
  testthat::local_mocked_bindings(
    .copula_rectangle_prob = function(u1, u2, l1, l2, ...) {
      evaluated <<- c(evaluated, length(u1))
      (u1 - l1) * (u2 - l2)
    },
    .package = "gamlss.longitudinal"
  )

  response <- c(0, 1, 2, 3)
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
    response = response,
    response_margin = rep(1:2, 2),
    response_subject = rep(1:2, each = 2)
  )
  out <- gamlss.longitudinal:::.gl_likelihood_evaluate_copula_pairs(
    eta_inv = list(theta = 0.4, zeta = 4),
    pair_cache = pair_cache,
    copula_dist = "t",
    margin_p = c(0.2, 0.4, NA_real_, 0.8),
    margin_d = c(0.1, 0.2, NA_real_, 0.2),
    margin_log_d = log(c(0.1, 0.2, NA_real_, 0.2)),
    margin_p_lower = c(0.1, 0.2, NA_real_, 0.6),
    likelihood_type = "discrete_rectangle"
  )

  expect_equal(out$par1, rep(0.4, 2))
  expect_equal(out$par2, rep(4, 2))
  expect_equal(evaluated, 1L)
  expect_true(is.finite(out$copula_log_d[[1L]]))
  expect_true(is.na(out$copula_log_d[[2L]]))
})

test_that("direct Student-t rectangles agree with four-CDF differences", {
  upper1 <- c(0.56, 0.82, 0.95, 0.4)
  upper2 <- c(0.46, 0.84, 0.92, 0.7)
  lower1 <- c(0.54, 0.81, 0.94, 0.35)
  lower2 <- c(0.43, 0.83, 0.90, 0.63)
  rho <- c(0.82, 0.5, -0.3, 0.1)
  degrees_freedom <- c(2.05, 5, 12, 30)

  corner <- gamlss.longitudinal:::.copula_t_cdf(
    upper1, upper2, rho, degrees_freedom
  ) - gamlss.longitudinal:::.copula_t_cdf(
    lower1, upper2, rho, degrees_freedom
  ) - gamlss.longitudinal:::.copula_t_cdf(
    upper1, lower2, rho, degrees_freedom
  ) + gamlss.longitudinal:::.copula_t_cdf(
    lower1, lower2, rho, degrees_freedom
  )
  direct <- gamlss.longitudinal:::.copula_rectangle_prob(
    upper1, upper2, lower1, lower2,
    family = "t", par = rho, par2 = degrees_freedom,
    floor_probability = FALSE
  )

  expect_equal(direct, corner, tolerance = 1e-9)
  expect_true(all(direct > 0))
})

test_that("joint selection isolates and checkpoints candidates", {
  skip_if_not_installed("callr")
  skip_if_not_installed("digest")

  dat <- simulate_longitudinal_dataset(
    n = 6,
    times = 1:3,
    margin_dist = gamlss.dist::PO(),
    margin_params = list(mu = 3),
    copula_dist = "N",
    copula_params = list(theta = 0.2),
    seed = 771
  )
  checkpoint_dir <- withr::local_tempdir()
  fit_args <- list(
    optimizer_control = gamlss_longitudinal_control(
      max_outer_iter = 1,
      rs = list(max_inner_iter = 1, warm_start_joint = FALSE)
    )
  )

  first <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    margin_families = "PO",
    copula_families = "N",
    progress = FALSE,
    fit_args = fit_args,
    isolate_candidates = "always",
    candidate_timeout = 60,
    checkpoint_dir = checkpoint_dir
  )
  messages <- character()
  second <- withCallingHandlers(
    select_joint_distribution(
      data = dat,
      response_var = "response",
      time_var = "time",
      subject_var = "subject",
      margin_families = "PO",
      copula_families = "N",
      progress = TRUE,
      fit_args = fit_args,
      isolate_candidates = "always",
      candidate_timeout = 60,
      checkpoint_dir = checkpoint_dir
    ),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    }
  )

  expect_equal(as.data.frame(second), as.data.frame(first))
  expect_true(any(grepl("Restored from checkpoint", messages, fixed = TRUE)))
  expect_length(list.files(checkpoint_dir, pattern = "[.]rds$"), 1L)
})

test_that("isolated joint candidates return structured timeout rows", {
  skip_if_not_installed("callr")

  dat <- simulate_longitudinal_dataset(
    n = 4,
    times = 1:2,
    margin_dist = gamlss.dist::PO(),
    margin_params = list(mu = 2),
    copula_dist = "N",
    copula_params = list(theta = 0.1),
    seed = 772
  )
  selected <- select_joint_distribution(
    data = dat,
    response_var = "response",
    time_var = "time",
    subject_var = "subject",
    margin_families = "PO",
    copula_families = "N",
    min_pairs = 1,
    progress = FALSE,
    isolate_candidates = "always",
    candidate_timeout = 1e-3
  )

  expect_identical(selected$stop_reason, "candidate_timeout")
  expect_match(selected$error, "exceeded")
  expect_false(selected$selection_eligible)
})
