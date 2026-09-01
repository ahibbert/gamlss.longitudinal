test_that("fit workflow preparation sequences setup phases and returns bundled state", {
  calls <- character()
  captured <- list()

  record <- function(label) {
    calls <<- c(calls, label)
  }

  dataset_in <- data.frame(y = c(1, 2), id = c(1, 1), t = c(1, 2))
  dataset_prepared <- data.frame(response = c(1, 2), subject = c(1, 1), time = c(1, 2))
  formulas_int <- list(
    mu = response ~ 1,
    sigma = ~ 1,
    nu = ~ 1,
    tau = ~ 1,
    theta = ~ 1,
    zeta = ~ 1
  )
  mm <- list(x = list(mu = matrix(1, nrow = 2, ncol = 1)), s = list(mu = list()))
  matrix_bundle <- list(mm = mm, copula_link = "identity")
  optimizer_context <- list(par_cov = c(`mu.(Intercept)` = 1), marker = "optimizer")
  optimizer_control <- gamlss_longitudinal_control(
    rs = list(
      warm_start_joint_iter = 3L,
      backtracking_max_halves = 7L,
      discrete_score_method = "finite",
      smooth_trust_radius = 13
    ),
    cg = list(
      max_stall = 8L,
      max_delta = 0.25,
      lambda_update_every = 9L,
      max_lambda_updates = 10L,
      raw_loglik_drop_tol = 11,
      line_search = "first",
      max_line_search_evals = 12L,
      gradient_method = "central",
      zeta_hessian = "finite",
      hessian_method = "auto"
    )
  )

  out <- gamlss.longitudinal:::.gl_prepare_fit_workflow(
    dataset = dataset_in,
    margin_dist = "NO",
    copula_dist = "N",
    time_var = "t",
    subject_var = "id",
    missingness = "error",
    mu.formula = y ~ x,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    tau.formula = ~ 1,
    theta.formula = ~ 1,
    zeta.formula = ~ 1,
    include_dlcopdpar = TRUE,
    start_from = NA,
    verbose = 0,
    true_val = NULL,
    method = "RS",
    optimizer_control = optimizer_control,
    vcov_method = "analytical",
    vcov_numderiv = FALSE,
    use_Rcpp = FALSE,
    lambda_start = 2,
    lambda_penalty_K = 2,
    data_fn = function(...) {
      record("data")
      captured$data <<- list(...)
      list(
        dataset_original = dataset_in,
        dataset = dataset_prepared,
        response_var = "y",
        var_map = list(response = "y"),
        formulas_int = formulas_int,
        miss_by_time = data.frame(),
        pair_summary = data.frame()
      )
    },
    matrix_fn = function(...) {
      record("matrix")
      captured$matrix <<- list(...)
      matrix_bundle
    },
    warm_start_fn = function(...) {
      record("warm-start")
      captured$warm_start <<- list(...)
      list(
        start_from = c(`mu.(Intercept)` = 4),
        warm_start_par_s = list(mu = list()),
        warm_start_info = list(used = TRUE)
      )
    },
    step_control_fn = function(...) {
      record("step-controls")
      captured$step_controls <<- list(...)
      list(start_step_size = 0.75, max_steps = 6L, step_adjustment = 0.5)
    },
    optimizer_context_fn = function(...) {
      record("optimizer-context")
      captured$optimizer_context <<- list(...)
      optimizer_context
    }
  )

  expect_equal(calls, c("data", "matrix", "warm-start", "step-controls", "optimizer-context"))
  expect_equal(captured$warm_start$method, "RS")
  expect_equal(captured$warm_start$warm_start_joint_iter, 3L)
  expect_identical(captured$warm_start$optimizer_control, optimizer_control)
  expect_equal(captured$warm_start$optimizer_control$rs$backtracking_max_halves, 7L)
  expect_equal(captured$warm_start$optimizer_control$cg$line_search, "first")
  expect_equal(captured$warm_start$optimizer_control$rs$discrete_score_method, "finite")
  expect_equal(captured$warm_start$vcov_method, "analytical")
  expect_identical(captured$matrix$formulas_int, formulas_int)
  expect_identical(captured$matrix$dataset, dataset_prepared)
  expect_equal(captured$step_controls$method, "RS")
  expect_equal(captured$optimizer_context$start_from, c(`mu.(Intercept)` = 4))
  expect_identical(captured$optimizer_context$warm_start_par_s, list(mu = list()))
  expect_identical(captured$optimizer_context$mm, mm)
  expect_equal(captured$optimizer_context$start_step_size, 0.75)
  expect_equal(captured$optimizer_context$copula_link, "identity")

  expect_equal(out$controls$method, "RS")
  expect_identical(out$fit_data$dataset, dataset_prepared)
  expect_identical(out$matrix_bundle, matrix_bundle)
  expect_true(out$warm_start$warm_start_info$used)
  expect_equal(out$step_controls$max_steps, 6L)
  expect_identical(out$optimizer_context, optimizer_context)
})
