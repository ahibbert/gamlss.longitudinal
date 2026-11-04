test_that("original formula bundle preserves names and formula objects", {
  out <- gamlss.longitudinal:::.gl_original_formula_bundle(
    mu.formula = y ~ x,
    sigma.formula = ~ z,
    nu.formula = ~ 1,
    tau.formula = ~ w,
    theta.formula = ~ a,
    zeta.formula = ~ b
  )

  expect_named(out, c("mu", "sigma", "nu", "tau", "theta", "zeta"))
  expect_equal(out$mu, y ~ x)
  expect_equal(out$sigma, ~ z)
  expect_equal(out$nu, ~ 1)
  expect_equal(out$tau, ~ w)
  expect_equal(out$theta, ~ a)
  expect_equal(out$zeta, ~ b)
})

test_that("fit entrypoint sequences public fit phases and preserves original fit arguments", {
  calls <- character()
  captured <- list()

  record <- function(label) {
    calls <<- c(calls, label)
  }

  fit_time <- as.POSIXct("2024-01-01 00:00:00", tz = "UTC")
  budget_checker <- function(stage) stage
  fit_workflow <- list(workflow = TRUE)
  optimizer_state <- list(optimizer = TRUE)
  final_fit <- list(class = "gamlss.longitudinal", result = TRUE)
  optimizer_control <- gamlss_longitudinal_control(
    max_outer_iter = 10L,
    max_elapsed_sec = 30,
    rs = list(
      max_inner_iter = 20L,
      max_negative_outer_streak = 3L,
      backtracking_max_halves = 7L
    )
  )

  out <- gamlss.longitudinal:::.gl_run_gamlss_longitudinal_entrypoint(
    dataset = data.frame(y = 1, id = 1, t = 1),
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
    check_dlcopdpar_gradient = FALSE,
    start_from = NA,
    verbose = 0,
    plot_results = FALSE,
    true_val = NA,
    method = "RS",
    optimizer_control = optimizer_control,
    compute_vcov = TRUE,
    vcov_method = "analytical",
    vcov_numderiv = FALSE,
    use_Rcpp = FALSE,
    lambda_start = NA,
    lambda_penalty_K = 2,
    time_fn = function() {
      record("time")
      fit_time
    },
    margin_normalizer_fn = function(margin_dist) {
      record("normalize")
      captured$margin_input <<- margin_dist
      "NO-normalized"
    },
    capability_preflight_fn = function(...) NULL,
    budget_checker_fn = function(fit_start_time, max_elapsed_sec) {
      record("budget")
      captured$budget <<- list(
        fit_start_time = fit_start_time,
        max_elapsed_sec = max_elapsed_sec
      )
      budget_checker
    },
    workflow_fn = function(...) {
      record("workflow")
      captured$workflow <<- list(...)
      fit_workflow
    },
    optimizer_fn = function(...) {
      record("optimizer")
      captured$optimizer <<- list(...)
      optimizer_state
    },
    finalize_fn = function(...) {
      record("finalize")
      captured$finalize <<- list(...)
      final_fit
    }
  )

  expect_equal(calls, c("time", "normalize", "budget", "workflow", "optimizer", "finalize"))
  expect_equal(captured$margin_input, "NO")
  expect_equal(captured$budget$fit_start_time, fit_time)
  expect_equal(captured$budget$max_elapsed_sec, 30)
  expect_equal(captured$workflow$margin_dist, "NO-normalized")
  expect_identical(captured$workflow$optimizer_control, optimizer_control)
  expect_equal(captured$workflow$optimizer_control$rs$max_inner_iter, 20L)
  expect_equal(captured$workflow$optimizer_control$cg$line_search, "best")
  expect_equal(captured$optimizer$fit_workflow, fit_workflow)
  expect_equal(captured$optimizer$check_elapsed_budget, budget_checker)
  expect_false("max_outer_iter" %in% names(captured$optimizer))
  expect_equal(captured$finalize$optimizer_state, optimizer_state)
  expect_equal(captured$finalize$fit_start_time, fit_time)
  expect_equal(captured$finalize$original_formulas$mu, y ~ x)
  expect_equal(captured$finalize$time_var, "t")
  expect_equal(captured$finalize$compute_vcov, TRUE)
  expect_equal(out, final_fit)
})
