test_that("RS design cache records fixed-only design metadata", {
  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1, x = c(0, 1, 2))),
    s = list(mu = list())
  )
  par_s <- list(mu = list())

  cache <- gamlss.longitudinal:::.gl_build_rs_design_cache(mm = mm, par_s = par_s)

  expect_named(cache, "mu")
  expect_equal(colnames(cache$mu$X), c("mu.(Intercept)", "mu.x"))
  expect_equal(cache$mu$fixed_names, c("mu.(Intercept)", "mu.x"))
  expect_equal(cache$mu$smooth_penalty_meta, list())
})

test_that("RS design cache appends smooth blocks and penalty metadata", {
  B <- matrix(c(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 1
  ), nrow = 4, byrow = TRUE)
  penalty <- diag(3)
  attr(B, "penalty") <- penalty

  mm <- list(
    x = list(mu = cbind(`(Intercept)` = 1, x = c(0, 1, 2, 3))),
    s = list(mu = list(`s(x)` = B))
  )
  par_s <- list(mu = list(`s(x)` = c(`mu.s(x).1` = 0, `mu.s(x).2` = 0, `mu.s(x).3` = 0)))

  cache <- gamlss.longitudinal:::.gl_build_rs_design_cache(mm = mm, par_s = par_s)
  meta <- cache$mu$smooth_penalty_meta$`s(x)`

  expect_equal(ncol(cache$mu$X), 5)
  expect_equal(colnames(cache$mu$X), c("mu.(Intercept)", "mu.x", names(par_s$mu$`s(x)`)))
  expect_equal(cache$mu$fixed_names, c("mu.(Intercept)", "mu.x"))
  expect_equal(meta$idx, 3:5)
  expect_equal(unname(meta$B), unname(B))
  expect_equal(colnames(meta$B), names(par_s$mu$`s(x)`))
  expect_equal(meta$S_base, penalty)
})

test_that("RS eta calculator uses cached path when fast option is enabled", {
  captured <- NULL
  rs_design_cache <- list(mu = list(X = matrix(1, nrow = 2)))
  mm <- list(x = list(mu = matrix(1, nrow = 2)))
  eta_out_current <- list(eta_inv = list(mu = c(1, 2)))

  rs_calc_eta <- gamlss.longitudinal:::.gl_build_rs_eta_calculator(
    rs_design_cache = rs_design_cache,
    mm = mm,
    margin_dist = "NO",
    copula_link = "identity",
    option_fn = function(name, default) {
      captured$option <<- c(name = name, default = default)
      TRUE
    },
    cached_eta_fn = function(rs_design_cache, par_cov, par_s, margin_dist, copula_link, update_only, eta_out) {
      captured$cached <<- list(
        rs_design_cache = rs_design_cache,
        par_cov = par_cov,
        par_s = par_s,
        margin_dist = margin_dist,
        copula_link = copula_link,
        update_only = update_only,
        eta_out = eta_out
      )
      list(path = "cached")
    },
    full_eta_fn = function(...) {
      stop("full eta path should not be used", call. = FALSE)
    }
  )

  out <- rs_calc_eta(
    par_cov_current = c(`mu.(Intercept)` = 1),
    par_s_current = list(mu = list()),
    update_only = "mu",
    eta_out_current = eta_out_current
  )

  expect_equal(out$path, "cached")
  expect_equal(unname(captured$option["name"]), "gamlss.longitudinal.fast_rs_eta")
  expect_true(as.logical(captured$option["default"]))
  expect_identical(captured$cached$rs_design_cache, rs_design_cache)
  expect_equal(captured$cached$par_cov, c(`mu.(Intercept)` = 1))
  expect_identical(captured$cached$par_s, list(mu = list()))
  expect_identical(captured$cached$margin_dist, "NO")
  expect_identical(captured$cached$copula_link, "identity")
  expect_identical(captured$cached$update_only, "mu")
  expect_identical(captured$cached$eta_out, eta_out_current)
})

test_that("RS eta calculator falls back to full calc_eta when fast option is disabled", {
  captured <- NULL
  rs_design_cache <- list(mu = list(X = matrix(1, nrow = 2)))
  mm <- list(x = list(mu = matrix(1, nrow = 2)))

  rs_calc_eta <- gamlss.longitudinal:::.gl_build_rs_eta_calculator(
    rs_design_cache = rs_design_cache,
    mm = mm,
    margin_dist = "NO",
    copula_link = "identity",
    option_fn = function(name, default) {
      captured$option <<- c(name = name, default = default)
      FALSE
    },
    cached_eta_fn = function(...) {
      stop("cached eta path should not be used", call. = FALSE)
    },
    full_eta_fn = function(par_cov, mm, margin_dist, copula_link, par_s) {
      captured$full <<- list(
        par_cov = par_cov,
        mm = mm,
        margin_dist = margin_dist,
        copula_link = copula_link,
        par_s = par_s
      )
      list(path = "full")
    }
  )

  out <- rs_calc_eta(
    par_cov_current = c(`mu.(Intercept)` = 1),
    par_s_current = list(mu = list()),
    update_only = "mu",
    eta_out_current = list(eta_inv = list(mu = c(1, 2)))
  )

  expect_equal(out$path, "full")
  expect_equal(unname(captured$option["name"]), "gamlss.longitudinal.fast_rs_eta")
  expect_true(as.logical(captured$option["default"]))
  expect_equal(captured$full$par_cov, c(`mu.(Intercept)` = 1))
  expect_identical(captured$full$mm, mm)
  expect_identical(captured$full$margin_dist, "NO")
  expect_identical(captured$full$copula_link, "identity")
  expect_identical(captured$full$par_s, list(mu = list()))
})

test_that("RS eta length validation catches silent row dropping", {
  mm <- list(
    x = list(
      mu = matrix(1, nrow = 3),
      sigma = matrix(1, nrow = 3),
      theta = matrix(1, nrow = 2)
    )
  )
  eta_inv <- list(
    mu = c(1, 2),
    sigma = c(1, 2, 3),
    theta = c(0.1, 0.2)
  )

  expect_error(
    gamlss.longitudinal:::.gl_validate_rs_eta_lengths(
      eta_inv = eta_inv,
      mm = mm,
      response = c(1, 2, 3)
    ),
    "mu=2 vs response=3",
    fixed = TRUE
  )

  eta_inv$mu <- c(1, 2, 3)
  expect_invisible(gamlss.longitudinal:::.gl_validate_rs_eta_lengths(
    eta_inv = eta_inv,
    mm = mm,
    response = c(1, 2, 3)
  ))
})

test_that("RS score input helper coerces and validates derivative lengths", {
  out <- gamlss.longitudinal:::.gl_prepare_rs_score_inputs(
    eta = list(mu = c(0, 1, 2)),
    eta_dr = list(mu = matrix(c(1, 1, 1), ncol = 1)),
    d1 = matrix(c(0.1, 0.2, 0.3), ncol = 1),
    par_name = "mu"
  )

  expect_equal(out$eta, c(0, 1, 2))
  expect_equal(out$d1, c(0.1, 0.2, 0.3))
  expect_equal(out$eta_dr, c(1, 1, 1))

  expect_error(
    gamlss.longitudinal:::.gl_prepare_rs_score_inputs(
      eta = list(mu = c(0, 1, 2)),
      eta_dr = list(mu = c(1, 1, 1)),
      d1 = c(0.1, 0.2),
      par_name = "mu"
    ),
    "Score derivative length mismatch for mu",
    fixed = TRUE
  )
  expect_error(
    gamlss.longitudinal:::.gl_prepare_rs_score_inputs(
      eta = list(mu = c(0, 1, 2)),
      eta_dr = list(mu = c(1, 1)),
      d1 = c(0.1, 0.2, 0.3),
      par_name = "mu"
    ),
    "Link-derivative length mismatch for mu",
    fixed = TRUE
  )
})

test_that("RS timer helper appends labelled injected elapsed seconds", {
  timer <- c(`Calc Lik` = 0.1)

  out <- gamlss.longitudinal:::.gl_record_rs_timer_step(
    timer = timer,
    timer_start = as.POSIXct("2024-01-01 00:00:00", tz = "UTC"),
    label = "Backfitting",
    elapsed_sec = 2.5
  )

  expect_equal(out, c(`Calc Lik` = 0.1, Backfitting = 2.5))
})

test_that("RS timer helper computes elapsed seconds when not supplied", {
  out <- gamlss.longitudinal:::.gl_record_rs_timer_step(
    timer = numeric(),
    timer_start = as.POSIXct("2024-01-01 00:00:00", tz = "UTC"),
    label = "Copula Derivatives",
    difftime_fn = function(time1, time2, units) {
      expect_equal(units, "secs")
      7
    },
    sys_time_fn = function() as.POSIXct("2024-01-01 00:00:07", tz = "UTC")
  )

  expect_equal(out, c(`Copula Derivatives` = 7))
})

test_that("RS beta-start helper combines fixed and smooth coefficients in design order", {
  design_info <- list(
    X = matrix(
      1,
      nrow = 2,
      ncol = 4,
      dimnames = list(NULL, c("mu.(Intercept)", "mu.x", "mu.s(x).1", "mu.s(x).2"))
    ),
    fixed_names = c("mu.(Intercept)", "mu.x")
  )
  par_cov <- c(`mu.(Intercept)` = 1, mu.x = 2)

  fixed_only <- gamlss.longitudinal:::.gl_rs_beta_start(
    par_name = "mu",
    par_cov = par_cov,
    par_s = list(mu = list()),
    design_info = design_info
  )
  expect_equal(fixed_only, par_cov)

  smooth_start <- gamlss.longitudinal:::.gl_rs_beta_start(
    par_name = "mu",
    par_cov = par_cov,
    par_s = list(mu = list(`s(x)` = c(`original name 1` = 3, `original name 2` = 4))),
    design_info = design_info
  )
  expect_equal(smooth_start, c(par_cov, `mu.s(x).1` = 3, `mu.s(x).2` = 4))
})

test_that("RS backfitting input helper prepares score and design context", {
  captured <- NULL
  design_info <- list(
    X = matrix(
      1,
      nrow = 3,
      ncol = 2,
      dimnames = list(NULL, c("mu.(Intercept)", "mu.x"))
    ),
    fixed_names = c("mu.(Intercept)", "mu.x"),
    smooth_penalty_meta = list()
  )
  par_cov <- c(`mu.(Intercept)` = 1, mu.x = 2)

  out <- gamlss.longitudinal:::.gl_rs_backfitting_inputs(
    eta = list(mu = c(0.1, 0.2, 0.3)),
    eta_dr = list(mu = c(1, 1.5, 2)),
    d1 = matrix(c(0.2, -0.4, 0.6), ncol = 1),
    par_name = "mu",
    rs_design_cache = list(mu = design_info),
    par_cov = par_cov,
    par_s = list(mu = list()),
    score_fn = function(eta, dldpar, d2ldpar, dpardeta) {
      captured <<- list(
        eta = eta,
        dldpar = dldpar,
        d2ldpar = d2ldpar,
        dpardeta = dpardeta
      )
      list(w_k = matrix(c(1, 2, 3), ncol = 1), z_k = c(4, 5, 6))
    }
  )

  expect_equal(captured$eta, c(0.1, 0.2, 0.3))
  expect_equal(captured$dldpar, c(0.2, -0.4, 0.6))
  expect_equal(captured$d2ldpar, -(c(0.2, -0.4, 0.6)^2))
  expect_equal(captured$dpardeta, c(1, 1.5, 2))
  expect_equal(out$d1, c(0.2, -0.4, 0.6))
  expect_equal(out$eta_dr_vec, c(1, 1.5, 2))
  expect_identical(out$design_info, design_info)
  expect_equal(out$w_k_vec, c(1, 2, 3))
  expect_equal(out$z_k, c(4, 5, 6))
  expect_equal(out$beta_start, par_cov)
})

test_that("RS backfitting input helper includes smooth starts in design order", {
  design_info <- list(
    X = matrix(
      1,
      nrow = 2,
      ncol = 3,
      dimnames = list(NULL, c("mu.(Intercept)", "mu.s(x).1", "mu.s(x).2"))
    ),
    fixed_names = "mu.(Intercept)",
    smooth_penalty_meta = list(`s(x)` = list(idx = 2:3))
  )

  out <- gamlss.longitudinal:::.gl_rs_backfitting_inputs(
    eta = list(mu = c(0.1, 0.2)),
    eta_dr = list(mu = c(1, 1)),
    d1 = c(0.2, 0.3),
    par_name = "mu",
    rs_design_cache = list(mu = design_info),
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(mu = list(`s(x)` = c(a = 2, b = 3))),
    score_fn = function(...) list(w_k = c(1, 1), z_k = c(0, 0))
  )

  expect_equal(out$beta_start, c(`mu.(Intercept)` = 1, `mu.s(x).1` = 2, `mu.s(x).2` = 3))
})

test_that("RS likelihood context forwards model data and unpacks likelihood fields", {
  captured <- NULL
  eta_inv <- list(mu = 1:3)
  mm <- list(x = list(mu = matrix(1, nrow = 3)))
  dataset <- data.frame(response = c(1, 2, 3), time = c(1, 2, 3), subject = c("a", "a", "b"))
  pair_cache <- list(pair = 1L)
  margin_eval_cache <- new.env(parent = emptyenv())
  calc_lik <- list(
    log_lik = c(margin = -1, copula = -2, joint = -3),
    margin_d = c(0.1, 0.2, 0.3),
    margin_p = c(0.4, 0.5, 0.6),
    margin_deriv = list(dldm = matrix(1:3, ncol = 1)),
    copula_d = c(0.7, 0.8),
    copula_p = c(0.9, 1),
    Fx_1_2 = matrix(c(0.2, 0.3), ncol = 2),
    order_copula = c(2, 1)
  )

  out <- gamlss.longitudinal:::.gl_rs_likelihood_context(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    likelihood_fn = function(
      eta_inv,
      mm,
      margin_dist,
      copula_dist,
      calc_d2,
      response,
      response_margin,
      response_subject,
      pair_cache,
      margin_eval_cache
    ) {
      captured <<- list(
        eta_inv = eta_inv,
        mm = mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        calc_d2 = calc_d2,
        response = response,
        response_margin = response_margin,
        response_subject = response_subject,
        pair_cache = pair_cache,
        margin_eval_cache = margin_eval_cache
      )
      calc_lik
    }
  )

  expect_identical(captured$eta_inv, eta_inv)
  expect_identical(captured$mm, mm$x)
  expect_identical(captured$margin_dist, "NO")
  expect_identical(captured$copula_dist, "N")
  expect_false(captured$calc_d2)
  expect_equal(captured$response, dataset$response)
  expect_equal(captured$response_margin, dataset$time)
  expect_equal(captured$response_subject, dataset$subject)
  expect_identical(captured$pair_cache, pair_cache)
  expect_identical(captured$margin_eval_cache, margin_eval_cache)
  expect_identical(out$calc_lik_out, calc_lik)
  expect_identical(out$log_lik, calc_lik$log_lik)
  expect_identical(out$margin_d, calc_lik$margin_d)
  expect_identical(out$margin_p, calc_lik$margin_p)
  expect_identical(out$margin_deriv, calc_lik$margin_deriv)
  expect_identical(out$copula_d, calc_lik$copula_d)
  expect_identical(out$copula_p, calc_lik$copula_p)
  expect_identical(out$Fx_1_2, calc_lik$Fx_1_2)
  expect_identical(out$order_copula, calc_lik$order_copula)
})

test_that("RS likelihood context requests only active marginal derivatives", {
  captured <- list()
  eta_inv <- list(mu = 1:3)
  mm <- list(x = list(mu = matrix(1, nrow = 3)))
  dataset <- data.frame(response = c(1, 2, 3), time = c(1, 2, 3), subject = c("a", "a", "b"))
  pair_cache <- list(pair = 1L)
  margin_eval_cache <- new.env(parent = emptyenv())
  calc_lik <- list(
    log_lik = c(margin = -1, copula = -2, joint = -3),
    margin_d = c(0.1, 0.2, 0.3),
    margin_p = c(0.4, 0.5, 0.6),
    margin_deriv = list(dldm = matrix(1:3, ncol = 1)),
    copula_d = c(0.7, 0.8),
    copula_p = c(0.9, 1),
    Fx_1_2 = matrix(c(0.2, 0.3), ncol = 2),
    order_copula = c(2, 1)
  )
  likelihood_fn <- function(
    eta_inv,
    mm,
    margin_dist,
    copula_dist,
    calc_d2,
    response,
    response_margin,
    response_subject,
    pair_cache,
    margin_eval_cache,
    margin_deriv_names = NULL
  ) {
    captured[[length(captured) + 1L]] <<- margin_deriv_names
    calc_lik
  }

  gamlss.longitudinal:::.gl_rs_likelihood_context(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    par_name = "mu",
    likelihood_fn = likelihood_fn
  )
  gamlss.longitudinal:::.gl_rs_likelihood_context(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    par_name = "theta",
    likelihood_fn = likelihood_fn
  )

  expect_identical(captured[[1]], "dldm")
  expect_identical(captured[[2]], character(0))
})

test_that("RS likelihood context reuses current margins for copula parameters", {
  likelihood_called <- FALSE
  updater_called <- FALSE
  eta_inv <- list(theta = c(0.2, 0.3))
  mm <- list(x = list(theta = matrix(1, nrow = 2)))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  pair_cache <- list(pair = 1L)
  margin_eval_cache <- new.env(parent = emptyenv())
  current_lik <- list(
    log_lik = c(marginal = -1, copula = -2, joint = -3),
    margin_d = c(0.1, 0.2),
    margin_p = c(0.4, 0.5),
    margin_deriv = list(),
    copula_d = c(0.7),
    copula_p = c(NA_real_),
    Fx_1_2 = matrix(c(0.2, 0.3), ncol = 2),
    order_copula = 1L
  )
  updated_lik <- current_lik
  updated_lik$log_lik <- c(marginal = -1, copula = -1.5, joint = -2.5)

  out <- gamlss.longitudinal:::.gl_rs_likelihood_context(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    current_calc_lik_out = current_lik,
    par_name = "theta",
    likelihood_fn = function(...) {
      likelihood_called <<- TRUE
      current_lik
    },
    copula_likelihood_update_fn = function(eta_inv, base_lik, copula_dist, pair_cache) {
      updater_called <<- TRUE
      expect_identical(base_lik, current_lik)
      expect_identical(copula_dist, "N")
      expect_identical(pair_cache, list(pair = 1L))
      updated_lik
    }
  )

  expect_false(likelihood_called)
  expect_true(updater_called)
  expect_identical(out$calc_lik_out, updated_lik)
  expect_identical(out$log_lik, updated_lik$log_lik)
  expect_identical(out$margin_d, updated_lik$margin_d)
  expect_identical(out$margin_p, updated_lik$margin_p)
})

test_that("RS likelihood context reuses current likelihood for margin score parameters", {
  likelihood_called <- FALSE
  eta_inv <- list(mu = c(1.1, 1.2))
  mm <- list(x = list(mu = matrix(1, nrow = 2)))
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  pair_cache <- list(pair = 1L)
  margin_eval_cache <- list(
    margin_deriv_cache = list(list(
      name = "dldm",
      FUN = function(y, mu) y + mu,
      args = c("y", "mu")
    ))
  )
  current_lik <- list(
    log_lik = c(marginal = -1, copula = -2, joint = -3),
    margin_d = c(0.1, 0.2),
    margin_p = c(0.4, 0.5),
    margin_deriv = list(),
    copula_d = c(0.7),
    copula_p = c(NA_real_),
    Fx_1_2 = matrix(c(0.2, 0.3), ncol = 2),
    order_copula = 1L
  )

  out <- gamlss.longitudinal:::.gl_rs_likelihood_context(
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    current_calc_lik_out = current_lik,
    par_name = "mu",
    likelihood_fn = function(...) {
      likelihood_called <<- TRUE
      current_lik
    }
  )

  expected_lik <- current_lik
  expected_lik$margin_deriv <- list(dldm = dataset$response + eta_inv$mu)
  expect_false(likelihood_called)
  expect_identical(out$calc_lik_out, expected_lik)
  expect_identical(out$margin_deriv, expected_lik$margin_deriv)
  expect_identical(out$log_lik, current_lik$log_lik)
  expect_identical(out$margin_d, current_lik$margin_d)
  expect_identical(out$margin_p, current_lik$margin_p)
})

test_that("RS iteration likelihood state sequences eta likelihood timers and history", {
  calls <- character()
  captured <- list()
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  mm <- list(x = list(mu = matrix(1, nrow = 2)))
  par_cov <- c(`mu.(Intercept)` = 0.2)
  par_s <- list(mu = list())
  pair_cache <- list(pair = TRUE)
  margin_eval_cache <- new.env(parent = emptyenv())
  log_lik_history <- matrix(ncol = 3, nrow = 0)
  par_history <- matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)"))
  eta_out <- list(
    eta = list(mu = c(0.1, 0.2)),
    eta_dr = list(mu = c(1, 1)),
    eta_inv = list(mu = c(1.1, 1.2))
  )
  calc_lik_out <- list(
    log_lik = c(margin = -1, copula = -2, joint = -3),
    margin_d = c(0.2, 0.3),
    margin_p = c(0.4, 0.5),
    margin_deriv = list(mu = c(1, 2)),
    copula_d = c(0.6),
    copula_p = c(0.7),
    Fx_1_2 = matrix(c(0.1, 0.2), ncol = 2),
    order_copula = 1L
  )

  out <- gamlss.longitudinal:::.gl_evaluate_rs_iteration_likelihood_state(
    rs_calc_eta = function(par_cov_current, par_s_current) {
      calls <<- c(calls, "eta")
      captured$eta_par_cov <<- par_cov_current
      captured$eta_par_s <<- par_s_current
      eta_out
    },
    par_cov = par_cov,
    par_s = par_s,
    mm = mm,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = dataset,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    first_outer_run = TRUE,
    outer_start_log_lik = 0,
    timer = c(),
    timer_start = as.POSIXct("2024-01-01", tz = "UTC"),
    log_lik_history = log_lik_history,
    par_history = par_history,
    validate_eta_fn = function(eta_inv, mm, response) {
      calls <<- c(calls, "validate")
      captured$validate_eta_inv <<- eta_inv
      captured$validate_mm <<- mm
      captured$validate_response <<- response
      invisible(TRUE)
    },
    likelihood_context_fn = function(eta_inv, mm, margin_dist, copula_dist, dataset, pair_cache, margin_eval_cache) {
      calls <<- c(calls, "likelihood")
      captured$likelihood_eta_inv <<- eta_inv
      captured$likelihood_mm <<- mm
      captured$likelihood_margin_dist <<- margin_dist
      captured$likelihood_copula_dist <<- copula_dist
      captured$likelihood_dataset <<- dataset
      captured$likelihood_pair_cache <<- pair_cache
      captured$likelihood_margin_eval_cache <<- margin_eval_cache
      c(list(calc_lik_out = calc_lik_out), calc_lik_out)
    },
    outer_start_fn = function(first_outer_run, outer_start_log_lik, log_lik) {
      calls <<- c(calls, "outer_start")
      captured$outer_start_first <<- first_outer_run
      captured$outer_start_previous <<- outer_start_log_lik
      captured$outer_start_log_lik <<- log_lik
      list(first_outer_run = FALSE, outer_start_log_lik = log_lik["joint"])
    },
    timer_fn = function(timer, timer_start, label) {
      calls <<- c(calls, paste0("timer:", label))
      timer <- c(timer, length(timer) + 1)
      names(timer)[length(timer)] <- label
      timer
    },
    history_fn = function(log_lik_history, par_history, calc_lik_out, par_cov) {
      calls <<- c(calls, "history")
      captured$history_log_lik_history <<- log_lik_history
      captured$history_par_history <<- par_history
      captured$history_calc_lik_out <<- calc_lik_out
      captured$history_par_cov <<- par_cov
      list(
        log_lik_history = rbind(log_lik_history, calc_lik_out$log_lik),
        par_history = rbind(par_history, par_cov[colnames(par_history)])
      )
    },
    state_builder_fn = function(...) {
      calls <<- c(calls, "state")
      gamlss.longitudinal:::.gl_build_rs_iteration_likelihood_state(...)
    }
  )

  expect_equal(calls, c(
    "eta",
    "validate",
    "likelihood",
    "outer_start",
    "timer:Calc Lik",
    "timer:Numerical Derivatives",
    "history",
    "state"
  ))
  expect_identical(captured$eta_par_cov, par_cov)
  expect_identical(captured$eta_par_s, par_s)
  expect_identical(captured$validate_eta_inv, eta_out$eta_inv)
  expect_identical(captured$validate_mm, mm)
  expect_equal(captured$validate_response, dataset$response)
  expect_identical(captured$likelihood_eta_inv, eta_out$eta_inv)
  expect_identical(captured$likelihood_pair_cache, pair_cache)
  expect_identical(captured$likelihood_margin_eval_cache, margin_eval_cache)
  expect_true(captured$outer_start_first)
  expect_equal(captured$outer_start_previous, 0)
  expect_equal(captured$outer_start_log_lik, calc_lik_out$log_lik)
  expect_identical(captured$history_calc_lik_out, calc_lik_out)
  expect_identical(captured$history_par_cov, par_cov)
  expect_identical(out$eta_out, eta_out)
  expect_identical(out$eta, eta_out$eta)
  expect_identical(out$eta_dr, eta_out$eta_dr)
  expect_identical(out$eta_inv, eta_out$eta_inv)
  expect_identical(out$calc_lik_out, calc_lik_out)
  expect_equal(out$log_lik, calc_lik_out$log_lik)
  expect_false(out$first_outer_run)
  expect_equal(out$outer_start_log_lik, c(joint = -3))
  expect_equal(names(out$timer), c("Calc Lik", "Numerical Derivatives"))
  expect_equal(out$log_lik_history[1, ], calc_lik_out$log_lik)
  expect_equal(out$par_history[1, ], par_cov)
})

test_that("RS iteration likelihood state builder preserves returned-state contract", {
  eta_out <- list(
    eta = list(mu = c(0.1, 0.2)),
    eta_dr = list(mu = c(1, 1)),
    eta_inv = list(mu = c(1.1, 1.2))
  )
  calc_lik_out <- list(log_lik = c(margin = -1, copula = -2, joint = -3))
  likelihood_context <- c(
    list(calc_lik_out = list(calc = TRUE)),
    calc_lik_out,
    list(Fx_1_2 = matrix(c(0.1, 0.2), ncol = 2))
  )
  history_state <- list(
    log_lik_history = matrix(c(-1, -2, -3), nrow = 1),
    par_history = matrix(1, nrow = 1)
  )

  out <- gamlss.longitudinal:::.gl_build_rs_iteration_likelihood_state(
    eta_out = eta_out,
    eta = eta_out$eta,
    eta_dr = eta_out$eta_dr,
    eta_inv = eta_out$eta_inv,
    likelihood_context = likelihood_context,
    outer_start_state = list(first_outer_run = FALSE, outer_start_log_lik = c(joint = -3)),
    timer = c(`Calc Lik` = 0.1, `Numerical Derivatives` = 0.2),
    history_state = history_state
  )

  expect_identical(out$eta_out, eta_out)
  expect_identical(out$eta, eta_out$eta)
  expect_identical(out$eta_dr, eta_out$eta_dr)
  expect_identical(out$eta_inv, eta_out$eta_inv)
  expect_identical(out$calc_lik_out, list(calc = TRUE))
  expect_equal(out$log_lik, calc_lik_out$log_lik)
  expect_false(out$first_outer_run)
  expect_equal(out$outer_start_log_lik, c(joint = -3))
  expect_equal(names(out$timer), c("Calc Lik", "Numerical Derivatives"))
  expect_identical(out$log_lik_history, history_state$log_lik_history)
  expect_identical(out$par_history, history_state$par_history)
})

test_that("RS discrete score helper skips non-discrete likelihoods", {
  called <- FALSE

  out <- gamlss.longitudinal:::.gl_rs_discrete_scores(
    calc_lik_out = list(likelihood_type = "continuous"),
    eta_inv = list(mu = 1:2),
    mm = list(x = list(mu = matrix(1, nrow = 2))),
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = 1:2, time = 1:2, subject = 1L),
    pair_cache = list(),
    discrete_score_method = "analytic",
    score_fn = function(...) {
      called <<- TRUE
      list(mu = 1:2)
    }
  )

  expect_null(out)
  expect_false(called)
})

test_that("RS discrete score helper forwards rectangle likelihood context", {
  captured <- NULL
  calc_lik <- list(likelihood_type = "discrete_rectangle", log_lik = c(joint = -1))
  eta_inv <- list(mu = 1:3)
  mm <- list(x = list(mu = matrix(1, nrow = 3)))
  dataset <- data.frame(response = 1:3, time = c(1, 2, 3), subject = c(1, 1, 2))
  pair_cache <- list(row_id1 = 1L)

  out <- gamlss.longitudinal:::.gl_rs_discrete_scores(
    calc_lik_out = calc_lik,
    eta_inv = eta_inv,
    mm = mm,
    margin_dist = "DEL",
    copula_dist = "C",
    dataset = dataset,
    pair_cache = pair_cache,
    discrete_score_method = "finite",
    score_fn = function(
      eta_inv,
      mm,
      margin_dist,
      copula_dist,
      response,
      response_margin,
      response_subject,
      pair_cache,
      calc_lik,
      method
    ) {
      captured <<- list(
        eta_inv = eta_inv,
        mm = mm,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        response = response,
        response_margin = response_margin,
        response_subject = response_subject,
        pair_cache = pair_cache,
        calc_lik = calc_lik,
        method = method
      )
      list(mu = c(0.1, 0.2, 0.3), theta = c(0.4, 0.5, 0.6))
    }
  )

  expect_equal(out$mu, c(0.1, 0.2, 0.3))
  expect_identical(captured$eta_inv, eta_inv)
  expect_identical(captured$mm, mm$x)
  expect_identical(captured$margin_dist, "DEL")
  expect_identical(captured$copula_dist, "C")
  expect_equal(captured$response, dataset$response)
  expect_equal(captured$response_margin, dataset$time)
  expect_equal(captured$response_subject, dataset$subject)
  expect_identical(captured$pair_cache, pair_cache)
  expect_identical(captured$calc_lik, calc_lik)
  expect_identical(captured$method, "finite")
})

test_that("RS copula derivative context clamps probabilities and forwards likelihood metadata", {
  captured <- NULL
  calc_lik <- list(
    copula_par1 = c(0.1, 0.2),
    copula_par2 = c(3, 4),
    pair_complete = c(TRUE, FALSE)
  )

  out <- gamlss.longitudinal:::.gl_rs_copula_derivative_context(
    eta_inv = list(theta = c(0.1, 0.2)),
    Fx_1_2 = matrix(c(-0.2, 0.4, 1.2, 0.8), ncol = 2),
    copula_dist = "N",
    calc_lik_out = calc_lik,
    derivative_fn = function(eta_inv, Fx_1_2, copula_dist, par1, par2, pair_complete) {
      captured <<- list(
        eta_inv = eta_inv,
        Fx_1_2 = Fx_1_2,
        copula_dist = copula_dist,
        par1 = par1,
        par2 = par2,
        pair_complete = pair_complete
      )
      list(
        dldth = c(1, 2),
        dcdth = c(3, 4),
        dcdu1 = c(5, 6),
        dcdu2 = c(7, 8)
      )
    }
  )

  expect_equal(captured$Fx_1_2, matrix(c(0, 0.4, 1, 0.8), ncol = 2))
  expect_identical(captured$eta_inv, list(theta = c(0.1, 0.2)))
  expect_identical(captured$copula_dist, "N")
  expect_identical(captured$par1, calc_lik$copula_par1)
  expect_identical(captured$par2, calc_lik$copula_par2)
  expect_identical(captured$pair_complete, calc_lik$pair_complete)
  expect_equal(out$Fx_1_2, captured$Fx_1_2)
  expect_equal(out$dldth, c(1, 2))
  expect_equal(out$dcdth, c(3, 4))
  expect_equal(out$dcdu1, c(5, 6))
  expect_equal(out$dcdu2, c(7, 8))
  expect_null(out$dldz)
  expect_null(out$dcdz)
})

test_that("RS copula derivative context returns zeta derivatives when present", {
  out <- gamlss.longitudinal:::.gl_rs_copula_derivative_context(
    eta_inv = list(theta = 0.1, zeta = 2),
    Fx_1_2 = matrix(c(0.2, 0.4), ncol = 2),
    copula_dist = "t",
    calc_lik_out = list(copula_par1 = 0.1, copula_par2 = 2, pair_complete = TRUE),
    derivative_fn = function(...) {
      list(
        dldth = 1,
        dcdth = 2,
        dcdu1 = 3,
        dcdu2 = 4,
        dldz = 5,
        dcdz = 6
      )
    }
  )

  expect_equal(out$dldz, 5)
  expect_equal(out$dcdz, 6)
})

test_that("RS copula derivative context requests only active block derivatives", {
  requested <- list()
  calc_lik <- list(copula_par1 = c(0.1, 0.2), copula_par2 = c(3, 4), pair_complete = c(TRUE, TRUE))
  derivative_fn <- function(eta_inv, Fx_1_2, copula_dist, par1, par2, pair_complete, derivatives) {
    requested[[length(requested) + 1L]] <<- derivatives
    out <- lapply(derivatives, function(x) c(1, 2))
    names(out) <- derivatives
    out
  }

  gamlss.longitudinal:::.gl_rs_copula_derivative_context(
    eta_inv = list(mu = c(1, 2), theta = c(0.1, 0.2), zeta = c(3, 4)),
    Fx_1_2 = matrix(c(0.2, 0.3, 0.4, 0.5), ncol = 2),
    copula_dist = "t",
    calc_lik_out = calc_lik,
    par_name = "mu",
    include_dlcopdpar = TRUE,
    derivative_fn = derivative_fn
  )
  gamlss.longitudinal:::.gl_rs_copula_derivative_context(
    eta_inv = list(theta = c(0.1, 0.2), zeta = c(3, 4)),
    Fx_1_2 = matrix(c(0.2, 0.3, 0.4, 0.5), ncol = 2),
    copula_dist = "t",
    calc_lik_out = calc_lik,
    par_name = "theta",
    derivative_fn = derivative_fn
  )
  gamlss.longitudinal:::.gl_rs_copula_derivative_context(
    eta_inv = list(theta = c(0.1, 0.2), zeta = c(3, 4)),
    Fx_1_2 = matrix(c(0.2, 0.3, 0.4, 0.5), ncol = 2),
    copula_dist = "t",
    calc_lik_out = calc_lik,
    par_name = "zeta",
    derivative_fn = derivative_fn
  )

  expect_equal(requested, list(c("dcdu1", "dcdu2"), "dldth", "dldz"))
})

test_that("RS parameter score state sequences copula discrete and score helpers", {
  calls <- character()
  captured <- list()
  printed <- NULL
  dataset <- data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a"))
  mm <- list(x = list(mu = matrix(1, nrow = 2)))
  eta <- list(mu = c(0.1, 0.2))
  eta_inv <- list(mu = c(1.1, 1.2), theta = c(0.3, 0.4))
  Fx_1_2 <- matrix(c(0.2, 0.3), ncol = 2)
  calc_lik_out <- list(likelihood_type = "continuous")
  pair_cache <- list(pair = TRUE)
  copula_context <- list(
    Fx_1_2 = matrix(c(0.4, 0.5), ncol = 2),
    copula_derivatives = list(dldth = c(1, 2)),
    dldth = c(1, 2),
    dcdth = c(3, 4),
    dcdu1 = c(5, 6),
    dcdu2 = c(7, 8),
    dldz = NULL,
    dcdz = NULL
  )
  discrete_scores <- list(mu = c(0.1, 0.2))
  score_assembly <- list(
    path = "margin",
    d1 = matrix(c(9, 10), ncol = 1),
    d1_m = matrix(c(11, 12), ncol = 1),
    d1_cop = matrix(c(13, 14), ncol = 1)
  )

  out <- gamlss.longitudinal:::.gl_evaluate_rs_parameter_score_state(
    par_name = "mu",
    eta = eta,
    eta_inv = eta_inv,
    Fx_1_2 = Fx_1_2,
    copula_dist = "N",
    calc_lik_out = calc_lik_out,
    mm = mm,
    margin_dist = "NO",
    dataset = dataset,
    pair_cache = pair_cache,
    discrete_score_method = "finite",
    include_dlcopdpar = TRUE,
    margin_deriv = list(dldm = c(1, 2)),
    copula_d = c(0.7),
    log_lik = c(joint = -3),
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 2L,
    verbose = 4,
    timer = c(),
    timer_start = as.POSIXct("2024-01-01", tz = "UTC"),
    copula_context_fn = function(eta_inv, Fx_1_2, copula_dist, calc_lik_out) {
      calls <<- c(calls, "copula_context")
      captured$copula_eta_inv <<- eta_inv
      captured$copula_Fx_1_2 <<- Fx_1_2
      captured$copula_dist <<- copula_dist
      captured$copula_calc_lik_out <<- calc_lik_out
      copula_context
    },
    timer_fn = function(timer, timer_start, label) {
      calls <<- c(calls, paste0("timer:", label))
      timer <- c(timer, length(timer) + 1)
      names(timer)[length(timer)] <- label
      timer
    },
    discrete_scores_fn = function(calc_lik_out, eta_inv, mm, margin_dist, copula_dist, dataset, pair_cache, discrete_score_method) {
      calls <<- c(calls, "discrete_scores")
      captured$discrete_calc_lik_out <<- calc_lik_out
      captured$discrete_eta_inv <<- eta_inv
      captured$discrete_mm <<- mm
      captured$discrete_margin_dist <<- margin_dist
      captured$discrete_copula_dist <<- copula_dist
      captured$discrete_dataset <<- dataset
      captured$discrete_pair_cache <<- pair_cache
      captured$discrete_method <<- discrete_score_method
      discrete_scores
    },
    parameter_score_fn = function(
      par_name,
      discrete_scores,
      include_dlcopdpar,
      eta,
      eta_inv,
      response,
      margin_deriv,
      margin_dist,
      copula_dist,
      dataset,
      mm,
      calc_lik_out,
      copula_derivatives,
      dcdu1,
      dcdu2,
      copula_d,
      log_lik,
      pair_cache,
      check_dlcopdpar_gradient,
      outer_only_run_counter,
      verbose
    ) {
      calls <<- c(calls, "parameter_score")
      captured$score_par_name <<- par_name
      captured$score_discrete_scores <<- discrete_scores
      captured$score_include_dlcopdpar <<- include_dlcopdpar
      captured$score_eta <<- eta
      captured$score_eta_inv <<- eta_inv
      captured$score_response <<- response
      captured$score_margin_deriv <<- margin_deriv
      captured$score_copula_derivatives <<- copula_derivatives
      captured$score_dcdu1 <<- dcdu1
      captured$score_dcdu2 <<- dcdu2
      captured$score_copula_d <<- copula_d
      captured$score_log_lik <<- log_lik
      captured$score_pair_cache <<- pair_cache
      captured$score_gradient_check <<- check_dlcopdpar_gradient
      captured$score_outer_only_run_counter <<- outer_only_run_counter
      captured$score_verbose <<- verbose
      score_assembly
    },
    print_fn = function(x) {
      calls <<- c(calls, "print")
      printed <<- x
    },
    state_builder_fn = function(...) {
      calls <<- c(calls, "state")
      gamlss.longitudinal:::.gl_build_rs_parameter_score_state(...)
    }
  )

  expect_equal(calls, c(
    "copula_context",
    "timer:Copula Derivatives",
    "discrete_scores",
    "parameter_score",
    "timer:Margin Derivatives",
    "print",
    "state"
  ))
  expect_identical(captured$copula_eta_inv, eta_inv)
  expect_equal(captured$copula_Fx_1_2, Fx_1_2)
  expect_identical(captured$discrete_pair_cache, pair_cache)
  expect_identical(captured$discrete_method, "finite")
  expect_identical(captured$score_discrete_scores, discrete_scores)
  expect_true(captured$score_include_dlcopdpar)
  expect_equal(captured$score_response, dataset$response)
  expect_identical(captured$score_copula_derivatives, copula_context$copula_derivatives)
  expect_equal(captured$score_dcdu1, copula_context$dcdu1)
  expect_equal(captured$score_dcdu2, copula_context$dcdu2)
  expect_false(captured$score_gradient_check)
  expect_equal(captured$score_outer_only_run_counter, 2L)
  expect_equal(captured$score_verbose, 4)
  expect_identical(out$score_assembly, score_assembly)
  expect_identical(out$discrete_scores, discrete_scores)
  expect_equal(out$Fx_1_2, copula_context$Fx_1_2)
  expect_equal(out$dldth, copula_context$dldth)
  expect_equal(out$d1, score_assembly$d1)
  expect_equal(out$d1_m, score_assembly$d1_m)
  expect_equal(out$d1_cop, score_assembly$d1_cop)
  expect_equal(names(out$timer), c("Copula Derivatives", "Margin Derivatives"))
  expect_identical(printed, out$timer)
})

test_that("RS parameter score state builder preserves returned-state contract", {
  copula_context <- list(
    Fx_1_2 = matrix(c(0.4, 0.5), ncol = 2),
    copula_derivatives = list(dldth = c(1, 2)),
    dldth = c(1, 2),
    dcdth = c(3, 4),
    dcdu1 = c(5, 6),
    dcdu2 = c(7, 8),
    dldz = NULL,
    dcdz = NULL
  )
  discrete_scores <- list(mu = c(0.1, 0.2))
  score_assembly <- list(
    path = "margin",
    d1 = matrix(c(9, 10), ncol = 1),
    d1_m = matrix(c(11, 12), ncol = 1),
    d1_cop = matrix(c(13, 14), ncol = 1)
  )
  timer <- c(`Copula Derivatives` = 0.1, `Margin Derivatives` = 0.2)

  out <- gamlss.longitudinal:::.gl_build_rs_parameter_score_state(
    copula_context = copula_context,
    discrete_scores = discrete_scores,
    score_assembly = score_assembly,
    timer = timer
  )

  expect_equal(out$Fx_1_2, copula_context$Fx_1_2)
  expect_identical(out$copula_derivatives, copula_context$copula_derivatives)
  expect_equal(out$dldth, copula_context$dldth)
  expect_equal(out$dcdth, copula_context$dcdth)
  expect_equal(out$dcdu1, copula_context$dcdu1)
  expect_equal(out$dcdu2, copula_context$dcdu2)
  expect_null(out$dldz)
  expect_null(out$dcdz)
  expect_identical(out$discrete_scores, discrete_scores)
  expect_identical(out$score_assembly, score_assembly)
  expect_equal(out$d1, score_assembly$d1)
  expect_equal(out$d1_m, score_assembly$d1_m)
  expect_equal(out$d1_cop, score_assembly$d1_cop)
  expect_identical(out$timer, timer)
})

test_that("RS copula parameter score aggregates theta pair scores by response row", {
  out <- gamlss.longitudinal:::.gl_rs_copula_parameter_score(
    par_name = "theta",
    eta = list(theta = rep(0, 4)),
    response = 1:4,
    calc_lik_out = list(copula_row_id1 = c(1, 2, 2, 4)),
    copula_derivatives = list(dldth = c(1, 2, 3, 4))
  )

  expect_equal(drop(out), c(1, 5, 0, 4))
  expect_equal(colnames(out), "dldtheta")
})

test_that("RS copula parameter score aggregates zeta pair scores by copula parameter index", {
  out <- gamlss.longitudinal:::.gl_rs_copula_parameter_score(
    par_name = "zeta",
    eta = list(zeta = rep(0, 2)),
    response = 1:4,
    calc_lik_out = list(
      copula_row_id1 = c(1, 2, 3, 4),
      copula_theta_index_map = c(1, 1, 2, NA)
    ),
    copula_derivatives = list(dldz = c(1, 2, 3, 4))
  )

  expect_equal(drop(out), c(3, 3))
  expect_equal(colnames(out), "dldzeta")
})

test_that("RS copula parameter score ignores invalid mapped indices and handles no pairs", {
  out <- gamlss.longitudinal:::.gl_rs_copula_parameter_score(
    par_name = "theta",
    eta = list(theta = rep(0, 2)),
    response = 1:4,
    calc_lik_out = list(
      copula_row_id1 = c(1, 2, 3, 4),
      copula_theta_index_map = c(1, 3, 2, NA)
    ),
    copula_derivatives = list(dldth = c(1, 2, 3, 4))
  )

  expect_equal(drop(out), c(1, 3))

  empty <- gamlss.longitudinal:::.gl_rs_copula_parameter_score(
    par_name = "theta",
    eta = list(theta = rep(0, 2)),
    response = 1:4,
    calc_lik_out = list(copula_row_id1 = integer(), copula_theta_index_map = integer()),
    copula_derivatives = list(dldth = numeric())
  )
  expect_equal(drop(empty), c(0, 0))
})

test_that("RS copula parameter score preserves unexpected parameter error", {
  expect_error(
    gamlss.longitudinal:::.gl_rs_copula_parameter_score(
      par_name = "rho",
      eta = list(rho = 0),
      response = 1,
      calc_lik_out = list(copula_row_id1 = 1),
      copula_derivatives = list(dldth = 1)
    ),
    "Unexpected copula parameter in optimisation: rho",
    fixed = TRUE
  )
})

test_that("RS margin parameter score returns the base margin derivative when dlcopdpar is off", {
  out <- gamlss.longitudinal:::.gl_rs_margin_parameter_score(
    par_name = "mu",
    margin_deriv = list(dldm = matrix(c(1, 2, 3), ncol = 1)),
    include_dlcopdpar = FALSE,
    eta = list(mu = 1:3),
    eta_inv = list(mu = 1:3),
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:3),
    mm = list(x = list(mu = matrix(1, nrow = 3))),
    calc_lik_out = list(),
    dcdu1 = numeric(),
    dcdu2 = numeric(),
    copula_d = numeric(),
    log_lik = c(joint = -1),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0
  )

  expect_equal(out$d1, c(1, 2, 3))
  expect_equal(drop(out$d1_m), c(1, 2, 3))
  expect_equal(drop(out$d1_cop), c(0, 0, 0))
})

test_that("RS margin parameter score adds dlcopdpar contribution when requested", {
  fx_args <- NULL
  dl_args <- NULL

  out <- gamlss.longitudinal:::.gl_rs_margin_parameter_score(
    par_name = "sigma",
    margin_deriv = list(dldd = matrix(c(1, 2, 3), ncol = 1)),
    include_dlcopdpar = TRUE,
    eta = list(sigma = 1:3),
    eta_inv = list(sigma = 1:3),
    margin_dist = "margin",
    copula_dist = "N",
    dataset = data.frame(response = 1:3),
    mm = list(x = list(sigma = matrix(1, nrow = 3))),
    calc_lik_out = list(
      copula_row_id1 = c(1, 2),
      copula_row_id2 = c(2, 3),
      pair_complete = c(TRUE, TRUE)
    ),
    dcdu1 = c(0.1, 0.2),
    dcdu2 = c(0.3, 0.4),
    copula_d = c(2, 2),
    log_lik = c(joint = -10),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0,
    fx_deriv_fn = function(eta_inv, mm, margin_dist, response, par_names) {
      fx_args <<- list(
        eta_inv = eta_inv,
        mm = mm,
        margin_dist = margin_dist,
        response = response,
        par_names = par_names
      )
      list(sigma = c(10, 20, 30))
    },
    dlcopdpar_fn = function(row_id1, row_id2, dcdu1, dcdu2, copula_d, F_nd, n_obs, pair_complete) {
      dl_args <<- list(
        row_id1 = row_id1,
        row_id2 = row_id2,
        dcdu1 = dcdu1,
        dcdu2 = dcdu2,
        copula_d = copula_d,
        F_nd = F_nd,
        n_obs = n_obs,
        pair_complete = pair_complete
      )
      matrix(c(0.5, 1.5, 2.5), ncol = 1)
    }
  )

  expect_equal(out$d1, c(1.5, 3.5, 5.5))
  expect_equal(drop(out$d1_m), c(1, 2, 3))
  expect_equal(drop(out$d1_cop), c(0.5, 1.5, 2.5))
  expect_equal(fx_args$par_names, "sigma")
  expect_equal(dl_args$F_nd, c(10, 20, 30))
  expect_equal(dl_args$n_obs, 3L)
})

test_that("RS margin parameter score preserves dlcopdpar gradient warning hook", {
  expect_warning(
    out <- gamlss.longitudinal:::.gl_rs_margin_parameter_score(
      par_name = "mu",
      margin_deriv = list(dldm = matrix(c(1, 2), ncol = 1)),
      include_dlcopdpar = TRUE,
      eta = list(mu = 1:2),
      eta_inv = list(mu = 1:2),
      margin_dist = NULL,
      copula_dist = "N",
      dataset = data.frame(response = 1:2),
      mm = list(x = list(mu = matrix(1, nrow = 2))),
      calc_lik_out = list(copula_row_id1 = 1L, copula_row_id2 = 2L, pair_complete = TRUE),
      dcdu1 = 0.1,
      dcdu2 = 0.2,
      copula_d = 1,
      log_lik = c(joint = -3),
      pair_cache = list(),
      check_dlcopdpar_gradient = TRUE,
      outer_only_run_counter = 1L,
      verbose = 0,
      fx_deriv_fn = function(...) list(mu = c(1, 1)),
      dlcopdpar_fn = function(...) matrix(c(0, 0), ncol = 1),
      gradient_check_fn = function(...) list(warned = TRUE, message = "gradient check warning")
    ),
    "gradient check warning",
    fixed = TRUE
  )
  expect_equal(out$d1, c(1, 2))
})

test_that("RS parameter score dispatch uses discrete scores for discrete joint paths", {
  out_margin <- gamlss.longitudinal:::.gl_rs_parameter_score(
    par_name = "mu",
    discrete_scores = list(mu = c(1, 2), theta = c(3, 4)),
    include_dlcopdpar = TRUE,
    eta = list(mu = 1:2),
    eta_inv = list(mu = 1:2),
    response = 1:2,
    margin_deriv = list(),
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:2),
    mm = list(x = list(mu = matrix(1, nrow = 2))),
    calc_lik_out = list(),
    copula_derivatives = list(),
    dcdu1 = numeric(),
    dcdu2 = numeric(),
    copula_d = numeric(),
    log_lik = c(joint = -1),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0
  )
  expect_equal(drop(out_margin$d1), c(1, 2))
  expect_equal(colnames(out_margin$d1), "dldmu")
  expect_identical(out_margin$path, "discrete")

  out_copula <- gamlss.longitudinal:::.gl_rs_parameter_score(
    par_name = "theta",
    discrete_scores = list(mu = c(1, 2), theta = c(3, 4)),
    include_dlcopdpar = FALSE,
    eta = list(theta = 1:2),
    eta_inv = list(theta = 1:2),
    response = 1:2,
    margin_deriv = list(),
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:2),
    mm = list(x = list(theta = matrix(1, nrow = 2))),
    calc_lik_out = list(),
    copula_derivatives = list(),
    dcdu1 = numeric(),
    dcdu2 = numeric(),
    copula_d = numeric(),
    log_lik = c(joint = -1),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0
  )
  expect_equal(drop(out_copula$d1), c(3, 4))
  expect_equal(colnames(out_copula$d1), "dldtheta")
  expect_identical(out_copula$path, "discrete")
})

test_that("RS parameter score dispatch routes continuous margin and copula paths", {
  margin_called <- FALSE
  copula_called <- FALSE

  margin_out <- gamlss.longitudinal:::.gl_rs_parameter_score(
    par_name = "mu",
    discrete_scores = NULL,
    include_dlcopdpar = FALSE,
    eta = list(mu = 1:2),
    eta_inv = list(mu = 1:2),
    response = 1:2,
    margin_deriv = list(dldm = c(1, 2)),
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:2),
    mm = list(x = list(mu = matrix(1, nrow = 2))),
    calc_lik_out = list(),
    copula_derivatives = list(),
    dcdu1 = numeric(),
    dcdu2 = numeric(),
    copula_d = numeric(),
    log_lik = c(joint = -1),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0,
    margin_score_fn = function(...) {
      margin_called <<- TRUE
      list(d1 = matrix(c(5, 6), ncol = 1), d1_m = matrix(c(5, 6), ncol = 1), d1_cop = matrix(c(0, 0), ncol = 1))
    }
  )
  expect_true(margin_called)
  expect_equal(drop(margin_out$d1), c(5, 6))
  expect_identical(margin_out$path, "margin")

  copula_out <- gamlss.longitudinal:::.gl_rs_parameter_score(
    par_name = "theta",
    discrete_scores = NULL,
    include_dlcopdpar = FALSE,
    eta = list(theta = 1:2),
    eta_inv = list(theta = 1:2),
    response = 1:2,
    margin_deriv = list(),
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:2),
    mm = list(x = list(theta = matrix(1, nrow = 2))),
    calc_lik_out = list(copula_row_id1 = integer()),
    copula_derivatives = list(dldth = numeric()),
    dcdu1 = numeric(),
    dcdu2 = numeric(),
    copula_d = numeric(),
    log_lik = c(joint = -1),
    pair_cache = NULL,
    check_dlcopdpar_gradient = FALSE,
    outer_only_run_counter = 1L,
    verbose = 0,
    copula_score_fn = function(...) {
      copula_called <<- TRUE
      matrix(c(7, 8), ncol = 1)
    }
  )
  expect_true(copula_called)
  expect_equal(drop(copula_out$d1), c(7, 8))
  expect_identical(copula_out$path, "copula")
})

test_that("RS backfitting proposal updates fixed effects and evaluates likelihood", {
  X <- cbind(`mu.(Intercept)` = 1, `mu.x` = c(0, 1, 2))
  design_info <- list(
    X = X,
    fixed_names = colnames(X),
    smooth_penalty_meta = list()
  )
  par_cov <- c(`mu.(Intercept)` = 0, `mu.x` = 0)
  z_k <- c(1, 2, 4)
  w_k_vec <- c(1, 1, 1)
  expected_beta <- as.vector(solve(t(X) %*% X, t(X) %*% z_k))
  names(expected_beta) <- colnames(X)

  rs_calc_eta <- function(par_cov_current, par_s_current) {
    list(eta_inv = list(mu = as.numeric(X %*% par_cov_current[colnames(X)])))
  }
  likelihood_calls <- 0L
  likelihood_fn <- function(...) {
    likelihood_calls <<- likelihood_calls + 1L
    list(log_lik = c(margin = -1, copula = -2, joint = -3))
  }

  out <- gamlss.longitudinal:::.gl_rs_backfitting_iteration(
    par_s = list(mu = list()),
    par_cov = par_cov,
    beta_start = par_cov,
    lambda_s = list(mu = list()),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = z_k, time = seq_along(z_k), subject = 1L),
    mm = list(x = list(mu = X)),
    df_s = list(mu = list()),
    step_size = 1,
    par_name = "mu",
    design_info = design_info,
    w_k_vec = w_k_vec,
    z_k = z_k,
    rs_smooth_trust_radius = Inf,
    rs_calc_eta = rs_calc_eta,
    calc_lik_out = list(log_lik = c(margin = -2, copula = -2, joint = -4)),
    likelihood_fn = likelihood_fn
  )

  expect_equal(out$par_cov[colnames(X)], expected_beta, tolerance = 1e-12)
  expect_equal(unname(out$calc_lik_out_end$log_lik["joint"]), -3)
  expect_equal(unname(out$GAIC_lambda_k), 6)
  expect_equal(likelihood_calls, 1L)
})

test_that("RS backfitting proposal applies smooth trust radius and records EDF", {
  X <- cbind(`mu.(Intercept)` = 1, `mu.s(x).1` = c(-1, 0, 1))
  B <- matrix(X[, "mu.s(x).1"], ncol = 1)
  colnames(B) <- "mu.s(x).1"
  design_info <- list(
    X = X,
    fixed_names = "mu.(Intercept)",
    smooth_penalty_meta = list(`s(x)` = list(idx = 2L, B = B, S_base = matrix(1)))
  )
  beta_start <- c(`mu.(Intercept)` = 0, `mu.s(x).1` = 0)
  z_k <- c(-10, 0, 10)

  rs_calc_eta <- function(par_cov_current, par_s_current) {
    beta <- c(par_cov_current["mu.(Intercept)"], par_s_current$mu$`s(x)`)
    list(eta_inv = list(mu = as.numeric(X %*% beta)))
  }
  likelihood_fn <- function(...) {
    list(log_lik = c(margin = -1, copula = -1, joint = -2))
  }

  out <- gamlss.longitudinal:::.gl_rs_backfitting_iteration(
    par_s = list(mu = list(`s(x)` = c(`mu.s(x).1` = 0))),
    par_cov = c(`mu.(Intercept)` = 0),
    beta_start = beta_start,
    lambda_s = list(mu = list(`s(x)` = 0)),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = z_k, time = seq_along(z_k), subject = 1L),
    mm = list(x = list(mu = X)),
    df_s = list(mu = list(`s(x)` = 0)),
    step_size = 1,
    par_name = "mu",
    design_info = design_info,
    w_k_vec = c(1, 1, 1),
    z_k = z_k,
    rs_smooth_trust_radius = 0.25,
    rs_calc_eta = rs_calc_eta,
    calc_lik_out = list(log_lik = c(margin = -2, copula = -2, joint = -4)),
    likelihood_fn = likelihood_fn
  )

  expect_lte(abs(unname(out$par_s$mu$`s(x)`)), 0.25 + 1e-12)
  expect_equal(unname(out$par_s$mu$`s(x)`), 0.25, tolerance = 1e-12)
  expect_equal(out$df_s$mu$`s(x)`, 1, tolerance = 1e-12)
  expect_equal(unname(out$GAIC_lambda_k), 6)
})

test_that("RS backfitting proposal uses copula likelihood update for theta fast path", {
  old_fast_copula <- getOption("gamlss.longitudinal.fast_copula_lik")
  on.exit(options(gamlss.longitudinal.fast_copula_lik = old_fast_copula), add = TRUE)
  options(gamlss.longitudinal.fast_copula_lik = TRUE)

  X <- matrix(1, nrow = 3, ncol = 1, dimnames = list(NULL, "theta.(Intercept)"))
  design_info <- list(
    X = X,
    fixed_names = colnames(X),
    smooth_penalty_meta = list()
  )
  par_cov <- c(`theta.(Intercept)` = 0)

  rs_calc_eta <- function(par_cov_current, par_s_current) {
    list(eta_inv = list(theta = rep(unname(par_cov_current), 3)))
  }
  likelihood_fn <- function(...) {
    stop("full likelihood path should not be called", call. = FALSE)
  }
  copula_update_calls <- 0L
  copula_likelihood_update_fn <- function(...) {
    copula_update_calls <<- copula_update_calls + 1L
    list(log_lik = c(margin = -1, copula = -1, joint = -2))
  }

  out <- gamlss.longitudinal:::.gl_rs_backfitting_iteration(
    par_s = list(theta = list()),
    par_cov = par_cov,
    beta_start = par_cov,
    lambda_s = list(theta = list()),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1:3, time = 1:3, subject = 1L),
    mm = list(x = list(theta = X)),
    df_s = list(theta = list()),
    step_size = 1,
    par_name = "theta",
    design_info = design_info,
    w_k_vec = c(1, 1, 1),
    z_k = c(0.2, 0.2, 0.2),
    rs_smooth_trust_radius = Inf,
    rs_calc_eta = rs_calc_eta,
    calc_lik_out = list(log_lik = c(margin = -2, copula = -2, joint = -4)),
    likelihood_fn = likelihood_fn,
    copula_likelihood_update_fn = copula_likelihood_update_fn
  )

  expect_equal(copula_update_calls, 1L)
  expect_equal(unname(out$calc_lik_out_end$log_lik["joint"]), -2)
})

test_that("RS inner backfitting step sequences preparation runner acceptance and timers", {
  calls <- character()
  captured <- list(timer_elapsed_sec = list())
  runner_token <- function(...) list()
  timer_start <- as.POSIXct("2024-01-01", tz = "UTC")
  backfitting_inputs <- list(
    score_inputs = list(mu = TRUE),
    score = list(w_k = c(1, 2)),
    d1 = c(0.1, 0.2),
    eta_dr_vec = c(1, 1),
    design_info = list(X = matrix(1, nrow = 2, ncol = 1)),
    w_k_vec = c(3, 4),
    z_k = c(5, 6),
    beta_start = c(`mu.(Intercept)` = 0)
  )
  rs_step <- list(
    lambda_s = list(mu = list(`s(x)` = 2)),
    par_cov = c(`mu.(Intercept)` = 1),
    par_s = list(mu = list()),
    calc_lik_out_end = list(log_lik = c(joint = -1)),
    df_s = list(mu = list(`s(x)` = 1)),
    rs_block_trace = list(data.frame(parameter = "mu")),
    change_log_lik = 0.25,
    run_counter = 3L,
    outer_run_counter = 4L,
    inner_run_counter = 5L,
    timer_label = "Plotting",
    elapsed_sec = 9
  )

  out <- gamlss.longitudinal:::.gl_run_rs_inner_backfitting_step(
    eta = list(mu = c(0.1, 0.2)),
    eta_dr = list(mu = c(1, 1)),
    d1 = c(0.3, 0.4),
    par_name = "mu",
    rs_design_cache = list(mu = list(X = matrix(1, nrow = 2))),
    par_cov = c(`mu.(Intercept)` = 0),
    par_s = list(mu = list()),
    timer = c(),
    timer_start = timer_start,
    rs_smooth_trust_radius = 0.75,
    rs_calc_eta = function(...) list(eta_inv = list(mu = 1)),
    calc_lik_out = list(log_lik = c(joint = -2)),
    pair_cache = list(pair = TRUE),
    margin_eval_cache = list(margin = TRUE),
    lambda_penalty_K = 2,
    lambda_s = list(mu = list()),
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a")),
    mm = list(x = list(mu = matrix(1, nrow = 2))),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 1L,
    outer_run_counter = 2L,
    outer_only_run_counter = 3L,
    verbose = 2,
    rs_block_trace = list(),
    run_counter = 1L,
    use_backtracking = TRUE,
    backtracking_max_halves = 4L,
    plot_results = FALSE,
    log_lik_history = matrix(ncol = 3, nrow = 0),
    par_history = matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "mu.(Intercept)")),
    true_val = NULL,
    backfitting_inputs_fn = function(eta, eta_dr, d1, par_name, rs_design_cache, par_cov, par_s) {
      calls <<- c(calls, "inputs")
      captured$inputs_d1 <<- d1
      captured$inputs_par_name <<- par_name
      captured$inputs_cache <<- rs_design_cache
      captured$inputs_par_cov <<- par_cov
      captured$inputs_par_s <<- par_s
      backfitting_inputs
    },
    timer_fn = function(timer, timer_start, label, elapsed_sec = NULL) {
      calls <<- c(calls, paste0("timer:", label))
      timer_elapsed_sec <- captured$timer_elapsed_sec
      timer_elapsed_sec[[label]] <- elapsed_sec
      captured$timer_elapsed_sec <<- timer_elapsed_sec
      timer <- c(timer, if (is.null(elapsed_sec)) length(timer) + 1 else elapsed_sec)
      names(timer)[length(timer)] <- label
      timer
    },
    runner_fn = function(design_info, w_k_vec, z_k, rs_smooth_trust_radius,
                         rs_calc_eta, calc_lik_out, pair_cache, margin_eval_cache) {
      calls <<- c(calls, "runner")
      captured$runner_design_info <<- design_info
      captured$runner_w_k_vec <<- w_k_vec
      captured$runner_z_k <<- z_k
      captured$runner_trust_radius <<- rs_smooth_trust_radius
      captured$runner_calc_lik_out <<- calc_lik_out
      captured$runner_pair_cache <<- pair_cache
      captured$runner_margin_eval_cache <<- margin_eval_cache
      runner_token
    },
    acceptance_fn = function(lambda_s, par_name, par_s, par_cov, beta_start, K,
                             margin_dist, copula_dist, dataset, mm, copula_link,
                             df_s, step_size, rs_update_lambda,
                             inner_run_counter, outer_run_counter,
                             outer_only_run_counter, verbose, backfitting_fn,
                             calc_lik_out, rs_block_trace, run_counter,
                             timer_start, use_backtracking,
                             backtracking_max_halves, plot_results,
                             log_lik_history, par_history, true_val) {
      calls <<- c(calls, "acceptance")
      captured$acceptance_lambda_s <<- lambda_s
      captured$acceptance_par_name <<- par_name
      captured$acceptance_par_s <<- par_s
      captured$acceptance_par_cov <<- par_cov
      captured$acceptance_beta_start <<- beta_start
      captured$acceptance_K <<- K
      captured$acceptance_margin_dist <<- margin_dist
      captured$acceptance_copula_dist <<- copula_dist
      captured$acceptance_step_size <<- step_size
      captured$acceptance_rs_update_lambda <<- rs_update_lambda
      captured$acceptance_inner_run_counter <<- inner_run_counter
      captured$acceptance_outer_run_counter <<- outer_run_counter
      captured$acceptance_outer_only_run_counter <<- outer_only_run_counter
      captured$acceptance_verbose <<- verbose
      captured$acceptance_backfitting_fn <<- backfitting_fn
      captured$acceptance_backtracking <<- use_backtracking
      captured$acceptance_max_halves <<- backtracking_max_halves
      captured$acceptance_plot_results <<- plot_results
      captured$acceptance_true_val <<- true_val
      rs_step
    }
  )

  expect_equal(calls, c("inputs", "timer:Backfitting", "runner", "acceptance", "timer:Plotting"))
  expect_equal(captured$inputs_d1, c(0.3, 0.4))
  expect_identical(captured$inputs_par_name, "mu")
  expect_equal(captured$runner_w_k_vec, backfitting_inputs$w_k_vec)
  expect_equal(captured$runner_z_k, backfitting_inputs$z_k)
  expect_equal(captured$runner_trust_radius, 0.75)
  expect_identical(captured$runner_pair_cache, list(pair = TRUE))
  expect_equal(captured$acceptance_beta_start, backfitting_inputs$beta_start)
  expect_equal(captured$acceptance_K, 2)
  expect_identical(captured$acceptance_backfitting_fn, runner_token)
  expect_true(captured$acceptance_rs_update_lambda)
  expect_true(captured$acceptance_backtracking)
  expect_equal(captured$acceptance_max_halves, 4L)
  expect_false(captured$acceptance_plot_results)
  expect_equal(captured$timer_elapsed_sec$Plotting, 9)
  expect_identical(out$backfitting_inputs, backfitting_inputs)
  expect_identical(out$score, backfitting_inputs$score)
  expect_equal(out$d1, backfitting_inputs$d1)
  expect_equal(out$eta_dr_vec, backfitting_inputs$eta_dr_vec)
  expect_identical(out$lambda_s, rs_step$lambda_s)
  expect_equal(out$par_cov, rs_step$par_cov)
  expect_identical(out$calc_lik_out_end, rs_step$calc_lik_out_end)
  expect_equal(out$change_log_lik, rs_step$change_log_lik)
  expect_equal(out$run_counter, rs_step$run_counter)
  expect_equal(out$outer_run_counter, rs_step$outer_run_counter)
  expect_equal(out$inner_run_counter, rs_step$inner_run_counter)
  expect_identical(out$rs_step, rs_step)
  expect_equal(names(out$timer), c("Backfitting", "Plotting"))
})

test_that("RS backfitting runner binds derivative state and forwards iteration arguments", {
  design_info <- list(X = matrix(1, nrow = 2, ncol = 1))
  w_k_vec <- c(1, 2)
  z_k <- c(3, 4)
  calc_lik_out <- list(log_lik = c(joint = -1))
  pair_cache <- list(pair = TRUE)
  margin_eval_cache <- list(margin = TRUE)

  runner <- gamlss.longitudinal:::.gl_build_rs_backfitting_runner(
    design_info = design_info,
    w_k_vec = w_k_vec,
    z_k = z_k,
    rs_smooth_trust_radius = 0.25,
    rs_calc_eta = function(...) list(eta_inv = list(mu = 1)),
    calc_lik_out = calc_lik_out,
    pair_cache = pair_cache,
    margin_eval_cache = margin_eval_cache,
    backfitting_iteration_fn = function(par_s, par_cov, beta_start, lambda_s,
                                        K, margin_dist, copula_dist, dataset,
                                        mm, df_s, step_size, par_name,
                                        design_info, w_k_vec, z_k,
                                        rs_smooth_trust_radius, rs_calc_eta,
                                        calc_lik_out, pair_cache,
                                        margin_eval_cache) {
      expect_equal(par_s, list(mu = list()))
      expect_equal(par_cov, c(`mu.(Intercept)` = 1))
      expect_equal(beta_start, c(`mu.(Intercept)` = 1))
      expect_equal(lambda_s, list(mu = list()))
      expect_equal(K, 2)
      expect_equal(copula_dist, "N")
      expect_equal(dataset$response, c(1, 2))
      expect_equal(mm$x$mu, matrix(1, nrow = 2, ncol = 1))
      expect_equal(df_s, list(mu = list()))
      expect_equal(step_size, 0.5)
      expect_equal(par_name, "mu")
      expect_equal(design_info, list(X = matrix(1, nrow = 2, ncol = 1)))
      expect_equal(w_k_vec, c(1, 2))
      expect_equal(z_k, c(3, 4))
      expect_equal(rs_smooth_trust_radius, 0.25)
      expect_type(rs_calc_eta, "closure")
      expect_equal(calc_lik_out, list(log_lik = c(joint = -1)))
      expect_equal(pair_cache, list(pair = TRUE))
      expect_equal(margin_eval_cache, list(margin = TRUE))
      list(
        par_cov = par_cov,
        par_s = par_s,
        calc_lik_out_end = calc_lik_out,
        GAIC_lambda_k = 1,
        df_s = df_s
      )
    }
  )

  out <- runner(
    par_s = list(mu = list()),
    par_cov = c(`mu.(Intercept)` = 1),
    beta_start = c(`mu.(Intercept)` = 1),
    lambda_s = list(mu = list()),
    first_inner_run = FALSE,
    K = 2,
    margin_dist = gamlss.dist::NO(),
    copula_dist = "N",
    dataset = data.frame(response = c(1, 2), time = c(1, 2), subject = c("a", "a")),
    mm = list(x = list(mu = matrix(1, nrow = 2, ncol = 1))),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    par_name = "mu"
  )

  expect_equal(out$GAIC_lambda_k, 1)
  expect_equal(out$calc_lik_out_end, calc_lik_out)
})

test_that("RS acceptance helper keeps improving proposal without backtracking", {
  make_result <- function(loglik, id) {
    list(
      par_cov = c(value = id),
      par_s = list(mu = list()),
      calc_lik_out_end = list(log_lik = c(margin = loglik - 1, copula = -1, joint = loglik)),
      GAIC_lambda_k = NA_real_,
      df_s = list(mu = list())
    )
  }

  proposal_calls <- 0L
  out <- gamlss.longitudinal:::.gl_rs_accept_backfitting_step(
    proposed_results = make_result(11, 2),
    current_results = make_result(10, 1),
    nominal_step_size = 0.5,
    use_backtracking = TRUE,
    backtracking_max_halves = 3,
    proposal_fn = function(step_size) {
      proposal_calls <<- proposal_calls + 1L
      make_result(99, step_size)
    },
    outer_iteration = 4,
    inner_iteration = 2,
    global_inner_iteration = 9,
    parameter = "mu",
    elapsed_sec = 1.25
  )

  expect_equal(out$accepted_results$par_cov[["value"]], 2)
  expect_equal(out$accepted_step_size, 0.5)
  expect_false(out$rejected)
  expect_false(out$backtracking_applied)
  expect_equal(out$backtracking_attempts, 0L)
  expect_equal(out$max_backtracking_attempts, 0L)
  expect_equal(proposal_calls, 0L)
  expect_equal(out$trace_row$accepted_logLik, 11)
  expect_equal(out$trace_row$accepted_change, 1)
  expect_false(out$trace_row$rejected)
})

test_that("RS acceptance helper accepts first non-downhill backtracking proposal", {
  make_result <- function(loglik, id) {
    list(
      par_cov = c(value = id),
      par_s = list(mu = list()),
      calc_lik_out_end = list(log_lik = c(margin = loglik - 1, copula = -1, joint = loglik)),
      GAIC_lambda_k = NA_real_,
      df_s = list(mu = list())
    )
  }

  trial_loglik <- c(9.5, 10.2)
  trial_steps <- numeric(0)
  out <- gamlss.longitudinal:::.gl_rs_accept_backfitting_step(
    proposed_results = make_result(9, 2),
    current_results = make_result(10, 1),
    nominal_step_size = 0.8,
    use_backtracking = TRUE,
    backtracking_max_halves = 3,
    proposal_fn = function(step_size) {
      trial_steps <<- c(trial_steps, step_size)
      make_result(trial_loglik[length(trial_steps)], step_size)
    },
    outer_iteration = 4,
    inner_iteration = 2,
    global_inner_iteration = 9,
    parameter = "theta",
    elapsed_sec = 1.25
  )

  expect_equal(trial_steps, c(0.4, 0.2))
  expect_equal(unname(out$accepted_results$par_cov[["value"]]), 0.2)
  expect_equal(out$accepted_step_size, 0.2)
  expect_false(out$rejected)
  expect_true(out$backtracking_applied)
  expect_equal(out$backtracking_attempts, 2L)
  expect_equal(out$max_backtracking_attempts, 3L)
  expect_equal(out$trace_row$proposed_logLik, 9)
  expect_equal(out$trace_row$accepted_logLik, 10.2)
  expect_equal(out$trace_row$accepted_step_size, 0.2)
  expect_false(out$trace_row$rejected)
})

test_that("RS acceptance helper rejects downhill proposal after exhausted backtracking", {
  make_result <- function(loglik, id) {
    list(
      par_cov = c(value = id),
      par_s = list(mu = list()),
      calc_lik_out_end = list(log_lik = c(margin = loglik - 1, copula = -1, joint = loglik)),
      GAIC_lambda_k = NA_real_,
      df_s = list(mu = list())
    )
  }

  proposal_calls <- 0L
  out <- gamlss.longitudinal:::.gl_rs_accept_backfitting_step(
    proposed_results = make_result(8, 2),
    current_results = make_result(10, 1),
    nominal_step_size = 1,
    use_backtracking = TRUE,
    backtracking_max_halves = 2,
    proposal_fn = function(step_size) {
      proposal_calls <<- proposal_calls + 1L
      make_result(9, step_size)
    },
    outer_iteration = 4,
    inner_iteration = 2,
    global_inner_iteration = 9,
    parameter = "sigma",
    elapsed_sec = 1.25
  )

  expect_equal(proposal_calls, 2L)
  expect_equal(out$accepted_results$par_cov[["value"]], 1)
  expect_equal(out$accepted_step_size, 0)
  expect_true(out$rejected)
  expect_true(out$backtracking_applied)
  expect_equal(out$backtracking_attempts, 2L)
  expect_equal(out$max_backtracking_attempts, 2L)
  expect_equal(out$trace_row$accepted_logLik, 10)
  expect_equal(out$trace_row$accepted_change, 0)
  expect_true(out$trace_row$rejected)
})

test_that("RS acceptance reporting respects verbosity thresholds", {
  cat_calls <- character()
  print_calls <- list()
  rs_acceptance <- list(
    backtracking_applied = TRUE,
    accepted_step_size = 0.25,
    backtracking_attempts = 2L,
    max_backtracking_attempts = 3L,
    start_joint_loglik = -10,
    proposed_joint_loglik = -11,
    accepted_joint_loglik = -9,
    rejected = FALSE
  )
  calc_lik_out_end <- list(log_lik = c(margin = -3, copula = -2, joint = -5))

  out <- gamlss.longitudinal:::.gl_report_rs_acceptance(
    par_name = "theta",
    step_size = 1,
    use_backtracking = TRUE,
    rs_acceptance = rs_acceptance,
    calc_lik_out_end = calc_lik_out_end,
    verbose = 1,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )

  expect_identical(out, TRUE)
  expect_equal(cat_calls, character())
  expect_equal(print_calls, list())

  gamlss.longitudinal:::.gl_report_rs_acceptance(
    par_name = "mu",
    step_size = 1,
    use_backtracking = TRUE,
    rs_acceptance = rs_acceptance,
    calc_lik_out_end = calc_lik_out_end,
    verbose = 2,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )

  expect_equal(length(cat_calls), 1)
  expect_match(cat_calls[[1]], "Backtracking applied for mu", fixed = TRUE)
  expect_match(cat_calls[[1]], "step_size 1 -> 0.25", fixed = TRUE)
  expect_equal(print_calls, list())
})

test_that("RS acceptance reporting includes theta diagnostics and log-likelihood at high verbosity", {
  cat_calls <- character()
  print_calls <- list()
  rs_acceptance <- list(
    backtracking_applied = TRUE,
    accepted_step_size = 0,
    backtracking_attempts = 2L,
    max_backtracking_attempts = 2L,
    start_joint_loglik = -10,
    proposed_joint_loglik = -11,
    accepted_joint_loglik = -10,
    rejected = TRUE
  )
  calc_lik_out_end <- list(log_lik = c(margin = -3, copula = -2, joint = -5))

  gamlss.longitudinal:::.gl_report_rs_acceptance(
    par_name = "theta",
    step_size = 1,
    use_backtracking = FALSE,
    rs_acceptance = rs_acceptance,
    calc_lik_out_end = calc_lik_out_end,
    verbose = 3,
    cat_fn = function(...) cat_calls <<- c(cat_calls, paste0(...)),
    print_fn = function(x, ...) print_calls[[length(print_calls) + 1L]] <<- x
  )

  expect_equal(length(cat_calls), 3)
  expect_match(cat_calls[[1]], "Backtracking applied for theta", fixed = TRUE)
  expect_match(cat_calls[[2]], "Theta step diagnostics: start=-10", fixed = TRUE)
  expect_match(cat_calls[[2]], "backtracking=off", fixed = TRUE)
  expect_match(cat_calls[[2]], "rejected=yes", fixed = TRUE)
  expect_equal(cat_calls[[3]], "\nLogLik:\n")
  expect_equal(print_calls[[1]], calc_lik_out_end$log_lik)
})

test_that("RS acceptance state helper applies accepted results and increments counters", {
  existing_trace <- list(data.frame(parameter = "mu", accepted_logLik = -4))
  new_trace <- data.frame(parameter = "sigma", accepted_logLik = -2)
  accepted_calc_lik <- list(log_lik = c(margin = -1, copula = -1, joint = -2))
  base_calc_lik <- list(log_lik = c(margin = -2, copula = -2, joint = -4))
  accepted_results <- list(
    par_cov = c(mu = 1, sigma = 2),
    par_s = list(mu = list(), sigma = list()),
    calc_lik_out_end = accepted_calc_lik,
    df_s = list(mu = list(), sigma = list())
  )

  out <- gamlss.longitudinal:::.gl_apply_rs_acceptance_state(
    rs_acceptance = list(
      accepted_results = accepted_results,
      trace_row = new_trace
    ),
    rs_block_trace = existing_trace,
    calc_lik_out = base_calc_lik,
    run_counter = 3,
    outer_run_counter = 7,
    inner_run_counter = 2
  )

  expect_identical(out$par_cov, accepted_results$par_cov)
  expect_identical(out$par_s, accepted_results$par_s)
  expect_identical(out$calc_lik_out_end, accepted_calc_lik)
  expect_identical(out$df_s, accepted_results$df_s)
  expect_equal(length(out$rs_block_trace), 2)
  expect_identical(out$rs_block_trace[[1]], existing_trace[[1]])
  expect_identical(out$rs_block_trace[[2]], new_trace)
  expect_equal(out$change_log_lik, c(joint = 2))
  expect_equal(out$run_counter, 4)
  expect_equal(out$outer_run_counter, 8)
  expect_equal(out$inner_run_counter, 3)
})

test_that("RS backfitting acceptance runner sequences proposal acceptance reporting and plotting", {
  call_order <- character()
  lambda_initial <- list(mu = list(`s(x)` = 1))
  lambda_updated <- list(mu = list(`s(x)` = 2))
  accepted_calc_lik <- list(log_lik = c(margin = -1, copula = -1, joint = -2))
  accepted_results <- list(
    par_cov = c(mu = 3),
    par_s = list(mu = list(`s(x)` = 4)),
    calc_lik_out_end = accepted_calc_lik,
    df_s = list(mu = list(`s(x)` = 1.5))
  )
  trace_row <- data.frame(parameter = "mu", accepted_logLik = -2)

  out <- gamlss.longitudinal:::.gl_run_rs_backfitting_acceptance(
    lambda_s = lambda_initial,
    par_name = "mu",
    par_s = list(mu = list(`s(x)` = 0)),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(x = list(mu = matrix(1, nrow = 1))),
    copula_link = "identity",
    df_s = list(mu = list(`s(x)` = 1)),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 2L,
    outer_only_run_counter = 3L,
    outer_run_counter = 4L,
    verbose = 2,
    backfitting_fn = function(...) list(),
    calc_lik_out = list(log_lik = c(margin = -3, copula = -2, joint = -5)),
    rs_block_trace = list(data.frame(parameter = "old")),
    run_counter = 7L,
    timer_start = as.POSIXct("2024-01-01 00:00:00", tz = "UTC"),
    use_backtracking = TRUE,
    backtracking_max_halves = 3L,
    plot_results = TRUE,
    log_lik_history = matrix(1, nrow = 1, ncol = 3),
    par_history = matrix(1, nrow = 1, ncol = 1),
    true_val = NA_real_,
    prepare_proposal_fn = function(lambda_s, par_name, par_s, par_cov,
                                   beta_start, K, margin_dist, copula_dist,
                                   dataset, mm, copula_link, df_s, step_size,
                                   rs_update_lambda, inner_run_counter,
                                   outer_only_run_counter, verbose,
                                   backfitting_fn) {
      call_order <<- c(call_order, "prepare")
      expect_equal(lambda_s, lambda_initial)
      expect_equal(par_name, "mu")
      expect_equal(inner_run_counter, 2L)
      expect_equal(outer_only_run_counter, 3L)
      expect_true(rs_update_lambda)
      list(lambda_s = lambda_updated, proposed_results = list(id = "proposal"))
    },
    accept_step_fn = function(proposed_results, current_results, nominal_step_size,
                              use_backtracking, backtracking_max_halves,
                              proposal_fn, outer_iteration, inner_iteration,
                              global_inner_iteration, parameter, elapsed_sec) {
      call_order <<- c(call_order, "accept")
      expect_equal(proposed_results, list(id = "proposal"))
      expect_equal(current_results$par_cov, c(mu = 0))
      expect_equal(current_results$calc_lik_out_end$log_lik["joint"], c(joint = -5))
      expect_equal(nominal_step_size, 0.5)
      expect_true(use_backtracking)
      expect_equal(backtracking_max_halves, 3L)
      expect_equal(outer_iteration, 3L)
      expect_equal(inner_iteration, 2L)
      expect_equal(global_inner_iteration, 4L)
      expect_equal(parameter, "mu")
      expect_equal(elapsed_sec, 10)
      trial <- proposal_fn(0.25)
      expect_equal(trial, list())
      list(accepted_results = accepted_results, trace_row = trace_row)
    },
    apply_acceptance_fn = function(rs_acceptance, rs_block_trace, calc_lik_out,
                                   run_counter, outer_run_counter,
                                   inner_run_counter) {
      call_order <<- c(call_order, "apply")
      expect_equal(length(rs_block_trace), 1L)
      expect_equal(run_counter, 7L)
      expect_equal(outer_run_counter, 4L)
      expect_equal(inner_run_counter, 2L)
      list(
        par_cov = rs_acceptance$accepted_results$par_cov,
        par_s = rs_acceptance$accepted_results$par_s,
        calc_lik_out_end = rs_acceptance$accepted_results$calc_lik_out_end,
        df_s = rs_acceptance$accepted_results$df_s,
        rs_block_trace = c(rs_block_trace, list(rs_acceptance$trace_row)),
        change_log_lik = c(joint = 3),
        run_counter = 8L,
        outer_run_counter = 5L,
        inner_run_counter = 3L
      )
    },
    report_acceptance_fn = function(par_name, step_size, use_backtracking,
                                    rs_acceptance, calc_lik_out_end, verbose) {
      call_order <<- c(call_order, "report")
      expect_equal(par_name, "mu")
      expect_equal(step_size, 0.5)
      expect_true(use_backtracking)
      expect_equal(calc_lik_out_end, accepted_calc_lik)
      expect_equal(verbose, 2)
      invisible(TRUE)
    },
    plot_progress_fn = function(log_lik_history, par_history, par_count, true_val) {
      call_order <<- c(call_order, "plot")
      expect_equal(par_count, 1L)
      invisible(TRUE)
    },
    difftime_fn = function(time1, time2, units) 10,
    sys_time_fn = function() as.POSIXct("2024-01-01 00:00:10", tz = "UTC")
  )

  expect_equal(call_order, c("prepare", "accept", "apply", "report", "plot"))
  expect_equal(out$lambda_s, lambda_updated)
  expect_equal(out$par_cov, c(mu = 3))
  expect_equal(out$calc_lik_out_end, accepted_calc_lik)
  expect_equal(out$change_log_lik, c(joint = 3))
  expect_equal(out$run_counter, 8L)
  expect_equal(out$outer_run_counter, 5L)
  expect_equal(out$inner_run_counter, 3L)
  expect_equal(out$timer_label, "Plotting")
  expect_equal(out$elapsed_sec, 10)
})

test_that("RS proposal helper evaluates fixed-only proposal once", {
  calls <- list()
  lambda_s <- list(mu = list())

  out <- gamlss.longitudinal:::.gl_rs_prepare_backfitting_proposal(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list()),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 1L,
    outer_only_run_counter = 1L,
    verbose = 0,
    backfitting_fn = function(
      par_s,
      par_cov,
      beta_start,
      lambda_s,
      first_inner_run,
      K,
      margin_dist,
      copula_dist,
      dataset,
      mm,
      copula_link,
      df_s,
      step_size,
      par_name
    ) {
      calls[[length(calls) + 1L]] <<- list(
        first_inner_run = first_inner_run,
        lambda_s = lambda_s,
        step_size = step_size,
        par_name = par_name
      )
      list(id = length(calls), first_inner_run = first_inner_run)
    },
    lambda_update_fn = function(...) {
      stop("lambda update should not be called for fixed-only terms")
    }
  )

  expect_equal(length(calls), 1)
  expect_false(calls[[1]]$first_inner_run)
  expect_identical(out$lambda_s, lambda_s)
  expect_equal(out$num_smooths, 0)
  expect_null(out$initial_results)
  expect_equal(out$proposed_results$id, 1)
})

test_that("RS proposal helper updates lambda before final smooth proposal after first outer run", {
  calls <- list()
  lambda_update_calls <- list()
  lambda_s <- list(mu = list(`s(x)` = 2))
  lambda_s_updated <- list(mu = list(`s(x)` = 5))

  out <- gamlss.longitudinal:::.gl_rs_prepare_backfitting_proposal(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list(`s(x)` = c(a = 0))),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list(`s(x)` = 1)),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 1L,
    outer_only_run_counter = 2L,
    verbose = 0,
    backfitting_fn = function(
      par_s,
      par_cov,
      beta_start,
      lambda_s,
      first_inner_run,
      K,
      margin_dist,
      copula_dist,
      dataset,
      mm,
      copula_link,
      df_s,
      step_size,
      par_name
    ) {
      calls[[length(calls) + 1L]] <<- list(
        first_inner_run = first_inner_run,
        lambda_s = lambda_s
      )
      list(id = length(calls), lambda_s = lambda_s, first_inner_run = first_inner_run)
    },
    lambda_update_fn = function(
      lambda_s,
      par_name,
      par_s,
      par_cov,
      beta_start,
      K,
      margin_dist,
      copula_dist,
      dataset,
      mm,
      copula_link,
      df_s,
      step_size,
      rs_update_lambda,
      inner_run_counter,
      verbose,
      backfitting_fn
    ) {
      lambda_update_calls[[length(lambda_update_calls) + 1L]] <<- list(
        lambda_s = lambda_s,
        par_name = par_name,
        rs_update_lambda = rs_update_lambda,
        inner_run_counter = inner_run_counter
      )
      lambda_s_updated
    }
  )

  expect_equal(length(lambda_update_calls), 1)
  expect_identical(lambda_update_calls[[1]]$lambda_s, lambda_s)
  expect_identical(lambda_update_calls[[1]]$par_name, "mu")
  expect_true(lambda_update_calls[[1]]$rs_update_lambda)
  expect_equal(lambda_update_calls[[1]]$inner_run_counter, 1L)
  expect_equal(length(calls), 1)
  expect_false(calls[[1]]$first_inner_run)
  expect_identical(calls[[1]]$lambda_s, lambda_s_updated)
  expect_identical(out$lambda_s, lambda_s_updated)
  expect_null(out$initial_results)
  expect_equal(out$num_smooths, 1)
  expect_equal(out$proposed_results$id, 1)
})

test_that("RS proposal helper evaluates first-outer smooth proposal once", {
  calls <- list()
  lambda_s <- list(mu = list(`s(x)` = 2))

  out <- gamlss.longitudinal:::.gl_rs_prepare_backfitting_proposal(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list(`s(x)` = c(a = 0))),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = "NO",
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list(`s(x)` = 1)),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 1L,
    outer_only_run_counter = 1L,
    verbose = 0,
    backfitting_fn = function(
      par_s,
      par_cov,
      beta_start,
      lambda_s,
      first_inner_run,
      K,
      margin_dist,
      copula_dist,
      dataset,
      mm,
      copula_link,
      df_s,
      step_size,
      par_name
    ) {
      calls[[length(calls) + 1L]] <<- first_inner_run
      list(id = length(calls), first_inner_run = first_inner_run)
    },
    lambda_update_fn = function(...) {
      stop("lambda update should not be called on first outer run")
    }
  )

  expect_equal(calls, list(FALSE))
  expect_identical(out$lambda_s, lambda_s)
  expect_equal(out$num_smooths, 1)
  expect_null(out$initial_results)
  expect_equal(out$proposed_results$id, 1)
})

test_that("RS progress plotting draws likelihood and parameter traces", {
  par_calls <- list()
  plot_calls <- list()
  abline_calls <- list()

  log_lik_history <- matrix(
    c(-3, -2, -5, -2, -1, -3),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("margin", "copula", "joint"))
  )
  par_history <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu", "sigma"))
  )

  out <- gamlss.longitudinal:::.gl_plot_rs_progress(
    log_lik_history = log_lik_history,
    par_history = par_history,
    par_count = 2L,
    true_val = NA_real_,
    par_fn = function(...) par_calls[[length(par_calls) + 1L]] <<- list(...),
    plot_fn = function(x, ...) plot_calls[[length(plot_calls) + 1L]] <<- c(list(x = x), list(...)),
    abline_fn = function(...) abline_calls[[length(abline_calls) + 1L]] <<- list(...)
  )

  expect_identical(out, TRUE)
  expect_equal(par_calls[[1]]$mfrow, c(3, 2))
  expect_equal(length(plot_calls), 5)
  expect_equal(plot_calls[[1]]$main, "LogLik - Overall")
  expect_equal(plot_calls[[2]]$main, "LogLik - Margin")
  expect_equal(plot_calls[[3]]$main, "LogLik - Copula")
  expect_equal(plot_calls[[4]]$main, "mu")
  expect_equal(plot_calls[[5]]$main, "sigma")
  expect_equal(length(abline_calls), 0)
})

test_that("RS progress plotting adds true-value reference lines when supplied", {
  plot_calls <- list()
  abline_calls <- list()

  log_lik_history <- matrix(1:6, nrow = 2, dimnames = list(NULL, c("margin", "copula", "joint")))
  par_history <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu", "sigma"))
  )

  gamlss.longitudinal:::.gl_plot_rs_progress(
    log_lik_history = log_lik_history,
    par_history = par_history,
    par_count = 2L,
    true_val = c(0.5, 0.6),
    par_fn = function(...) NULL,
    plot_fn = function(x, ...) plot_calls[[length(plot_calls) + 1L]] <<- c(list(x = x), list(...)),
    abline_fn = function(...) abline_calls[[length(abline_calls) + 1L]] <<- list(...)
  )

  expect_equal(length(abline_calls), 2)
  expect_equal(abline_calls[[1]], list(h = 0.5, col = "red"))
  expect_equal(abline_calls[[2]], list(h = 0.6, col = "red"))
  expect_equal(plot_calls[[4]]$ylim, range(c(par_history[, 1], 0.5)))
  expect_equal(plot_calls[[5]]$ylim, range(c(par_history[, 2], 0.6)))
})

test_that("RS lambda helper updates active smooths with optimizer results", {
  lambda_s <- list(mu = list(`s(x)` = 2, `s(z)` = 4))
  optim_calls <- list()
  objective_values <- numeric(0)
  backfitting_lambdas <- list()

  backfitting_fn <- function(par_s, par_cov, beta_start, lambda_s, first_inner_run,
                             K, margin_dist, copula_dist, dataset, mm, copula_link,
                             df_s, step_size, par_name) {
    backfitting_lambdas[[length(backfitting_lambdas) + 1L]] <<- lambda_s[[par_name]]
    active_lambda <- tail(unlist(lambda_s[[par_name]], use.names = FALSE), 1)
    list(
      calc_lik_out_end = list(log_lik = c(joint = -active_lambda)),
      df_s = df_s,
      GAIC_lambda_k = active_lambda + K
    )
  }
  optim_fn <- function(par, fn, method, lower, upper, control) {
    optim_calls[[length(optim_calls) + 1L]] <<- list(
      par = par,
      method = method,
      lower = lower,
      upper = upper,
      control = control
    )
    objective_values <<- c(objective_values, fn(par + 0.5))
    list(par = par + 1)
  }

  output <- capture.output({
    out <- gamlss.longitudinal:::.gl_rs_update_smoothing_parameters(
      lambda_s = lambda_s,
      par_name = "mu",
      par_s = list(mu = list()),
      par_cov = c(mu = 0),
      beta_start = c(mu = 0),
      K = 2,
      margin_dist = NULL,
      copula_dist = "N",
      dataset = data.frame(response = 1),
      mm = list(),
      copula_link = "identity",
      df_s = list(mu = list()),
      step_size = 0.5,
      rs_update_lambda = TRUE,
      inner_run_counter = 1L,
      verbose = 0,
      backfitting_fn = backfitting_fn,
      optim_fn = optim_fn
    )
  })

  expect_equal(out$mu$`s(x)`, 3)
  expect_equal(out$mu$`s(z)`, 5)
  expect_equal(length(optim_calls), 2L)
  expect_equal(vapply(optim_calls, `[[`, numeric(1), "par"), c(2, 4))
  expect_true(all(vapply(optim_calls, `[[`, character(1), "method") == "L-BFGS-B"))
  expect_true(all(vapply(optim_calls, `[[`, numeric(1), "lower") == 0.01))
  expect_true(all(vapply(optim_calls, `[[`, numeric(1), "upper") == 1e6))
  expect_equal(objective_values, c(6, 6.5))
  expect_equal(length(backfitting_lambdas), 2L)
  expect_equal(backfitting_lambdas[[1]]$`s(x)`, 2.5)
  expect_equal(backfitting_lambdas[[2]]$`s(z)`, 4.5)
  expect_true(any(grepl("Optimising smoothing parameter for mu - s\\(x\\)", output)))
  expect_true(any(grepl("Optimising smoothing parameter for mu - s\\(z\\)", output)))
})

test_that("RS lambda helper skips updates when disabled or not first inner iteration", {
  lambda_s <- list(mu = list(`s(x)` = 2))
  optim_calls <- 0L
  optim_fn <- function(...) {
    optim_calls <<- optim_calls + 1L
    list(par = 99)
  }
  backfitting_fn <- function(...) {
    stop("backfitting objective should not be called", call. = FALSE)
  }

  disabled <- gamlss.longitudinal:::.gl_rs_update_smoothing_parameters(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list()),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    rs_update_lambda = FALSE,
    inner_run_counter = 1L,
    verbose = 0,
    backfitting_fn = backfitting_fn,
    optim_fn = optim_fn
  )

  later_inner <- gamlss.longitudinal:::.gl_rs_update_smoothing_parameters(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list()),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 2L,
    verbose = 0,
    backfitting_fn = backfitting_fn,
    optim_fn = optim_fn
  )

  expect_equal(disabled, lambda_s)
  expect_equal(later_inner, lambda_s)
  expect_equal(optim_calls, 0L)
})

test_that("RS lambda helper handles parameters without smooths", {
  lambda_s <- list(mu = list())
  optim_calls <- 0L
  out <- gamlss.longitudinal:::.gl_rs_update_smoothing_parameters(
    lambda_s = lambda_s,
    par_name = "mu",
    par_s = list(mu = list()),
    par_cov = c(mu = 0),
    beta_start = c(mu = 0),
    K = 2,
    margin_dist = NULL,
    copula_dist = "N",
    dataset = data.frame(response = 1),
    mm = list(),
    copula_link = "identity",
    df_s = list(mu = list()),
    step_size = 0.5,
    rs_update_lambda = TRUE,
    inner_run_counter = 1L,
    verbose = 0,
    backfitting_fn = function(...) stop("should not be called", call. = FALSE),
    optim_fn = function(...) {
      optim_calls <<- optim_calls + 1L
      list(par = 99)
    }
  )

  expect_equal(out, lambda_s)
  expect_equal(optim_calls, 0L)
})
