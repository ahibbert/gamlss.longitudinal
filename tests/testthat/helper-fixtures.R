make_fixture_factor_time_interaction <- function(n_subject = 24L) {
  set.seed(42)

  subject_tbl <- data.frame(
    id = seq_len(n_subject),
    gender = factor(sample(c("F", "M"), n_subject, replace = TRUE)),
    age = round(runif(n_subject, min = 20, max = 70), 1),
    stringsAsFactors = FALSE
  )

  time_levels <- c("t1", "t2", "t3")
  grid <- expand.grid(
    id = subject_tbl$id,
    time_raw = factor(time_levels, levels = time_levels, ordered = TRUE),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  dat <- merge(grid, subject_tbl, by = "id", sort = FALSE)
  dat <- dat[order(dat$id, dat$time_raw), ]
  rownames(dat) <- NULL

  t_num <- as.integer(dat$time_raw)
  g_num <- ifelse(dat$gender == "M", 1, 0)

  # Include a time-gender interaction signal plus small noise.
  dat$y <- 1.5 + 0.25 * t_num + 0.4 * g_num + 0.3 * t_num * g_num + 0.01 * dat$age + rnorm(nrow(dat), sd = 0.15)

  dat
}

make_fixture_numeric_time <- function(n_subject = 24L) {
  dat <- make_fixture_factor_time_interaction(n_subject = n_subject)
  dat$time_raw <- as.numeric(dat$time_raw)
  dat
}

make_fixture_with_duplicate_subject_time <- function() {
  dat <- make_fixture_factor_time_interaction(n_subject = 12L)
  dup <- dat[1, , drop = FALSE]
  rbind(dat, dup)
}

make_fixture_missingness_margin_all_missing <- function() {
  dat <- make_fixture_factor_time_interaction(n_subject = 18L)
  t_last <- levels(dat$time_raw)[length(levels(dat$time_raw))]
  dat$y[dat$time_raw == t_last] <- NA_real_
  dat
}

make_fixture_missingness_zero_complete_pairs <- function() {
  dat <- make_fixture_factor_time_interaction(n_subject = 18L)
  t_levels <- levels(dat$time_raw)
  t1 <- t_levels[1]
  t2 <- t_levels[2]
  # Force no complete pairs for the first consecutive pair by splitting missingness.
  ids <- sort(unique(dat$id))
  half_ids <- ids[seq_len(floor(length(ids) / 2))]
  dat$y[dat$time_raw == t1 & dat$id %in% half_ids] <- NA_real_
  dat$y[dat$time_raw == t2 & !(dat$id %in% half_ids)] <- NA_real_
  dat
}

make_fixture_with_structural_missing_rows <- function() {
  dat <- make_fixture_factor_time_interaction(n_subject = 16L)
  # Remove a small number of subject-time combinations to test grid expansion
  # while keeping fitting numerically stable.
  drop_idx <- dat$id %in% c(4L, 12L) & dat$time_raw == levels(dat$time_raw)[2]
  dat[!drop_idx, , drop = FALSE]
}

make_fixture_single_level_factor <- function(n_subject = 12L) {
  dat <- make_fixture_factor_time_interaction(n_subject = n_subject)
  dat$gender <- factor("F", levels = c("F"))
  dat
}

fit_fixture_model <- function(
  dataset,
  include_dlcopdpar = TRUE,
  mu_formula = "y ~ time_raw * gender + age",
  sigma_formula = "~ time_raw + gender",
  nu_formula = "~ 1",
  tau_formula = "~ 1",
  theta_formula = "~ time_raw",
  zeta_formula = "~ 1",
  time_var = "time_raw",
  subject_var = "id",
  start_from = NA,
  max_outer_iter = 2,
  max_inner_iter = 2,
  outer_stop_crit = 1,
  inner_stop_crit = 1,
  method = "RS",
  missingness = "error",
  muffle_segment_warning = TRUE,
  use_backtracking = TRUE,
  compute_vcov = FALSE,
  verbose = 0
) {
  suppressPackageStartupMessages({
    library(gamlss)
    library(gamlss.dist)
  })

  withCallingHandlers(
    gamlss.longitudinal::gamlss_longitudinal(
      dataset = dataset,
      margin_dist = gamlss.dist::NO(),
      copula_dist = "N",
      time_var = time_var,
      subject_var = subject_var,
      missingness = missingness,
      mu.formula = mu_formula,
      sigma.formula = sigma_formula,
      nu.formula = nu_formula,
      tau.formula = tau_formula,
      theta.formula = theta_formula,
      zeta.formula = zeta_formula,
      start_from = start_from,
      include_dlcopdpar = include_dlcopdpar,
      method = method,
      optimizer_control = if (identical(toupper(method), "RS")) {
        gamlss.longitudinal::gamlss_longitudinal_control(
          outer_tol = outer_stop_crit,
          max_outer_iter = max_outer_iter,
          rs = list(
            inner_tol = inner_stop_crit,
            max_inner_iter = max_inner_iter,
            use_backtracking = use_backtracking
          )
        )
      } else {
        gamlss.longitudinal::gamlss_longitudinal_control(
          outer_tol = outer_stop_crit,
          max_outer_iter = max_outer_iter,
          cg = list(use_backtracking = use_backtracking)
        )
      },
      compute_vcov = compute_vcov,
      verbose = verbose
    ),
    warning = function(w) {
      if (
        inherits(w, "gamlss.longitudinal_nonconvergence_warning") ||
          (isTRUE(muffle_segment_warning) && inherits(w, "gamlss.longitudinal_segment_warning")) ||
          grepl("Model stopped at max_outer_iter", conditionMessage(w), fixed = TRUE)
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

get_valid_inference_fixture <- local({
  fixture <- NULL

  function() {
    if (is.null(fixture)) {
      dat <- make_fixture_factor_time_interaction(n_subject = 50L)
      fit <- fit_fixture_model(
        dat,
        include_dlcopdpar = FALSE,
        mu_formula = "y ~ time_raw + gender + age",
        sigma_formula = "~ 1",
        theta_formula = "~ 1",
        max_outer_iter = 25L,
        max_inner_iter = 25L,
        outer_stop_crit = 1e-5,
        inner_stop_crit = 1e-5
      )
      vcov_out <- stats::vcov(fit, method = "analytical")

      stopifnot(
        isTRUE(fit$convergence$converged),
        identical(vcov_out$hessian_diagnostics$status, "available")
      )
      fixture <- list(data = dat, fit = fit, vcov = vcov_out)
    }

    fixture
  }
})

get_valid_plot_fixture <- local({
  fixtures <- list()

  function(ordered = TRUE) {
    key <- if (isTRUE(ordered)) "ordered" else "unordered"
    if (is.null(fixtures[[key]])) {
      dat <- make_fixture_factor_time_interaction(n_subject = 50L)
      dat$time_raw <- factor(
        as.character(dat$time_raw),
        levels = levels(dat$time_raw),
        ordered = ordered
      )
      fit <- fit_fixture_model(
        dat,
        include_dlcopdpar = TRUE,
        max_outer_iter = 25L,
        max_inner_iter = 25L,
        outer_stop_crit = 1e-5,
        inner_stop_crit = 1e-5
      )
      vcov_out <- stats::vcov(fit, method = "analytical")

      stopifnot(
        isTRUE(fit$convergence$converged),
        identical(vcov_out$hessian_diagnostics$status, "available")
      )
      fixtures[[key]] <- list(data = dat, fit = fit, vcov = vcov_out)
    }

    fixtures[[key]]
  }
})

get_valid_smooth_plot_fixture <- local({
  fixture <- NULL

  function() {
    if (is.null(fixture)) {
      dat <- make_fixture_factor_time_interaction(n_subject = 50L)
      fit <- fit_fixture_model(
        dat,
        include_dlcopdpar = TRUE,
        mu_formula = "y ~ time_raw * gender + s(log(age), bs = 'ps')",
        max_outer_iter = 25L,
        max_inner_iter = 25L,
        outer_stop_crit = 1e-5,
        inner_stop_crit = 1e-5
      )
      vcov_out <- stats::vcov(fit, method = "analytical")

      stopifnot(
        isTRUE(fit$convergence$converged),
        identical(vcov_out$hessian_diagnostics$status, "available")
      )
      fixture <- list(data = dat, fit = fit, vcov = vcov_out)
    }

    fixture
  }
})

capture_warnings <- function(expr) {
  warnings <- character(0)
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}
