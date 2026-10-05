#!/usr/bin/env Rscript

# Bounded, reproducible diagnostics for ACTG zero-inflated fitting. This script
# never runs an unrestricted optimizer and is intentionally outside _targets.R.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", script_argument[[1L]]), winslash = "/", mustWork = TRUE
)
root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."),
                      winslash = "/", mustWork = TRUE)
setwd(root)

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[1L]], fixed = TRUE)
}

mode <- match.arg(argument_value("mode", "bounded"), c("bounded", "primitives"))
output <- argument_value(
  "output",
  file.path(dirname(script_path), paste0("benchmark-", mode, "-latest.csv"))
)

if (!requireNamespace("devtools", quietly = TRUE)) stop("Install devtools.")
suppressPackageStartupMessages(devtools::load_all(root, quiet = TRUE))
load(file.path(root, "data", "actg_cd4.rda"))
data <- actg_cd4
data$week_factor <- factor(data$week)

git_sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
session <- list(
  git_sha = git_sha[[1L]],
  r_version = paste(R.version$major, R.version$minor, sep = "."),
  gamlss_version = as.character(utils::packageVersion("gamlss")),
  gamlss_dist_version = as.character(utils::packageVersion("gamlss.dist")),
  rows = nrow(data),
  observed = sum(is.finite(data$cd4)),
  subjects = length(unique(data$subject))
)

elapsed_call <- function(fn) {
  gc(FALSE)
  value <- NULL
  elapsed <- system.time(value <- tryCatch(fn(), error = identity))[["elapsed"]]
  list(elapsed = unname(elapsed), value = value)
}

family_object <- function(name) {
  do.call(get(name, envir = asNamespace("gamlss.dist")), list())
}

models <- data.frame(
  margin = c("ZISICHEL", "ZISICHEL", "ZIBNB", "ZIBNB", "BCPE"),
  copula = c("G", "t", "G", "t", "t"),
  response = c(rep("cd4", 4L), "cd4_plus_1"),
  stringsAsFactors = FALSE
)

run_bounded <- function(spec) {
  fit_data <- data
  if (identical(spec$response, "cd4_plus_1")) {
    fit_data$cd4 <- ifelse(is.na(fit_data$cd4), NA_real_, fit_data$cd4 + 1)
  }
  control <- gamlss_longitudinal_control(
    outer_tol = 1e-8,
    max_outer_iter = 1L,
    max_elapsed_sec = 180,
    stop_on_convergence = FALSE,
    rs = list(
      max_inner_iter = 1L,
      max_steps = 1L,
      warm_start_joint = FALSE,
      use_backtracking = FALSE,
      update_lambda = FALSE,
      discrete_score_method = "analytical"
    )
  )
  warnings <- character()
  timed <- elapsed_call(function() withCallingHandlers(
    gamlss_longitudinal(
      data = fit_data,
      time_var = "week_factor",
      subject_var = "subject",
      margin_dist = family_object(spec$margin),
      copula_dist = spec$copula,
      mu.formula = cd4 ~ 1,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      tau.formula = ~ 1,
      theta.formula = ~ 1,
      zeta.formula = ~ 1,
      missingness = "segment",
      method = "RS",
      optimization = "joint",
      optimizer_control = control,
      compute_vcov = FALSE,
      verbose = 0
    ),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  ))
  fit <- timed$value
  data.frame(
    margin = spec$margin,
    copula = spec$copula,
    response = spec$response,
    operation = "bounded_one_rs_sweep",
    elapsed_sec = timed$elapsed,
    returned_fit = inherits(fit, "gamlss.longitudinal"),
    logLik = if (inherits(fit, "gamlss.longitudinal")) as.numeric(logLik(fit)) else NA_real_,
    warning_count = length(unique(warnings)),
    error = if (inherits(fit, "error")) conditionMessage(fit) else "",
    stringsAsFactors = FALSE
  )
}

make_eta <- function(family, n) {
  common <- list(mu = rep(25, n), sigma = rep(0.5, n))
  margin <- switch(
    family,
    ZISICHEL = list(nu = rep(-0.5, n), tau = rep(0.04, n)),
    ZIBNB = list(nu = rep(2, n), tau = rep(0.04, n)),
    BCPE = list(nu = rep(1, n), tau = rep(2, n))
  )
  c(common, margin, list(theta = rep(0.35, n), zeta = rep(6, n)))
}

run_primitives <- function(spec) {
  response <- data$cd4
  if (identical(spec$response, "cd4_plus_1")) {
    response <- ifelse(is.na(response), NA_real_, response + 1)
  }
  n <- length(response)
  margin <- family_object(spec$margin)
  eta <- make_eta(spec$margin, n)
  mm <- setNames(
    lapply(c("mu", "sigma", "nu", "tau", "theta", "zeta"),
           function(unused) matrix(1, n, 1L)),
    c("mu", "sigma", "nu", "tau", "theta", "zeta")
  )
  pair_cache <- gamlss.longitudinal:::build_copula_pair_cache(
    response, data$week_factor, data$subject
  )

  base_args <- list(
    eta_inv = eta,
    mm = mm,
    margin_dist = margin,
    copula_dist = spec$copula,
    calc_d2 = FALSE,
    response = response,
    response_margin = data$week_factor,
    response_subject = data$subject,
    pair_cache = pair_cache
  )
  no_deriv <- elapsed_call(function() do.call(
    gamlss.longitudinal:::calc_likelihood_minimal,
    c(base_args, list(calc_margin_deriv = FALSE))
  ))
  with_deriv <- elapsed_call(function() do.call(
    gamlss.longitudinal:::calc_likelihood_minimal,
    c(base_args, list(calc_margin_deriv = TRUE))
  ))

  rows <- list(data.frame(
    margin = spec$margin,
    copula = spec$copula,
    response = spec$response,
    operation = "likelihood_no_margin_derivatives",
    elapsed_sec = no_deriv$elapsed,
    returned_fit = !inherits(no_deriv$value, "error"),
    logLik = if (!inherits(no_deriv$value, "error"))
      unname(no_deriv$value$log_lik[["joint"]]) else NA_real_,
    warning_count = NA_integer_,
    error = if (inherits(no_deriv$value, "error")) conditionMessage(no_deriv$value) else "",
    stringsAsFactors = FALSE
  ))

  rows[[2L]] <- data.frame(
    margin = spec$margin,
    copula = spec$copula,
    response = spec$response,
    operation = "likelihood_with_margin_derivatives",
    elapsed_sec = with_deriv$elapsed,
    returned_fit = !inherits(with_deriv$value, "error"),
    logLik = if (!inherits(with_deriv$value, "error"))
      unname(with_deriv$value$log_lik[["joint"]]) else NA_real_,
    warning_count = NA_integer_,
    error = if (inherits(with_deriv$value, "error")) conditionMessage(with_deriv$value) else "",
    stringsAsFactors = FALSE
  )

  if (!identical(spec$margin, "BCPE") && !inherits(with_deriv$value, "error")) {
    scores <- elapsed_call(function() {
      gamlss.longitudinal:::.calc_discrete_rectangle_scores(
        eta, mm, margin, spec$copula, response, data$week_factor,
        data$subject, pair_cache = pair_cache, calc_lik = with_deriv$value,
        method = "analytical"
      )
    })
    rows[[3L]] <- data.frame(
      margin = spec$margin,
      copula = spec$copula,
      response = spec$response,
      operation = "all_discrete_scores",
      elapsed_sec = scores$elapsed,
      returned_fit = !inherits(scores$value, "error"),
      logLik = NA_real_,
      warning_count = NA_integer_,
      error = if (inherits(scores$value, "error")) conditionMessage(scores$value) else "",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

runner <- if (identical(mode, "bounded")) run_bounded else run_primitives
results <- do.call(rbind, lapply(seq_len(nrow(models)), function(index) {
  runner(as.list(models[index, , drop = FALSE]))
}))
for (name in names(session)) results[[name]] <- session[[name]]
results$completed_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(results, output, row.names = FALSE)
print(results)
message("Wrote: ", normalizePath(output, winslash = "/", mustWork = TRUE))
