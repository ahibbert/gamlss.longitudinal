#!/usr/bin/env Rscript

# First-pass simulation comparison for discrete GAMLSS margins in
# gamlss.longitudinal. The script intentionally targets the currently supported
# integer-valued families BI, PO, and NBI with a Gaussian pair-copula.

safe_source <- function(path) {
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(txt) > 0) {
    bom <- intToUtf8(65279)
    txt[1] <- sub(paste0("^", bom), "", txt[1], useBytes = FALSE)
  }
  eval(parse(text = txt), envir = .GlobalEnv)
}

safe_source("R/common_functions.R")
safe_source("R/link_functions.R")

suppressPackageStartupMessages({
  library(VineCopula)
  library(gamlss)
  library(gamlss.dist)
})

set.seed(20260513)

out_dir <- file.path("results", "discrete_gaussian_copula_comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n <- as.integer(Sys.getenv("DISCRETE_SIM_N", unset = "300"))
d <- as.integer(Sys.getenv("DISCRETE_SIM_D", unset = "4"))
n_rep <- as.integer(Sys.getenv("DISCRETE_SIM_REP", unset = "1"))
max_outer_iter <- as.integer(Sys.getenv("DISCRETE_MAX_OUTER", unset = "80"))
max_inner_iter <- as.integer(Sys.getenv("DISCRETE_MAX_INNER", unset = "20"))
compute_vcov <- as.logical(as.integer(Sys.getenv("DISCRETE_COMPUTE_VCOV", unset = "0")))
verbose <- as.integer(Sys.getenv("DISCRETE_VERBOSE", unset = "0"))

methods <- c("rectangle", "midpoint", "distributional_transform")
families <- c("PO", "NBI", "BI")
rho_true <- 0.45
copula_family <- as.numeric(VineCopula::BiCopName("N"))

clip_u <- function(u) pmin(pmax(u, 1e-8), 1 - 1e-8)

family_object <- function(family_name) {
  switch(
    family_name,
    PO = gamlss.dist::PO(mu.link = "log"),
    NBI = gamlss.dist::NBI(mu.link = "log", sigma.link = "log"),
    BI = gamlss.dist::BI(mu.link = "logit"),
    stop("Unsupported family: ", family_name)
  )
}

q_family <- function(family_name, u, mu, sigma = NULL) {
  u <- clip_u(u)
  switch(
    family_name,
    PO = gamlss.dist::qPO(u, mu = mu),
    NBI = gamlss.dist::qNBI(u, mu = mu, sigma = sigma),
    BI = gamlss.dist::qBI(u, mu = mu),
    stop("Unsupported family: ", family_name)
  )
}

true_eta <- function(family_name, x, tr) {
  if (family_name == "BI") {
    eta_mu <- -0.35 + 0.70 * x + 0.45 * tr
    return(list(eta_mu = eta_mu, eta_sigma = NULL))
  }
  eta_mu <- 1.00 + 0.35 * x + 0.25 * tr
  eta_sigma <- if (family_name == "NBI") log(0.70) + 0.15 * x else NULL
  list(eta_mu = eta_mu, eta_sigma = eta_sigma)
}

simulate_discrete_dataset <- function(family_name, seed) {
  set.seed(seed)
  subjects <- data.frame(
    id = seq_len(n),
    x = rnorm(n),
    stringsAsFactors = FALSE
  )

  dat <- merge(subjects, data.frame(time_index = seq_len(d)), by = NULL, all = TRUE)
  dat <- dat[order(dat$id, dat$time_index), ]
  dat$t <- if (d > 1) (dat$time_index - 1) / (d - 1) else 0
  dat$time <- dat$time_index

  eta <- true_eta(family_name, dat$x, dat$t)
  fam <- family_object(family_name)
  mu <- fam$mu.linkinv(eta$eta_mu)
  sigma <- if (!is.null(eta$eta_sigma)) fam$sigma.linkinv(eta$eta_sigma) else NULL

  U <- matrix(NA_real_, nrow = n, ncol = d)
  for (i in seq_len(n)) {
    u_prev <- runif(1)
    U[i, 1] <- u_prev
    if (d > 1) {
      for (j in 2:d) {
        u_new <- VineCopula::BiCopCondSim(
          N = 1,
          cond.val = u_prev,
          cond.var = 1,
          family = copula_family,
          par = rho_true,
          par2 = 0
        )
        u_prev <- clip_u(as.numeric(u_new))
        U[i, j] <- u_prev
      }
    }
  }

  dat$u <- as.vector(t(U))
  dat$response <- q_family(family_name, dat$u, mu = mu, sigma = sigma)
  dat$family <- family_name
  dat[, c("id", "time", "response", "x", "t", "family")]
}

coef_or_na <- function(x, name) {
  if (is.null(x) || !name %in% names(x)) return(NA_real_)
  as.numeric(x[[name]])
}

coef_first_available <- function(x, names_try) {
  for (nm in names_try) {
    val <- coef_or_na(x, nm)
    if (is.finite(val)) return(val)
  }
  NA_real_
}

fit_gamlss_baseline <- function(dat, family_name) {
  fam <- family_object(family_name)
  sigma_formula <- if (family_name == "NBI") ~ x else ~ 1
  fit <- tryCatch(
    gamlss(
      response ~ x + t,
      sigma.formula = sigma_formula,
      family = fam,
      data = dat,
      trace = FALSE
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      model = "gamlss",
      method = NA_character_,
      converged = FALSE,
      logLik = NA_real_,
      mu_intercept = NA_real_,
      mu_x = NA_real_,
      mu_t = NA_real_,
      sigma_intercept = NA_real_,
      sigma_x = NA_real_,
      theta = NA_real_,
      message = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  mu_coef <- coef(fit, what = "mu")
  sigma_coef <- tryCatch(coef(fit, what = "sigma"), error = function(e) NULL)
  data.frame(
    model = "gamlss",
    method = NA_character_,
    converged = TRUE,
    logLik = as.numeric(logLik(fit)),
    mu_intercept = coef_or_na(mu_coef, "(Intercept)"),
    mu_x = coef_or_na(mu_coef, "x"),
    mu_t = coef_or_na(mu_coef, "t"),
    sigma_intercept = coef_or_na(sigma_coef, "(Intercept)"),
    sigma_x = coef_or_na(sigma_coef, "x"),
    theta = NA_real_,
    message = "",
    stringsAsFactors = FALSE
  )
}

fit_longitudinal <- function(dat, family_name, method_name) {
  fam <- family_object(family_name)
  sigma_formula <- if (family_name == "NBI") ~ x else ~ 1
  fit <- tryCatch(
    {
      fit_expr <- quote(gamlss.longitudinal(
        dataset = dat,
        margin_dist = fam,
        copula_dist = "N",
        time_var = "time",
        subject_var = "id",
        mu.formula = response ~ x + t,
        sigma.formula = sigma_formula,
        theta.formula = ~ 1,
        discrete_copula_method = method_name,
        method = "RS",
        include_dlcopdpar = FALSE,
        compute_vcov = compute_vcov,
        vcov_method = "numderiv",
        max_outer_iter = max_outer_iter,
        max_inner_iter = max_inner_iter,
        verbose = verbose,
        plot_results = FALSE
      ))
      if (verbose <= 0) {
        out <- NULL
        capture.output(out <- eval(fit_expr))
        out
      } else {
        eval(fit_expr)
      }
    },
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(data.frame(
      model = "gamlss.longitudinal",
      method = method_name,
      converged = FALSE,
      logLik = NA_real_,
      mu_intercept = NA_real_,
      mu_x = NA_real_,
      mu_t = NA_real_,
      sigma_intercept = NA_real_,
      sigma_x = NA_real_,
      theta = NA_real_,
      message = conditionMessage(fit),
      stringsAsFactors = FALSE
    ))
  }

  p <- fit$par
  data.frame(
    model = "gamlss.longitudinal",
    method = method_name,
    converged = TRUE,
    logLik = as.numeric(fit$calc_lik_out_end$log_lik["joint"]),
    mu_intercept = coef_first_available(p, c("mu.(Intercept)", "mu.intercept")),
    mu_x = coef_or_na(p, "mu.x"),
    mu_t = coef_or_na(p, "mu.t"),
    sigma_intercept = coef_first_available(p, c("sigma.(Intercept)", "sigma.intercept")),
    sigma_x = coef_or_na(p, "sigma.x"),
    theta = coef_first_available(p, c("theta.(Intercept)", "theta.intercept")),
    message = "",
    stringsAsFactors = FALSE
  )
}

truth_row <- function(family_name) {
  if (family_name == "BI") {
    return(data.frame(
      model = "truth",
      method = NA_character_,
      converged = TRUE,
      logLik = NA_real_,
      mu_intercept = -0.35,
      mu_x = 0.70,
      mu_t = 0.45,
      sigma_intercept = NA_real_,
      sigma_x = NA_real_,
      theta = rho_true,
      message = "",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    model = "truth",
    method = NA_character_,
    converged = TRUE,
    logLik = NA_real_,
    mu_intercept = 1.00,
    mu_x = 0.35,
    mu_t = 0.25,
    sigma_intercept = if (family_name == "NBI") log(0.70) else NA_real_,
    sigma_x = if (family_name == "NBI") 0.15 else NA_real_,
    theta = rho_true,
    message = "",
    stringsAsFactors = FALSE
  )
}

all_results <- list()
counter <- 1L
partial_csv <- file.path(out_dir, "discrete_gaussian_copula_comparison_results_partial.csv")

for (rep_id in seq_len(n_rep)) {
  for (family_name in families) {
    cat("\nSimulating ", family_name, " replicate ", rep_id, "\n", sep = "")
    dat <- simulate_discrete_dataset(family_name, seed = 20260513 + 1000 * rep_id + match(family_name, families))

    write.csv(
      dat,
      file.path(out_dir, paste0("sim_data_", tolower(family_name), "_rep", rep_id, ".csv")),
      row.names = FALSE
    )

    rows <- list(truth_row(family_name), fit_gamlss_baseline(dat, family_name))
    for (method_name in methods) {
      cat("  Fitting ", method_name, "\n", sep = "")
      rows[[length(rows) + 1L]] <- fit_longitudinal(dat, family_name, method_name)
    }

    res <- do.call(rbind, rows)
    res$family <- family_name
    res$replicate <- rep_id
    res$n <- n
    res$d <- d
    all_results[[counter]] <- res
    counter <- counter + 1L
    partial_results <- do.call(rbind, all_results)
    partial_results <- partial_results[, c("family", "replicate", "n", "d", names(partial_results)[!names(partial_results) %in% c("family", "replicate", "n", "d")])]
    write.csv(partial_results, partial_csv, row.names = FALSE)
  }
}

results <- do.call(rbind, all_results)
results <- results[, c("family", "replicate", "n", "d", names(results)[!names(results) %in% c("family", "replicate", "n", "d")])]

out_csv <- file.path(out_dir, "discrete_gaussian_copula_comparison_results.csv")
write.csv(results, out_csv, row.names = FALSE)

cat("\nResults written to: ", out_csv, "\n", sep = "")
print(results)
