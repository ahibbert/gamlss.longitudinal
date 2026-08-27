# Sanitized private-analysis contract. This file contains no LIPID records,
# pseudo-observations, estimates, fitted objects, aggregates, or output assets.
jss_lipid_contract <- function() {
  list(
    required_columns = c("subject", "time", "response", "treatment"),
    subject = "subject", time = "time", response = "response",
    preparation = c("one row per subject/time", "finite non-negative time", "documented treatment coding"),
    candidate_margins = c("NO", "GA", "LOGNO", "BCPE"),
    candidate_copulas = c("N", "C", "F", "G", "J", "t")
  )
}

jss_validate_lipid_input <- function(data) {
  contract <- jss_lipid_contract()
  missing <- setdiff(contract$required_columns, names(data))
  if (length(missing)) stop("Private LIPID input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(data[c("subject", "time")])) stop("LIPID input must have one row per subject/time.", call. = FALSE)
  if (any(!is.finite(data$time)) || any(data$time < 0)) stop("LIPID time must be finite and non-negative.", call. = FALSE)
  invisible(TRUE)
}

jss_lipid_formulas <- function() {
  list(
    mu = response ~ treatment + time + treatment:time,
    sigma = ~ treatment + time,
    nu = ~ 1,
    tau = ~ 1,
    theta = ~ treatment + time,
    zeta = ~ 1
  )
}

jss_lipid_analysis_recipe <- function(data) {
  jss_validate_lipid_input(data)
  list(
    contract = jss_lipid_contract(), formulas = jss_lipid_formulas(),
    selection = c("screen_margin", "select_copula", "select_joint_distribution"),
    fit = "gamlss_longitudinal", diagnostics = c("check_model", "plot_dist", "plot_copula_diagnostics"),
    publication = c("publication_table", "plot_terms", "plot_margin_fit", "plot_copula_fit")
  )
}

jss_fit_lipid_private <- function(data, margin_dist, copula_dist,
                                  max_outer_iter = 100L, max_inner_iter = 50L) {
  jss_validate_lipid_input(data)
  f <- jss_lipid_formulas()
  gamlss.longitudinal::gamlss_longitudinal(
    dataset = data, margin_dist = margin_dist, copula_dist = copula_dist,
    subject_var = "subject", time_var = "time",
    mu.formula = f$mu, sigma.formula = f$sigma, nu.formula = f$nu,
    tau.formula = f$tau, theta.formula = f$theta, zeta.formula = f$zeta,
    method = "RS", max_outer_iter = max_outer_iter,
    max_inner_iter = max_inner_iter, compute_vcov = TRUE, verbose = 0
  )
}

jss_lipid_private_outputs <- function(fit) {
  list(
    model_check = gamlss.longitudinal::check_model(fit),
    coefficient_table = gamlss.longitudinal::publication_table(fit),
    margin_diagnostics = gamlss.longitudinal::plot_dist(fit),
    copula_diagnostics = gamlss.longitudinal::plot_copula_diagnostics(fit),
    term_plot = gamlss.longitudinal::plot_terms(fit)
  )
}
