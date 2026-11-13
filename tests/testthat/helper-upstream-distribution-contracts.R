zisichel_distribution_contract <- function(
    d_sichel, p_sichel, d_zisichel, logd_zisichel, p_zisichel,
    mu = 2, sigma = 0.5, nu = -0.5, tau = 0.15,
    tolerance = 1e-10) {
  x <- 0:3
  base_mass <- vapply(x, function(value) {
    as.numeric(d_sichel(value, mu = mu, sigma = sigma, nu = nu))
  }, numeric(1))
  base_cdf <- vapply(x, function(value) {
    as.numeric(p_sichel(value, mu = mu, sigma = sigma, nu = nu))
  }, numeric(1))
  expected_mass <- tau * (x == 0) + (1 - tau) * base_mass
  expected_cdf <- tau + (1 - tau) * base_cdf
  observed_mass <- vapply(x, function(value) {
    as.numeric(d_zisichel(value, mu = mu, sigma = sigma, nu = nu, tau = tau))
  }, numeric(1))
  observed_log_mass <- vapply(x, function(value) {
    as.numeric(logd_zisichel(value, mu = mu, sigma = sigma, nu = nu, tau = tau))
  }, numeric(1))
  observed_cdf <- vapply(x, function(value) {
    as.numeric(p_zisichel(value, mu = mu, sigma = sigma, nu = nu, tau = tau))
  }, numeric(1))

  valid_reference <- isTRUE(all(
    is.finite(base_mass), base_mass > 0,
    is.finite(base_cdf), base_cdf >= 0, base_cdf <= 1,
    is.finite(expected_mass), expected_mass > 0, expected_mass < 1,
    is.finite(expected_cdf), expected_cdf >= 0, expected_cdf <= 1,
    diff(c(0, expected_cdf)) >= 0,
    abs(expected_cdf - cumsum(expected_mass)) <= tolerance
  ))
  cdf_matches <- all(is.finite(observed_cdf)) &&
    all(abs(observed_cdf - expected_cdf) <= tolerance)
  zero_mass_matches <- is.finite(observed_mass[[1L]]) &&
    abs(observed_mass[[1L]] - expected_mass[[1L]]) <= tolerance &&
    is.finite(observed_log_mass[[1L]]) &&
    abs(observed_log_mass[[1L]] - log(expected_mass[[1L]])) <= tolerance
  expected_log_mass <- log(expected_mass)
  density_matches <- all(is.finite(observed_mass)) &&
    all(abs(observed_mass - expected_mass) <= tolerance) &&
    all(is.finite(observed_log_mass)) &&
    all(abs(observed_log_mass - expected_log_mass) <= tolerance)
  issue14_signature <- isTRUE(
    valid_reference && cdf_matches && zero_mass_matches &&
      all(observed_mass[x > 0] == 1) && all(observed_log_mass[x > 0] == 0)
  )

  status <- if (!valid_reference) {
    "invalid_reference_contract"
  } else if (density_matches && cdf_matches) {
    "valid"
  } else if (issue14_signature) {
    "known_issue14_density_defect"
  } else {
    "unexpected_dependency_contract_failure"
  }
  list(
    status = status, x = x, expected_mass = expected_mass,
    expected_cdf = expected_cdf, observed_mass = observed_mass,
    observed_log_mass = observed_log_mass, observed_cdf = observed_cdf,
    tolerance = tolerance
  )
}

zisichel_installed_distribution_contract <- function() {
  namespace <- asNamespace("gamlss.dist")
  zisichel_distribution_contract(
    d_sichel = get("dSICHEL", envir = namespace),
    p_sichel = get("pSICHEL", envir = namespace),
    d_zisichel = get("dZISICHEL", envir = namespace),
    logd_zisichel = function(x, ...) {
      get("dZISICHEL", envir = namespace)(x, ..., log = TRUE)
    },
    p_zisichel = get("pZISICHEL", envir = namespace)
  )
}

skip_if_known_zisichel_density_defect <- function(
    probe = zisichel_installed_distribution_contract) {
  contract <- withCallingHandlers(
    probe(),
    warning = function(w) {
      stop(paste("Upstream distribution probe warned:", conditionMessage(w)),
           call. = FALSE)
    }
  )
  if (identical(contract$status, "known_issue14_density_defect")) {
    version <- as.character(utils::packageVersion("gamlss.dist"))
    testthat::skip(paste0(
      "ZISICHEL density contract fails under gamlss.dist ", version,
      " (https://github.com/gamlss-dev/gamlss.dist/issues/14); ",
      "ZISICHEL integration was not exercised."
    ))
  }
  if (!identical(contract$status, "valid")) {
    stop(paste(
      "Unexpected gamlss.dist ZISICHEL distribution contract:", contract$status
    ), call. = FALSE)
  }
  invisible(contract)
}
