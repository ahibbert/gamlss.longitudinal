make_test_zisichel_contract_functions <- function(density_mode = "valid") {
  base_mass <- c(0.2, 0.3, 0.3, 0.15)
  base_cdf <- cumsum(base_mass)
  d_sichel <- function(x, mu, sigma, nu) base_mass[x + 1L]
  p_sichel <- function(x, mu, sigma, nu) if (x < 0) 0 else base_cdf[x + 1L]
  mixture_mass <- function(x, mu, sigma, nu, tau) {
    tau * (x == 0) + (1 - tau) * base_mass[x + 1L]
  }
  mixture_cdf <- function(x, mu, sigma, nu, tau) {
    tau + (1 - tau) * base_cdf[x + 1L]
  }
  d_zisichel <- switch(
    density_mode,
    valid = mixture_mass,
    issue14 = function(x, mu, sigma, nu, tau) {
      if (x > 0) 1 else mixture_mass(x, mu, sigma, nu, tau)
    },
    unexpected = function(x, mu, sigma, nu, tau) {
      if (x > 0) 0.5 else mixture_mass(x, mu, sigma, nu, tau)
    }
  )
  logd_zisichel <- function(x, mu, sigma, nu, tau) {
    value <- d_zisichel(x, mu, sigma, nu, tau)
    if (density_mode == "issue14" && x > 0) return(0)
    log(value)
  }
  list(
    d_sichel = d_sichel, p_sichel = p_sichel,
    d_zisichel = d_zisichel, logd_zisichel = logd_zisichel,
    p_zisichel = mixture_cdf
  )
}

evaluate_test_zisichel_contract <- function(functions) {
  do.call(zisichel_distribution_contract, functions)
}

test_that("ZISICHEL dependency guard recognizes only the exact density defect", {
  known_broken <- evaluate_test_zisichel_contract(
    make_test_zisichel_contract_functions("issue14")
  )
  expect_identical(known_broken$status, "known_issue14_density_defect")

  fixed <- evaluate_test_zisichel_contract(
    make_test_zisichel_contract_functions("valid")
  )
  expect_identical(fixed$status, "valid")

  unexpected <- evaluate_test_zisichel_contract(
    make_test_zisichel_contract_functions("unexpected")
  )
  expect_identical(unexpected$status, "unexpected_dependency_contract_failure")

  wrong_cdf <- make_test_zisichel_contract_functions("issue14")
  original_cdf <- wrong_cdf$p_zisichel
  wrong_cdf$p_zisichel <- function(...) original_cdf(...) + 0.01
  expect_identical(evaluate_test_zisichel_contract(wrong_cdf)$status,
                   "unexpected_dependency_contract_failure")

  wrong_zero_log <- make_test_zisichel_contract_functions("issue14")
  wrong_zero_log$logd_zisichel <- function(...) 0
  expect_identical(evaluate_test_zisichel_contract(wrong_zero_log)$status,
                   "unexpected_dependency_contract_failure")

  expect_error(
    skip_if_known_zisichel_density_defect(probe = function() {
      warning("unexpected provider condition", call. = FALSE)
      known_broken
    }),
    "Upstream distribution probe warned: unexpected provider condition"
  )
})
