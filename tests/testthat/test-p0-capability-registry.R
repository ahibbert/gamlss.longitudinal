test_that("capability registry is versioned, complete, and documentation-ready", {
  margins <- longitudinal_capabilities("margins")
  copulas <- longitudinal_capabilities("copulas")
  routes <- longitudinal_capabilities("routes")

  expect_s3_class(margins, "data.frame")
  expect_s3_class(copulas, "data.frame")
  expect_s3_class(routes, "data.frame")
  expect_true(all(nzchar(margins$registry_version)))
  expect_equal(unique(margins$registry_version), unique(routes$registry_version))
  expect_equal(unique(copulas$registry_version), unique(routes$registry_version))
  expect_named(margins, c(
    "registry_version", "family", "status", "family_type", "response_domain",
    "response_requirement", "parameters", "default_links", "inverse_links",
    "likelihood_route", "compatible_copulas", "hessian", "randomized_pit",
    "diagnostics", "paper_route", "limitations"
  ))
  expect_named(copulas, c(
    "registry_version", "copula", "status", "parameters", "links",
    "inverse_links", "parameter_domain", "limitations"
  ))
  expect_true(all(c("NO", "GA", "GG", "BCPE", "LOGNO", "PO", "NBI", "DEL", "ZIP", "ZAP", "ZINBI") %in% margins$family))
  expect_equal(copulas$copula, c("N", "C", "F", "G", "J", "t"))
  expect_true(all(nzchar(margins$response_requirement)))
  expect_true(all(nzchar(margins$limitations)))
  expect_true(all(nzchar(copulas$parameter_domain)))
})

test_that("registry identifies both mandatory public paper routes", {
  routes <- longitudinal_capabilities("routes")
  paper <- routes[routes$paper_route, c("margin_family", "copula", "likelihood_route")]

  expect_equal(nrow(paper), 2L)
  expect_true(any(paper$margin_family == "BCPE" & paper$copula == "t" & paper$likelihood_route == "continuous_density"))
  expect_true(any(paper$margin_family == "NBI" & paper$copula == "C" & paper$likelihood_route == "exact_discrete_rectangle"))
})

test_that("capability tables are direct projections of the registry specs", {
  margins <- longitudinal_capabilities("margins")
  routes <- longitudinal_capabilities("routes")
  supported_specs <- gamlss.longitudinal:::.gl_capability_margin_specs()

  expect_equal(sort(margins$family[margins$status == "supported"]), sort(names(supported_specs)))
  for (family in names(supported_specs)) {
    expected <- supported_specs[[family]]$compatible_copulas
    observed <- routes$copula[routes$margin_family == family]
    expect_equal(observed, expected, info = family)
  }

  rd <- readLines(test_path("..", "..", "man", "longitudinal_capabilities.Rd"), warn = FALSE)
  expect_true(any(grepl("BCPE", rd, fixed = TRUE)))
  expect_true(any(grepl("NBI", rd, fixed = TRUE)))
  expect_true(any(grepl("one homogeneous marginal family", tolower(rd), fixed = TRUE)))
})

test_that("registered routes and likelihood routing agree", {
  skip_if_not_installed("gamlss.dist")

  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::BCPE(), "t", response = c(1, 2)
  ))
  expect_invisible(gamlss.longitudinal:::.gl_validate_capability_route(
    gamlss.dist::NBI(), "C", response = c(0, 2)
  ))
  expect_identical(
    gamlss.longitudinal:::.gl_capability_likelihood_route(gamlss.dist::BCPE()),
    "continuous_density"
  )
  expect_identical(
    gamlss.longitudinal:::.gl_capability_likelihood_route(gamlss.dist::NBI()),
    "exact_discrete_rectangle"
  )
})

test_that("unsupported, mixed, and invalid-domain inputs fail with named conditions", {
  skip_if_not_installed("gamlss.dist")

  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::BI(), "N"),
    "denominator", class = "gamlss_longitudinal_unsupported_margin_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::BCPE(), "C"),
    "tested allowlist", class = "gamlss_longitudinal_unsupported_route_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(list(gamlss.dist::NO(), gamlss.dist::PO()), "N"),
    "one homogeneous", class = "gamlss_longitudinal_mixed_margin_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::GA(), "N", response = c(0, 1)),
    "strictly greater", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NBI(), "C", response = c(0, 1.5)),
    "non-negative integers", class = "gamlss_longitudinal_response_domain_error"
  )
})

test_that("fit preflight rejects unsupported margins before the workflow", {
  skip_if_not_installed("gamlss.dist")

  dat <- data.frame(id = c(1, 1), time = c(1, 2), y = c(0, 1))
  expect_error(
    gamlss.longitudinal:::.gl_preflight_fit_capabilities(
      dat, gamlss.dist::BI(), "N", y ~ 1
    ),
    class = "gamlss_longitudinal_unsupported_margin_error"
  )
})

test_that("coverage catalog is audit evidence constrained by the registry", {
  skip_if_not_installed("gamlss.dist")

  catalog <- gamlss.longitudinal:::.coverage_family_catalog()
  supported <- sort(catalog$family[catalog$supported])
  registered <- sort(longitudinal_capabilities("margins")$family[
    longitudinal_capabilities("margins")$status == "supported"
  ])
  expect_true(all(supported %in% registered))
  expect_false(catalog$supported[match("BI", catalog$family)])
})
