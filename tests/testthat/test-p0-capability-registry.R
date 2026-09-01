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
    "likelihood_route", "compatible_copulas", "hessian_path", "hessian_validity_guaranteed",
    "curvature_validation", "hessian_evidence", "randomized_pit",
    "diagnostics", "paper_route", "limitations"
  ))
  expect_named(copulas, c(
    "registry_version", "copula", "status", "parameters", "links",
    "inverse_links", "parameter_domain", "limitations"
  ))
  expect_equal(
    margins$family[margins$status == "supported"],
    c("NO", "GA", "GG", "BCPE", "LOGNO", "PO", "NBI", "DEL", "ZIP", "ZAP", "ZINBI")
  )
  expect_equal(nrow(routes), 19L)
  expect_equal(copulas$copula, c("N", "C", "F", "G", "J", "t"))
  expect_true(all(nzchar(margins$response_requirement)))
  expect_true(all(nzchar(margins$limitations)))
  expect_true(all(nzchar(copulas$parameter_domain)))
  expect_true(all(routes$curvature_validation == "required_per_fit_jss002"))
  expect_false(any(routes$hessian_validity_guaranteed))
  expect_true(all(nzchar(routes$route_evidence)))
  expect_true(all(nzchar(routes$hessian_evidence)))
  expect_true(all(nzchar(routes$diagnostics)))
  expect_true(all(nzchar(routes$diagnostic_evidence)))
  expect_true(all(grepl("not guaranteed valid inference", routes$limitations, fixed = TRUE)))
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
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(1, Inf)),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_error(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(-Inf, NA)),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
  )
  expect_invisible(
    gamlss.longitudinal:::.gl_validate_capability_route(gamlss.dist::NO(), "N", response = c(NA, NA))
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
  expect_error(
    gamlss.longitudinal:::.gl_preflight_fit_capabilities(
      transform(dat, y = c(1, Inf)), gamlss.dist::NO(), "N", y ~ 1
    ),
    "Inf", class = "gamlss_longitudinal_response_domain_error"
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
