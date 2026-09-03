test_that("public candidate datasets retain their audited structure", {
  data("patents_panel", package = "gamlss.longitudinal")
  data("pbc_prothrombin", package = "gamlss.longitudinal")
  data("vietnam_steps", package = "gamlss.longitudinal")

  expect_equal(nrow(patents_panel), 3460L)
  expect_equal(length(unique(patents_panel$firm)), 346L)
  expect_true(all(table(patents_panel$firm) == 10L))
  expect_true(all(patents_panel$patents >= 0L))
  expect_true(all(patents_panel$patents == floor(patents_panel$patents)))

  expect_equal(nrow(pbc_prothrombin), 1945L)
  expect_equal(length(unique(pbc_prothrombin$subject)), 312L)
  expect_false(anyNA(pbc_prothrombin$prothrombin))
  expect_true(all(pbc_prothrombin$prothrombin > 0))

  expect_equal(nrow(vietnam_steps), 3325L)
  expect_equal(length(unique(vietnam_steps$subject)), 475L)
  expect_equal(sum(is.na(vietnam_steps$steps)), 5L)
  expect_equal(sum(vietnam_steps$step_status == "fractional"), 3L)
  expect_equal(sum(vietnam_steps$step_status == "out_of_range"), 2L)
  observed <- !is.na(vietnam_steps$steps)
  expect_true(all(vietnam_steps$steps[observed] == floor(vietnam_steps$steps[observed])))
})

test_that("public candidate sources and licenses are installed", {
  notice <- system.file("DATA-LICENSES.md", package = "gamlss.longitudinal")
  expect_true(nzchar(notice))
  text <- paste(readLines(notice, warn = FALSE), collapse = "\n")
  expect_match(text, "GPL >= 2", fixed = TRUE)
  expect_match(text, "LGPL >= 2", fixed = TRUE)
  expect_match(text, "CC BY 4.0", fixed = TRUE)
})

test_that("public application specifications use supported routes", {
  root <- test_path("..", "..")
  script <- file.path(root, "paper", "R", "10-public-application-candidates.R")
  skip_if_not(file.exists(script), "Paper source is excluded from the installed package.")
  environment <- new.env(parent = globalenv())
  source(script, local = environment)

  for (name in environment$jss_public_application_names()) {
    spec <- environment$jss_public_application_spec(name)
    routes <- longitudinal_capabilities("routes")
    selected <- routes$margin_family == spec$margin & routes$copula == spec$copula
    expect_true(any(selected), info = name)

    prepared <- environment$jss_prepare_public_application(name)
    expect_identical(anyDuplicated(prepared[c(spec$subject, spec$time)]), 0L)
  }

  expect_equal(nrow(environment$jss_prepare_public_application("patents")), 3460L)
  expect_equal(length(unique(environment$jss_prepare_public_application("pbc")$subject)), 227L)
  expect_equal(nrow(environment$jss_prepare_public_application("steps")), 3325L)
})
