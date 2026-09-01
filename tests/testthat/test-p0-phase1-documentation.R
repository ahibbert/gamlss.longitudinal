test_that("Phase 1 documentation retains the inference and diagnostic contracts", {
  root <- testthat::test_path("..", "..")
  manuscript <- paste(readLines(
    file.path(root, "paper", "manuscript", "main.tex"),
    warn = FALSE
  ), collapse = "\n")
  standard <- paste(readLines(
    file.path(root, "vignettes", "standard-workflow.Rmd"),
    warn = FALSE
  ), collapse = "\n")
  site_guide <- paste(readLines(
    file.path(root, "vignettes", "site-guide.Rmd"),
    warn = FALSE
  ), collapse = "\n")
  fit_docs <- paste(readLines(
    file.path(root, "R", "model-fit.R"),
    warn = FALSE
  ), collapse = "\n")

  expect_false(grepl("sqrt{\\\\left|\\\\left(H", manuscript, fixed = TRUE))
  expect_false(grepl("flags failure if", manuscript, fixed = TRUE))
  expect_false(grepl("failed-check warnings", paste(standard, site_guide)))
  expect_false(grepl("residual dependence warnings", standard, fixed = TRUE))
  expect_false(grepl("tail fit warnings", standard, fixed = TRUE))
  expect_false(grepl("dependence_cor_cutoff = 0.25", standard, fixed = TRUE))
  expect_false(grepl("leading unobserved visits do not require", fit_docs, fixed = TRUE))
  expect_false(grepl("filling the entire joint Hessian matrix", manuscript, fixed = TRUE))
  expect_false(grepl("variance/covariance calculation for any fit", manuscript, fixed = TRUE))
  expect_false(grepl("based on the full hessian", manuscript, fixed = TRUE))

  expect_match(manuscript, "conditional on fitted smooth coefficients", fixed = TRUE)
  expect_match(manuscript, "conditional covariance for the fixed margin and copula coefficient", fixed = TRUE)
  expect_match(manuscript, "Discrete margins use seeded randomized PIT values", fixed = TRUE)
  expect_match(manuscript, "missingness = \"segment\"", fixed = TRUE)
  expect_match(fit_docs, "observed prefixes", fixed = TRUE)
})
