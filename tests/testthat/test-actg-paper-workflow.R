actg_workflow_environment <- function() {
  environment <- new.env(parent = globalenv())
  sys.source(
    test_path("..", "..", "paper", "R", "05-actg-workflow.R"),
    envir = environment
  )
  environment
}

test_that("ACTG manuscript summaries reproduce the public cohort", {
  environment <- actg_workflow_environment()
  data <- environment$jss_actg_load_data()
  summary <- environment$jss_actg_data_summary(data)

  expect_equal(nrow(data), 4708L)
  expect_equal(length(unique(data$subject)), 1177L)
  expect_equal(sum(is.finite(data$cd4)), 3352L)
  expect_equal(summary$time_summary$week, c(8, 16, 24, 32))
  expect_equal(summary$time_summary$n, c(906, 982, 695, 769))
  expect_equal(summary$time_summary$n_zero, c(31, 31, 30, 31))

  cohort <- stats::setNames(
    summary$cohort_summary$value,
    summary$cohort_summary$characteristic
  )
  expect_equal(unname(cohort["Participants with all four CD4 counts"]), 439)
  expect_equal(
    unname(cohort["Participants with a leading or intermittent gap"]), 475
  )
  expect_equal(unname(cohort["Contiguous observed response segments"]), 1503)
  expect_equal(unname(cohort["Participants with an intermittent gap"]), 326)
  expect_equal(unname(cohort["Participants with leading unobserved visits"]), 271)
  expect_equal(
    unname(cohort["Participants with terminal missing observations"]), 408
  )

  expect_equal(round(unname(summary$kendall["week_8", "week_16"]), 3), 0.605)
  expect_equal(round(unname(summary$upper_tail["week_8", "week_16"]), 3), 0.492)
  expect_equal(round(unname(summary$lower_tail["week_8", "week_16"]), 3), 0.274)
})

test_that("ACTG selection ledger is compact, portable, and complete", {
  root <- test_path("..", "..")
  directory <- file.path(root, "paper", "data", "public-derived", "actg")
  accepted <- utils::read.csv(
    file.path(directory, "covariate-selection-accepted.csv"),
    stringsAsFactors = FALSE
  )
  targeted <- utils::read.csv(
    file.path(directory, "covariate-selection-targeted.csv"),
    stringsAsFactors = FALSE
  )
  final <- utils::read.csv(
    file.path(directory, "final-formulas.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(accepted), 53L)
  expect_equal(nrow(targeted), 18L)
  expect_setequal(final$parameter, c("mu", "sigma", "nu", "tau", "theta"))
  expect_true(all(accepted$converged))
  expect_true(all(accepted$accepted))
  expect_false(any(grepl("[A-Za-z]:[/\\\\]", unlist(accepted), perl = TRUE)))
  expect_false(any(grepl("[A-Za-z]:[/\\\\]", unlist(targeted), perl = TRUE)))
})

test_that("ACTG paper screening ledger retains the manuscript shortlist", {
  environment <- actg_workflow_environment()
  root <- test_path("..", "..")
  screen <- environment$jss_actg_load_screen_ledger(list(
    public_data_dir = file.path(root, "paper", "data", "public-derived")
  ))

  expect_identical(screen$source, "tracked_public_screen_ledger")
  expect_equal(screen$time_varying$family, c("ZIBNB", "ZIPIG", "ZISICHEL"))
  expect_equal(screen$copula$family_code, c("t", "G", "N", "F", "J", "C"))
  expect_equal(nrow(screen$joint), 6L)
  expect_equal(screen$joint$margin_family[[1L]], "ZIBNB")
  expect_equal(screen$joint$copula_family[[1L]], "G")
  expect_true(all(screen$joint$converged))
})

test_that("ACTG warning audit rejects unclassified events", {
  source <- paste(
    readLines(test_path("..", "..", "paper", "R", "05-actg-workflow.R")),
    collapse = "\n"
  )
  expect_match(source, '"segmented_likelihood_notice"', fixed = TRUE)
  expect_match(source, '"analytical_hessian_fallback"', fixed = TRUE)
  expect_match(source, '"fitted_gradient_review"', fixed = TRUE)
  expect_match(source, 'warning_table$warning_code != "unexpected"', fixed = TRUE)
})

test_that("ACTG zero-inflated controls enforce the production time envelope", {
  environment <- actg_workflow_environment()
  control <- environment$jss_actg_optimizer_control(
    margin_family = "ZISICHEL",
    copula_family = "t",
    screening = TRUE
  )

  expect_equal(control$shared$max_elapsed_sec, 1190)
  expect_identical(control$rs$initial_inner_iter, 5L)
  expect_identical(control$rs$max_inner_iter, 20L)
  expect_identical(control$rs$backtracking_max_halves, 2L)
  expect_identical(control$rs$bfgs_max_iter, 100L)
})

test_that("ACTG four-route resume keys safely reuse stricter time caps", {
  path <- test_path(
    "..", "..", "paper", "reviews", "2026-09-05-medical-candidates",
    "run-zero-inflated-four-route-build.R"
  )
  expect_silent(parse(file = path))
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(
    source,
    "compatible_run_control(completion$run_control_signature)",
    fixed = TRUE
  )
  expect_match(
    source,
    "manifest_control$rs$initial_inner_iter",
    fixed = TRUE
  )
  expect_match(
    source,
    "manifest_control$shared$max_elapsed_sec",
    fixed = TRUE
  )
  expect_match(
    source,
    "as.numeric(stored$max_elapsed_sec) <= requested$max_elapsed_sec",
    fixed = TRUE
  )

  expressions <- parse(file = path)
  helper_environment <- new.env(parent = baseenv())
  for (expression in expressions) {
    assignment <- paste(deparse(expression, width.cutoff = 500L), collapse = " ")
    if (grepl("^compatible_run_control <-", assignment)) {
      eval(expression, envir = helper_environment)
    }
  }
  requested <- list(initial_inner_iter = 2L, max_elapsed_sec = 1500)
  expect_true(helper_environment$compatible_run_control(
    list(initial_inner_iter = 2L, max_elapsed_sec = 1190), requested
  ))
  expect_true(helper_environment$compatible_run_control(requested, requested))
  expect_false(helper_environment$compatible_run_control(
    list(initial_inner_iter = 2L, max_elapsed_sec = 1501), requested
  ))
  expect_false(helper_environment$compatible_run_control(
    list(initial_inner_iter = 1L, max_elapsed_sec = 1190), requested
  ))
})

test_that("ACTG four-route audit keeps the t runtime allowance isolated", {
  path <- test_path(
    "..", "..", "paper", "reviews", "2026-09-05-medical-candidates",
    "audit-zero-inflated-four-route-build.R"
  )
  expect_silent(parse(file = path))
  source <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_match(source, '"t-limit-sec"', fixed = TRUE)
  expect_match(source, 'endsWith(route, "-t")', fixed = TRUE)
  expect_match(source, "all(finite_elapsed < route_limit_sec)", fixed = TRUE)
})

test_that("ACTG manuscript table renderers use only promoted inputs", {
  environment <- actg_workflow_environment()
  data <- environment$jss_actg_load_data()
  summary <- environment$jss_actg_data_summary(data)
  output <- withr::local_tempdir()
  root <- test_path("..", "..")

  cohort <- environment$jss_actg_write_cohort_tex(
    summary$cohort_summary, file.path(output, "cohort.tex")
  )
  selection <- environment$jss_actg_write_selection_tex(
    root, file.path(output, "selection.tex")
  )
  specification <- environment$jss_actg_write_model_specification_tex(
    root, file.path(output, "specification.tex")
  )

  expect_true(all(file.exists(c(cohort, selection, specification))))
  expect_match(paste(readLines(cohort), collapse = "\n"), "1,177", fixed = TRUE)
  expect_match(
    paste(readLines(selection), collapse = "\n"),
    "Panel C: final retained predictors",
    fixed = TRUE
  )
  expect_match(
    paste(readLines(specification), collapse = "\n"),
    "ZISICHEL--Gumbel",
    fixed = TRUE
  )
})

test_that("ACTG workflow is registered as a public paper module", {
  root <- test_path("..", "..")
  manifest <- utils::read.csv(
    file.path(root, "paper", "manifest.csv"),
    stringsAsFactors = FALSE
  )
  seeds <- utils::read.csv(
    file.path(root, "paper", "seeds.csv"),
    stringsAsFactors = FALSE
  )
  actg <- manifest[manifest$producer == "module_05_actg", , drop = FALSE]
  workflow_table <- manifest[
    manifest$artifact_id == "inline_package_workflow",
    , drop = FALSE
  ]
  target_code <- paste(
    readLines(file.path(root, "paper", "_targets.R"), warn = FALSE),
    collapse = "\n"
  )
  entry_scripts <- paste(
    unlist(lapply(
      file.path(root, "paper", "scripts", "actg", c(
        "data-summary.R", "fit-final-model.R"
      )),
      readLines, warn = FALSE
    )),
    collapse = "\n"
  )

  expect_equal(nrow(actg), 18L)
  expect_equal(nrow(workflow_table), 1L)
  expect_identical(workflow_table$access, "non-data-static")
  expect_identical(workflow_table$publication_status, "active")
  expect_true(all(c(
    "actg-screening-pit",
    "actg-family-calibration",
    "actg-joint-simulation",
    "actg-conditional-continuation"
  ) %in% seeds$study))
  expect_true(all(actg$access == "public"))
  expect_true(all(actg$publication_status == "active"))
  expect_true(all(grepl("paper", actg$profiles, fixed = TRUE)))
  expect_match(target_code, "tar_target(actg_primary_fit", fixed = TRUE)
  expect_match(target_code, 'settings$profile, "full"', fixed = TRUE)
  expect_match(target_code, "module_05_actg", fixed = TRUE)
  expect_false(grepl("devtools::load_all", entry_scripts, fixed = TRUE))
  expect_false(grepl("paper/reviews/actg-shortlist", entry_scripts, fixed = TRUE))
})

test_that("manuscript audit records failures and rejects retired references", {
  root <- test_path("..", "..")
  audit_code <- readLines(
    file.path(root, "paper", "audit-manuscript.R"),
    warn = FALSE
  )
  write_line <- grep("utils::write.csv(audit, out", audit_code, fixed = TRUE)
  unclassified_line <- grep("if (any(!audit$classified))", audit_code, fixed = TRUE)

  expect_length(write_line, 1L)
  expect_length(unclassified_line, 1L)
  expect_lt(write_line, unclassified_line)
  expect_true(any(grepl(
    "Retired artifacts are still referenced by the manuscript",
    audit_code,
    fixed = TRUE
  )))
})
