local_jss_repo_root <- function() {
  repo_root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(
    file.exists(file.path(repo_root, "paper", "manifest.csv")),
    "paper replication sources are excluded from source-package checks"
  )
  repo_root
}

test_that("JSS replication skeleton declares eight paper modules", {
  repo_root <- local_jss_repo_root()
  module_files <- file.path(
    repo_root, "paper", "R",
    c(
      "01-simulation-bcpe-t.R",
      "02-simulation-delaporte-clayton.R",
      "03-joint-vs-separate-optimization.R",
      "04-missingness-dropout-sensitivity.R",
      "05-application-lipid.R",
      "06-application-rand-doctor-visits.R",
      "07-gamma-copula-misspecification.R",
      "08-simulation-sensitivity-correlation-misspecification.R"
    )
  )
  expect_true(all(file.exists(module_files)))

  manifest <- utils::read.csv(file.path(repo_root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  expected_modules <- sub("[.]R$", "", basename(module_files))
  manifest_modules <- setdiff(unique(manifest$module_id), "cross-module")
  expect_setequal(manifest_modules, expected_modules)

  module_rows <- manifest[manifest$module_id %in% expected_modules, , drop = FALSE]
  for (module_id in expected_modules) {
    module_types <- module_rows$result_type[module_rows$module_id == module_id]
    expect_true(all(c("data", "table", "figure") %in% module_types))
  }
  expect_true(all(module_rows$analysis_state %in% c("stub", "current")))
  expect_true(all(grepl("results/jss-replication/<profile>/", module_rows$output_path, fixed = TRUE)))
})

test_that("gamma copula mis-specification grid has smoke, pilot, and full stages", {
  repo_root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)

  settings <- list(profile = "expanded", seed = 20260528L)
  smoke <- env$jss_misspec_config(settings, stage = "smoke")
  pilot <- env$jss_misspec_config(settings, stage = "pilot")
  full <- env$jss_misspec_config(settings, stage = "full")

  expect_equal(nrow(env$jss_misspec_grid(smoke)), 36L)
  expect_equal(nrow(env$jss_misspec_grid(pilot)), 6L * 6L * 2L * 3L * 10L)
  expect_equal(nrow(env$jss_misspec_grid(full)), 6L * 6L * 2L * 3L * 100L)
  expect_equal(pilot$sample_sizes, c(50L, 150L, 500L))
  expect_equal(pilot$tau_levels$target_tau, c(0.25, 0.55))
})

test_that("gamma copula mis-specification full stage can resume from pilot checkpoints", {
  repo_root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "paper", "R", "07-gamma-copula-misspecification.R"), local = env)

  settings <- list(profile = "expanded", seed = 20260528L)
  pilot <- env$jss_misspec_config(settings, stage = "pilot")
  full <- env$jss_misspec_config(settings, stage = "full")
  pilot_grid <- env$jss_misspec_grid(pilot)
  full_grid <- env$jss_misspec_grid(full)

  checkpoint_dir <- file.path(tempdir(), paste0("misspec-checkpoints-", sample.int(1e6, 1)))
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  completed_id <- pilot_grid$fit_id[[1L]]
  utils::write.csv(
    data.frame(fit_id = completed_id, success = TRUE, stringsAsFactors = FALSE),
    env$jss_misspec_checkpoint_path(checkpoint_dir, completed_id),
    row.names = FALSE
  )

  pending <- env$jss_misspec_pending_grid(full_grid, checkpoint_dir)

  expect_false(completed_id %in% pending$fit_id)
  expect_equal(nrow(pending), nrow(full_grid) - 1L)
  expect_true(completed_id %in% full_grid$fit_id)
})

test_that("private JSS application data are represented by environment contracts", {
  repo_root <- local_jss_repo_root()
  source(file.path(repo_root, "paper", "R", "replication-helpers.R"), local = TRUE)

  expect_equal(jss_external_data_status("GAMLSS_LONGITUDINAL_LIPID_DATA")$envvar, "GAMLSS_LONGITUDINAL_LIPID_DATA")
  expect_equal(jss_external_data_status("GAMLSS_LONGITUDINAL_RAND_DATA")$envvar, "GAMLSS_LONGITUDINAL_RAND_DATA")

  manifest <- utils::read.csv(file.path(repo_root, "paper", "manifest.csv"), stringsAsFactors = FALSE)
  private_modules <- unique(manifest$module_id[manifest$private_data == "yes"])
  expect_setequal(private_modules, c("05-application-lipid", "06-application-rand-doctor-visits"))
})

test_that("joint-vs-separate module declares the fixed seven-case design", {
  repo_root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  cases <- env$jss_joint_case_definitions()

  expect_equal(nrow(cases), 7L)
  expect_identical(cases$case_id, sprintf("JVS%02d", 1:7))
  expect_setequal(cases$family, c("NO", "NBI"))
  expect_true(all(cases$copula == "N"))
  expect_equal(cases$total_observations, cases$n * cases$time_points)
  expect_equal(cases$hypothesis_role[[3L]], "central_joint_win")
  expect_equal(cases$hypothesis_role[[7L]], "discrete_analogue")
})

test_that("joint-vs-separate deltas are joint minus separate with expected signs", {
  repo_root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  base <- env$jss_joint_case_definitions()[1, , drop = FALSE]
  row_for <- function(method, train_ll, test_score, vario_p05, vario_p2, theta_rmse, tau_rmse) {
    data.frame(
      case_id = base$case_id,
      hypothesis_role = base$hypothesis_role,
      joint_review_rep = 1L,
      family = base$family,
      copula = base$copula,
      design = "seven_case",
      n = base$n,
      time_points = base$time_points,
      n_subject = base$n,
      n_time = base$time_points,
      total_observations = base$total_observations,
      mu_strength = base$mu_strength,
      sigma_strength = base$sigma_strength,
      theta_strength = base$theta_strength,
      time_shape = base$time_shape,
      dependence = "case_specific",
      missingness = "none",
      start_mode = "default",
      method = method,
      success = TRUE,
      train_joint_loglik = train_ll,
      test_log_score_per_obs = test_score,
      heldout_variogram_score_p05 = vario_p05,
      heldout_variogram_score_p2 = vario_p2,
      train_rmse_theta = theta_rmse,
      train_rmse_tau = tau_rmse,
      elapsed_sec = 1,
      stringsAsFactors = FALSE
    )
  }
  results <- rbind(
    row_for("rs_separate", 10, -0.4, 0.30, 0.20, 0.20, 0.10),
    row_for("rs_joint", 13, -0.3, 0.25, 0.12, 0.15, 0.08)
  )

  deltas <- env$jss_joint_delta_table(results)

  expect_equal(nrow(deltas), 1L)
  expect_equal(deltas$delta_train_joint_loglik, 3)
  expect_equal(deltas$delta_test_log_score_per_obs, 0.1)
  expect_equal(deltas$delta_heldout_variogram_score_p05, -0.05)
  expect_equal(deltas$delta_heldout_variogram_score_p2, -0.08)
  expect_equal(deltas$delta_train_rmse_theta, -0.05)
  expect_equal(deltas$delta_train_rmse_tau, -0.02)
})

test_that("joint-vs-separate hypothesis summary classifies gate outcomes", {
  repo_root <- local_jss_repo_root()
  env <- new.env(parent = globalenv())
  source(file.path(repo_root, "paper", "R", "03-joint-vs-separate-optimization.R"), local = env)

  cases <- env$jss_joint_case_definitions()
  delta_for <- function(case_id, train, test, vario_p05, vario_p2, theta, tau, success = TRUE) {
    case <- cases[cases$case_id == case_id, , drop = FALSE]
    data.frame(
      case_id = case$case_id,
      hypothesis_role = case$hypothesis_role,
      family = case$family,
      copula = case$copula,
      design = "seven_case",
      n_subject = case$n,
      n_time = case$time_points,
      total_observations = case$total_observations,
      mu_strength = case$mu_strength,
      sigma_strength = case$sigma_strength,
      theta_strength = case$theta_strength,
      time_shape = case$time_shape,
      rs_joint_success = success,
      rs_separate_success = success,
      delta_train_joint_loglik = train,
      delta_test_log_score_per_obs = test,
      delta_heldout_variogram_score_p05 = vario_p05,
      delta_heldout_variogram_score_p2 = vario_p2,
      delta_train_rmse_theta = theta,
      delta_train_rmse_tau = tau,
      stringsAsFactors = FALSE
    )
  }
  deltas <- do.call(rbind, list(
    delta_for("JVS01", 0.1, 0.0000, 0.001, 0.001, 0.001, 0.001),
    delta_for("JVS02", 1.0, 0.0010, -0.010, -0.005, 0.010, 0.005),
    delta_for("JVS03", 5.0, 0.0040, -0.050, -0.030, 0.050, 0.030),
    delta_for("JVS04", 2.0, 0.0030, -0.010, -0.004, 0.010, 0.004),
    delta_for("JVS05", 2.0, 0.0020, -0.015, -0.006, 0.015, 0.006),
    delta_for("JVS06", 2.0, 0.0010, -0.010, -0.004, 0.010, 0.004),
    delta_for("JVS07", 0.4, 0.0002, -0.002, -0.001, 0.002, 0.001)
  ))

  summary <- env$jss_joint_hypothesis_summary(deltas)

  expect_equal(nrow(summary), 6L)
  expect_true(all(summary$decision %in% c("confirmed", "disconfirmed", "inconclusive")))
  expect_true(summary$focal_gate_pass[summary$focal_case == "JVS03"][[1L]])
  expect_equal(
    summary$decision[summary$hypothesis == "Stronger dependence improves joint-vs-separate performance"],
    "confirmed"
  )
})
