# Fit-free correctness fixtures; these have NOT been executed in this candidate.
testthat::test_that("newdata smooth predictions retain fitted coordinates", {
  data <- data.frame(response = seq_len(80), x = seq(0, 1, length.out = 80))
  env <- new.env(parent = baseenv()); env$s <- mgcv::s
  form <- stats::as.formula('response ~ s(x, bs = "ps", k = 8)', env = env)
  retained <- .gl_build_smooth_model_matrices(form, data, env)
  term <- names(retained)[[1L]]; B <- retained[[term]]
  reference <- list(x = list(mu = matrix(1, nrow(data), 1L)), s = list(mu = retained))
  before_kind <- RNGkind(); before_seed <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
  whole <- .gl_build_prediction_smooth_model_matrices("mu", form, data, reference)[[term]]
  testthat::expect_equal(as.vector(whole), as.vector(B), tolerance = 1e-12)
  groups <- list(1:11, 12:46, 47:80)
  partition <- do.call(rbind, lapply(groups, function(index)
    .gl_build_prediction_smooth_model_matrices("mu", form, data[index,,drop = FALSE], reference)[[term]]))
  testthat::expect_equal(as.vector(partition), as.vector(whole), tolerance = 1e-12)
  single <- .gl_build_prediction_smooth_model_matrices("mu", form, data[37,,drop = FALSE], reference)[[term]]
  testthat::expect_equal(as.vector(single), as.vector(whole[37,,drop = FALSE]), tolerance = 1e-12)
  testthat::expect_identical(RNGkind(), before_kind)
  testthat::expect_identical(if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL, before_seed)
  legacy <- reference; attr(legacy$s$mu[[term]], "gl_smooth_spec") <- NULL
  testthat::expect_error(.gl_build_prediction_smooth_model_matrices("mu", form, data, legacy),
    class = "gamlss_longitudinal_prediction_unavailable")
  changed <- reference; changed$s$mu[[term]][1,1] <- changed$s$mu[[term]][1,1] + 1
  testthat::expect_error(.gl_build_prediction_smooth_model_matrices("mu", form, data, changed),
    class = "gamlss_longitudinal_prediction_unavailable")
  wrong_terms <- reference; names(wrong_terms$s$mu) <- "s(other)"
  testthat::expect_error(.gl_build_prediction_smooth_model_matrices("mu", form, data, wrong_terms),
    class = "gamlss_longitudinal_prediction_unavailable")
})


testthat::test_that("public newdata prediction preserves a complete finite design and refuses missing reference", {
  # A constructed coefficient/design fixture, not an estimated model or a simulation study.
  # This checks marginal outputs only; general theta/zeta panels and imputation are not batch-invariant claims.
  data <- data.frame(subject = rep(1:20, each = 4), time = rep(1:4, 20),
    response = seq_len(80) / 10, x = seq(0, 1, length.out = 80),
    treatment = factor(rep(c("a", "b"), each = 40)))
  env <- new.env(parent = baseenv()); env$s <- mgcv::s
  formulas <- list(mu = stats::as.formula('response ~ treatment + s(x, bs = "ps", k = 8)', env = env),
    sigma = ~ 1, nu = ~ 1, tau = ~ 1, theta = ~ 1, zeta = ~ 1)
  arguments <- c(setNames(formulas, paste0(names(formulas), ".formula")),
    list(margin.family = gamlss.dist::NO(), copula.family = "G",
      copula.link = get_copula_dist("G")$copula_link, dataset = data))
  mm <- do.call(create_model_matrices, arguments)
  fixed <- unlist(lapply(names(mm$x), function(parameter)
    setNames(rep(0.2, ncol(mm$x[[parameter]])), paste(parameter, colnames(mm$x[[parameter]]), sep = "."))))
  # unlist(lapply(...)) leaves the explicitly constructed coefficient names unchanged.
  smooth <- lapply(mm$s, function(terms) lapply(terms, function(B) seq_len(ncol(B)) / 10))
  object <- structure(list(model_matrix = mm, formulas_int = formulas, var_map = character(),
    dataset = data, dataset_original = data, margin_dist = gamlss.dist::NO(), copula_dist = "G",
    par = fixed, par_s = smooth, response = data$response, response_margin = data$time,
    response_subject = data$subject), class = "gamlss.longitudinal")
  whole <- stats::predict(object, newdata = data, type = "parameters", se.fit = FALSE, interval = "none")
  order <- c(seq(80, 2, by = -2), seq(79, 1, by = -2))
  groups <- split(order, rep(1:8, each = 10))
  parts <- lapply(groups, function(index) stats::predict(object, newdata = data[index,,drop = FALSE],
    type = "parameters", se.fit = FALSE, interval = "none"))
  combined <- do.call(rbind, parts)
  testthat::expect_equal(as.matrix(combined[c("mu", "sigma")]),
    as.matrix(whole[order,c("mu", "sigma")]), tolerance = 1e-12, ignore_attr = TRUE)
  testthat::expect_identical(as.character(combined$subject), as.character(data$subject[order]))
  testthat::expect_equal(stats::predict(object, newdata = data, type = "mean"),
    stats::predict(object, newdata = data, type = "mu"), tolerance = 1e-12)
  quantiles <- stats::predict(object, newdata = data, type = "quantile", probs = c(0.1, 0.5, 0.9))
  testthat::expect_equal(quantiles$q05, whole$mu, tolerance = 1e-12)
  broken <- object; broken$model_matrix <- NULL
  fresh_calls <- 0L
  testthat::local_mocked_bindings(.gl_build_smooth_model_matrices = function(...) {
    fresh_calls <<- fresh_calls + 1L; stop("Fresh smooth construction must not occur")
  }, .package = "gamlss.longitudinal")
  testthat::expect_error(stats::predict(broken, newdata = data, type = "mu"),
    class = "gamlss_longitudinal_prediction_unavailable")
  testthat::expect_identical(fresh_calls, 0L)
})
