test_that("fixed-term y limits pad finite ranges", {
  y_lim <- .plot_fixed_terms_y_limits(c(1, 2), c(0, 1), c(2, 3))

  expect_equal(y_lim, c(-0.15, 3.15))
})

test_that("fixed-term y limits handle constant finite values", {
  y_lim <- .plot_fixed_terms_y_limits(c(1, 1), c(1, 1))

  expect_equal(y_lim, c(1 - 5e-10, 1 + 5e-10))
})

test_that("fixed-term y limits ignore non-finite values and return NULL when empty", {
  expect_equal(.plot_fixed_terms_y_limits(c(NA, -Inf), c(2, Inf)), c(2 - 5e-10, 2 + 5e-10))
  expect_null(.plot_fixed_terms_y_limits(c(NA_real_, Inf, -Inf)))
})

test_that("fixed-term dashboard helper computes default layout", {
  plots <- list("p1", "p2", "p3")

  dashboard <- .plot_fixed_terms_dashboard(plots)

  expect_identical(dashboard$plotlist, plots)
  expect_identical(dashboard$ncol, 2)
  expect_equal(dashboard$nrow, 2)
})

test_that("fixed-term dashboard helper respects explicit layout", {
  dashboard <- .plot_fixed_terms_dashboard(list("p1", "p2", "p3", "p4"), ncol = 3)

  expect_identical(dashboard$ncol, 3)
  expect_equal(dashboard$nrow, 2)
  expect_null(.plot_fixed_terms_dashboard(list()))
})

test_that("fixed-term coefficient predicate requires parameter and vcov names", {
  object <- list(par = c(mu.age = 0.1, mu.sexB = 0.2))
  V <- diag(2)
  rownames(V) <- colnames(V) <- names(object$par)

  expect_true(.plot_fixed_terms_has_coef("mu.age", object, V))
  expect_false(.plot_fixed_terms_has_coef("mu.missing", object, V))

  missing_row <- V[-1, , drop = FALSE]
  expect_false(.plot_fixed_terms_has_coef("mu.age", object, missing_row))

  missing_col <- V[, -1, drop = FALSE]
  expect_false(.plot_fixed_terms_has_coef("mu.age", object, missing_col))
})

test_that("fixed-term coefficient info returns estimate variance and guarded SE", {
  object <- list(par = c(mu.age = 0.1, mu.sexB = 0.2))
  V <- diag(c(0.04, -0.09))
  rownames(V) <- colnames(V) <- names(object$par)

  age <- .plot_fixed_terms_coef_info("mu.age", object, V)
  missing <- .plot_fixed_terms_coef_info("mu.missing", object, V)

  expect_equal(age$estimate, 0.1)
  expect_equal(age$variance, 0.04)
  expect_equal(age$se, 0.2)
  expect_error(
    .plot_fixed_terms_coef_info("mu.sexB", object, V),
    "derived_variance_nonpositive",
    class = "gamlss_longitudinal_inference_unavailable"
  )
  expect_true(all(is.na(unlist(missing))))
})

test_that("fixed-term factor data builder maps reference and coefficient levels", {
  object <- list(par = c(mu.sexB = 0.2))
  V <- matrix(0.01, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  fg <- list(
    levels = c("A", "B"),
    ref_level = "A",
    level_col_map = list(B = "sexB")
  )

  term_data <- .plot_fixed_terms_factor_data(fg, "mu", object, V, z = 2)

  expect_equal(term_data$x_plot, 1:2)
  expect_equal(term_data$fitted, c(0, 0.2))
  expect_equal(term_data$se, c(0, 0.1))
  expect_equal(term_data$ci_lower, c(0, 0))
  expect_equal(term_data$ci_upper, c(0, 0.4))
  expect_true(all(term_data$keep))
  expect_equal(term_data$plot_df$fitted, term_data$fitted)
})

test_that("fixed-term factor-factor data builder returns grouped interaction rows", {
  object <- list(par = c(`mu.aB:bY` = 0.3))
  V <- matrix(0.04, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  ig <- list(
    panel_group = list(var_name = "a", levels = c("A", "B"), ref_level = "A"),
    other_group = list(var_name = "b", levels = c("X", "Y"), ref_level = "X"),
    interaction_col_map = list(B = list(Y = "aB:bY"))
  )

  term_data <- .plot_fixed_terms_factor_factor_data(ig, "mu", object, V, z = 2)

  expect_equal(term_data$x_plot, 1:2)
  expect_equal(term_data$x_labels, c("A", "B"))
  expect_equal(as.character(unique(term_data$plot_df$group)), c("X", "Y"))
  expect_equal(term_data$plot_df$fitted, c(0, 0, 0, 0.3))
  expect_equal(term_data$plot_df$se, c(0, 0, 0, 0.2))
  expect_true(all(term_data$plot_df$keep))
})

test_that("fixed-term continuous data builder sorts and preserves returned vectors", {
  X <- matrix(c(3, 1, 2), ncol = 1)
  colnames(X) <- "age"
  object <- list(
    model_matrix = list(x = list(mu = X)),
    par = c(mu.age = 0.5)
  )
  V <- matrix(0.04, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  spec <- list(col_name = "age", coef_name = "mu.age")

  term_data <- .plot_fixed_terms_continuous_data(spec, "mu", object, V, z = 2, sort_x = TRUE)

  expect_equal(term_data$x_plot, c(3, 1, 2))
  expect_equal(term_data$fitted, c(1.5, 0.5, 1))
  expect_equal(term_data$se, c(0.6, 0.2, 0.4))
  expect_equal(term_data$plot_df$x, c(1, 2, 3))
  expect_equal(term_data$plot_df$fitted, c(0.5, 1, 1.5))
  expect_equal(term_data$xlab_text, "age")
})

test_that("fixed-term factor plot builder returns point-interval panel", {
  object <- list(par = c(mu.sexB = 0.2))
  V <- matrix(0.01, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  fg <- list(
    var_name = "sex",
    levels = c("A", "B"),
    ref_level = "A",
    level_col_map = list(B = "sexB")
  )
  term_data <- .plot_fixed_terms_factor_data(fg, "mu", object, V, z = 2)

  p <- .plot_fixed_terms_factor_plot(
    term_data,
    fg,
    par_name = "mu",
    ci_col = "red",
    fit_col = "black",
    factor_cex = 1,
    ci_level = 0.95,
    show_legend = TRUE
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "mu: sex")
  expect_equal(p$labels$caption, "estimate / 95 % CI")
  expect_true("GeomPoint" %in% layer_classes)
  expect_true("GeomErrorbar" %in% layer_classes)
  expect_false("GeomRibbon" %in% layer_classes)
})

test_that("fixed-term factor-factor plot builder returns grouped line panel", {
  object <- list(par = c(`mu.aB:bY` = 0.3))
  V <- matrix(0.04, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  ig <- list(
    interaction_name = "a:b",
    panel_group = list(var_name = "a", levels = c("A", "B"), ref_level = "A"),
    other_group = list(var_name = "b", levels = c("X", "Y"), ref_level = "X"),
    interaction_col_map = list(B = list(Y = "aB:bY"))
  )
  term_data <- .plot_fixed_terms_factor_factor_data(ig, "mu", object, V, z = 2)

  p <- .plot_fixed_terms_factor_factor_plot(
    term_data,
    ig,
    par_name = "mu",
    factor_cex = 1,
    ci_level = 0.95,
    show_legend = TRUE
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "mu: a:b")
  expect_equal(p$labels$x, "a")
  expect_true("GeomLine" %in% layer_classes)
  expect_true("GeomPoint" %in% layer_classes)
  expect_true("GeomErrorbar" %in% layer_classes)
})

test_that("fixed-term continuous plot builder returns ribbon-line panel", {
  X <- matrix(c(3, 1, 2), ncol = 1)
  colnames(X) <- "age"
  object <- list(
    model_matrix = list(x = list(mu = X)),
    par = c(mu.age = 0.5)
  )
  V <- matrix(0.04, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  spec <- list(col_name = "age", coef_name = "mu.age")
  term_data <- .plot_fixed_terms_continuous_data(spec, "mu", object, V, z = 2, sort_x = TRUE)

  p <- .plot_fixed_terms_continuous_plot(
    term_data,
    par_name = "mu",
    ci_col = "red",
    fit_col = "black",
    fit_lwd = 2,
    ci_level = 0.95,
    show_legend = TRUE
  )

  layer_classes <- vapply(p$layers, function(layer) class(layer$geom)[1], character(1))
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "mu: age")
  expect_equal(p$labels$caption, "fit / 95 % CI")
  expect_true("GeomRibbon" %in% layer_classes)
  expect_true("GeomLine" %in% layer_classes)
  expect_false("GeomErrorbar" %in% layer_classes)
})

test_that("fixed-term CI caption helper preserves panel wording", {
  expect_equal(.plot_fixed_terms_ci_caption("estimate", 0.95), "estimate / 95 % CI")
  expect_equal(.plot_fixed_terms_ci_caption("fit", 0.954), "fit / 95 % CI")
  expect_equal(.plot_fixed_terms_ci_caption("fit", 0.995), "fit / 100 % CI")
})

test_that("fixed-term plot spec helper separates continuous and factor terms", {
  X <- matrix(
    c(
      1, 20, 0,
      1, 30, 1,
      1, 40, 0
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(X) <- c("intercept", "age", "sexB")
  attr(X, "assign") <- c(0L, 1L, 2L)
  attr(X, "term.labels") <- c("age", "sex")

  object <- list(
    model_matrix = list(x = list(mu = X)),
    par = c(mu.intercept = 1, mu.age = 0.1, mu.sexB = 0.2)
  )
  V <- diag(3)
  rownames(V) <- colnames(V) <- names(object$par)
  data_for_terms <- data.frame(
    age = c(20, 30, 40),
    sex = factor(c("A", "B", "A"), levels = c("A", "B"))
  )

  specs <- .plot_fixed_terms_plot_specs(
    object = object,
    V = V,
    data_for_terms = data_for_terms,
    include_intercept = FALSE,
    plot_interactions = FALSE
  )

  expect_equal(vapply(specs, `[[`, character(1), "type"), c("factor", "continuous"))
  expect_equal(specs[[1]]$var_name, "sex")
  expect_equal(specs[[2]]$col_name, "age")
})

test_that("fixed-term factor group helpers clean expressions and map levels", {
  expect_equal(.plot_fixed_terms_clean_expr_name(" factor(`visit`) "), "visit")
  expect_equal(.plot_fixed_terms_clean_expr_name("as.factor(group)"), "group")

  fg <- .plot_fixed_terms_make_factor_group_from_cols(
    term_name = "group",
    term_cols = c("groupB", "groupC"),
    levs = c("A", "B", "C")
  )

  expect_equal(fg$ref_level, "A")
  expect_equal(unlist(fg$level_col_map, use.names = FALSE), c("groupB", "groupC"))
  expect_equal(fg$matched_cols, c("groupB", "groupC"))
})

test_that("fixed-term factor group helpers recover levels from data and expressions", {
  data <- data.frame(
    visit_label = factor(c("baseline", "week1"), levels = c("baseline", "week1")),
    raw_group = c("control", "active")
  )

  expect_equal(
    .plot_fixed_terms_levels_from_data("visit_label", "visit_label", data),
    c("baseline", "week1")
  )
  expect_equal(
    .plot_fixed_terms_levels_from_data("factor(raw_group)", "raw_group", data),
    c("active", "control")
  )
  expect_equal(
    .plot_fixed_terms_levels_from_factor_expr("factor(group)", c("factor(group)B", "factor(group)C")),
    c("factor(group)", "B", "C")
  )
})

test_that("fixed-term fallback factor group helper matches time covariate aliases", {
  data <- data.frame(
    visit_time = factor(c("t0", "t1", "t2"), levels = c("t0", "t1", "t2"))
  )
  x_cols <- c("time_covariatet1", "time_covariatet2")

  fg <- .plot_fixed_terms_fallback_factor_group("visit_time", data, x_cols, groups = list())

  expect_equal(fg$var_name, "visit_time")
  expect_equal(fg$matched_cols, x_cols)
  expect_equal(unlist(fg$level_col_map, use.names = FALSE), x_cols)
})

test_that("fixed-term factor level matcher supports underscore column names", {
  matched <- .plot_fixed_terms_factor_level_col_map(
    levs = c("control", "active group"),
    var_prefixes = "arm",
    x_cols = c("arm_active.group")
  )

  expect_equal(matched$matched_cols, "arm_active.group")
  expect_equal(unlist(matched$level_col_map, use.names = FALSE), "arm_active.group")
})

test_that("fixed-term plot spec helper respects intercept interaction and vcov filters", {
  X <- matrix(
    c(
      1, 20, 1,
      1, 30, 2,
      1, 40, 3
    ),
    ncol = 3,
    byrow = TRUE
  )
  colnames(X) <- c("intercept", "age", "age:time")
  object <- list(
    model_matrix = list(x = list(mu = X)),
    par = c(mu.intercept = 1, mu.age = 0.1, `mu.age:time` = 0.2)
  )
  V <- diag(3)
  rownames(V) <- colnames(V) <- names(object$par)

  default_specs <- .plot_fixed_terms_plot_specs(
    object = object,
    V = V,
    include_intercept = FALSE,
    plot_interactions = FALSE
  )
  full_specs <- .plot_fixed_terms_plot_specs(
    object = object,
    V = V,
    include_intercept = TRUE,
    plot_interactions = TRUE
  )
  missing_vcov <- V[-3, -3, drop = FALSE]
  filtered_specs <- .plot_fixed_terms_plot_specs(
    object = object,
    V = missing_vcov,
    include_intercept = TRUE,
    plot_interactions = TRUE
  )

  expect_equal(vapply(default_specs, `[[`, character(1), "col_name"), "age")
  expect_equal(
    vapply(full_specs, `[[`, character(1), "type"),
    c("continuous", "continuous", "interaction_factor")
  )
  expect_false(any(vapply(filtered_specs, function(x) identical(x$col_name, "age:time"), logical(1))))
})

test_that("fixed-term renderer assembles continuous outputs and plot list", {
  X <- matrix(c(3, 1, 2), ncol = 1)
  colnames(X) <- "age"
  object <- list(
    model_matrix = list(x = list(mu = X)),
    par = c(mu.age = 0.5)
  )
  V <- matrix(0.04, nrow = 1, dimnames = list(names(object$par), names(object$par)))
  plot_specs <- list(list(type = "continuous", par_name = "mu", col_name = "age", coef_name = "mu.age"))

  rendered <- .plot_fixed_terms_render_specs(
    plot_specs = plot_specs,
    object = object,
    V = V,
    z = 2,
    ci_col = "red",
    fit_col = "black",
    ci_level = 0.95,
    factor_cex = 1,
    show_legend = TRUE,
    gg_add = .plot_fixed_terms_gg_add,
    sort_x = TRUE,
    fallback_to_index = TRUE,
    fit_lwd = 2
  )

  expect_named(rendered, c("out", "plot_objects"))
  expect_equal(rendered$out$mu$age$coefficient, "mu.age")
  expect_equal(rendered$out$mu$age$fitted, c(1.5, 0.5, 1))
  expect_equal(length(rendered$plot_objects), 1)
  expect_s3_class(rendered$plot_objects[[1]], "ggplot")
})
