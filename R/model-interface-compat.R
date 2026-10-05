.gl_warn_interface_deprecation <- function(key, message, class) {
  if (identical(getOption("gamlss.longitudinal.lifecycle_verbosity", "warning"), "quiet")) {
    return(invisible())
  }
  if (isTRUE(.gl_deprecation_state[[key]])) {
    return(invisible())
  }
  .gl_deprecation_state[[key]] <- TRUE
  warning(structure(
    list(message = message, call = NULL),
    class = c(class, "warning", "condition")
  ))
  invisible()
}

.gl_legacy_fit_defaults <- function() {
  list(
    dataset = NULL,
    include_dlcopdpar = TRUE,
    check_dlcopdpar_gradient = FALSE,
    inner_stop_crit = NA,
    outer_stop_crit = NA,
    start_step_size = 2,
    step_adjustment = NA,
    max_steps = 5L,
    start_from = NA,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 2L,
    plot_results = FALSE,
    true_val = NA,
    max_outer_iter = 100L,
    max_inner_iter = 20L,
    max_negative_outer_streak = 10L,
    max_elapsed_sec = Inf,
    use_backtracking = TRUE,
    backtracking_max_halves = 2L,
    cg_max_stall = 5L,
    cg_max_delta = 0.5,
    cg_armijo_c1 = 1e-4,
    cg_grad_tol = NA,
    cg_step_tol = NA,
    cg_update_lambda = TRUE,
    cg_lambda_update_every = 10L,
    cg_max_lambda_updates = NA,
    cg_raw_loglik_drop_tol = 10,
    cg_line_search = "best",
    cg_max_line_search_evals = 60L,
    cg_gradient_method = "forward",
    discrete_score_method = "analytical",
    cg_zeta_hessian = "analytical",
    cg_hessian_method = "analytical",
    vcov_numderiv = FALSE,
    use_Rcpp = FALSE,
    lambda_start = NA,
    lambda_penalty_K = 2,
    rs_update_lambda = TRUE,
    rs_smooth_trust_radius = Inf
  )
}

.gl_resolve_legacy_fit_dots <- function(dots) {
  defaults <- .gl_legacy_fit_defaults()
  if (!length(dots)) {
    return(list(values = defaults, supplied = setNames(rep(FALSE, length(defaults)), names(defaults))))
  }
  nms <- names(dots)
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms))) {
    stop("All arguments supplied through `...` must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("Arguments supplied through `...` must not be duplicated.", call. = FALSE)
  }
  unknown <- setdiff(nms, names(defaults))
  if (length(unknown)) {
    stop(
      "Unused argument", if (length(unknown) > 1L) "s" else "", ": ",
      paste0("`", unknown, "`", collapse = ", "), ".",
      call. = FALSE
    )
  }
  defaults[nms] <- dots
  supplied <- setNames(names(defaults) %in% nms, names(defaults))
  list(values = defaults, supplied = supplied)
}

.gl_resolve_fit_data_alias <- function(data, data_supplied, legacy) {
  old_supplied <- isTRUE(legacy$supplied[["dataset"]])
  if (isTRUE(data_supplied) && old_supplied) {
    stop("Supply only `data`; do not also supply its deprecated alias `dataset`.", call. = FALSE)
  }
  if (old_supplied) {
    .gl_warn_interface_deprecation(
      "dataset_argument",
      "`dataset` is deprecated; use `data` instead.",
      "gamlss.longitudinal_deprecated_dataset_argument"
    )
    data <- legacy$values$dataset
  }
  if (is.null(data)) {
    stop("`data` is required.", call. = FALSE)
  }
  data
}

.gl_resolve_fit_start_alias <- function(start, start_supplied, legacy) {
  old_supplied <- isTRUE(legacy$supplied[["start_from"]])
  if (isTRUE(start_supplied) && old_supplied) {
    stop("Supply only `start`; do not also supply its deprecated alias `start_from`.", call. = FALSE)
  }
  if (old_supplied) {
    .gl_warn_interface_deprecation(
      "start_from_argument",
      "`start_from` is deprecated; use `start` instead.",
      "gamlss.longitudinal_deprecated_start_argument"
    )
    return(legacy$values$start_from)
  }
  if (is.null(start)) NA else start
}

.gl_resolve_optimization <- function(method, optimization, optimization_supplied, legacy) {
  old_supplied <- isTRUE(legacy$supplied[["include_dlcopdpar"]])
  if (isTRUE(optimization_supplied) && old_supplied) {
    stop(
      "Supply only `optimization`; do not also supply the deprecated `include_dlcopdpar`.",
      call. = FALSE
    )
  }
  if (old_supplied) {
    old_value <- legacy$values$include_dlcopdpar
    if (!is.logical(old_value) || length(old_value) != 1L || is.na(old_value)) {
      stop("Deprecated `include_dlcopdpar` must be either TRUE or FALSE.", call. = FALSE)
    }
    .gl_warn_interface_deprecation(
      "include_dlcopdpar_argument",
      paste0(
        "`include_dlcopdpar` is deprecated; use `optimization = \"joint\"` ",
        "for TRUE or `optimization = \"separate\"` for FALSE."
      ),
      "gamlss.longitudinal_deprecated_optimization_argument"
    )
    optimization <- if (old_value) "joint" else "separate"
  } else {
    optimization <- match.arg(optimization, c("joint", "separate"))
  }
  if (identical(method, "CG") && identical(optimization, "separate")) {
    stop(
      "`optimization = \"separate\"` is available only with `method = \"RS\"`; ",
      "CG always optimizes the joint objective.",
      call. = FALSE
    )
  }
  optimization
}

.gl_normalize_verbose <- function(verbose) {
  if (is.logical(verbose) && length(verbose) == 1L && !is.na(verbose)) {
    return(if (isTRUE(verbose)) 1L else 0L)
  }
  if (is.numeric(verbose) && length(verbose) == 1L && is.finite(verbose) &&
      verbose %in% 0:3) {
    return(as.integer(verbose))
  }
  stop("`verbose` must be TRUE, FALSE, or an integer from 0 to 3.", call. = FALSE)
}
