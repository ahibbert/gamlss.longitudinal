#' Optimizer controls for longitudinal GAMLSS fits
#'
#' Create a validated, named control object for [gamlss_longitudinal()]. Common
#' controls apply to both optimizers. Method-specific controls live in the `rs`
#' and `cg` lists so that optimizer settings cannot silently leak across routes.
#'
#' @param outer_tol Positive outer objective-change tolerance, or `NA` for the
#'   data-adaptive default.
#' @param max_outer_iter Maximum number of outer optimizer iterations.
#' @param max_elapsed_sec Maximum elapsed fitting time in seconds; `Inf` means
#'   no time limit.
#' @param stop_on_convergence If `TRUE`, stop as soon as the method's convergence
#'   contract is satisfied. If `FALSE`, continue to a safeguard or iteration
#'   limit while still reporting whether the final state satisfies the contract.
#' @param rs Named list of RS-only controls. See Details.
#' @param cg Named list of CG-only controls. See Details.
#'
#' @details
#' RS controls are `inner_tol`, `max_inner_iter`, `max_negative_outer_streak`,
#' `start_step_size`, `step_adjustment`, `max_steps`, `warm_start_joint`,
#' `warm_start_joint_iter`, `use_backtracking`, `backtracking_max_halves`,
#' `update_lambda`, `smooth_trust_radius`, and `discrete_score_method`.
#'
#' CG controls are `max_stall`, `max_delta`, `armijo_c1`, `grad_tol`,
#' `step_tol`, `update_lambda`, `lambda_update_every`, `max_lambda_updates`,
#' `raw_loglik_drop_tol`, `line_search`, `max_line_search_evals`,
#' `gradient_method`, `zeta_hessian`, `hessian_method`, `use_backtracking`, and
#' `backtracking_max_halves`. `grad_tol` is an optimizer stopping tolerance and
#' is distinct from `inference_control(gradient_tol = ...)`, which validates
#' post-fit curvature.
#'
#' @return A `gamlss.longitudinal.control` object.
#' @examples
#' rs_control <- gamlss_longitudinal_control(
#'   outer_tol = 1e-3,
#'   rs = list(inner_tol = 1e-4, max_inner_iter = 50)
#' )
#' cg_control <- gamlss_longitudinal_control(
#'   outer_tol = 1e-4,
#'   cg = list(grad_tol = 1e-3, step_tol = 1e-5)
#' )
#' @export
gamlss_longitudinal_control <- function(
    outer_tol = NA,
    max_outer_iter = 100L,
    max_elapsed_sec = Inf,
    stop_on_convergence = TRUE,
    rs = list(),
    cg = list()) {
  specified <- list(
    shared = c(
      outer_tol = !missing(outer_tol),
      max_outer_iter = !missing(max_outer_iter),
      max_elapsed_sec = !missing(max_elapsed_sec),
      stop_on_convergence = !missing(stop_on_convergence)
    ),
    rs = .gl_control_list_names(rs, "rs"),
    cg = .gl_control_list_names(cg, "cg")
  )

  rs_values <- .gl_merge_control_section(.gl_rs_control_defaults(), rs, "rs")
  cg_values <- .gl_merge_control_section(.gl_cg_control_defaults(), cg, "cg")
  shared <- list(
    outer_tol = outer_tol,
    max_outer_iter = max_outer_iter,
    max_elapsed_sec = max_elapsed_sec,
    stop_on_convergence = stop_on_convergence
  )
  shared <- .gl_validate_shared_optimizer_controls(shared)

  structure(
    list(shared = shared, rs = rs_values, cg = cg_values),
    specified = specified,
    class = c("gamlss.longitudinal.control", "list")
  )
}

.gl_rs_control_defaults <- function() {
  list(
    inner_tol = NA,
    max_inner_iter = 100L,
    max_negative_outer_streak = 10L,
    start_step_size = 0.5,
    step_adjustment = NA,
    max_steps = 5L,
    warm_start_joint = TRUE,
    warm_start_joint_iter = 5L,
    use_backtracking = TRUE,
    backtracking_max_halves = 50L,
    update_lambda = TRUE,
    smooth_trust_radius = Inf,
    discrete_score_method = "analytical"
  )
}

.gl_cg_control_defaults <- function() {
  list(
    max_stall = 5L,
    max_delta = 0.5,
    armijo_c1 = 1e-4,
    grad_tol = NA,
    step_tol = NA,
    update_lambda = TRUE,
    lambda_update_every = 10L,
    max_lambda_updates = NA,
    raw_loglik_drop_tol = 10,
    line_search = "best",
    max_line_search_evals = 60L,
    gradient_method = "forward",
    zeta_hessian = "analytical",
    hessian_method = "analytical",
    use_backtracking = TRUE,
    backtracking_max_halves = 50L
  )
}

.gl_control_list_names <- function(x, section) {
  if (!is.list(x)) {
    stop("`", section, "` must be a named list.", call. = FALSE)
  }
  nms <- names(x)
  if (length(x) && (is.null(nms) || anyNA(nms) || any(!nzchar(nms)))) {
    stop("Every `", section, "` control must be named.", call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("Duplicate `", section, "` control names are not allowed: ",
      paste(unique(nms[duplicated(nms)]), collapse = ", "), ".", call. = FALSE)
  }
  nms
}

.gl_merge_control_section <- function(defaults, supplied, section) {
  nms <- .gl_control_list_names(supplied, section)
  unknown <- setdiff(nms, names(defaults))
  if (length(unknown)) {
    stop("Unknown `", section, "` control", if (length(unknown) > 1L) "s" else "",
      ": ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  defaults[nms] <- supplied
  defaults
}

.gl_validate_shared_optimizer_controls <- function(x) {
  if (!.gl_is_auto_stop_crit(x$outer_tol)) {
    x$outer_tol <- .gl_validate_stop_crit(x$outer_tol, "outer_tol")
  }
  if (!is.numeric(x$max_outer_iter) || length(x$max_outer_iter) != 1L ||
      is.na(x$max_outer_iter) || x$max_outer_iter < 1 ||
      x$max_outer_iter != as.integer(x$max_outer_iter)) {
    stop("max_outer_iter must be a single positive integer.", call. = FALSE)
  }
  x$max_outer_iter <- as.integer(x$max_outer_iter)
  if (!is.numeric(x$max_elapsed_sec) || length(x$max_elapsed_sec) != 1L ||
      is.na(x$max_elapsed_sec) || x$max_elapsed_sec <= 0) {
    stop("max_elapsed_sec must be a single positive number or Inf.", call. = FALSE)
  }
  if (!is.logical(x$stop_on_convergence) || length(x$stop_on_convergence) != 1L ||
      is.na(x$stop_on_convergence)) {
    stop("stop_on_convergence must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.gl_validate_selected_optimizer_controls <- function(control, method) {
  specified <- attr(control, "specified")
  inactive <- if (identical(method, "RS")) "cg" else "rs"
  if (length(specified[[inactive]])) {
    stop("Controls for `", inactive, "` were supplied, but method = \"", method,
      "\". Remove the wrong-method controls or change `method`.", call. = FALSE)
  }

  if (identical(method, "RS")) {
    x <- control$rs
    if (!.gl_is_auto_stop_crit(x$inner_tol)) x$inner_tol <- .gl_validate_stop_crit(x$inner_tol, "rs$inner_tol")
    x$max_inner_iter <- .gl_positive_integer(x$max_inner_iter, "rs$max_inner_iter")
    if (!(is.numeric(x$max_negative_outer_streak) && length(x$max_negative_outer_streak) == 1L &&
        !is.na(x$max_negative_outer_streak) && x$max_negative_outer_streak > 0)) {
      stop("rs$max_negative_outer_streak must be positive or Inf.", call. = FALSE)
    }
    x$max_negative_outer_streak <- if (is.infinite(x$max_negative_outer_streak)) Inf else as.integer(x$max_negative_outer_streak)
    x$start_step_size <- .gl_positive_number(x$start_step_size, "rs$start_step_size")
    if (!.gl_is_auto_stop_crit(x$step_adjustment)) x$step_adjustment <- .gl_positive_number(x$step_adjustment, "rs$step_adjustment")
    x$max_steps <- .gl_nonnegative_integer(x$max_steps, "rs$max_steps")
    warm <- .gl_normalize_warm_start_controls(x$warm_start_joint, x$warm_start_joint_iter)
    x$warm_start_joint <- warm$warm_start_joint
    x$warm_start_joint_iter <- warm$warm_start_joint_iter
    x$use_backtracking <- .gl_scalar_logical(x$use_backtracking, "rs$use_backtracking")
    x$backtracking_max_halves <- .gl_normalize_backtracking_halves(x$backtracking_max_halves)
    x$update_lambda <- .gl_scalar_logical(x$update_lambda, "rs$update_lambda")
    x$smooth_trust_radius <- .gl_validate_rs_smooth_trust_radius(x$smooth_trust_radius)
    x$discrete_score_method <- match.arg(x$discrete_score_method, c("analytical", "finite"))
    control$rs <- x
  } else {
    x <- control$cg
    x$max_stall <- .gl_positive_integer(x$max_stall, "cg$max_stall")
    x$max_delta <- .gl_positive_number(x$max_delta, "cg$max_delta")
    x$armijo_c1 <- .gl_positive_number(x$armijo_c1, "cg$armijo_c1")
    if (!.gl_is_auto_stop_crit(x$grad_tol)) x$grad_tol <- .gl_validate_stop_crit(x$grad_tol, "cg$grad_tol")
    if (!.gl_is_auto_stop_crit(x$step_tol)) x$step_tol <- .gl_validate_stop_crit(x$step_tol, "cg$step_tol")
    x$update_lambda <- .gl_scalar_logical(x$update_lambda, "cg$update_lambda")
    lambda <- .gl_normalize_cg_lambda_controls(x$lambda_update_every, x$max_lambda_updates, x$raw_loglik_drop_tol)
    x$lambda_update_every <- lambda$cg_lambda_update_every
    x$max_lambda_updates <- lambda$cg_max_lambda_updates
    x$raw_loglik_drop_tol <- lambda$cg_raw_loglik_drop_tol
    x$line_search <- match.arg(x$line_search, c("first", "best"))
    x$max_line_search_evals <- .gl_normalize_cg_line_search_evals(x$max_line_search_evals)
    x$gradient_method <- match.arg(x$gradient_method, c("analytical", "forward", "central"))
    x$zeta_hessian <- match.arg(x$zeta_hessian, c("analytical", "finite"))
    x$hessian_method <- match.arg(x$hessian_method, c("analytical", "finite", "auto"))
    x$use_backtracking <- .gl_scalar_logical(x$use_backtracking, "cg$use_backtracking")
    x$backtracking_max_halves <- .gl_normalize_backtracking_halves(x$backtracking_max_halves)
    control$cg <- x
  }
  control
}

.gl_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be a single positive integer.", call. = FALSE)
  }
  as.integer(x)
}

.gl_nonnegative_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x < 0 || x != as.integer(x)) {
    stop(name, " must be a single non-negative integer.", call. = FALSE)
  }
  as.integer(x)
}

.gl_positive_number <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(name, " must be a single positive finite number.", call. = FALSE)
  }
  as.numeric(x)
}

.gl_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.gl_deprecation_state <- new.env(parent = emptyenv())

.gl_reset_optimizer_control_deprecation <- function() {
  .gl_deprecation_state$optimizer_controls <- FALSE
  invisible()
}

.gl_warn_legacy_optimizer_controls <- function() {
  if (isTRUE(.gl_deprecation_state$optimizer_controls)) return(invisible())
  .gl_deprecation_state$optimizer_controls <- TRUE
  warning(structure(
    list(message = paste0(
      "Flat optimizer arguments are deprecated and will be removed after one release; ",
      "use optimizer_control = gamlss_longitudinal_control(...)."
    ), call = NULL),
    class = c("gamlss.longitudinal_deprecated_optimizer_control", "warning", "condition")
  ))
  invisible()
}

.gl_legacy_optimizer_map <- function() {
  c(
    inner_stop_crit = "rs.inner_tol", outer_stop_crit = "shared.outer_tol",
    start_step_size = "rs.start_step_size", step_adjustment = "rs.step_adjustment",
    max_steps = "rs.max_steps", warm_start_joint = "rs.warm_start_joint",
    warm_start_joint_iter = "rs.warm_start_joint_iter", max_outer_iter = "shared.max_outer_iter",
    max_inner_iter = "rs.max_inner_iter", max_negative_outer_streak = "rs.max_negative_outer_streak",
    max_elapsed_sec = "shared.max_elapsed_sec", cg_max_stall = "cg.max_stall",
    cg_max_delta = "cg.max_delta", cg_armijo_c1 = "cg.armijo_c1",
    cg_grad_tol = "cg.grad_tol", cg_step_tol = "cg.step_tol",
    cg_update_lambda = "cg.update_lambda", cg_lambda_update_every = "cg.lambda_update_every",
    cg_max_lambda_updates = "cg.max_lambda_updates", cg_raw_loglik_drop_tol = "cg.raw_loglik_drop_tol",
    cg_line_search = "cg.line_search", cg_max_line_search_evals = "cg.max_line_search_evals",
    cg_gradient_method = "cg.gradient_method", cg_zeta_hessian = "cg.zeta_hessian",
    cg_hessian_method = "cg.hessian_method", discrete_score_method = "rs.discrete_score_method",
    rs_update_lambda = "rs.update_lambda", rs_smooth_trust_radius = "rs.smooth_trust_radius"
  )
}

.gl_resolve_optimizer_control <- function(method, optimizer_control, legacy_values, legacy_supplied) {
  method <- toupper(as.character(method)[1L])
  if (!method %in% c("RS", "CG")) stop("ERROR: method must be one of 'RS' or 'CG'.", call. = FALSE)
  supplied_new <- !is.null(optimizer_control)
  if (!supplied_new) optimizer_control <- gamlss_longitudinal_control()
  if (!inherits(optimizer_control, "gamlss.longitudinal.control")) {
    stop("optimizer_control must be created by gamlss_longitudinal_control().", call. = FALSE)
  }

  map <- .gl_legacy_optimizer_map()
  legacy_names <- names(legacy_supplied)[legacy_supplied]
  if (length(legacy_names)) .gl_warn_legacy_optimizer_controls()
  specified <- attr(optimizer_control, "specified")
  for (old in legacy_names) {
    target <- unname(map[old])
    if (!length(target) || is.na(target)) next
    bits <- strsplit(target, ".", fixed = TRUE)[[1L]]
    section <- bits[1L]
    key <- bits[2L]
    new_has_key <- if (identical(section, "shared")) isTRUE(specified$shared[[key]]) else key %in% specified[[section]]
    if (supplied_new && new_has_key) {
      stop("Both legacy `", old, "` and optimizer_control's `", target,
        "` specify the same setting.", call. = FALSE)
    }
    optimizer_control[[section]][[key]] <- legacy_values[[old]]
    if (identical(section, "shared")) specified$shared[[key]] <- TRUE else specified[[section]] <- unique(c(specified[[section]], key))
  }

  # The old common backtracking controls apply to the selected method.
  for (old in intersect(c("use_backtracking", "backtracking_max_halves"), legacy_names)) {
    key <- if (identical(old, "use_backtracking")) "use_backtracking" else "backtracking_max_halves"
    section <- tolower(method)
    if (supplied_new && key %in% specified[[section]]) {
      stop("Both legacy `", old, "` and optimizer_control's `", section, ".", key,
        "` specify the same setting.", call. = FALSE)
    }
    optimizer_control[[section]][[key]] <- legacy_values[[old]]
    specified[[section]] <- unique(c(specified[[section]], key))
  }
  attr(optimizer_control, "specified") <- specified
  optimizer_control$shared <- .gl_validate_shared_optimizer_controls(optimizer_control$shared)
  optimizer_control <- .gl_validate_selected_optimizer_controls(optimizer_control, method)
  attr(optimizer_control, "method") <- method
  optimizer_control
}

.gl_optimizer_control_flat <- function(control, method = attr(control, "method")) {
  shared <- control$shared
  if (identical(method, "RS")) {
    rs <- control$rs
    return(c(list(
      inner_stop_crit = rs$inner_tol, outer_stop_crit = shared$outer_tol,
      max_outer_iter = shared$max_outer_iter, max_elapsed_sec = shared$max_elapsed_sec,
      stop_on_convergence = shared$stop_on_convergence
    ), rs))
  }
  cg <- control$cg
  c(list(
    inner_stop_crit = NA, outer_stop_crit = shared$outer_tol,
    max_outer_iter = shared$max_outer_iter, max_elapsed_sec = shared$max_elapsed_sec,
    stop_on_convergence = shared$stop_on_convergence
  ), cg)
}

.gl_effective_optimizer_control <- function(control, optimizer_context) {
  effective <- control
  effective$shared$outer_tol <- optimizer_context$outer_stop_crit
  effective$rs$inner_tol <- optimizer_context$inner_stop_crit
  effective$cg$grad_tol <- optimizer_context$cg_grad_tol_eff
  effective$cg$step_tol <- optimizer_context$cg_step_tol_eff
  attr(effective, "specified") <- NULL
  attr(effective, "method") <- attr(control, "method")
  effective
}
