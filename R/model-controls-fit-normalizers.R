.gl_normalize_backtracking_halves <- function(backtracking_max_halves) {
  if (!is.numeric(backtracking_max_halves) || length(backtracking_max_halves) != 1 ||
      is.na(backtracking_max_halves) || !is.finite(backtracking_max_halves) ||
      backtracking_max_halves < 0 ||
      backtracking_max_halves != as.integer(backtracking_max_halves)) {
    stop("ERROR: backtracking_max_halves must be a single non-negative integer.")
  }
  backtracking_max_halves <- as.integer(backtracking_max_halves)

  backtracking_max_halves
}

.gl_validate_rs_smooth_trust_radius <- function(rs_smooth_trust_radius) {
  if (length(rs_smooth_trust_radius) != 1 || !is.numeric(rs_smooth_trust_radius) ||
    is.na(rs_smooth_trust_radius) || rs_smooth_trust_radius <= 0) {
    stop("ERROR: rs_smooth_trust_radius must be a single positive numeric value or Inf.")
  }

  rs_smooth_trust_radius
}

.gl_normalize_warm_start_controls <- function(warm_start_joint, warm_start_joint_iter) {
  if (!is.logical(warm_start_joint) || length(warm_start_joint) != 1 || is.na(warm_start_joint)) {
    stop("ERROR: warm_start_joint must be TRUE or FALSE.")
  }
  if (!is.numeric(warm_start_joint_iter) || length(warm_start_joint_iter) != 1 || is.na(warm_start_joint_iter)) {
    stop("ERROR: warm_start_joint_iter must be a single non-negative integer.")
  }
  warm_start_joint_iter <- as.integer(warm_start_joint_iter)
  if (warm_start_joint_iter < 0) {
    stop("ERROR: warm_start_joint_iter must be a single non-negative integer.")
  }

  list(
    warm_start_joint = warm_start_joint,
    warm_start_joint_iter = warm_start_joint_iter
  )
}

.gl_normalize_vcov_controls <- function(vcov_method, vcov_numderiv) {
  vcov_method <- match.arg(vcov_method, c("analytical", "numderiv"))
  if (isTRUE(vcov_numderiv)) {
    vcov_method <- "numderiv"
  }
  vcov_numderiv <- identical(vcov_method, "numderiv")

  list(vcov_method = vcov_method, vcov_numderiv = vcov_numderiv)
}

.gl_normalize_cg_line_search_evals <- function(cg_max_line_search_evals) {
  if (length(cg_max_line_search_evals) != 1 || is.null(cg_max_line_search_evals)) {
    stop("ERROR: cg_max_line_search_evals must be a single non-negative integer or NA.")
  }
  if (is.na(cg_max_line_search_evals)) {
    cg_max_line_search_evals <- Inf
  } else {
    if (!is.numeric(cg_max_line_search_evals) ||
        !is.finite(cg_max_line_search_evals) || cg_max_line_search_evals < 0 ||
        cg_max_line_search_evals != as.integer(cg_max_line_search_evals)) {
      stop("ERROR: cg_max_line_search_evals must be a single non-negative integer or NA.")
    }
    cg_max_line_search_evals <- as.integer(cg_max_line_search_evals)
  }

  cg_max_line_search_evals
}

.gl_normalize_cg_lambda_controls <- function(
    cg_lambda_update_every,
    cg_max_lambda_updates,
    cg_raw_loglik_drop_tol) {
  if (!is.numeric(cg_lambda_update_every) || length(cg_lambda_update_every) != 1L ||
      is.na(cg_lambda_update_every) || !is.finite(cg_lambda_update_every) ||
      cg_lambda_update_every < 1L ||
      cg_lambda_update_every != as.integer(cg_lambda_update_every)) {
    stop("cg_lambda_update_every must be a positive integer.")
  }
  cg_lambda_update_every <- as.integer(cg_lambda_update_every)

  if (length(cg_max_lambda_updates) != 1 || is.null(cg_max_lambda_updates)) {
    stop("cg_max_lambda_updates must be a single non-negative integer or NA.")
  }
  if (is.na(cg_max_lambda_updates)) {
    cg_max_lambda_updates <- Inf
  } else {
    if (!is.numeric(cg_max_lambda_updates) || !is.finite(cg_max_lambda_updates) ||
        cg_max_lambda_updates < 0L ||
        cg_max_lambda_updates != as.integer(cg_max_lambda_updates)) {
      stop("cg_max_lambda_updates must be a single non-negative integer or NA.")
    }
    cg_max_lambda_updates <- as.integer(cg_max_lambda_updates)
  }

  if (length(cg_raw_loglik_drop_tol) != 1 || is.null(cg_raw_loglik_drop_tol)) {
    stop("cg_raw_loglik_drop_tol must be a single non-negative numeric value or NA.")
  }
  if (is.na(cg_raw_loglik_drop_tol)) {
    cg_raw_loglik_drop_tol <- Inf
  } else {
    cg_raw_loglik_drop_tol <- as.numeric(cg_raw_loglik_drop_tol)
    if (!is.finite(cg_raw_loglik_drop_tol) || cg_raw_loglik_drop_tol < 0) {
      stop("cg_raw_loglik_drop_tol must be a single non-negative numeric value or NA.")
    }
  }

  list(
    cg_lambda_update_every = cg_lambda_update_every,
    cg_max_lambda_updates = cg_max_lambda_updates,
    cg_raw_loglik_drop_tol = cg_raw_loglik_drop_tol
  )
}

.gl_normalize_cg_fallback_controls <- function(cg_max_stall, cg_max_delta) {
  cg_max_stall <- as.integer(cg_max_stall)
  if (!is.finite(cg_max_stall) || cg_max_stall < 1L) cg_max_stall <- 5L
  if (!is.finite(cg_max_delta) || cg_max_delta <= 0) cg_max_delta <- 0.5

  list(cg_max_stall = cg_max_stall, cg_max_delta = cg_max_delta)
}
