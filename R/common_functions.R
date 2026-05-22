#' @importFrom rlang .data
###########NEW SIMPLIFIED FUNCTIONS

# Null-coalescing operator (base R does not provide one)
`%||%` <- function(a, b) if (!is.null(a)) a else b

utils::globalVariables(c(
  "u1", "u2", "quartile", "tau_emp", "tau_fit", "density", "x_id", "time_pair", "split_group",
  "time", "z", "z_prev", "z_curr", "empirical", "fitted", "threshold", "tail", "probability",
  "emp_copula", "fit_copula", "lag", "cor_z", "n_pairs", "source", "cut_group"
))

.solve_linear_system <- function(A, b = NULL) {
  A <- as.matrix(A)

  if (is.null(b)) {
    b_mat <- diag(nrow(A))
  } else {
    b_mat <- as.matrix(b)
  }

  # Cholesky is fastest and most stable for positive-definite systems.
  chol_A <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(chol_A)) {
    sol <- backsolve(chol_A, forwardsolve(t(chol_A), b_mat))
    if (is.null(b)) {
      return(sol)
    }
    return(sol)
  }

  # General full-rank solve.
  sol <- tryCatch(solve(A, b_mat), error = function(e) NULL)
  if (!is.null(sol)) {
    return(sol)
  }

  # Rank-revealing fallback.
  sol <- tryCatch(qr.solve(A, b_mat), error = function(e) NULL)
  if (!is.null(sol)) {
    return(sol)
  }

  # Last-resort pseudo-inverse for near-singular systems.
  MASS::ginv(A) %*% b_mat
}

.cg_analytical_gradient <- function(
  beta_vec,
  mm_cg,
  eta_out,
  calc_lik,
  margin_dist,
  copula_dist,
  include_dlcopdpar,
  response,
  response_margin,
  response_subject
) {
  eta <- eta_out$eta
  eta_inv <- eta_out$eta_inv
  eta_dr <- eta_out$eta_dr

  grad <- rep(0, length(beta_vec))
  names(grad) <- names(beta_vec)

  margin_par <- intersect(names(mm_cg$x), c("mu", "sigma", "nu", "tau"))
  copula_par <- intersect(names(mm_cg$x), c("theta", "zeta"))

  copula_derivatives <- calc_copula_derivatives(
    eta_inv,
    calc_lik$Fx_1_2,
    copula_dist,
    par1 = calc_lik$copula_par1,
    par2 = calc_lik$copula_par2,
    pair_complete = calc_lik$pair_complete
  )

  margin_score_natural <- list()
  if(length(margin_par) > 0) {
    margin_deriv_subnames <- c("m", "d", "v", "t")
    names(margin_deriv_subnames) <- c("mu", "sigma", "nu", "tau")

    for(par_name in margin_par) {
      d_name <- paste0("dld", margin_deriv_subnames[par_name])
      hit <- grep(paste0("^", d_name, "$"), names(calc_lik$margin_deriv))
      if(length(hit) == 0) {
        hit <- grep(d_name, names(calc_lik$margin_deriv))
      }
      if(length(hit) == 0) {
        margin_score_natural[[par_name]] <- rep(0, length(eta[[par_name]]))
      } else {
        margin_score_natural[[par_name]] <- as.numeric(calc_lik$margin_deriv[[hit[1]]])
      }
    }

    if(isTRUE(include_dlcopdpar)) {
      nd_impact_F <- calc_Fx_derivatives(eta_inv, mm_cg$x, margin_dist, response = response)

      order_margin <- data.frame(time = response_margin, subject = response_subject)
      margin_deriv_1 <- matrix(0, ncol = length(margin_par), nrow = length(response))
      colnames(margin_deriv_1) <- paste0("dld", margin_par)
      for(par_name in margin_par) {
        margin_deriv_1[, paste0("dld", par_name)] <- margin_score_natural[[par_name]]
      }

      margin_components_base <- cbind(
        order_margin,
        response = response,
        margin_p = calc_lik$margin_p,
        margin_d = calc_lik$margin_d,
        margin_deriv_1,
        mu = eta_inv[["mu"]]
      )
      names(margin_components_base)[seq_len(ncol(order_margin))] <- c("time", "subject")

      copula_components <- cbind(
        calc_lik$order_copula,
        row_id1 = calc_lik$copula_row_id1,
        row_id2 = calc_lik$copula_row_id2,
        dcdu1 = copula_derivatives$dcdu1,
        dcdu2 = copula_derivatives$dcdu2,
        copula_d = calc_lik$copula_d
      )

      for(par_name in margin_par) {
        margin_components <- cbind(
          margin_components_base,
          F_nd = nd_impact_F[[par_name]]
        )
        margin_components_Ft_plus <- margin_components
        margin_components_Ft_plus$time <- normalize_lag_time(margin_components_Ft_plus$time)
        margin_plus <- merge(
          margin_components,
          margin_components_Ft_plus,
          by = c("time", "subject"),
          all.x = TRUE
        )
        copula_merged <- merge(
          copula_components,
          margin_plus,
          by.x = c("time1", "subject1"),
          by.y = c("time", "subject"),
          all.x = TRUE
        )
        d1_cop <- calc_deriv_copula_wrt_margin(
          copula_merged,
          margin_par,
          par_name,
          calc_d2 = FALSE
        )[, which(margin_par == par_name)]
        n_score <- length(margin_score_natural[[par_name]])
        if(length(d1_cop) >= n_score) {
          margin_score_natural[[par_name]] <- margin_score_natural[[par_name]] + d1_cop[seq_len(n_score)]
        }
      }
    }
  }

  for(par_name in margin_par) {
    score_eta <- as.numeric(margin_score_natural[[par_name]]) * as.numeric(eta_dr[[par_name]])
    score_eta[!is.finite(score_eta)] <- 0
    X <- as.matrix(mm_cg$x[[par_name]])
    par_grad <- as.numeric(crossprod(X, score_eta))
    x_names <- colnames(X)
    names(par_grad) <- ifelse(
      startsWith(x_names, paste0(par_name, ".")),
      x_names,
      paste(par_name, x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  if("theta" %in% copula_par) {
    n_par <- length(eta[["theta"]])
    d1_full <- rep(0, n_par)
    row_id1 <- calc_lik$copula_row_id1
    if(length(row_id1) > 0) {
      if(n_par == length(response)) {
        par_idx <- row_id1
      } else {
        par_idx <- calc_lik$copula_theta_index_map[row_id1]
      }
      valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
      if(any(valid_idx)) {
        d1_sum <- rowsum(copula_derivatives$dldth[valid_idx], par_idx[valid_idx], reorder = FALSE)
        d1_full[as.integer(rownames(d1_sum))] <- d1_sum[, 1]
      }
    }
    score_eta <- d1_full * as.numeric(eta_dr[["theta"]])
    score_eta[!is.finite(score_eta)] <- 0
    par_grad <- as.numeric(crossprod(as.matrix(mm_cg$x[["theta"]]), score_eta))
    x_names <- colnames(mm_cg$x[["theta"]])
    names(par_grad) <- ifelse(
      startsWith(x_names, "theta."),
      x_names,
      paste("theta", x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  if("zeta" %in% copula_par && "dldz" %in% names(copula_derivatives)) {
    n_par <- length(eta[["zeta"]])
    d1_full <- rep(0, n_par)
    row_id1 <- calc_lik$copula_row_id1
    if(length(row_id1) > 0) {
      if(n_par == length(response)) {
        par_idx <- row_id1
      } else {
        par_idx <- calc_lik$copula_theta_index_map[row_id1]
      }
      valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= n_par
      if(any(valid_idx)) {
        d1_sum <- rowsum(copula_derivatives$dldz[valid_idx], par_idx[valid_idx], reorder = FALSE)
        d1_full[as.integer(rownames(d1_sum))] <- d1_sum[, 1]
      }
    }
    score_eta <- d1_full * as.numeric(eta_dr[["zeta"]])
    score_eta[!is.finite(score_eta)] <- 0
    par_grad <- as.numeric(crossprod(as.matrix(mm_cg$x[["zeta"]]), score_eta))
    x_names <- colnames(mm_cg$x[["zeta"]])
    names(par_grad) <- ifelse(
      startsWith(x_names, "zeta."),
      x_names,
      paste("zeta", x_names, sep = ".")
    )
    idx <- match(names(par_grad), names(grad))
    valid <- !is.na(idx)
    grad[idx[valid]] <- par_grad[valid]
  }

  grad
}

#' Fit a longitudinal joint regression model
#'
#' This function fits a longitudinal model to a dataset with gamlss margins
#' and copula fit to dependence. Any linear or factor covariates can be fit
#' to any parameters of the copula or margin distributions. The model is fit
#' using RS() optimisation and the joint likelihood by default. Select
#' use dlcopdpar=FALSE to fit separately optimised models for the margin and
#' copula likelihoods which can be quicker with a slight loss to overall fit.
#'
#' @param margin_dist Marginal distribution specified as a gamlss family object,
#' e.g. GA(), NO(), PO(), NBI(), etc.
#' @param copula_dist Copula distribution code, one of "N", "C", "F", "G", "J", or "t".
#' @param mu.formula Formula for the mean parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param include_dlcopdpar Include the derivative of the copula likelihood with respect
#' to the margin parameters in the joint likelihood.
#' @param check_dlcopdpar_gradient If `TRUE`, run an optional finite-difference
#' diagnostic for the margin score contribution when `include_dlcopdpar = TRUE`.
#' @param inner_stop_crit Stopping criterion for the inner loop. If `NA` or
#' `NULL`, an automatic data-adaptive value is used.
#' @param outer_stop_crit Stopping criterion for the outer loop. If `NA` or
#' `NULL`, an automatic data-adaptive value is used.
#' @param start_step_size Initial step size for the backfitting algorithm
#' @param step_adjustment Step size adjustment factor
#' @param max_steps Maximum number of times for reducing the step size
#' @param start_from Starting values for the parameters if needed
#' @param warm_start_joint Logical; if `TRUE` (default), RS joint fits started
#' without explicit `start_from` first run a short separate RS stabilisation
#' phase and use those coefficients as the joint starting values.
#' @param warm_start_joint_iter Integer; number of separate RS outer iterations
#' used for the default joint warm start.
#' @param verbose Level of output to the console 3 = ALL, 0 = Minimal
#' @param plot_results Plot the results of the optimisation
#' @param true_val True values for the parameters if known for plotting
#' @param method Optimisation method to use, RS() is the default
#' @param max_outer_iter Maximum number of outer iterations
#' @param max_inner_iter Maximum number of inner iterations
#' @param max_negative_outer_streak Maximum number of consecutive negative outer
#' log-likelihood changes allowed before stopping.
#' @param max_elapsed_sec Optional maximum elapsed fitting time in seconds.
#' If finite, the optimiser stops with an error once this budget is exceeded.
#' @param use_backtracking Logical; if `TRUE` (default), apply step-halving
#' backtracking to reject downhill inner updates.
#' @param backtracking_max_halves Integer; maximum number of consecutive
#' step halvings attempted after a rejected update before taking no step.
#' @param cg_max_stall Integer; for `method = "CG"` only. Maximum number of
#' consecutive outer iterations where no improving step is found before CG stops.
#' @param cg_max_delta Numeric; for `method = "CG"` only. Maximum absolute
#' coefficient step size used to limit Newton/trust-region updates.
#' @param cg_armijo_c1 Numeric; for `method = "CG"` only. Minimum improvement
#' threshold used by the line-search acceptance rule.
#' @param cg_grad_tol Numeric; for `method = "CG"` only. Penalized-gradient
#' infinity-norm convergence tolerance. If `NA`, selected from `outer_stop_crit`.
#' @param cg_step_tol Numeric; for `method = "CG"` only. Accepted-step L2
#' convergence tolerance. If `NA`, selected from `outer_stop_crit`.
#' @param cg_update_lambda Logical; for `method = "CG"` only. If `TRUE`, update
#' smoother penalties during CG iterations.
#' @param cg_lambda_update_every Integer; for `method = "CG"` only. When
#' `cg_update_lambda = TRUE`, update each smoother's lambda every this many
#' outer iterations. Use `1` to update every CG iteration.
#' @param cg_max_lambda_updates Integer; for `method = "CG"` only. Maximum
#' number of smoother penalty update rounds. Use `NA` for no cap.
#' @param cg_raw_loglik_drop_tol Numeric; for `method = "CG"` only. Stop CG as
#' not converged if the raw joint log-likelihood drops this far below the best
#' raw joint log-likelihood seen after at least one lambda update. Use `NA` to
#' disable.
#' @param cg_line_search Character; for `method = "CG"` only. `"best"` evaluates
#' candidate steps up to `cg_max_line_search_evals` before taking the largest
#' improvement, while `"first"` accepts the first improving candidate step.
#' @param cg_max_line_search_evals Integer; for `method = "CG"` only. Optional
#' cap on the number of candidate likelihood evaluations per outer iteration.
#' @param cg_gradient_method Character; for `method = "CG"` only.
#' `"analytical"` uses the same score components as RS, `"forward"` uses
#' one-sided finite differences, and `"central"` uses two-sided finite
#' differences.
#' @param cg_zeta_hessian Character; for `method = "CG"` only. `"analytical"`
#' uses the analytical Hessian for the zeta block, while `"finite"` replaces
#' the zeta-zeta block with central finite differences of the raw joint
#' log-likelihood.
#' @param compute_vcov Logical; if `TRUE` (default), compute and store the
#' model variance-covariance output at the end of fitting.
#' @param vcov_method Character; fit-time vcov method when `compute_vcov = TRUE`.
#' One of `"numderiv"` or `"analytical"`.
#' @param vcov_numderiv Logical; passed to `vcov.gamlss.longitudinal()` when
#' `compute_vcov = TRUE`.
#' @param use_Rcpp Use Rcpp for matrix operations
#'
#' @export
gamlss.longitudinal=function(dataset,
                        margin_dist,
                        copula_dist,
                        time_var=NA,
                        subject_var=NA,
                        mu.formula = ("response ~ 1"),
                        sigma.formula = ("~ 1"),
                        nu.formula = ("~ 1"),
                        tau.formula = ("~ 1"),
                        theta.formula=("~ 1"),
                        zeta.formula=("~ 1"),
                        include_dlcopdpar=TRUE,
                        check_dlcopdpar_gradient=FALSE,
                        inner_stop_crit=NA,
                        outer_stop_crit=NA,
                        start_step_size=.5,
                        step_adjustment=NA,
                        max_steps=5,
                        start_from=NA,
                        warm_start_joint=TRUE,
                        warm_start_joint_iter=5,
                        verbose=1,
                        plot_results=FALSE,
                        true_val=NA,
                        method="RS",
                        max_outer_iter=100,
                        max_inner_iter=100,
                        max_negative_outer_streak=10,
                        max_elapsed_sec=Inf,
                        use_backtracking=TRUE,
                        backtracking_max_halves=50,
                        cg_max_stall=5,
                        cg_max_delta=0.5,
                        cg_armijo_c1=1e-4,
                        cg_grad_tol=NA,
                        cg_step_tol=NA,
                        cg_update_lambda=TRUE,
                        cg_lambda_update_every=10,
                        cg_max_lambda_updates=NA,
                        cg_raw_loglik_drop_tol=10,
                        cg_line_search="best",
                        cg_max_line_search_evals=60,
                        cg_gradient_method="forward",
                        cg_zeta_hessian="analytical",
                        compute_vcov=TRUE,
                        vcov_method=c("analytical","numderiv"),
                        vcov_numderiv=FALSE,
                        use_Rcpp=FALSE,
                        lambda_start=NA,
                        lambda_penalty_K=2
                      )
{
  fit_start_time <- Sys.time()
  check_elapsed_budget <- function(stage = "optimisation") {
    if (is.finite(max_elapsed_sec) && max_elapsed_sec > 0) {
      elapsed <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
      if (elapsed > max_elapsed_sec) {
        stop(
          sprintf(
            "Model exceeded max_elapsed_sec during %s (elapsed %.1f sec > %.1f sec).",
            stage, elapsed, max_elapsed_sec
          ),
          call. = FALSE
        )
      }
    }
    invisible(TRUE)
  }

  if (!is.numeric(backtracking_max_halves) || length(backtracking_max_halves) != 1 || is.na(backtracking_max_halves)) {
    stop("ERROR: backtracking_max_halves must be a single non-negative integer.")
  }
  backtracking_max_halves <- as.integer(backtracking_max_halves)
  if (backtracking_max_halves < 0) {
    stop("ERROR: backtracking_max_halves must be a single non-negative integer.")
  }

  method <- toupper(as.character(method)[1])
  cg_line_search <- match.arg(as.character(cg_line_search)[1], c("first", "best"))
  cg_gradient_method <- match.arg(as.character(cg_gradient_method)[1], c("analytical", "forward", "central"))
  cg_zeta_hessian <- match.arg(as.character(cg_zeta_hessian)[1], c("analytical", "finite"))
  if (length(cg_max_line_search_evals) != 1 || is.null(cg_max_line_search_evals)) {
    stop("ERROR: cg_max_line_search_evals must be a single non-negative integer or NA.")
  }
  if (is.na(cg_max_line_search_evals)) {
    cg_max_line_search_evals <- Inf
  } else {
    cg_max_line_search_evals <- as.integer(cg_max_line_search_evals)
    if (!is.finite(cg_max_line_search_evals) || cg_max_line_search_evals < 0) {
      stop("ERROR: cg_max_line_search_evals must be a single non-negative integer or NA.")
    }
  }
  if(!method %in% c("RS", "CG")) {
    stop("ERROR: method must be one of 'RS' or 'CG'.")
  }
  user_supplied_start <- !all(is.na(start_from))
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
  vcov_method <- match.arg(vcov_method)
  if (isTRUE(vcov_numderiv)) {
    vcov_method <- "numderiv"
  }
  vcov_numderiv <- identical(vcov_method, "numderiv")
  cg_lambda_update_every <- as.integer(cg_lambda_update_every)
  if(!is.finite(cg_lambda_update_every) || cg_lambda_update_every < 1L) {
    stop("cg_lambda_update_every must be a positive integer.")
  }
  if (length(cg_max_lambda_updates) != 1 || is.null(cg_max_lambda_updates)) {
    stop("cg_max_lambda_updates must be a single non-negative integer or NA.")
  }
  if (is.na(cg_max_lambda_updates)) {
    cg_max_lambda_updates <- Inf
  } else {
    cg_max_lambda_updates <- as.integer(cg_max_lambda_updates)
    if (!is.finite(cg_max_lambda_updates) || cg_max_lambda_updates < 0L) {
      stop("cg_max_lambda_updates must be a single non-negative integer or NA.")
    }
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
  cg_max_stall <- as.integer(cg_max_stall)
  if(!is.finite(cg_max_stall) || cg_max_stall < 1L) cg_max_stall <- 5L
  if(!is.finite(cg_max_delta) || cg_max_delta <= 0) cg_max_delta <- 0.5

  ##################### DATA CHECKS AND VALIDATION #####################

  # Save original dataset
  dataset_original <- dataset

  # Force plain data.frame (safe for tibble/data.table too)
  if (!is.data.frame(dataset)) {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  } else {
    dataset <- as.data.frame(dataset, stringsAsFactors = FALSE)
  }

  # Validate and prepare input data
  if(all(is.na(time_var)) || all(is.na(subject_var))) {
    stop("ERROR: Required input variables not specified.\n",
         "Please specify:\n",
         "  - time_var: column name for time/margin variable (e.g., 'time')\n",
         "  - subject_var: column name for subject ID variable (e.g., 'subject')\n",
         "Example: gamlss.longitudinal(..., time_var='time', subject_var='subject')")
  }

  # Validate dataset contains required columns
  if(!time_var %in% colnames(dataset)) {
    stop("ERROR: time_var='", time_var, "' not found in dataset.\n",
         "Available columns: ", paste(colnames(dataset), collapse=", "))
  }
  if(!subject_var %in% colnames(dataset)) {
    stop("ERROR: subject_var='", subject_var, "' not found in dataset.\n",
         "Available columns: ", paste(colnames(dataset), collapse=", "))
  }

  # Extract response variable name from mu formula
  mu_formula_obj <- as.formula(mu.formula)
  response_var <- all.vars(mu_formula_obj)[1]
  if (verbose > 1) print(paste("Identified response variable:", response_var))

  # Validate response variable exists
  if(!response_var %in% names(dataset)) {
    stop("ERROR: response variable '", response_var, "' not found in dataset.\n",
         "Available columns: ", paste(names(dataset), collapse=", "))
  }

  # Rename columns for internal use
  names(dataset)[names(dataset) == time_var] <- "time"
  names(dataset)[names(dataset) == subject_var] <- "subject"
  names(dataset)[names(dataset) == response_var] <- "response"

  # Preserve the user-facing time covariate (including factor type) for formulas,
  # while keeping an internal numeric time index for optimisation logic.
  dataset$time_covariate <- dataset$time
  time_covariate_is_factor <- is.factor(dataset$time_covariate)
  time_covariate_levels <- if (time_covariate_is_factor) levels(dataset$time_covariate) else NULL
  time_covariate_ordered <- if (time_covariate_is_factor) is.ordered(dataset$time_covariate) else FALSE

  if (is.factor(dataset$time_covariate)) {
    time_chr <- as.character(dataset$time_covariate)
    dataset$time <- match(time_chr, time_covariate_levels)
    if (anyNA(dataset$time)) {
      stop("ERROR: Failed to map factor time levels to internal numeric time index.")
    }
    dataset$time_covariate <- factor(time_chr, levels = time_covariate_levels, ordered = time_covariate_ordered)
    if (time_covariate_ordered) {
      time_contr <- contr.treatment(length(time_covariate_levels))
      if (length(time_covariate_levels) > 1) {
        colnames(time_contr) <- time_covariate_levels[-1]
      }
      contrasts(dataset$time_covariate) <- time_contr
    }
  } else if (is.numeric(dataset$time_covariate) || is.integer(dataset$time_covariate)) {
    dataset$time <- as.numeric(dataset$time_covariate)
  } else if (is.character(dataset$time_covariate)) {
    time_numeric <- suppressWarnings(as.numeric(dataset$time_covariate))
    if (anyNA(time_numeric)) {
      stop("ERROR: time must be numeric-like unless supplied as factor.\n",
           "If time is categorical for formulas/interactions, convert it to factor before fitting.")
    }
    dataset$time <- time_numeric
    dataset$time_covariate <- time_numeric
  } else {
    stop("ERROR: Unsupported time variable type: ", class(dataset$time_covariate)[1],
         ". Use numeric/integer, numeric-like character, or factor.")
  }

  if (is.factor(dataset$subject)) {
    dataset$subject <- as.character(dataset$subject)
  }

  if (any(is.na(dataset$time)) || any(is.na(dataset$subject))) {
    stop("ERROR: time and subject variables cannot contain NA values.")
  }

  # Robust formula normalizer:
  # - accepts formula objects
  # - accepts strings without "~" (e.g. "time + s(age)")
  # - for mu, adds response on LHS if missing
  normalize_formula <- function(fml, response_name = "response", require_lhs = FALSE) {
    if (inherits(fml, "formula")) {
      return(fml)
    }
    if (!is.character(fml) || length(fml) != 1 || is.na(fml) || nchar(trimws(fml)) == 0) {
      stop("ERROR: Invalid formula input: ", deparse(fml))
    }

    txt <- trimws(fml)

    if (!grepl("~", txt, fixed = TRUE)) {
      txt <- if (require_lhs) paste0(response_name, " ~ ", txt) else paste0("~ ", txt)
    } else if (require_lhs) {
      parts <- strsplit(txt, "~", fixed = TRUE)[[1]]
      lhs <- trimws(parts[1])
      rhs <- trimws(parts[2])
      if (nchar(lhs) == 0) txt <- paste0(response_name, " ~ ", rhs)
    }

    as.formula(txt, env = parent.frame())
  }

  # Variable-name translation from user names -> internal names
  translate_formula_vars <- function(fml, var_map, response_name = "response", require_lhs = FALSE) {
    f_obj <- normalize_formula(fml, response_name = response_name, require_lhs = require_lhs)
    f_txt <- paste(deparse(f_obj), collapse = " ")

    for (old_name in names(var_map)) {
      new_name <- var_map[[old_name]]
      if (!identical(old_name, new_name)) {
        old_esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", old_name)
        f_txt <- gsub(paste0("`", old_esc, "`"), new_name, f_txt, perl = TRUE)
        f_txt <- gsub(paste0("\\b", old_esc, "\\b"), new_name, f_txt, perl = TRUE)
      }
    }

    as.formula(f_txt, env = environment(f_obj))
  }

  var_map <- c()
  var_map[[time_var]] <- "time"
  var_map[[subject_var]] <- "subject"
  var_map[[response_var]] <- "response"

  # For formula parsing, map user time variable to preserved covariate column.
  formula_var_map <- var_map
  formula_var_map[[time_var]] <- "time_covariate"

  mu.formula.int    <- translate_formula_vars(mu.formula,    formula_var_map, response_name = "response", require_lhs = TRUE)
  sigma.formula.int <- translate_formula_vars(sigma.formula, formula_var_map, response_name = "response", require_lhs = FALSE)
  nu.formula.int    <- translate_formula_vars(nu.formula,    formula_var_map, response_name = "response", require_lhs = FALSE)
  tau.formula.int   <- translate_formula_vars(tau.formula,   formula_var_map, response_name = "response", require_lhs = FALSE)
  theta.formula.int <- translate_formula_vars(theta.formula, formula_var_map, response_name = "response", require_lhs = FALSE)
  zeta.formula.int  <- translate_formula_vars(zeta.formula,  formula_var_map, response_name = "response", require_lhs = FALSE)

  if(verbose > 1) {
    cat("Input validation successful.\n")
    cat("Data dimensions:", nrow(dataset), "x", ncol(dataset), "\n")
    cat("Response variable:", response_var, "-> renamed to 'response'\n")
    cat("Time variable:", time_var, "-> internal index 'time' and covariate 'time_covariate'\n")
    cat("Subject variable:", subject_var, "-> renamed to 'subject'\n")
    cat("Time points:", length(unique(dataset$time)), "\n")
    cat("Subjects:", length(unique(dataset$subject)), "\n")
  }

  # Validate that all subject/time combinations are unique
  subject_time_combo <- paste(dataset$subject, dataset$time, sep="_")
  if(length(subject_time_combo) != length(unique(subject_time_combo))) {
    duplicate_combos <- subject_time_combo[duplicated(subject_time_combo)]
    stop("ERROR: Duplicate subject/time combinations found.\n",
         "Each subject must have exactly one observation per time point.\n",
         "Duplicate combinations (first 10): ",
         paste(unique(duplicate_combos)[1:min(10, length(unique(duplicate_combos)))], collapse=", "))
  }

  if(verbose > 1) {
    cat("Subject/time uniqueness check passed.\n")
    cat("Unique subject/time combinations:", length(unique(subject_time_combo)), "\n\n")
  }

  # One-to-one map from internal time index back to preserved covariate values.
  time_lookup <- dataset[!duplicated(dataset$time), c("time", "time_covariate"), drop = FALSE]
  time_lookup <- time_lookup[order(time_lookup$time), , drop = FALSE]

  # Expand to full subject x time grid so structurally missing combinations
  # are represented explicitly as NA rows.
  observed_n <- nrow(dataset)
  full_grid <- expand.grid(
    subject = sort(unique(dataset$subject)),
    time = sort(unique(dataset$time)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  dataset <- merge(full_grid, dataset, by = c("subject", "time"), all.x = TRUE, sort = FALSE)
  dataset$time_covariate <- time_lookup$time_covariate[match(dataset$time, time_lookup$time)]
  if (time_covariate_is_factor) {
    dataset$time_covariate <- factor(as.character(dataset$time_covariate),
                                     levels = time_covariate_levels,
                                     ordered = time_covariate_ordered)
    if (time_covariate_ordered) {
      time_contr <- contr.treatment(length(time_covariate_levels))
      if (length(time_covariate_levels) > 1) {
        colnames(time_contr) <- time_covariate_levels[-1]
      }
      contrasts(dataset$time_covariate) <- time_contr
    }
  }
  dataset <- dataset[order(dataset$subject, dataset$time), , drop = FALSE]
  rownames(dataset) <- NULL

  inserted_n <- nrow(dataset) - observed_n
  if (verbose > 0 && inserted_n > 0) {
    cat("Inserted", inserted_n, "missing subject/time rows as NA entries.\n\n")
  }

  # Missingness summary by time and consecutive time pairs.
  time_levels <- sort(unique(dataset$time))
  n_time_levels <- length(time_levels)

  miss_by_time <- do.call(rbind, lapply(time_levels, function(ti) {
    idx <- dataset$time == ti
    n_total <- sum(idx)
    n_na <- sum(is.na(dataset$response[idx]))
    c(time = ti, n_total = n_total, n_na_response = n_na, n_observed_response = n_total - n_na)
  }))
  miss_by_time <- as.data.frame(miss_by_time)
  rownames(miss_by_time) <- NULL

  pair_summary <- data.frame(
    time1 = numeric(0),
    time2 = numeric(0),
    complete_pairs = integer(0),
    total_pairs = integer(0),
    total_observations = integer(0)
  )

  if (n_time_levels > 1) {
    for (i in seq_len(n_time_levels - 1)) {
      t1 <- time_levels[i]
      t2 <- time_levels[i + 1]
      d1 <- dataset[dataset$time == t1, c("subject", "response")]
      d2 <- dataset[dataset$time == t2, c("subject", "response")]
      names(d1) <- c("subject", "response_t1")
      names(d2) <- c("subject", "response_t2")
      merged_pair <- merge(d1, d2, by = "subject", all = FALSE)

      total_pairs <- nrow(merged_pair)
      complete_pairs <- sum(!is.na(merged_pair$response_t1) & !is.na(merged_pair$response_t2))

      pair_summary <- rbind(
        pair_summary,
        data.frame(
          time1 = t1,
          time2 = t2,
          complete_pairs = complete_pairs,
          total_pairs = total_pairs,
          total_observations = nrow(dataset)
        )
      )
    }
  }

  if (verbose > 0) {
    cat("Missingness Summary (by time):\n")
    print(miss_by_time)
    cat("\nConsecutive Pair Completeness:\n")
    if (nrow(pair_summary) > 0) {
      print(pair_summary)
    } else {
      cat("No consecutive time pairs available.\n")
    }
    cat("\n")
  }

  # Hard stop if any margin is 100% missing.
  margin_all_missing <- miss_by_time$n_observed_response == 0
  if (any(margin_all_missing)) {
    bad_times <- miss_by_time$time[margin_all_missing]
    stop(
      "ERROR: 100% missing response values detected for margin time point(s): ",
      paste(bad_times, collapse = ", "),
      "\nModel fitting stopped because at least one margin has no observed outcomes."
    )
  }

  # Hard stop if any consecutive copula pair has 0 complete pairs while pairs exist.
  if (nrow(pair_summary) > 0) {
    pair_all_missing <- pair_summary$total_pairs > 0 & pair_summary$complete_pairs == 0
    if (any(pair_all_missing)) {
      bad_pairs <- apply(pair_summary[pair_all_missing, c("time1", "time2"), drop = FALSE], 1, function(x) {
        paste0("(", x[1], ",", x[2], ")")
      })
      stop(
        "ERROR: 100% missing complete copula pairs detected for consecutive time pair(s): ",
        paste(bad_pairs, collapse = ", "),
        "\nModel fitting stopped because at least one copula pair contributes no complete observations."
      )
    }
  }

  ##################### END OF DATA CHECKS AND VALIDATION #####################

  ##################### MODEL SETUP #####################

  #Setup model matrix from given formulas
  copula_link=get_copula_dist(copula_dist)$copula_link
  mm=suppressWarnings(create_model_matrices(
    mu.formula.int,
    sigma.formula.int,
    nu.formula.int,
    tau.formula.int,
    theta.formula.int,
    zeta.formula.int,
    margin.family = margin_dist,
    copula.family = copula_dist,
    copula.link = copula_link,
    dataset = dataset
  ))

  warm_start_info <- list(
    used = FALSE,
    outer_iter = 0L,
    include_dlcopdpar = FALSE,
    log_lik = NULL
  )
  warm_start_par_s <- NULL

  if (
    method == "RS" &&
    isTRUE(include_dlcopdpar) &&
    isTRUE(warm_start_joint) &&
    warm_start_joint_iter > 0L &&
    !isTRUE(user_supplied_start)
  ) {
    if (verbose > 0) {
      cat(
        "\nRunning separate RS warm-start phase for ",
        warm_start_joint_iter,
        " outer iteration(s) before joint RS fit...\n",
        sep = ""
      )
    }

    warm_fit <- NULL
    warm_output <- NULL
    warm_err <- NULL
    tryCatch({
      warm_output <- capture.output({
        warm_fit <- gamlss.longitudinal(
          dataset = dataset_original,
          margin_dist = margin_dist,
          copula_dist = copula_dist,
          time_var = time_var,
          subject_var = subject_var,
          mu.formula = mu.formula,
          sigma.formula = sigma.formula,
          nu.formula = nu.formula,
          tau.formula = tau.formula,
          theta.formula = theta.formula,
          zeta.formula = zeta.formula,
          include_dlcopdpar = FALSE,
          check_dlcopdpar_gradient = FALSE,
          inner_stop_crit = inner_stop_crit,
          outer_stop_crit = outer_stop_crit,
          start_step_size = start_step_size,
          step_adjustment = step_adjustment,
          max_steps = max_steps,
          start_from = NA,
          warm_start_joint = FALSE,
          warm_start_joint_iter = 0L,
          verbose = 0,
          plot_results = FALSE,
          true_val = true_val,
          method = method,
          max_outer_iter = warm_start_joint_iter,
          max_inner_iter = max_inner_iter,
          max_negative_outer_streak = max_negative_outer_streak,
          max_elapsed_sec = max_elapsed_sec,
          use_backtracking = use_backtracking,
          backtracking_max_halves = backtracking_max_halves,
          cg_max_stall = cg_max_stall,
          cg_max_delta = cg_max_delta,
          cg_armijo_c1 = cg_armijo_c1,
          cg_grad_tol = cg_grad_tol,
          cg_step_tol = cg_step_tol,
          cg_update_lambda = cg_update_lambda,
          cg_lambda_update_every = cg_lambda_update_every,
          cg_line_search = cg_line_search,
          cg_max_line_search_evals = cg_max_line_search_evals,
          cg_gradient_method = cg_gradient_method,
          compute_vcov = FALSE,
          vcov_method = vcov_method,
          vcov_numderiv = vcov_numderiv,
          use_Rcpp = use_Rcpp,
          lambda_start = lambda_start,
          lambda_penalty_K = lambda_penalty_K
        )
      }, type = "output")
    }, error = function(e) {
      warm_err <<- e
    })

    if (!is.null(warm_err)) {
      stop(
        "Separate RS warm-start phase failed: ",
        conditionMessage(warm_err),
        "\nSet warm_start_joint = FALSE to force a cold-start joint fit.",
        call. = FALSE
      )
    }
    if (is.null(warm_fit) || is.null(warm_fit$par)) {
      stop(
        "Separate RS warm-start phase did not return coefficient starting values.\n",
        "Set warm_start_joint = FALSE to force a cold-start joint fit.",
        call. = FALSE
      )
    }

    start_from <- warm_fit$par
    warm_start_par_s <- warm_fit$par_s
    warm_start_info <- list(
      used = TRUE,
      outer_iter = warm_start_joint_iter,
      include_dlcopdpar = FALSE,
      log_lik = warm_fit$calc_lik_out_end$log_lik,
      carries_smooth = !is.null(warm_start_par_s) && any(vapply(warm_start_par_s, length, integer(1L)) > 0L),
      captured_output = warm_output
    )

    if (verbose > 1 && length(warm_output) > 0) {
      cat(paste(warm_output, collapse = "\n"), "\n")
    }
  }

  if (length(start_step_size) != 1 || !is.numeric(start_step_size) ||
      !is.finite(start_step_size) || start_step_size <= 0) {
    stop("ERROR: start_step_size must be a single positive finite numeric value.")
  }
  if (length(max_steps) != 1 || !is.numeric(max_steps) || is.na(max_steps)) {
    stop("ERROR: max_steps must be a single non-negative integer.")
  }
  max_steps <- as.integer(max_steps)
  if (max_steps < 0) {
    stop("ERROR: max_steps must be a single non-negative integer.")
  }
  if (length(step_adjustment) != 1 || is.null(step_adjustment)) {
    stop("ERROR: step_adjustment must be a single positive numeric value, or NA for the method-specific default.")
  }
  step_adjustment <- as.numeric(step_adjustment)
  if (is.na(step_adjustment)) {
    rs_joint_step_adjustment_default <- 1
    rs_separate_step_adjustment_default <- 1
    step_adjustment <- if (method == "RS" && isTRUE(include_dlcopdpar)) {
      rs_joint_step_adjustment_default
    } else if (method == "RS") {
      rs_separate_step_adjustment_default
    } else {
      1
    }
    if (verbose > 0) {
      cat(
        "\nUsing automatic step_adjustment=",
        signif(step_adjustment, 4),
        " for ",
        if (method == "RS" && isTRUE(include_dlcopdpar)) "joint RS" else if (method == "RS") "separate RS" else method,
        ".\n",
        sep = ""
      )
    }
  } else if (!is.finite(step_adjustment) || step_adjustment <= 0) {
    stop("ERROR: step_adjustment must be a single positive numeric value, or NA for the method-specific default.")
  }

  #Create vector of starting covariate values, currently starting at zero before first fit with the intercept as the mean
  if(all(is.na(start_from))) {
    par_eta=get_starting_values(copula_dist,margin_dist,dataset=dataset,eta_transform=TRUE)
    par_cov=as.numeric(vector())
    for (par_name in names(mm$x)) {
      par_cov_single=as.numeric(vector(length=length(colnames(mm$x[[par_name]]))))
      names(par_cov_single)=paste(par_name,colnames(mm$x[[par_name]]),sep=".")
      par_cov_single[1]=par_eta[par_name]
      if(length(par_cov_single)>1) {
        par_cov_single[2:length(par_cov_single)]=0
      }
      par_cov=c(par_cov,par_cov_single)
    }
  } else {
    par_cov=start_from
  }
  par_s=list()
  df_s=list()
  lambda_s=list()
  #names(par_s)=names(df_s)=names(lambda_s)=names(mm$x)
  for (par_name in names(mm$x)) {
    par_s[[par_name]]=list()
    df_s[[par_name]]=list()
    lambda_s[[par_name]]=list()
    for (s_name in names(mm$s[[par_name]])) {
      B=mm$s[[par_name]][[s_name]]
      par_s[[par_name]][[s_name]]=c(par_s[[par_name]][[s_name]],rep(0,ncol(B)))
      names(par_s[[par_name]][[s_name]])=paste(par_name,s_name,1:ncol(B),sep=".")
      df_s[[par_name]][[s_name]]=0
      # Data-adaptive starting lambda: tr(B'B) / tr(S) balances the penalty and
      # data terms regardless of n, k, or response scale. Used when the user has
      # not supplied an explicit lambda_start (i.e. lambda_start = NA).
      S_init <- attr(B, "penalty")
      if (is.na(lambda_start)) {
        if (!is.null(S_init) && is.matrix(S_init) && sum(diag(S_init)) > 0) {
          lambda_s[[par_name]][[s_name]] <- sum(diag(t(B) %*% B)) / sum(diag(S_init))
        } else {
          lambda_s[[par_name]][[s_name]] <- 10  # fallback if no penalty stored
        }
      } else {
        lambda_s[[par_name]][[s_name]] <- lambda_start
      }
      names(df_s[[par_name]][[s_name]])=names(lambda_s[[par_name]][[s_name]])=s_name
   }
  }
  if(!is.null(warm_start_par_s)) {
    for (par_name in intersect(names(par_s), names(warm_start_par_s))) {
      if(length(par_s[[par_name]]) == 0 || length(warm_start_par_s[[par_name]]) == 0) next
      for (s_name in intersect(names(par_s[[par_name]]), names(warm_start_par_s[[par_name]]))) {
        warm_beta <- warm_start_par_s[[par_name]][[s_name]]
        if(length(warm_beta) == length(par_s[[par_name]][[s_name]])) {
          par_s[[par_name]][[s_name]] <- warm_beta
        }
      }
    }
  }
  #Starting parameters for fixed parameters: par_cov

  rs_design_cache <- setNames(vector("list", length(names(mm$x))), names(mm$x))
  for (pn in names(mm$x)) {
    X_fixed <- as.matrix(mm$x[[pn]])
    fixed_names <- paste(pn, colnames(mm$x[[pn]]), sep = ".")
    X_parts <- list(X_fixed)
    smooth_penalty_meta <- list()

    if (length(mm$s[[pn]]) > 0) {
      start_idx <- ncol(X_fixed) + 1L
      for (s_name in names(mm$s[[pn]])) {
        B <- as.matrix(mm$s[[pn]][[s_name]])
        smooth_names <- names(par_s[[pn]][[s_name]])
        colnames(B) <- smooth_names
        X_parts[[length(X_parts) + 1L]] <- B

        n_B <- ncol(B)
        idx <- start_idx:(start_idx + n_B - 1L)
        pen_attr <- attr(mm$s[[pn]][[s_name]], "penalty")
        if (!is.null(pen_attr) && is.matrix(pen_attr) &&
            nrow(pen_attr) == n_B && ncol(pen_attr) == n_B) {
          S_base <- pen_attr
        } else {
          D <- diff(diag(n_B), differences = 2)
          S_base <- t(D) %*% D
        }
        smooth_penalty_meta[[s_name]] <- list(idx = idx, B = B, S_base = S_base)
        start_idx <- start_idx + n_B
      }
    }

    X_combined <- do.call(cbind, X_parts)
    colnames(X_combined)[seq_along(fixed_names)] <- fixed_names
    rs_design_cache[[pn]] <- list(
      X = X_combined,
      fixed_names = fixed_names,
      smooth_penalty_meta = smooth_penalty_meta
    )
  }

  #Parameters used in optimisation loops
  first_outer_run=TRUE
  outer_log_lik_change=outer_start_log_lik=outer_end_log_lik=0
  log_lik_history=matrix(ncol=3,nrow=0)
  par_history=matrix(ncol=length(par_cov),nrow=0); colnames(par_history)=names(par_cov)
  cg_stop_reason <- NA_character_
  cg_last_grad_inf <- NA_real_
  cg_last_step_l2 <- NA_real_
  cg_best_raw_loglik <- -Inf
  cg_best_iteration <- NA_integer_
  cg_raw_loglik_drop_from_best <- NA_real_
  rs_block_trace <- list()
  outer_run_counter=1; outer_only_run_counter=1
  outer_negative_streak=0
  step_size=start_step_size
  weights_final=list()
  pair_cache=build_copula_pair_cache(
    response=dataset$response,
    response_margin=dataset$time,
    response_subject=dataset$subject
  )
  margin_eval_cache=.build_margin_eval_cache(margin_dist, calc_d2 = FALSE)
  rs_calc_eta <- function(par_cov_current, par_s_current, update_only = NULL, eta_out_current = NULL) {
    if (isTRUE(getOption("gamlss.longitudinal.fast_rs_eta", TRUE))) {
      .calc_eta_rs_cached(
        rs_design_cache = rs_design_cache,
        par_cov = par_cov_current,
        par_s = par_s_current,
        margin_dist = margin_dist,
        copula_link = copula_link,
        update_only = update_only,
        eta_out = eta_out_current
      )
    } else {
      calc_eta(par_cov_current, mm, margin_dist, copula_link, par_s = par_s_current)
    }
  }

  .is_auto_stop_crit <- function(x) {
    is.null(x) || (length(x) == 1 && is.na(x))
  }

  .validate_stop_crit <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0) {
      stop(name, " must be a single positive finite number, or NA/NULL for automatic selection.")
    }
    as.numeric(x)
  }

  if (.is_auto_stop_crit(inner_stop_crit) || .is_auto_stop_crit(outer_stop_crit)) {
    eta_init_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)
    eta_init <- eta_init_out$eta_inv
    calc_lik_init <- calc_likelihood_minimal(
      eta_init,
      mm = mm$x,
      margin_dist,
      copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache,
      margin_eval_cache = margin_eval_cache
    )

    init_joint_ll <- as.numeric(calc_lik_init$log_lik["joint"])
    if (!is.finite(init_joint_ll)) {
      init_joint_ll <- 0
    }

    scale_base <- max(1, abs(init_joint_ll), nrow(dataset))
    auto_outer_stop_crit <- min(0.05, max(1e-4, 1e-6 * scale_base))
    if (identical(method, "CG")) {
      auto_outer_stop_crit <- auto_outer_stop_crit / 10
    }
    auto_inner_stop_crit <- min(0.01, max(1e-5, auto_outer_stop_crit / 5))

    if (.is_auto_stop_crit(outer_stop_crit)) {
      outer_stop_crit <- auto_outer_stop_crit
    } else {
      outer_stop_crit <- .validate_stop_crit(outer_stop_crit, "outer_stop_crit")
    }

    if (.is_auto_stop_crit(inner_stop_crit)) {
      inner_stop_crit <- auto_inner_stop_crit
    } else {
      inner_stop_crit <- .validate_stop_crit(inner_stop_crit, "inner_stop_crit")
    }

    if (verbose > 0) {
      cat("\nUsing stop criteria:",
          "inner_stop_crit=", format(inner_stop_crit, digits = 6),
          "| outer_stop_crit=", format(outer_stop_crit, digits = 6), "\n")
    }
  } else {
    inner_stop_crit <- .validate_stop_crit(inner_stop_crit, "inner_stop_crit")
    outer_stop_crit <- .validate_stop_crit(outer_stop_crit, "outer_stop_crit")
  }

  cg_grad_tol_eff <- if(.is_auto_stop_crit(cg_grad_tol)) {
    max(1e-3, 10 * outer_stop_crit)
  } else {
    .validate_stop_crit(cg_grad_tol, "cg_grad_tol")
  }
  cg_step_tol_eff <- if(.is_auto_stop_crit(cg_step_tol)) {
    max(1e-5, 0.1 * outer_stop_crit)
  } else {
    .validate_stop_crit(cg_step_tol, "cg_step_tol")
  }

  #OUTER ITERATION (MAIN LOOP)
  if(method == "CG") {
    if(!exists("calc_analytical_hessian", mode = "function")) {
      hess_path <- file.path(getwd(), "R", "analytical_hessian.R")
      if(file.exists(hess_path)) {
        source(hess_path, local = FALSE)
      }
    }
    if(!exists("calc_analytical_hessian", mode = "function")) {
      stop("CG requires calc_analytical_hessian(); source R/analytical_hessian.R first.")
    }
    if(verbose > 0) {
      cat("\nUsing optimization method: CG")
      cat(paste0(
        "\nCG controls: max_delta=", signif(cg_max_delta, 4),
        " | lambda_update_every=", cg_lambda_update_every,
        " | update_lambda=", isTRUE(cg_update_lambda),
        " | line_search=", cg_line_search,
        " | gradient=", cg_gradient_method,
        "\n"
      ))
    }

    build_cg_model <- function(mm, par_cov, par_s) {
      mm_cg <- mm
      beta <- par_cov
      for(pn in names(mm$x)) {
        if(length(mm$s[[pn]]) > 0) {
          for(sn in names(mm$s[[pn]])) {
            B <- mm$s[[pn]][[sn]]
            b <- par_s[[pn]][[sn]]
            colnames(B) <- sub(paste0("^", pn, "\\."), "", names(b))
            mm_cg$x[[pn]] <- cbind(mm_cg$x[[pn]], B)
            beta <- c(beta, b)
          }
        }
        mm_cg$s[[pn]] <- list()
      }
      list(mm = mm_cg, beta = beta)
    }

    unpack_cg_beta <- function(beta_vec) {
      par_cov_new <- beta_vec[names(par_cov)]
      par_s_new <- par_s
      for(pn in names(par_s_new)) {
        if(length(par_s_new[[pn]]) == 0) next
        for(sn in names(par_s_new[[pn]])) {
          b_names <- names(par_s_new[[pn]][[sn]])
          par_s_new[[pn]][[sn]] <- beta_vec[b_names]
        }
      }
      list(par_cov = par_cov_new, par_s = par_s_new)
    }

    build_cg_penalty <- function(beta_names, lambda_current) {
      P <- matrix(0, nrow = length(beta_names), ncol = length(beta_names),
                  dimnames = list(beta_names, beta_names))
      for(pn in names(par_s)) {
        if(length(par_s[[pn]]) == 0) next
        for(sn in names(par_s[[pn]])) {
          b_names <- names(par_s[[pn]][[sn]])
          idx <- match(b_names, beta_names)
          idx <- idx[!is.na(idx)]
          if(length(idx) == 0) next
          B <- mm$s[[pn]][[sn]]
          S <- attr(B, "penalty")
          if(is.null(S) || !is.matrix(S)) {
            D <- diff(diag(ncol(B)), differences = 2)
            S <- t(D) %*% D
          }
          P[idx, idx] <- P[idx, idx] + as.numeric(lambda_current[[pn]][[sn]]) * S
        }
      }
      P
    }

    cg_eval <- function(beta_vec, mm_cg) {
      unpacked <- unpack_cg_beta(beta_vec)
      eta_out <- calc_eta(beta_vec, mm_cg, margin_dist, copula_link,
                          par_s = setNames(lapply(names(mm_cg$x), function(x) list()), names(mm_cg$x)))
      eta_inv <- eta_out$eta_inv
      if(any(!is.finite(unlist(eta_inv, use.names = FALSE)))) return(NULL)
      if("theta" %in% names(eta_inv) && any(abs(eta_inv$theta) >= 0.999, na.rm = TRUE)) return(NULL)
      positive_names <- intersect(names(eta_inv), c("mu", "sigma", "tau", "zeta"))
      for(pn in positive_names) {
        if(any(eta_inv[[pn]] <= 1e-8, na.rm = TRUE)) return(NULL)
      }
      lik <- tryCatch(calc_likelihood_minimal(
        eta_inv, mm = mm_cg$x, margin_dist, copula_dist, calc_d2 = FALSE,
        response = dataset$response, response_margin = dataset$time,
        response_subject = dataset$subject, pair_cache = pair_cache,
        margin_eval_cache = margin_eval_cache
      ), error = function(e) NULL)
      if(is.null(lik)) return(NULL)
      list(loglik = as.numeric(lik$log_lik["joint"]), calc_lik = lik, eta_out = eta_out,
           par_cov = unpacked$par_cov, par_s = unpacked$par_s)
    }

    cg_objective <- function(beta_vec, loglik, penalty_current) {
      as.numeric(loglik) - 0.5 * sum(as.numeric(beta_vec) * as.numeric(penalty_current %*% beta_vec))
    }

    cg_gradient <- function(beta_vec, base_ll, mm_cg) {
      grad <- rep(0, length(beta_vec))
      names(grad) <- names(beta_vec)
      for(ii in seq_along(beta_vec)) {
        hk <- 1e-5 * max(1, abs(beta_vec[ii]))
        bp <- beta_vec
        bp[ii] <- bp[ii] + hk
        lp <- cg_eval(bp, mm_cg)
        lpv <- if(is.null(lp)) NA_real_ else lp$loglik
        if(identical(cg_gradient_method, "forward")) {
          if(is.finite(lpv) && is.finite(base_ll)) {
            grad[ii] <- (lpv - base_ll) / hk
          } else {
            bm <- beta_vec
            bm[ii] <- bm[ii] - hk
            lm <- cg_eval(bm, mm_cg)
            lmv <- if(is.null(lm)) NA_real_ else lm$loglik
            if(is.finite(lmv) && is.finite(base_ll)) grad[ii] <- (base_ll - lmv) / hk
          }
        } else {
          bm <- beta_vec
          bm[ii] <- bm[ii] - hk
          lm <- cg_eval(bm, mm_cg)
          lmv <- if(is.null(lm)) NA_real_ else lm$loglik
          if(is.finite(lpv) && is.finite(lmv)) grad[ii] <- (lpv - lmv) / (2 * hk)
          else if(is.finite(lpv) && is.finite(base_ll)) grad[ii] <- (lpv - base_ll) / hk
          else if(is.finite(lmv) && is.finite(base_ll)) grad[ii] <- (base_ll - lmv) / hk
        }
      }
      grad
    }

    cg_finite_hessian_block <- function(beta_vec, block_names, mm_cg, h = 1e-4) {
      block_names <- intersect(block_names, names(beta_vec))
      n_block <- length(block_names)
      H_block <- matrix(NA_real_, n_block, n_block,
                        dimnames = list(block_names, block_names))
      if(n_block == 0L) return(H_block)
      eval_base <- cg_eval(beta_vec, mm_cg)
      f0 <- if(is.null(eval_base)) NA_real_ else eval_base$loglik
      if(!is.finite(f0)) return(H_block)

      eval_ll <- function(beta_try) {
        out <- cg_eval(beta_try, mm_cg)
        if(is.null(out)) NA_real_ else out$loglik
      }

      for(ii in seq_len(n_block)) {
        ni <- block_names[ii]
        hi <- h * max(1, abs(beta_vec[ni]))
        bp <- beta_vec
        bm <- beta_vec
        bp[ni] <- bp[ni] + hi
        bm[ni] <- bm[ni] - hi
        fp <- eval_ll(bp)
        fm <- eval_ll(bm)
        if(is.finite(fp) && is.finite(fm)) {
          H_block[ii, ii] <- (fp - 2 * f0 + fm) / (hi^2)
        }

        if(ii < n_block) {
          for(jj in seq.int(ii + 1L, n_block)) {
            nj <- block_names[jj]
            hj <- h * max(1, abs(beta_vec[nj]))
            bpp <- beta_vec
            bpm <- beta_vec
            bmp <- beta_vec
            bmm <- beta_vec
            bpp[ni] <- bpp[ni] + hi
            bpp[nj] <- bpp[nj] + hj
            bpm[ni] <- bpm[ni] + hi
            bpm[nj] <- bpm[nj] - hj
            bmp[ni] <- bmp[ni] - hi
            bmp[nj] <- bmp[nj] + hj
            bmm[ni] <- bmm[ni] - hi
            bmm[nj] <- bmm[nj] - hj
            fpp <- eval_ll(bpp)
            fpm <- eval_ll(bpm)
            fmp <- eval_ll(bmp)
            fmm <- eval_ll(bmm)
            if(all(is.finite(c(fpp, fpm, fmp, fmm)))) {
              H_block[ii, jj] <- (fpp - fpm - fmp + fmm) / (4 * hi * hj)
              H_block[jj, ii] <- H_block[ii, jj]
            }
          }
        }
      }
      H_block
    }

    cg_smooth_edf_list <- function(H_obs_current, penalty_current, beta_names) {
      edf_out <- setNames(lapply(names(par_s), function(x) list()), names(par_s))
      for(pn in names(par_s)) {
        if(length(par_s[[pn]]) == 0) next
        for(sn in names(par_s[[pn]])) {
          idx <- match(names(par_s[[pn]][[sn]]), rownames(H_obs_current))
          idx <- idx[!is.na(idx)]
          if(length(idx) == 0) next
          H_block <- H_obs_current[idx, idx, drop = FALSE]
          P_block <- penalty_current[idx, idx, drop = FALSE]
          info_block <- -0.5 * (H_block + t(H_block))
          if(sum(diag(info_block), na.rm = TRUE) < 0) {
            info_block <- -info_block
          }
          info_block <- tryCatch({
            eg <- eigen(0.5 * (info_block + t(info_block)), symmetric = TRUE)
            eg$values[eg$values < 0] <- 0
            eg$vectors %*% diag(eg$values, nrow = length(eg$values)) %*% t(eg$vectors)
          }, error = function(e) info_block)
          P_block <- 0.5 * (P_block + t(P_block))
          edf_val <- tryCatch({
            k <- nrow(info_block)
            ridge <- max(1e-8, 1e-8 * max(1, max(abs(diag(info_block)), na.rm = TRUE)))
            sum(diag(.solve_linear_system(info_block + P_block + diag(ridge, k), info_block)))
          }, error = function(e) NA_real_)
          if(!is.finite(edf_val)) edf_val <- length(idx)
          edf_out[[pn]][[sn]] <- max(0, min(length(idx), as.numeric(edf_val)))
        }
      }
      edf_out
    }

    cg_update_lambda_once <- function(H_obs_current, beta_vec, grad_vec, lambda_current, mm_cg, trust_radius) {
      lambda_new <- lambda_current
      if(!isTRUE(cg_update_lambda)) return(lambda_new)
      for(pn in names(lambda_new)) {
        if(length(lambda_new[[pn]]) == 0) next
        for(sn in names(lambda_new[[pn]])) {
          lambda0 <- as.numeric(lambda_new[[pn]][[sn]])
          if(!is.finite(lambda0) || lambda0 <= 0) lambda0 <- 1
          candidates <- unique(pmax(0.01, pmin(1e6, lambda0 * c(0.1, 0.25, 0.5, 1, 2, 4, 10))))
          gaic_score <- rep(Inf, length(candidates))
          penalty_value <- rep(NA_real_, length(candidates))
          penalized_loglik <- rep(NA_real_, length(candidates))
          raw_loglik <- rep(NA_real_, length(candidates))
          edf_values <- rep(NA_real_, length(candidates))
          for(jj in seq_along(candidates)) {
            lambda_try <- lambda_new
            lambda_try[[pn]][[sn]] <- candidates[jj]
            P_try <- build_cg_penalty(names(beta_vec), lambda_try)
            g_try <- grad_vec - as.numeric(P_try %*% beta_vec)
            H_try <- H_obs_current - P_try
            delta <- tryCatch(-as.numeric(.solve_linear_system(H_try, g_try)), error = function(e) NULL)
            if(is.null(delta) || !all(is.finite(delta))) next
            dnorm <- sqrt(sum(delta^2))
            if(is.finite(dnorm) && dnorm > trust_radius) delta <- delta * trust_radius / dnorm
            dc <- max(abs(delta))
            if(is.finite(dc) && dc > cg_max_delta) delta <- delta * cg_max_delta / dc
            beta_try <- beta_vec + delta
            eval_try <- cg_eval(beta_try, mm_cg)
            if(is.null(eval_try) || !is.finite(eval_try$loglik)) next
            edf_try <- sum(unlist(cg_smooth_edf_list(H_obs_current, P_try, names(beta_vec))), na.rm = TRUE)
            penalty_try <- sum(as.numeric(beta_try) * as.numeric(P_try %*% beta_try))
            raw_loglik[jj] <- eval_try$loglik
            edf_values[jj] <- edf_try
            penalty_value[jj] <- penalty_try
            penalized_loglik[jj] <- cg_objective(beta_try, eval_try$loglik, P_try)
            gaic_score[jj] <- -2 * eval_try$loglik + lambda_penalty_K * edf_try
          }
          best <- which.max(penalized_loglik)
          if(length(best) == 1 && is.finite(penalized_loglik[best])) {
            trace_rows <- data.frame(
              outer_iteration = outer_only_run_counter,
              parameter = pn,
              smooth = sn,
              lambda_before = lambda0,
              lambda_candidate = candidates,
              raw_logLik_after_step = raw_loglik,
              smooth_penalty_after_step = penalty_value,
              penalized_logLik_after_step = penalized_loglik,
              edf_after_step = edf_values,
              gaic_score = gaic_score,
              chosen = seq_along(candidates) == best,
              row.names = NULL
            )
            cg_lambda_trace <<- rbind(cg_lambda_trace, trace_rows)
            lambda_new[[pn]][[sn]] <- candidates[best]
            if(verbose > 1) {
              cat(paste0("\nCG lambda update for ", pn, " - ", sn, ": ",
                         signif(lambda0, 4), " -> ", signif(candidates[best], 4)))
            }
          }
        }
      }
      lambda_new
    }

    cg_aug <- build_cg_model(mm, par_cov, par_s)
    mm_cg <- cg_aug$mm
    beta_all <- cg_aug$beta
    penalty_mat <- build_cg_penalty(names(beta_all), lambda_s)
    cg_trust_radius <- as.numeric(cg_max_delta)
    cg_stall_count <- 0L
    cg_converged <- FALSE
    cg_lambda_update_count <- 0L
    cg_has_smooths <- length(unlist(lambda_s, use.names = FALSE)) > 0L
    cg_lambda_trace <- data.frame()
    cg_step_trace <- list()

    while(!cg_converged && outer_only_run_counter < max_outer_iter) {
      check_elapsed_budget("CG outer iteration")
      cat(paste("\nOUTER ITERATION:", outer_only_run_counter))
      cg_trust_radius_start <- cg_trust_radius
      eval_start <- cg_eval(beta_all, mm_cg)
      if(is.null(eval_start) || !is.finite(eval_start$loglik)) stop("CG failed: current likelihood is not finite.")
      log_lik_history <- rbind(log_lik_history, eval_start$calc_lik$log_lik)
      par_history <- rbind(par_history, eval_start$par_cov[colnames(par_history)])
      outer_start_log_lik <- eval_start$loglik
      if(is.finite(outer_start_log_lik) && outer_start_log_lik > cg_best_raw_loglik) {
        cg_best_raw_loglik <- outer_start_log_lik
        cg_best_iteration <- outer_only_run_counter
      }
      obj_start <- cg_objective(beta_all, outer_start_log_lik, penalty_mat)
      grad <- if(identical(cg_gradient_method, "analytical")) {
        .cg_analytical_gradient(
          beta_all,
          mm_cg,
          eval_start$eta_out,
          eval_start$calc_lik,
          margin_dist,
          copula_dist,
          include_dlcopdpar,
          dataset$response,
          dataset$time,
          dataset$subject
        )
      } else {
        cg_gradient(beta_all, outer_start_log_lik, mm_cg)
      }

      tmp_obj <- list(
        response = dataset$response,
        response_margin = dataset$time,
        response_subject = dataset$subject,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        model_matrix = mm_cg,
        par = beta_all,
        par_s = setNames(lapply(names(mm_cg$x), function(x) list()), names(mm_cg$x))
      )
      H_obs <- calc_analytical_hessian(tmp_obj, progress = FALSE)
      H_zeta_fd <- NULL
      lambda_changed <- FALSE
      if(identical(cg_zeta_hessian, "finite")) {
        zeta_names <- grep("^zeta\\.", names(beta_all), value = TRUE)
        if(length(zeta_names) > 0L) {
          H_zeta_fd <- cg_finite_hessian_block(beta_all, zeta_names, mm_cg)
          if(all(is.finite(H_zeta_fd))) {
            H_obs[zeta_names, zeta_names] <- 0.5 * (H_zeta_fd + t(H_zeta_fd))
          } else if(verbose > 0) {
            cat("\nCG finite zeta Hessian skipped because the block was not finite.")
          }
        }
      }
      if(isTRUE(cg_update_lambda) && outer_only_run_counter > 1 &&
         cg_lambda_update_count < cg_max_lambda_updates &&
         (outer_only_run_counter %% cg_lambda_update_every == 0L)) {
        lambda_before <- lambda_s
        lambda_s <- cg_update_lambda_once(H_obs, beta_all, grad, lambda_s, mm_cg, cg_trust_radius)
        penalty_mat <- build_cg_penalty(names(beta_all), lambda_s)
        lambda_changed <- !isTRUE(all.equal(
          unlist(lambda_before, use.names = TRUE),
          unlist(lambda_s, use.names = TRUE),
          tolerance = 1e-12,
          check.attributes = FALSE
        ))
        if(isTRUE(lambda_changed)) {
          cg_trust_radius <- max(cg_step_tol_eff, cg_trust_radius / 2)
          if(verbose > 0) {
            cat(paste0("\nCG trust radius shrunk after lambda update to ", signif(cg_trust_radius, 4)))
          }
        }
        cg_lambda_update_count <- cg_lambda_update_count + 1L
      }
      df_s <- cg_smooth_edf_list(H_obs, penalty_mat, names(beta_all))

      g_pen <- grad - as.numeric(penalty_mat %*% beta_all)
      H_pen <- H_obs - penalty_mat
      candidate_steps <- list()
      grad_norm <- sqrt(sum(g_pen^2))
      if(is.finite(grad_norm) && grad_norm > 0) {
        candidate_steps[[length(candidate_steps) + 1L]] <- as.numeric(cg_trust_radius * g_pen / grad_norm)
      }
      for(ridge in c(0, 1e-8, 1e-6, 1e-4, 1e-2, 1, 10, 100)) {
        d <- tryCatch(-as.numeric(.solve_linear_system(H_pen - diag(ridge, nrow(H_pen)), g_pen)), error = function(e) NULL)
        if(!is.null(d) && all(is.finite(d))) {
          candidate_steps[[length(candidate_steps) + 1L]] <- d
          candidate_steps[[length(candidate_steps) + 1L]] <- -d
        }
      }

      best <- NULL
      line_eval_count <- 0L
      stop_line_search <- FALSE
      max_backtrack <- if(isTRUE(use_backtracking)) as.integer(backtracking_max_halves) else 0L
      for(delta0 in candidate_steps) {
        if(stop_line_search) break
        for(bt in seq_len(max_backtrack + 1L)) {
          if(line_eval_count >= cg_max_line_search_evals) {
            stop_line_search <- TRUE
            break
          }
          delta <- delta0 / (2 ^ (bt - 1L))
          dnorm <- sqrt(sum(delta^2))
          if(is.finite(dnorm) && dnorm > cg_trust_radius) delta <- delta * cg_trust_radius / dnorm
          dc <- max(abs(delta))
          if(is.finite(dc) && dc > cg_max_delta) delta <- delta * cg_max_delta / dc
          beta_try <- beta_all + delta
          line_eval_count <- line_eval_count + 1L
          eval_try <- cg_eval(beta_try, mm_cg)
          if(is.null(eval_try) || !is.finite(eval_try$loglik)) next
          obj_try <- cg_objective(beta_try, eval_try$loglik, penalty_mat)
          improvement <- obj_try - obj_start
          if(is.finite(improvement) && improvement > max(1e-8, cg_armijo_c1 * max(1, abs(obj_start)))) {
            if(is.null(best) || improvement > best$improvement) {
              best <- list(beta = beta_try, eval = eval_try, improvement = improvement,
                           step_l2 = sqrt(sum(delta^2)))
            }
            if(identical(cg_line_search, "first")) {
              stop_line_search <- TRUE
              break
            }
          }
        }
      }
      if(verbose > 1) {
        cat(paste0("\nCG line search likelihood evaluations: ", line_eval_count))
      }

      cg_prevented_deterioration <- FALSE
      cg_prevented_raw_loglik_drop <- NA_real_
      accepted_improvement <- NA_real_
      if(is.null(best)) {
        cg_stall_count <- cg_stall_count + 1L
        cg_trust_radius <- max(cg_trust_radius / 2, cg_step_tol_eff)
        calc_lik_out_end <- eval_start$calc_lik
        if(verbose > 0) cat(paste0("\nCG step rejected (stall ", cg_stall_count, "/", cg_max_stall, ")\n"))
      } else {
        prospective_best_raw_loglik <- max(cg_best_raw_loglik, outer_start_log_lik, na.rm = TRUE)
        prospective_raw_loglik_drop <- prospective_best_raw_loglik - best$eval$loglik
        cg_prevented_deterioration <- is.finite(cg_raw_loglik_drop_tol) &&
          cg_lambda_update_count > 0L &&
          is.finite(prospective_raw_loglik_drop) &&
          prospective_raw_loglik_drop >= cg_raw_loglik_drop_tol
        if(isTRUE(cg_prevented_deterioration)) {
          cg_prevented_raw_loglik_drop <- prospective_raw_loglik_drop
          calc_lik_out_end <- eval_start$calc_lik
          best <- NULL
        } else {
          beta_all <- best$beta
          unpacked <- unpack_cg_beta(beta_all)
          par_cov <- unpacked$par_cov
          par_s <- unpacked$par_s
          calc_lik_out_end <- best$eval$calc_lik
          cg_stall_count <- 0L
          accepted_improvement <- best$improvement
          if(is.finite(best$step_l2) && is.finite(cg_trust_radius) &&
             best$step_l2 >= 0.8 * cg_trust_radius) {
            cg_trust_radius <- min(as.numeric(cg_max_delta), max(cg_step_tol_eff, 1.5 * cg_trust_radius))
          }
        }
      }

      outer_end_log_lik <- as.numeric(calc_lik_out_end$log_lik["joint"])
      outer_log_lik_change <- outer_end_log_lik - outer_start_log_lik
      if(is.finite(outer_end_log_lik) && outer_end_log_lik > cg_best_raw_loglik) {
        cg_best_raw_loglik <- outer_end_log_lik
        cg_best_iteration <- outer_only_run_counter
      }
      cg_raw_loglik_drop_from_best <- cg_best_raw_loglik - outer_end_log_lik
      if(isTRUE(cg_prevented_deterioration) && is.finite(cg_prevented_raw_loglik_drop)) {
        cg_raw_loglik_drop_from_best <- max(cg_raw_loglik_drop_from_best, cg_prevented_raw_loglik_drop, na.rm = TRUE)
      }
      out_temp <- c(outer_start_log_lik, outer_end_log_lik, outer_log_lik_change)
      names(out_temp) <- c("Start LogLik", "End LogLik", "Change")
      cat("\n")
      print(out_temp)

      grad_inf <- max(abs(g_pen), na.rm = TRUE)
      step_l2 <- if(is.null(best)) 0 else best$step_l2
      cg_last_grad_inf <- grad_inf
      cg_last_step_l2 <- step_l2
      cg_tolerance_met <- abs(outer_log_lik_change) <= outer_stop_crit &&
        is.finite(grad_inf) && grad_inf <= cg_grad_tol_eff &&
        is.finite(step_l2) && step_l2 <= cg_step_tol_eff
      cg_max_stall_hit <- cg_stall_count >= cg_max_stall
      cg_deterioration_hit <- is.finite(cg_raw_loglik_drop_tol) &&
        cg_lambda_update_count > 0L &&
        is.finite(cg_raw_loglik_drop_from_best) &&
        cg_raw_loglik_drop_from_best >= cg_raw_loglik_drop_tol
      cg_deterioration_hit <- isTRUE(cg_deterioration_hit) || isTRUE(cg_prevented_deterioration)
      cg_stop_requested <- cg_max_stall_hit || cg_tolerance_met || cg_deterioration_hit

      cg_step_trace[[length(cg_step_trace) + 1L]] <- data.frame(
        outer_iteration = as.integer(outer_only_run_counter),
        start_logLik = as.numeric(outer_start_log_lik),
        end_logLik = as.numeric(outer_end_log_lik),
        raw_logLik_change = as.numeric(outer_log_lik_change),
        start_penalized_logLik = as.numeric(obj_start),
        accepted_penalized_improvement = as.numeric(accepted_improvement),
        grad_inf = as.numeric(grad_inf),
        step_l2 = as.numeric(step_l2),
        trust_radius_start = as.numeric(cg_trust_radius_start),
        trust_radius_end = as.numeric(cg_trust_radius),
        line_search_evals = as.integer(line_eval_count),
        accepted_step = !is.null(best),
        lambda_update_count = as.integer(cg_lambda_update_count),
        lambda_changed = isTRUE(lambda_changed),
        stall_count = as.integer(cg_stall_count),
        tolerance_met = isTRUE(cg_tolerance_met),
        max_stall_hit = isTRUE(cg_max_stall_hit),
        raw_deterioration_hit = isTRUE(cg_deterioration_hit),
        raw_loglik_drop_from_best = as.numeric(cg_raw_loglik_drop_from_best),
        row.names = NULL
      )

      if(isTRUE(cg_stop_requested)) {
        if(isTRUE(cg_update_lambda) && isTRUE(cg_has_smooths) && cg_lambda_update_count == 0L) {
          lambda_s <- cg_update_lambda_once(H_obs, beta_all, grad, lambda_s, mm_cg, cg_trust_radius)
          penalty_mat <- build_cg_penalty(names(beta_all), lambda_s)
          df_s <- cg_smooth_edf_list(H_obs, penalty_mat, names(beta_all))
          cg_lambda_update_count <- cg_lambda_update_count + 1L
          cg_stall_count <- 0L
          if(verbose > 0) {
            cat("\nCG convergence delayed for first smoother lambda update")
          }
        } else {
          cg_stop_reason <- if(isTRUE(cg_tolerance_met)) {
            "tolerance"
          } else if(isTRUE(cg_deterioration_hit)) {
            "raw_loglik_deterioration"
          } else {
            "max_stall"
          }
          if(identical(cg_stop_reason, "tolerance")) cat("\nOUTER CONVERGED")
          if(identical(cg_stop_reason, "raw_loglik_deterioration") && verbose > 0) {
            cat(paste0(
              "\nCG stopped after raw log-likelihood dropped ",
              signif(cg_raw_loglik_drop_from_best, 5),
              " below best seen value."
            ))
          }
          cg_converged <- TRUE
        }
      }

      outer_only_run_counter <- outer_only_run_counter + 1L
    }

    final_obj <- list(
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      model_matrix = mm_cg,
      par = beta_all,
      par_s = setNames(lapply(names(mm_cg$x), function(x) list()), names(mm_cg$x))
    )
    final_H <- tryCatch(calc_analytical_hessian(final_obj, progress = FALSE), error = function(e) NULL)
    if(!is.null(final_H)) {
      penalty_mat <- build_cg_penalty(names(beta_all), lambda_s)
      df_s <- cg_smooth_edf_list(final_H, penalty_mat, names(beta_all))
    }

    for(pn in names(mm$x)) weights_final[[pn]] <- rep(1, nrow(mm$x[[pn]]))
  } else {
  while ((first_outer_run==TRUE | (abs(outer_log_lik_change)>outer_stop_crit)) & outer_only_run_counter < max_outer_iter) {
    check_elapsed_budget("RS outer iteration")

    cat(paste("\nOUTER ITERATION:",outer_only_run_counter))
    first_outer_run=TRUE

    # RUN INNER ITERATION FOR EACH PARAMETER
    for (par_name in names(mm$x)) {

      if(verbose > 2) {
        cat(paste("\nINNER ITERATION: Parameter:",par_name))
      }

      first_inner_run=TRUE; change_log_lik=0; beta_change_inner=99
      run_counter=1
      inner_run_counter=1

      # INNER ITERATION (GLIM)
      while ( (first_inner_run==TRUE | abs(change_log_lik)>inner_stop_crit) & inner_run_counter<max_inner_iter) { #

        timer=c()
        timer_start=Sys.time()

        first_inner_run=FALSE

        eta_out=rs_calc_eta(par_cov_current = par_cov, par_s_current = par_s)
        eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv

        # Guard against silent row dropping in matrix construction when response has NAs.
        n_resp <- length(dataset$response)
        margin_params <- intersect(names(mm$x), c("mu", "sigma", "nu", "tau"))
        bad_lengths <- margin_params[sapply(margin_params, function(pn) length(eta_inv[[pn]]) != n_resp)]
        if (length(bad_lengths) > 0) {
          detail <- paste(sapply(bad_lengths, function(pn) {
            paste0(pn, "=", length(eta_inv[[pn]]), " vs response=", n_resp)
          }), collapse = ", ")
          stop(
            "ERROR: Parameter vector lengths do not match response length. ",
            detail,
            ".\nThis usually indicates model-matrix rows were dropped (often due to NA handling)."
          )
        }

        calc_lik_out=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=FALSE
          ,response=dataset$response,response_margin=(dataset$time),response_subject = dataset$subject
          ,pair_cache=pair_cache,margin_eval_cache=margin_eval_cache)
        log_lik=calc_lik_out$log_lik; margin_d=calc_lik_out$margin_d; margin_p=calc_lik_out$margin_p;
        margin_deriv=calc_lik_out$margin_deriv; copula_d=calc_lik_out$copula_d; copula_p=calc_lik_out$copula_p;
        Fx_1_2=calc_lik_out$Fx_1_2;order_copula=calc_lik_out$order_copula

        if(first_outer_run==TRUE) {
          outer_start_log_lik=log_lik["joint"]; first_outer_run=FALSE
        }

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Calc Lik")

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Numerical Derivatives")

        #Capturing log lik and parameter estimates for each iteration
        log_lik_history=rbind(log_lik_history,calc_lik_out$log_lik)
        par_history=rbind(par_history,par_cov[colnames(par_history)])

        #Fixing extreme values if they exist, though they shouldn't
        Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

        ########CALCULATE COPULA DERIVATIVES
        copula_derivatives=calc_copula_derivatives(
          eta_inv,
          Fx_1_2,
          copula_dist,
          par1 = calc_lik_out$copula_par1,
          par2 = calc_lik_out$copula_par2,
          pair_complete = calc_lik_out$pair_complete
        )
        dldth=copula_derivatives$dldth; dcdth=copula_derivatives$dcdth; dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2
        if("zeta" %in% names(eta_inv)) {dldz=copula_derivatives$dldz; dcdz=copula_derivatives$dcdz}

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Copula Derivatives")

        ### Calculate copula derivatives w.r.t margin parameters
        if(!par_name %in% c("mu","sigma","nu","tau")) {
          if(par_name == "theta") {
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)
            row_id1 <- calc_lik_out$copula_row_id1
            if(length(row_id1)>0) {
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                par_idx <- calc_lik_out$copula_theta_index_map[row_id1]
              }
              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              if(any(valid_idx)) {
                d1_sum <- rowsum(dldth[valid_idx], par_idx[valid_idx], reorder = FALSE)
                d1_full[as.integer(rownames(d1_sum)),1] <- d1_sum[,1]
              }
            }
            d1=as.matrix(d1_full)
            colnames(d1)="dldtheta"
            #d2=d2ldth2
          } else if(par_name == "zeta") {
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)
            row_id1 <- calc_lik_out$copula_row_id1
            if(length(row_id1)>0) {
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                par_idx <- calc_lik_out$copula_theta_index_map[row_id1]
              }
              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              if(any(valid_idx)) {
                d1_sum <- rowsum(dldz[valid_idx], par_idx[valid_idx], reorder = FALSE)
                d1_full[as.integer(rownames(d1_sum)),1] <- d1_sum[,1]
              }
            }
            d1=as.matrix(d1_full)
            colnames(d1)="dldzeta"
            #d2=d2ldz2
          } else {
            stop("Unexpected copula parameter in optimisation: ", par_name)
          }
        } else {

          ### MARGIN LIKELIHOOD DERIVATIVES
          margin_deriv_subnames=c("m","d","v","t")
          names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
          margin_par=names(mm$x)[names(mm$x) %in% c("mu","sigma","nu","tau")]
          response=dataset$response

          d1=as.matrix(margin_deriv[grepl(paste("dld",margin_deriv_subnames[par_name],sep=""),names(margin_deriv))][[1]])
          colnames(d1)=paste("dld",par_name,sep="")

          if(include_dlcopdpar==TRUE) {

            #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
            nd_impact_F=calc_Fx_derivatives(eta_inv,mm$x,margin_dist,response=dataset$response,par_names=par_name)

            # Calculate copula derivative with respect to marginal parameters.
            # The old path built a pair-expanded data frame using two merge()
            # calls. The likelihood path already has stable observation row
            # ids, so accumulate the endpoint contributions directly.
            d1_cop <- .calc_dlcopdpar_indexed(
              row_id1 = calc_lik_out$copula_row_id1,
              row_id2 = calc_lik_out$copula_row_id2,
              dcdu1 = dcdu1,
              dcdu2 = dcdu2,
              copula_d = copula_d,
              F_nd = nd_impact_F[[par_name]],
              n_obs = length(dataset$response),
              pair_complete = calc_lik_out$pair_complete
            )
            d1_m=d1
            d1=d1_m+d1_cop
            #d1=d1*0+(nd_impact[par_name]/nrow(d1))

            if (check_dlcopdpar_gradient && outer_only_run_counter == 1) {
              gradient_check <- check_dlcopdpar_gradient_margin_score(
                eta = eta,
                eta_inv = eta_inv,
                par_name = par_name,
                margin_dist = margin_dist,
                copula_dist = copula_dist,
                dataset = dataset,
                mm = mm$x,
                pair_cache = pair_cache,
                d1 = d1,
                base_loglik = log_lik["joint"],
                verbose = verbose
              )
              if (isTRUE(gradient_check$warned)) {
                warning(gradient_check$message, call. = FALSE)
              }
            }
          }

          timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
          names(timer)[length(timer)]=paste("Margin Derivatives")

          if(verbose>=4) {print(timer)}

          d1=d1[,grepl(par_name,colnames(d1))]

          if(include_dlcopdpar==FALSE) {d1_cop=d1*0; d1_m=d1}

          #nd=round(c(nd_impact[par_name],nd_impact_m[par_name],nd_impact_c[par_name],sum(d1),sum(d1_m),sum(d1_cop*0.5)),2)
          #names(nd)=c("joint_nd","marginal_nd","copula_nd","joint_calc","margin_calc","copula_calc")
          #print(nd)
        }

        ### INNER ITERATION / BACKFITTING STEP

        # Ensure score inputs have consistent lengths (prevents silent recycling).
        eta_len <- length(eta[[par_name]])
        d1 <- as.numeric(d1)
        eta_dr_vec <- as.numeric(eta_dr[[par_name]])

        if (length(d1) != eta_len) {
          stop(
            "Score derivative length mismatch for ", par_name,
            ": length(d1)=", length(d1),
            " but length(eta)=", eta_len,
            ". This indicates an index-alignment bug in derivative assembly."
          )
        }

        if (length(eta_dr_vec) != eta_len) {
          stop(
            "Link-derivative length mismatch for ", par_name,
            ": length(eta_dr)=", length(eta_dr_vec),
            " but length(eta)=", eta_len,
            "."
          )
        }

        # 1. Calculate y_k, w_k

        ########### FIRST ITERATION CALCULATES B_k without smooths
        score=score_function_v2(eta=eta[[par_name]],dldpar=d1,d2ldpar=-(d1*d1),dpardeta=eta_dr_vec)

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Backfitting")

        # Setup model matrices
        design_info <- rs_design_cache[[par_name]]
        X <- design_info$X
        fixed_names <- design_info$fixed_names
        smooth_penalty_meta <- design_info$smooth_penalty_meta
        w_k_vec=as.vector(score$w_k)

        z_k=score$z_k
        if(length(par_s[[par_name]])==0) {
          paste("No smooths found for parameter; running basic IRLS",par_name)
          beta_start=c(par_cov[fixed_names])
        } else {
            ############# UNPENALISED VERSION
          temp_par_s_unlisted=unlist(par_s[[par_name]],use.names=FALSE)
          names(temp_par_s_unlisted)=setdiff(colnames(X), fixed_names)
          beta_start=c(par_cov[fixed_names],temp_par_s_unlisted)
        }

        backfitting_iteration <- function(
          par_s,
          par_cov,
          beta_start,
          lambda_s,
          first_inner_run,
          K,
          margin_dist,
          copula_dist,
          dataset,
          mm,
          copula_link,
          df_s,
          step_size,
          par_name
        ) {

          ############# Backfitting with penalisation
          pen_mat=matrix(0,nrow=ncol(X),ncol=ncol(X))
          if(length(par_s[[par_name]])>0) {
            for (s_name in names(smooth_penalty_meta)) {
              meta=smooth_penalty_meta[[s_name]]
              B=meta$B
              idx=meta$idx
              S=meta$S_base
              # Always apply the current lambda penalty (lambda_s is initialised
              # to lambda_start, so the first outer iteration is not unpenalised).
              pen_mat[idx,idx]=S*lambda_s[[par_name]][[s_name]]
              # Weighted effective DF: tr((B'WB + lambda S)^(-1) B'WB)
              # Uses IRLS weights w_k_vec from the enclosing scope.
              # This is both correct and avoids building an n by n hat matrix.
              BtWB_s <- t(B) %*% (B * as.vector(w_k_vec))
              df_s[[par_name]][[s_name]] <- sum(.solve_linear_system(BtWB_s + pen_mat[idx,idx]) * BtWB_s)
            }
          }

          XtWX = t(X) %*% (X * w_k_vec)
          XtWz = t(X) %*% (z_k * w_k_vec)
          beta_update=as.vector(.solve_linear_system(XtWX + pen_mat, XtWz))
          beta_change_inner=beta_update-beta_start
          beta_new=beta_start*(1-step_size) + (step_size)*(beta_update)
          beta_new

          temp_par_cov_new=beta_new[fixed_names]
          par_cov_new=c(temp_par_cov_new,par_cov[!names(par_cov) %in% names(temp_par_cov_new)])
          #par_cov_new[names(beta)]=beta
          temp_par_s_new=beta_new[!names(beta_new) %in% names(par_cov_new)]
          #Select all beta_new which have names corresponding to s_name
          par_s_new=par_s
          for(s_name in names(par_s[[par_name]])) {
            smooth_col_names=colnames(X)[smooth_penalty_meta[[s_name]]$idx]
            par_s_new[[par_name]][[s_name]]=temp_par_s_new[smooth_col_names]
          }

          eta_out=rs_calc_eta(par_cov_current = par_cov_new, par_s_current = par_s_new)

          eta=eta_out$eta; eta_dr=eta_out$eta_dr; eta_inv=eta_out$eta_inv
          par_cov=par_cov_new
          par_s=par_s_new

          if(par_name %in% c("theta", "zeta") && isTRUE(getOption("gamlss.longitudinal.fast_copula_lik", TRUE))) {
            calc_lik_out_end=.calc_likelihood_update_copula(
              eta_inv = eta_inv,
              base_lik = calc_lik_out,
              copula_dist = copula_dist,
              pair_cache = pair_cache
            )
          } else {
            calc_lik_out_end=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=FALSE
              ,response=dataset$response,response_margin=(dataset$time),response_subject = dataset$subject
              ,pair_cache=pair_cache,margin_eval_cache=margin_eval_cache)
          }


          #print(sum(unlist(df_s[[par_name]])))
          GAIC_lambda_k=-2*calc_lik_out_end$log_lik["joint"]+ K*sum(unlist(df_s[[par_name]]))

          #print(paste("K*DF_S",K*sum(unlist(df_s[[par_name]]))))

          return_list=list(par_cov,par_s,calc_lik_out_end,GAIC_lambda_k,df_s)
          names(return_list)=c("par_cov","par_s","calc_lik_out_end","GAIC_lambda_k","df_s")
          return(return_list)
        }

        optim_lambda <- function(lambda_val,smooth_name,
        par_s,par_cov, beta_start, lambda_s, first_inner_run=FALSE,K=K,
                margin_dist, copula_dist, dataset, mm, copula_link,df_s,step_size,par_name) {
          lambda_s_temp=lambda_s
          lambda_s_temp[[par_name]][[smooth_name]]=lambda_val
          backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s_temp, first_inner_run=FALSE,K=K,
                margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                ,df_s=df_s,step_size=step_size,par_name=par_name)

          # Debug output
          loglik <- backfitting_iteration_results$calc_lik_out_end$log_lik["joint"]
          df_total <- sum(unlist(backfitting_iteration_results$df_s[[par_name]]))
          gaic_val <- backfitting_iteration_results$GAIC_lambda_k
          #print(sprintf("lambda=%.3f | LogLik=%.2f | DF=%.2f | GAIC=%.2f\n",
          #           lambda_val, loglik, df_total, gaic_val))

          return(backfitting_iteration_results$GAIC_lambda_k)
        }

        K=lambda_penalty_K
        num_smooths=length(lambda_s[[par_name]])
        if(num_smooths==0|outer_only_run_counter==1) {
          backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=TRUE,K=K,
                margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                ,df_s=df_s,step_size=step_size,par_name=par_name)
        } else {
          for (smooth_name in names(lambda_s[[par_name]])) {
           #Optimize lambda for each smooth
           if(inner_run_counter==1) {
             cat(paste("\nOptimising smoothing parameter for",par_name,"-",smooth_name))
              optim_lambda_out=optim(par=lambda_s[[par_name]][[smooth_name]],fn=optim_lambda,
                smooth_name=smooth_name,
                par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=FALSE,K=K,
                  margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link
                  ,df_s=df_s,step_size=step_size,par_name=par_name,
                  method="L-BFGS-B",lower=0.01,upper=1e6,control = list(factr=1,pgtol=.1)
              )
              lambda_s[[par_name]][[smooth_name]]=optim_lambda_out$par
              if(verbose>2) {
                print(paste("Chosen lambda:" ,round(lambda_s[[par_name]][[smooth_name]],2), "| Penalty K =", K))
              }
            } #end if inner_run_counter
          } #end for smooth_name
        } #end if num_smooths

        backfitting_iteration_results=backfitting_iteration(par_s=par_s,par_cov=par_cov, beta_start=beta_start, lambda_s=lambda_s, first_inner_run=FALSE,K=K,
          margin_dist=margin_dist, copula_dist=copula_dist, dataset=dataset, mm=mm, copula_link=copula_link,df_s=df_s,step_size=step_size,par_name=par_name)

        # Guard against downhill updates: if a proposed step lowers joint log-likelihood,
        # try smaller step sizes before accepting.
        start_joint_loglik <- as.numeric(calc_lik_out$log_lik["joint"])
        accepted_results <- backfitting_iteration_results
        accepted_step_size <- step_size
        proposed_joint_loglik <- as.numeric(backfitting_iteration_results$calc_lik_out_end$log_lik["joint"])
        theta_step_rejected <- FALSE
        backtracking_attempts_used <- 0L
        max_backtrack <- 0L

        if(isTRUE(use_backtracking) && is.finite(start_joint_loglik) && is.finite(proposed_joint_loglik) && proposed_joint_loglik < start_joint_loglik) {
          max_backtrack <- backtracking_max_halves
          trial_step <- step_size
          accepted <- FALSE

          for(bt in seq_len(max_backtrack)) {
            backtracking_attempts_used <- bt
            trial_step <- trial_step / 2
            trial_results <- backfitting_iteration(
              par_s=par_s,
              par_cov=par_cov,
              beta_start=beta_start,
              lambda_s=lambda_s,
              first_inner_run=FALSE,
              K=K,
              margin_dist=margin_dist,
              copula_dist=copula_dist,
              dataset=dataset,
              mm=mm,
              copula_link=copula_link,
              df_s=df_s,
              step_size=trial_step,
              par_name=par_name
            )

            trial_joint_loglik <- as.numeric(trial_results$calc_lik_out_end$log_lik["joint"])
            if(is.finite(trial_joint_loglik) && trial_joint_loglik >= start_joint_loglik) {
              accepted_results <- trial_results
              accepted_step_size <- trial_step
              accepted <- TRUE
              break
            }
          }

          if(!accepted) {
            accepted_results <- list(
              par_cov=par_cov,
              par_s=par_s,
              calc_lik_out_end=calc_lik_out,
              GAIC_lambda_k=NA_real_,
              df_s=df_s
            )
            accepted_step_size <- 0
            theta_step_rejected <- TRUE
          }

          if(verbose > 1) {
            cat(paste0(
              "\nBacktracking applied for ", par_name,
              ": step_size ", signif(step_size, 4),
              " -> ", signif(accepted_step_size, 4),
              " (halves tried=", backtracking_attempts_used,
              "/", max_backtrack, ")",
              "\n"
            ))
          }
        }

        accepted_joint_loglik <- as.numeric(accepted_results$calc_lik_out_end$log_lik["joint"])
        if(!identical(method, "CG")) {
          rs_block_trace[[length(rs_block_trace) + 1L]] <- data.frame(
            outer_iteration = as.integer(outer_only_run_counter),
            inner_iteration = as.integer(inner_run_counter),
            global_inner_iteration = as.integer(outer_run_counter),
            parameter = par_name,
            start_logLik = as.numeric(start_joint_loglik),
            proposed_logLik = as.numeric(proposed_joint_loglik),
            accepted_logLik = as.numeric(accepted_joint_loglik),
            proposed_change = as.numeric(proposed_joint_loglik - start_joint_loglik),
            accepted_change = as.numeric(accepted_joint_loglik - start_joint_loglik),
            nominal_step_size = as.numeric(step_size),
            accepted_step_size = as.numeric(accepted_step_size),
            backtracking_attempts = as.integer(backtracking_attempts_used),
            max_backtracking_attempts = as.integer(max_backtrack),
            rejected = isTRUE(theta_step_rejected),
            elapsed_sec = as.numeric(difftime(Sys.time(), timer_start, units = "secs")),
            stringsAsFactors = FALSE
          )
        }
        if(par_name == "theta" && verbose > 2) {
          cat(paste0(
            "\nTheta step diagnostics: start=", signif(start_joint_loglik, 8),
            ", proposed=", signif(proposed_joint_loglik, 8),
            ", accepted=", signif(accepted_joint_loglik, 8),
            ", backtracking=", if(isTRUE(use_backtracking)) "on" else "off",
            ", step=", signif(step_size, 4),
            ", accepted_step=", signif(accepted_step_size, 4),
            ", halves_tried=", backtracking_attempts_used,
            "/", max_backtrack,
            ", rejected=", if(theta_step_rejected) "yes" else "no",
            "\n"
          ))
        }

        par_cov=accepted_results$par_cov
        par_s=accepted_results$par_s
        calc_lik_out_end=accepted_results$calc_lik_out_end
        df_s=accepted_results$df_s

        if (verbose>2) {
          cat("\nLogLik:\n")
          print(calc_lik_out_end$log_lik)
        }

        if(plot_results==TRUE) {
          plot_count=3+length(par_cov)
          sides=round(sqrt(plot_count))

          par(mfrow=c(sides+1,sides))
          plot(log_lik_history[,3],type="l",main="LogLik - Overall")
          plot(log_lik_history[,1],type="l",main="LogLik - Margin")
          plot(log_lik_history[,2],type="l",main="LogLik - Copula")

          for(i in 1:length(colnames(par_history))) {


            if(!all(is.na(true_val))) {
              plot(par_history[,i],type="l",main=colnames(par_history)[i],xlab="Iteration",ylab="Parameter estimate",ylim=range(c(par_history[,i],true_val[i])))
              abline(h=true_val[i],col="red")
            } else {
              plot(par_history[,i],type="l",main=colnames(par_history)[i],xlab="Iteration",ylab="Parameter estimate")
            }
          }
        }

        timer=c(timer,difftime(Sys.time(),timer_start,units="secs"))
        names(timer)[length(timer)]=paste("Plotting")
        #print(timer)

        change_log_lik=calc_lik_out_end$log_lik["joint"]-calc_lik_out$log_lik["joint"]

        run_counter=run_counter+1
        outer_run_counter=outer_run_counter+1
        inner_run_counter=inner_run_counter+1

      }
      weights_final[[par_name]]=score$w_k
    }

    step_size = (step_adjustment^min(outer_only_run_counter,max_steps))*start_step_size
    outer_only_run_counter=outer_only_run_counter+1
    outer_end_log_lik=calc_lik_out_end$log_lik["joint"]
    outer_log_lik_change=outer_end_log_lik-outer_start_log_lik

    out_temp=c(outer_start_log_lik,outer_end_log_lik,outer_log_lik_change)
    names(out_temp) = c("Start LogLik","End LogLik","Change")
    cat("\n")
    print(out_temp)

    if(is.finite(outer_log_lik_change) && outer_log_lik_change < 0) {
      outer_negative_streak = outer_negative_streak + 1
    } else {
      outer_negative_streak = 0
    }

    if(outer_negative_streak >= max_negative_outer_streak) {
      msg = paste0(
        "Optimization stopped after ", max_negative_outer_streak, " consecutive negative outer log-likelihood changes. ",
        "We believe the model may be misspecified and the likelihood may be malformed. ",
        "Try different starting parameters or covariate combinations. Other options include switching between joint and separate optimisation. In general, joint optimisation provides more stable convergence.",
        "Alternatively, you can increase the max_negative_outer_streak parameter to allow more negative changes before stopping, but we recommend investigating the cause of the consecutive negative changes in likelihood."
      )
      warning(msg, call. = FALSE)
      stop(msg, call. = FALSE)
    }


    if(abs(outer_log_lik_change)<=outer_stop_crit) {
      print(c(outer_end_log_lik-outer_start_log_lik))
      cat("\nOUTER CONVERGED")
    }

  }
  }

  converged <- is.finite(outer_log_lik_change) && abs(outer_log_lik_change) <= outer_stop_crit
  if(identical(method, "CG")) {
    converged <- identical(cg_stop_reason, "tolerance")
  }
  hit_outer_limit <- outer_only_run_counter >= max_outer_iter && !isTRUE(converged)
  convergence_info <- list(
    converged = isTRUE(converged),
    hit_outer_limit = isTRUE(hit_outer_limit),
    hit_max_stall = isTRUE(identical(cg_stop_reason, "max_stall")),
    hit_raw_loglik_deterioration = isTRUE(identical(cg_stop_reason, "raw_loglik_deterioration")),
    stop_reason = if(identical(method, "CG")) cg_stop_reason else if(isTRUE(converged)) "tolerance" else NA_character_,
    grad_inf = as.numeric(cg_last_grad_inf),
    step_l2 = as.numeric(cg_last_step_l2),
    best_raw_loglik = as.numeric(cg_best_raw_loglik),
    best_raw_loglik_iteration = as.integer(cg_best_iteration),
    raw_loglik_drop_from_best = as.numeric(cg_raw_loglik_drop_from_best),
    raw_loglik_drop_tol = as.numeric(cg_raw_loglik_drop_tol),
    outer_iterations = max(0L, outer_only_run_counter - 1L),
    max_outer_iter = max_outer_iter,
    outer_log_lik_change = as.numeric(outer_log_lik_change),
    outer_stop_crit = outer_stop_crit,
    method = method,
    cg_gradient_method = if(identical(method, "CG")) cg_gradient_method else NA_character_,
    cg_zeta_hessian = if(identical(method, "CG")) cg_zeta_hessian else NA_character_
  )

  if (isTRUE(hit_outer_limit)) {
    warning(
      "Model stopped at max_outer_iter before satisfying outer_stop_crit; treat fit as not converged.",
      call. = FALSE
    )
  }

  cat("\n\n############ MODEL FIT ############\n")
  cat(paste("\nMargin distribution:",margin_dist$family[2]))
  cat(paste("\nCopula distribution:",copula_dist))
  cat("\n")
  cat(paste("\nParameter count:",length(par_cov)))
  cat(paste("\nObservations:",nrow(dataset)))
  cat(paste("\nMargins:",length(unique(dataset$time))))
  cat("\n")
  total_fit_time <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
  cat(paste("\nTotal time (seconds):",round(total_fit_time,2)))
  cat("\n\n")
  par_mat_out_temp=t(t((par_cov)))
  colnames(par_mat_out_temp) = c("estimate")
  print(par_mat_out_temp)
  cat("\n")
  cat("Model Selection Criteria:")
  cat("\n")

  p_cop=par_cov[grepl("theta",names(par_cov))|grepl("zeta",names(par_cov))]
  p_mar=par_cov[!(grepl("theta",names(par_cov))|grepl("zeta",names(par_cov)))]

  df_s_total=df_s_cop_total=df_s_margin_total=0
  for(par_name in names(par_s)) {
    if(par_name %in% c("theta","zeta"))
      df_s_cop_total=df_s_cop_total+sum(unlist(df_s[[par_name]]))
    else {
      df_s_margin_total=df_s_margin_total+sum(unlist(df_s[[par_name]]))
    }
    df_s_total=df_s_total+sum(unlist(df_s[[par_name]]))
  }

  aics=rbind(t(calc_lik_out_end$log_lik),
             t(-calc_lik_out_end$log_lik*2)+2*c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total),
             t(-calc_lik_out_end$log_lik*2)+c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total)*log(nrow(dataset)),
             t(c(length(p_mar)+df_s_margin_total,length(p_cop)+df_s_cop_total,length(par_cov)+df_s_total))
             )

  rownames(aics)=c("LogLik","AIC","BIC","EDF")
  print(aics)

  cat("\n####################################\n")

  return_list=list(par_cov,log_lik_history,par_history,calc_lik_out_end,mm,margin_dist,copula_dist,include_dlcopdpar,dataset$response,dataset$time,dataset$subject,par_s,lambda_s,df_s,weights_final)
  names(return_list)=c("par","log_lik_history","par_history","calc_lik_out_end","model_matrix","margin_dist","copula_dist","include_dlcopdpar","response","response_margin","response_subject","par_s","lambda_s","df_s","weights")
  return_list$formulas_int <- list(
    mu = mu.formula.int,
    sigma = sigma.formula.int,
    nu = nu.formula.int,
    tau = tau.formula.int,
    theta = theta.formula.int,
    zeta = zeta.formula.int
  )
  return_list$var_map <- var_map
  return_list$optim_method <- method
  return_list$warm_start_joint <- warm_start_info
  return_list$convergence <- convergence_info
  if(!identical(method, "CG")) {
    return_list$rs_block_trace <- if(length(rs_block_trace)) {
      do.call(rbind, rs_block_trace)
    } else {
      data.frame()
    }
  }
  if(identical(method, "CG")) {
    return_list$cg_lambda_trace <- cg_lambda_trace
    return_list$cg_step_trace <- if(length(cg_step_trace)) {
      do.call(rbind, cg_step_trace)
    } else {
      data.frame()
    }
  }

  # Store vcov metadata and optionally precompute vcov once at fit time.
  return_list$vcov <- NULL
  return_list$vcov_meta <- list(
    precomputed = FALSE,
    numderiv = isTRUE(vcov_numderiv),
    method = vcov_method
  )

  if(isTRUE(compute_vcov)) {
    if(verbose > 0) {
      cat("Calculating variance-covariance matrix at fit completion...\n")
    }
    vcov_cached <- NULL
    vcov_cached <- tryCatch({
      vcov.gamlss.longitudinal(
        return_list,
        numderiv = isTRUE(vcov_numderiv),
        method = vcov_method,
        progress = isTRUE(verbose > 0)
      )
    }, error = function(e) {
      warning(
        "Could not precompute variance-covariance matrix at fit completion: ",
        conditionMessage(e),
        call. = FALSE
      )
      NULL
    })

    if(!is.null(vcov_cached)) {
      return_list$vcov <- vcov_cached
      return_list$vcov_meta$precomputed <- TRUE
    }
  }

  class(return_list)="gamlss.longitudinal"
  return(return_list)
}

normalize_lag_time <- function(time) {
  if (is.factor(time)) {
    time <- as.character(time)
  }
  if (is.character(time)) {
    time_numeric <- suppressWarnings(as.numeric(time))
    if (anyNA(time_numeric)) {
      stop("ERROR: time must be numeric or numeric-like when use_dlcopdpar=TRUE.")
    }
    time <- time_numeric
  }
  time - 1
}

#' Create model matrices for model fitting
#'
#' This function takes the forumlas for each parameter mu,sigma,nu,tau,theta,zeta
#' and creates a list of model matrices mm with items mm$x and mm$s for
#' fixed and smooth terms respectively, with each of those lists being lists of each parameter
#' and their respective model matrices
#' @param mu.formula Formula for the mean parameter of the marginal distribution
#' @param sigma.formula Formula for the sigma parameter of the marginal distribution
#' @param nu.formula Formula for the nu parameter of the marginal distribution
#' @param tau.formula Formula for the tau parameter of the marginal distribution
#' @param theta.formula Formula for the theta parameter of the copula distribution
#' @param zeta.formula Formula for the zeta parameter of the copula distribution
#' @param margin.family Marginal distribution specified as a gamlss family object,
#' e.g. GA(), NO(), PO(), NBI(), etc.
#' @param copula.family Copula distribution code, one of "N", "C", "F", "G", "J", or "t".
#' @param copula.link List of link functions for the copula parameters
#' @return Returns a list mm with items mm$x and mm$s for fixed and smooth terms respectively,
#' with each of those lists being lists of each parameter and their respective model matrices
#'
#' @export
create_model_matrices<-function(
    mu.formula = ("response ~ 1"),
    sigma.formula = ("1"),
    nu.formula = ("1"),
    tau.formula = ("1"),
    theta.formula=("1"),
    zeta.formula=("1"),
    margin.family=NO(),
    copula.family="N",
    copula.link=NA,
    dataset=NA,
    quiet_gamlss2=TRUE
) {

  dataset_mm <- dataset
  normalize_ordered_factor <- function(col) {
    if (!is.factor(col)) return(col)
    if (!is.ordered(col)) return(col)

    levs <- levels(col)
    col_nom <- factor(as.character(col), levels = levs, ordered = FALSE)
    if (length(levs) > 1) {
      contr <- contr.treatment(length(levs))
      colnames(contr) <- levs[-1]
      contrasts(col_nom) <- contr
    }
    col_nom
  }

  mode_value <- function(x) {
    x_non_na <- x[!is.na(x)]
    if (length(x_non_na) == 0) return(NA)
    tab <- table(x_non_na)
    names(tab)[which.max(tab)]
  }

  # Build model matrices from an NA-free proxy dataset.
  # This does not alter likelihood calculations (which still use original dataset).
  for (nm in names(dataset_mm)) {
    if (!any(is.na(dataset_mm[[nm]]))) next
    if (nm %in% c("time", "subject")) next

    col <- dataset_mm[[nm]]
    if (is.numeric(col) || is.integer(col)) {
      obs <- col[!is.na(col)]
      fill_val <- if (length(obs) > 0) mean(obs) else 0
      col[is.na(col)] <- fill_val
      dataset_mm[[nm]] <- col
    } else if (is.factor(col)) {
      col <- normalize_ordered_factor(col)
      fill_val <- mode_value(col)
      if (is.na(fill_val)) {
        fill_val <- if (length(levels(col)) > 0) levels(col)[1] else "missing"
      }
      col_chr <- as.character(col)
      col_chr[is.na(col_chr)] <- fill_val
      dataset_mm[[nm]] <- factor(col_chr, levels = levels(col), ordered = FALSE)
    } else {
      fill_val <- mode_value(col)
      if (is.na(fill_val)) fill_val <- "missing"
      col[is.na(col)] <- fill_val
      dataset_mm[[nm]] <- col
    }
  }

  if ("response" %in% names(dataset_mm) && any(is.na(dataset_mm$response))) {
    obs_resp <- dataset_mm$response[!is.na(dataset_mm$response)]
    fill_val <- if (length(obs_resp) > 0) mean(obs_resp) else 0
    dataset_mm$response[is.na(dataset_mm$response)] <- fill_val
  }

  run_gamlss2 <- function(...) {
    if (isTRUE(quiet_gamlss2)) {
      fit <- NULL
      invisible(utils::capture.output({
        fit <- suppressMessages(suppressWarnings(gamlss2::gamlss2(...)))
      }, type = "output"))
      return(fit)
    }
    gamlss2::gamlss2(...)
  }

  normalize_time_covariate_colnames <- function(nms) {
    if (length(nms) == 0) return(nms)

    out <- nms
    suffix_map <- c(L = "1", Q = "2", C = "3")

    for (sx in names(suffix_map)) {
      out <- gsub(
        paste0("time_covariate\\.", sx, "\\b"),
        paste0("time_covariate.", suffix_map[[sx]]),
        out,
        perl = TRUE
      )
    }

    # contr.poly names can appear as ^4, ^5, ...; normalize to .4, .5, ...
    out <- gsub("time_covariate\\^([0-9]+)", "time_covariate.\\1", out, perl = TRUE)

    out
  }

  sanitize_for_gamlss2 <- function(data_in, fml) {
    vars_needed <- unique(all.vars(stats::as.formula(fml)))
    vars_needed <- vars_needed[vars_needed %in% names(data_in)]
    data_out <- data_in[, vars_needed, drop = FALSE]

    for (nm in names(data_out)) {
      col <- data_out[[nm]]
      if (is.factor(col)) {
        col <- droplevels(col)
        col <- normalize_ordered_factor(col)
        data_out[[nm]] <- col
      } else if (is.numeric(col) || is.integer(col)) {
        col[!is.finite(col)] <- NA
        if (any(is.na(col))) {
          obs <- col[!is.na(col)]
          fill_val <- if (length(obs) > 0) mean(obs) else 0
          col[is.na(col)] <- fill_val
        }
        data_out[[nm]] <- col
      } else {
        if (any(is.na(col))) {
          x_non_na <- col[!is.na(col)]
          fill_val <- if (length(x_non_na) > 0) {
            tab <- table(x_non_na)
            names(tab)[which.max(tab)]
          } else {
            "missing"
          }
          col[is.na(col)] <- fill_val
        }
        data_out[[nm]] <- col
      }
    }

    data_out
  }

  to_response_formula <- function(fml, response_name = "response") {
    if (inherits(fml, "formula")) {
      rhs_txt <- if (length(fml) == 3L) {
        paste(deparse(fml[[3]]), collapse = " ")
      } else {
        paste(deparse(fml[[2]]), collapse = " ")
      }
    } else if (is.character(fml) && length(fml) == 1L) {
      txt <- trimws(fml)
      if (grepl("~", txt, fixed = TRUE)) {
        parts <- strsplit(txt, "~", fixed = TRUE)[[1]]
        rhs_txt <- trimws(parts[length(parts)])
      } else {
        rhs_txt <- txt
      }
    } else {
      stop("Invalid formula input: ", deparse(fml))
    }

    as.formula(paste(response_name, "~", rhs_txt), env = parent.frame())
  }

  if(copula.family %in% c("t", "T", "Student")){two_par_cop=TRUE} else {two_par_cop=FALSE}
  included_parameters <- c(names(margin.family$parameters), if(two_par_cop) c("theta","zeta") else c("theta"))

  formulas=list()
  for (parameter in included_parameters) {
    formulas[[parameter]]=get(paste(parameter,"formula",sep="."))
  }

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required to construct smooth-term model matrices.")
  }

  formulas[["mu"]] <- as.formula(mu.formula)
  for (parameter in included_parameters[2:length(included_parameters)]) {
    formulas[[parameter]] <- to_response_formula(formulas[[parameter]], response_name = "response")
  }

  mm_x=list()
  mm_s=list()

  smooth_eval_env <- new.env(parent = baseenv())
  smooth_eval_env$s <- mgcv::s

  for(parameter in included_parameters) {
    data_for_par <- if(parameter %in% c("theta","zeta")) {
      dataset_mm[dataset_mm$time %in% unique(dataset_mm$time)[1:(length(unique(dataset_mm$time))-1)], , drop = FALSE]
    } else {
      dataset_mm
    }

    data_for_par <- sanitize_for_gamlss2(data_for_par, formulas[[parameter]])
    formula_terms <- stats::terms(formulas[[parameter]])
    has_intercept <- as.integer(attr(formula_terms, "intercept")) == 1L
    term_labels <- attr(formula_terms, "term.labels")

    # Keep all non-smooth RHS terms so model.matrix can expand interactions
    # like `time*gender` into main effects + interaction columns.
    fixed_terms <- term_labels[!grepl("^\\s*s\\(", term_labels)]

    if(length(fixed_terms) > 0 || has_intercept) {
      fixed_formula <- if (length(fixed_terms) == 0L && has_intercept) {
        stats::as.formula("~ 1")
      } else {
        stats::reformulate(termlabels = fixed_terms, intercept = has_intercept)
      }
      X_fixed <- stats::model.matrix(fixed_formula, data = data_for_par)
      colnames(X_fixed) <- sub("^\\(Intercept\\)$", "intercept", colnames(X_fixed))
      colnames(X_fixed) <- normalize_time_covariate_colnames(colnames(X_fixed))
      mm_x[[parameter]] <- as.data.frame(X_fixed, check.names = FALSE)
    } else {
      mm_x[[parameter]] <- data.frame(row.names = seq_len(nrow(data_for_par)))
    }

    smooth_terms <- term_labels[grepl("^\\s*s\\(", term_labels)]
    if(length(smooth_terms) == 0) {
      mm_s[[parameter]] <- NULL
    } else {
      mm_s[[parameter]] <- list()
      for (s_label in smooth_terms) {
        s_txt <- trimws(s_label)
        s_call <- tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)
        s_obj <- eval(parse(text = s_txt), envir = smooth_eval_env)
        s_con <- mgcv::smoothCon(s_obj, data = data_for_par, knots = NULL)
        if (length(s_con) > 0 && !is.null(s_con[[1]]$X)) {
          B_s <- s_con[[1]]$X
          # Store the basis-specific penalty matrix returned by smoothCon so the
          # optimizer can use it instead of a generic second-difference fallback.
          if (!is.null(s_con[[1]]$S) && length(s_con[[1]]$S) > 0) {
            attr(B_s, "penalty") <- s_con[[1]]$S[[1]]
          }
          if (!is.null(s_call) && length(s_call) >= 2) {
            x_expr <- s_call[[2]]
            x_var <- trimws(gsub("`", "", paste(deparse(x_expr), collapse = " "), fixed = TRUE))
            x_value <- tryCatch(eval(x_expr, envir = data_for_par, enclos = parent.frame()), error = function(e) NULL)
            if (!is.null(x_value) && length(x_value) == nrow(B_s)) {
              attr(B_s, "smooth_x") <- as.numeric(x_value)
              attr(B_s, "smooth_var") <- x_var
            } else if (nzchar(x_var)) {
              attr(B_s, "smooth_var") <- x_var
            }
          }
          mm_s[[parameter]][[s_txt]] <- B_s
        }
      }
      if(length(mm_s[[parameter]]) == 0) {
        mm_s[[parameter]] <- NULL
      }
    }
  }


  mm=list(mm_x,mm_s)
  names(mm)=c("x","s")
  return(mm)
}

#' Calculate eta, eta inverse and eta derivative based on the given parameters and model matrices
#'
#' This function calculates the linear predictors (eta) for each parameter
#' based on the given covariate parameters (par_cov) and model matrices (mm).
#' It also computes the inverse link function (eta_inv) and the derivative of the link function (eta_dr)
#' for each parameter using the specified marginal distribution and copula link functions.
#'
#' @param par_cov A named vector of covariate parameters for each model term.
#' @param mm A list containing model matrices for fixed effects (mm$x) and smooth terms (mm$s).
#' @param margin_dist A list of functions for the marginal distribution, including link inverse and derivative functions.
#' @param copula_link A list of functions for the copula link, including link inverse and derivative functions.
#' @param par_s A list of smooth term parameters for each model parameter (optional).
#'
#' @return A list containing:
#' \item{eta}{A list of linear predictors for each parameter.}
#' \item{eta_inv}{A list of inverse link function values for each parameter.}
#' \item{eta_dr}{A list of derivatives of the link function for each parameter.}
#'
#' @export
calc_eta=function(par_cov,mm,margin_dist,copula_link,par_s=NA) {
  eta=list()
  #par_s=list(); par_s[["mu"]]=list(); par_s[["sigma"]]=list(); par_s[["nu"]]=list(); par_s[["tau"]]=list()
  #par_s[["mu"]][["s(age)"]]=matrix(rep(1,ncol(mm$s[["mu"]][["s(age)"]])),ncol=1)
  for (par_name in names(mm$x)) {
    par_cov_single=par_cov[grepl(par_name,names(par_cov))]
    mm_temp=mm$x[[par_name]]
    #If there are no smooth terms for the parameter then just do standard calculation
    if(all(is.na(par_s[[par_name]]))) {
      eta[[par_name]]=rowSums(mm_temp * matrix(rep(par_cov_single,each=nrow(mm_temp)),ncol=length(par_cov_single),dimnames=list(NULL,c(names(par_cov_single)))))
    } else {
      eta[[par_name]]=
        rowSums(mm_temp * matrix(rep(par_cov_single,each=nrow(mm_temp)),ncol=length(par_cov_single),dimnames=list(NULL,c(names(par_cov_single)))))
      for (s_name in names(mm$s[[par_name]])) {
        eta[[par_name]]=eta[[par_name]] + mm$s[[par_name]][[s_name]] %*% par_s[[par_name]][[s_name]]
      }
    }
  }
  #Get link transforms (eta) and derivatives w.r.t to link for parameters
  eta_dr=eta_inv=list()
  for (par_name in names(mm$x)) {
    if(par_name %in% c("mu","sigma","nu","tau")) {
      eta_inv[[par_name]]=margin_dist[[paste(par_name,".linkinv",sep="")]](eta[[par_name]])
      eta_dr[[par_name]]=margin_dist[[paste(par_name,".dr",sep="")]](eta[[par_name]])
    }
    if(par_name %in% c("theta","zeta")) {
      eta_inv[[par_name]]=copula_link[[paste(par_name,".linkinv",sep="")]](eta[[par_name]])
      eta_dr[[par_name]]=copula_link[[paste(par_name,".dr",sep="")]](eta[[par_name]])
    }
  }
  return(list(eta=eta,eta_inv=eta_inv,eta_dr=eta_dr))
}

.calc_eta_rs_cached <- function(
  rs_design_cache,
  par_cov,
  par_s,
  margin_dist,
  copula_link,
  update_only = NULL,
  eta_out = NULL
) {
  if (is.null(eta_out)) {
    eta_out <- list(eta = list(), eta_inv = list(), eta_dr = list())
  }

  par_names <- names(rs_design_cache)
  if (!is.null(update_only)) {
    par_names <- intersect(update_only, par_names)
  }

  for (par_name in par_names) {
    design_info <- rs_design_cache[[par_name]]
    X <- design_info$X
    beta <- numeric(ncol(X))
    names(beta) <- colnames(X)

    fixed_names <- design_info$fixed_names
    beta[fixed_names] <- par_cov[fixed_names]

    if (length(par_s[[par_name]]) > 0) {
      smooth_beta <- unlist(par_s[[par_name]], use.names = FALSE)
      smooth_names <- setdiff(colnames(X), fixed_names)
      if (length(smooth_beta) != length(smooth_names)) {
        stop("Smooth coefficient length does not match cached design columns for ", par_name, ".", call. = FALSE)
      }
      beta[smooth_names] <- smooth_beta
    }

    eta_vec <- as.numeric(X %*% beta)
    eta_out$eta[[par_name]] <- eta_vec

    if (par_name %in% c("mu", "sigma", "nu", "tau")) {
      eta_out$eta_inv[[par_name]] <- margin_dist[[paste(par_name, ".linkinv", sep = "")]](eta_vec)
      eta_out$eta_dr[[par_name]] <- margin_dist[[paste(par_name, ".dr", sep = "")]](eta_vec)
    } else if (par_name %in% c("theta", "zeta")) {
      eta_out$eta_inv[[par_name]] <- copula_link[[paste(par_name, ".linkinv", sep = "")]](eta_vec)
      eta_out$eta_dr[[par_name]] <- copula_link[[paste(par_name, ".dr", sep = "")]](eta_vec)
    }
  }

  eta_out
}

#' Calculate the likelihood components for the joint model
#'
#' This function calculates the marginal and copula log likelihoods and components
#' for the joint model by organizing response data by margin and subject for
#' efficient pair-based copula calculations.
#'
#' @param response A numeric vector of response values.
#' @param response_margin A numeric vector indicating the margin (time) for each response.
#' @param response_subject A numeric vector indicating the subject for each response.
#' @return A list containing:
#' \item{log_lik}{A named vector with marginal, copula, and joint log-likelihoods.}
#' \item{margin_d}{A numeric vector of marginal densities.}
#' \item{copula_d}{A numeric vector of copula densities.}
#' \item{margin_p}{A numeric vector of marginal distribution function values.}
#' \item{Fx_1_2}{A matrix of marginal distribution function values for pairs of margins.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#' \item{margin_deriv}{A list of marginal derivatives.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#'
#' @export
build_copula_pair_cache <- function(response, response_margin, response_subject) {
  margin_names=sort(unique(response_margin))
  num_margins=length(margin_names)
  n_obs=length(response)
  obs_response=!is.na(response)

  base_df=data.frame(
    row_id=seq_len(n_obs),
    time=response_margin,
    subject=response_subject,
    observed=obs_response,
    stringsAsFactors = FALSE
  )

  pair_df_all=list()
  if(num_margins>1) {
    for (i in seq_len(num_margins-1)) {
      t1=margin_names[i]
      t2=margin_names[i+1]

      left=base_df[base_df$time==t1,c("row_id","subject","time","observed")]
      right=base_df[base_df$time==t2,c("row_id","subject","time","observed")]
      names(left)=c("row_id1","subject","time1","observed1")
      names(right)=c("row_id2","subject","time2","observed2")

      pair_i=merge(left,right,by="subject",all=FALSE)
      if(nrow(pair_i)>0) {
        pair_df_all[[length(pair_df_all)+1]]=pair_i
      }
    }
  }

  if(length(pair_df_all)==0) {
    pair_df=data.frame(
      subject=response_subject[0],
      row_id1=integer(0),
      time1=response_margin[0],
      observed1=logical(0),
      row_id2=integer(0),
      time2=response_margin[0],
      observed2=logical(0)
    )
  } else {
    pair_df=do.call(rbind,pair_df_all)
  }

  order_copula=as.matrix(pair_df[,c("time1","subject","time2","subject")])
  colnames(order_copula)=c("time1","subject1","time2","subject2")

  observed_pair_base=rep(FALSE,nrow(pair_df))
  if(nrow(pair_df)>0) {
    observed_pair_base=pair_df$observed1 & pair_df$observed2
  }

  Fx_1_2_template=matrix(NA_real_,nrow=nrow(pair_df),ncol=2)
  colnames(Fx_1_2_template)=c("u1","u2")

  theta_rows=which(response_margin %in% margin_names[seq_len(max(1, num_margins-1))])
  theta_index_map=rep(NA_integer_,n_obs)
  theta_index_map[theta_rows]=seq_along(theta_rows)

  cache=list(
    row_id1=pair_df$row_id1,
    row_id2=pair_df$row_id2,
    Fx_1_2_template=Fx_1_2_template,
    order_copula=order_copula,
    observed_pair_base=observed_pair_base,
    theta_index_map=theta_index_map,
    margin_names=margin_names,
    num_margins=num_margins,
    n_obs=n_obs
  )
  cache
}

.build_margin_eval_cache <- function(margin_dist, calc_d2=FALSE) {
  if(calc_d2==TRUE) {
    to_include=grepl("dld",names(margin_dist))|grepl("d2ld",names(margin_dist))
  } else {
    to_include=grepl("dld",names(margin_dist))
  }

  margin_deriv_names=names(margin_dist)[to_include]
  margin_deriv_cache=lapply(margin_deriv_names, function(deriv_name) {
    FUN=margin_dist[[deriv_name]]
    list(name=deriv_name,FUN=FUN,args=formalArgs(FUN))
  })

  margin_pFUN=get(paste("p",margin_dist$family[1],sep=""),mode="function",inherits=TRUE)
  margin_dFUN=get(paste("d",margin_dist$family[1],sep=""),mode="function",inherits=TRUE)

  list(
    calc_d2=calc_d2,
    margin_deriv_cache=margin_deriv_cache,
    margin_pFUN=margin_pFUN,
    margin_p_args=formalArgs(margin_pFUN),
    margin_dFUN=margin_dFUN,
    margin_d_args=formalArgs(margin_dFUN)
  )
}

.bcpe_FT <- function(t, tau, log_c=NULL) {
  if(is.null(log_c)) {
    log_c=0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
  }
  c_val=exp(log_c)
  s=0.5*((abs(t/c_val))^tau)
  F_s=pgamma(s, shape=1/tau, scale=1)
  0.5*(1+F_s*sign(t))
}

.bcpe_fT_log <- function(t, tau, log_c=NULL) {
  if(is.null(log_c)) {
    log_c=0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
  }
  c_val=exp(log_c)
  log(tau)-log_c-(0.5*(abs(t/c_val)^tau))-(1+(1/tau))*log(2)-lgamma(1/tau)
}

.eval_bcpe_margin_pd <- function(y, mu, sigma, nu, tau) {
  if(any(mu < 0)) stop(paste("mu must be positive", "\n", ""))
  if(any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if(any(tau < 0)) stop(paste("tau must be positive", "\n", ""))

  z=if(length(nu)>1) {
    ifelse(nu != 0, (((y/mu)^nu-1)/(nu*sigma)), log(y/mu)/sigma)
  } else if(nu != 0) {
    (((y/mu)^nu-1)/(nu*sigma))
  } else {
    log(y/mu)/sigma
  }

  log_c=0.5*(-(2/tau)*log(2)+lgamma(1/tau)-lgamma(3/tau))
  F_z=.bcpe_FT(z, tau, log_c=log_c)
  F_upper=.bcpe_FT(1/(sigma*abs(nu)), tau, log_c=log_c)
  F_lower=if(length(nu)>1) {
    ifelse(nu > 0, .bcpe_FT(-1/(sigma*abs(nu)), tau, log_c=log_c), 0)
  } else if(nu > 0) {
    .bcpe_FT(-1/(sigma*abs(nu)), tau, log_c=log_c)
  } else {
    0
  }

  margin_p=(F_z-F_lower)/F_upper
  log_fz=.bcpe_fT_log(z, tau, log_c=log_c)-log(F_upper)
  log_der=(nu-1)*log(y)-nu*log(mu)-log(sigma)
  margin_d=exp(log_der+log_fz)
  margin_d=ifelse(y <= 0, 0, margin_d)

  list(p=margin_p, d=margin_d)
}

.eval_margin_pd <- function(margin_deriv_input, margin_eval_cache) {
  if(identical(margin_eval_cache$family, "BCPE") &&
     all(c("y", "mu", "sigma", "nu", "tau") %in% names(margin_deriv_input))) {
    return(.eval_bcpe_margin_pd(
      y=margin_deriv_input[["y"]],
      mu=margin_deriv_input[["mu"]],
      sigma=margin_deriv_input[["sigma"]],
      nu=margin_deriv_input[["nu"]],
      tau=margin_deriv_input[["tau"]]
    ))
  }

  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%margin_eval_cache$margin_p_args]
  margin_p=do.call(margin_eval_cache$margin_pFUN,args=margin_deriv_input[FUN_args])

  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%margin_eval_cache$margin_d_args]
  margin_d=do.call(margin_eval_cache$margin_dFUN,args=margin_deriv_input[FUN_args])

  list(p=margin_p, d=margin_d)
}

#' @param pair_cache Optional cache built by build_copula_pair_cache to reuse pair indexing across repeated likelihood calls.
calc_likelihood_minimal <- function(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject,penalize_smooth=FALSE,par_s=NA,pair_cache=NULL,margin_eval_cache=NULL) {
  #Setup input matrix of response and parameters
  #response=dataset$response; response_subject=dataset$subject; response_margin=dataset$time; dataset=NA
  if(is.null(pair_cache)) {
    pair_cache=build_copula_pair_cache(response,response_margin,response_subject)
  }
  if(is.null(margin_eval_cache) || !identical(margin_eval_cache$calc_d2, calc_d2)) {
    margin_eval_cache=.build_margin_eval_cache(margin_dist, calc_d2=calc_d2)
  }

  margin_names=pair_cache$margin_names
  num_margins=pair_cache$num_margins
  n_obs=pair_cache$n_obs

  obs_response=!is.na(response)

  order_margin=cbind(response_margin,response_subject)
  colnames(order_margin)=c("time","subject")

  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu","sigma","nu","tau")) {
      margin_deriv_input[[par_name]]=eta_inv[[par_name]]
    }
  }

  #Calculate all derivatives

  ################## MARGIN DERIVATIVES
  margin_deriv=list()
  for (deriv_info in margin_eval_cache$margin_deriv_cache) {
    FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%deriv_info$args]
    deriv_val=do.call(deriv_info$FUN,args=margin_deriv_input[FUN_args])
    if(length(deriv_val)==n_obs) {
      deriv_val[!obs_response]=0
      deriv_val[!is.finite(deriv_val)]=0
    }
    margin_deriv[[deriv_info$name]]=deriv_val
  }

  margin_pd=.eval_margin_pd(margin_deriv_input, margin_eval_cache)
  margin_p=margin_pd$p
  margin_p[!obs_response]=NA
  margin_p[!is.finite(margin_p)]=NA

  margin_d=margin_pd$d
  margin_d[!obs_response]=NA
  margin_d[!is.finite(margin_d) | margin_d<=0]=NA

  ################COPULA DERIVATIVES
  #First calculate margin F(x1), F(x2) as inputs to copula

  row_id1=pair_cache$row_id1
  row_id2=pair_cache$row_id2
  order_copula=pair_cache$order_copula

  Fx_1_2=pair_cache$Fx_1_2_template
  if(is.null(Fx_1_2)) {
    Fx_1_2=matrix(NA_real_,nrow=length(row_id1),ncol=2)
    colnames(Fx_1_2)=c("u1","u2")
  }
  if(length(row_id1)>0) {
    Fx_1_2[,1]=margin_p[row_id1]
    Fx_1_2[,2]=margin_p[row_id2]
  }

  pair_complete=pair_cache$observed_pair_base & is.finite(Fx_1_2[,1]) & is.finite(Fx_1_2[,2])

  par1=rep(NA_real_,length(row_id1))
  par2=rep(NA_real_,length(row_id1))
  if(length(row_id1)>0) {
    theta_len <- length(eta_inv[["theta"]])
    if(theta_len == n_obs) {
      theta_idx <- row_id1
    } else {
      theta_idx=pair_cache$theta_index_map[row_id1]
    }

    par1=eta_inv[["theta"]][theta_idx]
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]][theta_idx]
    } else {
      par2=rep(0,length(par1))
    }
  }

  pair_complete=pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval=Fx_1_2
  if(nrow(Fx_eval)>0) {
    Fx_eval[!is.finite(Fx_eval)]=0.5
    Fx_eval[Fx_eval>1]=1
    Fx_eval[Fx_eval<0]=0
  }

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  if(length(par1_eval)==0) {
    copula_d=numeric(0)
  } else {
    copula_d=.copula_pdf(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval)
  }
  if(length(copula_d)>0) {
    copula_d[!is.finite(copula_d) | copula_d<=0]=1
    copula_d[!pair_complete]=1
  }

  ########COMBINE MARGINS AND COPULA DERVIATIVES

  margin_loglik_terms=log(margin_d[!is.na(margin_d)])
  margin_loglik_terms=margin_loglik_terms[is.finite(margin_loglik_terms)]
  copula_loglik_terms=log(copula_d[pair_complete])
  copula_loglik_terms=copula_loglik_terms[is.finite(copula_loglik_terms)]

  log_lik=c(sum(margin_loglik_terms),sum(copula_loglik_terms),sum(margin_loglik_terms)+sum(copula_loglik_terms))
  names(log_lik)=c("marginal","copula","joint")

  copula_p=rep(NA_real_,length(copula_d))

  return_list=list(log_lik,margin_d,copula_d,margin_p,copula_p,Fx_1_2,order_copula,margin_deriv,pair_complete,par1,par2,row_id1,row_id2,pair_cache$theta_index_map)
  names(return_list)=c("log_lik","margin_d","copula_d","margin_p","copula_p","Fx_1_2","order_copula","margin_deriv","pair_complete","copula_par1","copula_par2","copula_row_id1","copula_row_id2","copula_theta_index_map")
  return(return_list)
}

.calc_likelihood_update_copula <- function(eta_inv, base_lik, copula_dist, pair_cache) {
  row_id1 <- pair_cache$row_id1
  Fx_1_2 <- base_lik$Fx_1_2
  n_obs <- pair_cache$n_obs

  pair_complete <- pair_cache$observed_pair_base & is.finite(Fx_1_2[, 1]) & is.finite(Fx_1_2[, 2])

  par1 <- rep(NA_real_, length(row_id1))
  par2 <- rep(NA_real_, length(row_id1))
  if(length(row_id1) > 0) {
    theta_len <- length(eta_inv[["theta"]])
    if(theta_len == n_obs) {
      theta_idx <- row_id1
    } else {
      theta_idx <- pair_cache$theta_index_map[row_id1]
    }

    par1 <- eta_inv[["theta"]][theta_idx]
    if("zeta" %in% names(eta_inv)) {
      par2 <- eta_inv[["zeta"]][theta_idx]
    } else {
      par2 <- rep(0, length(par1))
    }
  }

  pair_complete <- pair_complete & is.finite(par1) & is.finite(par2)

  Fx_eval <- Fx_1_2
  if(nrow(Fx_eval) > 0) {
    Fx_eval[!is.finite(Fx_eval)] <- 0.5
    Fx_eval[Fx_eval > 1] <- 1
    Fx_eval[Fx_eval < 0] <- 0
  }

  par1_eval <- par1
  par2_eval <- par2
  par1_eval[!is.finite(par1_eval)] <- 0
  par2_eval[!is.finite(par2_eval)] <- 0

  if(copula_dist == "C") {
    par1_eval[par1_eval >= 28] <- 27.9
  }

  if(length(par1_eval) == 0) {
    copula_d <- numeric(0)
  } else {
    copula_d <- .copula_pdf(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval)
  }
  if(length(copula_d) > 0) {
    copula_d[!is.finite(copula_d) | copula_d <= 0] <- 1
    copula_d[!pair_complete] <- 1
  }

  copula_loglik_terms <- log(copula_d[pair_complete])
  copula_loglik_terms <- copula_loglik_terms[is.finite(copula_loglik_terms)]
  copula_loglik <- sum(copula_loglik_terms)
  marginal_loglik <- as.numeric(base_lik$log_lik["marginal"])
  log_lik <- c(marginal_loglik, copula_loglik, marginal_loglik + copula_loglik)
  names(log_lik) <- c("marginal", "copula", "joint")

  base_lik$log_lik <- log_lik
  base_lik$copula_d <- copula_d
  base_lik$pair_complete <- pair_complete
  base_lik$copula_par1 <- par1
  base_lik$copula_par2 <- par2
  base_lik
}

#'
#'
#' @export
score_function_v2 <- function(eta,dldpar,d2ldpar,dpardeta,response=NA,phi=1,step_size=1,verbose=FALSE,crit_wk=0.0000001) {

  u_k=dldeta = dldpar * dpardeta
  f_k=d2ldpar
  w_k=-f_k*(dpardeta*dpardeta)

  #Stop if weights are too small
  w_k[abs(w_k)<crit_wk]=1
  u_k[abs(w_k)<crit_wk]=0

  w_k[abs(u_k)<crit_wk]=1
  u_k[abs(u_k)<crit_wk]=0

  z_k=(1-phi)*eta+phi*(eta+step_size*(u_k/w_k))

  if(verbose==TRUE) {
    steps_mean=round(rbind(colMeans(as.matrix(eta))
                           ,colMeans(as.matrix(dldpar-dlcopdpar))
                           ,colMeans(as.matrix(dlcopdpar))
                           ,colMeans(as.matrix(dpardeta))
                           ,colMeans(as.matrix(dpardeta*dpardeta))
                           ,colMeans(as.matrix(f_k))
                           ,colMeans(as.matrix(w_k))
                           ,colMeans(as.matrix(u_k))
                           ,colMeans(as.matrix(u_k/w_k))
                           ,colMeans(as.matrix(z_k))
    ),8)
    rownames(steps_mean)=c("eta","dldpar","dlcopdpar","dpardeta","dpardeta2","f_k","w_k","u_k","(1/w_k)*u_k","z_k")
    print(steps_mean)
  }
  return_list=list(colMeans(as.matrix(z_k)),as.matrix(u_k),as.matrix(f_k),as.matrix(w_k),as.matrix(z_k))
  names(return_list)=c("par","u_k","f_k","w_k","z_k")
  return(return_list)
}

calc_deriv_copula_wrt_margin = function(input,margin_par,par_name,calc_d2=FALSE) {

    #Calculate copula derivative with respect to marginal parameters
    #input=copula_merged
    num_col <- function(nm) {
      val <- if (is.data.frame(input)) input[[nm]] else input[, nm]
      suppressWarnings(as.numeric(val))
    }

    if(calc_d2==FALSE) {

      dlcopdpar_row1=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      dlcopdpar_row2=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=num_col("dcdu1")
          dc_tplus_du_tplus=num_col("dcdu2")
          l_t=num_col(paste(paste("dld",inner_par_name,sep=""),".x",sep=""))
          l_t_plus=num_col(paste(paste("dld",inner_par_name,sep=""),".y",sep=""))
          x_t=num_col("response.x")
          x_t_plus=num_col("response.y")
          f_t=num_col("margin_d.x")
          f_t_plus=num_col("margin_d.y")
          c_tplus=num_col("copula_d")
          mu_t=num_col("mu.x")
          mu_t_plus=num_col("mu.y")

          F_nd_t=num_col("F_nd.x")
          F_nd_t_plus=num_col("F_nd.y")

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          # Exact endpoint attribution for pair log-copula derivative:
          # row_id1 gets (dc/du1)*(du1/dpar)/c, row_id2 gets (dc/du2)*(du2/dpar)/c.
          dlogc_row1=(dc_tplus_du_t * du_t_dmu) / c_tplus
          dlogc_row2=(dc_tplus_du_tplus * du_t_plus_dmu) / c_tplus
          dlogc_row1[!is.finite(dlogc_row1)] = 0
          dlogc_row2[!is.finite(dlogc_row2)] = 0

          dlcopdpar_row1[,i]=dlogc_row1
          dlcopdpar_row2[,i]=dlogc_row2

        }
        i=i+1
      }
      colnames(dlcopdpar_row1)=paste("dlcopd",margin_par,sep="")
      colnames(dlcopdpar_row2)=paste("dlcopd",margin_par,sep="")

      # Prefer explicit row-index accumulation to avoid merge-order instability.
      if(all(c("row_id1","row_id2") %in% colnames(input))) {
        n_obs <- suppressWarnings(max(c(num_col("row_id1"), num_col("row_id2")), na.rm = TRUE))
        if(!is.finite(n_obs) || n_obs < 1) {
          stop("Invalid row ids in copula-to-margin derivative assembly.")
        }
        d1_cop <- matrix(0, nrow = as.integer(n_obs), ncol = length(margin_par))
        colnames(d1_cop) <- margin_par

        row_id1 <- as.integer(num_col("row_id1"))
        row_id2 <- as.integer(num_col("row_id2"))
        for(j in seq_along(margin_par)) {
          contrib1 <- as.numeric(dlcopdpar_row1[, j])
          contrib2 <- as.numeric(dlcopdpar_row2[, j])
          valid1 <- is.finite(contrib1) & is.finite(row_id1) & row_id1 >= 1 & row_id1 <= n_obs
          valid2 <- is.finite(contrib2) & is.finite(row_id2) & row_id2 >= 1 & row_id2 <= n_obs
          d1_cop[row_id1[valid1], j] <- d1_cop[row_id1[valid1], j] + contrib1[valid1]
          d1_cop[row_id2[valid2], j] <- d1_cop[row_id2[valid2], j] + contrib2[valid2]
        }
        return(d1_cop)
      }

      dlcopdpar <- dlcopdpar_row1 + dlcopdpar_row2

      par_dlcopdpar=dlcopdpar[,paste("dlcopd",margin_par,sep="")]
      merged_dlcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),cbind(input[,c("time1","time2","subject1","subject2")],par_dlcopdpar),by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_dlcopdpar[is.na(merged_dlcopdpar)]=0

      x_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".x",colnames(merged_dlcopdpar))
      y_comp=grepl("dlcopd",colnames(merged_dlcopdpar))&grepl(".y",colnames(merged_dlcopdpar))

      d1_cop=0.5*(merged_dlcopdpar[,x_comp]+merged_dlcopdpar[,y_comp])

      return(d1_cop)

    } else {

      d2lcopdpar2=matrix(0,nrow=nrow(input),ncol=length(margin_par))
      i=1
      for (inner_par_name in margin_par) {

        if(inner_par_name==par_name) {
          #Take parameters from input for clarity
          dc_tplus_du_t=num_col("dcdu1")
          dc_tplus_du_tplus=num_col("dcdu2")
          #l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
          #l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
          #x_t=input[,"response.x"]
          #x_t_plus=input[,"response.y"]
          #f_t=input[,"margin_d.x"]
          #f_t_plus=input[,"margin_d.y"]
          c_tplus=num_col("copula_d")
          mu_t=num_col("mu.x")
          mu_t_plus=num_col("mu.y")

          F_nd_t=num_col("F_nd.x")
          F_nd_t_plus=num_col("F_nd.y")

          du_t_dmu=F_nd_t
          du_t_plus_dmu=F_nd_t_plus

          dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
          dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
          dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
          dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
          dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
          dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

          #dlcopdpar[,i]=dcdmu_tplus

          #######NOW FOR SECOND DERIVATIVE OF COPULA TERM

          F_nd2=num_col("F_nd2.x")
          F_nd2_plus=num_col("F_nd2.y")

          d2u_t_dmu2=F_nd2
          d2u_t_plus_dmu2=F_nd2_plus

          d2cdu_t2=num_col("d2cdu12")
          d2cdu_t_plus2=num_col("d2cdu22")
          d2cdu_t2[is.nan(d2cdu_t2)]=0
          d2cdu_t_plus2[is.nan(d2cdu_t_plus2)]=0

          d2cdmu2=  d2cdu_t2*du_t_dmu^2 +
                    dc_tplus_du_t * d2u_t_dmu2 +
                    d2cdu_t_plus2*du_t_plus_dmu^2 +
                    dc_tplus_du_tplus * d2u_t_plus_dmu2

          d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))
          d2lcdmu2=num_col("c_nd2")

          d2lcopdpar2[,i]=d2lcdmu2

        }
        i=i+1
      }
      colnames(d2lcopdpar2)=paste("d2lcopd",margin_par,sep="")

      par_d2lcopdpar=d2lcopdpar2[,paste("d2lcopd",margin_par,sep="")]
      merged_d2lcopdpar=merge(cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,cbind(input[,c("time1","time2","subject1","subject2")],par_d2lcopdpar)
                              ,by.x=c("time2","subject2"),by.y=c("time1","subject1"),all=TRUE)
      merged_d2lcopdpar[is.na(merged_d2lcopdpar)]=0

      x_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".x",colnames(merged_d2lcopdpar))
      y_comp=grepl("d2lcopd",colnames(merged_d2lcopdpar))&grepl(".y",colnames(merged_d2lcopdpar))

      d2_cop=0.5*(merged_d2lcopdpar[,x_comp]+merged_d2lcopdpar[,y_comp])

      #plot(d2lcopdpar2[,paste("d2lcopd",par_name,sep="")],input[,"c_nd2"],main="d2",ylab="numerical")

      return(d2_cop)
    }


}

.calc_dlcopdpar_indexed <- function(
  row_id1,
  row_id2,
  dcdu1,
  dcdu2,
  copula_d,
  F_nd,
  n_obs,
  pair_complete = NULL
) {
  row_id1 <- as.integer(row_id1)
  row_id2 <- as.integer(row_id2)
  n_pair <- length(row_id1)

  if (length(row_id2) != n_pair || length(dcdu1) != n_pair || length(dcdu2) != n_pair ||
      length(copula_d) != n_pair) {
    stop("Copula derivative inputs have inconsistent pair lengths.", call. = FALSE)
  }
  if (length(F_nd) != n_obs) {
    stop("F derivative length does not match the number of observations.", call. = FALSE)
  }

  if (is.null(pair_complete)) {
    pair_complete <- rep(TRUE, n_pair)
  } else {
    pair_complete <- as.logical(pair_complete)
    if (length(pair_complete) != n_pair) {
      stop("pair_complete length does not match copula pair length.", call. = FALSE)
    }
  }

  dlogc_row1 <- (as.numeric(dcdu1) * as.numeric(F_nd[row_id1])) / as.numeric(copula_d)
  dlogc_row2 <- (as.numeric(dcdu2) * as.numeric(F_nd[row_id2])) / as.numeric(copula_d)
  dlogc_row1[!pair_complete | !is.finite(dlogc_row1)] <- 0
  dlogc_row2[!pair_complete | !is.finite(dlogc_row2)] <- 0

  out <- numeric(n_obs)
  valid1 <- is.finite(row_id1) & row_id1 >= 1L & row_id1 <= n_obs
  valid2 <- is.finite(row_id2) & row_id2 >= 1L & row_id2 <= n_obs
  if (any(valid1)) {
    sum1 <- rowsum(dlogc_row1[valid1], row_id1[valid1], reorder = FALSE)
    out[as.integer(rownames(sum1))] <- out[as.integer(rownames(sum1))] + sum1[, 1]
  }
  if (any(valid2)) {
    sum2 <- rowsum(dlogc_row2[valid2], row_id2[valid2], reorder = FALSE)
    out[as.integer(rownames(sum2))] <- out[as.integer(rownames(sum2))] + sum2[, 1]
  }
  out
}

check_dlcopdpar_gradient_margin_score <- function(
  eta,
  eta_inv,
  par_name,
  margin_dist,
  copula_dist,
  dataset,
  mm,
  pair_cache,
  d1,
  base_loglik,
  verbose=FALSE
) {
  if (!par_name %in% c("mu", "sigma", "nu", "tau")) {
    return(list(warned = FALSE, message = NULL))
  }

  eta_vec <- as.numeric(eta[[par_name]])
  score_vec <- as.numeric(d1)
  finite_idx <- which(is.finite(eta_vec) & is.finite(score_vec))
  if (length(finite_idx) == 0) {
    return(list(warned = FALSE, message = NULL))
  }

  probe_candidates <- unique(round(seq(1, length(finite_idx), length.out = min(3, length(finite_idx)))))
  probe_idx <- finite_idx[probe_candidates]
  eps <- 1e-6
  tolerance <- 1e-3
  diffs <- numeric(length(probe_idx))

  linkinv_fun <- margin_dist[[paste(par_name, ".linkinv", sep = "")]]
  if (is.null(linkinv_fun)) {
    return(list(warned = FALSE, message = NULL))
  }

  for (k in seq_along(probe_idx)) {
    idx <- probe_idx[k]
    eta_pert <- eta
    eta_pert[[par_name]][idx] <- eta_pert[[par_name]][idx] + eps

    eta_inv_pert <- eta_inv
    eta_inv_pert[[par_name]] <- linkinv_fun(eta_pert[[par_name]])

    lik_pert <- calc_likelihood_minimal(
      eta_inv_pert,
      mm = mm,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      calc_d2 = FALSE,
      response = dataset$response,
      response_margin = dataset$time,
      response_subject = dataset$subject,
      pair_cache = pair_cache
    )$log_lik["joint"]

    finite_diff <- (lik_pert - base_loglik) / eps
    diffs[k] <- finite_diff - score_vec[idx]
  }

  max_abs_diff <- max(abs(diffs), na.rm = TRUE)
  if (!is.finite(max_abs_diff)) {
    return(list(warned = FALSE, message = NULL))
  }

  message_text <- paste0(
    "DLCOPDGRAD check for ", par_name,
    ": max abs(score - finite_diff) = ", signif(max_abs_diff, 4),
    " over ", length(probe_idx), " probe row(s)."
  )

  if (max_abs_diff > tolerance) {
    return(list(warned = TRUE, message = paste0(message_text, " Potential score mismatch detected.")))
  }

  if (!is.null(verbose) && verbose > 1) {
    message(message_text)
  }

  return(list(warned = FALSE, message = message_text))
}

calc_copula_derivatives = function(eta_inv, Fx_1_2, copula_dist, calc_d2=FALSE, calc_d2_marginal=FALSE, par1=NULL, par2=NULL, pair_complete=NULL) {

  if(is.null(par1)) {
    par1=eta_inv[["theta"]]
  }

  if(is.null(par2)) {
    if("zeta" %in% names(eta_inv)) {
      par2=eta_inv[["zeta"]]
    } else {
      par2=eta_inv[["theta"]]*0
    }
  }

  if(is.null(pair_complete)) {
    pair_complete=rep(TRUE,length(par1))
  }

  if(length(par1)==0) {
    if("zeta" %in% names(eta_inv)) {
      if(calc_d2==TRUE) {
        return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2ldz2=numeric(0),d2ldthdz=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
      }
      return(list(dldth=numeric(0),dcdth=numeric(0),dldz=numeric(0),dcdz=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
    }
    if(calc_d2==TRUE) {
      return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0),d2ldth2=numeric(0),d2cdu12=numeric(0),d2cdu22=numeric(0)))
    }
    return(list(dldth=numeric(0),dcdth=numeric(0),dcdu1=numeric(0),dcdu2=numeric(0)))
  }

  Fx_eval=as.matrix(Fx_1_2)
  Fx_eval[!is.finite(Fx_eval)]=0.5
  Fx_eval[Fx_eval>1]=1
  Fx_eval[Fx_eval<0]=0

  par1_eval=par1
  par2_eval=par2
  par1_eval[!is.finite(par1_eval)]=0
  par2_eval[!is.finite(par2_eval)]=0

  if(copula_dist=="C") {
    par1_eval[par1_eval>=28]=27.9
  }

  copula_d=.copula_pdf(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval)
  copula_d[!is.finite(copula_d) | copula_d<=0]=1
  copula_d[!pair_complete]=1

  dldth=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par",log=TRUE)
  dcdth=copula_d*dldth

  if(calc_d2==TRUE) {
    d2cdth=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par")
    d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  }

  if("zeta" %in% names(eta_inv)) {
    dldz=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par2",log=TRUE)
    dcdz=copula_d*dldz

    if(calc_d2==TRUE) {
      d2cdz=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par2")
      d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

      d2cdthdz=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="par1par2")
      d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
    }

  }
  dcdu1=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u1",log=FALSE)
  dcdu2=.copula_deriv(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u2",log=FALSE)

  dldth[!is.finite(dldth)]=0; if(calc_d2==TRUE) {d2ldth2[!is.finite(d2ldth2)]=0  }
  dldth[!pair_complete]=0
  dcdth[!pair_complete]=0
  dcdu1[!pair_complete]=0
  dcdu2[!pair_complete]=0

  if(calc_d2==TRUE) {
    d2cdu12=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u1")
    d2cdu22=.copula_deriv2(Fx_eval[,1],Fx_eval[,2],family = copula_dist,par=par1_eval,par2=par2_eval,deriv="u2")
    d2cdu12[!is.finite(d2cdu12)]=0
    d2cdu22[!is.finite(d2cdu22)]=0
    d2cdu12[!pair_complete]=0
    d2cdu22[!pair_complete]=0
  }

  if("zeta" %in% names(eta_inv)) {
    dldz[!is.finite(dldz)]=0
    dcdz[!is.finite(dcdz)]=0
    dldz[!pair_complete]=0
    dcdz[!pair_complete]=0
    if(calc_d2==TRUE) {
      d2ldz2[!is.finite(d2ldz2)]=0
      d2ldthdz[!is.finite(d2ldthdz)]=0
      d2ldz2[!pair_complete]=0
      d2ldthdz[!pair_complete]=0
    }
  }

  if(calc_d2==TRUE) {
    d2ldth2[!pair_complete]=0
  }

  ############# RETURN LIST

  if("zeta" %in% names(eta_inv)) {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2,d2ldth2,d2ldz2,d2ldthdz, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2","d2ldth2","d2ldz2","d2ldthdz","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dldz,dcdz,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dldz","dcdz","dcdu1","dcdu2")
    }
  } else {
    if(calc_d2==TRUE) {
      return_list=list(dldth,dcdth,dcdu1,dcdu2,d2ldth2, d2cdu12, d2cdu22)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2","d2ldth2","d2cdu12","d2cdu22")
    } else {
      return_list=list(dldth,dcdth,dcdu1,dcdu2)
      names(return_list)=c("dldth","dcdth","dcdu1","dcdu2")
    }
  }

  return(return_list)
}

calc_Fx_derivatives = function(eta_inv, mm, margin_dist,response,par_names=NULL) {
  # Allow callers to pass full model matrix object; we only need fixed-effect blocks.
  if (is.list(mm) && all(c("x", "s") %in% names(mm))) {
    mm = mm$x
  }

  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]
  if(!is.null(par_names)) {
    margin_par_names=intersect(margin_par_names, par_names)
  }

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change=change_m=change_c=c(0,0)
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=2)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      #change_c[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["copula"]
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    #nd_impact_c[eta_par_names_nd]=(change_c[2]-change_c[1])/(2*adj_fac)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]-change_F[,1])/(2*adj_fac)
  }
  return(nd_impact_F)
}

calc_Fx2_derivatives = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {
  nd_impact=nd_impact_m=nd_impact_c=rep(0,length(names(eta_inv)))
  nd_impact_F=nd_impact_c=list() ###### WE DON"T NEED TO BE CALCULATING THIS FOR ALL PARAMETERS EVERY TIME

  margin_par_names=names(eta_inv)[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {
    adj_fac=.0001
    change_F=matrix(0,nrow=length(eta_inv[[1]]),ncol=3)
    change_c=matrix(0,nrow=length(eta_inv[["theta"]]),ncol=3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      #change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["joint"]
      #change_m[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$log_lik["marginal"]
      if(testing==TRUE) {
        change_c[,i]=log(calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
      }
      change_F[,i]=calc_F_x(eta_inv_adj,mm,margin_dist,response)#calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist)$margin_p
      #print(change_F[,i]==calc_F_x(eta_inv_adj,mm,margin_dist))
      i=i+1
    }
    change_F[,3]=calc_F_x(eta_inv,mm,margin_dist,response)
    if(testing==TRUE) {
      change_c[,3]=log(calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2=FALSE,response,response_margin,response_subject)$copula_d)
    }
    #nd_impact[eta_par_names_nd]=(change[2]-change[1])/(2*adj_fac)
    #nd_impact_m[eta_par_names_nd]=(change_m[2]-change_m[1])/(2*adj_fac)
    nd_impact_c[[eta_par_names_nd]]=(change_c[,2]+change_c[,1]-2*change_c[,3])/(adj_fac^2)
    nd_impact_F[[eta_par_names_nd]]=(change_F[,2]+change_F[,1]-2*change_F[,3])/(adj_fac^2)
  }
  if(testing==FALSE) {
    return(nd_impact_F)
  } else {
    return(list(nd_impact_F,nd_impact_c))
  }
}

calc_true_SE_numderiv_only = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.001
  nd_impact=rep(0,length(names(eta_inv)))
  names(nd_impact)=margin_par_names=names(eta_inv)#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  for (eta_par_names_nd in margin_par_names) {

    change=rep(0,length(names(eta_inv)))
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj
      change[i]=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
      i=i+1
    }
    change[3]=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
    nd_impact[eta_par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    #print(c(change,nd_impact[eta_par_names_nd]))
  }

  nd_cross=matrix(0,nrow=length(names(eta_inv)),ncol=length(names(eta_inv)))
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_rowwise = function(eta_inv, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA) {

  adj_fac=.00001
  nd_impact_m=nd_impact_c=list()

  for (eta_par_names_nd in margin_par_names) {

    change_m=change_c=list()
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {
      eta_inv_adj=eta_inv
      eta_inv_adj[[eta_par_names_nd]]=eta_inv_adj[[eta_par_names_nd]]+adj

      calc_lik_temp=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)

      change_m[[i]]=calc_lik_temp$margin_d
      change_c[[i]]=calc_lik_temp$copula_d
      i=i+1
    }
    calc_lik_temp=calc_likelihood_minimal(eta_inv,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)
    change_m[[3]]=calc_lik_temp$margin_d
    change_c[[3]]=calc_lik_temp$copula_d
    nd_impact_m[[eta_par_names_nd]]=(change_m[[3]]+change_m[[1]]-2*change_m[[3]])/(adj_fac^2)
    nd_impact_c[[eta_par_names_nd]]=(change_c[[3]]+change_c[[1]]-2*change_c[[3]])/(adj_fac^2)
  }
  names(nd_impact_m)=names(nd_impact_c)=margin_par_names=names(eta_inv)

  nd_cross_m
  colnames(nd_cross)=rownames(nd_cross)=names(eta_inv)
  for (name1 in margin_par_names) {
    for (name2 in margin_par_names) {
      for(adj1 in c(-1*adj_fac,adj_fac)) {
        for(adj2 in c(-1*adj_fac,adj_fac)) {
          if(name1!=name2) {
            eta_inv_adj=eta_inv
            eta_inv_adj[[name1]]=eta_inv_adj[[name1]]+adj1
            eta_inv_adj[[name2]]=eta_inv_adj[[name2]]+adj2
            change=calc_likelihood_minimal(eta_inv_adj,mm,margin_dist,copula_dist,calc_d2 = FALSE,response,response_margin,response_subject)$log_lik["joint"]
            nd_cross[name1,name2]=nd_cross[name1,name2]+change*if(adj1==adj2){1} else {-1}
          }
        }
      }
    }
  }
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

calc_true_SE_numderiv_only_covariates = function(object, par, mm, margin_dist,response,testing=FALSE,response_margin=NA,response_subject=NA,h=.0001, progress=interactive()) {

  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  response=object$response
  response_margin=object$response_margin
  response_subject=object$response_subject

  #margin_names=unique(object$response_margin)
  #num_margins=length(margin_names)

  #se_out=object$par*0;
  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  par_cov=object$par
  par_s=object$par_s

  input_par=par_cov
  progress=isTRUE(progress)

  adj_fac=h
  par_names=names(input_par)
  nd_impact=rep(0,length(par_names))
  names(nd_impact)=par_names#[names(eta_inv) %in% c("mu","sigma","nu","tau")]

  # Reuse fixed copula pairing metadata across all numerical derivative evaluations.
  pair_cache=build_copula_pair_cache(response,response_margin,response_subject)
  eval_joint_loglik <- function(par_vec) {
    eta_out=calc_eta(par_cov=par_vec,mm=mm,margin_dist,copula_link,par_s)
    eta_inv=eta_out$eta_inv
    calc_likelihood_minimal(
      eta_inv,
      mm=mm$x,
      margin_dist,
      copula_dist,
      calc_d2=FALSE,
      response=response,
      response_margin=response_margin,
      response_subject=response_subject,
      pair_cache=pair_cache
    )$log_lik["joint"]
  }

  base_loglik=eval_joint_loglik(input_par)

  cat("Calculating numerical first derivates for Hessian matrix...\n")
  pb_first <- NULL
  if(progress) {
    pb_first <- utils::txtProgressBar(min=0, max=length(par_names), style=3)
    on.exit(close(pb_first), add=TRUE)
  }
  first_counter <- 0L
  for (par_names_nd in par_names) {

    #print(par_names_nd)
    change=rep(0,3)
    i=1
    for (adj in c(-1*adj_fac,adj_fac)) {

      par_cov=input_par
      par_cov[[par_names_nd]]=par_cov[[par_names_nd]]+adj

      change[i]=eval_joint_loglik(par_cov)
      i=i+1
    }
    change[3]=base_loglik
    nd_impact[par_names_nd]=(change[2]+change[1]-2*change[3])/(adj_fac^2)

    first_counter <- first_counter + 1L
    if(progress) {
      utils::setTxtProgressBar(pb_first, first_counter)
    }

    #print(c(change,nd_impact[eta_par_names_nd]))
  }
  if(progress) cat("\n")

  cat("Calculating numerical second derivates for Hessian matrix... this may take a while\n")
  p=length(par_names)
  nd_cross=matrix(0,nrow=p,ncol=p)
  colnames(nd_cross)=rownames(nd_cross)=par_names
  second_total <- if(p > 1) choose(p,2) else 0
  pb_second <- NULL
  if(progress && second_total > 0) {
    pb_second <- utils::txtProgressBar(min=0, max=second_total, style=3)
    on.exit(close(pb_second), add=TRUE)
  }
  second_counter <- 0L
  if(p > 1) {
    for (i in 1:(p-1)) {
      name1=par_names[i]
      for (j in (i+1):p) {
        name2=par_names[j]
        cross_sum=0
        for(adj1 in c(-1*adj_fac,adj_fac)) {
          for(adj2 in c(-1*adj_fac,adj_fac)) {
            par=input_par
            par[[name1]]=par[[name1]]+adj1
            par[[name2]]=par[[name2]]+adj2

            change=eval_joint_loglik(par)
            cross_sum=cross_sum+change*if(adj1==adj2){1} else {-1}
          }
        }
        nd_cross[name1,name2]=cross_sum
        nd_cross[name2,name1]=cross_sum

        second_counter <- second_counter + 1L
        if(progress && !is.null(pb_second)) {
          utils::setTxtProgressBar(pb_second, second_counter)
        }
      }
    }
  }
  if(progress && !is.null(pb_second)) cat("\n")
  nd_cross=nd_cross/(4*(adj_fac^2))

  nd2=(diag(nd_impact)+nd_cross)

  return(nd2)
}

#' This function returns the log likelihood for a fitted gamlss.longitudinal object
#' @export
logLik.gamlss.longitudinal=function(object, ...) {
  return(object$calc_lik_out$log_lik)
}

# This function returns the coefficients for a fitted gamlss.longitudinal object
#' @export
coef.gamlss.longitudinal=function(object, ...) {
  return(object$par)
}

# This function returns the variance-covariance matrix for a given gamlss longitudinal object
#' @param method Character; Hessian method to use. \code{"analytical"} (default)
#'   uses the semi-analytical Hessian from \code{R/analytical_hessian.R}.
#'   \code{"numderiv"} uses full finite-difference numerical second
#'   derivatives as a slower reference path. The legacy
#'   \code{numderiv} logical argument is still accepted and maps to
#'   \code{method = "numderiv"} when \code{TRUE}.
#' @export
vcov.gamlss.longitudinal=function(object,par=NA,sep_d2=TRUE,numderiv=FALSE,
                                   method=c("analytical","numderiv","analytical_only"),
                                   progress=interactive(), h=1e-4, ...) {

  #object=fit; par=NA; numderiv=TRUE; sep_d2=TRUE

  method <- match.arg(method)
  # Legacy: numderiv=TRUE overrides method selection
  if (isTRUE(numderiv)) method <- "numderiv"

  progress = isTRUE(progress)

  include_dlcopdpar=TRUE
  response=object$response
  response_margin=object$response_margin
  response_subject=object$response_subject

  margin_names=unique(object$response_margin)
  num_margins=length(margin_names)

  #se_out=object$par*0;
  margin_dist=object$margin_dist; copula_dist=object$copula_dist; copula_link=get_copula_dist(copula_dist)$copula_link
  mm=object$model_matrix

  if(all(is.na(par))) {
    par_cov=object$par
    par_s=object$par_s
  } else {
    par_cov=par$par
    par_s=par$par_s
  }

  eta_out=calc_eta(par_cov,mm,margin_dist,copula_link,par_s=par_s)
  eta_inv=eta_out$eta_inv; eta_dr=eta_out$eta_dr; eta=eta_out$eta; eta_dr=eta_out$eta_dr

  #if(!all(is.na(par))) {response=eta_inv[["mu"]]}
  calc_lik_out=calc_likelihood_minimal(eta_inv,mm=mm$x,margin_dist,copula_dist,calc_d2=TRUE
            ,response=response,response_margin=response_margin,response_subject = response_subject)

  Fx_1_2=calc_lik_out$Fx_1_2; margin_p=calc_lik_out$margin_p; margin_d=calc_lik_out$margin_d; copula_d=calc_lik_out$copula_d

  ###Calculate derivaties: margin and copula d1 and d2
  margin_derivatives=calc_lik_out$margin_deriv
  copula_derivatives=calc_copula_derivatives(eta_inv, Fx_1_2, copula_dist,calc_d2 = TRUE)

  #Calculate numerical derivatives for F(x) - also used as a reference for overall likelihood derivative
  nd_impact_F=calc_Fx_derivatives(eta_inv,mm$x,margin_dist,response)
  nd_impact_F2=calc_Fx2_derivatives(eta_inv,mm$x,margin_dist,response)

  if (method == "numderiv") {
    #nd2_joint_lik=calc_true_SE_numderiv_only(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)
    hessian_nd=calc_true_SE_numderiv_only_covariates(object=object,par=par_cov,mm=mm$x,margin_dist=margin_dist,response=response,testing=FALSE,response_margin=response_margin,response_subject=response_subject,progress=progress)
  } else if (method %in% c("analytical", "analytical_only")) {
    # Source the analytical Hessian helpers if not already loaded.
    if (!exists("calc_analytical_hessian", mode = "function")) {
      # Try relative to this file, then working directory
      candidate_paths <- c(
        file.path(dirname(attr(body(vcov.gamlss.longitudinal), "srcfile")$filename %||% ""), "analytical_hessian.R"),
        "R/analytical_hessian.R",
        file.path(getwd(), "R", "analytical_hessian.R")
      )
      loaded <- FALSE
      for (cp in candidate_paths) {
        if (file.exists(cp)) { source(cp, local = FALSE); loaded <- TRUE; break }
      }
      if (!loaded) stop("Cannot locate analytical_hessian.R. Source it manually or ensure the working directory is the package root.")
    }
    hessian_nd <- calc_analytical_hessian(object, progress = progress, h = h)
  } else {

    #######to delete############
    nd_impact_C2=calc_Fx2_derivatives(eta_inv,mm,margin_dist,response,testing=TRUE,response_margin,response_subject)[[2]]

    ### MARGIN LIKELIHOOD DERIVATIVES
    margin_deriv_subnames=c("m","d","v","t")
    names(margin_deriv_subnames)=c("mu","sigma","nu","tau")
    margin_par=names(mm)[names(mm) %in% c("mu","sigma","nu","tau")]
    order_margin=cbind(object$response_margin,object$response_subject)
    colnames(order_margin)=c("time","subject")

    dcdu1=copula_derivatives$dcdu1; dcdu2=copula_derivatives$dcdu2;d2cdu12=copula_derivatives$d2cdu12;d2cdu22=copula_derivatives$d2cdu22

    #####################For each parameter...
    d1_cop=d2_cop=matrix(0,nrow=length(response),ncol=length(margin_par))
    colnames(d1_cop)=colnames(d2_cop)=margin_par
    pair_cache_diag <- build_copula_pair_cache(response, response_margin, response_subject)
    for (par_name in margin_par) {
      if(object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) {

        order_copula=data.frame()
        for (i in 1:(num_margins-1)) {
          order_copula=rbind(order_copula,cbind(order_margin[response_margin == margin_names[i],c("time","subject")],order_margin[response_margin == margin_names[i+1],c("time","subject")]))
        }
        colnames(order_copula)=c("time1","subject1","time2","subject2")

        margin_deriv_1=matrix(0,ncol=length(margin_par),nrow=length(response))
        colnames(margin_deriv_1)=paste("dld",margin_par,sep="")
        margin_deriv_1[,paste("dld",par_name,sep="")]=margin_derivatives[grepl("dld",names(margin_derivatives))][[which(margin_par==par_name)]]

        #COPULA DERIVS WITH RESPECT TO

        mu=eta_inv[["mu"]]
        F_nd=nd_impact_F[[par_name]]
        F_nd2=nd_impact_F2[[par_name]]
        c_nd2=nd_impact_C2[[par_name]]

        margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,mu,F_nd,F_nd2)
        margin_components_Ft_plus=margin_components
        margin_components_Ft_plus[,"time"]=normalize_lag_time(margin_components_Ft_plus[,"time"])
        margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

        copula_components=cbind(
          order_copula,
          row_id1=pair_cache_diag$row_id1,
          row_id2=pair_cache_diag$row_id2,
          dcdu1,
          dcdu2,
          copula_d,
          d2cdu12,
          d2cdu22,
          c_nd2
        )
        copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)

        #Calculate copula derivative with respect to marginal parameters
        input=copula_merged
        d1_cop[,par_name]=calc_deriv_copula_wrt_margin(input,margin_par,par_name,calc_d2=FALSE)[,which(margin_par==par_name)]

        #OK so let's calcute the numerical d2lcopdpar and pass it through input
        d2_cop[,par_name]=calc_deriv_copula_wrt_margin(input,margin_par,par_name,calc_d2=TRUE)[,which(margin_par==par_name)]
      }
    }

    ###########Need d1 and d2 for score function

    m_d1_names=names(margin_derivatives)[grepl("dld",names(margin_derivatives))]
    c_d1_names=names(copula_derivatives)[grepl("dld",names(copula_derivatives))]

    m_d2_names=names(margin_derivatives)[grepl("d2ld",names(margin_derivatives))]
    c_d2_names=names(copula_derivatives)[grepl("d2ld",names(copula_derivatives))]

    d1_all=list(); d2_all=list()

    i=1
    for(par_name in c(m_d1_names)) {
      d1_all[[par_name]]=c(margin_derivatives[[par_name]])
      if((object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE)) {
        d1_all[[par_name]]=c(margin_derivatives[[par_name]]+d1_cop[,i])
        i=i+1
      }
    }

    for(par_name in c(c_d1_names)) {
      d1_all[[par_name]]=c(copula_derivatives[[par_name]])
    }

    names(d1_all)=names(mm)

    i=1
    for(par_name in c(m_d2_names)) {
      d2_all[[par_name]]=c(margin_derivatives[[par_name]])
      if((object$include_dlcopdpar==TRUE | include_dlcopdpar==TRUE) & endsWith(par_name,"2")) {
        d2_all[[par_name]]=c(margin_derivatives[[par_name]]+d2_cop[,i]*(if(sep_d2==TRUE) {0} else {1}))
        i=i+1
      }
    }
    for(par_name in c(c_d2_names)) {
      d2_all[[par_name]]=c(copula_derivatives[[par_name]])
    }

    d2_all_mean=rep(0,length=length(d2_all))
    names(d2_all_mean)=names(d2_all)
    for (deriv_name in names(d2_all)) {
      d2_all_mean[deriv_name]=mean(d2_all[[deriv_name]])
    }

    d2_mat_diag=d2_all_mean[endsWith(names(d2_all_mean),"2")]
    d2_mat_cross=d2_all_mean[!endsWith(names(d2_all_mean),"2")]
    d2_mat=matrix(nrow=length(eta),ncol=length(eta))
    #print(d2_mat);print(names(eta))
    colnames(d2_mat)=rownames(d2_mat)=names(eta)

    copula_deriv_subnames=c("th","z")
    names(copula_deriv_subnames)=c("theta","zeta")
    all_names=c(margin_deriv_subnames,copula_deriv_subnames)
    sub_names_in=all_names[names(eta)]
    print(sub_names_in)

    for (row_name in rownames(d2_mat)) {
      for (col_name in colnames(d2_mat)) {
        if(!row_name==col_name) {
          deriv_name_temp=paste("d2ld",sub_names_in[row_name],"d",sub_names_in[col_name],sep="")

          if(is.na(d2_all_mean[deriv_name_temp])) {
            deriv_name_temp=paste("d2ld",sub_names_in[col_name],"d",sub_names_in[row_name],sep="")
          }

          deriv_val_temp=d2_all_mean[deriv_name_temp]
          d2_mat[row_name,col_name]=deriv_val_temp
        }
      }
    }

    #cop_row=(grepl("theta",rownames(d2_mat))|grepl("zeta",rownames(d2_mat)))
    #d2_mat[!cop_row,!cop_row][upper.tri(d2_mat[!cop_row,!cop_row])]=d2_mat_cross
    diag(d2_mat)=d2_mat_diag

    d2_mat[is.na(d2_mat)]=0

  }

  if (method %in% c("numderiv", "analytical", "analytical_only")) {
    vcov_final=-solve(hessian_nd)
    se_final=sqrt(abs(diag(solve(hessian_nd))))
  } else {
    vcov_final = -(solve((d2_mat)))/(length(response))
    se_final=sqrt(abs(diag(vcov_final)))
  }

  ###########TESTING APPROACH FOR ESTIMATING SMOOTHER VARIANCE

  # Method 1: Bayesian/Mixed Model Approach for Smoother Variance
  # Calculate variance-covariance matrix for smooth terms using: Var(beta) = (X'WX + lambda*P)^(-1) * sigma^2

  smooth_vcov_list = list()
  smooth_se_list = list()

  # Extract residual variance estimate (using reciprocal of mean weights as proxy for sigma^2)
  if(!is.null(object$weights) && length(object$weights) > 0 && is.numeric(object$weights)) {
    sigma2_est = 1 / mean(object$weights, na.rm = TRUE)
  } else {
    # Fallback: estimate from residuals if weights not available
    fitted_response = eta_inv[["mu"]]
    residuals = response - fitted_response
    sigma2_est = var(residuals, na.rm = TRUE)
  }

  # Ensure scalar numeric to avoid deprecated array recycling warnings.
  sigma2_est = as.numeric(sigma2_est)[1]
  if(!is.finite(sigma2_est)) {
    sigma2_est = 1
  }

  # Process each parameter that has smooth terms
  for(par_name in names(object$par_s)) {
    if(length(object$par_s[[par_name]]) > 0) {

      smooth_vcov_list[[par_name]] = list()
      smooth_se_list[[par_name]] = list()

      # Process each smooth term for this parameter
      for(s_name in names(object$par_s[[par_name]])) {

        # Get the B-spline basis matrix
        B = object$model_matrix$s[[par_name]][[s_name]]

        # Get the smoothing parameter
        lambda = object$lambda_s[[par_name]][[s_name]]

        # Use the mgcv-generated penalty stored on the basis matrix; fall back to
        # a generic second-difference penalty only when unavailable.
        k = ncol(B)
        pen_attr = attr(B, "penalty")
        if (!is.null(pen_attr) && is.matrix(pen_attr) &&
            nrow(pen_attr) == k && ncol(pen_attr) == k) {
          P = pen_attr
        } else if (k > 2) {
          D2 = diff(diag(k), differences = 2)
          P = t(D2) %*% D2
        } else {
          P = diag(k)
        }

        # Get per-parameter IRLS working weights. object$weights is a named list
        # keyed by parameter name; fall back to unit weights if not available.
        w_par = object$weights[[par_name]]
        if (!is.null(w_par) && is.numeric(w_par) && length(w_par) == nrow(B)) {
          w_diag = as.vector(w_par)
        } else {
          w_diag = rep(1, nrow(B))
        }
        W = diag(w_diag)

        # Per-parameter sigma2: scale consistent with IRLS, 1/mean(w)
        sigma2_par = if (all(w_diag > 0)) 1 / mean(w_diag) else sigma2_est

        # Calculate the penalized precision matrix: X'WX + lambda*P
        XWX = t(B) %*% W %*% B
        penalized_precision = XWX + lambda * P

        # Variance-covariance matrix for this smooth: (X'WX + lambda*P)^(-1) * sigma^2
        tryCatch({
          smooth_vcov = solve(penalized_precision) * sigma2_par
          smooth_se = sqrt((diag(smooth_vcov)))

          # Store results
          smooth_vcov_list[[par_name]][[s_name]] = smooth_vcov
          smooth_se_list[[par_name]][[s_name]] = smooth_se

          # Also calculate the smoother matrix for fitted values variance
          # A = X(X'WX + lambda*P)^(-1)X'W
          smoother_matrix = B %*% solve(penalized_precision) %*% t(B) %*% W
          fitted_se = sqrt(abs(as.vector(diag(smoother_matrix))) * sigma2_par)

          cat(sprintf("\nSmooth term variance estimates for %s:%s\n", par_name, s_name))
          cat(sprintf("  Basis coefficients SE: min=%.4f, max=%.4f, mean=%.4f\n",
                     min(smooth_se), max(smooth_se), mean(smooth_se)))
          cat(sprintf("  Fitted values SE: min=%.4f, max=%.4f, mean=%.4f\n",
                     min(fitted_se), max(fitted_se), mean(fitted_se)))
          cat(sprintf("  Effective DF: %.2f (trace of smoother matrix)\n",
                     sum(diag(smoother_matrix))))
          cat(sprintf("  Smoothing parameter lambda: %.4f\n", lambda))

        }, error = function(e) {
          warning(sprintf("Could not calculate variance for smooth %s:%s - %s",
                         par_name, s_name, e$message))
          smooth_vcov_list[[par_name]][[s_name]] = NULL
          smooth_se_list[[par_name]][[s_name]] = NULL
        })
      }
    }
  }

  # Add smooth variance results to return list
  vcov_final_with_smooth = list(
    overall = vcov_final,
    smooth_vcov = smooth_vcov_list,
    smooth_se = smooth_se_list
  )

  se_final_with_smooth = list(
    overall = se_final,
    smooth_se = smooth_se_list
  )

  return(list(vcov=vcov_final_with_smooth, se=se_final_with_smooth))

}

.can_use_cached_vcov <- function(object, numderiv = FALSE, method = NULL, extra_args = list()) {
  if(!inherits(object, "gamlss.longitudinal")) return(FALSE)
  extra_args_cache <- extra_args
  extra_args_cache$method <- NULL
  if(!is.null(extra_args_cache) && length(extra_args_cache) > 0) return(FALSE)
  if(is.null(object$vcov) || !is.list(object$vcov)) return(FALSE)
  if(is.null(object$vcov$vcov) || is.null(object$vcov$vcov$overall)) return(FALSE)

  if(!is.null(object$vcov_meta) && !is.null(object$vcov_meta$numderiv)) {
    numderiv_ok <- identical(isTRUE(object$vcov_meta$numderiv), isTRUE(numderiv))
    method_ok <- is.null(method) ||
      is.null(object$vcov_meta$method) ||
      identical(as.character(object$vcov_meta$method)[1], as.character(method)[1])
    return(numderiv_ok && method_ok)
  }

  TRUE
}

.resolve_vcov <- function(object, numderiv = FALSE, extra_args = list()) {
  vcov_method <- extra_args$method %||% if (isTRUE(numderiv)) "numderiv" else "analytical"
  if(.can_use_cached_vcov(object, numderiv = numderiv, method = vcov_method, extra_args = extra_args)) {
    return(object$vcov)
  }

  if(is.null(extra_args$method)) {
    extra_args$method <- vcov_method
  }

  do.call(
    vcov.gamlss.longitudinal,
    c(list(object = object, numderiv = numderiv), extra_args)
  )
}

#' Summarize a fitted gamlss.longitudinal model
#'
#' Creates a compact summary of key model diagnostics and coefficient estimates
#' for a fitted `gamlss.longitudinal` object. By default it wraps
#' `vcov.gamlss.longitudinal()` to provide standard errors and confidence
#' intervals for fixed effects.
#'
#' @param object A fitted object of class `gamlss.longitudinal`.
#' @param include_vcov Logical; if `TRUE`, compute and include variance-covariance
#'   output via `vcov.gamlss.longitudinal()`.
#' @param numderiv Logical passed to `vcov.gamlss.longitudinal()`.
#' @param ci_level Confidence level for coefficient intervals.
#' @param ... Additional arguments passed to `vcov.gamlss.longitudinal()`.
#'
#' @return An object of class `summary.gamlss.longitudinal` containing:
#' - model dimensions and parameter counts,
#' - likelihood and information criteria,
#' - fixed-effect coefficient table,
#' - optional `vcov` output.
#' @export
summary.gamlss.longitudinal = function(
  object,
  include_vcov = TRUE,
  numderiv = FALSE,
  ci_level = 0.95,
  ...
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  n_obs = length(object$response)
  n_subjects = length(unique(object$response_subject))
  n_timepoints = length(unique(object$response_margin))

  n_fixed = length(object$par)
  n_smooth_terms = 0
  if(!is.null(object$par_s)) {
    n_smooth_terms = sum(vapply(object$par_s, length, integer(1)))
  }

  edf_smooth = NA_real_
  if(!is.null(object$df_s) && length(object$df_s) > 0) {
    df_vals = suppressWarnings(as.numeric(unlist(object$df_s, use.names = FALSE)))
    df_vals = df_vals[is.finite(df_vals)]
    if(length(df_vals) > 0) {
      edf_smooth = sum(df_vals)
    }
  }

  loglik_vec = c(marginal = NA_real_, copula = NA_real_, joint = NA_real_)
  if(!is.null(object$calc_lik_out_end) && !is.null(object$calc_lik_out_end$log_lik)) {
    ll_in = object$calc_lik_out_end$log_lik
    for(nm in c("marginal", "copula", "joint")) {
      if(nm %in% names(ll_in)) {
        loglik_vec[nm] = as.numeric(ll_in[[nm]])
      }
    }
  }
  loglik_joint = as.numeric(loglik_vec["joint"])

  p_cop = object$par[grepl("theta", names(object$par)) | grepl("zeta", names(object$par))]
  p_mar = object$par[!(grepl("theta", names(object$par)) | grepl("zeta", names(object$par)))]

  df_s_total = 0
  df_s_cop_total = 0
  df_s_margin_total = 0
  if(!is.null(object$df_s) && length(object$df_s) > 0) {
    for(par_name in names(object$df_s)) {
      df_val = suppressWarnings(sum(as.numeric(unlist(object$df_s[[par_name]])), na.rm = TRUE))
      if(!is.finite(df_val)) df_val = 0

      if(par_name %in% c("theta", "zeta")) {
        df_s_cop_total = df_s_cop_total + df_val
      } else {
        df_s_margin_total = df_s_margin_total + df_val
      }
      df_s_total = df_s_total + df_val
    }
  }

  edf_vec = c(
    marginal = length(p_mar) + df_s_margin_total,
    copula = length(p_cop) + df_s_cop_total,
    joint = length(object$par) + df_s_total
  )

  aic_vec = -2 * loglik_vec + 2 * edf_vec
  bic_vec = -2 * loglik_vec + log(max(1, n_obs)) * edf_vec
  model_selection = rbind(LogLik = loglik_vec, AIC = aic_vec, BIC = bic_vec, EDF = edf_vec)

  coef_tbl = data.frame(
    term = names(object$par),
    estimate = as.numeric(object$par),
    std_error = NA_real_,
    p_value = NA_real_,
    signif = NA_character_,
    stringsAsFactors = FALSE
  )

  coef_tbl$.original_order = seq_len(nrow(coef_tbl))
  coef_tbl$parameter = sub("\\..*$", "", coef_tbl$term)

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  coef_tbl$.param_rank = match(coef_tbl$parameter, param_order)
  coef_tbl$.param_rank[is.na(coef_tbl$.param_rank)] = length(param_order) + 1L

  vcov_out = NULL
  if(isTRUE(include_vcov)) {
    vcov_out = .resolve_vcov(
      object = object,
      numderiv = numderiv,
      extra_args = list(...)
    )

    if(!is.null(vcov_out$vcov) && !is.null(vcov_out$vcov$overall)) {
      V = vcov_out$vcov$overall
      se = NULL
      if(!is.null(vcov_out$se) && !is.null(vcov_out$se$overall)) {
        se = as.numeric(vcov_out$se$overall)
        se_names = names(vcov_out$se$overall)
      } else {
        se = sqrt(pmax(0, diag(V)))
        se_names = names(diag(V))
      }

      if(is.null(se_names) && !is.null(rownames(V)) && length(rownames(V)) == length(se)) {
        se_names = rownames(V)
      }

      if(!is.null(se_names)) {
        names(se) = se_names
      }

      if(!is.null(names(se))) {
        idx = match(coef_tbl$term, names(se))
        coef_tbl$std_error = se[idx]
      } else if(length(se) == nrow(coef_tbl)) {
        coef_tbl$std_error = se
      }

      z_abs = abs(coef_tbl$estimate / coef_tbl$std_error)
      coef_tbl$p_value = 2 * stats::pnorm(z_abs, lower.tail = FALSE)
      coef_tbl$signif = ifelse(
        is.na(coef_tbl$p_value),
        NA_character_,
        ifelse(coef_tbl$p_value < 0.001, "***",
               ifelse(coef_tbl$p_value < 0.01, "**",
                      ifelse(coef_tbl$p_value < 0.05, "*",
                             ifelse(coef_tbl$p_value < 0.1, ".", " "))))
      )
    }
  }

  coef_tbl = coef_tbl[order(coef_tbl$.param_rank, coef_tbl$.original_order), , drop = FALSE]
  rownames(coef_tbl) = NULL

  out = list(
    model = list(
      margin_dist = if(!is.null(object$margin_dist$family[1])) as.character(object$margin_dist$family[1]) else NA_character_,
      copula_dist = object$copula_dist,
      n_obs = n_obs,
      n_subjects = n_subjects,
      n_timepoints = n_timepoints,
      n_fixed = n_fixed,
      n_smooth_terms = n_smooth_terms,
      edf_smooth = edf_smooth
    ),
    fit = list(
      logLik = loglik_joint,
      AIC = as.numeric(aic_vec["joint"]),
      BIC = as.numeric(bic_vec["joint"]),
      ci_level = ci_level,
      vcov_included = isTRUE(include_vcov),
      vcov_numderiv = isTRUE(numderiv),
      model_selection = model_selection
    ),
    smooth_terms = {
      st = list()
      if(!is.null(object$par_s) && length(object$par_s) > 0) {
        for(par_name in names(object$par_s)) {
          if(length(object$par_s[[par_name]]) == 0) next
          for(s_name in names(object$par_s[[par_name]])) {
            st[[length(st) + 1]] = data.frame(
              parameter = par_name,
              smooth_term = s_name,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      if(length(st) == 0) {
        data.frame(parameter = character(0), smooth_term = character(0), stringsAsFactors = FALSE)
      } else {
        do.call(rbind, st)
      }
    },
    coefficients = within(coef_tbl, {
      .original_order = NULL
      .param_rank = NULL
    }),
    vcov = vcov_out
  )
  class(out) = "summary.gamlss.longitudinal"
  out
}

.copula_v2_gaussian_limit_warning <- function(object, threshold = 7.5) {
  if (is.null(object) || is.null(object$model$copula_dist)) {
    return(NULL)
  }

  copula_dist <- object$model$copula_dist
  if (!identical(copula_dist, "t")) {
    return(NULL)
  }

  coef_tbl <- object$coefficients
  if (is.null(coef_tbl) || nrow(coef_tbl) == 0 || !("parameter" %in% names(coef_tbl))) {
    return(NULL)
  }

  zeta_rows <- coef_tbl[coef_tbl$parameter == "zeta" & is.finite(coef_tbl$estimate), , drop = FALSE]
  if (nrow(zeta_rows) == 0) {
    return(NULL)
  }

  zeta_link <- zeta_rows$estimate
  near_limit <- is.finite(zeta_link) & zeta_link >= threshold
  if (!any(near_limit)) {
    return(NULL)
  }

  zeta_nat <- exp(zeta_link[near_limit]) + 2
  zeta_label <- paste0(formatC(zeta_link[near_limit], format = "f", digits = 2), collapse = ", ")
  df_label <- paste0(formatC(zeta_nat, format = "f", digits = 1), collapse = ", ")

  paste0(
    "WARNING: t-copula zeta is near the Gaussian-limit regime (link-scale zeta = ",
    zeta_label,
    "; implied degrees of freedom ",
    df_label,
    "). The t-copula may be collapsing toward Gaussian dependence."
  )
}

#' @export
print.summary.gamlss.longitudinal = function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nGAMLSS Longitudinal Model Summary\n")
  cat("--------------------------------\n")
  cat("Margin distribution:", x$model$margin_dist, "\n")
  cat("Copula distribution:", x$model$copula_dist, "\n")

  copula_warning <- .copula_v2_gaussian_limit_warning(x)
  if (!is.null(copula_warning)) {
    cat(copula_warning, "\n")
  }

  cat("Observations:", x$model$n_obs,
      " | Subjects:", x$model$n_subjects,
      " | Time points:", x$model$n_timepoints, "\n")
  cat("Fixed coefficients:", x$model$n_fixed,
      " | Smooth terms:", x$model$n_smooth_terms,
      " | Smooth EDF:", format(round(x$model$edf_smooth, digits), nsmall = 2), "\n")

  cat("\nFixed coefficients:\n")
  cat("--------------------\n")
  coef_tbl = x$coefficients
  p_value_raw = coef_tbl$p_value
  coef_tbl$estimate = round(coef_tbl$estimate, digits)
  coef_tbl$std_error = round(coef_tbl$std_error, digits)
  coef_tbl$p_value = round(coef_tbl$p_value, digits + 1)

  fmt_num = function(v, d) ifelse(is.na(v), "NA", formatC(v, format = "f", digits = d))
  fmt_p_value = function(v, v_raw, d) {
    ifelse(
      is.na(v_raw),
      "NA",
      ifelse(v_raw > 0 & v_raw < 10^(-d), paste0("<", formatC(10^(-d), format = "f", digits = d)), fmt_num(v, d))
    )
  }
  coef_disp = data.frame(
    term = as.character(coef_tbl$term),
    estimate = fmt_num(coef_tbl$estimate, digits),
    std_error = fmt_num(coef_tbl$std_error, digits),
    p_value = fmt_p_value(coef_tbl$p_value, p_value_raw, digits + 1),
    signif = ifelse(is.na(coef_tbl$signif), "", as.character(coef_tbl$signif)),
    parameter = as.character(coef_tbl$parameter),
    stringsAsFactors = FALSE
  )

  w_term = max(nchar("term"), nchar(coef_disp$term, type = "width"), na.rm = TRUE)
  w_est = max(nchar("estimate"), nchar(coef_disp$estimate, type = "width"), na.rm = TRUE)
  w_se = max(nchar("std_error"), nchar(coef_disp$std_error, type = "width"), na.rm = TRUE)
  w_p = max(nchar("p_value"), nchar(coef_disp$p_value, type = "width"), na.rm = TRUE)
  w_sig = max(nchar("signif"), nchar(coef_disp$signif, type = "width"), na.rm = TRUE)

  format_row = function(term, estimate, std_error, p_value, signif) {
    sprintf(
      "%-*s  %*s  %*s  %*s  %-*s",
      w_term, term,
      w_est, estimate,
      w_se, std_error,
      w_p, p_value,
      w_sig, signif
    )
  }

  print_coef_block = function(block, prefix = "    ") {
    hdr = format_row("term", "estimate", "std_error", "p_value", "signif")
    cat(prefix, hdr, "\n", sep = "")
    for(ii in seq_len(nrow(block))) {
      row_txt = format_row(
        block$term[ii],
        block$estimate[ii],
        block$std_error[ii],
        block$p_value[ii],
        block$signif[ii]
      )
      cat(prefix, row_txt, "\n", sep = "")
    }
  }

  param_order = c("mu", "sigma", "nu", "tau", "theta", "zeta")
  params_present = unique(coef_tbl$parameter)
  params_print = c(param_order[param_order %in% params_present], setdiff(params_present, param_order))

  for(k in seq_along(params_print)) {
    p = params_print[k]
    block = coef_disp[coef_disp$parameter == p, c("term", "estimate", "std_error", "p_value", "signif"), drop = FALSE]
    cat(sprintf("  [%s]\n", p))
    print_coef_block(block, prefix = "    ")
    if(k < length(params_print)) cat("  --------------------\n")
  }

  cat("  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")


  cat("\nSmooth terms:\n")
  cat("--------------------\n")
  if(!is.null(x$smooth_terms) && nrow(x$smooth_terms) > 0) {
    print(x$smooth_terms, row.names = FALSE)
    cat("Use plot(object) to visualize smooth and fixed terms with confidence bands.\n")
  } else {
    cat("None\n")
  }

  cat("\nModel Selection Criteria:\n")
  cat("--------------------\n")
  if(!is.null(x$fit$model_selection)) {
    print(round(x$fit$model_selection, digits))
  } else {
    fit_tbl = data.frame(
      metric = c("logLik", "AIC", "BIC"),
      value = c(x$fit$logLik, x$fit$AIC, x$fit$BIC),
      stringsAsFactors = FALSE
    )
    print(fit_tbl, row.names = FALSE, digits = digits)
  }
  cat("--------------------------------\n")

  invisible(x)
}

#' Plot all smooth terms with confidence bands
#'
#' This utility plots every smooth term in a fitted `gamlss.longitudinal` object
#' and computes pointwise confidence bands using the smooth coefficient
#' covariance matrices returned by `vcov.gamlss.longitudinal()`.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with the analytical vcov path.
#' @param data Optional data frame containing original covariates used for the
#' x-axis variable of each smooth.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted smooth line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted smooth.
#' @param sort_x Logical; sort points by x before plotting lines.
#' @param even_grid Logical; if TRUE, plot smooths on an evenly spaced x-grid
#' built over observed x-range.
#' @param grid_n Number of grid points when `even_grid = TRUE`.
#' @param fallback_to_index Logical; if x variable cannot be inferred, plot
#' against row index.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each smooth term.
#' @export
plot_smooth_terms = function(
  object,
  vcov_obj = NULL,
  data = NULL,
  ci_level = 0.95,
  ncol = NULL,
  ci_col = "red",
  fit_col = "black",
  ci_lty = 2,
  fit_lwd = 2,
  sort_x = TRUE,
  even_grid = TRUE,
  grid_n = 200,
  fallback_to_index = TRUE,
  setup_mfrow = TRUE,
  show_legend = TRUE
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if(is.null(vcov_obj)) {
    vcov_obj = .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))
  }

  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal().")
  }

  smooth_vcov_list = vcov_obj$vcov$smooth_vcov
  smooth_se_list = vcov_obj$vcov$smooth_se

  extract_smooth_var = function(s_name) {
    s_txt = trimws(s_name)
    s_call = tryCatch(parse(text = s_txt)[[1]], error = function(e) NULL)
    if(!is.null(s_call) && length(s_call) >= 2) {
      out = paste(deparse(s_call[[2]]), collapse = " ")
    } else {
      out = sub("^s\\((.*)\\)$", "\\1", s_txt)
    }
    out = trimws(gsub("`", "", out, fixed = TRUE))
    out
  }

  eval_smooth_x = function(x_expr, data_frame) {
    if(is.null(data_frame) || !is.data.frame(data_frame)) {
      return(NULL)
    }
    tryCatch(eval(parse(text = x_expr)[[1]], envir = data_frame), error = function(e) NULL)
  }

  get_x_for_smooth = function(par_name, s_name, B) {
    x_var = extract_smooth_var(s_name)
    x = NULL

    # Prefer x saved with the smooth basis because it is guaranteed row-aligned.
    x_basis = attr(B, "smooth_x")
    x_basis_var = attr(B, "smooth_var")
    if(!is.null(x_basis) && length(x_basis) == nrow(B)) {
      x = x_basis
      if(!is.null(x_basis_var) && nzchar(x_basis_var)) {
        x_var = x_basis_var
      }
    }

    if(is.null(x) && !is.null(data) && is.data.frame(data)) {
      data_names = names(data)
      # Prefer exact match, then case-insensitive match, then make.names match.
      idx_exact = which(data_names == x_var)
      idx_ci = which(tolower(data_names) == tolower(x_var))
      idx_mn = which(make.names(data_names) == make.names(x_var))
      idx = c(idx_exact, idx_ci, idx_mn)
      idx = idx[!duplicated(idx)]
      if(length(idx) > 0) {
        matched_name = data_names[idx[1]]
        x_candidate = data[[matched_name]]
        if(length(x_candidate) == nrow(B)) {
          x = x_candidate
        } else if(!is.null(rownames(B)) && !is.null(rownames(data))) {
          row_idx = match(rownames(B), rownames(data))
          if(all(!is.na(row_idx))) {
            x = x_candidate[row_idx]
          }
        }
        x_var = matched_name
      }

      if(is.null(x)) {
        x_candidate = eval_smooth_x(x_var, data)
        if(!is.null(x_candidate) && length(x_candidate) == nrow(B)) {
          x = x_candidate
        }
      }
    }

    if(is.null(x) && !is.null(object$model_matrix$x[[par_name]]) && x_var %in% colnames(object$model_matrix$x[[par_name]])) {
      x = object$model_matrix$x[[par_name]][, x_var]
    }

    if(is.null(x) && fallback_to_index) {
      x = seq_len(nrow(B))
      x_var = "index"
      warning("Falling back to index for smooth term '", s_name,
              "' because covariate was not found in supplied data/model matrix.")
    }

    if(is.null(x)) {
      stop("Could not infer x-axis for smooth term '", s_name, "'. Provide 'data' with the smooth covariate columns.")
    }

    if(length(x) != nrow(B)) {
      stop("Length mismatch for smooth term '", s_name, "': length(x)=", length(x), " but nrow(B)=", nrow(B), ".")
    }

    list(x = as.numeric(x), x_var = x_var)
  }

  z = qnorm((1 + ci_level) / 2)
  smooth_index = list()

  for(par_name in names(object$par_s)) {
    if(length(object$par_s[[par_name]]) == 0) next
    for(s_name in names(object$par_s[[par_name]])) {
      B = object$model_matrix$s[[par_name]][[s_name]]
      beta_s = object$par_s[[par_name]][[s_name]]
      if(is.null(B) || is.null(beta_s)) next

      smooth_index[[length(smooth_index) + 1]] = list(par_name = par_name, s_name = s_name)
    }
  }

  n_plots = length(smooth_index)
  if(n_plots == 0) {
    warning("No smooth terms found to plot.")
    return(invisible(list()))
  }

  out = list()
  plot_objects = list()
  for(i in seq_len(n_plots)) {
    par_name = smooth_index[[i]]$par_name
    s_name = smooth_index[[i]]$s_name

    B = object$model_matrix$s[[par_name]][[s_name]]
    beta_s = object$par_s[[par_name]][[s_name]]
    x_info = get_x_for_smooth(par_name, s_name, B)
    x = x_info$x

    fitted_smooth = as.numeric(B %*% beta_s)

    smooth_vcov = NULL
    smooth_se = NULL
    if(!is.null(smooth_vcov_list) && !is.null(smooth_vcov_list[[par_name]])) {
      smooth_vcov = smooth_vcov_list[[par_name]][[s_name]]
    }
    if(!is.null(smooth_se_list) && !is.null(smooth_se_list[[par_name]])) {
      smooth_se = smooth_se_list[[par_name]][[s_name]]
    }

    if(!is.null(smooth_vcov) && all(dim(smooth_vcov) == c(ncol(B), ncol(B)))) {
      smooth_fit_se = sqrt(pmax(0, diag(B %*% smooth_vcov %*% t(B))))
    } else if(!is.null(smooth_se) && length(smooth_se) == ncol(B)) {
      beta_var_diag = as.numeric(smooth_se)^2
      smooth_fit_se = sqrt(pmax(0, rowSums((B^2) * rep(beta_var_diag, each = nrow(B)))))
    } else {
      smooth_fit_se = rep(NA_real_, nrow(B))
    }

    ci_lower = fitted_smooth - z * smooth_fit_se
    ci_upper = fitted_smooth + z * smooth_fit_se
    main_title = paste(par_name, s_name, sep = ": ")
    ylab_text = paste("smooth(", x_info$x_var, ")", sep = "")

    if(isTRUE(even_grid)) {
      x_ok = is.finite(x)
      df_obs = data.frame(
        x = x[x_ok],
        fitted = fitted_smooth[x_ok],
        ci_lower = ci_lower[x_ok],
        ci_upper = ci_upper[x_ok]
      )

      if(nrow(df_obs) >= 2 && length(unique(df_obs$x)) >= 2) {
        agg_df = stats::aggregate(df_obs[, c("fitted", "ci_lower", "ci_upper")], by = list(x = df_obs$x), FUN = mean)
        agg_df = agg_df[order(agg_df$x), , drop = FALSE]
        n_grid_use = max(20, as.integer(grid_n))
        x_grid = seq(min(agg_df$x), max(agg_df$x), length.out = n_grid_use)

        plot_df = data.frame(
          x = x_grid,
          fitted = stats::approx(agg_df$x, agg_df$fitted, xout = x_grid, method = "linear", rule = 2)$y,
          ci_lower = stats::approx(agg_df$x, agg_df$ci_lower, xout = x_grid, method = "linear", rule = 2)$y,
          ci_upper = stats::approx(agg_df$x, agg_df$ci_upper, xout = x_grid, method = "linear", rule = 2)$y
        )
      } else {
        ord = if(sort_x) order(x) else seq_along(x)
        plot_df = data.frame(
          x = x[ord],
          fitted = fitted_smooth[ord],
          ci_lower = ci_lower[ord],
          ci_upper = ci_upper[ord]
        )
      }
    } else {
      ord = if(sort_x) order(x) else seq_along(x)
      plot_df = data.frame(
        x = x[ord],
        fitted = fitted_smooth[ord],
        ci_lower = ci_lower[ord],
        ci_upper = ci_upper[ord]
      )
    }

    y_vals = c(plot_df$fitted, plot_df$ci_lower, plot_df$ci_upper)
    y_vals = y_vals[is.finite(y_vals)]
    y_lim = NULL
    if(length(y_vals) > 0) {
      y_rng = range(y_vals)
      y_pad = 0.05 * max(1e-8, diff(y_rng))
      y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
    }

    p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16) +
      ggplot2::geom_line(color = fit_col, linewidth = fit_lwd) +
      ggplot2::labs(title = main_title, x = x_info$x_var, y = ylab_text)

    if(!is.null(y_lim)) {
      p = p + ggplot2::coord_cartesian(ylim = y_lim)
    }

    if(show_legend) {
      p = p + ggplot2::labs(caption = paste("fit /", round(ci_level * 100), "% CI"))
    }

    p = p + ggplot2::theme_minimal()

    plot_objects[[length(plot_objects) + 1]] = p

    if(is.null(out[[par_name]])) out[[par_name]] = list()
    out[[par_name]][[s_name]] = list(
      x = x,
      fitted = fitted_smooth,
      se = smooth_fit_se,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      plot = p
    )
  }

  if(length(plot_objects) > 0) {
    if(is.null(ncol)) {
      ncol = min(2, n_plots)
    }
    nrow = ceiling(length(plot_objects) / ncol)
    dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow)
    if(setup_mfrow) {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
      for(i_plot in seq_along(plot_objects)) {
        r = ((i_plot - 1) %/% ncol) + 1
        c = ((i_plot - 1) %% ncol) + 1
        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
      }
      grid::popViewport()
    }
    out$plots = plot_objects
    out$dashboard = dashboard
  }

  invisible(out)
}

#' Plot all fixed terms with confidence bands
#'
#' This utility plots fixed-effect term contributions for a fitted
#' `gamlss.longitudinal` object using coefficient uncertainty from
#' `vcov.gamlss.longitudinal()`.
#'
#' For each fixed-effect design-matrix column \eqn{x_j}, it plots
#' \eqn{x_j \hat{\beta}_j} with pointwise confidence bands
#' \eqn{x_j \hat{\beta}_j \pm z_{\alpha/2}\sqrt{x_j^2 \mathrm{Var}(\hat{\beta}_j)}}.
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param vcov_obj Optional output from `vcov(object, ...)`. If `NULL`, this is
#' computed internally with the analytical vcov path.
#' @param ci_level Confidence level for pointwise intervals.
#' @param ncol Number of plot columns (defaults to 2 or fewer if needed).
#' @param include_intercept Logical; include intercept columns in plots.
#' @param ci_col Color for confidence bands.
#' @param fit_col Color for fitted fixed-term line.
#' @param ci_lty Line type for confidence bands.
#' @param fit_lwd Line width for fitted fixed-term line.
#' @param sort_x Logical; sort x-values before drawing lines.
#' @param fallback_to_index Logical; if x has one unique value, use index on x-axis.
#' @param setup_mfrow Logical; if TRUE (default), configure par(mfrow) inside
#' this function. Set FALSE when caller configures layout.
#' @param data Optional data frame used to detect factor columns and show
#' factor levels on x-axis for categorical fixed terms.
#' @param factor_pch Point symbol for factor-level estimates.
#' @param factor_cex Point size for factor-level estimates.
#' @param show_legend Logical; if TRUE, draw a small legend in each panel.
#'
#' @return Invisibly returns a nested list with x, fitted values, standard
#' errors, and confidence limits for each fixed term.
#' @export
plot_fixed_terms = function(
  object,
  vcov_obj = NULL,
  ci_level = 0.95,
  ncol = NULL,
  include_intercept = FALSE,
  plot_interactions = FALSE,
  ci_col = "red",
  fit_col = "black",
  ci_lty = 2,
  fit_lwd = 2,
  sort_x = TRUE,
  fallback_to_index = TRUE,
  setup_mfrow = TRUE,
  data = NULL,
  factor_pch = 16,
  factor_cex = 1.2,
  show_legend = TRUE
) {
  if(!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be of class 'gamlss.longitudinal'.")
  }

  if(is.null(vcov_obj)) {
    vcov_obj = .resolve_vcov(object, numderiv = FALSE, extra_args = list(method = "analytical"))
  }

  if(!is.list(vcov_obj) || is.null(vcov_obj$vcov) || is.null(vcov_obj$vcov$overall)) {
    stop("'vcov_obj' must be the output of vcov.gamlss.longitudinal() with vcov$overall present.")
  }

  V = vcov_obj$vcov$overall
  if(is.null(rownames(V)) || is.null(colnames(V))) {
    stop("vcov$overall must have row and column names matching fixed coefficients.")
  }

  z = qnorm((1 + ci_level) / 2)
  gg_add = function(plot, object, object_name = "") {
    ggplot2::ggplot_add(object, plot, object_name)
  }

  build_factor_groups = function(X, data) {
    groups = list()
    if(is.null(data) || !is.data.frame(data) || is.null(X) || ncol(X) == 0) return(groups)

    x_cols = colnames(X)
    for(var_name in names(data)) {
      v = data[[var_name]]
      if(!is.factor(v)) next
      levs = levels(v)
      if(length(levs) < 2) next

      var_tokens = strsplit(var_name, "_", fixed = TRUE)[[1]]
      var_prefixes = unique(c(
        var_name,
        make.names(var_name),
        if(length(var_tokens) > 0) var_tokens[1] else character(0),
        if(length(var_tokens) > 0) make.names(var_tokens[1]) else character(0)
      ))
      # Internal fitting often renames the user time variable to time_covariate.
      if(grepl("time", var_name, ignore.case = TRUE) && any(grepl("^time_covariate", x_cols))) {
        var_prefixes = unique(c(var_prefixes, "time_covariate", make.names("time_covariate")))
      }
      var_prefixes = var_prefixes[nzchar(var_prefixes)]

      level_col_map = list()
      matched_cols = character(0)
      for(lev in levs[-1]) {
        lev_plain = as.character(lev)
        lev_mn = make.names(lev_plain)
        candidates = unique(unlist(lapply(var_prefixes, function(pref) {
          c(
            paste0(pref, lev_plain),
            paste0(pref, lev_mn),
            paste0(pref, "_", lev_plain),
            paste0(pref, "_", lev_mn)
          )
        }), use.names = FALSE))
        hit = candidates[candidates %in% x_cols]
        if(length(hit) > 0) {
          level_col_map[[lev]] = hit[1]
          matched_cols = c(matched_cols, hit[1])
        }
      }

      if(length(matched_cols) == 0 && length(levs) == 2 && var_name %in% x_cols) {
        level_col_map[[levs[2]]] = var_name
        matched_cols = var_name
      }

      if(length(level_col_map) > 0) {
        groups[[var_name]] = list(
          var_name = var_name,
          levels = levs,
          ref_level = levs[1],
          level_col_map = level_col_map,
          matched_cols = unique(matched_cols)
        )
      }
    }
    groups
  }

  build_factor_interaction_groups = function(X, factor_groups) {
    groups = list()
    if(is.null(X) || ncol(X) == 0 || length(factor_groups) == 0) return(groups)

    x_cols = colnames(X)
    fg_names = names(factor_groups)
    if(length(fg_names) < 2) return(groups)

    for(i in seq_len(length(fg_names) - 1)) {
      for(j in (i + 1):length(fg_names)) {
        g1 = factor_groups[[fg_names[i]]]
        g2 = factor_groups[[fg_names[j]]]

        if(length(g1$level_col_map) == 0 || length(g2$level_col_map) == 0) next

        n1 = length(g1$levels)
        n2 = length(g2$levels)
        is_gender_1 = grepl("gender|sex", g1$var_name, ignore.case = TRUE)
        is_gender_2 = grepl("gender|sex", g2$var_name, ignore.case = TRUE)
        is_time_1 = grepl("time", g1$var_name, ignore.case = TRUE)
        is_time_2 = grepl("time", g2$var_name, ignore.case = TRUE)

        if(n1 > n2 || (n1 == n2 && is_time_1 && !is_time_2) || (n1 == n2 && !is_gender_1 && is_gender_2)) {
          panel_group = g1
          other_group = g2
        } else {
          panel_group = g2
          other_group = g1
        }

        panel_level_col_map = panel_group$level_col_map
        other_level_col_map = other_group$level_col_map

        interaction_col_map = list()
        matched_cols = character(0)
        for(panel_lev in names(panel_level_col_map)) {
          panel_col = panel_level_col_map[[panel_lev]]
          interaction_col_map[[panel_lev]] = list()

          for(other_lev in names(other_level_col_map)) {
            other_col = other_level_col_map[[other_lev]]
            candidates = c(
              paste0(other_col, ":", panel_col),
              paste0(panel_col, ":", other_col)
            )
            hit = candidates[candidates %in% x_cols]
            if(length(hit) > 0) {
              interaction_col_map[[panel_lev]][[other_lev]] = hit[1]
              matched_cols = c(matched_cols, hit[1])
            }
          }
        }

        if(length(matched_cols) > 0) {
          interaction_name = paste(other_group$var_name, panel_group$var_name, sep = ":")
          groups[[interaction_name]] = list(
            interaction_name = interaction_name,
            panel_group = panel_group,
            other_group = other_group,
            interaction_col_map = interaction_col_map,
            matched_cols = unique(matched_cols)
          )
        }
      }
    }

    groups
  }

  plot_specs = list()
  for(par_name in names(object$model_matrix$x)) {
    X = object$model_matrix$x[[par_name]]
    if(is.null(X) || ncol(X) == 0) next

    factor_groups = build_factor_groups(X, data)
    interaction_groups = if(plot_interactions) build_factor_interaction_groups(X, factor_groups) else list()
    grouped_cols = unique(unlist(lapply(factor_groups, function(g) g$matched_cols), use.names = FALSE))
    if(length(grouped_cols) == 0) grouped_cols = character(0)
    grouped_interaction_cols = unique(unlist(lapply(interaction_groups, function(g) g$matched_cols), use.names = FALSE))
    if(length(grouped_interaction_cols) == 0) grouped_interaction_cols = character(0)

    for(var_name in names(factor_groups)) {
      fg = factor_groups[[var_name]]
      has_valid_coef = FALSE
      for(lev in names(fg$level_col_map)) {
        coef_name = paste(par_name, fg$level_col_map[[lev]], sep = ".")
        if(coef_name %in% names(object$par) && coef_name %in% rownames(V) && coef_name %in% colnames(V)) {
          has_valid_coef = TRUE
          break
        }
      }
      if(has_valid_coef) {
        plot_specs[[length(plot_specs) + 1]] = list(
          type = "factor",
          par_name = par_name,
          var_name = var_name,
          group = fg
        )
      }
    }

    for(inter_name in names(interaction_groups)) {
      ig = interaction_groups[[inter_name]]
      has_valid_coef = FALSE
      for(panel_lev in names(ig$interaction_col_map)) {
        for(other_lev in names(ig$interaction_col_map[[panel_lev]])) {
          coef_name = paste(par_name, ig$interaction_col_map[[panel_lev]][[other_lev]], sep = ".")
          if(coef_name %in% names(object$par) && coef_name %in% rownames(V) && coef_name %in% colnames(V)) {
            has_valid_coef = TRUE
            break
          }
        }
        if(has_valid_coef) break
      }
      if(has_valid_coef) {
        plot_specs[[length(plot_specs) + 1]] = list(
          type = "interaction_factor_factor",
          par_name = par_name,
          group = ig
        )
      }
    }

    for(col_name in colnames(X)) {
      if(col_name %in% grouped_cols) next
      if(col_name %in% grouped_interaction_cols) next
      if(!include_intercept && col_name == "intercept") next
      if(!plot_interactions && grepl(":", col_name, fixed = TRUE)) next
      coef_name = paste(par_name, col_name, sep = ".")
      if(!coef_name %in% names(object$par)) next
      if(!coef_name %in% rownames(V) || !coef_name %in% colnames(V)) next

      term_type = if(grepl(":", col_name, fixed = TRUE)) "interaction_factor" else "continuous"

      plot_specs[[length(plot_specs) + 1]] = list(
        type = term_type,
        par_name = par_name,
        col_name = col_name,
        coef_name = coef_name
      )
    }
  }

  n_plots = length(plot_specs)
  if(n_plots == 0) {
    warning("No fixed terms found to plot with matching vcov entries.")
    return(invisible(list()))
  }

  out = list()
  plot_objects = list()
  for(i in seq_len(n_plots)) {
    spec = plot_specs[[i]]
    par_name = spec$par_name

    if(identical(spec$type, "factor")) {
      fg = spec$group
      levs = fg$levels
      x_plot = seq_along(levs)
      fitted_term = rep(NA_real_, length(levs))
      term_se = rep(NA_real_, length(levs))

      for(j in seq_along(levs)) {
        lev = levs[j]
        if(identical(lev, fg$ref_level)) {
          fitted_term[j] = 0
          term_se[j] = 0
        } else if(lev %in% names(fg$level_col_map)) {
          col_name_lev = fg$level_col_map[[lev]]
          coef_name_lev = paste(par_name, col_name_lev, sep = ".")
          if(coef_name_lev %in% names(object$par) && coef_name_lev %in% rownames(V) && coef_name_lev %in% colnames(V)) {
            fitted_term[j] = as.numeric(object$par[coef_name_lev])
            term_se[j] = sqrt(pmax(0, as.numeric(V[coef_name_lev, coef_name_lev])))
          }
        }
      }

      keep = is.finite(fitted_term) & is.finite(term_se)
      ci_lower = fitted_term - z * term_se
      ci_upper = fitted_term + z * term_se

      y_vals = c(fitted_term[keep], ci_lower[keep], ci_upper[keep])
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      plot_df = data.frame(
        x = x_plot,
        fitted = fitted_term,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        keep = keep
      )

      p = ggplot2::ggplot(plot_df[plot_df$keep, , drop = FALSE], ggplot2::aes(x = x, y = fitted))
      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
      p = gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")
      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")
      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = levs), "scale_x_continuous")
      p = gg_add(
        p,
        ggplot2::labs(
          title = paste(par_name, fg$var_name, sep = ": "),
          x = fg$var_name,
          y = paste("fixed contribution:", paste(par_name, fg$var_name, sep = "."))
        ),
        "labs"
      )
      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")

      if(!is.null(y_lim)) {
        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")
      }

      if(show_legend) {
        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")
      }

      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[fg$var_name]] = list(
        coefficient = paste(par_name, fg$var_name, sep = "."),
        x = x_plot,
        levels = levs,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )
      plot_objects[[length(plot_objects) + 1]] = p
    } else if(identical(spec$type, "interaction_factor_factor")) {
      ig = spec$group
      pg = ig$panel_group
      og = ig$other_group

      panel_levels = pg$levels
      other_levels = og$levels
      x_plot = seq_along(panel_levels)
      x_labels = panel_levels

      plot_rows = list()
      for(other_lev in other_levels) {
        fitted_term = rep(NA_real_, length(panel_levels))
        term_se = rep(NA_real_, length(panel_levels))

        for(j in seq_along(panel_levels)) {
          panel_lev = panel_levels[j]
          if(identical(panel_lev, pg$ref_level) || identical(other_lev, og$ref_level)) {
            fitted_term[j] = 0
            term_se[j] = 0
          } else if(panel_lev %in% names(ig$interaction_col_map) && other_lev %in% names(ig$interaction_col_map[[panel_lev]])) {
            col_name_lev = ig$interaction_col_map[[panel_lev]][[other_lev]]
            coef_name_lev = paste(par_name, col_name_lev, sep = ".")
            if(coef_name_lev %in% names(object$par) && coef_name_lev %in% rownames(V) && coef_name_lev %in% colnames(V)) {
              fitted_term[j] = as.numeric(object$par[coef_name_lev])
              term_se[j] = sqrt(pmax(0, as.numeric(V[coef_name_lev, coef_name_lev])))
            }
          }
        }

        keep = is.finite(fitted_term) & is.finite(term_se)
        ci_lower = fitted_term - z * term_se
        ci_upper = fitted_term + z * term_se

        plot_rows[[length(plot_rows) + 1]] = data.frame(
          x = x_plot,
          group = factor(rep(other_lev, length(panel_levels)), levels = other_levels),
          fitted = fitted_term,
          se = term_se,
          ci_lower = ci_lower,
          ci_upper = ci_upper,
          keep = keep,
          stringsAsFactors = FALSE
        )
      }

      plot_df = do.call(rbind, plot_rows)
      y_vals = c(plot_df$fitted[plot_df$keep], plot_df$ci_lower[plot_df$keep], plot_df$ci_upper[plot_df$keep])
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      p = ggplot2::ggplot(plot_df[plot_df$keep, , drop = FALSE], ggplot2::aes(x = x, y = fitted, color = group, group = group))
      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
      p = gg_add(p, ggplot2::geom_line(linewidth = 0.8), "geom_line")
      p = gg_add(p, ggplot2::geom_point(size = factor_cex), "geom_point")
      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), width = 0.15), "geom_errorbar")
      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = x_labels), "scale_x_continuous")
      p = gg_add(p, ggplot2::scale_color_discrete(name = og$var_name), "scale_color_discrete")
      p = gg_add(
        p,
        ggplot2::labs(
          title = paste(par_name, ig$interaction_name, sep = ": "),
          x = pg$var_name,
          y = paste("fixed contribution:", paste(par_name, ig$interaction_name, sep = "."))
        ),
        "labs"
      )
      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")

      if(!is.null(y_lim)) {
        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")
      }

      if(show_legend) {
        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")
      }

      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[ig$interaction_name]] = list(
        coefficient = paste(par_name, ig$interaction_name, sep = "."),
        x = x_plot,
        levels = panel_levels,
        series = other_levels,
        fitted = plot_df$fitted,
        se = plot_df$se,
        ci_lower = plot_df$ci_lower,
        ci_upper = plot_df$ci_upper,
        plot_data = plot_df,
        plot = p
      )
      plot_objects[[length(plot_objects) + 1]] = p
    } else if(identical(spec$type, "interaction_factor")) {
      col_name = spec$col_name
      coef_name = spec$coef_name
      X = object$model_matrix$x[[par_name]]
      x_raw = as.numeric(X[, col_name])
      beta_hat = as.numeric(object$par[coef_name])
      var_beta = as.numeric(V[coef_name, coef_name])

      x_levels = sort(unique(x_raw[is.finite(x_raw)]))
      if(length(x_levels) == 0) next

      fitted_term = x_levels * beta_hat
      term_se = abs(x_levels) * sqrt(pmax(0, var_beta))
      ci_lower = fitted_term - z * term_se
      ci_upper = fitted_term + z * term_se

      x_labels = as.character(signif(x_levels, 6))
      x_plot = seq_along(x_levels)

      y_vals = c(fitted_term, ci_lower, ci_upper)
      y_vals = y_vals[is.finite(y_vals)]
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      plot_df = data.frame(
        x = x_plot,
        fitted = fitted_term,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        stringsAsFactors = FALSE
      )

      p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted))
      p = gg_add(p, ggplot2::geom_hline(yintercept = 0, color = "grey70", linetype = 3), "geom_hline")
      p = gg_add(p, ggplot2::geom_point(color = fit_col, size = factor_cex), "geom_point")
      p = gg_add(p, ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), color = ci_col, width = 0.15), "geom_errorbar")
      p = gg_add(p, ggplot2::scale_x_continuous(breaks = x_plot, labels = x_labels), "scale_x_continuous")
      p = gg_add(
        p,
        ggplot2::labs(
          title = paste(par_name, col_name, sep = ": "),
          x = paste(col_name, "(interaction level)"),
          y = paste("fixed contribution:", coef_name)
        ),
        "labs"
      )
      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")

      if(!is.null(y_lim)) {
        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")
      }

      if(show_legend) {
        p = gg_add(p, ggplot2::labs(caption = paste("estimate /", round(ci_level * 100), "% CI")), "labs")
      }

      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[col_name]] = list(
        coefficient = coef_name,
        x = x_levels,
        levels = x_labels,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )
      plot_objects[[length(plot_objects) + 1]] = p
    } else {
      col_name = spec$col_name
      coef_name = spec$coef_name
      X = object$model_matrix$x[[par_name]]
      x_raw = as.numeric(X[, col_name])
      beta_hat = as.numeric(object$par[coef_name])
      var_beta = as.numeric(V[coef_name, coef_name])

      fitted_term = x_raw * beta_hat
      term_se = sqrt(pmax(0, (x_raw^2) * var_beta))
      ci_lower = fitted_term - z * term_se
      ci_upper = fitted_term + z * term_se

      if(length(unique(x_raw)) <= 1 && fallback_to_index) {
        x_plot = seq_along(x_raw)
        xlab_text = paste(col_name, "(index)")
        ord = seq_along(x_plot)
      } else {
        x_plot = x_raw
        xlab_text = col_name
        ord = if(sort_x) order(x_plot) else seq_along(x_plot)
      }

      main_title = paste(par_name, col_name, sep = ": ")
      ylab_text = paste("fixed contribution:", coef_name)

      y_vals = c(fitted_term, ci_lower, ci_upper)
      y_vals = y_vals[is.finite(y_vals)]
      if(length(y_vals) > 0) {
        y_rng = range(y_vals)
        y_pad = 0.05 * max(1e-8, diff(y_rng))
        y_lim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
      } else {
        y_lim = NULL
      }

      plot_df = data.frame(
        x = x_plot[ord],
        fitted = fitted_term[ord],
        ci_lower = ci_lower[ord],
        ci_upper = ci_upper[ord]
      )

      p = ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = fitted))
      p = gg_add(p, ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper), fill = ci_col, alpha = 0.16), "geom_ribbon")
      p = gg_add(p, ggplot2::geom_line(color = fit_col, linewidth = fit_lwd), "geom_line")
      p = gg_add(p, ggplot2::labs(title = main_title, x = xlab_text, y = ylab_text), "labs")
      p = gg_add(p, ggplot2::theme_minimal(), "theme_minimal")

      if(!is.null(y_lim)) {
        p = gg_add(p, ggplot2::coord_cartesian(ylim = y_lim), "coord_cartesian")
      }

      if(show_legend) {
        p = gg_add(p, ggplot2::labs(caption = paste("fit /", round(ci_level * 100), "% CI")), "labs")
      }

      if(is.null(out[[par_name]])) out[[par_name]] = list()
      out[[par_name]][[col_name]] = list(
        coefficient = coef_name,
        x = x_plot,
        fitted = fitted_term,
        se = term_se,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        plot = p
      )
      plot_objects[[length(plot_objects) + 1]] = p
    }
  }

  if(length(plot_objects) > 0) {
    if(is.null(ncol)) {
      ncol = min(2, n_plots)
    }
    nrow = ceiling(length(plot_objects) / ncol)
    dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow)
    if(setup_mfrow) {
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
      for(i_plot in seq_along(plot_objects)) {
        r = ((i_plot - 1) %/% ncol) + 1
        c = ((i_plot - 1) %% ncol) + 1
        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
      }
      grid::popViewport()
    }
    out$plots = plot_objects
    out$dashboard = dashboard
  }

  invisible(out)
}

#' Plot term effects for a fitted `gamlss.longitudinal` object
#'
#' This is the original term-wise plotting behavior (smooth and fixed effects)
#' that was previously exposed through `plot.gamlss.longitudinal()`.
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param y Unused; included for S3 generic compatibility.
#' @param data Optional data frame containing original covariates used to infer
#' smooth-term x-axis variables. If omitted, smooth terms may fall back to index.
#' @param ci_level Confidence level for pointwise intervals (default: 0.95).
#' @param ncol Number of columns in each plot frame (default: 2).
#' @param include_intercept Logical; include intercept terms in fixed plots.
#' @param ci_col Color for confidence bands (default: "red").
#' @param fit_col Color for fitted lines (default: "black").
#' @param show_legend Logical; if TRUE, draw a legend in each panel.
#' @param smooth_even_grid Logical; if TRUE, draw smooth terms on an evenly
#' spaced x-grid.
#' @param smooth_grid_n Number of x-grid points for smooth-term plots when
#' `smooth_even_grid = TRUE`.
#' @param paginate Logical; if `TRUE`, render one chart at a time and prompt
#' for Enter between plots (interactive sessions). If `FALSE` (default), render
#' all plots in a multi-panel dashboard.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list with `smooth_terms` and `fixed_terms`.
#' @export
plot.terms = function(x, ...) {
  UseMethod("plot.terms")
}

#' @method plot.terms gamlss.longitudinal
#' @export
plot.terms.gamlss.longitudinal = function(
  x,
  y,
  data = NULL,
  ci_level = 0.95,
  ncol = 4,
  include_intercept = FALSE,
  plot_interactions = FALSE,
  ci_col = "red",
  fit_col = "black",
  show_legend = TRUE,
  smooth_even_grid = TRUE,
  smooth_grid_n = 200,
  paginate = FALSE,
  ...
) {
  cat("\n=== Plotting term effects for gamlss.longitudinal object ===\n")

  count_plot_terms = function(obj) {
    n_smooth = 0
    n_fixed = 0

    for(par_name in names(obj$par_s)) {
      if(length(obj$par_s[[par_name]]) > 0) {
        n_smooth = n_smooth + length(obj$par_s[[par_name]])
      }
    }

    for(par_name in names(obj$model_matrix$x)) {
      X = obj$model_matrix$x[[par_name]]
      if(!is.null(X) && ncol(X) > 0) {
        coef_names = paste(par_name, colnames(X), sep = ".")
        coef_names = coef_names[!(colnames(X) == "intercept" & !include_intercept)]
        if(!plot_interactions) {
          coef_names = coef_names[!grepl(":", coef_names, fixed = TRUE)]
        }
        n_fixed = n_fixed + sum(coef_names %in% names(obj$par))
      }
    }

    list(smooth = n_smooth, fixed = n_fixed, total = n_smooth + n_fixed)
  }

  counts = count_plot_terms(x)
  cat(sprintf("Found %d smooth terms and %d fixed terms (total: %d plots).\n\n",
              counts$smooth, counts$fixed, counts$total))

  if(counts$total == 0) {
    warning("No term plots to display.")
    return(invisible(list(smooth_terms = list(), fixed_terms = list())))
  }

  vcov_obj = .resolve_vcov(x, numderiv = FALSE, extra_args = list(method = "analytical"))

  smooth_results = list()
  fixed_results = list()
  plot_objects = list()

  if(counts$smooth > 0) {
    smooth_results = plot_smooth_terms(
      object = x,
      vcov_obj = vcov_obj,
      data = data,
      ci_level = ci_level,
      ncol = ncol,
      ci_col = ci_col,
      fit_col = fit_col,
      even_grid = smooth_even_grid,
      grid_n = smooth_grid_n,
      setup_mfrow = FALSE,
      show_legend = show_legend
    )
    if(!is.null(smooth_results$plots)) {
      plot_objects = c(plot_objects, smooth_results$plots)
    }
  }

  if(counts$fixed > 0) {
    fixed_results = plot_fixed_terms(
      object = x,
      vcov_obj = vcov_obj,
      ci_level = ci_level,
      ncol = ncol,
      include_intercept = include_intercept,
      plot_interactions = plot_interactions,
      ci_col = ci_col,
      fit_col = fit_col,
      setup_mfrow = FALSE,
      data = data,
      show_legend = show_legend
    )
    if(!is.null(fixed_results$plots)) {
      plot_objects = c(plot_objects, fixed_results$plots)
    }
  }

  dashboard = NULL
  if(length(plot_objects) > 0) {
    if(length(plot_objects) > 16 && !isTRUE(paginate)) {
      warning(
        "More than 16 charts (", length(plot_objects), ") were generated. ",
        "Rendering all charts at once may fail in some environments. ",
        "Use paginate=TRUE to view one chart at a time.",
        call. = FALSE
      )
    }

    if(isTRUE(paginate)) {
      dashboard = list(plotlist = plot_objects, paginate = TRUE)
      for(i_plot in seq_along(plot_objects)) {
        grid::grid.newpage()
        print(plot_objects[[i_plot]])
        if(i_plot < length(plot_objects) && interactive()) {
          invisible(readline(prompt = sprintf("Press [Enter] for next chart (%d/%d)... ", i_plot, length(plot_objects))))
        }
      }
    } else {
      if(is.null(ncol)) {
        ncol = min(2, length(plot_objects))
      }
      nrow = ceiling(length(plot_objects) / ncol)
      dashboard = list(plotlist = plot_objects, ncol = ncol, nrow = nrow, paginate = FALSE)

      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow, ncol)))
      for(i_plot in seq_along(plot_objects)) {
        r = ((i_plot - 1) %/% ncol) + 1
        c = ((i_plot - 1) %% ncol) + 1
        print(plot_objects[[i_plot]], vp = grid::viewport(layout.pos.row = r, layout.pos.col = c))
      }
      grid::popViewport()
    }
  }

  invisible(list(
    smooth_terms = smooth_results,
    fixed_terms = fixed_results,
    dashboard = dashboard
  ))
}

#' Plot diagnostics dashboard for fitted `gamlss.longitudinal` objects
#'
#' Displays six ggplot-based panels by default:
#' 1) PIT histogram
#' 2) QQ residual plot
#' 3) Worm plot
#' 4) Rootogram
#' 5) Fitted-data quantile forecast plot
#' 6) Newdata quantile forecast plot (if `newdata` or `data` supplied)
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param y Unused; included for S3 generic compatibility.
#' @param data Optional data frame used as fallback for `newdata` plotting.
#' @param newdata Optional data frame for the newdata forecast panel.
#' @param newdata_n Number of rows to use from `data` when `newdata` is NULL.
#' @param quantiles Quantiles for forecast panels.
#' @param randomize Logical; randomized PIT/residual diagnostics.
#' @param time_stratified Logical; if TRUE, show time-stratified PIT and worm plots.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly returns a list of generated plot/data objects.
#' @export
plot.gamlss.longitudinal = function(
  x,
  y,
  data = NULL,
  newdata = NULL,
  newdata_n = 8,
  quantiles = c(0.1, 0.5, 0.9),
  randomize = TRUE,
  time_stratified = FALSE,
  ...
) {
  if(!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be of class 'gamlss.longitudinal'.")
  }

  q_col_name = function(prob) {
    paste0("q", gsub("^0\\.", "", format(prob, trim = TRUE)))
  }

  make_empty_plot = function(title_txt, msg_txt) {
    ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_text(label = msg_txt) +
      ggplot2::xlim(-1, 1) +
      ggplot2::ylim(-1, 1) +
      ggplot2::labs(title = title_txt) +
      ggplot2::theme_void()
  }

  # 1-4: Standard diagnostics
  p_diag1 = pithist(x, bins = 20, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag2 = qqrplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag3 = wormplot(x, randomize = randomize, plot = TRUE, by_time = time_stratified)
  p_diag4 = rootogram(x, bins = 20, plot = TRUE)

  # 5: Fitted-data forecast quantiles
  fc_fit_q = procast(x, type = "quantile", at = quantiles)
  fc_fit_q$idx = seq_len(nrow(fc_fit_q))

  q_low = q_col_name(min(quantiles))
  q_high = q_col_name(max(quantiles))
  q_mid = q_col_name(if(0.5 %in% quantiles) 0.5 else quantiles[ceiling(length(quantiles) / 2)])

  p_fit_quant = ggplot2::ggplot(fc_fit_q, ggplot2::aes(x = idx)) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_low, ymax = q_high), fill = "#4e79a7", alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes_string(y = q_mid), color = "#1f4e79", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(y = response), color = "black", alpha = 0.35, size = 0.9) +
    ggplot2::labs(
      title = "Fitted Forecast Quantiles",
      x = "Observation Index",
      y = "Response"
    ) +
    ggplot2::theme_minimal()

  # 6: Newdata forecast quantiles
  nd_use = newdata
  if(is.null(nd_use) && !is.null(data) && is.data.frame(data)) {
    nd_use = utils::head(data, newdata_n)
    # ensure quantile-only mode (response optional)
    if(is.null(x$var_map) || !"response" %in% x$var_map) {
      nd_use$response = NA_real_
    } else {
      response_orig = names(x$var_map)[x$var_map == "response"][1]
      if(!is.na(response_orig) && !response_orig %in% names(nd_use) && !"response" %in% names(nd_use)) {
        nd_use[[response_orig]] = NA_real_
      }
    }
  }

  p_new_quant = NULL
  fc_new_q = NULL
  if(!is.null(nd_use) && is.data.frame(nd_use) && nrow(nd_use) > 0) {
    fc_new_q = tryCatch(
      procast.gamlss.longitudinal(x, type = "quantile", at = quantiles, newdata = nd_use),
      error = function(e) NULL
    )

    if(!is.null(fc_new_q)) {
      time_candidates = c("time", if(!is.null(x$var_map)) names(x$var_map)[x$var_map == "time"] else character(0))
      person_candidates = c("subject", if(!is.null(x$var_map)) names(x$var_map)[x$var_map == "subject"] else character(0))
      time_col = time_candidates[time_candidates %in% names(nd_use)][1]
      person_col = person_candidates[person_candidates %in% names(nd_use)][1]

      if(is.na(time_col) || is.null(time_col) || nchar(time_col) == 0) {
        fc_new_q$time_plot = seq_len(nrow(fc_new_q))
      } else {
        fc_new_q$time_plot = nd_use[[time_col]]
      }
      if(is.na(person_col) || is.null(person_col) || nchar(person_col) == 0) {
        fc_new_q$person_plot = factor(seq_len(nrow(fc_new_q)))
      } else {
        fc_new_q$person_plot = as.factor(nd_use[[person_col]])
      }

      p_new_quant = ggplot2::ggplot(fc_new_q, ggplot2::aes(x = time_plot, color = person_plot, group = person_plot)) +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = q_low, ymax = q_high, fill = "person_plot"), alpha = 0.14, color = NA, show.legend = FALSE) +
        ggplot2::geom_line(ggplot2::aes_string(y = q_mid), linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes_string(y = q_mid), size = 1.7) +
        ggplot2::labs(
          title = "Newdata Forecast Quantiles",
          x = "Time",
          y = paste0("Predicted ", q_mid)
        ) +
        ggplot2::theme_minimal()
    }
  }

  if(is.null(p_new_quant)) {
    p_new_quant = make_empty_plot("Newdata Forecast Quantiles", "Provide 'newdata' or 'data' for this panel")
  }

  dashboard = ggpubr::ggarrange(
    p_diag1, p_diag2, p_diag3, p_diag4, p_fit_quant, p_new_quant,
    ncol = 2,
    nrow = 3
  )
  print(dashboard)

  invisible(list(
    diagnostics = list(pithist = p_diag1, qqrplot = p_diag2, wormplot = p_diag3, rootogram = p_diag4),
    forecasts = list(fitted_quantiles = p_fit_quant, newdata_quantiles = p_new_quant),
    fitted_data = fc_fit_q,
    newdata_data = fc_new_q,
    dashboard = dashboard
  ))
}

.copula_v2_clamp01 <- function(x) {
  pmin(pmax(x, 0.001), 0.999)
}

.copula_v2_tau_from_par <- function(family_num, par, par2 = NA_real_) {
  # Handle NA inputs immediately
  if (!is.finite(par)) {
    return(NA_real_)
  }

  tau <- tryCatch({
    if (is.finite(par2)) {
      suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = par2))
    } else {
      suppressWarnings(.copula_par_to_tau(family = family_num, par = par, par2 = 0))
    }
  }, error = function(e) NA_real_)

  if (is.finite(tau)) {
    return(as.numeric(tau))
  }

  # Fallback for Gaussian copulas using formula.
  if (identical(family_num, "N") && is.finite(par)) {
    return(2 / pi * asin(max(min(par, 0.999999), -0.999999)))
  }

  NA_real_
}

.copula_v2_bicop_cdf <- function(u1, u2, family_num, par, par2 = NA_real_) {
  n <- max(length(u1), length(u2), length(par), length(par2))
  if (!is.character(family_num) || length(family_num) != 1L || n < 1) return(rep(NA_real_, length(u1)))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  par <- rep(par, length.out = n)
  par2 <- rep(par2, length.out = n)
  vapply(seq_len(n), function(i) {
    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) return(NA_real_)
    tryCatch({
      .copula_cdf(
        u1[i],
        u2[i],
        family = family_num,
        par = par[i],
        par2 = if (is.finite(par2[i])) par2[i] else 0
      )
    }, error = function(e) NA_real_)
  }, numeric(1), USE.NAMES = FALSE)
}

.copula_v2_bicop_cond_u2_given_u1 <- function(u1, u2, family_num, par, par2 = NA_real_) {
  n <- max(length(u1), length(u2), length(par), length(par2))
  if (!is.character(family_num) || length(family_num) != 1L || n < 1) return(rep(NA_real_, length(u1)))
  u1 <- rep(u1, length.out = n)
  u2 <- rep(u2, length.out = n)
  par <- rep(par, length.out = n)
  par2 <- rep(par2, length.out = n)
  out <- vapply(seq_len(n), function(i) {
    if (!is.finite(u1[i]) || !is.finite(u2[i]) || !is.finite(par[i])) return(NA_real_)
    tryCatch({
      # BiCopHfunc1 gives dC(u1, u2) / du1, i.e. F(U2 <= u2 | U1 = u1).
      .copula_hfunc1(
        u1[i],
        u2[i],
        family = family_num,
        par = par[i],
        par2 = if (is.finite(par2[i])) par2[i] else 0
      )
    }, error = function(e) NA_real_)
  }, numeric(1), USE.NAMES = FALSE)
  .copula_v2_clamp01(as.numeric(out))
}

.copula_v2_rosenblatt_pair_data <- function(pair_data, family_num) {
  pair_data$rosenblatt <- .copula_v2_bicop_cond_u2_given_u1(
    pair_data$u1,
    pair_data$u2,
    family_num = family_num,
    par = pair_data$theta_pair,
    par2 = pair_data$zeta_pair
  )
  pair_data$z <- stats::qnorm(.copula_v2_clamp01(pair_data$rosenblatt))
  pair_data$z_prev <- stats::qnorm(.copula_v2_clamp01(pair_data$u1))
  pair_data$z_curr <- pair_data$z
  pair_data
}

.copula_v2_rosenblatt_series <- function(fit_data, family_num) {
  pair_data <- .copula_v2_pair_data(fit_data, lags = 1)
  pair_data <- .copula_v2_rosenblatt_pair_data(pair_data, family_num)

  out <- fit_data[, c("subject", "time", "u"), drop = FALSE]
  out$key <- paste(out$subject, as.character(out$time), sep = "::")
  out$rosenblatt <- NA_real_

  first_idx <- ave(seq_len(nrow(out)), out$subject, FUN = function(x) x == min(x))
  out$rosenblatt[as.logical(first_idx)] <- out$u[as.logical(first_idx)]

  pair_key <- paste(pair_data$subject, as.character(pair_data$time_right), sep = "::")
  out$rosenblatt <- ifelse(
    is.na(out$rosenblatt),
    pair_data$rosenblatt[match(out$key, pair_key)],
    out$rosenblatt
  )
  out$rosenblatt <- .copula_v2_clamp01(out$rosenblatt)
  out$z <- stats::qnorm(out$rosenblatt)
  out[is.finite(out$z), c("subject", "time", "rosenblatt", "z"), drop = FALSE]
}

.copula_v2_kendall_diagnostic <- function(pair_data, family_num) {
  if (nrow(pair_data) < 2) return(data.frame())

  emp_copula <- vapply(seq_len(nrow(pair_data)), function(i) {
    mean(pair_data$u1 <= pair_data$u1[i] & pair_data$u2 <= pair_data$u2[i], na.rm = TRUE)
  }, numeric(1))

  fit_copula <- .copula_v2_bicop_cdf(
    pair_data$u1,
    pair_data$u2,
    family_num = family_num,
    par = pair_data$theta_pair,
    par2 = pair_data$zeta_pair
  )

  ok <- is.finite(emp_copula) & is.finite(fit_copula)
  data.frame(
    empirical = sort(emp_copula[ok]),
    fitted = sort(fit_copula[ok]),
    stringsAsFactors = FALSE
  )
}

.copula_v2_tail_diagnostics <- function(pair_data, family_num, thresholds = c(0.05, 0.10, 0.20)) {
  rows <- lapply(thresholds, function(alpha) {
    c_lower <- .copula_v2_bicop_cdf(
      rep(alpha, nrow(pair_data)),
      rep(alpha, nrow(pair_data)),
      family_num = family_num,
      par = pair_data$theta_pair,
      par2 = pair_data$zeta_pair
    )
    upper_cut <- 1 - alpha
    c_upper_cut <- .copula_v2_bicop_cdf(
      rep(upper_cut, nrow(pair_data)),
      rep(upper_cut, nrow(pair_data)),
      family_num = family_num,
      par = pair_data$theta_pair,
      par2 = pair_data$zeta_pair
    )
    lower_fit <- mean(c_lower, na.rm = TRUE)
    upper_fit <- mean(1 - 2 * upper_cut + c_upper_cut, na.rm = TRUE)

    data.frame(
      threshold = alpha,
      tail = c("Lower", "Upper"),
      empirical = c(
        mean(pair_data$u1 <= alpha & pair_data$u2 <= alpha, na.rm = TRUE),
        mean(pair_data$u1 >= upper_cut & pair_data$u2 >= upper_cut, na.rm = TRUE)
      ),
      fitted = c(lower_fit, upper_fit),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.copula_v2_conditional_tail_diagnostics <- function(tail_df) {
  if (nrow(tail_df) == 0) return(tail_df)
  out <- tail_df
  out$empirical <- out$empirical / out$threshold
  out$fitted <- out$fitted / out$threshold
  out$empirical <- pmin(pmax(out$empirical, 0), 1)
  out$fitted <- pmin(pmax(out$fitted, 0), 1)
  out
}

.copula_v2_rosenblatt_lag_summary <- function(rosenblatt_df, lag_values = 1:3) {
  lag_values <- sort(unique(as.integer(lag_values)))
  lag_values <- lag_values[lag_values > 0]
  rows <- lapply(lag_values, function(lag_value) {
    pair_list <- lapply(split(rosenblatt_df, rosenblatt_df$subject), function(x) {
      x <- x[order(x$time), , drop = FALSE]
      if (nrow(x) <= lag_value) return(NULL)
      data.frame(
        z_prev = x$z[seq_len(nrow(x) - lag_value)],
        z_curr = x$z[(lag_value + 1):nrow(x)]
      )
    })
    pair_list <- pair_list[!vapply(pair_list, is.null, logical(1))]
    if (length(pair_list) == 0) {
      return(data.frame(lag = lag_value, cor_z = NA_real_, n_pairs = 0L))
    }
    pairs <- do.call(rbind, pair_list)
    data.frame(
      lag = lag_value,
      cor_z = suppressWarnings(stats::cor(pairs$z_prev, pairs$z_curr, use = "complete.obs")),
      n_pairs = nrow(pairs)
    )
  })
  do.call(rbind, rows)
}

.copula_v2_message_plot <- function(title, subtitle, message) {
  ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_text(label = message, size = 4) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::theme_void()
}

.copula_v2_fit_data <- function(object) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  eta_out <- calc_eta(
    par_cov = object$par,
    mm = object$model_matrix,
    margin_dist = object$margin_dist,
    copula_link = copula_spec$copula_link,
    par_s = object$par_s
  )

  response <- object$response
  subject <- object$response_subject
  time <- object$response_margin

  # Extract only margin parameters that are actually in eta_out$eta_inv
  margin_param_names <- names(object$margin_dist$parameters)
  margin_params <- list()
  for (param_name in margin_param_names) {
    if (param_name %in% names(eta_out$eta_inv)) {
      margin_params[[param_name]] <- eta_out$eta_inv[[param_name]]
    }
  }

  theta_fit <- if ("theta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$theta else numeric(0)
  zeta_fit <- if ("zeta" %in% names(eta_out$eta_inv)) eta_out$eta_inv$zeta else numeric(0)

  # Align response-side vectors to a common leading length.
  margin_min_n <- if (length(margin_params) > 0) {
    min(vapply(margin_params, length, integer(1)))
  } else {
    length(response)
  }

  common_n <- min(length(response), length(subject), length(time), margin_min_n)
  if (!is.finite(common_n) || common_n < 1) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  response <- response[seq_len(common_n)]
  subject <- subject[seq_len(common_n)]
  time <- time[seq_len(common_n)]
  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  align_copula_param <- function(param_vec) {
    n_resp <- common_n
    if (length(param_vec) == 0) {
      return(rep(NA_real_, n_resp))
    }

    # Full-row parameterization.
    if (length(param_vec) == n_resp) {
      return(param_vec)
    }

    # Pair-row parameterization: parameters correspond to times 1:(T-1) only.
    margin_names <- sort(unique(time))
    left_time_rows <- which(time %in% margin_names[seq_len(max(1, length(margin_names) - 1))])
    if (length(param_vec) == length(left_time_rows)) {
      out <- rep(NA_real_, n_resp)
      out[left_time_rows] <- param_vec
      return(out)
    }

    # Fallback for unexpected lengths.
    rep(param_vec, length.out = n_resp)
  }

  theta_fit <- align_copula_param(theta_fit)
  zeta_fit <- align_copula_param(zeta_fit)

  # Filter by finite values
  keep <- is.finite(response)
  for (param_name in names(margin_params)) {
    keep <- keep & is.finite(margin_params[[param_name]])
  }

  response <- response[keep]
  subject <- subject[keep]
  time <- time[keep]
  margin_params <- lapply(margin_params, function(x) x[keep])
  theta_fit <- theta_fit[keep]
  zeta_fit <- zeta_fit[keep]

  if (length(response) == 0) {
    stop("No finite fitted observations are available for copula diagnostics.")
  }

  # Convert margin_dist$family to family name if needed
  family_name <- object$margin_dist$family[1]
  if (!is.character(family_name)) {
    family_name <- object$margin_dist$family[1]$family
  }

  u <- .gl_call_family_fun("p", family_name, response, margin_params)
  u <- .copula_v2_clamp01(u)

  family_num <- tryCatch({
    .copula_family_code(copula_family_name)
  }, error = function(e) NA_character_)

  # Compute tau_fit, suppressing coercion warnings
  tau_fit <- suppressWarnings(
    vapply(seq_along(theta_fit), function(i) {
      .copula_v2_tau_from_par(family_num, theta_fit[i], zeta_fit[i])
    }, numeric(1), USE.NAMES = FALSE)
  )

  data.frame(
    subject = subject,
    time = time,
    response = response,
    u = u,
    theta_fit = theta_fit,
    zeta_fit = zeta_fit,
    tau_fit = tau_fit,
    stringsAsFactors = FALSE
  )
}

.copula_v2_pair_data <- function(fit_data, lags = 1) {
  time_vec <- fit_data$time
  time_levels <- if (is.factor(time_vec)) {
    lev <- levels(time_vec)
    lev[lev %in% as.character(unique(time_vec))]
  } else {
    u <- unique(time_vec)
    if (is.numeric(u) || is.integer(u)) sort(u) else sort(as.character(u))
  }
  if (length(time_levels) < 2) {
    stop("Need at least two time points to build copula pair diagnostics.")
  }

  time_lookup <- setNames(seq_along(time_levels), as.character(time_levels))
  fit_data$time_idx <- unname(time_lookup[as.character(fit_data$time)])
  if (any(!is.finite(fit_data$time_idx))) {
    stop("Could not map time values to an ordered index for copula pair diagnostics.")
  }

  lag_values <- sort(unique(as.integer(lags)))
  lag_values <- lag_values[lag_values > 0]
  if (length(lag_values) == 0) {
    lag_values <- 1L
  }

  pair_list <- list()
  idx <- 1L

  for (lag_value in lag_values) {
    for (subject_id in unique(fit_data$subject)) {
      subject_rows <- fit_data[fit_data$subject == subject_id, , drop = FALSE]
      subject_rows <- subject_rows[order(subject_rows$time_idx), , drop = FALSE]
      if (nrow(subject_rows) < 2) next

      for (j in seq_len(nrow(subject_rows) - lag_value)) {
        k <- j + lag_value
        if (k > nrow(subject_rows)) next

        t1 <- subject_rows$time[j]
        t2 <- subject_rows$time[k]
        t1_idx <- subject_rows$time_idx[j]
        t2_idx <- subject_rows$time_idx[k]
        if ((t2_idx - t1_idx) != lag_value) next

        row1 <- subject_rows[j, , drop = FALSE]
        row2 <- subject_rows[k, , drop = FALSE]

        # Match likelihood indexing: pair (t, t+lag) uses the left-row copula parameter.
        theta_pair <- as.numeric(row1$theta_fit)
        zeta_pair <- as.numeric(row1$zeta_fit)
        tau_pair <- as.numeric(row1$tau_fit)
        if (!is.finite(theta_pair)) theta_pair <- NA_real_
        if (!is.finite(zeta_pair)) zeta_pair <- NA_real_
        if (!is.finite(tau_pair)) tau_pair <- NA_real_

        pair_list[[idx]] <- data.frame(
          subject = subject_id,
          time_left = as.character(t1),
          time_right = as.character(t2),
          time_pair = paste0("T", as.character(t1), " vs T", as.character(t2)),
          lag = lag_value,
          u1 = row1$u,
          u2 = row2$u,
          theta_pair = theta_pair,
          zeta_pair = zeta_pair,
          tau_fit = tau_pair,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }

  if (length(pair_list) == 0) {
    stop("No complete subject-time pairs were found for copula diagnostics.")
  }

  do.call(rbind, pair_list)
}

.copula_v2_attach_group <- function(pair_data, object, by, data = NULL) {
  if (is.null(by) || (is.character(by) && length(by) == 1 && !nzchar(by))) {
    pair_data$split_group <- factor(pair_data$time_pair)
    return(pair_data)
  }

  if (!is.character(by) || length(by) != 1) {
    stop("'by' must be NULL or a single column name as a character string.")
  }

  if (by %in% c("time", "time_pair")) {
    pair_data$split_group <- factor(pair_data$time_pair)
    return(pair_data)
  }
  if (by %in% c("subject", "lag") && by %in% names(pair_data)) {
    pair_data$split_group <- factor(pair_data[[by]])
    return(pair_data)
  }

  if (is.null(data)) {
    stop("To split plot.copula by '", by, "', provide data= containing that column.")
  }

  df <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!is.null(object$var_map)) {
    for (old_name in names(object$var_map)) {
      new_name <- object$var_map[[old_name]]
      if (old_name %in% names(df) && !new_name %in% names(df)) {
        names(df)[names(df) == old_name] <- new_name
      }
    }
  }

  by_col <- by
  if (!by_col %in% names(df) && !is.null(object$var_map) && by %in% names(object$var_map)) {
    mapped_col <- object$var_map[[by]]
    if (mapped_col %in% names(df)) {
      by_col <- mapped_col
    }
  }
  if (!by_col %in% names(df)) {
    stop("Column '", by, "' not found in provided data after internal name mapping.")
  }
  if (!all(c("subject", "time") %in% names(df))) {
    stop("Provided data must contain subject and time columns (or names mappable via object$var_map) to split by '", by, "'.")
  }

  key_df <- paste(df$subject, as.character(df$time), sep = "::")
  key_pair <- paste(pair_data$subject, as.character(pair_data$time_left), sep = "::")
  matched <- df[[by_col]][match(key_pair, key_df)]

  if (all(is.na(matched))) {
    by_subj <- tapply(df[[by_col]], as.character(df$subject), function(v) {
      vv <- unique(v[!is.na(v)])
      if (length(vv) == 1) vv else NA
    })
    matched <- by_subj[as.character(pair_data$subject)]
  }

  pair_data$split_group <- factor(matched)
  pair_data <- pair_data[!is.na(pair_data$split_group), , drop = FALSE]

  if (nrow(pair_data) == 0) {
    stop("No valid paired rows remained after grouping by '", by, "'.")
  }

  pair_data
}

.copula_v2_transform_data <- function(data, transform = "uniform") {
  # Transform uniform [0,1] data to normal scale or other scales
  if (transform == "normal") {
    # Clamp to avoid infinite values from qnorm at 0 or 1
    data$u1 <- stats::qnorm(.copula_v2_clamp01(data$u1))
    data$u2 <- stats::qnorm(.copula_v2_clamp01(data$u2))
  }
  data
}

.copula_v2_average_density_grid <- function(family_num, pair_data, grid_n = 35, max_pairs_overlay = 300) {
  grid <- seq(0.02, 0.98, length.out = grid_n)
  grid_df <- expand.grid(u1 = grid, u2 = grid)

  pair_data <- pair_data[is.finite(pair_data$theta_pair), , drop = FALSE]
  if (nrow(pair_data) == 0) {
    grid_df$density <- NA_real_
    return(grid_df)
  }

  if (nrow(pair_data) > max_pairs_overlay) {
    set.seed(1)
    pair_data <- pair_data[sample(seq_len(nrow(pair_data)), max_pairs_overlay), , drop = FALSE]
  }

  density_sum <- rep(0, nrow(grid_df))
  density_count <- 0L

  for (i in seq_len(nrow(pair_data))) {
    par <- pair_data$theta_pair[i]
    par2 <- pair_data$zeta_pair[i]
    density_i <- tryCatch({
      if (is.finite(par2)) {
        .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = par2)
      } else {
        .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = 0)
      }
    }, error = function(e) rep(NA_real_, nrow(grid_df)))

    if (all(!is.finite(density_i))) next
    density_i[!is.finite(density_i)] <- 0
    density_sum <- density_sum + density_i
    density_count <- density_count + 1L
  }

  if (density_count == 0L) {
    grid_df$density <- NA_real_
  } else {
    grid_df$density <- density_sum / density_count
  }

  grid_df
}

.copula_v2_empirical_density_grid <- function(pair_data, grid_n = 35, lims = NULL) {
  if (is.null(lims)) {
    x_rng <- range(pair_data$u1, na.rm = TRUE)
    y_rng <- range(pair_data$u2, na.rm = TRUE)
    x_pad <- max(0.001, 0.025 * diff(x_rng))
    y_pad <- max(0.001, 0.025 * diff(y_rng))
    lims <- c(x_rng[1] - x_pad, x_rng[2] + x_pad, y_rng[1] - y_pad, y_rng[2] + y_pad)
  }
  kde <- MASS::kde2d(pair_data$u1, pair_data$u2, n = grid_n, lims = lims)
  data.frame(
    u1 = rep(kde$x, each = length(kde$y)),
    u2 = rep(kde$y, times = length(kde$x)),
    density = as.vector(kde$z),
    stringsAsFactors = FALSE
  )
}

.copula_v2_surface_metrics <- function(emp_density, fit_density, overlap_probs = c(0.7, 0.85, 0.95)) {
  emp <- as.numeric(emp_density)
  fit <- as.numeric(fit_density)
  ok <- is.finite(emp) & is.finite(fit)
  emp <- emp[ok]
  fit <- fit[ok]

  if (length(emp) == 0) {
    return(list(summary = data.frame(), overlap = data.frame()))
  }

  # Scale both surfaces to unit mass before computing distance metrics.
  emp <- pmax(emp, 0)
  fit <- pmax(fit, 0)
  emp <- emp / max(sum(emp), .Machine$double.eps)
  fit <- fit / max(sum(fit), .Machine$double.eps)

  summary_df <- data.frame(
    rmse = sqrt(mean((fit - emp)^2)),
    mae = mean(abs(fit - emp)),
    surface_cor = suppressWarnings(stats::cor(emp, fit, use = "complete.obs")),
    stringsAsFactors = FALSE
  )

  overlap_df <- do.call(rbind, lapply(overlap_probs, function(p) {
    thr_emp <- stats::quantile(emp, probs = p, na.rm = TRUE, type = 7)
    thr_fit <- stats::quantile(fit, probs = p, na.rm = TRUE, type = 7)
    mask_emp <- emp >= thr_emp
    mask_fit <- fit >= thr_fit
    union_n <- sum(mask_emp | mask_fit)
    iou <- if (union_n == 0) NA_real_ else sum(mask_emp & mask_fit) / union_n
    data.frame(level_prob = p, contour_iou = iou, stringsAsFactors = FALSE)
  }))

  list(summary = summary_df, overlap = overlap_df)
}

#' Compare fitted and empirical copula contour surfaces
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param grid_n Grid size used for density surfaces.
#' @param max_pairs_overlay Maximum number of paired observations used for fitted surface averaging.
#' @param contour_bins Number of contour levels to draw in the surface panels.
#' @param transform Character; "uniform" compares surfaces on copula scale, "normal" compares them on z-scale.
#' @param diff_scale_limit Positive numeric; fixed symmetric color scale limit for the difference panel.
#' @param time_stratified Logical; if TRUE, compare surfaces by time pair.
#' @param plot Logical; if TRUE, print the dashboard.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns plots, grid-level surfaces, and numeric similarity metrics.
#' @export
plot.copula_contour_compare <- function(x, lags = 1, grid_n = 45, max_pairs_overlay = 300, contour_bins = 10, transform = "uniform", diff_scale_limit = 0.05, time_stratified = FALSE, plot = TRUE, ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.")
  }
  object <- x

  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!is.numeric(diff_scale_limit) || length(diff_scale_limit) != 1 || !is.finite(diff_scale_limit) || diff_scale_limit <= 0) {
    stop("'diff_scale_limit' must be a single positive numeric value.")
  }

  fit_data <- .copula_v2_fit_data(object)
  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)

  copula_spec <- get_copula_dist(object$copula_dist)
  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch({
    .copula_family_code(copula_family_name)
  }, error = function(e) NA_character_)

  split_data <- if (isTRUE(time_stratified)) split(pair_data, pair_data$time_pair) else list(All = pair_data)

  grid_list <- lapply(names(split_data), function(nm) {
    pd <- split_data[[nm]]

    fit_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pd,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )

    # Build empirical surface on the same copula grid as fit_grid, then transform both
    # together if requested. This avoids grid mismatch artifacts in contouring.
    emp_grid <- .copula_v2_empirical_density_grid(pd, grid_n = grid_n, lims = c(0.02, 0.98, 0.02, 0.98))

    if (transform == "normal") {
      z1 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u1))
      z2 <- stats::qnorm(.copula_v2_clamp01(fit_grid$u2))
      jacobian <- stats::dnorm(z1) * stats::dnorm(z2)
      fit_grid$u1 <- z1
      fit_grid$u2 <- z2
      fit_grid$density <- fit_grid$density * jacobian

      emp_grid$u1 <- z1
      emp_grid$u2 <- z2
      emp_grid$density <- emp_grid$density * jacobian
    } else {
      emp_grid <- emp_grid
    }

    # Merge on grid coordinates to ensure pointwise comparisons.
    g <- merge(
      emp_grid,
      fit_grid,
      by = c("u1", "u2"),
      suffixes = c("_emp", "_fit"),
      all = FALSE
    )
    g$density_diff <- g$density_fit - g$density_emp
    g$time_pair <- nm
    g
  })

  grid_df <- do.call(rbind, grid_list)

  metric_list <- lapply(split(grid_df, grid_df$time_pair), function(g) {
    m <- .copula_v2_surface_metrics(g$density_emp, g$density_fit)
    out <- cbind(data.frame(time_pair = unique(g$time_pair), stringsAsFactors = FALSE), m$summary)
    if (nrow(m$overlap) > 0) {
      overlap <- cbind(data.frame(time_pair = unique(g$time_pair), stringsAsFactors = FALSE), m$overlap)
    } else {
      overlap <- data.frame()
    }
    list(summary = out, overlap = overlap)
  })

  metric_summary <- do.call(rbind, lapply(metric_list, function(x) x$summary))
  metric_overlap <- do.call(rbind, lapply(metric_list, function(x) x$overlap))

  x_label <- if (transform == "normal") expression(Phi^-1 * (U[t])) else expression(U[t])
  y_label <- if (transform == "normal") expression(Phi^-1 * (U[t + 1])) else expression(U[t + 1])

  p_emp <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_emp)) +
    ggplot2::geom_contour(color = "#4d4d4d", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Empirical Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_fit <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, z = density_fit)) +
    ggplot2::geom_contour(color = "#e41a1c", bins = contour_bins, linewidth = 0.9) +
    ggplot2::labs(title = "Fitted Copula Contours", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  p_diff <- ggplot2::ggplot(grid_df, ggplot2::aes(x = u1, y = u2, fill = density_diff)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-diff_scale_limit, diff_scale_limit),
      oob = scales::squish,
      name = "Fit - Emp"
    ) +
    ggplot2::labs(title = "Contour Difference Surface", x = x_label, y = y_label) +
    ggplot2::theme_minimal()

  if (isTRUE(time_stratified)) {
    p_emp <- p_emp + ggplot2::facet_wrap(~time_pair)
    p_fit <- p_fit + ggplot2::facet_wrap(~time_pair)
    p_diff <- p_diff + ggplot2::facet_wrap(~time_pair)
  }

  dashboard <- ggpubr::ggarrange(p_emp, p_fit, p_diff, ncol = 1, nrow = 3)

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    plots = list(empirical_contours = p_emp, fitted_contours = p_fit, difference_surface = p_diff),
    dashboard = dashboard,
    grid = grid_df,
    metrics = list(summary = metric_summary, overlap = metric_overlap)
  ))
}

#' Plot copula diagnostics for a fitted gamlss.longitudinal object
#'
#' @param x A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param grid_n Grid size used for contour averaging.
#' @param max_pairs_overlay Maximum number of paired observations used for the fitted overlay.
#' @param transform Character; "uniform" (default) shows empirical copula on [0,1], "normal" transforms to standard normal scale.
#' @param plot1_style Character; "bins" (default) draws a binned empirical layer, "scatter" draws points.
#' @param contour_bins Integer number of contour levels for the fitted copula overlay in plot 1.
#' @param time_stratified Logical; if TRUE, facet both plots by time pair.
#' @param by Optional grouping variable name for stratified plots. Defaults to
#'   time-pair grouping when NULL. Use `data` for covariates not stored on the
#'   fitted pair object (for example gender).
#' @param data Optional data frame used when grouping by a covariate via `by`.
#' @param tau_ylim Optional numeric vector of length 2 specifying y-axis limits
#'   for Kendall's tau chart(s). If `NULL` (default), y-axis scales are automatic.
#' @param plot2_cuts Integer number of quantile-based cuts used in plot 2 (default 10).
#' @param tail_thresholds Numeric vector of lower-tail probabilities used for
#'   tail co-occurrence and conditional exceedance diagnostics.
#' @param residual_lags Integer lags used for Rosenblatt normal-score
#'   autocorrelation diagnostics.
#' @param dashboard_ncol Number of columns in the combined diagnostic dashboard.
#' @param plot Logical; if TRUE, print the dashboard.
#' @param ... Additional arguments reserved for future methods.
#'
#' @return Invisibly returns a list with plot objects and summaries.
#' @export
plot.copula <- function(x, lags = 1, grid_n = 35, max_pairs_overlay = 300, transform = "normal", plot1_style = "bins", contour_bins = 8, time_stratified = FALSE, by = NULL, data = NULL, tau_ylim = NULL, plot2_cuts = 10, tail_thresholds = c(0.05, 0.10, 0.20), residual_lags = 1:3, dashboard_ncol = 2, plot = TRUE, ...) {
  if (!inherits(x, "gamlss.longitudinal")) {
    stop("'x' must be a fitted 'gamlss.longitudinal' object.")
  }
  object <- x

  # Validate and apply transformation
  if (!transform %in% c("uniform", "normal")) {
    stop("'transform' must be either 'uniform' or 'normal'.")
  }

  if (!plot1_style %in% c("bins", "scatter")) {
    stop("'plot1_style' must be either 'bins' or 'scatter'.")
  }

  if (!is.numeric(contour_bins) || length(contour_bins) != 1 || !is.finite(contour_bins) || contour_bins < 1) {
    stop("'contour_bins' must be a single finite number >= 1.")
  }
  contour_bins <- as.integer(round(contour_bins))

  if (!is.logical(time_stratified) || length(time_stratified) != 1 || is.na(time_stratified)) {
    stop("'time_stratified' must be TRUE or FALSE.")
  }

  if (!is.numeric(plot2_cuts) || length(plot2_cuts) != 1 || !is.finite(plot2_cuts) || plot2_cuts < 2) {
    stop("'plot2_cuts' must be a single finite number >= 2.")
  }
  plot2_cuts <- as.integer(round(plot2_cuts))

  if (!is.null(tau_ylim)) {
    if (!is.numeric(tau_ylim) || length(tau_ylim) != 2 || any(!is.finite(tau_ylim)) || tau_ylim[1] >= tau_ylim[2]) {
      stop("'tau_ylim' must be NULL or a numeric vector of length 2 with tau_ylim[1] < tau_ylim[2].")
    }
    tau_ylim <- as.numeric(tau_ylim)
  }

  tail_thresholds <- sort(unique(as.numeric(tail_thresholds)))
  tail_thresholds <- tail_thresholds[is.finite(tail_thresholds) & tail_thresholds > 0 & tail_thresholds < 0.5]
  if (length(tail_thresholds) == 0) {
    tail_thresholds <- c(0.05, 0.10, 0.20)
  }

  residual_lags <- sort(unique(as.integer(residual_lags)))
  residual_lags <- residual_lags[residual_lags > 0]
  if (length(residual_lags) == 0) {
    residual_lags <- 1:3
  }

  if (!is.numeric(dashboard_ncol) || length(dashboard_ncol) != 1 || !is.finite(dashboard_ncol) || dashboard_ncol < 1) {
    stop("'dashboard_ncol' must be a single positive integer.")
  }
  dashboard_ncol <- as.integer(round(dashboard_ncol))

  fit_data <- .copula_v2_fit_data(object)
  pair_data_uniform <- .copula_v2_pair_data(fit_data, lags = lags)

  if (isTRUE(time_stratified) && is.null(by)) {
    by <- "time_pair"
  } else if (isTRUE(time_stratified) && !is.null(by)) {
    warning("Both time_stratified and by were supplied; using by='", by, "'.", call. = FALSE)
  }

  pair_data_uniform <- .copula_v2_attach_group(pair_data_uniform, object = object, by = by, data = data)
  pair_data_plot <- pair_data_uniform

  # Apply transform to pair data if requested
  if (transform == "normal") {
    pair_data_plot <- .copula_v2_transform_data(pair_data_plot, transform = "normal")
  }

  copula_spec <- get_copula_dist(object$copula_dist)

  copula_family_name <- .copula_family_code(copula_spec$copula_dist)

  family_num <- tryCatch({
    .copula_family_code(copula_family_name)
  }, error = function(e) NA_character_)

  is_grouped <- !is.null(by) || isTRUE(time_stratified)

  if (is_grouped) {
    density_list <- lapply(split(pair_data_uniform, pair_data_uniform$split_group), function(x) {
      grid_i <- .copula_v2_average_density_grid(
        family_num = family_num,
        pair_data = x,
        grid_n = grid_n,
        max_pairs_overlay = max_pairs_overlay
      )
      grid_i$split_group <- as.character(x$split_group[1])
      grid_i
    })
    density_grid <- do.call(rbind, density_list)
  } else {
    density_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pair_data_uniform,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )
  }

  # Apply transform to density grid if requested
  if (transform == "normal") {
    # Transform coordinates to normal scale
    z1 <- stats::qnorm(.copula_v2_clamp01(density_grid$u1))
    z2 <- stats::qnorm(.copula_v2_clamp01(density_grid$u2))

    # Apply Jacobian correction: multiply by phi(z1) * phi(z2)
    # where phi is the standard normal PDF
    jacobian_correction <- stats::dnorm(z1) * stats::dnorm(z2)

    density_grid$u1 <- z1
    density_grid$u2 <- z2
    density_grid$density <- density_grid$density * jacobian_correction
  }

  # Set axis labels based on transform
  x_label <- if (transform == "normal") {
    expression(Phi^-1 * (U[t]))
  } else {
    expression(U[t])
  }

  y_label <- if (transform == "normal") {
    expression(Phi^-1 * (U[t + 1]))
  } else {
    expression(U[t + 1])
  }

  p1 <- ggplot2::ggplot(pair_data_plot, ggplot2::aes_string(x = "u1", y = "u2"))

  if (plot1_style == "scatter") {
    p1 <- p1 +
      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.45, size = 1.2)
  } else {
    p1 <- p1 +
      ggplot2::geom_bin2d(bins = 25, alpha = 0.8) +
      ggplot2::scale_fill_gradient(low = "white", high = "black", name = "Count")
  }

  p1 <- p1 +
    ggplot2::geom_contour(
      data = density_grid,
      ggplot2::aes(x = u1, y = u2, z = density),
      inherit.aes = FALSE,
      color = "#e41a1c",
      linewidth = 1.2,
      bins = contour_bins
    ) +
    ggplot2::labs(
      title = "Empirical Copula with Fitted Overlay",
      subtitle = paste0("Red contours show the average fitted copula across paired observations; family = ", copula_family_name),
      x = x_label,
      y = y_label
    ) +
    ggplot2::theme_minimal()

  if (is_grouped) {
    p1 <- p1 + ggplot2::facet_wrap(~split_group)
  }

  if (all(!is.finite(density_grid$density))) {
    p1 <- .copula_v2_message_plot(
      title = "Empirical Copula with Fitted Overlay",
      subtitle = paste0("Red contours show the average fitted copula across paired observations; family = ", copula_family_name),
      message = "No finite fitted copula density"
    )
  }

  build_cut_summary <- function(df, split_name = NULL) {
    if (nrow(df) < 1) {
      return(data.frame())
    }

    # Use rank-based bins to avoid collapsed quantile cuts when many fitted tau values are tied.
    df <- df[is.finite(df$tau_fit), , drop = FALSE]
    if (nrow(df) < 1) {
      return(data.frame())
    }

    effective_cuts <- min(plot2_cuts, nrow(df))
    cut_labels <- paste0("C", seq_len(effective_cuts))
    tau_rank <- rank(df$tau_fit, ties.method = "first", na.last = "keep")
    df$cut_group <- cut(tau_rank, breaks = effective_cuts, include.lowest = TRUE, labels = cut_labels)

    out <- do.call(rbind, lapply(split(df, df$cut_group), function(x) {
      tau_emp <- suppressWarnings(stats::cor(x$u1, x$u2, method = "kendall", use = "complete.obs"))
      tau_fit <- mean(x$tau_fit, na.rm = TRUE)
      data.frame(
        cut_group = as.character(x$cut_group[1]),
        tau_emp = tau_emp,
        tau_fit = tau_fit,
        n_pairs = nrow(x),
        stringsAsFactors = FALSE
      )
    }))

    if (!is.null(split_name)) {
      out$split_group <- split_name
    }
    out
  }

  if (is_grouped) {
    quartile_list <- lapply(split(pair_data_plot, pair_data_plot$split_group), function(x) {
      build_cut_summary(x, split_name = as.character(x$split_group[1]))
    })
    quartile_df <- do.call(rbind, quartile_list)
  } else {
    quartile_df <- build_cut_summary(pair_data_plot)
  }

  if (nrow(quartile_df) == 0 || all(!is.finite(quartile_df$tau_emp)) || all(!is.finite(quartile_df$tau_fit))) {
    p2 <- .copula_v2_message_plot(
      title = "Observed vs Fitted Correlation by Quantile Bin",
      subtitle = "Bins are formed from fitted copula strength",
      message = "No finite cut summaries"
    )
  } else {
    cut_levels <- paste0("C", sort(unique(as.integer(sub("^C", "", quartile_df$cut_group)))))
    quartile_df$cut_group <- factor(quartile_df$cut_group, levels = cut_levels)
    p2 <- ggplot2::ggplot(quartile_df, ggplot2::aes(x = cut_group)) +
      ggplot2::geom_point(ggplot2::aes(y = tau_emp), color = "#4d4d4d", size = 2.8) +
      ggplot2::geom_line(ggplot2::aes(y = tau_emp, group = 1), color = "#4d4d4d", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(y = tau_fit), color = "#e41a1c", size = 2.8, shape = 4, stroke = 1.1) +
      ggplot2::geom_line(ggplot2::aes(y = tau_fit, group = 1), color = "#e41a1c", linewidth = 0.8, linetype = "dashed") +
      ggplot2::labs(
        title = "Observed vs Fitted Correlation by Quantile Bin",
        subtitle = paste0("", plot2_cuts, " cuts formed from fitted copula strength"),
        x = "Cut",
        y = "Kendall's tau"
      ) +
      ggplot2::theme_minimal()

    if (is_grouped) {
      p2 <- p2 + ggplot2::facet_wrap(~split_group, scales = if (is.null(tau_ylim)) "free_y" else "fixed")
    }

    if (!is.null(tau_ylim)) {
      p2 <- p2 + ggplot2::coord_cartesian(ylim = tau_ylim)
    }
  }

  rosenblatt_df <- tryCatch(
    .copula_v2_rosenblatt_series(fit_data, family_num),
    error = function(e) data.frame()
  )

  rosenblatt_pair_df <- tryCatch(
    .copula_v2_rosenblatt_pair_data(pair_data_uniform, family_num),
    error = function(e) data.frame()
  )

  if (nrow(rosenblatt_df) == 0 || all(!is.finite(rosenblatt_df$z))) {
    p_ros_time <- .copula_v2_message_plot(
      title = "Rosenblatt Normal Scores by Time",
      subtitle = "Scores are qnorm of pairwise conditional Rosenblatt residuals",
      message = "No finite Rosenblatt residuals"
    )
  } else {
    p_ros_time <- ggplot2::ggplot(rosenblatt_df, ggplot2::aes(x = factor(time), y = z)) +
      ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
      ggplot2::geom_boxplot(fill = "#9ecae1", color = "#4d4d4d", outlier.alpha = 0.35) +
      ggplot2::labs(
        title = "Rosenblatt Normal Scores by Time",
        subtitle = "Each time point should be centered near zero with similar spread",
        x = "Time",
        y = "Normal score"
      ) +
      ggplot2::theme_minimal()
  }

  if (nrow(rosenblatt_pair_df) == 0 || all(!is.finite(rosenblatt_pair_df$z_prev)) || all(!is.finite(rosenblatt_pair_df$z_curr))) {
    p_ros_lag <- .copula_v2_message_plot(
      title = "Rosenblatt Lag Plot",
      subtitle = "Current conditional score against previous marginal score",
      message = "No finite Rosenblatt lag pairs"
    )
  } else {
    p_ros_lag <- ggplot2::ggplot(rosenblatt_pair_df, ggplot2::aes(x = z_prev, y = z_curr)) +
      ggplot2::geom_hline(yintercept = 0, color = "#d9d9d9") +
      ggplot2::geom_vline(xintercept = 0, color = "#d9d9d9") +
      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.35, size = 1.1) +
      ggplot2::geom_smooth(method = "loess", se = FALSE, color = "#e41a1c", linewidth = 0.7) +
      ggplot2::labs(
        title = "Rosenblatt Lag Plot",
        subtitle = "The smooth should be approximately flat at zero",
        x = expression(Phi^-1 * (U[t])),
        y = expression(Phi^-1 * (R[t + 1] ~ "|" ~ U[t]))
      ) +
      ggplot2::theme_minimal()
  }

  kendall_df <- tryCatch(
    .copula_v2_kendall_diagnostic(pair_data_uniform, family_num),
    error = function(e) data.frame()
  )

  if (nrow(kendall_df) == 0) {
    p_kendall <- .copula_v2_message_plot(
      title = "Kendall Function Diagnostic",
      subtitle = "Empirical copula values compared with fitted copula values at observed pairs",
      message = "No finite Kendall diagnostic values"
    )
  } else {
    p_kendall <- ggplot2::ggplot(kendall_df, ggplot2::aes(x = fitted, y = empirical)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "#666666", linetype = "dashed") +
      ggplot2::geom_point(color = "#4d4d4d", alpha = 0.55, size = 1.2) +
      ggplot2::labs(
        title = "Kendall Function Diagnostic",
        subtitle = "Sorted empirical copula probabilities should track sorted fitted probabilities",
        x = "Fitted copula probability",
        y = "Empirical copula probability"
      ) +
      ggplot2::theme_minimal()
  }

  tail_df <- tryCatch(
    .copula_v2_tail_diagnostics(pair_data_uniform, family_num, thresholds = tail_thresholds),
    error = function(e) data.frame()
  )
  cond_tail_df <- .copula_v2_conditional_tail_diagnostics(tail_df)

  tail_long <- if (nrow(tail_df) > 0) {
    rbind(
      data.frame(threshold = tail_df$threshold, tail = tail_df$tail, source = "Empirical", probability = tail_df$empirical),
      data.frame(threshold = tail_df$threshold, tail = tail_df$tail, source = "Fitted", probability = tail_df$fitted)
    )
  } else {
    data.frame()
  }

  if (nrow(tail_long) == 0 || all(!is.finite(tail_long$probability))) {
    p_tail <- .copula_v2_message_plot(
      title = "Tail Co-occurrence",
      subtitle = "Observed joint tail probability against fitted copula probability",
      message = "No finite tail diagnostics"
    )
  } else {
    p_tail <- ggplot2::ggplot(tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~tail) +
      ggplot2::scale_color_manual(values = c(Empirical = "#4d4d4d", Fitted = "#e41a1c")) +
      ggplot2::labs(
        title = "Tail Co-occurrence",
        subtitle = "Lower: P(Ut <= a, Ut+1 <= a); Upper: P(Ut >= 1-a, Ut+1 >= 1-a)",
        x = "Tail probability a",
        y = "Joint probability",
        color = NULL
      ) +
      ggplot2::theme_minimal()
  }

  cond_tail_long <- if (nrow(cond_tail_df) > 0) {
    rbind(
      data.frame(threshold = cond_tail_df$threshold, tail = cond_tail_df$tail, source = "Empirical", probability = cond_tail_df$empirical),
      data.frame(threshold = cond_tail_df$threshold, tail = cond_tail_df$tail, source = "Fitted", probability = cond_tail_df$fitted)
    )
  } else {
    data.frame()
  }

  if (nrow(cond_tail_long) == 0 || all(!is.finite(cond_tail_long$probability))) {
    p_cond_tail <- .copula_v2_message_plot(
      title = "Conditional Tail Exceedance",
      subtitle = "Observed conditional tail probability against fitted copula probability",
      message = "No finite conditional tail diagnostics"
    )
  } else {
    p_cond_tail <- ggplot2::ggplot(cond_tail_long, ggplot2::aes(x = threshold, y = probability, color = source, group = source)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~tail) +
      ggplot2::scale_color_manual(values = c(Empirical = "#4d4d4d", Fitted = "#e41a1c")) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(
        title = "Conditional Tail Exceedance",
        subtitle = "Lower: P(Ut+1 <= a | Ut <= a); Upper: P(Ut+1 >= 1-a | Ut >= 1-a)",
        x = "Tail probability a",
        y = "Conditional probability",
        color = NULL
      ) +
      ggplot2::theme_minimal()
  }

  lag_summary_df <- tryCatch(
    .copula_v2_rosenblatt_lag_summary(rosenblatt_df, lag_values = residual_lags),
    error = function(e) data.frame()
  )

  if (nrow(lag_summary_df) == 0 || all(!is.finite(lag_summary_df$cor_z))) {
    p_lag_summary <- .copula_v2_message_plot(
      title = "Residual Dependence by Lag",
      subtitle = "Correlation of Rosenblatt normal scores within subject",
      message = "No finite residual lag correlations"
    )
  } else {
    p_lag_summary <- ggplot2::ggplot(lag_summary_df, ggplot2::aes(x = factor(lag), y = cor_z)) +
      ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed") +
      ggplot2::geom_col(fill = "#4d4d4d", alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = paste0("n=", n_pairs)), vjust = -0.35, size = 3) +
      ggplot2::labs(
        title = "Residual Dependence by Lag",
        subtitle = "Correlations should be close to zero after the Rosenblatt transform",
        x = "Lag",
        y = "Correlation"
      ) +
      ggplot2::theme_minimal()
  }

  dashboard_plots <- list(p1, p2, p_ros_time, p_ros_lag, p_kendall, p_tail, p_cond_tail, p_lag_summary)
  dashboard <- do.call(
    ggpubr::ggarrange,
    c(
      dashboard_plots,
      list(
        ncol = min(dashboard_ncol, length(dashboard_plots)),
        nrow = ceiling(length(dashboard_plots) / dashboard_ncol)
      )
    )
  )

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    plots = list(
      empirical_overlay = p1,
      quartile_correlation = p2,
      rosenblatt_by_time = p_ros_time,
      rosenblatt_lag = p_ros_lag,
      kendall_function = p_kendall,
      tail_cooccurrence = p_tail,
      conditional_tail_exceedance = p_cond_tail,
      residual_lag_correlation = p_lag_summary
    ),
    dashboard = dashboard,
    fit_data = fit_data,
    pair_data = pair_data_plot,
    pair_data_uniform = pair_data_uniform,
    rosenblatt = rosenblatt_df,
    rosenblatt_pairs = rosenblatt_pair_df,
    quartile_summary = quartile_df,
    kendall_summary = kendall_df,
    tail_summary = tail_df,
    conditional_tail_summary = cond_tail_df,
    residual_lag_summary = lag_summary_df
  ))
}

#' @export
eta_to_par=function(eta,margin_dist,copula_dist) {
  par=eta*0
  for (par_name in names(eta)) {
    if(par_name %in% names(margin_dist$parameters)) {
      FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkinv",sep="."))))
      par[par_name]=FUN(eta[par_name])
    }
    if(par_name %in% names(copula_dist$parameters)) {
      FUN = eval(parse(text=paste(paste(paste("copula_dist$",par_name,sep=""),"linkinv",sep="."))))
      par[par_name]=FUN(eta[par_name])
    }
  }
  return(par)
}

.select_t_copula_zeta_start <- function(dataset, margin_dist, copula_dist, margin_par, theta_start) {
  # Grid is ordered low-to-high: return the first candidate that yields a finite
  # joint log-likelihood. Starting as low as possible avoids the optimizer being
  # trapped at high df (Gaussian limit) where the link-scale step sizes collapse.
  zeta_grid <- c(2.05, 2.2, 2.5, 3, 4, 5, 8, 12, 20, 35)
  fallback_zeta <- 3

  param_names <- c(names(margin_dist$parameters), get_copula_dist(copula_dist)$parameters)
  mm_stub <- as.list(setNames(rep(1, length(param_names)), param_names))
  pair_cache <- build_copula_pair_cache(dataset$response, dataset$time, dataset$subject)

  base_eta_inv <- c(as.list(margin_par), list(theta = as.numeric(theta_start)[1]))

  for (candidate_zeta in zeta_grid) {
    eta_inv <- base_eta_inv
    eta_inv$zeta <- candidate_zeta

    candidate_fit <- tryCatch(
      calc_likelihood_minimal(
        eta_inv = eta_inv,
        mm = mm_stub,
        margin_dist = margin_dist,
        copula_dist = copula_dist,
        calc_d2 = FALSE,
        response = dataset$response,
        response_margin = dataset$time,
        response_subject = dataset$subject,
        pair_cache = pair_cache
      ),
      error = function(e) NULL
    )

    if (!is.null(candidate_fit) && is.finite(candidate_fit$log_lik["joint"])) {
      return(candidate_zeta)
    }
  }

  fallback_zeta
}

#' @export
get_starting_values = function(copula_dist,margin_dist,dataset,eta_transform=FALSE) {

  margin_names=unique(dataset$time)
  num_margins=length(margin_names)
  finite_response <- dataset$response[is.finite(dataset$response)]
  margin_par_already_eta <- FALSE
  moment_skewness <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 3) return(0)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) return(0)
    mean(((x - mean(x)) / s)^3)
  }
  moment_kurtosis <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 4) return(3)
    s <- stats::sd(x)
    if (!is.finite(s) || s <= 0) return(3)
    mean(((x - mean(x)) / s)^4)
  }

  tau_start=cor(dataset[dataset$time%in%(margin_names[1:(num_margins-1)]),"response"]
                ,dataset[dataset$time%in%(margin_names[2:(num_margins)]),"response"],method="kendall",use="complete.obs")
  if(!is.finite(tau_start)) {
    warning("Non-finite Kendall tau in get_starting_values(); using tau = 0 for copula initialisation.")
    tau_start=0
  }
  tau_start=max(min(tau_start,0.9999),-0.9999)

  copula_spec=get_copula_dist(copula_dist)
  theta_start=.copula_tau_to_par(
    family=copula_dist,
    tau=tau_start
  )

  if(margin_dist$family[1]=="GA" | margin_dist$family[1]=="EXP") {
    margin_par=c(
      mean(finite_response)
      , stats::sd(finite_response)/mean(finite_response)
      , moment_skewness(finite_response)
      , moment_kurtosis(finite_response)
    )
  } else if (margin_dist$family[1]=="NO") {
    margin_par=c(
      mean(finite_response)
      , stats::sd(finite_response)
    )
  } else if (margin_dist$family[1]=="PO") {
    margin_par=c(
      mean(finite_response)
    )
  } else if (margin_dist$family[1]=="NBI") {
    margin_par=c(
      mean(finite_response),
      stats::sd(finite_response)/mean(finite_response)
    )
  } else {
    cat("Fitting initial GAMLSS model for margin to obtain starting values...\n")
    # Deliberately low-iteration startup fit; silence expected convergence warnings.
    start_fit=suppressWarnings(suppressMessages(
      gamlss(dataset$response~1, family=margin_dist)
    ))
    margin_par=unlist(coefAll(start_fit))
    names(margin_par)=names(margin_dist$parameters)
    #margin_par=eta_to_par(margin_par_temp,margin_dist,get_copula_dist(copula_dist))
    margin_par_already_eta <- TRUE
  }

  names(margin_par)=names(margin_dist$parameters)
  margin_par=margin_par[!is.na(names(margin_par))]

  if("zeta" %in% copula_spec$parameters) {
    # .copula_tau_to_par() returns only theta for t-copula; select zeta by a small grid search.
    zeta_start <- .select_t_copula_zeta_start(
      dataset = dataset,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      margin_par = margin_par,
      theta_start = theta_start
    )
    cop_par=c(theta=as.numeric(theta_start)[1], zeta=as.numeric(zeta_start))
  } else {
    cop_par=c(theta=as.numeric(theta_start)[1])
  }

  if(eta_transform==TRUE) {
    margin_par_eta=margin_par
    cop_par_eta=cop_par

    if(!isTRUE(margin_par_already_eta)) {
      for (par_name in names(margin_par)) {
        FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkfun",sep="."))))
        margin_par_eta[par_name]=FUN(margin_par[par_name])
      }
    }

    for (par_name in names(cop_par)) {
      cop_par_eta[par_name]=get_copula_dist(copula_dist)$copula_link[[paste(par_name,".linkfun",sep="")]](cop_par[par_name])
    }

    return_list=c(margin_par_eta,cop_par_eta)
  } else {
    return_list=c(margin_par,cop_par)
  }

  return(return_list)
}
#' @export
par_to_eta = function(par,copula_dist,margin_dist) {

  margin_par=par[names(margin_dist$parameters)]
  names(margin_par)=names(margin_dist$parameters)

  cop_par=par[get_copula_dist(copula_dist)$parameters]
  names(cop_par)=get_copula_dist(copula_dist)$parameters

    margin_par_eta=margin_par
    cop_par_eta=cop_par

    for (par_name in names(margin_par)) {
      FUN = eval(parse(text=paste(paste(paste("margin_dist$",par_name,sep=""),"linkfun",sep="."))))
      margin_par_eta[par_name]=FUN(margin_par[par_name])
    }

    for (par_name in names(cop_par)) {
      cop_par_eta[par_name]=get_copula_dist(copula_dist)$copula_link[[paste(par_name,".linkfun",sep="")]](cop_par[par_name])
      names(cop_par_eta)=names(cop_par)
    }

    return_list=c(margin_par_eta,cop_par_eta)

  return(return_list)
}

#' @export
fit_jointreg_nocov <- function(input_par,margin_dist,copula_dist,data
                               , use_dlcopdpar=TRUE, verbose=TRUE, plot_results=TRUE
                               , crit_lik_change=0.05, start_step_size=.5, step_adjustment=.9, max_steps=5
                               , true_val = NA) {

  log_lik_history=matrix(ncol=3+2,nrow=0)
  par_history=matrix(ncol=length(input_par)+2,nrow=0)

  ### Run fit for separate and joint optimisation
    copula_deriv=if(use_dlcopdpar==TRUE){1}else{0}
    ### CORE ITERATION
    change=1;log_lik_start=0;log_lik_change=1000;run_counter=1;step_size=start_step_size;
    while (abs(log_lik_change)>crit_lik_change) {
      step_size=step_size*(step_adjustment^min(max_steps,run_counter))
      par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))

      #Run optimisation
      outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,use_dlcopdpar=use_dlcopdpar,verbose=FALSE,step_size=step_size)

      #Capture outputs
      input_par=outer_optim_output$par_end
      change=sum(outer_optim_output$par_change)
      #print(outer_optim_output$log_lik)

      log_lik=outer_optim_output$log_lik["joint"]
      log_lik_change=log_lik-log_lik_start
      log_lik_start=log_lik

      #Capture changes in parameters
      log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))
      run_counter=run_counter+1

    }
    par_history=rbind(par_history,c(copula_deriv,run_counter,input_par))

    outer_optim_output=optim_outer(par=input_par,dataset,margin_dist,copula_dist,use_dlcopdpar=use_dlcopdpar,verbose=FALSE,step_size=step_size)
    log_lik=outer_optim_output$log_lik["joint"]
    log_lik_change=log_lik-log_lik_start
    log_lik_start=log_lik
    log_lik_history=rbind(log_lik_history,c(copula_deriv,run_counter,outer_optim_output$log_lik))

    colnames(log_lik_history)[1:2]=colnames(par_history)[1:2]=c("use_dlcopdpar","run_counter")

  #Plot likelihood and parameters
  if(plot_results==TRUE) {

    plot.new()
    par_count=round(sqrt((ncol(par_history)+1)),0)+1
    par(mfrow=c(par_count,par_count))

    for (i in colnames(log_lik_history)[3:5]) {
      plot( log_lik_history[,i],xlab="Iteration",ylab="LogLik",main=i,type = "l",col="blue",xlim=c(1,max(log_lik_history[,"run_counter"])))
      #lines(log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="red",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))
      #legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))

    }

    for (i in 1:(ncol(par_history)-2)) {
      #lines(par_nodlcop[,i+2],col="red",type="l")
      if (!all(is.na(true_val))) {
        plot(par_history[,i+2],col="blue",type="l",main=colnames(par_history)[i+2],ylab="Parameter estimate",ylim=range(c(par_history[,i+2],true_val[i])))
        abline(h=true_val[i])
      } else {
        plot(par_history[,i+2],col="blue",type="l",main=colnames(par_history)[i+2],ylab="Parameter estimate")
      }
      #legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    }
  }

  return_list=list(par_history,log_lik_history)
  names(return_list)=c("par_history","log_lik_history")
  return(return_list)
}




#' @keywords internal
#' @noRd
calc_F_x <- function(eta_inv,mm,margin_dist,response) {
  #Setup input matrix of response and parameters
  #margin_names=unique(response_margin)
  #num_margins=length(margin_names)

  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(mm)) {
    if (par_name %in% c("mu","sigma","nu","tau")) {
      margin_deriv_input[[par_name]]=eta_inv[[par_name]]
    }
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])

  return(margin_p)
}

#' @export
get_copula_dist=function(copula_dist) {

  copula_dist <- .copula_family_code(copula_dist)

  if(copula_dist=="C") {
    copula_link=list(log,exp,dloginv=exp); two_par_cop=FALSE
    parameters=c("theta")
  }
  else if(copula_dist=="F") {
    copula_link=list(identity,identity,function(x) rep(1, length(x))); two_par_cop=FALSE
    parameters=c("theta")
  }
  else if(copula_dist=="J") {
    copula_link=list(log_1plus,log_1plus_inv,dlog_1plus_inv); two_par_cop=FALSE
    parameters=c("theta")
  }
  else if(copula_dist=="G") {
    copula_link=list(gumbel_linkfun,gumbel_linkinv,dgumbel_linkinv); two_par_cop=FALSE
    parameters=c("theta")
  }
  else if(copula_dist=="N") {
    copula_link=list(fisher_z,fisher_z_inv,dfisher_z_inv); two_par_cop=FALSE
    parameters=c("theta")
  } else if(copula_dist=="t") {
    copula_link=list(fisher_z,fisher_z_inv,dfisher_z_inv,log_2plus,log_2plus_inv,dlog_2plus_inv); two_par_cop=TRUE
    parameters=c("theta","zeta")
  } else {
    stop("ERROR: COPULA DIST LINK FUNCTIONS NOT YET IMPLEMENTED.")
  }

  if(two_par_cop) {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr","zeta.linkfun","zeta.linkinv","zeta.dr")} else {names(copula_link)=c("theta.linkfun","theta.linkinv","theta.dr")}

  return_list=list()
  return_list[["copula_link"]]=copula_link
  return_list[["copula_dist"]]=copula_dist
  return_list[["parameters"]]=parameters

  return(return_list)
}

#' Summarise fitted copula parameters by time
#'
#' @param object A fitted `gamlss.longitudinal` object.
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param stat Character summary statistic for fitted values, one of "mean" or "median".
#'
#' @return A data frame with fitted theta and tau summaries by time.
#' @export
copula_time_summary <- function(object, lags = 1, stat = c("mean", "median")) {
  if (!inherits(object, "gamlss.longitudinal")) {
    stop("'object' must be a fitted 'gamlss.longitudinal' object.")
  }

  stat <- match.arg(stat)
  copula_info <- get_copula_dist(object$copula_dist)
  has_zeta <- "zeta" %in% copula_info$parameters

  fit_data <- .copula_v2_fit_data(object)
  pair_data <- .copula_v2_pair_data(fit_data, lags = lags)

  agg_fun <- if (stat == "median") stats::median else mean

  time_summary <- do.call(rbind, lapply(split(fit_data, fit_data$time), function(x) {
    out <- data.frame(
      time = x$time[1],
      n_obs = nrow(x),
      theta_fit = agg_fun(x$theta_fit, na.rm = TRUE),
      tau_fit = agg_fun(x$tau_fit, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_zeta) {
      out$zeta_fit <- agg_fun(x$zeta_fit, na.rm = TRUE)
    }
    out
  }))

  pair_summary <- do.call(rbind, lapply(split(pair_data, pair_data$time_pair), function(x) {
    out <- data.frame(
      time_pair = x$time_pair[1],
      n_pairs = nrow(x),
      theta_pair = agg_fun(x$theta_pair, na.rm = TRUE),
      tau_pair = agg_fun(x$tau_fit, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if (has_zeta) {
      out$zeta_pair <- agg_fun(x$zeta_pair, na.rm = TRUE)
    }
    out
  }))

  time_summary <- time_summary[order(time_summary$time), , drop = FALSE]

  if (!has_zeta) {
    # Keep the returned data tidy for one-parameter copulas.
    if ("zeta_fit" %in% names(fit_data)) {
      fit_data$zeta_fit <- NULL
    }
    if ("zeta_pair" %in% names(pair_data)) {
      pair_data$zeta_pair <- NULL
    }
  }

  out <- list(
    time_summary = time_summary,
    pair_summary = pair_summary,
    fit_data = fit_data,
    pair_data = pair_data
  )
  class(out) <- "copula_time_summary"
  out
}

#' Plot fitted copula trends by time
#'
#' @param x A `copula_time_summary` object or a fitted `gamlss.longitudinal` object.
#' @param ... Additional arguments (currently unused).
#' @param lags Integer lags to assess, measured in ordered time steps.
#' @param stat Character summary statistic for fitted values, one of "mean" or "median".
#' @param plot Logical; if TRUE, print the plot.
#'
#' @return Invisibly returns a list with the summary data and plot objects.
#' @export
plot.copula_time_summary <- function(x, ..., lags = 1, stat = c("mean", "median"), plot = TRUE) {
  summary_out <- if (inherits(x, "copula_time_summary")) {
    x
  } else {
    copula_time_summary(object = x, lags = lags, stat = stat)
  }
  time_summary <- summary_out$time_summary

  if (nrow(time_summary) == 0) {
    stop("No fitted copula summaries are available for plotting.")
  }

  time_summary$time <- as.factor(time_summary$time)

  p_theta <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = theta_fit, group = 1)) +
    ggplot2::geom_line(color = "#1f4e79", linewidth = 0.8) +
    ggplot2::geom_point(color = "#1f4e79", size = 2.5) +
    ggplot2::labs(
      title = "Fitted Copula Theta by Time",
      x = "Time",
      y = "Theta"
    ) +
    ggplot2::theme_minimal()

  p_tau <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = tau_fit, group = 1)) +
    ggplot2::geom_line(color = "#e41a1c", linewidth = 0.8) +
    ggplot2::geom_point(color = "#e41a1c", size = 2.5) +
    ggplot2::labs(
      title = "Fitted Copula Kendall's Tau by Time",
      x = "Time",
      y = "Tau"
    ) +
    ggplot2::theme_minimal()

  if ("zeta_fit" %in% names(time_summary)) {
    p_zeta <- ggplot2::ggplot(time_summary, ggplot2::aes(x = time, y = zeta_fit, group = 1)) +
      ggplot2::geom_line(color = "#4d4d4d", linewidth = 0.8) +
      ggplot2::geom_point(color = "#4d4d4d", size = 2.5) +
      ggplot2::labs(
        title = "Fitted Copula Zeta by Time",
        x = "Time",
        y = "Zeta"
      ) +
      ggplot2::theme_minimal()
  } else {
    p_zeta <- NULL
  }

  dashboard <- if (is.null(p_zeta)) {
    ggpubr::ggarrange(p_theta, p_tau, ncol = 1, nrow = 2)
  } else {
    ggpubr::ggarrange(p_theta, p_tau, p_zeta, ncol = 1, nrow = 3)
  }

  if (isTRUE(plot)) {
    print(dashboard)
  }

  invisible(list(
    summary = summary_out,
    p_theta = p_theta,
    p_tau = p_tau,
    p_zeta = p_zeta,
    dashboard = dashboard
  ))
}
#' @export
plotDist <- function (dataset,dist,offdiag_scale=c("response","pseudo"),show_cor_stats=TRUE) {

  offdiag_scale <- match.arg(offdiag_scale)

  time_values <- sort(unique(dataset[, "time"]))
  num_margins=length(time_values)

  margin_data=list()
  margin_pseudo=list()
  for (i in seq_len(num_margins)) {
    margin_data[[i]] <- dataset[dataset[,"time"] == time_values[i], c("subject", "response")]

    r <- rank(margin_data[[i]]$response, ties.method = "average", na.last = "keep")
    n_obs <- sum(!is.na(margin_data[[i]]$response))
    u <- r / (n_obs + 1)
    margin_pseudo[[i]] <- data.frame(subject = margin_data[[i]]$subject, u = u)
  }

  ##plot.new()
  #par(mfrow=c(1,num_margins))

  #for (i in 1:num_margins) {histDist(margin_data[[i]],family=dist,xlab=TeX(paste("$Y_",i,"$")),main=paste("Histogram of margin",i,"and fitted",dist))}
  #invisible(readline(prompt="Press [enter] to continue"))

  plots=list()

  z=1
  for (i in seq_len(num_margins)) {
    for (j in seq_len(num_margins)) {
      if(i==j) {
        input_data=data.frame(X1 = margin_data[[i]]$response)
        x_lab <- latex2exp::TeX(paste("$Y_",i,"$"))

        p <- ggplot2::ggplot(input_data, ggplot2::aes(x=X1)) +
          ggplot2::geom_histogram(bins=30, na.rm=TRUE) +
          ggplot2::labs(x = x_lab)
      }
      if(i!=j) {
        if (offdiag_scale == "pseudo") {
          input_data <- merge(
            margin_pseudo[[i]],
            margin_pseudo[[j]],
            by = "subject",
            suffixes = c(".i", ".j"),
            all = FALSE
          )
          input_data <- input_data[complete.cases(input_data$u.i, input_data$u.j), c("u.i", "u.j")]
          names(input_data) <- c("X1", "X2")
          x_lab <- latex2exp::TeX(paste("$U_",i,"$"))
          y_lab <- latex2exp::TeX(paste("$U_",j,"$"))
        } else {
          input_data <- merge(
            margin_data[[i]],
            margin_data[[j]],
            by = "subject",
            suffixes = c(".i", ".j"),
            all = FALSE
          )
          input_data <- input_data[complete.cases(input_data$response.i, input_data$response.j), c("response.i", "response.j")]
          names(input_data) <- c("X1", "X2")
          x_lab <- latex2exp::TeX(paste("$Y_",i,"$"))
          y_lab <- latex2exp::TeX(paste("$Y_",j,"$"))
        }

        p=ggplot2::ggplot(data=input_data,ggplot2::aes(x=X1,y=X2)) +
          ggplot2::geom_point(size=0.4, alpha=0.25, color="black", na.rm=TRUE) +
          ggplot2::geom_density_2d(contour_var="density",bins=10,color="black") +
          ggplot2::labs(x = x_lab, y = y_lab)

        if (show_cor_stats) {
          if (nrow(input_data) >= 3) {
            pearson_r <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "pearson", use = "complete.obs"))
            kendall_tau <- suppressWarnings(cor(input_data$X1, input_data$X2, method = "kendall", use = "complete.obs"))
            stats_lab <- sprintf("Pearson r = %.3f | Kendall tau = %.3f", pearson_r, kendall_tau)
          } else {
            stats_lab <- "Pearson r = NA | Kendall tau = NA"
          }

          p <- p + ggplot2::labs(subtitle = stats_lab)
        }
      }

      plots[[z]]=p
      z=z+1
    }
  }
  ggpubr::ggarrange(plotlist=plots,ncol=num_margins,nrow=num_margins)

}
#' @keywords internal
#' @noRd
create_longitudinal_dataset <- function(response,covariates,labels=NA) {
  num_time_points=ncol(response)
  if(num_time_points <=1) {print('Not enough time points')}

  dataset<-matrix(data=NA,ncol=2+length(covariates),nrow=0)
  subject<-as.factor(seq(1:nrow(response)))

  for (t in 1:ncol(response)) {

    dataset_temp<-cbind(subject,t,response[,t])

    for (i in 1:length(covariates)) {
      if (ncol(covariates[[i]]) == 1 ) {
        covariate_for_time=covariates[[i]]
      } else {
        covariate_for_time=covariates[[i]][,t]
      }
      dataset_temp<-cbind(dataset_temp,covariate_for_time)
    }

    ###Add dataset temp to full table
    dataset <- rbind(dataset,dataset_temp)
  }

  if(!all(is.na(labels))) {
    colnames(dataset) <- labels
  }

  dataset=dataset[order(dataset$time,dataset$subject),] ###NOTE THIS WILL BREAK GLMM
  rownames(dataset)=1:nrow(dataset)

  return(dataset)
}

#' @export
loadDataset <- function(simOption=5,plot_dist=FALSE,n=100,d=3,copula_dist=NA, margin_dist,copula.link=NA,par.copula,par.margin,covariates_input=NA) {

  if (simOption==1) {
    load("Data/rand_mvt.rds")
    head(rand_mvt)

    # Basic data setup
    response = rand_mvt[,4:18]#[,4:18](4+2) ####Currently limiting to just 5 margins for simplicity
    covariates=list()
    covariates[[1]] = as.data.frame(rand_mvt[,19]) #Age 19:33 - changed to age at start to avoid correlation with time
    covariates[[2]] = as.data.frame(rand_mvt[,34:48]) #Time 34:48
    covariates[[3]] = as.data.frame(rand_mvt[,3]) #Gender

    # Setup data as longitudinal file
    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))
  }
  else if (simOption==2) {

    # set up D-vine copula model with mixed pair-copulas
    d <- 3
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(2, 2, 0)
    par <- c(logit_inv(.8), logit_inv(.8), logit_inv(.8))
    par2 <- c(log_2plus_inv(2.1),log_2plus_inv(2.1),log_2plus_inv(2.1))

    # transform to R-vine matrix notation
    RVM <- .copula_dvine(order, family, par, par2)
    contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n*t,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {
      margin[covariates[[3]]==0,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i),sigma=exp(0.3+0.2*i),nu=-0.8,tau=0.05)[covariates[[3]]==0]#Update to i*mu/sigma as needed
      margin[covariates[[3]]==1,i]=qZISICHEL(copsim[,i],mu=exp(0.3+0.2*i+0.1),sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.05)[covariates[[3]]==1]#Update to i*mu/sigma as needed
    }

    response = as.data.frame(margin)

    dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  }
  else if (simOption==3) {

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(copula.link$theta.linkinv(par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(copula.link$theta.linkinv(par.copula[1:(length(par.copula)/2)]), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(copula.link$zeta.linkinv(par.copula[(length(par.copula)/2+1):(length(par.copula))]),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==4) {

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==5) {

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==d-1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==6) {

    t=d
    margin_sim=matrix(0,ncol=d,nrow=n)

    for (i in 1:d) {
      margin_sim[,i]=rnorm(n,1,3)
    }

    W=rnorm(n,0,3)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin_sim_out=matrix(0,ncol=d,nrow=n)
    for (i in 1:d) {
      margin_sim_out[,i]=margin_sim[,i]+W
    }

    response=as.data.frame(margin_sim_out)

  }
  else if (simOption==7) {

    t=d
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))

    #par.copula=c(.3); names(par.copula)=c("theta")
    theta_intercept=unlist(par.copula["theta"])
    theta_out=theta_intercept+matrix(rep(covariates_input$theta.time*1:(d-1),n),ncol=d-1,byrow=TRUE) +
    matrix(rep(as.matrix(covariates_input$theta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1)

    theta_inv=copula_input$copula_link$theta.linkinv(theta_out)

    if(length(par.copula)==d-1){
      par=c((par.copula),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(par.copula[1:(length(par.copula)/2)], rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(par.copula[(length(par.copula)/2+1):(length(par.copula))],rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notatio

    RVM=list()

    for (i in 1:n) {
      RVM[[i]] = .copula_dvine(order, c(rep(copula.family,length(theta_inv[i,])),rep(0,dd-(length(theta_inv[i,])))), par=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))), par2=c(theta_inv[i,],rep(0,dd-(length(theta_inv[i,])))))
    }
    #RVM <- .copula_dvine(order, rep(family[1],nrow(theta_inv)), theta_inv, theta_inv*0)
    #contour(RVM)

    copsim=.copula_rvine_sim(n,RVM)


    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }
  else if (simOption==8) {
    #Multivariate Gamma

    U=margin_sim=matrix(0,ncol=d,nrow=n)

    a=.25;b=1.75;mu=rep(1,d)
    W=rbeta(n,shape1=a,shape2=b)

    for (i in 1:d) {
      U[,i]=rgamma(n,shape=a+b,rate=1/mu[i])
      margin_sim[,i]=U[,i]*W
    }

    #Fake covariates
    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=d,nrow=n))*(1:d)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    response=margin_sim
  } else if (simOption==9) { ########TIME VARIANT MU

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)


    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  } else if (simOption==9) { ########TIME VARIANT MU

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))


    if(length(par.copula)==1){
      par=c(rep(par.copula,d-1),rep(0,dd-(d-1)))
      par2=par*0
    } else {
      par <- c(rep(par.copula["theta"],d-1), rep(0,dd-(d-1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"],d-1),rep(0,dd-(d-1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    #contour(RVM)

    t=d
    copsim=.copula_rvine_sim(n,RVM)

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,1),0)) #Gender

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=par.margin[1],sigma=par.margin[2],nu=par.margin[3],tau=par.margin[4])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }  else if (simOption==10) { ########TIME VARIANT SIGMA AND MU

    print("WARNING: SIMULATION MAPS MARGIN AND COPULA PARAMETERS THROUGH LINK-INVERSE FUNCTIONS AFTER ADDING COVARIATE EFFECTS.")

    t=d

    # Setup covariates from covariates_input

    covariates=list()
    covariates[[1]] = as.data.frame(round(runif(n,0,100),0)) #Age
    covariates[[2]] = t(t(matrix(1,ncol=t,nrow=n))*(1:t)) #Time
    covariates[[3]] = as.data.frame(round(runif(n,0,2),0)) #Gender

    # Build treatment-coded factor effects for simulation inputs.
    # Supported coefficient formats for each *.time / *.gender entry:
    # - scalar: same effect for all non-reference levels
    # - vector length (L-1): explicit non-reference effects
    # - vector length L: full per-level effects
    resolve_factor_effect <- function(levels, coef_input, label) {
      lvl <- as.character(levels)
      n_lvl <- length(lvl)

      if (is.null(coef_input) || any(is.na(coef_input))) {
        return(stats::setNames(rep(0, n_lvl), lvl))
      }

      coef_vec <- as.numeric(coef_input)
      if (length(coef_vec) == 1) {
        out <- c(0, rep(coef_vec, max(0, n_lvl - 1)))
      } else if (length(coef_vec) == (n_lvl - 1)) {
        out <- c(0, coef_vec)
      } else if (length(coef_vec) == n_lvl) {
        out <- coef_vec
      } else {
        stop(
          "simOption 10 factor effect '", label, "' has invalid length ", length(coef_vec),
          ". Expected 1, ", n_lvl - 1, ", or ", n_lvl, " for levels: ",
          paste(lvl, collapse = ", "),
          "."
        )
      }

      stats::setNames(out, lvl)
    }

    make_time_factor_component <- function(coef_input, n_cols, label) {
      levels <- as.character(seq_len(n_cols))
      level_effects <- resolve_factor_effect(levels, coef_input, label)
      matrix(rep(level_effects[levels], n), ncol = n_cols, byrow = TRUE)
    }

    make_gender_factor_component <- function(coef_input, n_cols, label) {
      gender_vals <- as.character(as.vector(covariates[[3]][, 1]))
      levels <- sort(unique(gender_vals))
      level_effects <- resolve_factor_effect(levels, coef_input, label)
      subj_effect <- as.numeric(level_effects[gender_vals])
      matrix(rep(subj_effect, n_cols), ncol = n_cols)
    }

    copula_input=get_copula_dist(copula_dist)
    copula.family=copula_input$copula_dist

    apply_margin_link <- function(par_name, par_value, eta_component) {
      if (par_name %in% names(margin_dist$parameters)) {
        linkfun_name <- paste0(par_name, ".linkfun")
        linkinv_name <- paste0(par_name, ".linkinv")
        par_eta_base <- eval(parse(text=paste0("margin_dist$", linkfun_name)))(par_value)
        par_eta <- par_eta_base + eta_component
        return(eval(parse(text=paste0("margin_dist$", linkinv_name)))(par_eta))
      }

      return(NULL)
    }

    mu_out = NULL
    sigma_out = NULL
    nu_out = NULL
    tau_out = NULL

    if ("mu" %in% names(margin_dist$parameters)) {
      mu_eta = make_time_factor_component(covariates_input$mu.time, d, "mu.time") +
        matrix(rep(as.matrix(covariates_input$mu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$mu.gender, d, "mu.gender")
      mu_out = apply_margin_link("mu", par.margin[1], mu_eta)
    }

    if ("sigma" %in% names(margin_dist$parameters)) {
      sigma_eta = make_time_factor_component(covariates_input$sigma.time, d, "sigma.time") +
        matrix(rep(as.matrix(covariates_input$sigma.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$sigma.gender, d, "sigma.gender")
      sigma_out = apply_margin_link("sigma", par.margin[2], sigma_eta)
    }

    if ("nu" %in% names(margin_dist$parameters)) {
      nu_eta = make_time_factor_component(covariates_input$nu.time, d, "nu.time") +
        matrix(rep(as.matrix(covariates_input$nu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$nu.gender, d, "nu.gender")
      nu_out = apply_margin_link("nu", par.margin[3], nu_eta)
    }

    if ("tau" %in% names(margin_dist$parameters)) {
      tau_eta = make_time_factor_component(covariates_input$tau.time, d, "tau.time") +
        matrix(rep(as.matrix(covariates_input$tau.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$tau.gender, d, "tau.gender")
      tau_out = apply_margin_link("tau", par.margin[4], tau_eta)
    }
    theta_eta_out=par.copula[1]+make_time_factor_component(covariates_input$theta.time, d - 1, "theta.time") +
      matrix(rep(as.matrix(covariates_input$theta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) +
      make_gender_factor_component(covariates_input$theta.gender, d - 1, "theta.gender")
    theta_out = copula_input$copula_link$theta.linkinv(theta_eta_out)

    if ("zeta" %in% copula_input$parameters) {
      zeta_eta_out=par.copula[2]+make_time_factor_component(covariates_input$zeta.time, d - 1, "zeta.time") +
        matrix(rep(as.matrix(covariates_input$zeta.age*((covariates[[1]]-50)/100)^2),d-1),ncol=d-1) +
        make_gender_factor_component(covariates_input$zeta.gender, d - 1, "zeta.gender")
      zeta_out = copula_input$copula_link$zeta.linkinv(zeta_eta_out)
    } else {
      zeta_out = matrix(0, nrow = n, ncol = d - 1)
    }

    if (!is.null(mu_out) && any(!is.finite(mu_out))) {
      stop("simOption 10 generated non-finite mu values after link inverse transformation.")
    }
    if (!is.null(sigma_out) && any(!is.finite(sigma_out))) {
      stop("simOption 10 generated non-finite sigma values after link inverse transformation.")
    }
    if (!is.null(nu_out) && any(!is.finite(nu_out))) {
      stop("simOption 10 generated non-finite nu values after link inverse transformation.")
    }
    if (!is.null(tau_out) && any(!is.finite(tau_out))) {
      stop("simOption 10 generated non-finite tau values after link inverse transformation.")
    }

    if (any(!is.finite(theta_out))) {
      stop("simOption 10 generated non-finite theta values after link inverse transformation.")
    }
    if ("zeta" %in% copula_input$parameters && any(!is.finite(zeta_out))) {
      stop("simOption 10 generated non-finite zeta values after link inverse transformation.")
    }

    # Print parameter ranges for quick simulation diagnostics.
    range_str <- function(label, x) {
      if (is.null(x)) return(paste0(label, ": [NA, NA]"))
      sprintf("%s: [%.2f, %.2f]", label, min(x), max(x))
    }
    margin_range_msg <- paste(
      c(
        range_str("MU", mu_out),
        range_str("SIGMA", sigma_out),
        range_str("NU", nu_out),
        range_str("TAU", tau_out)
      ),
      collapse = " | "
    )
    copula_range_msg <- paste0(
      sprintf("THETA: [%.2f, %.2f]", min(theta_out), max(theta_out)),
      if ("zeta" %in% copula_input$parameters) sprintf(" | ZETA: [%.2f, %.2f]", min(zeta_out), max(zeta_out)) else ""
    )
    print(paste("MARGIN RANGES ->", margin_range_msg, "| COPULA RANGES ->", copula_range_msg))

    #Define margin distribution
    qFUN=paste("q",margin_dist$family[1],sep="")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d*(d-1)/2
    order <- 1:d
    family <- c(rep(copula.family,d-1), rep(0,dd-(d-1)))

    # OK now for each row in theta_out and zeta_out, we need to create a new RVM and simulate from it, then apply the qFUN with the appropriate parameters to get the margin values for that row. This is going to be computationally intensive but should work.

    # row-specific copula simulation from theta_out / zeta_out
    copsim <- matrix(NA_real_, nrow = n, ncol = d)
    for (r in 1:n) {
      par_r  <- c(as.numeric(theta_out[r, ]), rep(0, dd - (d - 1)))
      par2_r <- c(
        if ("zeta" %in% copula_input$parameters) as.numeric(zeta_out[r, ]) else rep(0, d - 1),
        rep(0, dd - (d - 1))
      )
      copsim[r, ] <- as.numeric(.copula_rvine_sim(1, .copula_dvine(order, family, par_r, par2_r)))
    }

    margin=matrix(0,ncol=ncol(copsim),nrow=nrow(copsim))
    for ( i in 1:ncol(copsim)) {

      input_list=list(p=copsim[,i]
                      ,mu=mu_out[,i],sigma=sigma_out[,i],nu=nu_out[,i],tau=tau_out[,i])
      args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      qFunOutput_1=do.call(qFUN,args=(input_list[args]))
      #input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      #args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      #qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[,i]=qFunOutput_1#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      #margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response = as.data.frame(margin)
    }
  }

  dataset<-create_longitudinal_dataset(response,covariates,labels=c("subject","time","response","age","year","gender"))

  if(plot_dist==TRUE) {plotDist(dataset,margin_dist)}

  return(dataset)
}
#' @keywords internal
#' @noRd
bvt_norm_true_SE_B0_Bt <- function(sigma_x,sigma_y,rho,n,d){

  #sigma_x=2
  #sigma_y=2
  #rho=.75

  #((1)/(1-rho^2))
  #(1/(sigma_x^2))
  #(1/(sigma_y^2))
  #(1/(sigma_x))
  #(1/(sigma_y))
  #(rho/(sigma_x*sigma_y))

  hessian=((-1)/(1-(rho^2)))*matrix(c( (1/(sigma_x^2)) +(1/(sigma_y^2)) - 2*(rho/(sigma_x*sigma_y)),
                                       ((rho/(sigma_x*sigma_y))-(1/(sigma_y^2))),
                                       ((rho/(sigma_x*sigma_y))-(1/(sigma_y^2))),
                                       (1/(sigma_y^2))),nrow=2)

  vcov_matrix=-solve(hessian)

  true_SE=(diag(vcov_matrix))
  names(true_SE)=c("B0","Bt")

  return(true_SE)
}

########## ARCHIVE ###########

#Given a parameter vector starting values par = (mu,sigma,nu,tau,theta,zeta), return best fit parameters
#' @export
optim_outer <- function(par,dataset,margin_dist,copula_dist,
                        step_size=0.1,verbose=TRUE,use_dlcopdpar=TRUE) {

  #print("THIS FUNCTION ASSUMES RESPONSE IS ORDERED AS TIME, SUBJECT | PAR INPUT MUST BE NAMED")

  copula_input=get_copula_dist(copula_dist)
  copula_number=copula_input$copula_dist
  copula_link=copula_input$copula_link

  num_margins=length(unique(dataset$time))
  margin_names=unique(dataset$time)
  response=dataset$response

  #Set up parameter vector so names are consistent with the distributions

  if(all(is.null(names(par))|is.na(names(par)))) {stop("ERROR: par vector must be named")}
  margin_par=par[names(par)%in%c("mu","sigma","nu","tau")]
  copula_par=par[!names(par)%in%c("mu","sigma","nu","tau")]

  ##### Calculate all relevant derivatives / CG method with first and second derivatives

  ### Calculate margin derivatives w.r.t. margin parameters

  #Get names for margin derivatives from margin_dist
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)

  #Get link transforms (eta) and derivatives w.r.t to link for parameters
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                margin_names = sort(unique(dataset$time))
                theta_rows = which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names)-1))])
                theta_index_map=rep(NA_integer_,length(dataset$response))
                theta_index_map[theta_rows]=seq_along(theta_rows)
                par_idx <- theta_index_map[row_id1]
              }

              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              d1_full[par_idx[valid_idx],1] <- dldth[valid_idx]
    if(par_name %in% names(margin_par)) {
      par_eta[par_name]=margin_dist[[paste(par_name,".linkfun",sep="")]](par[par_name])
      par_eta_dr[par_name]=margin_dist[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
    if(par_name %in% names(copula_par)) {
            n_par <- length(eta[[par_name]])
            d1_full=matrix(0,nrow=n_par,ncol=1)
      par_eta_dr[par_name]=copula_link[[paste(par_name,".dr",sep="")]](par_eta[par_name])
    }
              if(n_par == length(dataset$response)) {
                par_idx <- row_id1
              } else {
                margin_names = sort(unique(dataset$time))
                theta_rows = which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names)-1))])
                theta_index_map=rep(NA_integer_,length(dataset$response))
                theta_index_map[theta_rows]=seq_along(theta_rows)
                par_idx <- theta_index_map[row_id1]
              }

              valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
              d1_full[par_idx[valid_idx],1] <- dldz[valid_idx]
  #Setup input matrix of response and parameters
  margin_deriv_input=list()
  margin_deriv_input[["y"]]=response
  margin_deriv_input[["q"]]=response
  margin_deriv_input[["x"]]=response
  for (par_name in names(margin_par)) {
    margin_deriv_input[[par_name]]=rep(margin_par[par_name],length(response))
  }

  #Calculate all derivatives
  margin_deriv=list()
  for (deriv_name in margin_deriv_names) {
    FUN=margin_dist[[deriv_name]]
    FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
    margin_deriv[[deriv_name]]=do.call(FUN,args=margin_deriv_input[FUN_args])
    margin_deriv[[deriv_name]][!is.finite(margin_deriv[[deriv_name]])]=0
  }

  margin_pFUN=eval(parse( text=paste("p",margin_dist$family[1],sep="") ))
  FUN=margin_pFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_p=do.call(FUN,args=margin_deriv_input[FUN_args])

  margin_dFUN=eval(parse( text=paste("d",margin_dist$family[1],sep="") ))
  FUN=margin_dFUN
  FUN_args=names(margin_deriv_input)[names(margin_deriv_input)%in%formalArgs(FUN)]
  margin_d=do.call(FUN,args=margin_deriv_input[FUN_args])

  ### Calculate copula derivatives w.r.t. copula parameters

  #First calculate margin F(x1), F(x2) as inputs to copula

  Fx_1_2=matrix(ncol=2,nrow=0)
  order_copula=data.frame()
  for (i in 1:(num_margins-1)) {
    Fx_1_2=rbind(Fx_1_2,cbind(margin_p[dataset$time == margin_names[i]],margin_p[dataset$time == margin_names[i+1]]))
    order_copula=rbind(order_copula,cbind(dataset[dataset$time == margin_names[i],c("time","subject")],dataset[dataset$time == margin_names[i+1],c("time","subject")]))
  }
  names(order_copula)=c("time1","subject1","time2","subject2")

  par1=copula_par["theta"]
  if(is.na(copula_par["zeta"])) {par2=0} else {par2=copula_par["zeta"]}

  #Handling extreme values
  Fx_1_2[Fx_1_2>1]=1;Fx_1_2[Fx_1_2<0]=0

  if(copula_number==3) {
    if(par1>28){par1=28}
  }

  copula_d=.copula_pdf(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2)
  dldth=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=TRUE)
  dcdth=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par",log=FALSE)
  d2cdth=.copula_deriv2( Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par")
  d2ldth2=(1/(copula_d^2))*(copula_d*d2cdth-dcdth^2)
  if(!is.na(copula_par["zeta"])) {
    dldz=.copula_deriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=TRUE)
    dcdz=.copula_deriv(    Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2",log=FALSE)
    d2cdz=.copula_deriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par2")
    d2ldz2=(1/(copula_d^2))*(copula_d*d2cdz-dcdz^2)

    d2cdthdz=.copula_deriv2(  Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="par1par2")
    d2ldthdz=(d2cdthdz*copula_d-dcdth*dcdz)/(copula_d^2)
  }
  dcdu1=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  dcdu2=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

  d2cdu12=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u1",log=FALSE)
  d2cdu22=.copula_deriv(   Fx_1_2[,1],Fx_1_2[,2],family = copula_number,par=par1,par2=par2,deriv="u2",log=FALSE)

  d2ldth2[!is.finite(d2ldth2)]=0

  ### Calculate copula derivatives w.r.t margin parameters

  #Extract margin calculations for F(x), f(x), response and derivatives at time 1 and time 2, join to copula values for time 1 and time 2
  margin_deriv_1=margin_deriv_2=margin_deriv_2cross=matrix(ncol=length(margin_par),nrow=length(response))
  for (i in 1:length(margin_par)) {
    margin_deriv_1[,i]=margin_deriv[grepl("dld",names(margin_deriv))][[i]]
    margin_deriv_2[,i]=margin_deriv[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")][[i]]
  }
  colnames(margin_deriv_1)=paste("dld",names(margin_par),sep="")
  colnames(margin_deriv_2)=paste(paste("d2ld",names(margin_par),sep=""),"2",sep="")

  #colnames(margin_deriv_2)=paste("d2ld",names(margin_par),sep="")

  order_margin=dataset[,c("time","subject")]
  margin_components=cbind(order_margin,response,margin_p,margin_d,margin_deriv_1,margin_deriv_2)
  margin_components_Ft_plus=margin_components
  margin_components_Ft_plus$time=normalize_lag_time(margin_components_Ft_plus$time)
  margin_plus=merge(margin_components,margin_components_Ft_plus,by=c("time","subject"),all.x=TRUE)

  copula_components=cbind(order_copula,dcdu1,dcdu2,copula_d,d2cdu12,d2cdu22)
  copula_merged=merge(copula_components,margin_plus,by.x=c("time1","subject1"),by.y=c("time","subject"),all.x=TRUE)

  #Calculate copula derivative with respect to marginal parameters
  input=copula_merged
  dlcopdpar=matrix(0,nrow=nrow(input),ncol=length(margin_par))
  d2lcopdpar2=matrix(0,nrow=nrow(input),ncol=length(margin_par))

  i=1
  for (par_name in names(margin_par)) {

    #Take parameters from input for clarity
    dc_tplus_du_t=input[,"dcdu1"]
    dc_tplus_du_tplus=input[,"dcdu2"]
    l_t=input[,paste(paste("dld",par_name,sep=""),".x",sep="")]
    l_t_plus=input[,paste(paste("dld",par_name,sep=""),".y",sep="")]
    x_t=input[,"response.x"]
    x_t_plus=input[,"response.y"]
    f_t=input[,"margin_d.x"]
    f_t_plus=input[,"margin_d.y"]
    du_t_dmu=x_t*f_t*l_t
    du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus
    c_tplus=input[,"copula_d"]

    du_t_dmu=x_t*f_t*l_t
    du_t_plus_dmu=x_t_plus*f_t_plus*l_t_plus

    dc_plus_dt_dmu=dc_tplus_du_t * du_t_dmu
    dc_plus_dt_plus_dmu=dc_tplus_du_tplus * du_t_plus_dmu
    dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)]=0
    dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)]=0
    dcdmu_tplus=((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
    dcdmu_tplus[is.nan(dcdmu_tplus)|is.na(dcdmu_tplus)]=0

    dlcopdpar[,i]=dcdmu_tplus


    #######NOW FOR SECOND DERIVATIVE OF COPULA TERM

    l2_t=input[,paste(paste(paste("d2ld",par_name,sep=""),"2",sep=""),".x",sep="")]
    l2_tplus=input[,paste(paste(paste("d2ld",par_name,sep=""),"2",sep=""),".y",sep="")]

    df_t_dmu=f_t*l_t
    df_t_plus_dmu=f_t_plus*l_t_plus

    d2f_t_dmu=df_t_dmu*l_t + f_t*l2_t
    d2f_t_plus_dmu=df_t_plus_dmu*l_t_plus + f_t_plus*l2_tplus

    d2u_t_dmu2=x_t*d2f_t_dmu
    d2u_t_plus_dmu2=x_t_plus*d2f_t_plus_dmu

    d2cdu_t2=input[,"d2cdu12"]
    d2cdu_t_plus2=input[,"d2cdu22"]
    d2cdu_t2[is.nan(d2cdu_t2)]=0
    d2cdu_t_plus2[is.nan(d2cdu_t_plus2)]=0

    d2cdmu2=d2cdu_t2*du_t_dmu^2 + dc_tplus_du_t * d2u_t_dmu2 + d2cdu_t_plus2*du_t_plus_dmu^2 + dc_tplus_du_tplus * d2u_t_plus_dmu2

    d2lcdmu2=as.matrix((d2cdmu2*c_tplus-(dcdmu_tplus^2))/(c_tplus^2))

    d2lcopdpar2[,i]=d2lcdmu2
    #num_deriv=margin_copula_merged_2[,"num_dlcopdpar_ordered.Ft"]
    #num_deriv_nolog=margin_copula_merged_2[,"num_dlcopdpar_nolog_ordered.Ft"]

    i=i+1
  }
  colnames(dlcopdpar)=paste("dlcopd",names(margin_par),sep="")
  colnames(d2lcopdpar2)=paste(paste("d2lcd",names(margin_par),sep=""),"2",sep="")

  dlcopdpar[!is.finite(dlcopdpar)]=0
  d2lcopdpar2[!is.finite(d2lcopdpar2)]=0

  #### Define score and hessian

  score=par*0
  hessian=matrix(0,nrow=length(par),ncol=length(par))
  colnames(hessian)=names(par);rownames(hessian)=names(par)
  names(score)=names(par)

  margin_deriv_sum=vector()
  for (i in 1:length(margin_deriv)) {
    margin_deriv[[i]][!is.finite(margin_deriv[[i]])]=0
    margin_deriv_sum[i]=sum(margin_deriv[[i]])
  }
  names(margin_deriv_sum)=names(margin_deriv)

  margin_d1=margin_deriv_sum[grepl("dld",names(margin_deriv))]
  margin_d2=margin_deriv_sum[grepl("d2ld",names(margin_deriv))&endsWith(names(margin_deriv),"2")]
  margin_d2d=margin_deriv_sum[grepl("d2ld",names(margin_deriv))&!endsWith(names(margin_deriv),"2")]

  if(is.na(copula_par["zeta"])) {
    copula_d1=sum(dldth)
    copula_d2=sum(d2ldth2)
  } else {
    copula_d1=colSums(cbind(dldth,dldz))
    copula_d2=colSums(cbind(d2ldth2,d2ldz2))
  }
  margin_d1_dlcopdpar=margin_d1+if(use_dlcopdpar==TRUE){ colSums(dlcopdpar)} else {colSums(dlcopdpar)*0}
  margin_d2_dlcopdpar=margin_d2+if(use_dlcopdpar==TRUE){ colSums(d2lcopdpar2)*0} else {colSums(d2lcopdpar2)*0}
  score=c(margin_d1_dlcopdpar,copula_d1)

  ###CALCULATING HESSIAN USING D2
  diag(hessian)=c(margin_d2_dlcopdpar,copula_d2)
  hessian[1:length(margin_par),1:length(margin_par)][upper.tri(hessian[1:length(margin_par),1:length(margin_par)])]=margin_d2d
  hessian[1:length(margin_par),1:length(margin_par)][lower.tri(hessian[1:length(margin_par),1:length(margin_par)])]=margin_d2d

  #Why isn't d2 for copula negative?
  copula_hess=hessian[(length(margin_par)+1):(length(margin_par)+length(copula_par)),(length(margin_par)+1):(length(margin_par)+length(copula_par))]
  if(!is.na(copula_par["zeta"])) {
    copula_hess[upper.tri(copula_hess)]=sum(d2ldthdz)
    copula_hess[lower.tri(copula_hess)]=sum(d2ldthdz)
  }
  hessian[(length(margin_par)+1):(length(margin_par)+length(copula_par)),(length(margin_par)+1):(length(margin_par)+length(copula_par))]=copula_hess


  ###STILL NEED TO CALCULATE d2 for marginal parameters with respect to copula likelihood and add to hessian values

  #par_end=par-(solve(-hessian)%*%(score))

  #weights_eta=diag((1/(score_eta^2)))
  #weights=-diag(score*score)

  #score=score
  #weights=-solve(hessian)

  weights_eta=-solve(hessian*par_eta_dr*par_eta_dr)

  #weights_eta=diag(1/(score*score*par_eta_dr*par_eta_dr))
  score_eta=score*par_eta_dr
  #par_end=par*(1-step_size) + step_size*(par+par_change)

  par_end=par*0
  names(par_end)=names(par_eta_end)=names(par)
  #Get end paraemters re-transformed
  #for (par_name in names(par)) {
  #  if(par_name %in% names(margin_par)) {
  #    par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  #  if(par_name %in% names(copula_par)) {
  #    par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  #}
  ###If calculating for eta
  for (par_name in names(par)) {
    if(par_name %in% names(margin_par)) {
      par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
    }
    if(par_name %in% names(copula_par)) {
      par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
    }
  }


  sum_log_margin_p=sum(log(margin_d)[is.finite(log(margin_d))])
  sum_log_copula_d=sum(log(copula_d)[is.finite(log(copula_d))])

  log_lik=c(sum_log_copula_d,sum_log_margin_p,sum_log_copula_d+sum_log_margin_p)
  names(log_lik)=c("copula","margin","joint")

  if(verbose==TRUE) {
    print("Start Parameters")
    print(par)
    print("End Parameters:")
    print(par_end)
    print("Score:")
    print(score)
    print("Hessian:")
    print(hessian)
    print("Weights:")
    print(weights_eta)

    print(log_lik)
  }

  return(list(score=score,hessian=hessian,par_end=par_end,par_eta_end=par_eta_end,par_start=par,log_lik=log_lik))
}

# Load analytical Hessian helpers when sourcing this file directly (development workflow).
local({
  candidates <- c("R/analytical_hessian.R",
                  file.path(getwd(), "R", "analytical_hessian.R"))
  for (p in candidates) {
    if (file.exists(p)) { source(p, local = FALSE); break }
  }
})
