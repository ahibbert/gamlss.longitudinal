#' Fit a longitudinal joint regression model

#'

#' This function fits a longitudinal model to a dataset with gamlss margins

#' and copula fit to dependence. Any linear or factor covariates can be fit

#' to any parameters of the copula or margin distributions. The model is fit

#' using RS() optimisation and the joint likelihood by default. Select

#' use dlcopdpar=FALSE to fit separately optimised models for the margin and

#' copula likelihoods which can be quicker with a slight loss to overall fit.

#'

#' Formula inputs are converted to fixed-effect and smooth model matrices by

#' `create_model_matrices()`. The user-supplied `time_var` is preserved for

#' formulas as an internal `time_covariate` column, while an internal numeric

#' `time` index is used for ordering longitudinal margins and adjacent copula

#' pairs. Numeric and integer time inputs are used directly, numeric-like

#' character time inputs are converted with [as.numeric()], and categorical time

#' should be supplied as a factor. Ordered factors are treated with treatment

#' contrasts rather than polynomial contrasts so factor-level effects are

#' interpretable as level comparisons.

#'

#' The fitting data are converted to a plain data frame. Structurally missing

#' subject-time combinations are expanded to explicit rows with missing

#' responses. Response `NA` values may be present unless an entire margin or an

#' adjacent copula pair has no complete responses, in which case fitting stops.

#' Submitted predictor values must be observed and finite. Response `NA` values

#' represent missing outcomes, but response `NaN`, `Inf`, and `-Inf` values stop

#' before fitting. The package does not impute missing values; structurally

#' inserted rows are represented in `model.frame(fit, type = "expanded")`, while

#' observed response rows are available from `model.frame(fit, type = "observed")`.

#'

#' @param dataset Long-format data frame containing the response, subject, time,

#'   and covariate columns.

#' @param margin_dist Marginal distribution specified as a gamlss family object,

#' e.g. GA(), NO(), PO(), NBI(), etc.

#' @param copula_dist Copula distribution code, one of "N", "C", "F", "G", "J", or "t".

#' @param time_var Name of the time variable in `dataset`.

#' @param subject_var Name of the subject identifier in `dataset`.

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

#' @param discrete_score_method Character. For discrete margins using exact

#' rectangle likelihoods, choose `"analytical"` for vectorised rectangle-score

#' assembly or `"finite"` for slow row-wise finite-difference scores.

#' @param cg_zeta_hessian Character; for `method = "CG"` only. `"analytical"`

#' uses the analytical Hessian for the zeta block, while `"finite"` replaces

#' the zeta-zeta block with central finite differences of the raw joint

#' log-likelihood.

#' @param cg_hessian_method Character; for `method = "CG"` only. `"analytical"`

#' uses the semi-analytical Hessian for Newton steps, `"finite"` uses a full

#' finite-difference Hessian, and `"auto"` tries analytical then falls back to

#' finite differences when needed.

#' @param compute_vcov Logical; if `TRUE` (default), compute and store the

#' model variance-covariance output at the end of fitting.

#' @param vcov_method Character; fit-time vcov method when `compute_vcov = TRUE`.

#' One of `"analytical"` or `"numderiv"`. Analytical vcov falls back to the

#' numerical reference path if the analytical Hessian cannot be inverted.

#' @param vcov_numderiv Logical; passed to `vcov.gamlss.longitudinal()` when

#' `compute_vcov = TRUE`.

#' @param use_Rcpp Use Rcpp for matrix operations

#' @param lambda_start Optional starting value for smooth-term penalties.

#' @param lambda_penalty_K Penalty strength used when updating smooth-term

#'   smoothing parameters.

#' @param rs_update_lambda Logical; for `method = "RS"` only. If `TRUE`,

#' update smoothing parameters by the RS GAIC step; if `FALSE`, keep

#' `lambda_start` fixed.

#' @param rs_smooth_trust_radius Numeric; for `method = "RS"` only. Optional

#' L2 trust radius applied separately to each smooth coefficient block after

#' the RS weighted least-squares proposal. Use `Inf` to disable.

#'

#' @export

gamlss_longitudinal=function(dataset,

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

                        discrete_score_method=c("analytical","finite"),

                        cg_zeta_hessian="analytical",

                        cg_hessian_method=c("analytical","finite","auto"),

                        compute_vcov=TRUE,

                        vcov_method=c("analytical","numderiv"),

                        vcov_numderiv=FALSE,

                        use_Rcpp=FALSE,

                        lambda_start=NA,

                        lambda_penalty_K=2,

                        rs_update_lambda=TRUE,

                        rs_smooth_trust_radius=Inf

                      )

{

  fit_start_time <- Sys.time()

  margin_dist <- .normalise_margin_dist_links(margin_dist)

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

  if (length(rs_smooth_trust_radius) != 1 || !is.numeric(rs_smooth_trust_radius) ||

      is.na(rs_smooth_trust_radius) || rs_smooth_trust_radius <= 0) {

    stop("ERROR: rs_smooth_trust_radius must be a single positive numeric value or Inf.")

  }


  method <- toupper(as.character(method)[1])

  cg_line_search <- match.arg(as.character(cg_line_search)[1], c("first", "best"))

  cg_gradient_method <- match.arg(as.character(cg_gradient_method)[1], c("analytical", "forward", "central"))

  discrete_score_method <- match.arg(as.character(discrete_score_method)[1], c("analytical", "finite"))

  cg_zeta_hessian <- match.arg(as.character(cg_zeta_hessian)[1], c("analytical", "finite"))

  cg_hessian_method <- match.arg(as.character(cg_hessian_method)[1], c("analytical", "finite", "auto"))

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

  .gl_validate_tabular_shape(dataset, context = "dataset")


  # Validate and prepare input data

  if(all(is.na(time_var)) || all(is.na(subject_var))) {

    stop("ERROR: Required input variables not specified.\n",

         "Please specify:\n",

         "  - time_var: column name for time/margin variable (e.g., 'time')\n",

         "  - subject_var: column name for subject ID variable (e.g., 'subject')\n",

         "Example: gamlss_longitudinal(..., time_var='time', subject_var='subject')")

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

    warning(

      "Converted character time variable '", time_var,

      "' to numeric for fitting; convert it to factor before fitting if visits are categorical.",

      call. = FALSE

    )

    dataset$time <- time_numeric

    dataset$time_covariate <- time_numeric

  } else {

    stop("ERROR: Unsupported time variable type: ", class(dataset$time_covariate)[1],

         ". Use numeric/integer, numeric-like character, or factor.")

  }


  if (is.factor(dataset$subject)) {

    dataset$subject <- as.character(dataset$subject)

  }


  if (any(is.na(dataset$time)) || any(is.na(dataset$subject)) || any(!is.finite(dataset$time))) {

    stop("ERROR: time and subject variables cannot contain missing or non-finite values.")

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


  .gl_validate_fitting_data_policy(

    dataset,

    formulas = list(mu.formula.int, sigma.formula.int, nu.formula.int, tau.formula.int, theta.formula.int, zeta.formula.int),

    response_name = "response"

  )


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


  .warn_rank_deficient_model_matrices <- function(mm_x) {

    for (parameter in names(mm_x)) {

      X <- mm_x[[parameter]]

      if (!is.matrix(X) && !is.data.frame(X)) next

      X <- as.matrix(X)

      if (nrow(X) == 0L || ncol(X) < 2L) next

      finite_rows <- stats::complete.cases(X)

      if (sum(finite_rows) < 2L) next

      qr_x <- qr(X[finite_rows, , drop = FALSE])

      if (qr_x$rank < ncol(X)) {

        warning(

          "Fixed-effect model matrix for parameter '", parameter,

          "' is rank deficient; estimates may be non-identifiable.",

          call. = FALSE

        )

      }

    }

    invisible(NULL)

  }

  .warn_rank_deficient_model_matrices(mm$x)


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

    warm_warnings <- character(0)

    warm_err <- NULL

    tryCatch({

      warm_output <- capture.output(

        withCallingHandlers(

          {

            warm_fit <- gamlss_longitudinal(

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

              discrete_score_method = discrete_score_method,

              cg_zeta_hessian = cg_zeta_hessian,

              cg_hessian_method = cg_hessian_method,

              compute_vcov = FALSE,

              vcov_method = vcov_method,

              vcov_numderiv = vcov_numderiv,

              use_Rcpp = use_Rcpp,

              lambda_start = lambda_start,

              lambda_penalty_K = lambda_penalty_K

            )

          },

          warning = function(w) {

            warm_warnings <<- c(warm_warnings, conditionMessage(w))

            invokeRestart("muffleWarning")

          }

        ),

        type = "output"

      )

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

      captured_output = warm_output,

      captured_warnings = warm_warnings

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

      hess_paths <- file.path(getwd(), "R", c(
        "hessian-setup.R",
        "hessian-margin-cdf.R",
        "hessian-copula.R",
        "hessian-assembly.R",
        "hessian-analytical.R"
      ))

      if(all(file.exists(hess_paths))) {

        for(hess_path in hess_paths) source(hess_path, local = FALSE)

      }

    }

    if(!exists("calc_analytical_hessian", mode = "function")) {

      stop("CG requires calc_analytical_hessian(); source the R/hessian-*.R files first.")

    }

    if(verbose > 0) {

      cat("\nUsing optimization method: CG")

      cat(paste0(

        "\nCG controls: max_delta=", signif(cg_max_delta, 4),

        " | lambda_update_every=", cg_lambda_update_every,

        " | update_lambda=", isTRUE(cg_update_lambda),

        " | line_search=", cg_line_search,

        " | gradient=", cg_gradient_method,

        " | hessian=", cg_hessian_method,

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

      if(copula_dist %in% c("N", "t") &&

         "theta" %in% names(eta_inv) &&

         any(abs(eta_inv$theta) >= 0.999, na.rm = TRUE)) return(NULL)

      positive_names <- intersect(names(eta_inv), c("sigma", "tau", "zeta"))

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


    cg_hessian_ok <- function(H, beta_names) {

      is.matrix(H) &&

        identical(dim(H), c(length(beta_names), length(beta_names))) &&

        all(is.finite(H))

    }


    cg_observed_hessian <- function(tmp_obj, beta_vec, mm_cg, context = "CG iteration") {

      beta_names <- names(beta_vec)

      use_finite <- identical(cg_hessian_method, "finite")

      H <- NULL

      analytical_error <- NULL


      if(!use_finite) {

        H <- tryCatch(

          calc_analytical_hessian(tmp_obj, progress = FALSE),

          error = function(e) {

            analytical_error <<- conditionMessage(e)

            NULL

          }

        )

        if(cg_hessian_ok(H, beta_names)) {

          return(0.5 * (H + t(H)))

        }

        if(verbose > 0) {

          msg <- if(!is.null(analytical_error)) analytical_error else "non-finite analytical Hessian"

          warning(

            "CG analytical Hessian failed during ", context,

            "; falling back to finite-difference Hessian. Reason: ", msg,

            call. = FALSE

          )

        }

      }


      H_fd <- cg_finite_hessian_block(beta_vec, beta_names, mm_cg)

      if(!cg_hessian_ok(H_fd, beta_names)) {

        msg <- if(!is.null(analytical_error)) analytical_error else "finite-difference Hessian was non-finite"

        stop("CG failed to construct a usable Hessian during ", context, ". Reason: ", msg, call. = FALSE)

      }

      0.5 * (H_fd + t(H_fd))

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

      if(verbose > 0) cat(paste("\nOUTER ITERATION:", outer_only_run_counter))

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

      H_obs <- cg_observed_hessian(tmp_obj, beta_all, mm_cg, context = paste0("outer iteration ", outer_only_run_counter))

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

      if(verbose > 0) {

        cat("\n")

        print(out_temp)

      }


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

          if(identical(cg_stop_reason, "tolerance") && verbose > 0) cat("\nOUTER CONVERGED")

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

    final_H <- tryCatch(

      cg_observed_hessian(final_obj, beta_all, mm_cg, context = "final smooth EDF update"),

      error = function(e) NULL

    )

    if(!is.null(final_H)) {

      penalty_mat <- build_cg_penalty(names(beta_all), lambda_s)

      df_s <- cg_smooth_edf_list(final_H, penalty_mat, names(beta_all))

    }


    for(pn in names(mm$x)) weights_final[[pn]] <- rep(1, nrow(mm$x[[pn]]))

  } else {

  while ((first_outer_run==TRUE | (abs(outer_log_lik_change)>outer_stop_crit)) & outer_only_run_counter < max_outer_iter) {

    check_elapsed_budget("RS outer iteration")


    if(verbose > 0) cat(paste("\nOUTER ITERATION:",outer_only_run_counter))

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


        discrete_scores <- NULL

        if (identical(calc_lik_out$likelihood_type, "discrete_rectangle")) {

          discrete_scores <- .calc_discrete_rectangle_scores(

            eta_inv,

            mm$x,

            margin_dist,

            copula_dist,

            dataset$response,

            dataset$time,

            dataset$subject,

            pair_cache = pair_cache,

            calc_lik = calc_lik_out,

            method = discrete_score_method

          )

        }


        ### Calculate copula derivatives w.r.t margin parameters

        if (

          !is.null(discrete_scores) &&

          (!par_name %in% c("mu", "sigma", "nu", "tau") || isTRUE(include_dlcopdpar))

        ) {

          d1 <- as.matrix(discrete_scores[[par_name]])

          colnames(d1) <- paste0("dld", par_name)

        } else if(!par_name %in% c("mu","sigma","nu","tau")) {

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

          if (is.finite(rs_smooth_trust_radius) && length(par_s[[par_name]]) > 0) {

            for (s_name in names(smooth_penalty_meta)) {

              idx <- smooth_penalty_meta[[s_name]]$idx

              delta_s <- beta_new[idx] - beta_start[idx]

              delta_norm <- sqrt(sum(delta_s^2))

              if (is.finite(delta_norm) && delta_norm > rs_smooth_trust_radius) {

                beta_new[idx] <- beta_start[idx] + delta_s * (rs_smooth_trust_radius / delta_norm)

              }

            }

          }


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

              ,pair_cache=pair_cache,margin_eval_cache=margin_eval_cache,calc_margin_deriv=FALSE)

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

           if(isTRUE(rs_update_lambda) && inner_run_counter==1) {

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

    if(verbose > 0) {

      cat("\n")

      print(out_temp)

    }


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

      if(verbose > 0) {

        print(c(outer_end_log_lik-outer_start_log_lik))

        cat("\nOUTER CONVERGED")

      }

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

    cg_zeta_hessian = if(identical(method, "CG")) cg_zeta_hessian else NA_character_,

    cg_hessian_method = if(identical(method, "CG")) cg_hessian_method else NA_character_

  )


  if (isTRUE(hit_outer_limit)) {

    warning(

      "Model stopped at max_outer_iter before satisfying outer_stop_crit; treat fit as not converged.",

      call. = FALSE

    )

  }


  total_fit_time <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))


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

  if(verbose > 0) {

    cat("\n\n############ MODEL FIT ############\n")

    cat(paste("\nMargin distribution:",margin_dist$family[2]))

    cat(paste("\nCopula distribution:",copula_dist))

    cat("\n")

    cat(paste("\nParameter count:",length(par_cov)))

    cat(paste("\nObservations:",nrow(dataset)))

    cat(paste("\nMargins:",length(unique(dataset$time))))

    cat("\n")

    cat(paste("\nTotal time (seconds):",round(total_fit_time,2)))

    cat("\n\n")

    par_mat_out_temp=t(t((par_cov)))

    colnames(par_mat_out_temp) = c("estimate")

    print(par_mat_out_temp)

    cat("\n")

    cat("Model Selection Criteria:")

    cat("\n")

    print(aics)

    cat("\n####################################\n")

  }


  return_list=list(par_cov,log_lik_history,par_history,calc_lik_out_end,mm,margin_dist,copula_dist,include_dlcopdpar,dataset$response,dataset$time,dataset$subject,par_s,lambda_s,df_s,weights_final)

  names(return_list)=c("par","log_lik_history","par_history","calc_lik_out_end","model_matrix","margin_dist","copula_dist","include_dlcopdpar","response","response_margin","response_subject","par_s","lambda_s","df_s","weights")

  return_list$dataset <- dataset

  return_list$dataset_original <- dataset_original

  return_list$response_var <- response_var

  return_list$time_var <- time_var

  return_list$subject_var <- subject_var

  return_list$formulas <- list(

    mu = mu.formula,

    sigma = sigma.formula,

    nu = nu.formula,

    tau = tau.formula,

    theta = theta.formula,

    zeta = zeta.formula

  )

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

      if(!is.null(vcov_cached$method)) {

        return_list$vcov_meta$method_used <- vcov_cached$method

      }

      if(!is.null(vcov_cached$method_requested)) {

        return_list$vcov_meta$method <- vcov_cached$method_requested

      }

    }

  }


  class(return_list)="gamlss.longitudinal"

  return(return_list)

}


gamlss.longitudinal <- gamlss_longitudinal


