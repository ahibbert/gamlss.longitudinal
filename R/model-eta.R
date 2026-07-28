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
#' \item{eta_inv}{A list of inverse link function transforms of the linear predictor values for each parameter.}
#' \item{eta_dr}{A list of derivatives of the link function with respect to the linear predictor values for each parameter.}
#'
#' @keywords internal
#' @noRd
calc_eta <- function(par_cov, mm, margin_dist, copula_link, par_s = NA) {
  eta <- list()

  for (par_name in names(mm$x)) {
    par_cov_single <- par_cov[grepl(par_name, names(par_cov))]
    mm_temp <- mm$x[[par_name]]

    # If there are no smooth terms for the parameter then just do standard fixed calculation, otherwise add in the smooth term contributions
    if (all(is.na(par_s[[par_name]]))) {
      eta[[par_name]] <- rowSums(mm_temp * matrix(rep(par_cov_single, each = nrow(mm_temp)), ncol = length(par_cov_single), dimnames = list(NULL, c(names(par_cov_single)))))
    } else {
      eta[[par_name]] <-
        rowSums(mm_temp * matrix(rep(par_cov_single, each = nrow(mm_temp)), ncol = length(par_cov_single), dimnames = list(NULL, c(names(par_cov_single)))))
      for (s_name in names(mm$s[[par_name]])) {
        eta[[par_name]] <- eta[[par_name]] + mm$s[[par_name]][[s_name]] %*% par_s[[par_name]][[s_name]]
      }
    }
  }

  # Get link transforms (eta) and derivatives w.r.t to link for parameters
  eta_dr <- eta_inv <- list()
  for (par_name in names(mm$x)) {
    if (par_name %in% c("mu", "sigma", "nu", "tau")) {
      eta_inv[[par_name]] <- margin_dist[[paste(par_name, ".linkinv", sep = "")]](eta[[par_name]])
      eta_dr[[par_name]] <- margin_dist[[paste(par_name, ".dr", sep = "")]](eta[[par_name]])
    }
    if (par_name %in% c("theta", "zeta")) {
      eta_inv[[par_name]] <- copula_link[[paste(par_name, ".linkinv", sep = "")]](eta[[par_name]])
      eta_dr[[par_name]] <- copula_link[[paste(par_name, ".dr", sep = "")]](eta[[par_name]])
    }
  }
  return(list(eta = eta, eta_inv = eta_inv, eta_dr = eta_dr))
}

.calc_eta_rs_cached <- function(
    rs_design_cache,
    par_cov,
    par_s,
    margin_dist,
    copula_link,
    update_only = NULL,
    eta_out = NULL) {
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
