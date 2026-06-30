#' Validate and coerce RS score inputs
#'
#' @noRd
.gl_prepare_rs_score_inputs <- function(eta, eta_dr, d1, par_name) {
  eta_len <- length(eta[[par_name]])
  d1_vec <- as.numeric(d1)
  eta_dr_vec <- as.numeric(eta_dr[[par_name]])

  if (length(d1_vec) != eta_len) {
    stop(
      "Score derivative length mismatch for ", par_name,
      ": length(d1)=", length(d1_vec),
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

  list(
    eta = eta[[par_name]],
    d1 = d1_vec,
    eta_dr = eta_dr_vec
  )
}

#' Build the starting coefficient vector for one RS parameter block
#'
#' @noRd
.gl_rs_beta_start <- function(par_name, par_cov, par_s, design_info) {
  fixed_names <- design_info$fixed_names
  X <- design_info$X

  if (length(par_s[[par_name]]) == 0) {
    paste("No smooths found for parameter; running basic IRLS", par_name)
    return(c(par_cov[fixed_names]))
  }

  temp_par_s_unlisted <- unlist(par_s[[par_name]], use.names = FALSE)
  names(temp_par_s_unlisted) <- setdiff(colnames(X), fixed_names)
  c(par_cov[fixed_names], temp_par_s_unlisted)
}

#' Prepare RS backfitting inputs for one parameter block
#'
#' @noRd
.gl_rs_backfitting_inputs <- function(
    eta,
    eta_dr,
    d1,
    par_name,
    rs_design_cache,
    par_cov,
    par_s,
    score_fn = score_function_v2) {
  score_inputs <- .gl_prepare_rs_score_inputs(
    eta = eta,
    eta_dr = eta_dr,
    d1 = d1,
    par_name = par_name
  )

  score <- score_fn(
    eta = score_inputs$eta,
    dldpar = score_inputs$d1,
    d2ldpar = -(score_inputs$d1 * score_inputs$d1),
    dpardeta = score_inputs$eta_dr
  )

  design_info <- rs_design_cache[[par_name]]
  w_k_vec <- as.vector(score$w_k)
  z_k <- score$z_k
  beta_start <- .gl_rs_beta_start(
    par_name = par_name,
    par_cov = par_cov,
    par_s = par_s,
    design_info = design_info
  )

  list(
    score_inputs = score_inputs,
    score = score,
    d1 = score_inputs$d1,
    eta_dr_vec = score_inputs$eta_dr,
    design_info = design_info,
    w_k_vec = w_k_vec,
    z_k = z_k,
    beta_start = beta_start
  )
}
