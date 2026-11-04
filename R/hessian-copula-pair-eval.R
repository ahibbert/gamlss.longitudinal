#' @noRd
.copula_hessian_pair_derivatives <- function(u1, u2, fam_num, par1_e, par2_e,
                                             pair_ok, has_zeta,
                                             u1_active = TRUE, u2_active = TRUE,
                                             theta_active = TRUE, zeta_active = TRUE) {
  inputs <- .copula_recycle(u1, u2, par1_e, par2_e, pair_ok)
  u1 <- inputs[[1L]]; u2 <- inputs[[2L]]
  par1_e <- inputs[[3L]]; par2_e <- inputs[[4L]]; pair_ok <- inputs[[5L]]
  # Excluded pairs must not be sent as missing probabilities or invalid
  # parameters to a delegated backend. Their derivatives are masked below.
  u1[!pair_ok] <- u2[!pair_ok] <- 0.5
  par1_e[!pair_ok] <- if (fam_num %in% c("G", "J")) 1.5 else
    if (fam_num %in% c("C", "F")) 1 else 0
  par2_e[!pair_ok] <- if (has_zeta) 5 else 0
  c_val <- .copula_pdf(u1, u2, family = fam_num, par = par1_e, par2 = par2_e)
  # A backend's DBL_MIN floor cannot supply reliable density-scale curvature.
  # Include log/exp round-trip rounding around the normal floating-point floor.
  if (any(pair_ok & (!is.finite(c_val) |
      c_val <= .Machine$double.xmin * (1 + 1e-12)))) {
    stop("Continuous copula density is nonfinite, nonpositive or underflowed for an included pair.", call. = FALSE)
  }
  c_val[!is.finite(c_val) | c_val <= 0] <- 1

  dcdu1 <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u1", log = FALSE)
  dcdu2 <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u2", log = FALSE)
  parameter_log <- .copula_hessian_parameter_log_curvature(u1, u2,
    fam_num, par1_e, par2_e, c_val)
  dldth <- parameter_log$theta
  dcdth <- c_val * dldth

  argument_curvature <- .copula_density_argument_curvature(u1, u2, fam_num, par1_e, par2_e)
  if (identical(.copula_backend(), "native")) {
    dcdu1 <- argument_curvature$du
    dcdu2 <- argument_curvature$dv
  }
  d2cdu1_2 <- argument_curvature$duu
  d2cdu2_2 <- argument_curvature$dvv
  d2cdth2 <- c_val * (parameter_log$theta2 + dldth^2)

  d2cdu1u2 <- argument_curvature$duv

  parameter_mixed <- .copula_parameter_argument_mixed(u1, u2, fam_num, par1_e, par2_e)
  d2cdthu1 <- parameter_mixed$u1
  d2cdthu2 <- parameter_mixed$u2
  new_curvatures <- list(d2cdu1_2, d2cdu2_2, d2cdu1u2, d2cdthu1, d2cdthu2,
    dldth, parameter_log$theta2)
  if (any(vapply(new_curvatures, function(x) any(pair_ok & !is.finite(x)), logical(1)))) {
    stop("Continuous copula curvature is nonfinite for an included pair.", call. = FALSE)
  }

  dcdu1 <- .copula_hessian_zero_invalid(dcdu1, pair_ok)
  dcdu2 <- .copula_hessian_zero_invalid(dcdu2, pair_ok)
  dcdth <- .copula_hessian_zero_invalid(dcdth, pair_ok)
  dldth <- .copula_hessian_zero_invalid(dldth, pair_ok)
  d2cdu1_2 <- .copula_hessian_zero_invalid(d2cdu1_2, pair_ok)
  d2cdu2_2 <- .copula_hessian_zero_invalid(d2cdu2_2, pair_ok)
  d2cdth2 <- .copula_hessian_zero_invalid(d2cdth2, pair_ok)
  d2cdu1u2 <- .copula_hessian_zero_invalid(d2cdu1u2, pair_ok)
  d2cdthu1 <- .copula_hessian_zero_invalid(d2cdthu1, pair_ok)
  d2cdthu2 <- .copula_hessian_zero_invalid(d2cdthu2, pair_ok)
  c_val[!pair_ok] <- 1

  zeta_derivs <- NULL
  if (has_zeta) {
    dcdz <- c_val * parameter_log$zeta

    zeta_mixed <- .copula_parameter_argument_mixed(
      u1, u2, fam_num, par1_e, par2_e, parameter = "par2"
    )
    d2cdzu1 <- zeta_mixed$u1
    d2cdzu2 <- zeta_mixed$u2
    if (any(pair_ok & (!is.finite(d2cdzu1) | !is.finite(d2cdzu2)))) {
      stop("Continuous copula df/argument curvature is nonfinite for an included pair.", call. = FALSE)
    }

    dcdz <- .copula_hessian_zero_invalid(dcdz, pair_ok)
    d2cdzu1 <- .copula_hessian_zero_invalid(d2cdzu1, pair_ok)
    d2cdzu2 <- .copula_hessian_zero_invalid(d2cdzu2, pair_ok)

    if (any(pair_ok & (!is.finite(parameter_log$zeta) |
        !is.finite(parameter_log$zeta2) | !is.finite(parameter_log$thetazeta)))) {
      stop("Continuous copula parameter curvature is nonfinite for an included pair.", call. = FALSE)
    }
    dldz_pair <- .copula_hessian_zero_invalid(parameter_log$zeta, pair_ok)
    d2ldz2_pair <- parameter_log$zeta2
    d2ldthdz_pair <- parameter_log$thetazeta
    d2ldz2_pair <- .copula_hessian_zero_invalid(d2ldz2_pair, pair_ok)
    d2ldthdz_pair <- .copula_hessian_zero_invalid(d2ldthdz_pair, pair_ok)

    zeta_derivs <- list(
      dcdz = dcdz,
      d2cdzu1 = d2cdzu1,
      d2cdzu2 = d2cdzu2,
      dldz_pair = dldz_pair,
      d2ldz2_pair = d2ldz2_pair,
      d2ldthdz_pair = d2ldthdz_pair
    )
  }

  d2ldth2_pair <- parameter_log$theta2
  d2ldth2_pair <- .copula_hessian_zero_invalid(d2ldth2_pair, pair_ok)

  # The density is evaluated at effective backend arguments. Outside a clamp
  # its composition with the original margin/parameter is locally constant.
  # These masks are supplied by full-model assembly; direct density-derivative
  # callers retain derivatives with respect to the effective arguments.
  u1_active <- rep_len(u1_active, length(u1)) & pair_ok
  u2_active <- rep_len(u2_active, length(u2)) & pair_ok
  theta_active <- rep_len(theta_active, length(u1)) & pair_ok
  zeta_active <- rep_len(zeta_active, length(u1)) & pair_ok
  dcdu1[!u1_active] <- d2cdu1_2[!u1_active] <- 0
  dcdu2[!u2_active] <- d2cdu2_2[!u2_active] <- 0
  d2cdu1u2[!(u1_active & u2_active)] <- 0
  dcdth[!theta_active] <- dldth[!theta_active] <- d2ldth2_pair[!theta_active] <- 0
  d2cdthu1[!(theta_active & u1_active)] <- 0
  d2cdthu2[!(theta_active & u2_active)] <- 0
  if (has_zeta) {
    zeta_derivs$dcdz[!zeta_active] <- 0
    zeta_derivs$dldz_pair[!zeta_active] <- 0
    zeta_derivs$d2ldz2_pair[!zeta_active] <- 0
    zeta_derivs$d2cdzu1[!(zeta_active & u1_active)] <- 0
    zeta_derivs$d2cdzu2[!(zeta_active & u2_active)] <- 0
    zeta_derivs$d2ldthdz_pair[!(theta_active & zeta_active)] <- 0
  }

  list(
    c_val = c_val,
    dcdu1 = dcdu1,
    dcdu2 = dcdu2,
    dcdth = dcdth,
    dldth = dldth,
    d2cdu1_2 = d2cdu1_2,
    d2cdu2_2 = d2cdu2_2,
    d2cdu1u2 = d2cdu1u2,
    d2cdthu1 = d2cdthu1,
    d2cdthu2 = d2cdthu2,
    d2ldth2_pair = d2ldth2_pair,
    zeta = zeta_derivs
  )
}
