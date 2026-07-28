#' @noRd
.copula_hessian_pair_derivatives <- function(u1, u2, fam_num, par1_e, par2_e,
                                             pair_ok, has_zeta) {
  c_val <- .copula_pdf(u1, u2, family = fam_num, par = par1_e, par2 = par2_e)
  c_val[!is.finite(c_val) | c_val <= 0] <- 1

  dcdu1 <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u1", log = FALSE)
  dcdu2 <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "u2", log = FALSE)
  dcdth <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par", log = FALSE)
  dldth <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par", log = TRUE)

  d2cdu1_2 <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "u1")
  d2cdu2_2 <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "u2")
  d2cdth2 <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par")

  hu <- 1e-5
  u1c <- pmax(pmin(u1, 1 - hu), hu)
  u2c <- pmax(pmin(u2, 1 - hu), hu)
  c_pp <- .copula_pdf(u1c + hu, u2c + hu, fam_num, par1_e, par2_e)
  c_pm <- .copula_pdf(u1c + hu, u2c - hu, fam_num, par1_e, par2_e)
  c_mp <- .copula_pdf(u1c - hu, u2c + hu, fam_num, par1_e, par2_e)
  c_mm <- .copula_pdf(u1c - hu, u2c - hu, fam_num, par1_e, par2_e)
  d2cdu1u2 <- (c_pp - c_pm - c_mp + c_mm) / (4 * hu^2)
  d2cdu1u2[!is.finite(d2cdu1u2)] <- 0

  hth <- 1e-5
  dcdu1_p <- .copula_deriv(u1, u2, fam_num, par1_e + hth, par2_e, deriv = "u1", log = FALSE)
  dcdu1_m <- .copula_deriv(u1, u2, fam_num, par1_e - hth, par2_e, deriv = "u1", log = FALSE)
  d2cdthu1 <- (dcdu1_p - dcdu1_m) / (2 * hth)
  d2cdthu1[!is.finite(d2cdthu1)] <- 0

  dcdu2_p <- .copula_deriv(u1, u2, fam_num, par1_e + hth, par2_e, deriv = "u2", log = FALSE)
  dcdu2_m <- .copula_deriv(u1, u2, fam_num, par1_e - hth, par2_e, deriv = "u2", log = FALSE)
  d2cdthu2 <- (dcdu2_p - dcdu2_m) / (2 * hth)
  d2cdthu2[!is.finite(d2cdthu2)] <- 0

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
    dcdz <- .copula_deriv(u1, u2, fam_num, par1_e, par2_e, deriv = "par2", log = FALSE)
    d2cdz2 <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par2")
    d2cdthdz <- .copula_deriv2(u1, u2, fam_num, par1_e, par2_e, deriv = "par1par2")

    hz <- 1e-5 * pmax(1, abs(par2_e))
    par2_p <- par2_e + hz
    par2_m <- pmax(par2_e - hz, 2 + 1e-8)
    hz_eff <- par2_p - par2_m

    dcdu1_zp <- .copula_deriv(u1, u2, fam_num, par1_e, par2_p, deriv = "u1", log = FALSE)
    dcdu1_zm <- .copula_deriv(u1, u2, fam_num, par1_e, par2_m, deriv = "u1", log = FALSE)
    d2cdzu1 <- (dcdu1_zp - dcdu1_zm) / hz_eff
    dcdu2_zp <- .copula_deriv(u1, u2, fam_num, par1_e, par2_p, deriv = "u2", log = FALSE)
    dcdu2_zm <- .copula_deriv(u1, u2, fam_num, par1_e, par2_m, deriv = "u2", log = FALSE)
    d2cdzu2 <- (dcdu2_zp - dcdu2_zm) / hz_eff

    dcdz <- .copula_hessian_zero_invalid(dcdz, pair_ok)
    d2cdz2 <- .copula_hessian_zero_invalid(d2cdz2, pair_ok)
    d2cdthdz <- .copula_hessian_zero_invalid(d2cdthdz, pair_ok)
    d2cdzu1 <- .copula_hessian_zero_invalid(d2cdzu1, pair_ok)
    d2cdzu2 <- .copula_hessian_zero_invalid(d2cdzu2, pair_ok)

    dldz_pair <- .copula_hessian_zero_invalid(dcdz / c_val, pair_ok)
    d2ldz2_pair <- (c_val * d2cdz2 - dcdz^2) / c_val^2
    d2ldthdz_pair <- (d2cdthdz * c_val - dcdth * dcdz) / c_val^2
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

  d2ldth2_pair <- (c_val * d2cdth2 - dcdth^2) / c_val^2
  d2ldth2_pair <- .copula_hessian_zero_invalid(d2ldth2_pair, pair_ok)

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
