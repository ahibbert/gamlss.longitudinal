calc_copula_derivatives <- function(eta_inv, Fx_1_2, copula_dist, calc_d2 = FALSE, calc_d2_marginal = FALSE, par1 = NULL, par2 = NULL, pair_complete = NULL, derivatives = NULL) {
  if (is.null(par1)) {
    par1 <- eta_inv[["theta"]]
  }

  if (is.null(par2)) {
    if ("zeta" %in% names(eta_inv)) {
      par2 <- eta_inv[["zeta"]]
    } else {
      par2 <- eta_inv[["theta"]] * 0
    }
  }

  if (is.null(pair_complete)) {
    pair_complete <- rep(TRUE, length(par1))
  }

  has_zeta <- "zeta" %in% names(eta_inv)
  if (is.null(derivatives)) {
    derivatives <- c("dldth", "dcdth", "dcdu1", "dcdu2")
    if (has_zeta) {
      derivatives <- c(derivatives, "dldz", "dcdz")
    }
  }
  derivatives <- unique(as.character(derivatives))
  want <- function(name) name %in% derivatives
  zero_vec <- function() numeric(length(par1))

  if (length(par1) == 0) {
    if (has_zeta) {
      if (calc_d2 == TRUE) {
        return(list(dldth = numeric(0), dcdth = numeric(0), dldz = numeric(0), dcdz = numeric(0), dcdu1 = numeric(0), dcdu2 = numeric(0), d2ldth2 = numeric(0), d2ldz2 = numeric(0), d2ldthdz = numeric(0), d2cdu12 = numeric(0), d2cdu22 = numeric(0)))
      }
      return(as.list(stats::setNames(rep(list(numeric(0)), length(derivatives)), derivatives)))
    }
    if (calc_d2 == TRUE) {
      return(list(dldth = numeric(0), dcdth = numeric(0), dcdu1 = numeric(0), dcdu2 = numeric(0), d2ldth2 = numeric(0), d2cdu12 = numeric(0), d2cdu22 = numeric(0)))
    }
    return(as.list(stats::setNames(rep(list(numeric(0)), length(derivatives)), derivatives)))
  }

  Fx_eval <- as.matrix(Fx_1_2)
  Fx_eval[!is.finite(Fx_eval)] <- 0.5
  Fx_eval[Fx_eval > 1] <- 1
  Fx_eval[Fx_eval < 0] <- 0

  par1_eval <- par1
  par2_eval <- par2
  par1_eval[!is.finite(par1_eval)] <- 0
  par2_eval[!is.finite(par2_eval)] <- 0

  if (copula_dist == "C") {
    par1_eval[par1_eval >= 28] <- 27.9
  }

  if (identical(.copula_family_code(copula_dist), "t") && !isTRUE(calc_d2)) {
    return(.calc_copula_t_first_derivatives(
      Fx_eval = Fx_eval,
      par1_eval = par1_eval,
      par2_eval = par2_eval,
      pair_complete = pair_complete,
      derivatives = derivatives,
      has_zeta = has_zeta
    ))
  }

  need_copula_d <- want("dcdth") || want("dcdz") || calc_d2
  copula_d <- NULL
  if (need_copula_d) {
    copula_d <- .copula_pdf(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval)
    copula_d[!is.finite(copula_d) | copula_d <= 0] <- 1
    copula_d[!pair_complete] <- 1
  }

  dldth <- dcdth <- NULL
  if (want("dldth") || want("dcdth") || calc_d2) {
    dldth <- .copula_deriv(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "par", log = TRUE)
    dldth[!is.finite(dldth)] <- 0
    dldth[!pair_complete] <- 0
  }
  if (want("dcdth") || calc_d2) {
    dcdth <- copula_d * dldth
    dcdth[!pair_complete] <- 0
  }

  if (calc_d2 == TRUE) {
    d2cdth <- .copula_deriv2(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "par")
    d2ldth2 <- (1 / (copula_d^2)) * (copula_d * d2cdth - dcdth^2)
  }

  dldz <- dcdz <- NULL
  if (has_zeta && (want("dldz") || want("dcdz") || calc_d2)) {
    dldz <- .copula_deriv(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "par2", log = TRUE)
    dldz[!is.finite(dldz)] <- 0
    dldz[!pair_complete] <- 0
  }
  if (has_zeta && (want("dcdz") || calc_d2)) {
    dcdz <- copula_d * dldz
    dcdz[!is.finite(dcdz)] <- 0
    dcdz[!pair_complete] <- 0

    if (calc_d2 == TRUE) {
      d2cdz <- .copula_deriv2(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "par2")
      d2ldz2 <- (1 / (copula_d^2)) * (copula_d * d2cdz - dcdz^2)

      d2cdthdz <- .copula_deriv2(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "par1par2")
      d2ldthdz <- (d2cdthdz * copula_d - dcdth * dcdz) / (copula_d^2)
    }
  }

  dcdu1 <- dcdu2 <- NULL
  if (want("dcdu1") || calc_d2) {
    dcdu1 <- .copula_deriv(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "u1", log = FALSE)
    dcdu1[!pair_complete] <- 0
  }
  if (want("dcdu2") || calc_d2) {
    dcdu2 <- .copula_deriv(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "u2", log = FALSE)
    dcdu2[!pair_complete] <- 0
  }

  if (calc_d2 == TRUE) {
    d2ldth2[!is.finite(d2ldth2)] <- 0
  }

  if (calc_d2 == TRUE) {
    d2cdu12 <- .copula_deriv2(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "u1")
    d2cdu22 <- .copula_deriv2(Fx_eval[, 1], Fx_eval[, 2], family = copula_dist, par = par1_eval, par2 = par2_eval, deriv = "u2")
    d2cdu12[!is.finite(d2cdu12)] <- 0
    d2cdu22[!is.finite(d2cdu22)] <- 0
    d2cdu12[!pair_complete] <- 0
    d2cdu22[!pair_complete] <- 0
  }

  if (has_zeta) {
    if (calc_d2 == TRUE) {
      d2ldz2[!is.finite(d2ldz2)] <- 0
      d2ldthdz[!is.finite(d2ldthdz)] <- 0
      d2ldz2[!pair_complete] <- 0
      d2ldthdz[!pair_complete] <- 0
    }
  }

  if (calc_d2 == TRUE) {
    d2ldth2[!pair_complete] <- 0
  }

  ############# RETURN LIST

  return_list <- list()
  if (want("dldth") || calc_d2) return_list$dldth <- dldth %||% zero_vec()
  if (want("dcdth") || calc_d2) return_list$dcdth <- dcdth %||% zero_vec()
  if (has_zeta && (want("dldz") || calc_d2)) return_list$dldz <- dldz %||% zero_vec()
  if (has_zeta && (want("dcdz") || calc_d2)) return_list$dcdz <- dcdz %||% zero_vec()
  if (want("dcdu1") || calc_d2) return_list$dcdu1 <- dcdu1 %||% zero_vec()
  if (want("dcdu2") || calc_d2) return_list$dcdu2 <- dcdu2 %||% zero_vec()

  if (calc_d2 == TRUE) {
    return_list$d2ldth2 <- d2ldth2
    if (has_zeta) {
      return_list$d2ldz2 <- d2ldz2
      return_list$d2ldthdz <- d2ldthdz
    }
    return_list$d2cdu12 <- d2cdu12
    return_list$d2cdu22 <- d2cdu22
  }

  return(return_list)
}

.calc_copula_t_first_derivatives <- function(
    Fx_eval,
    par1_eval,
    par2_eval,
    pair_complete,
    derivatives,
    has_zeta) {
  want <- function(name) name %in% derivatives
  t_derivs <- character(0)
  if (want("dldth") || want("dcdth")) t_derivs <- c(t_derivs, "par")
  if (has_zeta && (want("dldz") || want("dcdz"))) t_derivs <- c(t_derivs, "par2")
  if (want("dcdu1")) t_derivs <- c(t_derivs, "u1")
  if (want("dcdu2")) t_derivs <- c(t_derivs, "u2")
  t_derivs <- unique(t_derivs)

  if (length(t_derivs) == 0L) {
    return(list())
  }

  log_flags <- c(u1 = FALSE, u2 = FALSE, par = TRUE, par2 = TRUE)
  t_out <- .copula_t_deriv_many(
    Fx_eval[, 1],
    Fx_eval[, 2],
    par1_eval,
    par2_eval,
    derivs = t_derivs,
    log = log_flags[t_derivs]
  )
  copula_d <- attr(t_out, "density")
  copula_d[!is.finite(copula_d) | copula_d <= 0] <- 1
  copula_d[!pair_complete] <- 1

  clean_pair <- function(x) {
    x[!is.finite(x)] <- 0
    x[!pair_complete] <- 0
    x
  }

  out <- list()
  if (want("dldth")) out$dldth <- clean_pair(t_out$par)
  if (want("dcdth")) out$dcdth <- clean_pair(copula_d * t_out$par)
  if (has_zeta && want("dldz")) out$dldz <- clean_pair(t_out$par2)
  if (has_zeta && want("dcdz")) out$dcdz <- clean_pair(copula_d * t_out$par2)
  if (want("dcdu1")) out$dcdu1 <- clean_pair(t_out$u1)
  if (want("dcdu2")) out$dcdu2 <- clean_pair(t_out$u2)

  out
}
