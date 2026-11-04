#' Legacy archived outer optimizer helper
#'
#' Retained for compatibility and internal historical support. This helper is
#' not part of the core reviewer path for fitting, prediction, or diagnostics.
#'
#' @noRd
NULL

########## ARCHIVE ###########

# Given a parameter vector starting values par = (mu,sigma,nu,tau,theta,zeta), return best fit parameters
#' @keywords internal
#' @noRd
optim_outer <- function(par, dataset, margin_dist, copula_dist,
                        step_size = 0.1, verbose = TRUE, use_dlcopdpar = TRUE) {
  # print("THIS FUNCTION ASSUMES RESPONSE IS ORDERED AS TIME, SUBJECT | PAR INPUT MUST BE NAMED")

  copula_input <- get_copula_dist(copula_dist)
  copula_number <- copula_input$copula_dist
  copula_link <- copula_input$copula_link

  num_margins <- length(unique(dataset$time))
  margin_names <- unique(dataset$time)
  response <- dataset$response

  # Set up parameter vector so names are consistent with the distributions

  if (all(is.null(names(par)) | is.na(names(par)))) {
    stop("ERROR: par vector must be named")
  }
  margin_par <- par[names(par) %in% c("mu", "sigma", "nu", "tau")]
  copula_par <- par[!names(par) %in% c("mu", "sigma", "nu", "tau")]

  ##### Calculate all relevant derivatives / CG method with first and second derivatives

  ### Calculate margin derivatives w.r.t. margin parameters

  # Get names for margin derivatives from margin_dist
  n_par <- length(eta[[par_name]])
  d1_full <- matrix(0, nrow = n_par, ncol = 1)

  # Get link transforms (eta) and derivatives w.r.t to link for parameters
  if (n_par == length(dataset$response)) {
    par_idx <- row_id1
  } else {
    margin_names <- sort(unique(dataset$time))
    theta_rows <- which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names) - 1))])
    theta_index_map <- rep(NA_integer_, length(dataset$response))
    theta_index_map[theta_rows] <- seq_along(theta_rows)
    par_idx <- theta_index_map[row_id1]
  }

  valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
  d1_full[par_idx[valid_idx], 1] <- dldth[valid_idx]
  if (par_name %in% names(margin_par)) {
    par_eta[par_name] <- margin_dist[[paste(par_name, ".linkfun", sep = "")]](par[par_name])
    par_eta_dr[par_name] <- margin_dist[[paste(par_name, ".dr", sep = "")]](par_eta[par_name])
  }
  if (par_name %in% names(copula_par)) {
    n_par <- length(eta[[par_name]])
    d1_full <- matrix(0, nrow = n_par, ncol = 1)
    par_eta_dr[par_name] <- copula_link[[paste(par_name, ".dr", sep = "")]](par_eta[par_name])
  }
  if (n_par == length(dataset$response)) {
    par_idx <- row_id1
  } else {
    margin_names <- sort(unique(dataset$time))
    theta_rows <- which(dataset$time %in% margin_names[seq_len(max(1, length(margin_names) - 1))])
    theta_index_map <- rep(NA_integer_, length(dataset$response))
    theta_index_map[theta_rows] <- seq_along(theta_rows)
    par_idx <- theta_index_map[row_id1]
  }

  valid_idx <- is.finite(par_idx) & par_idx >= 1 & par_idx <= nrow(d1_full)
  d1_full[par_idx[valid_idx], 1] <- dldz[valid_idx]
  # Setup input matrix of response and parameters
  margin_deriv_input <- list()
  margin_deriv_input[["y"]] <- response
  margin_deriv_input[["q"]] <- response
  margin_deriv_input[["x"]] <- response
  for (par_name in names(margin_par)) {
    margin_deriv_input[[par_name]] <- rep(margin_par[par_name], length(response))
  }

  # Calculate all derivatives
  margin_deriv <- list()
  for (deriv_name in margin_deriv_names) {
    FUN <- margin_dist[[deriv_name]]
    FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% formalArgs(FUN)]
    margin_deriv[[deriv_name]] <- do.call(FUN, args = margin_deriv_input[FUN_args])
    margin_deriv[[deriv_name]][!is.finite(margin_deriv[[deriv_name]])] <- 0
  }

  margin_pFUN <- eval(parse(text = paste("p", margin_dist$family[1], sep = "")))
  FUN <- margin_pFUN
  FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% formalArgs(FUN)]
  margin_p <- do.call(FUN, args = margin_deriv_input[FUN_args])

  margin_dFUN <- eval(parse(text = paste("d", margin_dist$family[1], sep = "")))
  FUN <- margin_dFUN
  FUN_args <- names(margin_deriv_input)[names(margin_deriv_input) %in% formalArgs(FUN)]
  margin_d <- do.call(FUN, args = margin_deriv_input[FUN_args])

  ### Calculate copula derivatives w.r.t. copula parameters

  # First calculate margin F(x1), F(x2) as inputs to copula

  Fx_1_2 <- matrix(ncol = 2, nrow = 0)
  order_copula <- data.frame()
  for (i in 1:(num_margins - 1)) {
    Fx_1_2 <- rbind(Fx_1_2, cbind(margin_p[dataset$time == margin_names[i]], margin_p[dataset$time == margin_names[i + 1]]))
    order_copula <- rbind(order_copula, cbind(dataset[dataset$time == margin_names[i], c("time", "subject")], dataset[dataset$time == margin_names[i + 1], c("time", "subject")]))
  }
  names(order_copula) <- c("time1", "subject1", "time2", "subject2")

  par1 <- copula_par["theta"]
  if (is.na(copula_par["zeta"])) {
    par2 <- 0
  } else {
    par2 <- copula_par["zeta"]
  }

  # Handling extreme values
  Fx_1_2[Fx_1_2 > 1] <- 1
  Fx_1_2[Fx_1_2 < 0] <- 0

  if (copula_number == 3) {
    if (par1 > 28) {
      par1 <- 28
    }
  }

  copula_d <- .copula_pdf(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2)
  dldth <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par", log = TRUE)
  dcdth <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par", log = FALSE)
  d2cdth <- .copula_deriv2(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par")
  d2ldth2 <- (1 / (copula_d^2)) * (copula_d * d2cdth - dcdth^2)
  if (!is.na(copula_par["zeta"])) {
    dldz <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par2", log = TRUE)
    dcdz <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par2", log = FALSE)
    d2cdz <- .copula_deriv2(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par2")
    d2ldz2 <- (1 / (copula_d^2)) * (copula_d * d2cdz - dcdz^2)

    d2cdthdz <- .copula_deriv2(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "par1par2")
    d2ldthdz <- (d2cdthdz * copula_d - dcdth * dcdz) / (copula_d^2)
  }
  dcdu1 <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "u1", log = FALSE)
  dcdu2 <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "u2", log = FALSE)

  d2cdu12 <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "u1", log = FALSE)
  d2cdu22 <- .copula_deriv(Fx_1_2[, 1], Fx_1_2[, 2], family = copula_number, par = par1, par2 = par2, deriv = "u2", log = FALSE)

  d2ldth2[!is.finite(d2ldth2)] <- 0

  ### Calculate copula derivatives w.r.t margin parameters

  # Extract margin calculations for F(x), f(x), response and derivatives at time 1 and time 2, join to copula values for time 1 and time 2
  margin_deriv_1 <- margin_deriv_2 <- margin_deriv_2cross <- matrix(ncol = length(margin_par), nrow = length(response))
  for (i in 1:length(margin_par)) {
    margin_deriv_1[, i] <- margin_deriv[grepl("dld", names(margin_deriv))][[i]]
    margin_deriv_2[, i] <- margin_deriv[grepl("d2ld", names(margin_deriv)) & endsWith(names(margin_deriv), "2")][[i]]
  }
  colnames(margin_deriv_1) <- paste("dld", names(margin_par), sep = "")
  colnames(margin_deriv_2) <- paste(paste("d2ld", names(margin_par), sep = ""), "2", sep = "")

  # colnames(margin_deriv_2)=paste("d2ld",names(margin_par),sep="")

  order_margin <- dataset[, c("time", "subject")]
  margin_components <- cbind(order_margin, response, margin_p, margin_d, margin_deriv_1, margin_deriv_2)
  margin_components_Ft_plus <- margin_components
  margin_components_Ft_plus$time <- normalize_lag_time(margin_components_Ft_plus$time)
  margin_plus <- merge(margin_components, margin_components_Ft_plus, by = c("time", "subject"), all.x = TRUE)

  copula_components <- cbind(order_copula, dcdu1, dcdu2, copula_d, d2cdu12, d2cdu22)
  copula_merged <- merge(copula_components, margin_plus, by.x = c("time1", "subject1"), by.y = c("time", "subject"), all.x = TRUE)

  # Calculate copula derivative with respect to marginal parameters
  input <- copula_merged
  dlcopdpar <- matrix(0, nrow = nrow(input), ncol = length(margin_par))
  d2lcopdpar2 <- matrix(0, nrow = nrow(input), ncol = length(margin_par))

  i <- 1
  for (par_name in names(margin_par)) {
    # Take parameters from input for clarity
    dc_tplus_du_t <- input[, "dcdu1"]
    dc_tplus_du_tplus <- input[, "dcdu2"]
    l_t <- input[, paste(paste("dld", par_name, sep = ""), ".x", sep = "")]
    l_t_plus <- input[, paste(paste("dld", par_name, sep = ""), ".y", sep = "")]
    x_t <- input[, "response.x"]
    x_t_plus <- input[, "response.y"]
    f_t <- input[, "margin_d.x"]
    f_t_plus <- input[, "margin_d.y"]
    du_t_dmu <- x_t * f_t * l_t
    du_t_plus_dmu <- x_t_plus * f_t_plus * l_t_plus
    c_tplus <- input[, "copula_d"]

    du_t_dmu <- x_t * f_t * l_t
    du_t_plus_dmu <- x_t_plus * f_t_plus * l_t_plus

    dc_plus_dt_dmu <- dc_tplus_du_t * du_t_dmu
    dc_plus_dt_plus_dmu <- dc_tplus_du_tplus * du_t_plus_dmu
    dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)] <- 0
    dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)] <- 0
    dcdmu_tplus <- ((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
    dcdmu_tplus[is.nan(dcdmu_tplus) | is.na(dcdmu_tplus)] <- 0

    dlcopdpar[, i] <- dcdmu_tplus

    ####### NOW FOR SECOND DERIVATIVE OF COPULA TERM

    l2_t <- input[, paste(paste(paste("d2ld", par_name, sep = ""), "2", sep = ""), ".x", sep = "")]
    l2_tplus <- input[, paste(paste(paste("d2ld", par_name, sep = ""), "2", sep = ""), ".y", sep = "")]

    df_t_dmu <- f_t * l_t
    df_t_plus_dmu <- f_t_plus * l_t_plus

    d2f_t_dmu <- df_t_dmu * l_t + f_t * l2_t
    d2f_t_plus_dmu <- df_t_plus_dmu * l_t_plus + f_t_plus * l2_tplus

    d2u_t_dmu2 <- x_t * d2f_t_dmu
    d2u_t_plus_dmu2 <- x_t_plus * d2f_t_plus_dmu

    d2cdu_t2 <- input[, "d2cdu12"]
    d2cdu_t_plus2 <- input[, "d2cdu22"]
    d2cdu_t2[is.nan(d2cdu_t2)] <- 0
    d2cdu_t_plus2[is.nan(d2cdu_t_plus2)] <- 0

    d2cdmu2 <- d2cdu_t2 * du_t_dmu^2 + dc_tplus_du_t * d2u_t_dmu2 + d2cdu_t_plus2 * du_t_plus_dmu^2 + dc_tplus_du_tplus * d2u_t_plus_dmu2

    d2lcdmu2 <- as.matrix((d2cdmu2 * c_tplus - (dcdmu_tplus^2)) / (c_tplus^2))

    d2lcopdpar2[, i] <- d2lcdmu2
    # num_deriv=margin_copula_merged_2[,"num_dlcopdpar_ordered.Ft"]
    # num_deriv_nolog=margin_copula_merged_2[,"num_dlcopdpar_nolog_ordered.Ft"]

    i <- i + 1
  }
  colnames(dlcopdpar) <- paste("dlcopd", names(margin_par), sep = "")
  colnames(d2lcopdpar2) <- paste(paste("d2lcd", names(margin_par), sep = ""), "2", sep = "")

  dlcopdpar[!is.finite(dlcopdpar)] <- 0
  d2lcopdpar2[!is.finite(d2lcopdpar2)] <- 0

  #### Define score and hessian

  score <- par * 0
  hessian <- matrix(0, nrow = length(par), ncol = length(par))
  colnames(hessian) <- names(par)
  rownames(hessian) <- names(par)
  names(score) <- names(par)

  margin_deriv_sum <- vector()
  for (i in 1:length(margin_deriv)) {
    margin_deriv[[i]][!is.finite(margin_deriv[[i]])] <- 0
    margin_deriv_sum[i] <- sum(margin_deriv[[i]])
  }
  names(margin_deriv_sum) <- names(margin_deriv)

  margin_d1 <- margin_deriv_sum[grepl("dld", names(margin_deriv))]
  margin_d2 <- margin_deriv_sum[grepl("d2ld", names(margin_deriv)) & endsWith(names(margin_deriv), "2")]
  margin_d2d <- margin_deriv_sum[grepl("d2ld", names(margin_deriv)) & !endsWith(names(margin_deriv), "2")]

  if (is.na(copula_par["zeta"])) {
    copula_d1 <- sum(dldth)
    copula_d2 <- sum(d2ldth2)
  } else {
    copula_d1 <- colSums(cbind(dldth, dldz))
    copula_d2 <- colSums(cbind(d2ldth2, d2ldz2))
  }
  margin_d1_dlcopdpar <- margin_d1 + if (use_dlcopdpar == TRUE) {
    colSums(dlcopdpar)
  } else {
    colSums(dlcopdpar) * 0
  }
  margin_d2_dlcopdpar <- margin_d2 + if (use_dlcopdpar == TRUE) {
    colSums(d2lcopdpar2) * 0
  } else {
    colSums(d2lcopdpar2) * 0
  }
  score <- c(margin_d1_dlcopdpar, copula_d1)

  ### CALCULATING HESSIAN USING D2
  diag(hessian) <- c(margin_d2_dlcopdpar, copula_d2)
  hessian[1:length(margin_par), 1:length(margin_par)][upper.tri(hessian[1:length(margin_par), 1:length(margin_par)])] <- margin_d2d
  hessian[1:length(margin_par), 1:length(margin_par)][lower.tri(hessian[1:length(margin_par), 1:length(margin_par)])] <- margin_d2d

  # Why isn't d2 for copula negative?
  copula_hess <- hessian[(length(margin_par) + 1):(length(margin_par) + length(copula_par)), (length(margin_par) + 1):(length(margin_par) + length(copula_par))]
  if (!is.na(copula_par["zeta"])) {
    copula_hess[upper.tri(copula_hess)] <- sum(d2ldthdz)
    copula_hess[lower.tri(copula_hess)] <- sum(d2ldthdz)
  }
  hessian[(length(margin_par) + 1):(length(margin_par) + length(copula_par)), (length(margin_par) + 1):(length(margin_par) + length(copula_par))] <- copula_hess

  ### STILL NEED TO CALCULATE d2 for marginal parameters with respect to copula likelihood and add to hessian values

  # par_end=par-(solve(-hessian)%*%(score))

  # weights_eta=diag((1/(score_eta^2)))
  # weights=-diag(score*score)

  # score=score
  # weights=-solve(hessian)

  weights_eta <- -solve(hessian * par_eta_dr * par_eta_dr)

  # weights_eta=diag(1/(score*score*par_eta_dr*par_eta_dr))
  score_eta <- score * par_eta_dr
  # par_end=par*(1-step_size) + step_size*(par+par_change)

  par_end <- par * 0
  names(par_end) <- names(par_eta_end) <- names(par)
  # Get end paraemters re-transformed
  # for (par_name in names(par)) {
  #  if(par_name %in% names(margin_par)) {
  #    par_end[par_name]=margin_dist[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  #  if(par_name %in% names(copula_par)) {
  #    par_end[par_name]=copula_link[[paste(par_name,".linkinv",sep="")]](par_end[par_name])
  #  }
  # }
  ### If calculating for eta
  for (par_name in names(par)) {
    if (par_name %in% names(margin_par)) {
      par_end[par_name] <- margin_dist[[paste(par_name, ".linkinv", sep = "")]](par_end[par_name])
    }
    if (par_name %in% names(copula_par)) {
      par_end[par_name] <- copula_link[[paste(par_name, ".linkinv", sep = "")]](par_end[par_name])
    }
  }

  sum_log_margin_p <- sum(log(margin_d)[is.finite(log(margin_d))])
  sum_log_copula_d <- sum(log(copula_d)[is.finite(log(copula_d))])

  log_lik <- c(sum_log_copula_d, sum_log_margin_p, sum_log_copula_d + sum_log_margin_p)
  names(log_lik) <- c("copula", "margin", "joint")

  if (verbose == TRUE) {
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

  return(list(score = score, hessian = hessian, par_end = par_end, par_eta_end = par_eta_end, par_start = par, log_lik = log_lik))
}

# Load analytical Hessian helpers when sourcing this file directly (development workflow).
local({
  hessian_files <- c(
    "hessian-linkinv-derivatives.R",
    "hessian-fd-step.R",
    "hessian-warnings.R",
    "hessian-margin-cdf.R",
    "hessian-copula.R",
    "hessian-assembly.R",
    "hessian-analytical.R"
  )
  candidate_dirs <- c("R", file.path(getwd(), "R"))
  for (dir in candidate_dirs) {
    candidates <- file.path(dir, hessian_files)
    if (all(file.exists(candidates))) {
      for (p in candidates) source(p, local = FALSE)
      break
    }
  }
})
