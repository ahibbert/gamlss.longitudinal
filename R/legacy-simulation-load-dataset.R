#' Legacy simulation and research-era helpers
#'
#' These helpers are retained for compatibility and internal historical support.
#' They are not part of the core reviewer path for fitting, prediction, or diagnostics.
#'
#' @noRd
NULL

loadDataset <- function(simOption = 5, plot_dist = FALSE, n = 100, d = 3, copula_dist = NA, margin_dist, copula.link = NA, par.copula, par.margin, covariates_input = NA) {
  if (simOption == 1) {
    stop(
      "Legacy RAND data loading has been retired and is not shipped. ",
      "Use a caller-supplied public dataset or one of the simulation options instead.",
      call. = FALSE
    )
  } else if (simOption == 2) {
    # set up D-vine copula model with mixed pair-copulas
    d <- 3
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(2, 2, 0)
    par <- c(logit_inv(.8), logit_inv(.8), logit_inv(.8))
    par2 <- c(log_2plus_inv(2.1), log_2plus_inv(2.1), log_2plus_inv(2.1))

    # transform to R-vine matrix notation
    RVM <- .copula_dvine(order, family, par, par2)
    contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n * t, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      margin[covariates[[3]] == 0, i] <- qZISICHEL(copsim[, i], mu = exp(0.3 + 0.2 * i), sigma = exp(0.3 + 0.2 * i), nu = -0.8, tau = 0.05)[covariates[[3]] == 0] # Update to i*mu/sigma as needed
      margin[covariates[[3]] == 1, i] <- qZISICHEL(copsim[, i], mu = exp(0.3 + 0.2 * i + 0.1), sigma = exp(0.3 + 0.2 * i + 0.1), nu = -0.8, tau = 0.05)[covariates[[3]] == 1] # Update to i*mu/sigma as needed
    }

    response <- as.data.frame(margin)

    dataset <- create_longitudinal_dataset(response, covariates, labels = c("subject", "time", "response", "age", "year", "gender"))
  } else if (simOption == 3) {
    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    if (length(par.copula) == d - 1) {
      par <- c(copula.link$theta.linkinv(par.copula), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(copula.link$theta.linkinv(par.copula[1:(length(par.copula) / 2)]), rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(copula.link$zeta.linkinv(par.copula[(length(par.copula) / 2 + 1):(length(par.copula))]), rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    # contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 4) {
    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    if (length(par.copula) == d - 1) {
      par <- c(rep(par.copula, d - 1), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(rep(par.copula["theta"], d - 1), rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"], d - 1), rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    # contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 5) {
    copula_input <- get_copula_dist(copula_dist)
    copula.family <- copula_input$copula_dist

    qFUN <- paste("q", margin_dist$family[1], sep = "")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    if (length(par.copula) == d - 1) {
      par <- c(rep(par.copula, d - 1), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(rep(par.copula["theta"], d - 1), rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"], d - 1), rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    # contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 6) {
    t <- d
    margin_sim <- matrix(0, ncol = d, nrow = n)

    for (i in 1:d) {
      margin_sim[, i] <- rnorm(n, 1, 3)
    }

    W <- rnorm(n, 0, 3)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin_sim_out <- matrix(0, ncol = d, nrow = n)
    for (i in 1:d) {
      margin_sim_out[, i] <- margin_sim[, i] + W
    }

    response <- as.data.frame(margin_sim_out)
  } else if (simOption == 7) {
    t <- d
    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    copula_input <- get_copula_dist(copula_dist)
    copula.family <- copula_input$copula_dist

    qFUN <- paste("q", margin_dist$family[1], sep = "")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    # par.copula=c(.3); names(par.copula)=c("theta")
    theta_intercept <- unlist(par.copula["theta"])
    theta_out <- theta_intercept + matrix(rep(covariates_input$theta.time * 1:(d - 1), n), ncol = d - 1, byrow = TRUE) +
      matrix(rep(as.matrix(covariates_input$theta.age * ((covariates[[1]] - 50) / 100)^2), d - 1), ncol = d - 1)

    theta_inv <- copula_input$copula_link$theta.linkinv(theta_out)

    if (length(par.copula) == d - 1) {
      par <- c((par.copula), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(par.copula[1:(length(par.copula) / 2)], rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(par.copula[(length(par.copula) / 2 + 1):(length(par.copula))], rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notatio

    RVM <- list()

    for (i in 1:n) {
      RVM[[i]] <- .copula_dvine(order, c(rep(copula.family, length(theta_inv[i, ])), rep(0, dd - (length(theta_inv[i, ])))), par = c(theta_inv[i, ], rep(0, dd - (length(theta_inv[i, ])))), par2 = c(theta_inv[i, ], rep(0, dd - (length(theta_inv[i, ])))))
    }
    # RVM <- .copula_dvine(order, rep(family[1],nrow(theta_inv)), theta_inv, theta_inv*0)
    # contour(RVM)

    copsim <- .copula_rvine_sim(n, RVM)

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 8) {
    # Multivariate Gamma

    U <- margin_sim <- matrix(0, ncol = d, nrow = n)

    a <- .25
    b <- 1.75
    mu <- rep(1, d)
    W <- rbeta(n, shape1 = a, shape2 = b)

    for (i in 1:d) {
      U[, i] <- rgamma(n, shape = a + b, rate = 1 / mu[i])
      margin_sim[, i] <- U[, i] * W
    }

    # Fake covariates
    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = d, nrow = n)) * (1:d)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    response <- margin_sim
  } else if (simOption == 9) { ######## TIME VARIANT MU

    copula_input <- get_copula_dist(copula_dist)
    copula.family <- copula_input$copula_dist

    qFUN <- paste("q", margin_dist$family[1], sep = "")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    if (length(par.copula) == 1) {
      par <- c(rep(par.copula, d - 1), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(rep(par.copula["theta"], d - 1), rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"], d - 1), rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    # contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 9) { ######## TIME VARIANT MU

    copula_input <- get_copula_dist(copula_dist)
    copula.family <- copula_input$copula_dist

    qFUN <- paste("q", margin_dist$family[1], sep = "")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    if (length(par.copula) == 1) {
      par <- c(rep(par.copula, d - 1), rep(0, dd - (d - 1)))
      par2 <- par * 0
    } else {
      par <- c(rep(par.copula["theta"], d - 1), rep(0, dd - (d - 1))) #+1*1:(d-1)
      par2 <- c(rep(par.copula["zeta"], d - 1), rep(0, dd - (d - 1))) #+0.5*1:(d-1)
    }

    # transform to R-vine matrix notation

    RVM <- .copula_dvine(order, family, par, par2)
    # contour(RVM)

    t <- d
    copsim <- .copula_rvine_sim(n, RVM)

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 1), 0)) # Gender

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = par.margin[1], sigma = par.margin[2], nu = par.margin[3], tau = par.margin[4]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  } else if (simOption == 10) { ######## TIME VARIANT SIGMA AND MU

    print("WARNING: SIMULATION MAPS MARGIN AND COPULA PARAMETERS THROUGH LINK-INVERSE FUNCTIONS AFTER ADDING COVARIATE EFFECTS.")

    t <- d

    # Setup covariates from covariates_input

    covariates <- list()
    covariates[[1]] <- as.data.frame(round(runif(n, 0, 100), 0)) # Age
    covariates[[2]] <- t(t(matrix(1, ncol = t, nrow = n)) * (1:t)) # Time
    covariates[[3]] <- as.data.frame(round(runif(n, 0, 2), 0)) # Gender

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

    copula_input <- get_copula_dist(copula_dist)
    copula.family <- copula_input$copula_dist

    apply_margin_link <- function(par_name, par_value, eta_component) {
      if (par_name %in% names(margin_dist$parameters)) {
        linkfun_name <- paste0(par_name, ".linkfun")
        linkinv_name <- paste0(par_name, ".linkinv")
        par_eta_base <- eval(parse(text = paste0("margin_dist$", linkfun_name)))(par_value)
        par_eta <- par_eta_base + eta_component
        return(eval(parse(text = paste0("margin_dist$", linkinv_name)))(par_eta))
      }

      return(NULL)
    }

    mu_out <- NULL
    sigma_out <- NULL
    nu_out <- NULL
    tau_out <- NULL

    if ("mu" %in% names(margin_dist$parameters)) {
      mu_eta <- make_time_factor_component(covariates_input$mu.time, d, "mu.time") +
        matrix(rep(as.matrix(covariates_input$mu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$mu.gender, d, "mu.gender")
      mu_out <- apply_margin_link("mu", par.margin[1], mu_eta)
    }

    if ("sigma" %in% names(margin_dist$parameters)) {
      sigma_eta <- make_time_factor_component(covariates_input$sigma.time, d, "sigma.time") +
        matrix(rep(as.matrix(covariates_input$sigma.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$sigma.gender, d, "sigma.gender")
      sigma_out <- apply_margin_link("sigma", par.margin[2], sigma_eta)
    }

    if ("nu" %in% names(margin_dist$parameters)) {
      nu_eta <- make_time_factor_component(covariates_input$nu.time, d, "nu.time") +
        matrix(rep(as.matrix(covariates_input$nu.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$nu.gender, d, "nu.gender")
      nu_out <- apply_margin_link("nu", par.margin[3], nu_eta)
    }

    if ("tau" %in% names(margin_dist$parameters)) {
      tau_eta <- make_time_factor_component(covariates_input$tau.time, d, "tau.time") +
        matrix(rep(as.matrix(covariates_input$tau.age * ((covariates[[1]] - 50) / 100)^2), d), ncol = d) +
        make_gender_factor_component(covariates_input$tau.gender, d, "tau.gender")
      tau_out <- apply_margin_link("tau", par.margin[4], tau_eta)
    }
    theta_eta_out <- par.copula[1] + make_time_factor_component(covariates_input$theta.time, d - 1, "theta.time") +
      matrix(rep(as.matrix(covariates_input$theta.age * ((covariates[[1]] - 50) / 100)^2), d - 1), ncol = d - 1) +
      make_gender_factor_component(covariates_input$theta.gender, d - 1, "theta.gender")
    theta_out <- copula_input$copula_link$theta.linkinv(theta_eta_out)

    if ("zeta" %in% copula_input$parameters) {
      zeta_eta_out <- par.copula[2] + make_time_factor_component(covariates_input$zeta.time, d - 1, "zeta.time") +
        matrix(rep(as.matrix(covariates_input$zeta.age * ((covariates[[1]] - 50) / 100)^2), d - 1), ncol = d - 1) +
        make_gender_factor_component(covariates_input$zeta.gender, d - 1, "zeta.gender")
      zeta_out <- copula_input$copula_link$zeta.linkinv(zeta_eta_out)
    } else {
      zeta_out <- matrix(0, nrow = n, ncol = d - 1)
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
      if (is.null(x)) {
        return(paste0(label, ": [NA, NA]"))
      }
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

    # Define margin distribution
    qFUN <- paste("q", margin_dist$family[1], sep = "")

    # set up D-vine copula model with mixed pair-copulas
    d <- d
    dd <- d * (d - 1) / 2
    order <- 1:d
    family <- c(rep(copula.family, d - 1), rep(0, dd - (d - 1)))

    # OK now for each row in theta_out and zeta_out, we need to create a new RVM and simulate from it, then apply the qFUN with the appropriate parameters to get the margin values for that row. This is going to be computationally intensive but should work.

    # row-specific copula simulation from theta_out / zeta_out
    copsim <- matrix(NA_real_, nrow = n, ncol = d)
    for (r in 1:n) {
      par_r <- c(as.numeric(theta_out[r, ]), rep(0, dd - (d - 1)))
      par2_r <- c(
        if ("zeta" %in% copula_input$parameters) as.numeric(zeta_out[r, ]) else rep(0, d - 1),
        rep(0, dd - (d - 1))
      )
      copsim[r, ] <- as.numeric(.copula_rvine_sim(1, .copula_dvine(order, family, par_r, par2_r)))
    }

    margin <- matrix(0, ncol = ncol(copsim), nrow = nrow(copsim))
    for (i in 1:ncol(copsim)) {
      input_list <- list(
        p = copsim[, i],
        mu = mu_out[, i], sigma = sigma_out[, i], nu = nu_out[, i], tau = tau_out[, i]
      )
      args <- names(input_list)[names(input_list) %in% formalArgs(qFUN)]
      qFunOutput_1 <- do.call(qFUN, args = (input_list[args]))
      # input_list=list(p=copsim[,i],mu=par.margin[1],sigma=exp(0.3+0.2*i+0.1),nu=-0.8,tau=0.1)
      # args=names(input_list)[names(input_list)%in%formalArgs(qFUN)]
      # qFunOutput_2=do.call(qFUN,args=(input_list[args]))

      margin[, i] <- qFunOutput_1 # Update to i*mu/sigma as needed
      # margin[covariates[[3]]==0,i]=qFunOutput_1[covariates[[3]]==0]#Update to i*mu/sigma as needed
      # margin[covariates[[3]]==1,i]=qFunOutput_2[covariates[[3]]==1]#Update to i*mu/sigma as needed

      response <- as.data.frame(margin)
    }
  }

  dataset <- create_longitudinal_dataset(response, covariates, labels = c("subject", "time", "response", "age", "year", "gender"))

  if (plot_dist == TRUE) {
    plot_dist(
      dataset,
      margin_dist = margin_dist,
      subject_var = "subject",
      time_var = "time",
      response_var = "response"
    )
  }

  return(dataset)
}
