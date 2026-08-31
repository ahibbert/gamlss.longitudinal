#' Clamp perturbed natural margin parameters inside common support boundaries
#'
#' @noRd
.discrete_hessian_clamp_margin_value <- function(value, par_name, margin_dist) {
  link_name <- margin_dist[[paste(par_name, "link", sep = ".")]]
  family <- as.character(margin_dist$family[1])

  if (is.character(link_name) && identical(link_name[1], "logit")) {
    return(pmin(pmax(value, 1e-8), 1 - 1e-8))
  }

  if (is.character(link_name) && identical(link_name[1], "log")) {
    return(pmax(value, 1e-8))
  }

  if (par_name %in% c("mu", "sigma", "tau") && !identical(family, "BI")) {
    return(pmax(value, 1e-8))
  }

  value
}

#' Construct plus/minus finite-difference values for a natural parameter
#'
#' @noRd
.discrete_hessian_margin_pm <- function(value, par_name, margin_dist, h) {
  step <- .natural_fd_step(value, par_name, margin_dist, h)
  plus <- .discrete_hessian_clamp_margin_value(value + step, par_name, margin_dist)
  minus <- .discrete_hessian_clamp_margin_value(value - step, par_name, margin_dist)
  if (!is.finite(plus) || !is.finite(minus) || abs(plus - minus) < .Machine$double.eps) {
    plus <- value + h
    minus <- value - h
  }
  c(plus = plus, minus = minus)
}

#' Safe scalar finite-difference calculations
#'
#' @noRd
.discrete_hessian_first <- function(gp, gm, xp, xm) {
  out <- (gp - gm) / (xp - xm)
  if (is.finite(out)) out else 0
}

#' @noRd
.discrete_hessian_second <- function(gp, g0, gm, xp, xm) {
  step <- (xp - xm) / 2
  out <- (gp - 2 * g0 + gm) / (step^2)
  if (is.finite(out)) out else 0
}

#' @noRd
.discrete_hessian_cross <- function(gpp, gpm, gmp, gmm, xp, xm, yp, ym) {
  out <- (gpp - gpm - gmp + gmm) / ((xp - xm) * (yp - ym))
  if (is.finite(out)) out else 0
}

#' Evaluate one endpoint's discrete CDF bounds and PMF
#'
#' @noRd
.discrete_hessian_endpoint_eval <- function(
    obs_id,
    eta_inv,
    margin_dist,
    response,
    pfun,
    dfun,
    override = NULL) {
  y <- response[[obs_id]]
  if (!is.finite(y)) {
    return(list(upper = NA_real_, lower = NA_real_, mass = NA_real_))
  }

  n <- length(response)
  args <- list(q = y, x = y, y = y)
  for (pn in names(eta_inv)[names(eta_inv) %in% c("mu", "sigma", "nu", "tau")]) {
    args[[pn]] <- eta_inv[[pn]][[obs_id]]
  }
  if (length(override) > 0L) {
    for (nm in names(override)) args[[nm]] <- override[[nm]]
  }
  args <- c(args, lapply(.gl_margin_fixed_family_args(margin_dist, n), function(x) x[[obs_id]]))

  upper_args <- args
  upper_args$q <- y
  lower_args <- args
  lower_args$q <- y - 1
  density_args <- args
  density_args$x <- y

  upper <- tryCatch(
    do.call(pfun, upper_args[names(upper_args) %in% formalArgs(pfun)]),
    error = function(e) NA_real_
  )
  lower <- if (is.finite(y) && y <= 0) {
    0
  } else {
    tryCatch(
      do.call(pfun, lower_args[names(lower_args) %in% formalArgs(pfun)]),
      error = function(e) NA_real_
    )
  }
  mass <- tryCatch(
    do.call(dfun, density_args[names(density_args) %in% formalArgs(dfun)]),
    error = function(e) NA_real_
  )

  list(
    upper = as.numeric(upper)[1],
    lower = as.numeric(lower)[1],
    mass = as.numeric(mass)[1]
  )
}

#' Evaluate one adjacent-pair discrete copula ratio log-likelihood
#'
#' @noRd
.discrete_hessian_pair_logratio <- function(
    pair_id,
    eta_inv,
    margin_dist,
    copula_dist,
    response,
    row_id1,
    row_id2,
    par1,
    par2,
    pfun,
    dfun,
    endpoint1_override = NULL,
    endpoint2_override = NULL,
    theta = NULL,
    zeta = NULL) {
  e1 <- .discrete_hessian_endpoint_eval(
    row_id1[[pair_id]], eta_inv, margin_dist, response, pfun, dfun,
    override = endpoint1_override
  )
  e2 <- .discrete_hessian_endpoint_eval(
    row_id2[[pair_id]], eta_inv, margin_dist, response, pfun, dfun,
    override = endpoint2_override
  )

  th <- theta %||% par1[[pair_id]]
  ze <- zeta %||% par2[[pair_id]]
  rect <- .copula_rectangle_prob(
    e1$upper, e2$upper, e1$lower, e2$lower,
    family = copula_dist, par = th, par2 = ze
  )

  out <- log(pmax(rect, 1e-300)) -
    log(pmax(e1$mass, 1e-300)) -
    log(pmax(e2$mass, 1e-300))
  if (is.finite(out)) as.numeric(out) else 0
}

#' Pair-level Hessian contributions for exact discrete rectangle likelihoods
#'
#' @noRd
.calc_discrete_rectangle_hessian_contributions <- function(
    eta_inv,
    eta_dr,
    eta_d2,
    pair_cache,
    margin_dist,
    copula_dist,
    response,
    calc_lik,
    mm,
    h = 1e-4) {
  n_obs <- length(response)
  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]
  row_id1 <- calc_lik$copula_row_id1
  row_id2 <- calc_lik$copula_row_id2
  n_pairs <- length(row_id1)
  has_zeta <- "zeta" %in% names(eta_inv)

  accum <- .copula_hessian_margin_accumulators(margin_pars, n_obs, n_pairs, has_zeta)
  cross_pair_contribs <- vector("list", length(margin_pars))
  names(cross_pair_contribs) <- margin_pars
  for (pn in margin_pars) {
    cross_pair_contribs[[pn]] <- vector("list", length(margin_pars))
    names(cross_pair_contribs[[pn]]) <- margin_pars
    for (pn2 in margin_pars) cross_pair_contribs[[pn]][[pn2]] <- numeric(n_pairs)
  }

  cop_d2l_theta <- numeric(n_pairs)
  cop_d2l_zeta <- if (has_zeta) numeric(n_pairs) else NULL
  cop_d2l_thetazeta <- if (has_zeta) numeric(n_pairs) else NULL

  pfun <- get(paste0("p", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), mode = "function")
  dfun <- get(paste0("d", margin_dist$family[1]), envir = asNamespace("gamlss.dist"), mode = "function")

  par1 <- calc_lik$copula_par1
  par2 <- calc_lik$copula_par2
  par1[!is.finite(par1)] <- 0
  par2[!is.finite(par2)] <- 0

  eval_pair <- function(k, e1 = NULL, e2 = NULL, theta = NULL, zeta = NULL) {
    .discrete_hessian_pair_logratio(
      pair_id = k,
      eta_inv = eta_inv,
      margin_dist = margin_dist,
      copula_dist = copula_dist,
      response = response,
      row_id1 = row_id1,
      row_id2 = row_id2,
      par1 = par1,
      par2 = par2,
      pfun = pfun,
      dfun = dfun,
      endpoint1_override = e1,
      endpoint2_override = e2,
      theta = theta,
      zeta = zeta
    )
  }

  pair_ok <- calc_lik$pair_complete
  for (k in seq_len(n_pairs)) {
    if (!isTRUE(pair_ok[[k]])) next
    i1 <- row_id1[[k]]
    i2 <- row_id2[[k]]
    g0 <- eval_pair(k)

    for (pn in margin_pars) {
      v1 <- eta_inv[[pn]][[i1]]
      pm1 <- .discrete_hessian_margin_pm(v1, pn, margin_dist, h)
      gp <- eval_pair(k, e1 = stats::setNames(list(pm1[["plus"]]), pn))
      gm <- eval_pair(k, e1 = stats::setNames(list(pm1[["minus"]]), pn))
      accum$cop_d1l_margin[[pn]][i1] <- accum$cop_d1l_margin[[pn]][i1] +
        .discrete_hessian_first(gp, gm, pm1[["plus"]], pm1[["minus"]])
      accum$cop_d2l_margin[[pn]][[pn]][i1] <- accum$cop_d2l_margin[[pn]][[pn]][i1] +
        .discrete_hessian_second(gp, g0, gm, pm1[["plus"]], pm1[["minus"]])

      v2 <- eta_inv[[pn]][[i2]]
      pm2 <- .discrete_hessian_margin_pm(v2, pn, margin_dist, h)
      gp <- eval_pair(k, e2 = stats::setNames(list(pm2[["plus"]]), pn))
      gm <- eval_pair(k, e2 = stats::setNames(list(pm2[["minus"]]), pn))
      accum$cop_d1l_margin[[pn]][i2] <- accum$cop_d1l_margin[[pn]][i2] +
        .discrete_hessian_first(gp, gm, pm2[["plus"]], pm2[["minus"]])
      accum$cop_d2l_margin[[pn]][[pn]][i2] <- accum$cop_d2l_margin[[pn]][[pn]][i2] +
        .discrete_hessian_second(gp, g0, gm, pm2[["plus"]], pm2[["minus"]])

      for (pn2 in margin_pars) {
        if (pn2 != pn) {
          v1b <- eta_inv[[pn2]][[i1]]
          pm1b <- .discrete_hessian_margin_pm(v1b, pn2, margin_dist, h)
          gpp <- eval_pair(k,
            e1 = list(
              stats::setNames(list(pm1[["plus"]]), pn),
              stats::setNames(list(pm1b[["plus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gpm <- eval_pair(k,
            e1 = list(
              stats::setNames(list(pm1[["plus"]]), pn),
              stats::setNames(list(pm1b[["minus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gmp <- eval_pair(k,
            e1 = list(
              stats::setNames(list(pm1[["minus"]]), pn),
              stats::setNames(list(pm1b[["plus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gmm <- eval_pair(k,
            e1 = list(
              stats::setNames(list(pm1[["minus"]]), pn),
              stats::setNames(list(pm1b[["minus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          same1 <- .discrete_hessian_cross(
            gpp, gpm, gmp, gmm,
            pm1[["plus"]], pm1[["minus"]],
            pm1b[["plus"]], pm1b[["minus"]]
          )
          accum$cop_d2l_margin[[pn]][[pn2]][i1] <- accum$cop_d2l_margin[[pn]][[pn2]][i1] + same1

          v2c <- eta_inv[[pn2]][[i2]]
          pm2c <- .discrete_hessian_margin_pm(v2c, pn2, margin_dist, h)
          gpp <- eval_pair(k,
            e2 = list(
              stats::setNames(list(pm2[["plus"]]), pn),
              stats::setNames(list(pm2c[["plus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gpm <- eval_pair(k,
            e2 = list(
              stats::setNames(list(pm2[["plus"]]), pn),
              stats::setNames(list(pm2c[["minus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gmp <- eval_pair(k,
            e2 = list(
              stats::setNames(list(pm2[["minus"]]), pn),
              stats::setNames(list(pm2c[["plus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          gmm <- eval_pair(k,
            e2 = list(
              stats::setNames(list(pm2[["minus"]]), pn),
              stats::setNames(list(pm2c[["minus"]]), pn2)
            ) |> unlist(recursive = FALSE)
          )
          same2 <- .discrete_hessian_cross(
            gpp, gpm, gmp, gmm,
            pm2[["plus"]], pm2[["minus"]],
            pm2c[["plus"]], pm2c[["minus"]]
          )
          accum$cop_d2l_margin[[pn]][[pn2]][i2] <- accum$cop_d2l_margin[[pn]][[pn2]][i2] + same2
        }

        v2b <- eta_inv[[pn2]][[i2]]
        pm2b <- .discrete_hessian_margin_pm(v2b, pn2, margin_dist, h)
        gpp <- eval_pair(k,
          e1 = stats::setNames(list(pm1[["plus"]]), pn),
          e2 = stats::setNames(list(pm2b[["plus"]]), pn2)
        )
        gpm <- eval_pair(k,
          e1 = stats::setNames(list(pm1[["plus"]]), pn),
          e2 = stats::setNames(list(pm2b[["minus"]]), pn2)
        )
        gmp <- eval_pair(k,
          e1 = stats::setNames(list(pm1[["minus"]]), pn),
          e2 = stats::setNames(list(pm2b[["plus"]]), pn2)
        )
        gmm <- eval_pair(k,
          e1 = stats::setNames(list(pm1[["minus"]]), pn),
          e2 = stats::setNames(list(pm2b[["minus"]]), pn2)
        )
        cross_pair_contribs[[pn]][[pn2]][k] <- .discrete_hessian_cross(
          gpp, gpm, gmp, gmm,
          pm1[["plus"]], pm1[["minus"]],
          pm2b[["plus"]], pm2b[["minus"]]
        )
      }

      th <- par1[[k]]
      th_pm <- c(plus = th + h * max(1, abs(th)), minus = th - h * max(1, abs(th)))
      if (identical(copula_dist, "C")) th_pm[["minus"]] <- max(th_pm[["minus"]], 1e-8)
      gp <- eval_pair(k, e1 = stats::setNames(list(pm1[["plus"]]), pn), theta = th_pm[["plus"]])
      gm <- eval_pair(k, e1 = stats::setNames(list(pm1[["minus"]]), pn), theta = th_pm[["plus"]])
      hp <- eval_pair(k, e1 = stats::setNames(list(pm1[["plus"]]), pn), theta = th_pm[["minus"]])
      hm <- eval_pair(k, e1 = stats::setNames(list(pm1[["minus"]]), pn), theta = th_pm[["minus"]])
      accum$cop_d2l_margin_theta_u1[[pn]][k] <- .discrete_hessian_cross(
        gp, hp, gm, hm,
        pm1[["plus"]], pm1[["minus"]],
        th_pm[["plus"]], th_pm[["minus"]]
      )

      gp <- eval_pair(k, e2 = stats::setNames(list(pm2[["plus"]]), pn), theta = th_pm[["plus"]])
      gm <- eval_pair(k, e2 = stats::setNames(list(pm2[["minus"]]), pn), theta = th_pm[["plus"]])
      hp <- eval_pair(k, e2 = stats::setNames(list(pm2[["plus"]]), pn), theta = th_pm[["minus"]])
      hm <- eval_pair(k, e2 = stats::setNames(list(pm2[["minus"]]), pn), theta = th_pm[["minus"]])
      accum$cop_d2l_margin_theta_u2[[pn]][k] <- .discrete_hessian_cross(
        gp, hp, gm, hm,
        pm2[["plus"]], pm2[["minus"]],
        th_pm[["plus"]], th_pm[["minus"]]
      )

      if (has_zeta) {
        ze <- par2[[k]]
        ze_pm <- c(plus = ze + h * max(1, abs(ze)), minus = max(ze - h * max(1, abs(ze)), 2 + 1e-8))
        gp <- eval_pair(k, e1 = stats::setNames(list(pm1[["plus"]]), pn), zeta = ze_pm[["plus"]])
        gm <- eval_pair(k, e1 = stats::setNames(list(pm1[["minus"]]), pn), zeta = ze_pm[["plus"]])
        hp <- eval_pair(k, e1 = stats::setNames(list(pm1[["plus"]]), pn), zeta = ze_pm[["minus"]])
        hm <- eval_pair(k, e1 = stats::setNames(list(pm1[["minus"]]), pn), zeta = ze_pm[["minus"]])
        accum$cop_d2l_margin_zeta_u1[[pn]][k] <- .discrete_hessian_cross(
          gp, hp, gm, hm,
          pm1[["plus"]], pm1[["minus"]],
          ze_pm[["plus"]], ze_pm[["minus"]]
        )

        gp <- eval_pair(k, e2 = stats::setNames(list(pm2[["plus"]]), pn), zeta = ze_pm[["plus"]])
        gm <- eval_pair(k, e2 = stats::setNames(list(pm2[["minus"]]), pn), zeta = ze_pm[["plus"]])
        hp <- eval_pair(k, e2 = stats::setNames(list(pm2[["plus"]]), pn), zeta = ze_pm[["minus"]])
        hm <- eval_pair(k, e2 = stats::setNames(list(pm2[["minus"]]), pn), zeta = ze_pm[["minus"]])
        accum$cop_d2l_margin_zeta_u2[[pn]][k] <- .discrete_hessian_cross(
          gp, hp, gm, hm,
          pm2[["plus"]], pm2[["minus"]],
          ze_pm[["plus"]], ze_pm[["minus"]]
        )
      }
    }

    th <- par1[[k]]
    th_pm <- c(plus = th + h * max(1, abs(th)), minus = th - h * max(1, abs(th)))
    if (identical(copula_dist, "C")) th_pm[["minus"]] <- max(th_pm[["minus"]], 1e-8)
    gp <- eval_pair(k, theta = th_pm[["plus"]])
    gm <- eval_pair(k, theta = th_pm[["minus"]])
    d1_th <- .discrete_hessian_first(gp, gm, th_pm[["plus"]], th_pm[["minus"]])
    d2_th <- .discrete_hessian_second(gp, g0, gm, th_pm[["plus"]], th_pm[["minus"]])
    th_idx <- if (length(eta_inv[["theta"]]) == n_obs) {
      i1
    } else {
      pair_cache$theta_index_map[[i1]]
    }
    if (is.finite(th_idx) && th_idx >= 1 && th_idx <= length(eta_dr$theta)) {
      cop_d2l_theta[[k]] <- d2_th * eta_dr$theta[[th_idx]]^2 + d1_th * eta_d2$theta[[th_idx]]
    }

    if (has_zeta) {
      ze <- par2[[k]]
      ze_pm <- c(plus = ze + h * max(1, abs(ze)), minus = max(ze - h * max(1, abs(ze)), 2 + 1e-8))
      gp <- eval_pair(k, zeta = ze_pm[["plus"]])
      gm <- eval_pair(k, zeta = ze_pm[["minus"]])
      d1_ze <- .discrete_hessian_first(gp, gm, ze_pm[["plus"]], ze_pm[["minus"]])
      d2_ze <- .discrete_hessian_second(gp, g0, gm, ze_pm[["plus"]], ze_pm[["minus"]])
      if (is.finite(th_idx) && th_idx >= 1 && th_idx <= length(eta_dr$zeta)) {
        cop_d2l_zeta[[k]] <- d2_ze * eta_dr$zeta[[th_idx]]^2 + d1_ze * eta_d2$zeta[[th_idx]]
      }

      gpp <- eval_pair(k, theta = th_pm[["plus"]], zeta = ze_pm[["plus"]])
      gpm <- eval_pair(k, theta = th_pm[["plus"]], zeta = ze_pm[["minus"]])
      gmp <- eval_pair(k, theta = th_pm[["minus"]], zeta = ze_pm[["plus"]])
      gmm <- eval_pair(k, theta = th_pm[["minus"]], zeta = ze_pm[["minus"]])
      cop_d2l_thetazeta[[k]] <- .discrete_hessian_cross(
        gpp, gpm, gmp, gmm,
        th_pm[["plus"]], th_pm[["minus"]],
        ze_pm[["plus"]], ze_pm[["minus"]]
      )
    }
  }

  list(
    cop_d1l_margin = accum$cop_d1l_margin,
    cop_d2l_margin = accum$cop_d2l_margin,
    cop_d2l_theta = cop_d2l_theta,
    cop_d2l_zeta = cop_d2l_zeta,
    cop_d2l_thetazeta = cop_d2l_thetazeta,
    cop_d2l_margin_theta_u1 = accum$cop_d2l_margin_theta_u1,
    cop_d2l_margin_theta_u2 = accum$cop_d2l_margin_theta_u2,
    cop_d2l_margin_zeta_u1 = accum$cop_d2l_margin_zeta_u1,
    cop_d2l_margin_zeta_u2 = accum$cop_d2l_margin_zeta_u2,
    cross_pair_contribs = cross_pair_contribs,
    row_id1 = row_id1,
    row_id2 = row_id2,
    pair_ok = pair_ok
  )
}
