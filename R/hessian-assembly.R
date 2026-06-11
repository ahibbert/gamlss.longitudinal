#' Assemble the full covariate Hessian from per-observation second derivatives

#'

#' `H[a,b] = X_a' diag(w_ab) X_b`, where `w_ab[i] = d2l/(deta_a deta_b)[i]`.

#'

#' @keywords internal

.assemble_covariate_hessian <- function(

    object,

    margin_d1l,           # list[[pn]]: per-obs d log margin / d par

    margin_d2l,           # list[[pn1]][[pn2]]: per-obs d2 log margin / d par1 d par2

    copula_hess,          # output of .calc_copula_hessian_contributions

    eta_dr,               # link function derivatives d par/d eta, named list

    eta_d2,               # link inverse second derivatives d2 par / d eta2

    mm,                   # model matrix list (mm$x, mm$s)

    pair_cache            # pair_cache for theta_index_map

) {

  par_names  <- names(object$par)

  n_par      <- length(par_names)

  H          <- matrix(0, nrow = n_par, ncol = n_par,

                       dimnames = list(par_names, par_names))


  mm_x <- mm$x

  # Map each covariate coefficient to its distribution parameter block

  # e.g. "mu.(Intercept)", "mu.time", "sigma.(Intercept)", ...

  # par_names has format "param.term" or just "param" for intercept.

  param_of <- function(cn) {

    # parameter block is everything before the first dot

    m <- regexpr("^[^.]+", cn)

    regmatches(cn, m)

  }


  param_block <- sapply(par_names, param_of)


  # Helper: get design matrix column for covariate coefficient a

  get_x_col <- function(a_name) {

    blk <- param_block[a_name]

    if (!blk %in% names(mm_x)) return(NULL)

    X <- mm_x[[blk]]

    # column name within the block is the part after "blk."

    col_name <- sub(paste0("^", blk, "\\."), "", a_name)

    if (col_name == a_name) col_name <- "(Intercept)"  # no dot = intercept

    if (!col_name %in% colnames(X)) {

      # Smooth columns can be appended to mm$x as numbered basis columns while

      # the coefficient names retain the full smooth label, e.g.

      # theta.s(s1, k = 10).3 -> design column "3".

      basis_col <- sub("^.*\\.([0-9]+)$", "\\1", col_name)

      if (!identical(basis_col, col_name) && basis_col %in% colnames(X)) {

        col_name <- basis_col

      }

    }

    if (!col_name %in% colnames(X)) return(NULL)

    X[, col_name, drop = TRUE]

  }


  # Helper: per-obs margin second derivative (margin x margin only)

  d2l_par_margin_obs <- function(pn1, pn2) {

    m_pars <- names(margin_d2l)

    n_obs  <- length(eta_dr[[if (pn1 %in% names(eta_dr)) pn1 else pn2]])

    val    <- numeric(n_obs)

    if (pn1 %in% m_pars && pn2 %in% m_pars) {

      if (!is.null(margin_d2l[[pn1]][[pn2]]))

        val <- val + margin_d2l[[pn1]][[pn2]]

      if (!is.null(copula_hess$cop_d2l_margin[[pn1]][[pn2]]))

        val <- val + copula_hess$cop_d2l_margin[[pn1]][[pn2]]

    }

    val

  }


  d1l_par_margin_obs <- function(pn) {

    n_obs <- length(eta_dr[[pn]])

    val <- numeric(n_obs)

    if (pn %in% names(margin_d1l)) val <- val + margin_d1l[[pn]]

    if (!is.null(copula_hess$cop_d1l_margin) && pn %in% names(copula_hess$cop_d1l_margin)) {

      val <- val + copula_hess$cop_d1l_margin[[pn]]

    }

    val

  }


  row_id1 <- copula_hess$row_id1

  row_id2 <- copula_hess$row_id2

  pair_ok <- copula_hess$pair_ok

  id1_ok  <- row_id1[pair_ok]

  id2_ok  <- row_id2[pair_ok]


  for (a_idx in seq_len(n_par)) {

    a_name <- par_names[a_idx]

    pa     <- param_block[a_name]

    xa     <- get_x_col(a_name)

    if (is.null(xa)) next

    dra    <- as.numeric(eta_dr[[pa]])


    for (b_idx in seq_len(a_idx)) {  # lower triangle + diagonal

      b_name <- par_names[b_idx]

      pb     <- param_block[b_name]

      xb     <- get_x_col(b_name)

      if (is.null(xb)) next

      drb    <- as.numeric(eta_dr[[pb]])


      H_ab <- 0


      # Note: theta/zeta design matrices are indexed by theta_index_map[row_id1]

      # (one row per time-1..T-1 observation = same as mm$x$theta row order).

      # Margin design matrices are obs-indexed (n_obs rows).

      if (pa == "theta" && pb == "theta") {

        # cop_d2l_theta is already in eta space (curvature correction applied);

        # do NOT multiply by dra*drb here.

        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]

        valid  <- !is.na(th_idx)

        d2_pairs <- copula_hess$cop_d2l_theta[pair_ok]

        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],

                    na.rm = TRUE)


      } else if (!is.null(copula_hess$cop_d2l_zeta) &&

                 pa == "zeta" && pb == "zeta") {

        # cop_d2l_zeta is already in eta space (curvature correction applied);

        # do NOT multiply by dra*drb here.

        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]

        valid  <- !is.na(th_idx)

        d2_pairs <- copula_hess$cop_d2l_zeta[pair_ok]

        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid],

                    na.rm = TRUE)


      } else if (!is.null(copula_hess$cop_d2l_thetazeta) &&

                 ((pa == "theta" && pb == "zeta") ||

                  (pa == "zeta"  && pb == "theta"))) {

        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]

        valid  <- !is.na(th_idx)

        d2_pairs <- copula_hess$cop_d2l_thetazeta[pair_ok]

        H_ab <- sum(xa[th_idx[valid]] * xb[th_idx[valid]] * d2_pairs[valid] *

                    dra[th_idx[valid]] * drb[th_idx[valid]], na.rm = TRUE)


      } else if ((pa %in% names(margin_d2l) && pb == "theta") ||

                 (pb %in% names(margin_d2l) && pa == "theta")) {

        # Margin x theta cross block: pair-indexed to avoid design-matrix length mismatch.

        # H[sigma.a, theta.b] = sum_k {

        #   X_sigma[i1_k] * dr_sigma[i1_k] * d2logc/(d par_i1 d theta_k) * dr_theta[th_k] * X_theta[th_k]

        # + X_sigma[i2_k] * dr_sigma[i2_k] * d2logc/(d par_i2 d theta_k) * dr_theta[th_k] * X_theta[th_k]

        # }

        mp  <- if (pa %in% names(margin_d2l)) pa else pb  # margin parameter name

        xm  <- if (pa %in% names(margin_d2l)) xa else xb  # margin design column (n_obs)

        drm <- if (pa %in% names(margin_d2l)) dra else drb # margin eta_dr (n_obs)

        # theta design column and eta_dr (n_theta_rows); index via th_idx

        xth <- if (pb == "theta") xb else xa

        drth <- if (pb == "theta") drb else dra


        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]

        valid  <- !is.na(th_idx)

        th_v   <- th_idx[valid]

        id1_v  <- id1_ok[valid]

        id2_v  <- id2_ok[valid]


        u1_vals <- copula_hess$cop_d2l_margin_theta_u1[[mp]][pair_ok][valid]

        u2_vals <- copula_hess$cop_d2l_margin_theta_u2[[mp]][pair_ok][valid]


        H_ab <- sum(

          xm[id1_v] * drm[id1_v] * u1_vals * drth[th_v] * xth[th_v] +

          xm[id2_v] * drm[id2_v] * u2_vals * drth[th_v] * xth[th_v],

          na.rm = TRUE

        )


      } else if (!is.null(copula_hess$cop_d2l_margin_zeta_u1) &&

                 ((pa %in% names(margin_d2l) && pb == "zeta") ||

                  (pb %in% names(margin_d2l) && pa == "zeta"))) {

        # Margin x zeta cross block. This is analogous to margin x theta,

        # but zeta uses the log_2plus inverse link, so include dzeta/deta.

        mp  <- if (pa %in% names(margin_d2l)) pa else pb

        xm  <- if (pa %in% names(margin_d2l)) xa else xb

        drm <- if (pa %in% names(margin_d2l)) dra else drb

        xze <- if (pb == "zeta") xb else xa

        drze <- if (pb == "zeta") drb else dra


        th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]

        valid  <- !is.na(th_idx)

        th_v   <- th_idx[valid]

        id1_v  <- id1_ok[valid]

        id2_v  <- id2_ok[valid]


        u1_vals <- copula_hess$cop_d2l_margin_zeta_u1[[mp]][pair_ok][valid]

        u2_vals <- copula_hess$cop_d2l_margin_zeta_u2[[mp]][pair_ok][valid]


        H_ab <- sum(

          xm[id1_v] * drm[id1_v] * u1_vals * drze[th_v] * xze[th_v] +

          xm[id2_v] * drm[id2_v] * u2_vals * drze[th_v] * xze[th_v],

          na.rm = TRUE

        )


      } else {

        # Margin x margin: obs-level

        d2l_obs <- d2l_par_margin_obs(pa, pb) * dra * drb

        if (identical(pa, pb) && pa %in% names(eta_d2)) {

          d2l_obs <- d2l_obs + d1l_par_margin_obs(pa) * as.numeric(eta_d2[[pa]])

        }

        H_ab <- sum(xa * xb * d2l_obs, na.rm = TRUE)


        # Cross-pair copula contributions (obs i1 in u1 role, obs i2 in u2 role)

        if (pa %in% names(copula_hess$cross_pair_contribs) &&

            pb %in% names(copula_hess$cross_pair_contribs[[pa]])) {

          # Term 1: d2logc / (d pa at i1) (d pb at i2)

          cross_ab <- copula_hess$cross_pair_contribs[[pa]][[pb]][pair_ok]

          # Term 2: d2logc / (d pa at i2) (d pb at i1) = cross_pair_contribs[[pb]][[pa]]

          cross_ba <- copula_hess$cross_pair_contribs[[pb]][[pa]][pair_ok]

          H_ab <- H_ab + sum(

            xa[id1_ok] * dra[id1_ok] * cross_ab * xb[id2_ok] * drb[id2_ok] +

            xa[id2_ok] * dra[id2_ok] * cross_ba * xb[id1_ok] * drb[id1_ok],

            na.rm = TRUE

          )

        }

      }


      H[a_idx, b_idx] <- H_ab

      H[b_idx, a_idx] <- H_ab

    }

  }

  H

}



# ------------------------------------------------------------

# 4.  Margin per-observation second derivative (observed, via FD on log-density)

# ------------------------------------------------------------


#' Compute observed per-observation margin log-density second derivatives

#' via finite differences in the natural (pre-link) parameter space.

#'

#' This replaces the gamlss family d2ld*2 functions which return the

#' expected Fisher information (y-independent) rather than the observed

#' second derivative.

#'

#' Returns a `list[[pn1]][[pn2]]` of length n_obs.

#' @keywords internal

.calc_margin_d2l_fd <- function(eta_inv, mm, margin_dist, response, h = 1e-4) {

  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]

  n_obs <- length(response)

  obs_valid <- is.finite(response)


  dfun_name <- paste0("d", margin_dist$family[1])

  dfun <- eval(parse(text = dfun_name))


  # Build base argument list

  base_args <- list(x = response, log = TRUE)

  for (pn in margin_pars) {

    base_args[[pn]] <- as.numeric(eta_inv[[pn]])

  }


  # log-density at base parameters (all observations)

  ll0 <- do.call(dfun, base_args)

  ll0[!obs_valid] <- 0


  result <- vector("list", length(margin_pars))

  names(result) <- margin_pars

  for (pni in margin_pars) {

    result[[pni]] <- vector("list", length(margin_pars))

    names(result[[pni]]) <- margin_pars

  }


  for (i in seq_along(margin_pars)) {

    pni <- margin_pars[i]

    par_i <- as.numeric(eta_inv[[pni]])

    hi <- .natural_fd_step(par_i, pni, margin_dist, h)


    # Log-density at +hi and -hi for pni

    args_ip <- base_args; args_ip[[pni]] <- par_i + hi

    args_im <- base_args; args_im[[pni]] <- par_i - hi

    ll_ip <- do.call(dfun, args_ip); ll_ip[!obs_valid] <- 0

    ll_im <- do.call(dfun, args_im); ll_im[!obs_valid] <- 0


    # Diagonal: d2l/d pni^2 via 3-point stencil

    d2_diag <- (ll_ip - 2 * ll0 + ll_im) / hi^2

    d2_diag[!is.finite(d2_diag)] <- 0

    result[[pni]][[pni]] <- d2_diag


    for (j in seq_along(margin_pars)) {

      if (j <= i) next  # fill upper triangle, then mirror

      pnj <- margin_pars[j]

      par_j <- as.numeric(eta_inv[[pnj]])

      hj <- .natural_fd_step(par_j, pnj, margin_dist, h)


      args_jp <- base_args; args_jp[[pnj]] <- par_j + hj

      args_jm <- base_args; args_jm[[pnj]] <- par_j - hj


      args_ipjp <- args_jp; args_ipjp[[pni]] <- par_i + hi

      args_ipjm <- args_jm; args_ipjm[[pni]] <- par_i + hi

      args_imjp <- args_jp; args_imjp[[pni]] <- par_i - hi

      args_imjm <- args_jm; args_imjm[[pni]] <- par_i - hi


      ll_pp <- do.call(dfun, args_ipjp); ll_pp[!obs_valid] <- 0

      ll_pm <- do.call(dfun, args_ipjm); ll_pm[!obs_valid] <- 0

      ll_mp <- do.call(dfun, args_imjp); ll_mp[!obs_valid] <- 0

      ll_mm <- do.call(dfun, args_imjm); ll_mm[!obs_valid] <- 0


      d2_cross <- (ll_pp - ll_pm - ll_mp + ll_mm) / (4 * hi * hj)

      d2_cross[!is.finite(d2_cross)] <- 0

      result[[pni]][[pnj]] <- d2_cross

      result[[pnj]][[pni]] <- d2_cross  # symmetric

    }

  }

  result

}


#' Extract per-observation margin log-likelihood second derivatives

#'

#' Returns a `list[[pn1]][[pn2]]`, each of length n_obs, sourced from the

#' gamlss family d2ldm2, d2ldmdd, ... functions already called in

#' calc_likelihood_minimal.

#'

#' @keywords internal

.get_margin_d1l <- function(margin_deriv, mm) {

  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]

  subs <- c(mu = "m", sigma = "d", nu = "v", tau = "t")

  result <- vector("list", length(margin_pars))

  names(result) <- margin_pars

  for (pn in margin_pars) {

    key <- paste0("dld", subs[[pn]])

    if (key %in% names(margin_deriv)) {

      result[[pn]] <- as.numeric(margin_deriv[[key]])

    } else {

      n <- if (length(margin_deriv)) length(margin_deriv[[1]]) else nrow(mm$x[[pn]])

      result[[pn]] <- numeric(n)

    }

  }

  result

}


.get_margin_d2l <- function(margin_deriv, mm) {

  margin_pars <- names(mm$x)[names(mm$x) %in% c("mu", "sigma", "nu", "tau")]

  n_par <- length(margin_pars)


  # Map gamlss derivative name endings to (par1, par2) pairs

  subs <- c(mu = "m", sigma = "d", nu = "v", tau = "t")


  result <- vector("list", n_par)

  names(result) <- margin_pars

  for (pn in margin_pars) {

    result[[pn]] <- vector("list", n_par)

    names(result[[pn]]) <- margin_pars

    for (pn2 in margin_pars) result[[pn]][[pn2]] <- NULL

  }


  for (i in seq_along(margin_pars)) {

    pni <- margin_pars[i]

    si  <- subs[pni]

    for (j in seq_along(margin_pars)) {

      pnj <- margin_pars[j]

      sj  <- subs[pnj]

      if (i == j) {

        key <- paste0("d2ld", si, "2")

      } else {

        # gamlss uses alphabetical order of the two sub-letters

        key1 <- paste0("d2ld", si, "d", sj)

        key2 <- paste0("d2ld", sj, "d", si)

        key  <- if (key1 %in% names(margin_deriv)) key1 else

                if (key2 %in% names(margin_deriv)) key2 else NA_character_

      }

      if (!is.na(key) && key %in% names(margin_deriv)) {

        result[[pni]][[pnj]] <- as.numeric(margin_deriv[[key]])

      } else {

        # Not provided by this family: assume zero cross term

        n <- length(margin_deriv[[1]])

        result[[pni]][[pnj]] <- numeric(n)

      }

    }

  }

  result

}



# ------------------------------------------------------------

# 5.  Top-level analytical Hessian function

# ------------------------------------------------------------


