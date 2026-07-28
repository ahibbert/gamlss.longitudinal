#' Assemble the full covariate Hessian from per-observation second derivatives

#'

#' `H[a,b] = X_a' diag(w_ab) X_b`, where `w_ab[i] = d2l/(deta_a deta_b)[i]`.

#'

#' @keywords internal

.assemble_covariate_hessian <- function(
    object,
    margin_d1l, # list[[pn]]: per-obs d log margin / d par
    margin_d2l, # list[[pn1]][[pn2]]: per-obs d2 log margin / d par1 d par2
    copula_hess, # output of .calc_copula_hessian_contributions
    eta_dr, # link function derivatives d par/d eta, named list
    eta_d2, # link inverse second derivatives d2 par / d eta2
    mm, # model matrix list (mm$x, mm$s)
    pair_cache # pair_cache for theta_index_map
    ) {
  par_names <- names(object$par)

  n_par <- length(par_names)

  H <- matrix(0,
    nrow = n_par, ncol = n_par,
    dimnames = list(par_names, par_names)
  )

  mm_x <- mm$x

  # Map each covariate coefficient to its distribution parameter block.

  param_block <- vapply(par_names, .hessian_assembly_param_of, character(1))

  row_id1 <- copula_hess$row_id1

  row_id2 <- copula_hess$row_id2

  pair_ok <- copula_hess$pair_ok

  id1_ok <- row_id1[pair_ok]

  id2_ok <- row_id2[pair_ok]

  for (a_idx in seq_len(n_par)) {
    a_name <- par_names[a_idx]

    pa <- param_block[a_name]

    xa <- .hessian_assembly_x_col(a_name, param_block, mm_x)

    if (is.null(xa)) next

    dra <- as.numeric(eta_dr[[pa]])

    for (b_idx in seq_len(a_idx)) { # lower triangle + diagonal

      b_name <- par_names[b_idx]

      pb <- param_block[b_name]

      xb <- .hessian_assembly_x_col(b_name, param_block, mm_x)

      if (is.null(xb)) next

      drb <- as.numeric(eta_dr[[pb]])

      H_ab <- 0

      # Note: theta/zeta design matrices are indexed by theta_index_map[row_id1]

      # (one row per time-1..T-1 observation = same as mm$x$theta row order).

      # Margin design matrices are obs-indexed (n_obs rows).

      copula_parameter_block <- .hessian_assembly_copula_parameter_block(
        pa,
        pb,
        xa,
        xb,
        dra,
        drb,
        copula_hess,
        pair_cache,
        row_id1,
        pair_ok
      )

      if (!is.null(copula_parameter_block)) {
        H_ab <- copula_parameter_block
      } else if ((pa %in% names(margin_d2l) && pb == "theta") ||

        (pb %in% names(margin_d2l) && pa == "theta")) {
        H_ab <- .hessian_assembly_margin_copula_block(
          pa = pa,
          pb = pb,
          target = "theta",
          xa = xa,
          xb = xb,
          dra = dra,
          drb = drb,
          margin_d2l = margin_d2l,
          copula_hess = copula_hess,
          pair_cache = pair_cache,
          row_id1 = row_id1,
          pair_ok = pair_ok,
          id1_ok = id1_ok,
          id2_ok = id2_ok
        )
      } else if (!is.null(copula_hess$cop_d2l_margin_zeta_u1) &&

        ((pa %in% names(margin_d2l) && pb == "zeta") ||

          (pb %in% names(margin_d2l) && pa == "zeta"))) {
        H_ab <- .hessian_assembly_margin_copula_block(
          pa = pa,
          pb = pb,
          target = "zeta",
          xa = xa,
          xb = xb,
          dra = dra,
          drb = drb,
          margin_d2l = margin_d2l,
          copula_hess = copula_hess,
          pair_cache = pair_cache,
          row_id1 = row_id1,
          pair_ok = pair_ok,
          id1_ok = id1_ok,
          id2_ok = id2_ok
        )
      } else {
        H_ab <- .hessian_assembly_margin_margin_block(
          pa = pa,
          pb = pb,
          xa = xa,
          xb = xb,
          dra = dra,
          drb = drb,
          margin_d1l = margin_d1l,
          margin_d2l = margin_d2l,
          copula_hess = copula_hess,
          eta_dr = eta_dr,
          eta_d2 = eta_d2,
          pair_ok = pair_ok,
          id1_ok = id1_ok,
          id2_ok = id2_ok
        )
      }

      H[a_idx, b_idx] <- H_ab

      H[b_idx, a_idx] <- H_ab
    }
  }

  H
}
