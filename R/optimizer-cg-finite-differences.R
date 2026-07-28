#' Finite-difference gradient for CG optimization
#'
#' @noRd
.gl_cg_finite_gradient <- function(beta_vec, base_ll, mm_cg, eval_fn, gradient_method) {
  grad <- rep(0, length(beta_vec))
  names(grad) <- names(beta_vec)
  for (ii in seq_along(beta_vec)) {
    hk <- 1e-5 * max(1, abs(beta_vec[ii]))
    bp <- beta_vec
    bp[ii] <- bp[ii] + hk
    lp <- eval_fn(bp, mm_cg)
    lpv <- if (is.null(lp)) NA_real_ else lp$loglik
    if (identical(gradient_method, "forward")) {
      if (is.finite(lpv) && is.finite(base_ll)) {
        grad[ii] <- (lpv - base_ll) / hk
      } else {
        bm <- beta_vec
        bm[ii] <- bm[ii] - hk
        lm <- eval_fn(bm, mm_cg)
        lmv <- if (is.null(lm)) NA_real_ else lm$loglik
        if (is.finite(lmv) && is.finite(base_ll)) grad[ii] <- (base_ll - lmv) / hk
      }
    } else {
      bm <- beta_vec
      bm[ii] <- bm[ii] - hk
      lm <- eval_fn(bm, mm_cg)
      lmv <- if (is.null(lm)) NA_real_ else lm$loglik
      if (is.finite(lpv) && is.finite(lmv)) {
        grad[ii] <- (lpv - lmv) / (2 * hk)
      } else if (is.finite(lpv) && is.finite(base_ll)) {
        grad[ii] <- (lpv - base_ll) / hk
      } else if (is.finite(lmv) && is.finite(base_ll)) grad[ii] <- (base_ll - lmv) / hk
    }
  }
  grad
}

#' Finite-difference Hessian block for CG optimization
#'
#' @noRd
.gl_cg_finite_hessian_block <- function(beta_vec, block_names, mm_cg, eval_fn, h = 1e-4) {
  block_names <- intersect(block_names, names(beta_vec))
  n_block <- length(block_names)
  H_block <- matrix(NA_real_, n_block, n_block,
    dimnames = list(block_names, block_names)
  )
  if (n_block == 0L) {
    return(H_block)
  }
  eval_base <- eval_fn(beta_vec, mm_cg)
  f0 <- if (is.null(eval_base)) NA_real_ else eval_base$loglik
  if (!is.finite(f0)) {
    return(H_block)
  }

  eval_ll <- function(beta_try) {
    out <- eval_fn(beta_try, mm_cg)
    if (is.null(out)) NA_real_ else out$loglik
  }

  for (ii in seq_len(n_block)) {
    ni <- block_names[ii]
    hi <- h * max(1, abs(beta_vec[ni]))
    bp <- beta_vec
    bm <- beta_vec
    bp[ni] <- bp[ni] + hi
    bm[ni] <- bm[ni] - hi
    fp <- eval_ll(bp)
    fm <- eval_ll(bm)
    if (is.finite(fp) && is.finite(fm)) {
      H_block[ii, ii] <- (fp - 2 * f0 + fm) / (hi^2)
    }

    if (ii < n_block) {
      for (jj in seq.int(ii + 1L, n_block)) {
        nj <- block_names[jj]
        hj <- h * max(1, abs(beta_vec[nj]))
        bpp <- beta_vec
        bpm <- beta_vec
        bmp <- beta_vec
        bmm <- beta_vec
        bpp[ni] <- bpp[ni] + hi
        bpp[nj] <- bpp[nj] + hj
        bpm[ni] <- bpm[ni] + hi
        bpm[nj] <- bpm[nj] - hj
        bmp[ni] <- bmp[ni] - hi
        bmp[nj] <- bmp[nj] + hj
        bmm[ni] <- bmm[ni] - hi
        bmm[nj] <- bmm[nj] - hj
        fpp <- eval_ll(bpp)
        fpm <- eval_ll(bpm)
        fmp <- eval_ll(bmp)
        fmm <- eval_ll(bmm)
        if (all(is.finite(c(fpp, fpm, fmp, fmm)))) {
          H_block[ii, jj] <- (fpp - fpm - fmp + fmm) / (4 * hi * hj)
          H_block[jj, ii] <- H_block[ii, jj]
        }
      }
    }
  }
  H_block
}

#' Check CG Hessian shape and finiteness
#'
