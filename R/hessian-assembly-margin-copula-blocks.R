#' Assemble margin-by-copula Hessian cross blocks
#'
#' Margin parameters are observation-indexed while theta/zeta parameters are
#' pair-indexed through `theta_index_map`. This helper keeps that index
#' translation visible and shared by the theta and zeta cross blocks.
#'
#' @keywords internal
#' @noRd
.hessian_assembly_margin_copula_block <- function(
    pa,
    pb,
    target,
    xa,
    xb,
    dra,
    drb,
    margin_d2l,
    copula_hess,
    pair_cache,
    row_id1,
    pair_ok,
    id1_ok,
    id2_ok) {
  mp <- if (pa %in% names(margin_d2l)) pa else pb
  xm <- if (pa %in% names(margin_d2l)) xa else xb
  drm <- if (pa %in% names(margin_d2l)) dra else drb
  xcp <- if (pb == target) xb else xa
  drcp <- if (pb == target) drb else dra

  th_idx <- pair_cache$theta_index_map[row_id1[pair_ok]]
  valid <- !is.na(th_idx)
  th_v <- th_idx[valid]
  id1_v <- id1_ok[valid]
  id2_v <- id2_ok[valid]

  u1_name <- paste0("cop_d2l_margin_", target, "_u1")
  u2_name <- paste0("cop_d2l_margin_", target, "_u2")
  u1_vals <- copula_hess[[u1_name]][[mp]][pair_ok][valid]
  u2_vals <- copula_hess[[u2_name]][[mp]][pair_ok][valid]

  sum(
    xm[id1_v] * drm[id1_v] * u1_vals * drcp[th_v] * xcp[th_v] +
      xm[id2_v] * drm[id2_v] * u2_vals * drcp[th_v] * xcp[th_v],
    na.rm = TRUE
  )
}
