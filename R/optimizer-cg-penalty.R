#' Build CG smooth penalty matrix
#'
#' @noRd
.gl_build_cg_penalty <- function(beta_names, lambda_current, par_s, mm) {
  P <- matrix(0,
    nrow = length(beta_names), ncol = length(beta_names),
    dimnames = list(beta_names, beta_names)
  )
  for (pn in names(par_s)) {
    if (length(par_s[[pn]]) == 0) next
    for (sn in names(par_s[[pn]])) {
      b_names <- names(par_s[[pn]][[sn]])
      idx <- match(b_names, beta_names)
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      B <- mm$s[[pn]][[sn]]
      S <- attr(B, "penalty")
      if (is.null(S) || !is.matrix(S)) {
        D <- diff(diag(ncol(B)), differences = 2)
        S <- t(D) %*% D
      }
      P[idx, idx] <- P[idx, idx] + as.numeric(lambda_current[[pn]][[sn]]) * S
    }
  }
  P
}
