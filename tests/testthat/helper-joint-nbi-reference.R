joint_nbi_reference <- function(q, object, B) {
  # NBI size=1/sigma, mean=mu; exact Clayton CDF rectangles. This reference
  # calls no package likelihood, link, pair-cache, derivative or rectangle helper.
  mu <- exp(q[["mu.intercept"]] + q[["mu.x"]] * object$dataset$x +
              as.numeric(B %*% q[colnames(B)]))
  sigma <- exp(q[["sigma.intercept"]])
  theta <- exp(q[["theta.intercept"]])
  y <- object$response
  upper <- stats::pnbinom(y, size = 1 / sigma, mu = mu)
  lower <- stats::pnbinom(y - 1, size = 1 / sigma, mu = mu)
  mass <- stats::dnbinom(y, size = 1 / sigma, mu = mu)
  cdf <- function(u, v) {
    out <- numeric(length(u)); positive <- u > 0 & v > 0
    out[positive] <- (u[positive]^(-theta) + v[positive]^(-theta) - 1)^(-1 / theta)
    out
  }
  left <- which(object$response_margin < 3L); right <- left + 1L
  rectangle <- cdf(upper[left], upper[right]) - cdf(lower[left], upper[right]) -
    cdf(upper[left], lower[right]) + cdf(lower[left], lower[right])
  stopifnot(all(is.finite(rectangle)), all(rectangle > 0), all(mass > 0))
  sum(log(mass)) + sum(log(rectangle) - log(mass[left]) - log(mass[right]))
}
