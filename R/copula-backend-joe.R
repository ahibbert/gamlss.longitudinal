.copula_joe_theta <- function(par) {
  pmax(as.numeric(par), 1)
}


.copula_joe_parts <- function(u1, u2, par) {
  vals <- .copula_recycle(.copula_clamp01(u1), .copula_clamp01(u2), .copula_joe_theta(par))

  u1 <- vals[[1]]

  u2 <- vals[[2]]

  theta <- vals[[3]]

  a <- (1 - u1)^theta

  b <- (1 - u2)^theta

  s <- a + b - a * b

  list(u1 = u1, u2 = u2, theta = theta, a = a, b = b, s = s, cdf = 1 - s^(1 / theta))
}


.copula_joe_cdf <- function(u1, u2, par) {
  .copula_clamp01(.copula_joe_parts(u1, u2, par)$cdf)
}


.copula_joe_pdf <- function(u1, u2, par) {
  p <- .copula_joe_parts(u1, u2, par)

  au <- (1 - p$u1)^(p$theta - 1)

  bu <- (1 - p$u2)^(p$theta - 1)

  du <- 1 - p$a

  dv <- 1 - p$b

  m <- 1 / p$theta - 1

  out <- au * bu * p$s^(m - 1) * (p$theta * p$s + (p$theta - 1) * du * dv)

  out[!is.finite(out)] <- 0

  out
}


.copula_joe_hfunc1 <- function(u1, u2, par) {
  p <- .copula_joe_parts(u1, u2, par)

  out <- (1 - p$u1)^(p$theta - 1) * (1 - p$b) * p$s^(1 / p$theta - 1)

  .copula_clamp01(out)
}
