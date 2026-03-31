suppressPackageStartupMessages({library(gamlss2); library(gamlss)})
set.seed(1)
n <- 200
d <- data.frame(
  response = rnorm(n),
  time = sample(1:4, n, TRUE),
  age_new_name = runif(n, 20, 80),
  gender = factor(sample(0:1, n, TRUE))
)
m <- gamlss2(response ~ time + gender + s(age_new_name, bs='ps'), family = NO(), data = d, control = gamlss2_control(maxit = 1))
cat('xterms:', paste(m$xterms$mu, collapse=' | '), '\n')
cat('model cols:', paste(names(m$model), collapse=' | '), '\n')
hasx <- ('x' %in% names(m)) && ('mu' %in% names(m$x))
cat('has x_mu:', hasx, '\n')
if (hasx) cat('x_mu cols:', paste(colnames(m$x$mu), collapse=' | '), '\n')
