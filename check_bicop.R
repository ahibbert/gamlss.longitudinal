library(VineCopula)
library(gamlss.dist)
fam <- 1; par <- 0.5
mu <- 5; sigma <- 0.5; nu <- 1; y1 <- 4.5; y2 <- 4.8
h <- 1e-5

u1v <- pGG(y1, mu=mu, sigma=sigma, nu=nu)
u2v <- pGG(y2, mu=mu, sigma=sigma, nu=nu)
cat("u1:", u1v, "u2:", u2v, "\n")

logc <- function(u1, u2) log(BiCopPDF(u1, u2, fam, par))

# True d2 log c(F(y1,s), F(y2,s)) / ds^2 at s=sigma
logc_s <- function(s) logc(pGG(y1, mu=mu, sigma=s, nu=nu), pGG(y2, mu=mu, sigma=s, nu=nu))
true_d2 <- (logc_s(sigma+h) - 2*logc_s(sigma) + logc_s(sigma-h)) / h^2
cat("True d2 log c / dsigma^2:", true_d2, "\n")

# Analytical components
dF1 <- (pGG(y1, mu=mu, sigma=sigma+h, nu=nu) - pGG(y1, mu=mu, sigma=sigma-h, nu=nu)) / (2*h)
dF2 <- (pGG(y2, mu=mu, sigma=sigma+h, nu=nu) - pGG(y2, mu=mu, sigma=sigma-h, nu=nu)) / (2*h)
d2F1 <- (pGG(y1, mu=mu, sigma=sigma+h, nu=nu) + pGG(y1, mu=mu, sigma=sigma-h, nu=nu) - 2*u1v) / h^2
d2F2 <- (pGG(y2, mu=mu, sigma=sigma+h, nu=nu) + pGG(y2, mu=mu, sigma=sigma-h, nu=nu) - 2*u2v) / h^2

cv <- BiCopPDF(u1v, u2v, fam, par)
dc_u1 <- BiCopDeriv(u1v, u2v, fam, par, deriv="u1")
dc_u2 <- BiCopDeriv(u1v, u2v, fam, par, deriv="u2")
d2c_u1 <- BiCopDeriv2(u1v, u2v, fam, par, deriv="u1")
d2c_u2 <- BiCopDeriv2(u1v, u2v, fam, par, deriv="u2")
hu <- 1e-5
d2c_u1u2 <- (BiCopPDF(u1v+hu, u2v+hu, fam, par) - BiCopPDF(u1v+hu, u2v-hu, fam, par) -
             BiCopPDF(u1v-hu, u2v+hu, fam, par) + BiCopPDF(u1v-hu, u2v-hu, fam, par)) / (4*hu^2)

d2logc_u1 <- d2c_u1/cv - (dc_u1/cv)^2
d2logc_u2 <- d2c_u2/cv - (dc_u2/cv)^2
d2logc_u1u2 <- d2c_u1u2/cv - (dc_u1/cv)*(dc_u2/cv)

term1 <- d2logc_u1 * dF1^2
term2 <- d2logc_u2 * dF2^2
term3 <- 2 * d2logc_u1u2 * dF1 * dF2
term4 <- (dc_u1/cv) * d2F1
term5 <- (dc_u2/cv) * d2F2

cat("term1 (d2logc/du1^2 * dF1^2):", term1, "\n")
cat("term2 (d2logc/du2^2 * dF2^2):", term2, "\n")
cat("term3 (2*d2logc/(du1u2) * dF1 * dF2):", term3, "\n")
cat("term4 (dlogc/du1 * d2F1):", term4, "\n")
cat("term5 (dlogc/du2 * d2F2):", term5, "\n")
cat("Sum of all 5 terms:", term1+term2+term3+term4+term5, "\n")
cat("True:", true_d2, "\n")

u1 <- 0.5; u2 <- 0.6; par <- 0.5; fam <- 1
c_val <- BiCopPDF(u1, u2, fam, par)
dcu1 <- BiCopDeriv(u1, u2, fam, par, deriv = "u1")
d2cu1 <- BiCopDeriv2(u1, u2, fam, par, deriv = "u1")
d2logc_du1 <- d2cu1/c_val - (dcu1/c_val)^2
cat("d2 log c/du1^2:", d2logc_du1, "\n")
cat("d log c/du1 (=dcu1/c):", dcu1/c_val, "\n")

# Now: what is dF/dsigma and d2F/dsigma^2 for GG distribution?
# Typical GG parameters: mu=5, sigma=0.5, nu=1
mu <- 5; sigma <- 0.5; nu <- 1; y <- 4.5
h <- 1e-4
F0 <- pGG(y, mu=mu, sigma=sigma, nu=nu)
Fp <- pGG(y, mu=mu, sigma=sigma+h, nu=nu)
Fm <- pGG(y, mu=mu, sigma=sigma-h, nu=nu)
dF_sigma <- (Fp - Fm) / (2*h)
d2F_sigma <- (Fp + Fm - 2*F0) / h^2
cat("\nGG CDF derivatives w.r.t. sigma:\n")
cat("F(y):", F0, "\n")
cat("dF/dsigma:", dF_sigma, "\n")
cat("d2F/dsigma^2:", d2F_sigma, "\n")

# Compute diag_i1 term:
diag_i1_term1 <- d2logc_du1 * dF_sigma^2
diag_i1_term2 <- (dcu1/c_val) * d2F_sigma  # d log c / du1 * d2F/dsigma^2
diag_i1 <- diag_i1_term1 + diag_i1_term2
cat("\ndiag_i1 from d2 log c/du1^2 * dF^2:", diag_i1_term1, "\n")
cat("diag_i1 from d log c/du1 * d2F/ds^2:", diag_i1_term2, "\n")
cat("Total diag_i1:", diag_i1, "\n")

# What does the true d2 log c(F(y1,s), F(y2,s)) / ds^2 look like?
# Via FD on log copula density
y2 <- 4.8
logc_s <- function(s) {
  u1 <- pGG(y, mu=mu, sigma=s, nu=nu)
  u2 <- pGG(y2, mu=mu, sigma=s, nu=nu)
  log(BiCopPDF(u1, u2, fam, par))
}
d2logc_ds2_FD <- (logc_s(sigma+h) - 2*logc_s(sigma) + logc_s(sigma-h)) / h^2
cat("\nTrue d2 log c(F_i1, F_i2) / dsigma^2 (FD on logc(F(y1,s), F(y2,s))):", d2logc_ds2_FD, "\n")

# i2 contribution
F0_2 <- pGG(y2, mu=mu, sigma=sigma, nu=nu)
Fp_2 <- pGG(y2, mu=mu, sigma=sigma+h, nu=nu)
Fm_2 <- pGG(y2, mu=mu, sigma=sigma-h, nu=nu)
dF_sigma_2 <- (Fp_2 - Fm_2) / (2*h)
d2F_sigma_2 <- (Fp_2 + Fm_2 - 2*F0_2) / h^2
dcu2 <- BiCopDeriv(u1, F0_2, fam, par, deriv = "u2")
d2cu2 <- BiCopDeriv2(u1, F0_2, fam, par, deriv = "u2")
d2logc_du2 <- d2cu2/BiCopPDF(u1, F0_2, fam, par) - (dcu2/BiCopPDF(u1, F0_2, fam, par))^2
diag_i2 <- d2logc_du2 * dF_sigma_2^2 + (dcu2/BiCopPDF(u1, F0_2, fam, par)) * d2F_sigma_2
cat("diag_i2 (i2 self-term):", diag_i2, "\n")
cat("Total diag_i1 + diag_i2:", diag_i1 + diag_i2, "\n")

# Cross term (d2 log c / d sigma_i1 d sigma_i2)
u1v <- F0; u2v <- F0_2
cv <- BiCopPDF(u1v, u2v, fam, par)
dc_u1 <- BiCopDeriv(u1v, u2v, fam, par, deriv="u1")
dc_u2 <- BiCopDeriv(u1v, u2v, fam, par, deriv="u2")
hu <- 1e-5
d2c_u1u2 <- (BiCopPDF(u1v+hu, u2v+hu, fam, par) - BiCopPDF(u1v+hu, u2v-hu, fam, par) -
             BiCopPDF(u1v-hu, u2v+hu, fam, par) + BiCopPDF(u1v-hu, u2v-hu, fam, par)) / (4*hu^2)
cross <- d2c_u1u2 * dF_sigma * dF_sigma_2 / cv - (dc_u1 * dF_sigma / cv) * (dc_u2 * dF_sigma_2 / cv)
cat("\nCross pair contribution (d2logc/dsigma_i1 dsigma_i2):", cross, "\n")
cat("Expected from FD: use separate sigma perturbations\n")
logc_ss <- function(s1, s2) {
  u1 <- pGG(y, mu=mu, sigma=s1, nu=nu)
  u2 <- pGG(y2, mu=mu, sigma=s2, nu=nu)
  log(BiCopPDF(u1, u2, fam, par))
}
cross_FD <- (logc_ss(sigma+h, sigma+h) - logc_ss(sigma+h, sigma-h) -
             logc_ss(sigma-h, sigma+h) + logc_ss(sigma-h, sigma-h)) / (4*h^2)
cat("Cross FD:", cross_FD, "\n")
