.copula_dvine <- function(order, family, par, par2 = 0) {
  .copula_require_vinecopula("D-vine simulation setup via .copula_dvine()")

  VineCopula::D2RVine(order, .copula_family_numbers(family), par, par2)
}

.copula_rvine_sim <- function(n, rvm) {
  .copula_require_vinecopula("R-vine simulation via .copula_rvine_sim()")

  VineCopula::RVineSim(n, rvm)
}
