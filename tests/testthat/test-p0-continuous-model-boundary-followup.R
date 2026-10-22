test_that("included invalid pairs fail while excluded pairs remain neutral", {
  old <- options(gamlss.longitudinal.copula_backend="native")
  on.exit(options(old),add=TRUE)
  cache <- list(row_id1=1L,row_id2=2L,theta_index_map=c(1L,NA_integer_),
    observed_pair_base=TRUE)
  base <- list(mu=c(0,0),theta=0.4,margin_p_cache=c(0.2,0.4))
  for(backend in c("native","vinecopula")) {
    options(gamlss.longitudinal.copula_backend=backend)
    for(value in c(NA_real_,Inf,-0.1,1.1)) {
      eta <- base;eta$margin_p_cache[1L]<-value
      expect_error(.copula_hessian_pair_inputs(eta,cache,"N"),"Invalid.*included")
    }
    for(value in c(NA_real_,Inf)) {
      eta <- base;eta$theta<-value
      expect_error(.copula_hessian_pair_inputs(eta,cache,"N"),"Invalid.*included")
      eta <- base;eta$zeta<-value
      expect_error(.copula_hessian_pair_inputs(eta,cache,"t"),"Invalid.*included")
    }
    for(value in c(1e-12,1-1e-12)) {
      eta<-base;eta$margin_p_cache[1L]<-value
      expect_error(.copula_hessian_pair_inputs(eta,cache,"N"),"nondifferentiable copula probability")
    }
    eta<-base;eta$margin_p_cache[1L]<-1e-13
    input<-.copula_hessian_pair_inputs(eta,cache,"N")
    expect_identical(input$u1,1e-12)
    expect_false(input$u1_active)
    expect_true(input$u2_active)
    excluded<-cache;excluded$observed_pair_base<-FALSE
    eta<-base;eta$margin_p_cache[]<-NA_real_;eta$theta<-Inf
    input<-.copula_hessian_pair_inputs(eta,excluded,"N")
    expect_false(input$pair_ok)
    expect_identical(input$par1_e,0)
  }
})

test_that("included nonfinite Hessian arithmetic is never replaced by zero", {
  expect_error(.copula_hessian_zero_invalid(c(1, Inf), c(TRUE, TRUE)), "Nonfinite.*included")
  expect_equal(.copula_hessian_zero_invalid(c(1, NaN), c(TRUE, FALSE)), c(1, 0))
  expect_error(.copula_hessian_finite_or_zero(NaN), "Nonfinite.*included")
  derivs <- list(c_val=1, dcdu1=Inf, dcdu2=1, d2cdu1u2=1)
  expect_error(.copula_hessian_cross_pair_terms("mu",1L,2L,TRUE,
    list(mu=c(1,1)),derivs), "Nonfinite.*included")
})

test_that("full Hessians cover meaningful dependence parameter boundaries", {
  old<-options(gamlss.longitudinal.copula_backend="native")
  on.exit(options(old),add=TRUE)
  for(backend in c("native","vinecopula"))
  for(family in c("N","C","F","G","J","t")) {
    options(gamlss.longitudinal.copula_backend=backend)
    object<-curvature_model_fixture(family)
    object$reference_backend<-backend
    object$par[["theta.intercept"]]<-switch(family,N=atanh(-0.85),
      t=atanh(-0.5),C=log(0.05),F=-3,G=log(0.01),J=log(0.01))
    if(family=="t") object$par[["zeta.intercept"]]<-log(0.01)
    boundary_followup_compare(object,paste("parameter",backend,family,sep="/"))
  }
})

test_that("Gaussian CDF derivatives retain relative accuracy in both tails", {
  y<-c(-6.6,6.6)
  eta<-list(mu=c(0,0),sigma=c(1,1))
  mm<-list(x=list(mu=matrix(1,2,1),sigma=matrix(1,2,1)))
  margin<-gamlss.dist::NO()
  d1<-.calc_dFdpar(eta,mm,margin,y)
  d2<-.calc_d2Fdpar2(eta,mm,margin,y)
  cross<-.calc_d2Fdpar_cross(eta,mm,margin,y)
  phi<-dnorm(y)
  references<-list(d1$mu/(-phi),d1$sigma/(-y*phi),
    d2$mu/(-y*phi),d2$sigma/(y*(2-y^2)*phi),
    cross$mu$sigma/((1-y^2)*phi))
  for(ratio in references) expect_equal(ratio,c(1,1),tolerance=1e-4)
})

test_that("Gumbel link plateaus have zero parameter curvature and exact joins refuse", {
  old<-options(gamlss.longitudinal.copula_backend="native")
  on.exit(options(old),add=TRUE)
  for(backend in c("native","vinecopula")) {
    options(gamlss.longitudinal.copula_backend=backend)
    object<-curvature_model_fixture("G")
    object$reference_backend<-backend
    object$par[["theta.intercept"]]<-log(17)
    boundary_followup_compare(object,paste("gumbel-plateau",backend,sep="/"))
    object$par[["theta.intercept"]]<-log(16)
    expect_error(calc_analytical_hessian(object,progress=FALSE),"nondifferentiable Gumbel")
  }
})
