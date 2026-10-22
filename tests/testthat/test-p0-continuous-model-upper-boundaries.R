test_that("full Hessians cover upper tails and upper probability saturation", {
  old <- options(gamlss.longitudinal.copula_backend="native")
  on.exit(options(old),add=TRUE)
  for(backend in c("native","vinecopula"))
  for(family in c("N","C","F","G","J","t"))
  for(residual in c(5.5,7.5)) {
    options(gamlss.longitudinal.copula_backend=backend)
    object<-curvature_model_fixture(family)
    object$reference_backend<-backend
    object$response[2L]<-0.3+0.7*object$dataset$x[2L]+exp(-0.2)*residual
    object$dataset$response<-object$response
    object$dataset_original<-object$dataset
    boundary_followup_compare(object,paste("upper",backend,family,residual,sep="/"))
  }
})
