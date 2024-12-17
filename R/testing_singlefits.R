source("R/common_functions.R");source("R/link_functions.R"); library("gamlss"); library("VineCopula");library("moments"); library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix")
#library(gamlss.longitudinal)

#########DATASET
n=1000; d=2

#copula_dist="C";margin_dist=GA(); mu=20; sigma=2;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NO(); mu=1; sigma=2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C";margin_dist=BE(); mu=0.4; sigma=0.6;nu=NA; tau=NA; theta=2; zeta=NA; simOption=7;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=PO(); mu=0.5; sigma=NA;nu=NA; tau=NA; theta=2; zeta=NA; simOption=7;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="N";margin_dist=PO(); mu=0.5; sigma=NA;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
#copula_dist="C"; margin_dist=ZISICHEL(); mu=3;sigma=exp(1);nu=-2;tau=.2;theta=5;zeta=NA;simOption=5;
#copula_dist="C"; margin_dist=ZISICHEL(); simOption
#copula_dist="C"; margin_dist=GA(); simOption=8; input_par=NA

copula_dist="N";margin_dist=NO(); mu=2; sigma=2;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=9; #MU TIME VARIANT
copula_dist="N";margin_dist=NO(); mu=2; sigma=2;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=10; #MU AND SIGMA TIME VARIANT


#copula_dist="C";margin_dist=GA(); mu=1; sigma=0.2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=7;
#covariates_input=list(mu.time=0,sigma.time=.2,nu.time=1,tau.time=1,theta.time=.5,zeta.time=1
#                      , mu.age=0,sigma.age=0,nu.age=1,tau.age=1,theta.age=0,zeta.age=1
#                      , mu.gender=0, sigma.gender=0, nu.gender=0,tau.gender=0,theta.gender=0, zeta.gender=0)
#    min_par=c(1,1,2);       max_par=c(10,10,10)


#simOption=6; margin_dist=JSU(); copula_dist="N"

mu_formula="response ~ (time)"
sigma_formula="~ 1"

#########Generate dataset

dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                    , par.margin=c(mu,sigma,nu,tau), par.copula=c(theta),covariates_input=NA)
plotDist(dataset,margin_dist)

##########FIT
#source("R/common_functions.R")
no_dl=fit_jointreg(dataset, margin_dist,copula_dist
                   , mu.formula = mu_formula
                   , sigma.formula = sigma_formula
                   , nu.formula = ("~ 1")
                   , tau.formula = ("~ 1")
                   , theta.formula=("~1")
                   , zeta.formula=("~1")
                   , include_dlcopdpar=FALSE
                   , verbose=3, plot_results=FALSE,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                   , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=0.01, outer_stop_crit=0.005
)

w_dl=fit_jointreg(dataset, margin_dist, copula_dist
                  , mu.formula = mu_formula
                  , sigma.formula = sigma_formula
                  , nu.formula = ("~ 1")
                  , tau.formula = ("~ 1")
                  , theta.formula=("~1")
                  , zeta.formula=("~1")
                  , include_dlcopdpar=TRUE
                  , verbose=3,plot_results=FALSE, true_val=par_to_eta(input_par,copula_dist,margin_dist)
                  , use_Rcpp=FALSE, start_step_size=.25, step_adjustment = 0.5, inner_stop_crit=0.01, outer_stop_crit=0.005
)


input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
#true_par=par_to_eta(input_par,margin_dist=margin_dist,copula_dist=copula_dist)
#names(true_par)=names((w_dl$par))

library(gee);
#gee_model=(gee(response~-1+as.factor(time),family=Gamma(link="log"),corstr="exchangeable",id=subject,data=dataset[order(dataset$subject,dataset$time),]))
gee_model=gee(mu_formula,family=gaussian
               ,corstr="AR-M",id=subject,data=dataset[order(dataset$subject,dataset$time),],Mv=1)

source("R/common_functions.R");
true_var_b0_bt=bvt_norm_true_SE_B0_Bt(sigma_x=sigma,sigma_y=sigma,rho=theta,n=n,d=d)
true_se_b0_bt=sqrt(true_var_b0_bt)/sqrt(n)

plot_true=TRUE
if(plot_true==TRUE) {
  results_table=cbind(true_par
                      ,unlist(coef(w_dl))
                      ,unlist(coef(no_dl))
                      ,c(gee_model$coefficients,log( sqrt(gee_model$scale)),NA,logit(gee_model$working.correlation[2]))
                      ,c(true_se_b0_bt,NA,NA,NA)
                      ,(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])
                      #,c(rep(sigma/sqrt(n),d),NA,NA)     #sqrt(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,(vcov(w_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,(vcov(no_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,c(sqrt(diag(gee_model$robust.variance)),0,0,0)
  )
  colnames(results_table)=c("True Est","Joint Est","Sep Est","GEE Est","True SE","ND True SE", "Joint SE", "Sep SE","GEE SE")
} else {
  results_table=cbind(unlist(coef(w_dl))
                      ,unlist(coef(no_dl))
                      ,c(gee_model$coefficients,log(sqrt( gee_model$scale)),logit(gee_model$working.correlation[2]))
                      ,sqrt(vcov(w_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,sqrt(vcov(no_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,c(sqrt(diag(gee_model$robust.variance)),0,0)
  )
  colnames(results_table)=c("Joint Est","Sep Est", "GEE Est","Joint SE","Sep SE","GEE SE")
}
print(round(results_table,4))
round(results_table,6)

