########### 0. Starting parameters and libraries
set.seed(1000)
source("R/common_functions.R"); source("R/link_functions.R")
library("gamlss"); library("VineCopula");library("moments"); library("parallel"); library("foreach"); library("doParallel")#library(gamlss.longitudinal)
options(scipen=999);set.seed(1000)

#1. Choose distribution and simulation parameters
n=100; d=2; sims=100 #simulation parameters
#copula_dist="C";margin_dist=EXP(); mu=10; sigma=NA;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=PO(); mu=10; sigma=NA;nu=NA; tau=NA; theta=2; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=GA(); mu=20; sigma=10;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C";margin_dist=NO(); mu=3; sigma=2;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NO(); mu=2; sigma=.2;nu=NA; tau=NA; theta=.25; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C";margin_dist=NBI(); mu=3; sigma=1;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NBI(); mu=3; sigma=1;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C"; margin_dist=ST1(); mu=5;sigma=exp(1);nu=10;tau=2;theta=5;zeta=NA; simOption=5;
#copula_dist="N"; margin_dist=ST1(); mu=1;sigma=exp(1);nu=10;tau=10;theta=5;zeta=NA; simOption=5
#copula_dist="C"; margin_dist=ZISICHEL(); mu=3;sigma=exp(1);nu=-2;tau=.2;theta=5;zeta=NA;simOption=5;
#copula_dist="N"; margin_dist=ZISICHEL(); simOption=1
#copula_dist="C"; margin_dist=GA(); simOption=8; input_par=NA; mu=sigma=nu=tau=theta=zeta=NA;
covariates_input=NA

###TIME VARIANT MU AND SIGMA (plus 1 per time in link form)
copula_dist="N";margin_dist=NO(); mu=2; sigma=2;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=9;

###TIME VARIANT MU AND SIGMA (plus 1 per time in link form)
#copula_dist="N";margin_dist=NO(); mu=2; sigma=2;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=10;

#copula_dist="C";margin_dist=GA(); mu=1; sigma=0.2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=7;
#covariates_input=list(mu.time=0,sigma.time=.2,nu.time=1,tau.time=1,theta.time=.5,zeta.time=1
#                      , mu.age=0,sigma.age=0,nu.age=1,tau.age=1,theta.age=0,zeta.age=1
#                      , mu.gender=0, sigma.gender=0, nu.gender=0,tau.gender=0,theta.gender=0, zeta.gender=0)
#    min_par=c(1,1,2);       max_par=c(10,10,10)


#simOption=6; margin_dist=JSU(); copula_dist="N"


############################# 2. CHOOSE MODEL FORMULAS

mu_formula="response ~ (time)"
sigma_formula="~ 1"

########################################## 3. SIMULATION LOOP - GENERATE, FIT, REPEAT ##########################################
input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
true_val=matrix(nrow=0,ncol=length(input_par[!is.na(input_par)]))
fit_gee=TRUE

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl, cores = detectCores() - 1)
out=foreach(run_counter=1:sims,.packages = c("VineCopula","gamlss","moments","gee"),.errorhandling = 'remove') %dopar% {

  plot_runs=FALSE;
  no_dl_outer=w_dl_outer=matrix(nrow=0,ncol=4+length(input_par[!is.na(input_par)]))

  #source("common_functions.R")

  #Generate dataset
  dataset=loadDataset(simOption=simOption,n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                      , par.margin=input_par[c("mu","sigma","nu","tau")], par.copula=c(theta=theta)
                      , covariates_input=covariates_input)

  #plotDist(dataset,margin_dist)

  #3. Fit model with and without joint component
  no_dl=fit_jointreg(dataset, margin_dist,copula_dist,
                     mu.formula = mu_formula, sigma.formula = sigma_formula, nu.formula = ("~ 1"), tau.formula = ("~ 1"),
                     theta.formula=("~1"), zeta.formula=("~1"),
                     include_dlcopdpar=FALSE,
                     verbose=3, plot_results=FALSE#,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                     , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=.01, outer_stop_crit=.005
  )
  w_dl=fit_jointreg(dataset, margin_dist, copula_dist,
                    mu.formula = mu_formula, sigma.formula = sigma_formula, nu.formula = ("~ 1"), tau.formula = ("~ 1"),
                    theta.formula=("~1"), zeta.formula=("~1"),
                    include_dlcopdpar=TRUE,
                    verbose=3,plot_results=FALSE#, true_val=par_to_eta(input_par,copula_dist,margin_dist)
                    , use_Rcpp=FALSE, start_step_size=.25, step_adjustment = 0.5, inner_stop_crit=.01, outer_stop_crit=.005
                    , start_from=no_dl$par
  )

  if(fit_gee==TRUE) {
    gee_fit=
    if(margin_dist$family[1]=="GA") {
      gee_fit=gee(formula=mu_formula,data=dataset[order(dataset$subject,dataset$time),],family=Gamma(link="log"),id=subject,corstr="AR-M",Mv=1)
    } else {
      gee_fit=gee(formula=mu_formula,data=dataset[order(dataset$subject,dataset$time),],family=gaussian,id=subject,corstr="AR-M",Mv=1)
    }
  }

  # PLOTTING
  #plot_runs=TRUE
  if(plot_runs==TRUE) {
    plot.new()
    num_plots=3+ncol(w_dl[[3]])
    #par(mfrow=c(round(sqrt(num_plots),0)+1,round(sqrt(num_plots),0)))
    par(mfrow=c(round((num_plots/3),0)+1,3))

    w_iter=nrow(w_dl[[2]])
    n_iter=nrow(no_dl[[2]])

    w_diff=round(w_dl[[2]][nrow(w_dl[[2]]),]-no_dl[[2]][nrow(no_dl[[2]]),],2)

    if(w_iter>n_iter) {
      for (type_plot in c("marginal","copula","joint")) {
        plot(1:nrow(w_dl[[2]]),w_dl[[2]][,type_plot],ylim=range(c(w_dl[[2]][,type_plot],no_dl[[2]][,type_plot])),type='l',col="red",main=paste(type_plot,w_diff[type_plot]),ylab="Likelihood",xlab="Iteration")
        lines(1:nrow(no_dl[[2]]),no_dl[[2]][,type_plot],col="black")
      }
    } else {
      for (type_plot in c("marginal","copula","joint")) {
        plot(1:nrow(no_dl[[2]]),no_dl[[2]][,type_plot],ylim=range(c(w_dl[[2]][,type_plot],no_dl[[2]][,type_plot])),type='l',col="black",main=paste(type_plot,w_diff[type_plot]),ylab="Likelihood",xlab="Iteration")
        lines(1:nrow(w_dl[[2]]),w_dl[[2]][,type_plot],col="red")
      }
    }

    par_names=colnames(w_dl$par_history)
    for (i in par_names) {
      if((all(exists("input_par")))) {
        true_val=par_to_eta(input_par,copula_dist,margin_dist)
        names(true_val)=par_names
        y_range=range(c(no_dl$par_history[,i],w_dl$par_history[,i],true_val[i]))
      } else {
        y_range=range(c(no_dl$par_history[,i],w_dl$par_history[,i]))
      }
      plot(w_dl$par_history[,i],type="l",main=i,xlab="Iteration",ylab='Parameter Value',ylim=y_range,col="red")
      lines(no_dl$par_history[,i],col="black")

      if((all(exists("input_par")))) {abline(h=true_val[i],col="blue")}
    }

  } # END PLOTTING

  list(no_dl,w_dl,gee_fit)
}
stopCluster(cl)

######################################### 4. PLOTTING ############################################
plot_true=FALSE
########TAKE OUT RESULTS
par_results_list=list()
#if(plot_true==TRUE) {
#  input=out[[1]][[1]]
#  true_par=par_to_eta(input_par,margin_dist=margin_dist,copula_dist=copula_dist)
#  names(true_par)=names((input$par))
#}

#true_par=c(rep(mu,d),log(sigma),logit(theta)) #Means
true_par=c(c(mu-1,1),log(sigma),logit(theta)) #Time parameter for mu
#true_par=c(c(mu-1,1),log(sigma)-1,0,logit(theta)) #Time parameter for mu and sigma

names(true_par)=names((coef(out[[1]][[2]])))

source("R/common_functions.R");
true_var_b0_bt=bvt_norm_true_SE_B0_Bt(sigma_x=sigma,sigma_y=sigma,rho=theta,n=n,d=d)
true_se_b0_bt=sqrt(true_var_b0_bt)/sqrt(n)

# Calculate variance / covariance for parameters to show
for (i in 1:length(out)) {

  #input=out[[i]][[1]]

  no_dl=out[[i]][[1]]
  w_dl=out[[i]][[2]]
  gee_model=out[[i]][[3]]
  #source("R/common_functions.R")
  #input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")

  if(plot_true==TRUE) {
    results_table=cbind(true_par
                        ,unlist(coef(w_dl))
                        ,unlist(coef(no_dl))
                        ,c(gee_model$coefficients,log( sqrt(gee_model$scale)),NA,logit(gee_model$working.correlation[2]))
                        ,c(true_se_b0_bt,NA,NA,NA)
                        ,(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])
                        ,(vcov(w_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                        ,(vcov(no_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                        ,c(sqrt(diag(gee_model$robust.variance)),0,0)
    )
    colnames(results_table)=c("True Est","Joint Est","Sep Est","GEE Est","True SE","ND True SE", "Joint SE", "Sep SE","GEE SE")
  } else {
    results_table=cbind(true_par
                        ,unlist(coef(w_dl))
                        ,unlist(coef(no_dl))
                        ,c(gee_model$coefficients,log( sqrt(gee_model$scale)),logit(gee_model$working.correlation[2]))
                        #,c(true_se_b0_bt,NA,NA,NA)
                        ,(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])
                        #,c(rep(sigma/sqrt(n),d),NA,NA)     #sqrt(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])*sqrt(n*d)
                        ,(vcov(w_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                        ,(vcov(no_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                        ,c(sqrt(diag(gee_model$robust.variance)),0,0)
    )
    colnames(results_table)=c("True Est","Joint Est","Sep Est","GEE Est","ND True SE", "Joint SE", "Sep SE","GEE SE")
  }
  print(round(results_table,4))
  round(results_table,6)

  for(par_temp in rownames(results_table)) {
    if(i==1) {
      par_results_list[[par_temp]]=results_table[par_temp,]
    } else {
      par_results_list[[par_temp]]=rbind(par_results_list[[par_temp]],results_table[par_temp,])
    }
  }
  print(i/length(out))

}



svg(file="results/charts/plot.svg",width=12, height=9)

# Plot means and standard errors for parameters
plot.new()
par(mfrow=c(2,length(par_results_list)))
for (i in 1:length(par_results_list)) {
  #if(grepl("mu",names(par_results_list)[i])) {
    boxplot(par_results_list[[i]][,c("Joint Est","Sep Est","GEE Est")],main=names(par_results_list)[i]);
  #} else {
    #boxplot(par_results_list[[i]][,c("Joint Est","Sep Est")],main=names(par_results_list)[i]);
  #}
  {abline(h=mean(par_results_list[[i]][,"True Est"]),col="red")}
}
for (i in 1:length(par_results_list)) {
  if(grepl("mu",names(par_results_list)[i])) {
    boxplot(par_results_list[[i]][,c("Joint SE","Sep SE","GEE SE")],main=paste(names(par_results_list)[i],"SE"))

  } else {
    boxplot(par_results_list[[i]][,c("Joint SE","Sep SE")],main=paste(names(par_results_list)[i],"SE"))
  }

  if(plot_true==TRUE){
    abline(h=mean(par_results_list[[i]][,"True SE"]),col="red")

    #legend("bottomright",legend=c("Exact","Numerical"),col=c("red","blue"),lty=c(1,2))
  }
  abline(h=mean(par_results_list[[i]][,"ND True SE"]),col="blue",lty=2)

}

dev.off()

#print(c(n=n, d=d, mu=mu, sigma=sigma,nu=nu, tau=tau, theta=theta, zeta=zeta))
