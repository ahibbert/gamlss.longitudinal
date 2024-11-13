source("R/common_functions.R");source("R/link_functions.R"); library("gamlss"); library("VineCopula");library("moments"); library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix")
#library(gamlss.longitudinal)

#########DATASET
n=100; d=2

copula_dist="C";margin_dist=GA(); mu=20; sigma=2;nu=NA; tau=NA; theta=5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="N";margin_dist=NO(); mu=1; sigma=2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=5;#    min_par=c(1,1,2);       max_par=c(10,10,10)
#copula_dist="C";margin_dist=BE(); mu=0.4; sigma=0.6;nu=NA; tau=NA; theta=2; zeta=NA; simOption=7;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="C";margin_dist=PO(); mu=0.5; sigma=NA;nu=NA; tau=NA; theta=2; zeta=NA; simOption=7;  #Note these can be times by 2 or half as start parameters so be careful
#copula_dist="N";margin_dist=PO(); mu=0.5; sigma=NA;nu=NA; tau=NA; theta=.75; zeta=NA; simOption=5;  #Note these can be times by 2 or half as start parameters so be careful
#input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
#copula_dist="C"; margin_dist=ZISICHEL(); mu=3;sigma=exp(1);nu=-2;tau=.2;theta=5;zeta=NA;simOption=5;
#copula_dist="C"; margin_dist=ZISICHEL(); simOption
#copula_dist="C"; margin_dist=GA(); simOption=8; input_par=NA

#########Generate dataset

dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                    , par.margin=c(mu,sigma,nu,tau), par.copula=c(theta),covariates_input=NA)
plotDist(dataset,margin_dist)

##########FIT
#source("R/common_functions.R")
no_dl=fit_jointreg(dataset, margin_dist,copula_dist
                   , mu.formula = ("response ~ -1+as.factor(time)")
                   , sigma.formula = ("~1")
                   , nu.formula = ("~ 1")
                   , tau.formula = ("~ 1")
                   , theta.formula=("~1")
                   , zeta.formula=("~1")
                   , include_dlcopdpar=FALSE
                   , verbose=3, plot_results=FALSE,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                   , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=0.02, outer_stop_crit=0.01
)

w_dl=fit_jointreg(dataset, margin_dist, copula_dist
                  , mu.formula = ("response ~ -1+as.factor(time)")
                  , sigma.formula = ("~ 1")
                  , nu.formula = ("~ 1")
                  , tau.formula = ("~ 1")
                  , theta.formula=("~1")
                  , zeta.formula=("~1")
                  , include_dlcopdpar=TRUE
                  , verbose=3,plot_results=FALSE, true_val=par_to_eta(input_par,copula_dist,margin_dist)
                  , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=0.02, outer_stop_crit=0.01
)


input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
#true_par=par_to_eta(input_par,margin_dist=margin_dist,copula_dist=copula_dist)
#names(true_par)=names((w_dl$par))

library(gee);
gee_model=(gee(response~-1+as.factor(time),family=Gamma(link="log"),corstr="exchangeable",id=subject,data=dataset[order(dataset$subject,dataset$time),]))
#gee_model=(gee(response~1,family=gaussian,corstr="exchangeable",id=subject,data=dataset[order(dataset$subject,dataset$time),]))

results_table=NA
plot_true=FALSE
if(plot_true==TRUE) {
  results_table=cbind(true_par
                      ,unlist(coef(w_dl))
                      ,unlist(coef(no_dl))
                      ,c(gee_model$coefficients,0,0)
                      ,sqrt(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,sqrt(vcov(w_dl,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,sqrt(vcov(no_dl,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,c(sqrt(gee_model$robust.variance),0,0)
  )
  colnames(results_table)=c("True Est","Joint Est","Sep Est","GEE Est","True SE", "Joint SE", "Sep SE","GEE SE")
} else {
  results_table=cbind(unlist(coef(w_dl))
                      ,unlist(coef(no_dl))
                      ,c(gee_model$coefficients,0,0)
                      ,sqrt(vcov(w_dl,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,sqrt(vcov(no_dl,numderiv = TRUE)[[2]])*sqrt(n*d)
                      ,c(sqrt(diag(gee_model$robust.variance)),0,0)
  )
  colnames(results_table)=c("Joint Est","Sep Est", "GEE Est","Joint SE","Sep SE","GEE SE")
  rownames(results_table)=names(w_dl$par)
}
print(round(results_table,4))


#######TESTING TO DELETE
source("R/common_functions.R")
vcov_out=vcov(w_dl,par=true_par,numderiv = TRUE)

sqrt(diag(vcov_out[[1]]))
sqrt(unlist(vcov_out[[2]]))*sqrt(n*d)


#####TRUE SE OVERALL FOR MU with LOG LINK
sqrt(vcov_out[[1]][1,1]/(margin_dist$mu.dr(results_table[1,"True Est"])^2)) #Variance covariance matrix without covariates

sqrt(vcov_out[[2]]$mu)/sqrt(n*d)



sqrt(vcov(w_dl,par=true_par,numderiv = TRUE)[[1]][1,1]
diag(vcov(w_dl,par=true_par,numderiv = TRUE))


vcov(w_dl,par=true_par,numderiv = TRUE)


source("R/common_functions.R")
vcov_true=vcov(w_dl,par=true_par,numderiv=TRUE)
vcov_joint=vcov(w_dl,par=w_dl$par,numderiv=TRUE)
vcov_sep=vcov(w_dl,par=no_dl$par,numderiv=TRUE)



vcov_est=(d2[[2]])

se_est=c(sigma
,d2[[1]][1] #Not adjusted for correlation
,gee_model$robust.variance #GEE
,vcov_est[1] #Numderiv with estimated paramters
,true_SE_in[1] #Numderiv with true parameters
)

names(se_est)=c("TRUE","","GEE","Numderiv w/ est","Numderiv w TRUE est")
se_est


sqrt(-solve(vcov(w_dl,par=w_dl$par,numderiv=TRUE)[[2]]))
solve(vcov(w_dl,par=w_dl$par,numderiv=TRUE)[[1]])
sqrt(-solve(vcov(no_dl,par=true_par,numderiv=TRUE)[[2]]))
sqrt(-solve(vcov(no_dl,par=true_par,numderiv=TRUE)[[1]]))

#sqrt(vcov(w_dl,par=w_dl$par,numderiv=TRUE)[[2]])
#sqrt(vcov(no_dl,par=no_dl$par,numderiv=TRUE)[[2]])
#sqrt(vcov(no_dl,par=true_par,numderiv=TRUE)[[2]])
#sqrt(vcov(no_dl,par=true_par,numderiv=TRUE)[[2]])
#sqrt(vcov(no_dl,par=true_par,numderiv=FALSE)[[1]])/sqrt(n*d)
#sqrt(vcov(no_dl,par=true_par,numderiv=FALSE)[[1]])/sqrt(n*d)





#SE for joint


true_SE=diag(true_SE_in[[1]])/(n*d)
sqrt(true_SE)

results_table=cbind(true_par
              ,unlist(coef(w_dl))
              ,unlist(coef(no_dl))
              ,true_SE
              ,diag(vcov(w_dl,sep_d2 = FALSE)/n*d)
              ,diag(vcov(w_dl,sep_d2 = TRUE)/n*d)
              ,diag(vcov(no_dl)/n*d)

              ,

              )
colnames(results_table)=c("True Est","w_dl","no_dl","True SE", "w_dl_SE (adj)","w_dl_SE","no_dl_SE")
round(results_table,6)



results_table=cbind(unlist(coef(w_dl))
                    ,unlist(coef(no_dl))
                    ,diag(vcov(w_dl,sep_d2 = FALSE)/n*d)
                    ,diag(vcov(w_dl,sep_d2 = TRUE)/n*d)
                    ,diag(vcov(no_dl)/n*d) )
colnames(results_table)=c("w_dl","no_dl", "w_dl_SE (adj)","w_dl_SE","no_dl_SE")
round(results_table,6)


plot_runs=TRUE
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
    if((all(exists("input_par")))&!all(is.na(input_par))) {
      true_val=par_to_eta(input_par,copula_dist,margin_dist)
      names(true_val)=par_names
      y_range=range(c(no_dl$par_history[,i],w_dl$par_history[,i],true_val[i]))
    } else {
      y_range=range(c(no_dl$par_history[,i],w_dl$par_history[,i]))
    }
    plot(w_dl$par_history[,i],type="l",main=i,xlab="Iteration",ylab='Parameter Value',ylim=y_range,col="red")
    lines(no_dl$par_history[,i],col="black")

    if((all(exists("input_par")))&!all(is.na(input_par))) {abline(h=true_val[i],col="blue")}
  }

} # END PLOTTING


norm_cov=matrix(c(sigma^2,theta*sigma^2,theta*sigma^2,sigma^2),2,2)
sqrt(norm_cov)/sqrt(n*d)


###DIFFERENCE DECREASES WITH SAMPLE SIZE???




