# gamlss.longitudinal: A simple to use, easy to inrepret  and highly flexible framework for regression of longitudinal data
  
## Overview of the approach
The purpose of the gamlss.longitudinal package is to provide a simple to use and easily interpretable modelling approach for longitudinal data. It allows users to fit longitudinal datasets with any gamlss marginal distributions with parameters relying on any parametric or smooth covariates, and simultaneously model the correlation structure of the data using copulas. Any copulas available in the VineCopula package can be used, and, unique to this longitudinal modelling approach, the correlation parameter can also be made to rely on any number of parametric or smooth covariates.

In the provided example code below, we show an example for a dataset with four time points, where the marginal distribution is the ST1 which is a 4-parameter gamlss marginal distribution with each parameter relying on a time-effect covariate and a smooth age effect. In addition, the example fits a copula which has increasing correlation with time and a fitted smooth age effect.

References:
Motivation for the approach and its performance compared to alternative methods in the bivariate case is available in the work in progress paper here: *'A comparison between copula-based, mixed model, and estimating equation methods for regression of bivariate correlated data'* Link: https://arxiv.org/abs/2410.11892.

## Installation

The gamlss.longitudinal package relies on the latest version of gamlss2 which is not yet available on CRAN, so to install gamlss.longitudinal, you must also install gamlss2 from source as below:

```R
# Install gamlss2 and gamlss.longitudinal
install.packages("gamlss2",
                 repos = c("https://gamlss-dev.R-universe.dev",
                           "https://cloud.R-project.org"))
devtools::install_github("ahibbert/gamlss.longitudinal")
```
  


## Example usage with a simulated dataset 

Example code for using gamlss.longitudinal to fit a longitudinal GAMLSS model

  

**R code:**

```R

  

# Install devtools if not already installed

if(!require(devtools)) install.packages("devtools")

  
  

# Install gamlss2 and gamlss.longitudinal

install.packages("gamlss2",

repos = c("https://gamlss-dev.R-universe.dev",

"https://cloud.R-project.org"))

devtools::install_github("ahibbert/gamlss.longitudinal")

  

# Load libraries & set seed

library(gamlss2)

library(gamlss.longitudinal)

set.seed(100)

  

# Set parameters

n=1000; d=4

copula_dist="N"; margin_dist=NO(); mu=0; sigma=1;nu=NA; tau=NA; theta=-1; zeta=NA; simOption=7;

#copula_dist="C";margin_dist=GA(); mu=1; sigma=0.2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=7;

#copula_dist="C"; margin_dist=PE(); mu=0.4; sigma=0.6;nu=1; tau=NA; theta=-0.5; zeta=NA; simOption=7;

#copula_dist="N"; margin_dist=ST1(); mu=1; sigma=1;nu=1; tau=1; theta=-0.5; zeta=NA; simOption=7;

covariates_input=list( mu.time=.1 ,sigma.time=.1 ,nu.time=.1 ,tau.time=.1 ,theta.time=.1 ,zeta.time=0

,mu.age=10 ,sigma.age=10 ,nu.age=10 ,tau.age=10 ,theta.age=10 ,zeta.age=0

,mu.gender=0 ,sigma.gender=0 ,nu.gender=0 ,tau.gender=0 ,theta.gender=0 ,zeta.gender=0)

  

# Generate and plot dataset

  

dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist

, par.margin=c(mu,sigma,nu,tau), par.copula=c(theta=theta),covariates_input=covariates_input)

plotDist(dataset,margin_dist)

  

fit=gamlss.longitudinal(dataset

, margin_dist

, copula_dist

, mu.formula="response ~ time + s(age,bs=\"ps\",k=10)"

, sigma.formula="time + s(age,bs=\"ps\")"

, nu.formula="time + s(age,bs=\"ps\")"

, tau.formula="time + s(age,bs=\"ps\")"

, theta.formula="time + s(age,bs=\"ps\")"

, zeta.formula="time"

, plot_results = FALSE

, verbose = 3

)

  

# Plot smoothers manually for now

par(mfrow=c(3,2))

plot(dataset$age,fit$model_matrix$s$mu$'s(age)'%*%fit$par_s$mu$'s(age)',main="mu age smooth",ylab="coefficient",xlab="age")

plot(dataset$age,fit$model_matrix$s$sigma$'s(age)'%*%fit$par_s$sigma$'s(age)', main="sigma age smooth",ylab="coefficient",xlab="age")

plot(dataset$age[1:(n*(d-1))],fit$model_matrix$s$theta$'s(age)'%*%fit$par_s$theta$'s(age)', main="theta age smooth",ylab="coefficient",xlab="age")

plot(dataset$age,fit$model_matrix$s$nu$'s(age)'%*%fit$par_s$nu$'s(age)', main="nu age smooth",ylab="coefficient",xlab="age")

plot(dataset$age,fit$model_matrix$s$tau$'s(age)'%*%fit$par_s$tau$'s(age)', main="tau age smooth",ylab="coefficient",xlab="age")

  

```