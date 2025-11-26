source("R/common_functions.R");source("R/link_functions.R"); library("VineCopula");library("moments"); 
library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix"); library("MASS")
#library(gamlss.longitudinal); 
library(gamlss2); library(gamlss)
set.seed(100)
#########DATASET
n=1000; d=4

#copula_dist="N"; margin_dist=NO(); mu=0; sigma=1;nu=NA; tau=NA; theta=-1; zeta=NA; simOption=7;
#copula_dist="C";margin_dist=GA(); mu=1; sigma=0.2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=7;
#copula_dist="C"; margin_dist=PE(); mu=0.4; sigma=0.6;nu=1; tau=NA; theta=-0.5; zeta=NA; simOption=7;
copula_dist="N"; margin_dist=ST1(); mu=1; sigma=1;nu=1; tau=1; theta=-0.5; zeta=NA; simOption=7;
# USE THIS WITH SIMOPTION 7
covariates_input=list( mu.time=.1   ,sigma.time=.1   ,nu.time=.1    ,tau.time=.1   ,theta.time=.1  ,zeta.time=0
                        ,mu.age=10    ,sigma.age=10     ,nu.age=10     ,tau.age=10    ,theta.age=10    ,zeta.age=0
                        ,mu.gender=0 ,sigma.gender=0  ,nu.gender=0  ,tau.gender=0 ,theta.gender=0 ,zeta.gender=0)

#simOption=6; margin_dist=JSU(); copula_dist="N"

#########Generate dataset

dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                    , par.margin=c(mu,sigma,nu,tau), par.copula=c(theta=theta),covariates_input=covariates_input)
plotDist(dataset,margin_dist)

#########PLOTTING#############
#Group dataset by age categories in buckets of 10 years then calculate the correlation in each bucket and plot
library(dplyr)
# Create age groups with extended range and filter out any remaining NAs
age_range = range(dataset$age, na.rm=TRUE)
age_breaks = seq(floor(age_range[1]/10)*10, ceiling(age_range[2]/10)*10 + 10, by=10)
dataset=dataset%>%
  mutate(age_group=cut(age, breaks=age_breaks, include.lowest=TRUE)) %>%
  filter(!is.na(age_group))

# Calculate correlations between consecutive time points
time_points = sort(unique(dataset$time))
age_corrs_list = list()

for(i in 1:(length(time_points)-1)) {
  t1 = time_points[i]
  t2 = time_points[i+1]
  temp_corrs = dataset %>% 
    group_by(age_group) %>% 
    summarise(cor = cor(response[time==t1], response[time==t2]), .groups='drop') %>%
    mutate(time_pair = paste0("T", t1, " vs T", t2)) %>%
    filter(!is.na(cor))
  age_corrs_list[[i]] = temp_corrs
}

age_corrs = bind_rows(age_corrs_list) %>% filter(!is.na(age_group))

p1 = ggplot(age_corrs, aes(x=age_group, y=cor, color=time_pair, group=time_pair)) + 
  geom_point() + 
  geom_line() + 
  ylim(0,1) +
  labs(title="Correlation by Age Group", 
       x="Age Group", 
       y="Correlation",
       color="Time Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate mean response by age group and time
age_response = dataset %>%
  group_by(age_group, time) %>%
  summarise(mean_response = mean(response), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p2 = ggplot(age_response, aes(x=age_group, y=mean_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Mean Response by Age Group",
       x="Age Group",
       y="Mean Response",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate standard deviation by age group and time
age_sd = dataset %>%
  group_by(age_group, time) %>%
  summarise(sd_response = sd(response), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p3 = ggplot(age_sd, aes(x=age_group, y=sd_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="SD of Response by Age Group",
       x="Age Group",
       y="Standard Deviation",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate kurtosis by age group and time
library(moments)
age_kurtosis = dataset %>%
  group_by(age_group, time) %>%
  summarise(kurtosis_response = kurtosis(response), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p4 = ggplot(age_kurtosis, aes(x=age_group, y=kurtosis_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Kurtosis of Response by Age Group",
       x="Age Group",
       y="Kurtosis",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate skewness by age group and time
age_skewness = dataset %>%
  group_by(age_group, time) %>%
  summarise(skewness_response = skewness(response), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p5 = ggplot(age_skewness, aes(x=age_group, y=skewness_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Skewness of Response by Age Group",
       x="Age Group",
       y="Skewness",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all five plots
library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3, common.legend=FALSE)

##########FIT
source("R/common_functions.R")
mu_formula="response ~ time + s(age,bs=\"ps\",k=10)"
sigma_formula="time + s(age,bs=\"ps\")"
nu_formula="time + s(age,bs=\"ps\")"
tau_formula="time + s(age,bs=\"ps\")"
theta_formula="time + s(age,bs=\"ps\")"
zeta_formula="time"

no_dl=gamlss.longitudinal(dataset, margin_dist,copula_dist
                   , mu.formula = mu_formula
                   , sigma.formula = sigma_formula
                   , nu.formula = nu_formula
                   , tau.formula = tau_formula
                   , theta.formula=theta_formula
                   , zeta.formula=zeta_formula
                   , include_dlcopdpar=FALSE
                   , verbose=3, plot_results=FALSE,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                   , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=.1, outer_stop_crit=.1
                   , lambda_start = 5
                   , lambda_penalty_K = 2 #Optimising for AIC
)
par(mfrow=c(3,2))
plot(dataset$age,no_dl$model_matrix$s$mu$'s(age)'%*%no_dl$par_s$mu$'s(age)',main="mu age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,no_dl$model_matrix$s$sigma$'s(age)'%*%no_dl$par_s$sigma$'s(age)', main="sigma age smooth",ylab="coefficient",xlab="age")
plot(dataset$age[1:(n*(d-1))],no_dl$model_matrix$s$theta$'s(age)'%*%no_dl$par_s$theta$'s(age)', main="theta age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,no_dl$model_matrix$s$nu$'s(age)'%*%no_dl$par_s$nu$'s(age)', main="nu age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,no_dl$model_matrix$s$tau$'s(age)'%*%no_dl$par_s$tau$'s(age)', main="tau age smooth",ylab="coefficient",xlab="age")


# gamlss no random effect
gamlss(response ~ time + ps(age), sigma.formula = ~ time + ps(age), nu.formula = ~ time + ps(age), data=dataset, family=margin_dist)


w_dl=gamlss.longitudinal(dataset, margin_dist, copula_dist
                  , mu.formula = mu_formula
                  , sigma.formula = sigma_formula
                  , nu.formula = nu_formula
                  , tau.formula = tau_formula
                  , theta.formula=theta_formula
                  , zeta.formula=zeta_formula
                  , include_dlcopdpar=TRUE
                  , verbose=3,plot_results=FALSE, true_val=par_to_eta(input_par,copula_dist,margin_dist)
                  , use_Rcpp=FALSE, start_step_size=1, step_adjustment = 0.5, inner_stop_crit=.1, outer_stop_crit=.1
                  , lambda_start = 5
                  , lambda_penalty_K = 2 #Optimising for AIC
)

par(mfrow=c(3,2))
plot(dataset$age,w_dl$model_matrix$s$mu$'s(age)'%*%w_dl$par_s$mu$'s(age)',main="mu age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,w_dl$model_matrix$s$sigma$'s(age)'%*%w_dl$par_s$sigma$'s(age)', main="sigma age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,w_dl$model_matrix$s$nu$'s(age)'%*%w_dl$par_s$nu$'s(age)', main="nu age smooth",ylab="coefficient",xlab="age")
plot(dataset$age,w_dl$model_matrix$s$tau$'s(age)'%*%w_dl$par_s$tau$'s(age)', main="tau age smooth",ylab="coefficient",xlab="age")
plot(dataset$age[1:(n*(d-1))],w_dl$model_matrix$s$theta$'s(age)'%*%w_dl$par_s$theta$'s(age)', main="theta age smooth",ylab="coefficient",xlab="age")


input_par=c(mu,sigma,nu,tau,theta,zeta); names(input_par)=c("mu","sigma","nu","tau","theta","zeta")
true_par=par_to_eta(input_par,margin_dist=margin_dist,copula_dist=copula_dist)
names(true_par)=names((w_dl$par))

library(gee);
#gee_model=(gee(response~-1+as.factor(time),family=Gamma(link="log"),corstr="exchangeable",id=subject,data=dataset[order(dataset$subject,dataset$time),]))
gee_model=gee(mu_formula,family=gaussian
               ,corstr="AR-M",id=subject,data=dataset[order(dataset$subject,dataset$time),],Mv=1)

source("R/common_functions.R");
true_var_b0_bt=bvt_norm_true_SE_B0_Bt(sigma_x=sigma,sigma_y=sigma,rho=theta,n=n,d=d)
true_se_b0_bt=sqrt(true_var_b0_bt)/sqrt(n)

plot_true=FALSE
if(plot_true==TRUE) {
  results_table=cbind(true_par
                      ,unlist(coef(w_dl))
                      ,unlist(coef(no_dl))
                      ,c(gee_model$coefficients,log( sqrt(gee_model$scale)),NA,logit(gee_model$working.correlation[2]))
                      ,c(true_se_b0_bt,NA,NA,NA)
                      #,(vcov(w_dl,par=true_par,numderiv = TRUE)[[2]])
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
                      ,(vcov(w_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,(vcov(no_dl,numderiv = TRUE)[[2]])#*sqrt(n*d)
                      ,c(sqrt(diag(gee_model$robust.variance)),0,0)
  )
  colnames(results_table)=c("Joint Est","Sep Est", "GEE Est","Joint SE","Sep SE","GEE SE")
}
print(round(results_table,4))
round(results_table,6)

