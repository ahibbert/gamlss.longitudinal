source("R/common_functions.R");source("R/link_functions.R"); library("VineCopula");library("moments"); 
library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix"); library("MASS")
#library(gamlss.longitudinal); 
library(gamlss2); library(gamlss)
set.seed(100)
#########DATASET
n=500; d=4

#copula_dist="N"; margin_dist=NO(); mu=0; sigma=1;nu=NA; tau=NA; theta=-1; zeta=NA; simOption=7;
#copula_dist="C";margin_dist=GA(); mu=1; sigma=0.2;nu=NA; tau=NA; theta=.5; zeta=NA; simOption=7;
#copula_dist="C"; margin_dist=PE(); mu=0.4; sigma=0.6;nu=1; tau=NA; theta=-0.5; zeta=NA; simOption=7;
#copula_dist="N"; margin_dist=ST1(); mu=1; sigma=1;nu=1; tau=1; theta=-0.5; zeta=NA; simOption=7;

copula_dist="N"; margin_dist=BCPEo(); mu=1; sigma=.1;nu=1; tau=1; theta=0; zeta=NA; simOption=10;

# USE THIS WITH SIMOPTION 10
covariates_input=list( mu.time=1   ,sigma.time=1   ,nu.time=1    ,tau.time=1   ,theta.time=.1  ,zeta.time=0
                        ,mu.age=5    ,sigma.age=5     ,nu.age=0     ,tau.age=0    ,theta.age=.5    ,zeta.age=0
                        ,mu.gender=1 ,sigma.gender=1  ,nu.gender=0  ,tau.gender=0 ,theta.gender=.25 ,zeta.gender=0)

#########Generate dataset

rm(dataset)
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

########## FIT ###########
source("R/common_functions.R")
mu_formula="random_name ~ time_of_observation_random_name + s(age_new_name,bs='ps') + gender"
sigma_formula="~ time_of_observation_random_name + gender + s(age_new_name,bs='ps')"
nu_formula="~ time_of_observation_random_name"
tau_formula="~ time_of_observation_random_name"
theta_formula="~ time_of_observation_random_name + s(age_new_name,bs='ps')"
zeta_formula="~ time_of_observation_random_name"

### FOR TESTING DATASETS WITH NON STANDARD NAMING
if("age_group" %in% names(dataset)) dataset$age_group = NULL
colnames(dataset)=c("person","time_of_observation_random_name","random_name","age_new_name","year","gender")
data_in=dataset; data_in$gender=as.factor(data_in$gender)
rm(dataset)

fit=gamlss.longitudinal(dataset=data_in
                   , margin_dist=margin_dist
                   , copula_dist=copula_dist
                   , time_var="time_of_observation_random_name"
                   , subject_var="person"
                   , mu.formula = mu_formula
                   , sigma.formula = sigma_formula
                   , nu.formula = nu_formula
                   , tau.formula = tau_formula
                   , theta.formula=theta_formula
                   , zeta.formula=zeta_formula
                   , include_dlcopdpar=FALSE
                   , verbose=1, plot_results=FALSE,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                   , use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=.1, outer_stop_crit=.1
                   , lambda_start = 5
                   , lambda_penalty_K = 2 #Optimising for AIC
)

#vcov_fit=vcov.gamlss.longitudinal(fit, numderiv=TRUE)

#################### PLOT METHOD ####################
source("R/common_functions.R")
summary(fit)
plot(
  fit,
  data = data_in,
  ci_level=0.90,
  max_plots_per_page=9,
  ncol=3,
  include_intercept=FALSE
)