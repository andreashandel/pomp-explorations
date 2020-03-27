#######################################################################
# doing a stochastic simulator with pomp
# no data fitting
# trying to do COVID
#######################################################################
#examples are here
#https://kingaa.github.io/pomp/vignettes/oaxaca.html#a_more_complex_example:_a_stochastic,_seasonal_sir_model
#https://kingaa.github.io/pomp/vignettes/He2010.html
#https://kingaa.github.io/sbied/stochsim/stochsim.html
#https://kingaa.github.io/sbied/index.html
#https://kingaa.github.io/sbied/stochsim/exercises.html#basic-exercise-the-seir-model


## ----prelims,cache=FALSE-------------------------------------------------
set.seed(123L)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pomp)


#######################################################################
#1 step for model
#R code for debugging/testing
#######################################################################
covid_step_R <- function(S, E1, E2, E3, E4, E5, E6, I1, I2, I3, I4, Iu1, Iu2, Iu3, Iu4, H, Ru, beta_d, beta_u, beta_e, beta_red_factor, t_int1, t_int2, t_int3, gamma_u, gamma_d, detect_frac_0,detect_frac_1,sigma,delta.t,...)
{
 
  Epresymptom = E1+E2+E3+E4+E5+E6;  #all pre-symptomatic
  Idetected = I1+I2+I3+I4;          #all symptomatic that wil be detected
  Iundetected = Iu1+Iu2+Iu3+Iu4;    #all symptomatic/asymptomatic that won't be detected

  dt = delta.t #using dt below
  
  trans = rep(0,17)

  detect_frac = detect_frac_0
  gamma = gamma_u
  foi = beta_d*Idetected + beta_u*Iundetected + beta_e*Epresymptom   
  
    
  #if (t_int1<t) {
  #     gamma = gamma_u
  #}
  #else {
  #     gamma = gamma_d
  #}
  # 
  # if (t<=t_int2) {
  #   detect_frac = detect_frac_0
  # }
  # else {
  #   detect_frac = detect_frac_1
  # }
  # 
  # if (t<=t_int3) {
  #       foi = beta_d*Idetected + beta_u*Iundetected + beta_e*Epresymptom   
  # }
  # else {
  #       foi = beta_red_factor*(beta_d*Idetected + beta_u*Iundetected + beta_e*Epresymptom)   
  # }

  
  trans[1] = rbinom(1,S,1-exp(-foi*dt));            
  trans[2] = rbinom(1,E1,1-exp(-sigma*dt));         
  trans[3] = rbinom(1,E2,1-exp(-sigma*dt));         
  trans[4] = rbinom(1,E3,1-exp(-sigma*dt));         
  trans[5] = rbinom(1,E4,1-exp(-sigma*dt));         
  trans[6] = rbinom(1,E5,1-exp(-sigma*dt));         
  trans[7] = rbinom(1,E6,(1-exp(-sigma*dt))*detect_frac);     
  trans[8] = rbinom(1,E6,(1-exp(-sigma*dt))*(1-detect_frac));
  trans[9] = rbinom(1,I1,1-exp(-gamma*dt));           
  trans[10] = rbinom(1,I2,1-exp(-gamma*dt));          
  trans[11] = rbinom(1,I3,1-exp(-gamma*dt));          
  trans[12] = rbinom(1,I4,1-exp(-gamma*dt));            
  trans[13] = rbinom(1,Iu1,1-exp(-gamma_u*dt));           
  trans[14] = rbinom(1,Iu2,1-exp(-gamma_u*dt));           
  trans[15] = rbinom(1,Iu3,1-exp(-gamma_u*dt));           
  trans[16] = rbinom(1,Iu4,1-exp(-gamma_u*dt));           

  S <- S - trans[1];
  E1 <- E1 + trans[1] - trans[2];
  E2 <- E2 + trans[2] - trans[3];
  E3 <- E3 + trans[3] - trans[4];
  E4 <- E4 + trans[4] - trans[5];
  E5 <- E5 + trans[5] - trans[6];
  E6 <- E6 + trans[6] - trans[7] - trans[8];
  I1 <- I1 + trans[7] - trans[9];
  I2 <- I2 +  trans[9] - trans[10];
  I3 <- I3 + trans[10] - trans[11];
  I4 <- I4 + trans[11] - trans[12];
  Iu1 <- Iu1 + trans[8] - trans[13];
  Iu2 <- Iu2 + trans[13] - trans[14];
  Iu3 <- Iu3 + trans[14] - trans[15];
  Iu4 <- Iu4 + trans[15] - trans[16];
  H <- H + trans[12]; 
  Ru <- Ru + trans[16]; 
  
  c(S = S, E1 = E1, E2 = E2, E3 =E3, E4=E4, E5=E5, E6=E6, I1=I1,I2=I2,I3=I3,I4=I4,Iu1=Iu1,Iu2=Iu2,Iu3=Iu3,Iu4=Iu4,H=H,Ru=Ru)
  
} #end 1-step function R code




#######################################################################
#initial conditions for model
#######################################################################
covid_init_R <- function (S_0, E1_0, E2_0, E3_0, E4_0, E5_0, E6_0, I1_0, I2_0, I3_0, I4_0, Iu1_0, Iu2_0, Iu3_0, Iu4_0, H_0, Ru_0,...) {
  c(S = S_0, E1 = E1_0, E2 = E2_0, E3 =E3_0, E4=E4_0, E5=E5_0, E6=E6_0, I1=I1_0,I2=I2_0,I3=I3_0,I4=I4_0,Iu1=Iu1_0,Iu2=Iu2_0,Iu3=Iu3_0,Iu4=Iu4_0,H=H_0,Ru=Ru_0)
}

#parameter names
#beta_d/u/e are transmission rates for those in the I (detected) Iu (undetected) and E (pre-symptomatic) classes
#beta_red_factor is reduction of transmission due to interventions
#t_int1 is day at which intervention kicks in that leads to faster diagnosis (gamma_u -> gamma_d)
#t_int2 is day at which intervention kicks in that leads to larger fraction diagnosed (detect_frac_0 -> detect_frac_1)
#t_int3 is day at which intervention kicks in that leads to reduction in transmssion by factor beta_red_factor
#gamma_u is rate of movement through Iu categories (inverse is duration of undetected infection). 
#gamma_d is rate of movement through I categories post intervention, i.e. shorter duration than Iu because of diagnosis and movement into H
#detect_frac_0 is the fraction detected before intervention
#detect_frac_1 is the fraction detected post intervention
#sigma is rate of movement through E categories (inverse is duration of latent period)
parnames1 = c("beta_d", "beta_u", "beta_e", "beta_red_factor", "t_int1", "t_int2", "t_int3", "gamma_u", "gamma_d", "detect_frac_0","detect_frac_1","sigma")

#variable names
varnames = c("S", "E1", "E2", "E3", "E4", "E5", "E6", "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4", "H", "Ru")

#initial conditions of state variables are also parameters
parnames2 = c("S_0", "E1_0", "E2_0", "E3_0", "E4_0", "E5_0", "E6_0", "I1_0", "I2_0", "I3_0", "I4_0", "Iu1_0", "Iu2_0", "Iu3_0", "Iu4_0", "H_0", "Ru_0")
parnames = c(parnames1,parnames2)


#######################################################################
#data loading
#note that we don't fit, but supply the data to the pomp object anyway since it wants data
#######################################################################
data <- read.csv('https://raw.githubusercontent.com/CEIDatUGA/COVID-19-DATA/master/GA_daily_status_report_GDPH.csv?token=ABQMQXYHNZXBU7W3UPEIJCS6QKTCE')
data$date <- as.Date(data$date, format='%m/%d/%y')

covid_ga_data <- data %>% dplyr::select(date, cases_cumulative, fatalities_cumulative) %>%
                          tidyr::replace_na(list(cases_cumulative = 0, fatalities_cumulative = 0)) %>%
                          dplyr::mutate(days = 1:nrow(data))


#create pomp object
covid_model_R <- pomp(data = covid_ga_data,
                      times="days",t0=0,
                      rprocess=euler(covid_step_R,delta.t=0.05),
                      rinit=covid_init_R,
                      statenames = varnames,
                      paramnames = parnames
                      ) 

inivals = c(S_0 = 10600000, E1_0 = 35, E2_0 = 35, E3_0 = 35, E4_0 = 35, E5_0 =35, E6_0 =35, I1_0 = 14, I2_0 = 14, I3_0 = 14, I4_0 = 14, Iu1_0 = 111, Iu2_0= 111, Iu3_0= 111, Iu4_0= 111, H_0  = 0, Ru_0 = 0 )
parvals = c(beta_d = 0.34/sum(inivals), beta_u = 0.5*0.34/sum(inivals), beta_e = 0.1*0.34/sum(inivals), beta_red_factor = 0.5, t_int1 = 12, t_int2 = 12, t_int3 = 12, gamma_u = 0.1, gamma_d = 0.2, detect_frac_0 = 0.1, detect_frac_1 = 0.9, sigma = 4*0.18, dt = 0.05)

#run simulation
sims <- pomp::simulate(covid_model_R,
                         params=c(parvals,inivals),
                         nsim=5, format="data.frame", include.data=FALSE)
  
pl <- sims %>%
  ggplot(aes(x=days,y=H,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

plot(pl)
