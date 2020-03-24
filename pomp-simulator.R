# doing a stochastic simulator with pomp
# no data fitting
# adopted from: https://kingaa.github.io/pomp/vignettes/He2010.html


## ----prelims,cache=FALSE-------------------------------------------------
set.seed(123L)
library(ggplot2)
library(dplyr)
library(pomp)

#######################################################################
# simple example from here
# https://kingaa.github.io/sbied/stochsim/stochsim.html
# with modifications
#######################################################################


read_csv("https://kingaa.github.io/sbied/stochsim/Measles_Consett_1948.csv") %>% select(week,reports=cases) -> meas
as.data.frame(meas)

seir_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")


covid_step <- Csnippet("
  double dN_SE = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_EI = rbinom(E,1-exp(-mu_EI*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  S -= dN_SE;
  E += dN_SE - dN_EI;
  I += dN_EI - dN_IR;
  R += dN_IR;
  H += dN_IR;
  
  
")




seir_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")


#r code equivalent
#sir_init <- function (N, eta, ...) {
#  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
#}


dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
  reports = rbinom(H,rho);
")

#create pomp object
measSEIR <- pomp(data = meas,
                times="week",t0=0,
                rprocess=euler(seir_step,delta.t=1/7),
                rinit=seir_init,
                rmeasure=rmeas,
                dmeasure=dmeas,
                accumvars="H",
                paramnames=c("N","Beta","mu_EI","mu_IR","rho","eta"),
                statenames=c("S","E","I","R","H")
              ) 

#run simulation
sims <- pomp::simulate(measSEIR,
                       params=c(Beta=40,mu_EI=0.8,mu_IR=1.3,rho=0.5,eta=0.06,N=38000),
                       nsim=20,format="data.frame",include.data=TRUE)

sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)




## ----rprocess------------------------------------------------------------
rproc <- Csnippet("
  double beta, br, seas, foi, dw, births;
  double rate[6], trans[6];
  
  // cohort effect
  if (fabs(t-floor(t)-251.0/365.0) < 0.5*dt) 
    br = cohort*birthrate/dt + (1-cohort)*birthrate;
  else 
  	br = (1.0-cohort)*birthrate;

  // term-time seasonality
  t = (t-floor(t))*365.25;
  if ((t>=7&&t<=100) || (t>=115&&t<=199) || (t>=252&&t<=300) || (t>=308&&t<=356))
      seas = 1.0+amplitude*0.2411/0.7589;
    else
      seas = 1.0-amplitude;

  // transmission rate
  beta = R0*seas*(1.0-exp(-(gamma+mu)*dt))/dt;
  // force of infection
  foi = beta*pow(I+iota,alpha)/pop;
  // white noise (extra-demographic stochasticity)
  dw = rgammawn(sigmaSE,dt);

  rate[0] = foi*dw/dt;  //infection rate (stochastic)
  rate[1] = mu;			    // natural S death
  rate[2] = sigma;		  // rate of ending of latent stage
  rate[3] = mu;			    // natural E death
  rate[4] = gamma;		  // recovery
  rate[5] = mu;			    // natural I death

  // Poisson births
  births = rpois(br*dt);
  
  // transitions between classes
  // more on this, see here: 
  //https://kingaa.github.io/pomp/FAQ.html#eulermultinomial-approximation
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
  C += trans[4];           // true incidence
")


## ----initializer---------------------------------------------------------
rinit <- Csnippet("
  double m = 100000/(S_0+E_0+I_0+R_0);
  S = nearbyint(m*S_0);
  E = nearbyint(m*E_0);
  I = nearbyint(m*I_0);
  W = 0;
  C = 0;
")

## ----create fake data to give to pomp---------------------------------------------------
dat = data.frame(time = 0:100, fakevar  = 100:200)

## ----pomp-construction---------------------------------------------------
dat %>%  pomp(t0=0,
    time="time",
    rprocess=euler(rproc,delta.t=1/20),
    rinit=rinit,
    accumvars=c("C","W"),
    statenames=c("S","E","I","C","W"),
    paramnames=c("R0","mu","sigma","gamma","alpha","iota",
      "rho","sigmaSE","psi","cohort","amplitude",
      "S_0","E_0","I_0","R_0")
  ) -> m1



## ----simulate model--------------------------------------------------
sim1 <- m1 %>% simulate(params=theta,nsim=9,format="d",include.data=TRUE) 


