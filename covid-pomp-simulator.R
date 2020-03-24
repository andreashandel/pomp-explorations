#######################################################################
# doing a stochastic simulator with pomp
# no data fitting
# trying to do COVID
#######################################################################

## ----prelims,cache=FALSE-------------------------------------------------
set.seed(123L)
library(ggplot2)
library(dplyr)
library(pomp)

#######################################################################
#1 step for model
#######################################################################
covid_step <- Csnippet("
  double trans[10];
  double rate[6];
  double Idetected;
  double Iundetected;
  double Itot;

  Idetected = I1+I2+I3+I4
  Iundetected = Iu1+Iu2+Iu3+Iu4
  Itot = Idetected + Iundetected  
  
  N = S + E1 + E2 + E3 + E4 + E5 + E6 + Itot + R;
  
  //time-dependent transmission rate gamma
  if (z<0)
    gamma = 1/(max(-0.57437*(t+12)+15.45651, a0)) 
    else
    gamma <- ifelse(t<=z, gamma <- b, gamma <- a0);
  
  
  sigmai <- 6*sigma  # multiplier 2 for pseudo stages
  etat <- eta(t,w)     # case notification rate
  betat <- beta(t,w)   # time dependent transmissibility, presymptomatic=1 causes this transmissibility to apply to late stage latent cases as well
  
  rate[0] = betat*Idetected/N+betat*c*Iundetected/N+presymptomatic*betat*c*E6/N;      // transmission
  rate[1] = 6*sigma;           // transition between E compartments 1/2
  rate[2] = 6*sigma;           // transition between E compartments
  rate[3] = 6*sigma;           // transition between E compartments
  rate[4] = 6*sigma;           // transition between E compartments
  rate[5] = 6*sigma;           // transition between E compartments 5/6
  rate[6] = 6*sigma;           // transition between E compartment and I or Iu
  rate[7] = 4*gamma;           // transition between I compartments 1/2
  rate[8] = 4*gamma;           // transition between I compartments 2/3
  rate[9] = 4*gamma;           // transition between I compartments 3/4
  rate[10] = 4*gamma;          // transition between I compartments and H
  rate[11] = b;           // transition between Iu compartments 1/2
  rate[12] = b;           // transition between Iu compartments 2/3
  rate[13] = b;           // transition between Iu compartments 3/4
  rate[14] = b;           // transition between Iu compartments and Ru
  rate[15] = etat;        // transition between H compartments and C
  
  
  // transitions between classes
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  
  S += trans[0];
  E1 += ;
  E2 += trans[2] - trans[3];
  E3 += trans[3] - trans[4];
  E4 += trans[4] - trans[5];
  E5 += trans[5] - trans[6];
  E6 += dN_SE - dN_EI;
  I1 += dN_EI - dN_IR;
  I2 += dN_EI - dN_IR;
  I3 += dN_EI - dN_IR;
  I4 += dN_EI - dN_IR;
  Iu1 += dN_EI - dN_IR;
  Iu2 += dN_EI - dN_IR;
  Iu3 += dN_EI - dN_IR;
  Iu4 += dN_EI - dN_IR;
  H += dN_IR;
  Ru += dN_IR;
  C += dN_IR;
  
  
")




params <- c( )


covid_init <- Csnippet("
  S = nearbyint(eta*N);
  E = 0;
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

#create pomp object
covid_model <- pomp(data = fake,
                times="days",t0=0,
                rprocess=euler(covid_step,delta.t=0.05),
                rinit=covid_init,
                rmeasure=rmeas,
                dmeasure=dmeas,
                accumvars="H",
                paramnames = ,
                statenames = c("S", "E1", "E2", "E3", "E4", "E5", "E6", "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4", "H", "Ru", "C") 
              ) 




# load scenarios, run model for each scenario
scenarios <- read.csv('Georgia scenarios - Sheet1.csv')
start=as.Date("2020-03-01")

scenarios.output <- list()
for(i in 1:dim(scenarios)[1]){
  parms <- as.list(scenarios[i,3:13])
  init <- as.list(scenarios[i,14:31])
  startdate = as.Date("2020-03-01")

  #run simulation
  sims <- pomp::simulate(covid_model,
                         params=c(Beta=40,mu_EI=0.8,mu_IR=1.3,rho=0.5,eta=0.06,N=38000),
                         nsim=25, format="data.frame", include.data=FALSE)
  
  #scenarios.output[[i]] <- evaluate.model(params=parms, init = init, nsims=25, nstep=NULL, start=startdate)
}


