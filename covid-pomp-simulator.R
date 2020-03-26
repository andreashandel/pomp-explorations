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
library(pomp)

#######################################################################
#1 step for model
#######################################################################
covid_step <- Csnippet("
  double trans[17]; //C indexes at 0, I hate that so I'm making things 1 bigger and start with index 1, i.e. leave trans[0] empty
  double Idetected;
  double Iundetected;
  double Itot;
  double N;

  Idetected = I1+I2+I3+I4
  Iundetected = Iu1+Iu2+Iu3+Iu4
  Itot = Idetected + Iundetected  
  N = S + E1 + E2 + E3 + E4 + E5 + E6 + Itot + H + R;
  
  //time-dependent rate of movement through infected and detected classes
  //z is the date at which the intervention starts - it shouldn't really be less than 0
  //z days after simulation start, the time at which individuals are diagnosed and thus the time spent in the I categories before ending up in H decreases
  if (z<0)
    gamma = 1/(max(-0.57437*(t+12)+15.45651, a0)) 
  else
    if (t<=z) //if time is less then intervention time, the 
      gamma = b
    else
      gamma = a0

  //time dependent fraction of those that move into detected category at the end of the E phase
  //w days after simulation start, the fraction detected (those that move into I instead of Iu after exiting E) increases from q0 to q1
  //note that both higher fraction detected and faster rate of detection speed up arrival of individuals in H and subsequently C
  //this parameter is called q in the original code
  if (t<=w)
    detect_frac = q0
  else
    detect_frac = q1

  //time-dependent case notification rate, i.e. rate at which individuals in H compartment (confirmed cases) move into C (notified cases)
  //currently no change after time w is assumed
  //also note that this leads to an exponentially distributed process from H to C, which is unrealistic
  //this rate is called etat in the original code
  if (t<=w)
    asc_rate = 0.2
  else
    asc_rate = 0.2


  //time-dependent transmission 
  //w days after simulation start, an intervention reduces transmission rate by some factor
  if (t<=w)
    betat = beta0
  else
    betat = beta0/beta.factor
  
  //force of infection
  //time dependent transmission, multiplied by different groups
  //undetected are assumed to be less transmissible by factor c
  //if variable presymptomatic is 1, then the last stage of presymptomatic can also transmit
  foi = betat*Idetected/N+betat*c*Iundetected/N+presymptomatic*betat*c*E6/N;      // transmission


  //need to multiply by 6 since we have 6 latent stages
  sigma = 6*sigma
  
  //gamma_u is called b in the original code
  
  
  // define all transmission rates
  trans[1] = rbinom(S,1-exp(-foi*dt));              //transition from S to E
  trans[2] = rbinom(E1,1-exp(-sigma*dt));           // transition between E compartments 1/2
  trans[3] = rbinom(E2,1-exp(-sigma*dt));           // transition between E compartments
  trans[4] = rbinom(E3,1-exp(-sigma*dt));           // transition between E compartments
  trans[5] = rbinom(E4,1-exp(-sigma*dt));           // transition between E compartments 4/5
  trans[6] = rbinom(E5,1-exp(-sigma*dt));           // transition between E compartments 5/6
  
  trans[7] = rbinom(E6,(1-exp(-sigma*dt))*detect_frac);           // transition between E6 compartment and I 
  trans[8] = rbinom(E6,(1-exp(-sigma*dt))*(1-detect_frac));           // transition between E6 compartment and I 
  
  trans[9] = rbinom(I1,1-exp(-gamma*dt));           // transition between I compartments 1/2
  trans[10] = rbinom(I2,1-exp(-gamma*dt));           // transition between I compartments 2/3
  trans[11] = rbinom(I3,1-exp(-gamma*dt));           // transition between I compartments 3/4
  trans[12] = rbinom(I4,1-exp(-gamma*dt));          // transition between I compartments and H
  
  trans[13] = rbinom(Iu1,1-exp(-gamma_u*dt));           // transition between Iu compartments 1/2
  trans[14] = rbinom(Iu2,1-exp(-gamma_u*dt));           // transition between Iu compartments 2/3
  trans[15] = rbinom(Iu3,1-exp(-gamma_u*dt));           // transition between Iu compartments 3/4
  trans[16] = rbinom(Iu4,1-exp(-gamma_u*dt));           // transition between Iu compartments and Ru
  
  trans[17] = rbinom(H,1-exp(-asc_rate*dt));        // transition between H compartments and C
  
  // define all transmissions for each compartment
  S += - trans[1];
  
  E1 += trans[1] - trans[2];
  E2 += trans[2] - trans[3];
  E3 += trans[3] - trans[4];
  E4 += trans[4] - trans[5];
  E5 += trans[5] - trans[6];
  E6 += trans[6] - trans[7] - trans[8];
  
  I1 += trans[7] - trans[9];
  I2 += trans[9] - trans[10];
  I3 += trans[10] - trans[11];
  I4 += trans[11] - trans[12];
  
  Iu1 += trans[8] - trans[13];
  Iu2 += trans[13] - trans[14];
  Iu3 += trans[14] - trans[15];
  Iu4 += trans[15] - trans[16];
  
  H += trans[12] - trans[17];
  Ru += trans[16];
  C += trans[17];
  
  
")

covid_init <- Csnippet("
  S = S;
  E1 = E1;
  E2 = E2;
  E3 = E3;
  E4 = E4;
  E5 = E5;
  E6 = E6;
  I1 = I1;
  I2 = I2;
  I3 = I3;
  I4 = I4;
  Iu1 = Iu1;
  Iu2 = Iu2;
  Iu3 = Iu3;
  Iu4 = Iu4;
  H = H;
  Ru = Ru;
  C = C;
")

#parameter names
#beta0 is base transmission rate
#sigma is rate of movement through E categories (inverse is duration of latent period)
#b is rate of movement through Iu categories (inverse is duration of undetected infection). Also called gamma_u in my code
#a0 is rate of movement through I categories post intervention, i.e. shorter duration than Iu because of diagnosis and movement into H
#z is day at which intervention kicks in that leads to faster diagnosis
#w is day at which intervention kicks in that leads to larger fraction diagnosed and reduction in transmission rate
#c is a factor by which transmission of Iu individuals is lowered compared to I individuals 
#presymptomatic is either 0 or 1 and indicates if some presymptomatic individuals transmit
#q0 is the fraction detected before intervention
#q1 is the fraction detected post intervention
parnames = c("beta0", "sigma", "z", "b", "a0", "w", "c", "presymptomatic","q0","q1")
#variable names
varnames = c("S", "E1", "E2", "E3", "E4", "E5", "E6", "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4", "H", "Ru", "C")

#create pomp object
covid_model <- pomp(data = fake,
                times="days",t0=0,
                rprocess=euler(covid_step,delta.t=0.05),
                rinit=covid_init,
                rmeasure=rmeas,
                dmeasure=dmeas,
                accumvars="H",
                paramnames = parnames,
                statenames = varnames  
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
                         params=parms,
                         nsim=25, format="data.frame", include.data=FALSE)
  
  #scenarios.output[[i]] <- evaluate.model(params=parms, init = init, nsims=25, nstep=NULL, start=startdate)
}


