library(ggplot2)
library(tidyr)
library(pomp)


#make process model
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
  reulermultinom(2,S,&rate[0],dt,&trans[0]);
  reulermultinom(2,E,&rate[2],dt,&trans[2]);
  reulermultinom(2,I,&rate[4],dt,&trans[4]);

  S += births   - trans[0] - trans[1];
  E += trans[0] - trans[2] - trans[3];
  I += trans[2] - trans[4] - trans[5];
  W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
  C += trans[4];           // true incidence
")





#simple model without data
simulate(t0=0, times=1:20,
         params=c(r=1.2,K=200,sigma=0.1,N_0=50),
         rinit=function (N_0, ...) {
           c(N=N_0)
         },
         rprocess=discrete_time(
           function (N, r, K, sigma, ...) {
             eps <- rnorm(n=1,mean=0,sd=sigma)
             c(N=r*N*exp(1-N/K+eps))
           },
           delta.t=1
         )
) -> sim1

simulate(
  sim1,
  params=c(r=1.2,K=200,sigma=0.1,N_0=50,b=0.05),
  rmeasure=function (N, b, ...) {
    c(Y=rpois(n=1,lambda=b*N))
  }
) -> sim2

p <- parmat(coef(sim2),3)
p["sigma",] <- c(0.05,0.25,1)
colnames(p) <- LETTERS[1:3]


simulate(sim2,params=p,
         times=seq(0,3),
         nsim=500,format="data.frame") -> sims

ggplot(data=tidyr::separate(sims,.id,c("parset","rep")),
       aes(x=N,fill=parset,group=parset,color=parset))+
  geom_density(alpha=0.5)+
  # geom_histogram(aes(y=..density..),position="dodge")+
  facet_grid(time~.,labeller=label_both,scales="free_y")+
  lims(x=c(NA,1000))



#pomp object containing data and simple model 

parus %>%
  pomp(
    times="year", t0=1960,
    rinit=function (N_0, ...) {
      c(N=N_0)
    },
    rprocess=discrete_time(
      function (N, r, K, sigma, ...) {
        eps <- rnorm(n=1,mean=0,sd=sigma)
        c(N=r*N*exp(1-N/K+eps))
      },
      delta.t=1
    ),
    rmeasure=function (N, b, ...) {
      c(pop=rpois(n=1,lambda=b*N))
    }
  ) -> rick

# continuous time stochastic SDE

#model
vpstep <- function (N, r, K, sigma, delta.t, ...) {
  dW <- rnorm(n=1,mean=0,sd=sqrt(delta.t))
  c(N = N + r*N*(1-N/K)*delta.t + sigma*N*dW)
}

# full pomp object
rick %>% pomp(rprocess=euler(vpstep,delta.t=1/365)) -> vp
