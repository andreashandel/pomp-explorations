library(tidyverse)
library(pomp)

read_csv("https://kingaa.github.io/sbied/stochsim/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) -> meas
as.data.frame(meas)


sir_step_R <- function (S, I, R, H, N, Beta, mu_IR, delta.t, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  
  if (t>=70) 
  {
    dN_SI = 0;
  }
  else 
  {
    dN_SI = dN_SI;
  }
  
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

sir_init_R <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}

sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
  double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
  double amplitude;
  
  if (t>=7)
      dN_SI = 0;
    else
      dN_SI = dN_SI;
  
  
  S -= dN_SI;
  I += dN_SI - dN_IR;
  R += dN_IR;
  H += dN_IR;
")

sir_init <- Csnippet("
  S = nearbyint(eta*N);
  I = 1;
  R = nearbyint((1-eta)*N);
  H = 0;
")

dmeas <- Csnippet("
  lik = dbinom(reports,H,rho,give_log);
")

rmeas <- Csnippet("
  reports = rbinom(H,rho);
")


measSIR <- meas %>%
  pomp(times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
       rinit=sir_init,
       rmeasure=rmeas,
       dmeasure=dmeas,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho")
  ) 

measSIR_R <- meas %>%
  pomp(times="week",t0=0,
       rprocess=euler(sir_step_R,delta.t=1/7),
       rinit=sir_init_R,
       rmeasure=rmeas,
       dmeasure=dmeas,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho")
  ) 


sims <- measSIR %>%
  simulate(params=c(Beta=40,mu_IR=0.5,rho=0.5,eta=0.03,N=38000),
           nsim=20,format="data.frame",include.data=TRUE) 

sims_R <- measSIR %>%
  simulate(params=c(Beta=40,mu_IR=0.5,rho=0.5,eta=0.03,N=38000),
           nsim=20,format="data.frame",include.data=TRUE) 


pl <- sims_R %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)

plot(pl)



## ----seir-diagram,echo=FALSE,cache=FALSE,eval=FALSE----------------------
library(DiagrammeR)
DiagrammeR("digraph SEIR {
  graph [rankdir=TD, overlap=false, fontsize = 10]
  node[shape=egg, label='B'] b;
  subgraph {
    rank=same;
    node[shape=oval, label='S'] S;
    node[shape=oval, label='E'] E;
    node[shape=oval, label='I'] I;
    node[shape=oval, label='R'] R;
    S->E E->I I->R
  }
  node[shape=diamond, label='dead'] d;
  b->S
  {S E I R}->d
   }",type="grViz",engine="dot",height=300,width=800)

#' 
#' ![model diagram](./model_diagram.png)

