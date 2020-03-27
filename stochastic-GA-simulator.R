##putting all functions in separate R file

#main function to perform single simulation step
onestep <- function (x, params) {  #function to calculate one step of stochastic SIR
  
  S <- x[2]                            #local variable for susceptibles
  
  E1 <- x[3]                           #exposed classes
  E2 <- x[4]
  E3 <- x[5]
  E4 <- x[6]
  E5 <- x[7]
  E6 <- x[8]
  
  I1 <- x[9]                           #detected infectious classes
  I2 <- x[10]
  I3 <- x[11]
  I4 <- x[12]
  
  Iu1 <- x[13]                           #undetected infectious classes
  Iu2 <- x[14]
  Iu3 <- x[15]
  Iu4 <- x[16]
  
  I.detected <-I1+I2+I3+I4
  I.undetected <- Iu1+Iu2+Iu3+Iu4
  I <- I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4     #total infectious
  
  H <- x[17]                           #local variable for Isolated
  Ru <- x[18]                          #local variable for undetected recovereds
  
  C <- x[19]    # local variable for notifications
  
  N <- S+E1+E2+E3+E4+E5+E6+I1+I2+I3+I4+Iu1+Iu2+Iu3+Iu4+H+Ru          # total size of population
  
  t <- x[20]                          #get current time
  
  with(                                #use with to simplify code
    as.list(params), 
    {
      gammai <- 4*gamma(z=z, b=b, a0=a0, t=as.numeric(t))  # multiplier 4 for pseudo stages
      sigmai <- 6*sigma  # multiplier 2 for pseudo stages
      etat <- eta(t,w)     # case notification rate
      betat <- beta(t,w)   # time dependent transmissibility, presymptomatic=1 causes this transmissibility to apply to late stage latent cases as well
      
      rates <- as.numeric(c(betat*I.detected/N+betat*c*I.undetected/N+presymptomatic*betat*c*E6/N,                           # movements out of S
                            sigmai, sigmai, sigmai, sigmai, sigmai, sigmai,               # movements out of E (left dummy compartments with really high rates)
                            gammai, gammai, gammai, gammai,                   # movements out of I (detected)
                            4*b, 4*b, 4*b, 4*b,                                       # movements out of I (undetected)
                            etat))
      
      states0 <- x[2:(length(x)-1)]
      
      # transition probabilities
      
      p <- matrix(0, nrow=length(rates),ncol=length(states0))                                     # matrix to hold transitions probs
      
      p[1,]  <- c(exp(-rates[1]*dt),1-exp(-rates[1]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # S-> E
      
      p[2,]  <- c(0, exp(-rates[2]*dt), 1-exp(-rates[2]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[3,]  <- c(0, 0, exp(-rates[3]*dt), 1-exp(-rates[3]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[4,]  <- c(0, 0, 0, exp(-rates[4]*dt), 1-exp(-rates[4]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[5,]  <- c(0, 0, 0, 0, exp(-rates[5]*dt), 1-exp(-rates[5]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[6,]  <- c(0, 0, 0, 0, 0, exp(-rates[6]*dt), 1-exp(-rates[6]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of E
      p[7,]  <- c(0, 0, 0, 0, 0, 0, exp(-rates[7]*dt), (1-exp(-rates[7]*dt))*q(t,w,q0,q1), 0, 0, 0, (1-exp(-rates[7]*dt))*(1-q(t,w,q0,q1)), 0, 0, 0, 0, 0, 0)    # Transitions out of E
      
      p[8,]  <- c(0, 0, 0, 0, 0, 0, 0, exp(-rates[8]*dt), 1-exp(-rates[8]*dt), 0, 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I (detected)
      p[9,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[9]*dt), 1-exp(-rates[9]*dt), 0, 0, 0, 0, 0, 0, 0, 0)    # Transitions out of I
      p[10,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[10]*dt), 1-exp(-rates[10]*dt), 0, 0, 0, 0, 0, 0, 0)  # Transitions out of I
      p[11,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[11]*dt), 0, 0, 0, 0, 1-exp(-rates[11]*dt), 0, 0)  # Transitions out of I -> H
      
      p[12,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[12]*dt), 1-exp(-rates[12]*dt), 0, 0, 0, 0, 0) # Transitions out of I (undetected)
      p[13,]  <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[13]*dt), 1-exp(-rates[13]*dt), 0, 0, 0, 0) # Transitions out of I
      p[14,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[14]*dt), 1-exp(-rates[14]*dt), 0, 0, 0)  # Transitions out of I
      p[15,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  exp(-rates[15]*dt), 0, 1-exp(-rates[15]*dt), 0)  # Transitions out of I -> R_u
      
      
      p[16,] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, exp(-rates[16]*dt), 0, 1-exp(-rates[16]*dt))  # Transitions R_d -> to C (notification)
      
      # update states
      
      states1 <- matrix(0, nrow=length(rates),ncol=length(states0))                                # matrix to hold future states
      
      for(i in 1:length(rates)){
        states1[i,] <- t(rmultinom(1, states0[i], p[i,])) 
      }
      
      states1 <- colSums(states1)
      states1[17] <- states1[17]+Ru  #add formerly Recovered undetected cases
      states1[18] <- states1[18]+C  #add formerly notified cases
      
      return(x <- c(dt, states1, tail(x,1)+dt))
    }
  )
}

#model function, calls the main onestep function multiple times
model <- function (x, params, nstep) {  #function to simulate stochastic SIR
  output <- array(dim=c(nstep+1,length(x)))         #set up array to store results
  colnames(output) <- c("time","S",
                        "E1", "E2", "E3", "E4", "E5", "E6",
                        "I1", "I2", "I3", "I4", "Iu1", "Iu2", "Iu3", "Iu4",
                        "H", "Ru", "C", "cum.time") #name variables
  output[1,] <- x                           #first record of output is initial condition
  for (k in 1:nstep) {                      #iterate for nstep steps
    output[k+1,] <- x <- as.numeric(onestep(x,params))
  }
  output                                    #return output
}


#wrapper that calls the 'model' function for a specific setting
evaluate.model <- function(params=list(beta0=0.6584, sigma=1/6.4, z=12, b=0.143, a0=1/1.5, w=12, c=1, presymptomatic=1, dt=0.05),
                           init = list(S=10600000, E1=0, E2=0, E3=0, E4=0, E5=6, E6=0,
                                       I1 = 1, I2= 0, I3=0, I4=0, Iu1=0, Iu2=0, Iu3=0, Iu4=0,
                                       H=0, Ru=0, C=0),
                           nsims=2, nstep=NULL, start=as.Date("2020-03-01"),today=Sys.Date()){
  
  if(is.null(nstep)) nstep <- (as.numeric(today-start)+1+28)/params$dt #run simulation from start to current time plus four weeks
  
  xstart <- c(time=0, unlist(init), cum.time = 0) #initial conditions
  
  data <- vector(mode='list',length=nsims) #initialize list to store the output
  
  for (k in 1:nsims) {              #simulate nsims times
    data[[k]] <- as.data.frame(model(xstart,params,nstep))
    data[[k]]$cum.time <- cumsum(data[[k]]$time)
    
  }
  
  return(data)
}



######
#functions that specify time-varying rates

# gamma <- function(z = 12, b=0.143, a0=1/1.5, t){
#   # piecewise function
#   # default parameters z = 12, b=1/7, a0=1/1.5
#   #    z: time at start of intervention (notionally March 12)
#   #    b: intercept (positive)
#   #    a0: post intervention isolation ratae
#   #    t: time in the model
#   
#   gamma <- ifelse(t<=z, gamma <- b, gamma <- a0)
#   return(gamma)
# }

gamma <- function(z = 12, b=0.143, a0=1/1.5, t){
  # regression equation
  # default parameters z = 12, b=1/7, a0=1/1.5
  #    z: time at start of intervention (notionally March 12) -- if z<0 then follow regression equation
  #    b: intercept (positive)
  #    a0: post intervention isolation ratae
  #    t: time in the model
  
  if(z<0){
    gamma <- 1/(max(-0.57437*(t+12)+15.45651, a0)) # 12 day shift because regression equation starts 18 Feb and model starts 1 March ; probably needs a somewhat smaller slope due to bias, but not too bad (many observations in the region of observability)
  }
  else {
    gamma <- ifelse(t<=z, gamma <- b, gamma <- a0)
  }
  return(gamma)
}


eta <- function(t, w=12) ifelse(t<=w,1/5,1/5)   # assume average five day reporting lag based on "4-6 days" for test to be returned
#eta <- function(t, w) 1/(max(-0.59549*(t+12)+18.79376, 1))   #based on regression equation from Georgia line list, minimum report time of one day, # 12 day shift because regression equation starts 18 Feb and model starts 1 March; NOTE: this biased by unobserved data

q <- function(t, w=12, q0=0.11, q1=0.5) ifelse(t<=w,q0,q1)   #reporting rate, serious testing starts in Georgia on 18 March (subtract a 5 day reporting lag)

beta <- function(t, w=12, beta0=0.6584, beta.factor=2) ifelse(t<=w,beta0,beta0/beta.factor)


