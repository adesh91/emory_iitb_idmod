# Session 4 Exercise

library(kableExtra)
library(deSolve)
library(reshape2)
library(ggplot2)

# A new respiratory pathogen has emerged. Here's what we know so far:
#   All individuals are susceptible and become infected by contact with an infected individual
#   People are infectious for a while and then become immune
#   Immunity eventually decays.
# We can partition our population of size N into susceptible (S), infected (I), and immune (R) people,
# the above-mentioned state variables, so N = S + I + R.

# Then, we write these as difference equations:
# In these equations,
#   lambda is the force (or rate) of infection per susceptible person. Which is: 
#       - contact rate, alpha = 5, multiplied by 
#       - probability of infection on contact with an infectious person beta = 0.1, multiplied by 
#       - the probability that a randomly encountered person is infectious I/N 
#   mu = 0.05 is the birth and death rate, maintaining constant population size 
#   omega = 0.02 is the rate at which immunity is lost 
#   sigma = 0.2 is the recovery rate 

# Define the model -------------------------------------------------------------

# Please define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha = 5,         # alpha = daily contacts
           beta = 0.1,        # beta = probability of infection on contact
           sigma = 0.2,       # sigma = rate of recovery per day
           mu = 0.05,         # mu =  per capita birth and death rate
           omega = 0.02)      # omega = rate of immune loss per day

# Please define the initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
init <- c(S = 100,          # number initially susceptible
          I = 1,            # number initially infectious
          R = 0)            # number initially immune or "recovered"

# Please write the model equations dS, dI, and dR.  They are written as a function called sir_ode.  
sir_ode <- function(times,init,parms){
  with(as.list(c(parms,init)), {
    # ODEs
    dS <- mu*(S+I+R) + omega*R - alpha*beta*(I/(S+I+R))*S - mu*S 
    dI <- alpha*beta*(I/(S+I+R))*S - sigma*I - mu*I
    dR <- sigma*I  - omega*R - mu*R
    list(c(dS,dI,dR))
  })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,100,length.out=100)
sir_out <- lsoda(init,times,sir_ode,parms)
sir_out_long <- melt(as.data.frame(sir_out),"time")

#Plotting the model output
ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable))+
  geom_line(lwd=2)+             #Add line
  xlab("Time")+ylab("Number")   #Add labels

