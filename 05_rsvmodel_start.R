# Session 4 Exercise

library(kableExtra)
library(deSolve)
library(reshape2)
library(ggplot2)

# Let's build off of this SIR model to code an RSV model:
# 1) All individuals are infectious for several days before becoming symptomatic (𝜼)
# 2) Some infectious individuals become symptomatic, and some do not. 


# Then, we write these as difference equations:
# In these equations,
#   lambda is the force (or rate) of infection per susceptible person. Which is: 
#       - contact rate, alpha = 5, multiplied by 
#       - probability of infection on contact with an infectious person beta = 0.1, multiplied by 
#       - the probability that a randomly encountered person is infectious I/N 
#   mu = 0.05 is the birth and death rate, maintaining constant population size 
#   omega = 0.02 is the rate at which immunity is lost 
#   sigma = 0.2 is the recovery rate 
#   phi = 0.2 is the rate of becoming symptomatically infectious 
#   epsilon = 0.35 is the symptomatic ratio

# Define the model -------------------------------------------------------------

# Please define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha = 5,         # alpha = daily contacts
           beta = 0.1,        # beta = probability of infection on contact
           sigma = 0.2,       # sigma = rate of recovery per day
           mu = 0.05,         # mu =  per capita birth and death rate
           omega = 0.02)      # omega = rate of immune loss per day

# Please define the initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
# Start with 3 million susceptible, 3 asymptomatically infectious, and 1 symptomatically infectious
init <- c(S = 100,          # number initially susceptible
          I = 1,            # number initially infectious
          R = 0)            # number initially immune or "recovered"

# Please write the model equations dS, dE, dA, dI, and dR.  They are written as a function called sir_ode.  
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
# Run this model for 10 years
times <- seq(0,100,length.out=100)
sir_out <- lsoda(init,times,sir_ode,parms)
sir_out_long <- melt(as.data.frame(sir_out),"time")

#Plotting the model output
ggplot(sir_out_long,aes(x=time,y=value,colour=variable,group=variable))+
  geom_line(lwd=2)+             #Add line
  xlab("Time")+ylab("Number")   #Add labels

# - What is the peak magnitude? 
peak_mag = max(sir_out_long %>% filter(variable == "I") %>% select(value))
peak_mag 

# - When does the peak occur (what timestep)?
# (Time at which I compartment is at maximum)
peak_time <- sir_out_long %>% filter(value == peak_mag) %>% select(time)
peak_time 

# - What is the cumulative incidence (better question than attack rate for this outbreak)?
sir_out_wide <- sir_out_long %>% pivot_wider(., id_cols = c("time"), names_from = c("variable"), values_from = c("value"))
sir_out_wide <- sir_out_wide %>%
  mutate(incid = parms["alpha"]*parms["beta"]*(I/(S+I+R))*S)
sum(sir_out_wide$incid) 