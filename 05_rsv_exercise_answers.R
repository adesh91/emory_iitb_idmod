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
#       - probability of infection on contact with a (symptomatically) infectious person beta = 0.3 OR an asymptomatically infectious person gamma = 0.1, multiplied by 
#       - the probability that a randomly encountered person is infectious I/N 
#   mu = 0.05 is the birth and death rate, maintaining constant population size 
#   omega = 0.02 is the rate at which immunity is lost 
#   sigma = 0.2 is the recovery rate 
#   phi = 0.2 is the rate of becoming infectious 
#   epsilon = 0.35 is the symptomatic ratio

# Define the model -------------------------------------------------------------

# Please define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha = 5,         # alpha = daily contacts
           beta = 0.3,        # beta = probability of infection on contact from symptomatic infectious
           beta1 = 0.3396,       # beta1 = amplitude of seasonal variation for symptomatic infectious
           gamma = 0.1,       # gamma = probability of infection on contact from asymptomatic infectious
           gamma1 = 0.3396,       # gamma1 = amplitude of seasonal variation from asymptomatic infectious
           sigma = 0.2,       # sigma = rate of recovery per day
           mu = 0.01,         # mu =  per capita birth and death rate
           omega = 0.02,      # omega = rate of immune loss per day
           phi = 0.2,         # phi = rate of becoming infectious 
           epsilon = 0.35)    # epsilon = symptomatic ratio

# Please define the initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
init <- c(S = 3000000,      # number initially susceptible
          E = 0,            # number initially pre-infectious
          A = 3,            # number initially asymptomatically infectious
          I = 1,            # number initially symptomatically infectious
          R = 0)            # number initially immune or "recovered"

# Please write the model equations.  They are written as a function called seair_ode.  
seair_ode <- function(times,init,parms){
  with(as.list(c(parms,init)), {
    # ODEs
    beta_t = beta * (1 + beta1 * cos(2 * pi * times / 365))
    gamma_t = gamma * (1 + gamma1 * cos(2 * pi * times / 365))
    
    dS <- mu*(S+E+A+I+R) + omega*R - alpha*beta_t*(I/(S+E+A+I+R))*S - alpha*gamma_t*(A/(S+E+A+I+R))*S - mu*S 
    dE <- alpha*beta_t*(I/(S+E+A+I+R))*S + alpha*gamma_t*(A/(S+E+A+I+R))*S - (1-epsilon)*phi*E - epsilon*phi*E - mu*E
    dA <- (1-epsilon)*phi*E - sigma*A - mu*A
    dI <- epsilon*phi*E - sigma*I - mu*I
    dR <- sigma*A + sigma*I  - omega*R - mu*R
    list(c(dS,dE,dA,dI,dR))
  })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,365*10,length.out=365*10)
seair_out <- lsoda(init,times,seair_ode,parms)
seair_out_long <- melt(as.data.frame(seair_out),"time")

#Plotting the model output
ggplot(seair_out_long,aes(x=time,y=value,colour=variable,group=variable))+
  geom_line(lwd=2)+  #Add line
  xlab("Time")+ylab("Number")                   #Add labels


# - What is the peak magnitude? 
peak_mag = max(seair_out_long %>% filter(variable == "I") %>% select(value))
peak_mag 

# - When does the peak occur (what timestep)?
# (Time at which I compartment is at maximum)
peak_time <- seair_out_long %>% filter(value == peak_mag) %>% select(time)
peak_time 

# - What is the cumulative incidence?
seair_out_wide <- seair_out_long %>% pivot_wider(., id_cols = c("time"), names_from = c("variable"), values_from = c("value"))
seair_out_wide <- seair_out_wide %>%
  mutate(incid = parms["alpha"]*parms["beta"]*(1 + parms["beta1"] * cos(2 * pi * times / 365))*(I/(S+E+A+I+R)*S) +
           parms["alpha"]*parms["gamma"]*(1 + parms["gamma1"] * cos(2 * pi * times / 365))*(A/(S+E+A+I+R)*S))
sum(seair_out_wide$incid) 

# - What is the symptomatic cumulative incidence?
seair_out_wide <- seair_out_wide %>%
  mutate(symp_incid = parms["alpha"]*parms["beta"]*(1 + parms["beta1"] * cos(2 * pi * times / 365))*(I/(S+E+A+I+R)*S))
sum(seair_out_wide$symp_incid) 

#Plot Incidence
ggplot(seair_out_wide,aes(x=time,y=incid))+
  geom_line(lwd=2)+  #Add line
  xlab("Time")+ylab("Incidence")                   #Add labels

#Plot Symptomatic Incidence
ggplot(seair_out_wide,aes(x=time,y=symp_incid))+
  geom_line(lwd=2)+  #Add line
  xlab("Time")+ylab("Symptomatic Incidence")                   #Add labels
