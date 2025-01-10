# Session 4 Exercise
library(kableExtra)
library(deSolve)
library(reshape2)
library(tidyverse)

####
# Rotavirus models:
## v1. Simple
## v2. Added age categories
## v3. Added vaccination

# Version 1
# Define the model -------------------------------------------------------------
# Please define parameters -- THESE ARE THE MODEL PARAMETERS
# we asume 1 month = 28 days
parms <- c(alpha = 14,         # alpha = daily contacts
           beta = 0.63,        # beta = probability of infection on contact
           sigma = 1/7,       # sigma = rate of recovery per day (1/duration of infection)
           mu = (16/1000),         # mu =  per capita birth and death rate (16 births per 1000)
           omega = 1/(28*12*250),      # omega = rate of waning immunity per day; Immunity lasts 250 years
           omega_m = 1/(6*28),  # omega_m = rate of waning for maternal immunity per day; lasts 6 months
           symp1 = 0.30, # proportion of primary infections that are symptomatic
           symp2 = 0.28, # proportion of secondary infections that are symptomatic
           symp3 = 0.18, # proportion of tertiary infections that are symptomatic
           imm_single_inf = 0.74, # relative immunity of having 1 prior infection vs. none
           imm_two_inf = 0.74, # relative immunity of having two prior infections vs. none 
           rel_inf_2 = 0.5,
           rel_inf_3 = 0.1)  

# Please define the initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
init <- c(M = 0,           # number with maternal immunity
          S1 = 98,          # number of completely susceptible
          S2 = 0,            # number of susceptible with 1 prior infection
          S3 = 0,            # number of susceptible with 2 prior infection
          I1 = 2,            # number of individuals with first infection
          I2 = 0,            # number of individuals with second infection
          I3 = 0,            # number of individuals with third infection
          R = 0)             # number of individuals with durable immunity

# Please write the model equations dS, dI, and dR.  They are written as a function called sir_ode.  
sir_ode <- function(times,init,parms){
  with(as.list(c(parms,init)), {
    # ODEs
    N <- sum(M, S1, S2, S3, I1, I2, I3, R)
    I_prev <- (symp1*I1 + (symp2*I2*rel_inf_2) + (symp3*I3*rel_inf_3))/N

    dM <- mu*N - omega_m*M - mu*M

    dS1 <- omega_m*M + omega*R - alpha*beta*I_prev*S1 - mu*S1
    dI1 <- alpha*beta*I_prev*S1 - sigma*I1 - mu*I1

    dS2 <- sigma*I1 - alpha*beta*I_prev*S2*imm_single_inf - mu*S2 
    dI2 <- alpha*beta*I_prev*S2*imm_single_inf - sigma*I2 - mu*I2

    dS3 <- sigma*I2 - alpha*beta*I_prev*S3*imm_two_inf - mu*S3 
    dI3 <- alpha*beta*I_prev*S3*imm_two_inf - sigma*I3 - mu*I3

    dR <- sigma*I3 - omega*R - mu*R

    incidence <- alpha*beta*I_prev*S1 + alpha*beta*I_prev*S2 +
      alpha*beta*I_prev*S3

    list(c(dM,dS1, dS2, dS3, dI1, dI2, dI3, dR), incidence = incidence)
  })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,365*100,1)
sir_out <- lsoda(init,times,sir_ode,parms)
sir_out_long <- melt(as.data.frame(sir_out),"time")
365*(filter(sir_out_long, time == max(sir_out_long$time), variable == 'incidence')$value)

#Plotting the model output
ggplot() +
  geom_line(data = filter(sir_out_long, variable %in% c('incidence'), time > 365*5),
    aes(x=time, y=value, colour=variable,group=variable)) +             #Add line
  xlab("Time")+ylab("Number") +  #Add labels
  theme_bw() +
  facet_wrap(. ~ variable) +
  ggtitle('Daily Incidence')


# Version 2
# Add Age Categories --------------------------------------------------------------
# Important age categories for rotavirus:
## 0-1 months: before eligible for vaccine
## 2-3 months: eligible for dose 1 of vaccine and maternal immunity
## 4-6 months: eligible for dose 2 of vaccine and maternal immunity
## 7 months to 5 years: no maternal immunity and under 5 (target demographic for rotavirus)
## Over 5 years: Older population that is not the focus of rotavirus epidemiology

## Therefore we need 5 strata for the model

# Please define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha = 14,         # alpha = daily contacts
           beta = 0.63,        # beta = probability of infection on contact
           sigma = 1/7,       # sigma = rate of recovery per day (1/duration of infection)
           mu = 16/1000,         # mu =  per capita birth and death rate (16 births per 1000)
           omega = 1/(28*12*250),      # omega = rate of waning immunity per day; Immunity lasts 250 years
           omega_m = 1/(6*28),  # omega_m = rate of waning for maternal immunity per day; lasts 6 months
           symp1 = 0.30, # proportion of primary infections that are symptomatic
           symp2 = 0.28, # proportion of secondary infections that are symptomatic
           symp3 = 0.18, # proportion of tertiary infections that are symptomatic
           imm_single_inf = 0.74, # relative immunity of having 1 prior infection vs. none
           imm_two_inf = 0.74, # relative immunity of having two prior infections vs. none 
           rel_inf_2 = 0.5,
           rel_inf_3 = 0.1,   
           age_2m = 1/56,
           age_3m = 1/84,
           age_53m = 1/1484
           )    

# Please define the initial conditions --  THESE ARE THE CONDITIONS AT THE START OF THE SIMULATION
## replicate initial states for each stratum of the model
init <- c(

          M_0_1 = 0,  S1_0_1 = 1,  S2_0_1 = 0,  S3_0_1 = 0,  I1_0_1 = 1,  I2_0_1 = 0,  I3_0_1 = 0,  R_0_1 = 0,

          M_2_3 = 0,  S1_2_3 = 1,  S2_2_3 = 0,  S3_2_3 = 0,  I1_2_3 = 1,  I2_2_3 = 0,  I3_2_3 = 0,  R_2_3 = 0,

          M_4_6 = 0,  S1_4_6 = 3,  S2_4_6 = 0,  S3_4_6 = 0,  I1_4_6 = 1,  I2_4_6 = 0,  I3_4_6 = 0,  R_4_6 = 0,

          M_7_59 = 0, S1_7_59 = 71, S2_7_59 = 0, S3_7_59 = 0, I1_7_59 = 1, I2_7_59 = 0, I3_7_59 = 0, R_7_59 = 0,

          M_o59 = 0,    S1_o59 = 1920,  S2_o59 = 0,  S3_o59 = 0,  I1_o59 = 0,  I2_o59 = 0,  I3_o59 = 0,  R_o59 = 0
)

# Please write the model equations dS, dI, and dR.  They are written as a function called sir_ode.
# We need to add aging rates to each compartment so individuals age into higher age categories at appropriate rates
# This can be done using parameters to change the time-scale of the simulation or as fixed numbers if the time-scale
# of the simulation will not change.  
sir_ode_v2 <- function(times,init,parms){
  with(as.list(c(parms,init)), {
    # ODEs
    N <- sum(
      M_0_1, S1_0_1, I1_0_1, S2_0_1, I2_0_1, S3_0_1, I3_0_1, R_0_1,
      M_2_3, S1_2_3, I1_2_3, S2_2_3, I2_2_3, S3_2_3, I3_2_3, R_2_3,
      M_4_6, S1_4_6, I1_4_6, S2_4_6, I2_4_6, S3_4_6, I3_4_6, R_4_6,
      M_7_59, S1_7_59, I1_7_59, S2_7_59, I2_7_59, S3_7_59, I3_7_59, R_7_59,
      M_o59, S1_o59, I1_o59, S2_o59, I2_o59, S3_o59, I3_o59, R_o59)
    I_sum <- (symp1*(I1_0_1 + I1_2_3 + I1_4_6 + I1_7_59 + I1_o59) + 
        (symp2*rel_inf_2*(I2_0_1 + I2_2_3 + I2_4_6 + I2_7_59 + I2_o59)) + 
        (symp3*rel_inf_3*(I3_0_1 + I3_2_3 + I3_4_6 + I3_7_59 + I3_o59))
      )
    I_prev <- I_sum/N

    dM_0_1 <- mu*N - omega_m*M_0_1 - age_2m*M_0_1 - mu*M_0_1
    dS1_0_1 <- omega_m*M_0_1 + omega*R_0_1 - alpha*beta*I_prev*S1_0_1 - mu*S1_0_1 - age_2m*S1_0_1
    dI1_0_1 <- alpha*beta*I_prev*S1_0_1 - sigma*I1_0_1 - mu*I1_0_1 - age_2m*I1_0_1
    dS2_0_1 <- sigma*I1_0_1 - alpha*beta*I_prev*S2_0_1*imm_single_inf - mu*S2_0_1 - age_2m*S2_0_1 
    dI2_0_1 <- alpha*beta*I_prev*S2_0_1*imm_single_inf - sigma*I2_0_1 - mu*I2_0_1 - age_2m*I2_0_1
    dS3_0_1 <- sigma*I2_0_1 - alpha*beta*I_prev*S3_0_1*imm_two_inf - mu*S3_0_1 - age_2m*S3_0_1 
    dI3_0_1 <- alpha*beta*I_prev*S3_0_1*imm_two_inf - sigma*I3_0_1 - mu*I3_0_1 - age_2m*I3_0_1
    dR_0_1 <- sigma*I3_0_1 - omega*R_0_1 - mu*R_0_1 - age_2m*R_0_1

    dM_2_3 <- - omega_m*M_2_3 - age_2m*M_2_3 + age_2m*M_0_1 - mu*M_2_3
    dS1_2_3 <- omega_m*M_2_3 + omega*R_2_3 - alpha*beta*I_prev*S1_2_3 - mu*S1_2_3 - age_2m*S1_2_3 + age_2m*S1_0_1 
    dI1_2_3 <- alpha*beta*I_prev*S1_2_3 - sigma*I1_2_3 - mu*I1_2_3 - age_2m*I1_2_3 + age_2m*I1_0_1
    dS2_2_3 <- sigma*I1_2_3 - alpha*beta*I_prev*S2_2_3*imm_single_inf - mu*S2_2_3 - age_2m*S2_2_3 + age_2m*S2_0_1
    dI2_2_3 <- alpha*beta*I_prev*S2_2_3*imm_single_inf - sigma*I2_2_3 - mu*I2_2_3 - age_2m*I2_2_3 + age_2m*I2_0_1
    dS3_2_3 <- sigma*I2_2_3 - alpha*beta*I_prev*S3_2_3*imm_two_inf - mu*S3_2_3 - age_2m*S3_2_3 + age_2m*S3_0_1
    dI3_2_3 <- alpha*beta*I_prev*S3_2_3*imm_two_inf - sigma*I3_2_3 - mu*I3_2_3 - age_2m*I3_2_3 + age_2m*I3_0_1
    dR_2_3 <- sigma*I3_2_3 - omega*R_2_3 - mu*R_2_3 - age_2m*R_2_3 + age_2m*R_0_1

    dM_4_6 <- - omega_m*M_4_6 - age_3m*M_4_6 + age_2m*M_2_3 - mu*M_4_6
    dS1_4_6 <- omega_m*M_4_6 + omega*R_4_6 - alpha*beta*I_prev*S1_4_6 - mu*S1_4_6 - age_3m*S1_4_6 + age_2m*S1_2_3
    dI1_4_6 <- alpha*beta*I_prev*S1_4_6 - sigma*I1_4_6 - mu*I1_4_6 - age_3m*I1_4_6 + age_2m*I1_2_3
    dS2_4_6 <- sigma*I1_4_6 - alpha*beta*I_prev*S2_4_6*imm_single_inf - mu*S2_4_6 - age_3m*S2_4_6 + age_2m*S2_2_3 
    dI2_4_6 <- alpha*beta*I_prev*S2_4_6*imm_single_inf - sigma*I2_4_6 - mu*I2_4_6 - age_3m*I2_4_6 + age_2m*I2_2_3
    dS3_4_6 <- sigma*I2_4_6 - alpha*beta*I_prev*S3_4_6*imm_two_inf - mu*S3_4_6 - age_3m*S3_4_6 + age_2m*S3_2_3 
    dI3_4_6 <- alpha*beta*I_prev*S3_4_6*imm_two_inf - sigma*I3_4_6 - mu*I3_4_6 - age_3m*I3_4_6 + age_2m*I3_2_3
    dR_4_6 <- sigma*I3_4_6 - omega*R_4_6 - mu*R_4_6 - age_3m*R_4_6 + age_2m*R_2_3

    dM_7_59 <- - omega_m*M_7_59 - age_53m*M_7_59 + age_3m*M_4_6 - mu*M_7_59
    dS1_7_59 <- omega_m*M_7_59 + omega*R_7_59 - alpha*beta*I_prev*S1_7_59 - mu*S1_7_59 - age_53m*S1_7_59 + age_3m*S1_4_6
    dI1_7_59 <- alpha*beta*I_prev*S1_7_59 - sigma*I1_7_59 - mu*I1_7_59 - age_53m*I1_7_59 + age_3m*I1_4_6
    dS2_7_59 <- sigma*I1_7_59 - alpha*beta*I_prev*S2_7_59*imm_single_inf - mu*S2_7_59 - age_53m*S2_7_59 + age_3m*S2_4_6 
    dI2_7_59 <- alpha*beta*I_prev*S2_7_59*imm_single_inf - sigma*I2_7_59 - mu*I2_7_59 - age_53m*I2_7_59 + age_3m*I2_4_6
    dS3_7_59 <- sigma*I2_7_59 - alpha*beta*I_prev*S3_7_59*imm_two_inf - mu*S3_7_59 - age_53m*S3_7_59 + age_3m*S3_4_6 
    dI3_7_59 <- alpha*beta*I_prev*S3_7_59*imm_two_inf - sigma*I3_7_59 - mu*I3_7_59 - age_53m*I3_7_59 + age_3m*I3_4_6
    dR_7_59 <- sigma*I3_7_59 - omega*R_7_59 - mu*R_7_59 - age_53m*R_7_59 + age_3m*R_4_6

    dM_o59 <- - omega_m*M_o59 + age_53m*M_7_59 - mu*M_o59
    dS1_o59 <- omega_m*M_o59 + omega*R_o59 - alpha*beta*I_prev*S1_o59 - mu*S1_o59 + age_53m*S1_7_59
    dI1_o59 <- alpha*beta*I_prev*S1_o59 - sigma*I1_o59 - mu*I1_o59 + age_53m*I1_7_59
    dS2_o59 <- sigma*I1_o59 - alpha*beta*I_prev*S2_o59*imm_single_inf - mu*S2_o59 + age_53m*S2_7_59
    dI2_o59 <- alpha*beta*I_prev*S2_o59*imm_single_inf - sigma*I2_o59 - mu*I2_o59 + age_53m*I2_7_59
    dS3_o59 <- sigma*I2_o59 - alpha*beta*I_prev*S3_o59*imm_two_inf - mu*S3_o59 + age_53m*S3_7_59
    dI3_o59 <- alpha*beta*I_prev*S3_o59*imm_two_inf - sigma*I3_o59 - mu*I3_o59 + age_53m*I3_7_59
    dR_o59 <- sigma*I3_o59 - omega*R_o59 - mu*R_o59 + age_53m*R_7_59

    incidence_0_1 = c(alpha*beta*I_prev*S1_0_1 + alpha*beta*I_prev*S2_0_1*imm_single_inf + alpha*beta*I_prev*S3_0_1*imm_two_inf)
    incidence_2_3 = c(alpha*beta*I_prev*S1_2_3 + alpha*beta*I_prev*S2_2_3*imm_single_inf + alpha*beta*I_prev*S3_2_3*imm_two_inf)
    incidence_4_6 = c(alpha*beta*I_prev*S1_4_6 + alpha*beta*I_prev*S2_4_6*imm_single_inf + alpha*beta*I_prev*S3_4_6*imm_two_inf)
    incidence_7_59 = c(alpha*beta*I_prev*S1_7_59 + alpha*beta*I_prev*S2_7_59*imm_single_inf + alpha*beta*I_prev*S3_7_59*imm_two_inf)
    incidence_o59 = c(alpha*beta*I_prev*S1_o59 + alpha*beta*I_prev*S2_o59*imm_single_inf + alpha*beta*I_prev*S3_o59*imm_two_inf)


    list(c(dM_0_1,  dS1_0_1, dS2_0_1,  dS3_0_1,  dI1_0_1,  dI2_0_1,  dI3_0_1,  dR_0_1,

          dM_2_3,  dS1_2_3,  dS2_2_3,  dS3_2_3,  dI1_2_3,  dI2_2_3,  dI3_2_3,  dR_2_3,

          dM_4_6,  dS1_4_6,  dS2_4_6,  dS3_4_6,  dI1_4_6,  dI2_4_6,  dI3_4_6,  dR_4_6,

          dM_7_59, dS1_7_59, dS2_7_59, dS3_7_59, dI1_7_59, dI2_7_59, dI3_7_59, dR_7_59,

          dM_o59, dS1_o59, dS2_o59, dS3_o59, dI1_o59, dI2_o59, dI3_o59, dR_o59), 

      incidence_0_1 = incidence_0_1, incidence_2_3 = incidence_2_3,
      incidence_4_6 = incidence_4_6, incidence_7_59 = incidence_7_59, 
      incidence_o59 = incidence_o59)
  })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,365*5,1)
sir_out_v2 <- lsoda(init,times,sir_ode_v2,parms)
sir_out_long_v2 <- melt(as.data.frame(sir_out_v2),"time")

#Plotting the model output
# Plot all compartments
ggplot() +
  geom_line(data = sir_out_long_v2, 
    aes(x=time, y=value, colour=variable,group=variable)) +             #Add line
  xlab("Time")+ylab("Number") +  #Add labels
  theme_bw() +
  facet_wrap(. ~ variable) +
  ggtitle('Daily Incidence')

# Focused plot of incidence  
ggplot() +
  geom_line(data = filter(sir_out_long_v2, time > 300, 
    variable %in% c('incidence_0_1', 'incidence_2_3', 'incidence_4_6', 'incidence_7_59', 'incidence_o59')),
    aes(x=time, y=value, colour=variable,group=variable)) +             #Add line
  xlab("Time")+ylab("Number") +  #Add labels
  theme_bw() +
  ggtitle('Daily Incidence')

