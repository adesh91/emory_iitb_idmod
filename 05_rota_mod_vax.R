# Version 3
# Add Vaccination --------------------------------------------------------------
## We have already created the age structure in version 2. We need to add
## relevant vaccine parameters and a tracker compartment to explicilty
## track number of vaccinated children

# Please define parameters -- THESE ARE THE MODEL PARAMETERS
parms <- c(alpha = 14,         # alpha = daily contacts
           beta = 0.63,        # beta = probability of infection on contact
           sigma = 1/7,       # sigma = rate of recovery per day (1/duration of infection)
           mu = (16/1000)/365,         # mu =  per capita birth and death rate (16 births per 1000)
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
           age_53m = 1/1484,
           chi = 0.1
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


sir_ode_v3 <- function(times,init,parms){
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

    print(N)

    dM_0_1 <- mu*N - omega_m*M_0_1 - age_2m*M_0_1 - mu*M_0_1
    dS1_0_1 <- omega_m*M_0_1 + omega*R_0_1 - alpha*beta*I_prev*S1_0_1 - mu*S1_0_1 - age_2m*S1_0_1
    dI1_0_1 <- alpha*beta*I_prev*S1_0_1 - sigma*I1_0_1 - mu*I1_0_1 - age_2m*I1_0_1
    dS2_0_1 <- sigma*I1_0_1 - alpha*beta*I_prev*S2_0_1*imm_single_inf - mu*S2_0_1 - age_2m*S2_0_1 
    dI2_0_1 <- alpha*beta*I_prev*S2_0_1*imm_single_inf - sigma*I2_0_1 - mu*I2_0_1 - age_2m*I2_0_1
    dS3_0_1 <- sigma*I2_0_1 - alpha*beta*I_prev*S3_0_1*imm_two_inf - mu*S3_0_1 - age_2m*S3_0_1 
    dI3_0_1 <- alpha*beta*I_prev*S3_0_1*imm_two_inf - sigma*I3_0_1 - mu*I3_0_1 - age_2m*I3_0_1
    dR_0_1 <- sigma*I3_0_1 - omega*R_0_1 - mu*R_0_1 - age_2m*R_0_1

    dM_2_3 <- - omega_m*M_2_3 - age_2m*M_2_3 + age_2m*M_0_1 - mu*M_2_3
    dS1_2_3 <- omega_m*M_2_3 + omega*R_2_3 - alpha*beta*I_prev*S1_2_3 - mu*S1_2_3 - age_2m*S1_2_3 + age_2m*S1_0_1 - chi*S1_2_3 
    dI1_2_3 <- alpha*beta*I_prev*S1_2_3 - sigma*I1_2_3 - mu*I1_2_3 - age_2m*I1_2_3 + age_2m*I1_0_1
    dS2_2_3 <- sigma*I1_2_3 - alpha*beta*I_prev*S2_2_3*imm_single_inf - mu*S2_2_3 - age_2m*S2_2_3 + age_2m*S2_0_1 + chi*S1_2_3
    dI2_2_3 <- alpha*beta*I_prev*S2_2_3*imm_single_inf - sigma*I2_2_3 - mu*I2_2_3 - age_2m*I2_2_3 + age_2m*I2_0_1
    dS3_2_3 <- sigma*I2_2_3 - alpha*beta*I_prev*S3_2_3*imm_two_inf - mu*S3_2_3 - age_2m*S3_2_3 + age_2m*S3_0_1
    dI3_2_3 <- alpha*beta*I_prev*S3_2_3*imm_two_inf - sigma*I3_2_3 - mu*I3_2_3 - age_2m*I3_2_3 + age_2m*I3_0_1
    dR_2_3 <- sigma*I3_2_3 - omega*R_2_3 - mu*R_2_3 - age_2m*R_2_3 + age_2m*R_0_1

    dM_4_6 <- - omega_m*M_4_6 - age_3m*M_4_6 + age_2m*M_2_3 - mu*M_4_6
    dS1_4_6 <- omega_m*M_4_6 + omega*R_4_6 - alpha*beta*I_prev*S1_4_6 - mu*S1_4_6 - age_3m*S1_4_6 + age_2m*S1_2_3 - chi*S1_4_6
    dI1_4_6 <- alpha*beta*I_prev*S1_4_6 - sigma*I1_4_6 - mu*I1_4_6 - age_3m*I1_4_6 + age_2m*I1_2_3
    dS2_4_6 <- sigma*I1_4_6 - alpha*beta*I_prev*S2_4_6*imm_single_inf - mu*S2_4_6 - age_3m*S2_4_6 + age_2m*S2_2_3 - chi*S2_4_6 + chi*S1_4_6
    dI2_4_6 <- alpha*beta*I_prev*S2_4_6*imm_single_inf - sigma*I2_4_6 - mu*I2_4_6 - age_3m*I2_4_6 + age_2m*I2_2_3
    dS3_4_6 <- sigma*I2_4_6 - alpha*beta*I_prev*S3_4_6*imm_two_inf - mu*S3_4_6 - age_3m*S3_4_6 + age_2m*S3_2_3 + chi*S2_4_6
    dI3_4_6 <- alpha*beta*I_prev*S3_4_6*imm_two_inf - sigma*I3_4_6 - mu*I3_4_6 - age_3m*I3_4_6 + age_2m*I3_2_3
    dR_4_6 <- sigma*I3_4_6 - omega*R_4_6 - mu*R_4_6 - age_3m*R_4_6 + age_2m*R_2_3

    dM_7_59 <- - omega_m*M_7_59 - age_53m*M_7_59 + age_3m*M_4_6 - mu*M_7_59
    dS1_7_59 <- omega_m*M_7_59 + omega*R_7_59 - alpha*beta*I_prev*S1_7_59 - mu*S1_7_59 - age_53m*S1_7_59 + age_3m*S1_4_6 - chi*S1_7_59
    dI1_7_59 <- alpha*beta*I_prev*S1_7_59 - sigma*I1_7_59 - mu*I1_7_59 - age_53m*I1_7_59 + age_3m*I1_4_6
    dS2_7_59 <- sigma*I1_7_59 - alpha*beta*I_prev*S2_7_59*imm_single_inf - mu*S2_7_59 - age_53m*S2_7_59 + age_3m*S2_4_6 + chi*S1_7_59 - chi*S2_7_59
    dI2_7_59 <- alpha*beta*I_prev*S2_7_59*imm_single_inf - sigma*I2_7_59 - mu*I2_7_59 - age_53m*I2_7_59 + age_3m*I2_4_6
    dS3_7_59 <- sigma*I2_7_59 - alpha*beta*I_prev*S3_7_59*imm_two_inf - mu*S3_7_59 - age_53m*S3_7_59 + age_3m*S3_4_6 + chi*S2_7_59
    dI3_7_59 <- alpha*beta*I_prev*S3_7_59*imm_two_inf - sigma*I3_7_59 - mu*I3_7_59 - age_53m*I3_7_59 + age_3m*I3_4_6
    dR_7_59 <- sigma*I3_7_59 - omega*R_7_59 - mu*R_7_59 - age_53m*R_7_59 + age_3m*R_4_6

    dM_o59 <- - omega_m*M_o59 + age_53m*M_7_59 - mu*M_o59
    dS1_o59 <- omega_m*M_o59 + omega*R_o59 - alpha*beta*I_prev*S1_o59 - mu*S1_o59 + age_53m*S1_7_59 - chi*S1_o59
    dI1_o59 <- alpha*beta*I_prev*S1_o59 - sigma*I1_o59 - mu*I1_o59 + age_53m*I1_7_59
    dS2_o59 <- sigma*I1_o59 - alpha*beta*I_prev*S2_o59*imm_single_inf - mu*S2_o59 + age_53m*S2_7_59 + chi*S1_o59 - chi*S2_o59
    dI2_o59 <- alpha*beta*I_prev*S2_o59*imm_single_inf - sigma*I2_o59 - mu*I2_o59 + age_53m*I2_7_59
    dS3_o59 <- sigma*I2_o59 - alpha*beta*I_prev*S3_o59*imm_two_inf - mu*S3_o59 + age_53m*S3_7_59 + chi*S2_o59
    dI3_o59 <- alpha*beta*I_prev*S3_o59*imm_two_inf - sigma*I3_o59 - mu*I3_o59 + age_53m*I3_7_59
    dR_o59 <- sigma*I3_o59 - omega*R_o59 - mu*R_o59 + age_53m*R_7_59

    vaccinated = c(chi*S1_2_3 + chi*S1_4_6 + chi*S2_4_6 + chi*S1_7_59 + chi*S2_7_59 + chi*S1_o59 + chi*S2_o59)

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
      incidence_o59 = incidence_o59, vaccinated = vaccinated)
  })
}

# This creates the output from model equations.  
#If you want to run the model for longer, change the second term eg: seq(0,200,...)
times <- seq(0,365*5,1)
sir_out_v3 <- lsoda(init,times,sir_ode_v3,parms)
sir_out_long_v3 <- melt(as.data.frame(sir_out_v3),"time")

#Plotting the model output
# Plot all compartments
ggplot() +
  geom_line(data = sir_out_long_v3,  
    aes(x=time, y=value, colour=variable,group=variable)) +             #Add line
  xlab("Time")+ylab("Number") +  #Add labels
  theme_bw() +
  facet_wrap(. ~ variable) +
  ggtitle('Daily Incidence')

# Focused plot of incidence 
# age distribution at equilibrium
age_0_1 <- filter(sir_out_long_v2, time == 365*40, 
  variable %in% c('M_0_1', 'S1_0_1', 'I1_0_1', 'S2_0_1', 'I2_0_1', 'S3_0_1', 'I3_0_1', 'R_0_1'))$value %>% sum()
age_2_3 <- filter(sir_out_long_v2, time == 365*40, 
  variable %in% c('M_2_3', 'S1_2_3', 'I1_2_3', 'S2_2_3', 'I2_2_3', 'S3_2_3', 'I3_2_3', 'R_2_3'))$value %>% sum()
age_4_6 <- filter(sir_out_long_v2, time == 365*40, 
  variable %in% c('M_4_6', 'S1_4_6', 'I1_4_6', 'S2_4_6', 'I2_4_6', 'S3_4_6', 'I3_4_6', 'R_4_6'))$value %>% sum()
age_7_59 <- filter(sir_out_long_v2, time == 365*40, 
  variable %in% c('M_7_59', 'S1_7_59', 'I1_7_59', 'S2_7_59', 'I2_7_59', 'S3_7_59', 'I3_7_59', 'R_7_59'))$value %>% sum()
age_o59 <- filter(sir_out_long_v2, time == 365*40, 
  variable %in% c('M_o59', 'S1_o59', 'I1_o59', 'S2_o59', 'I2_o59', 'S3_o59', 'I3_o59', 'R_o59'))$value %>% sum()

incidence_rate_df <- filter(sir_out_long_v2, 
  variable %in% c('incidence_0_1', 'incidence_2_3', 'incidence_4_6', 'incidence_7_59', 'incidence_o59'),
  time %in% seq(365*40, 365*45, 1)) %>%
  group_by(variable) %>%
  summarize(yearly_incidence = sum(value)/5)
incidence_rate_df$population <- c(age_0_1, age_2_3, age_4_6, age_7_59, age_o59)
incidence_rate_df <- mutate(incidence_rate_df, incidence_rate = yearly_incidence/population)

ggplot(incidence_rate_df) +
  geom_bar(aes(x = variable, y = incidence_rate, fill = variable), stat = 'identity') +
  theme_bw() +
  xlab('Age Category') +
  ylab('Rate')

# Focused plot of vaccination  
ggplot() +
  geom_line(data = filter(sir_out_long_v3, time > 300, 
    variable %in% c('vaccinated')),
    aes(x=time, y=value, colour=variable,group=variable)) +             #Add line
  xlab("Time")+ylab("Number") +  #Add labels
  theme_bw() +
  ggtitle('Daily Incidence')