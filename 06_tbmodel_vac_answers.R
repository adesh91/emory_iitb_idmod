library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

# Define parameters - based on those from Dye et al. J Roy Soc Inter 2008 (PMCID: PMC3226985) and Dowdy et al. IJTLD 2014 (PMCID: PMC4041555)
parms_vac <- c(
    
  # all time units are years
  
  # demography
  mu = 0.01,    # Birth rate 
  pi = 5.89,    # Effective per capita contact rate (contact rate, or alpha * probability of infection per contact, or beta)
  # force of infection, lambda = alpha * beta * I(t)/N
  
  # natural history 
  theta = 0.03,   # Rate of fast progression
  nu = 0.2,       # Rate of immune stabilization 
  sigma = 0.0005, # Rate of reactivation  
  p = 0.5,        # Reduced risk of secondary infection (E2 -> E1, R -> E2, R -> I)
  tau = 0.9,      # Rate of recovery 
  m = 0,          # Rate of mortality from TB 
  
  # vaccination
  # AT BIRTH:
  v_b = 0.8,      # coverage of vaccination at birth (BCG)
  ve = 0.5,       # vaccine efficacy against disease
  # RANDOM:
  v_rate = 0.05,  # rate of random vaccination 
  
  # waning immunity from vaccination
  omega = 0.2
)

# Define the differential equations model
tb_model_vac <- function(times, init_vac, parms_vac) {
  with(as.list(c(init_vac, parms_vac)), {
    # Differential equations
    N <- Su + Eu1 + Eu2 + Iu + Ru + Sv + Ev1 + Ev2 + Iv + Rv
    
    dSu <- mu*N*(1-v_b) - pi*((Iu+Iv)/N)*Su - mu*Su - v_rate*Su + omega*Sv
    dEu1 <- pi*((Iu+Iv)/N)*Su + p*pi*((Iu+Iv)/N)*Eu2 - nu*Eu1 -theta*Eu1 - mu*Eu1 - v_rate*Eu1 + omega*Ev1
    dEu2 <- nu*Eu1 - p*pi*((Iu+Iv)/N)*Eu2 - sigma*Eu2 - mu*Eu2 - v_rate*Eu2 +  omega*Ev2
    dIu <- theta*Eu1 + sigma*Eu2 + p*pi*((Iu+Iv)/N)*Ru - tau*Iu - mu*Iu - v_rate*Iu + omega*Iv
    dRu <- tau*Iu - p*pi*((Iu+Iv)/N)*Ru  - mu*Ru - v_rate*Ru + omega*Rv
    
    dSv <-  mu*N*v_b - pi*((Iu+Iv)/N)*Sv - mu*Sv + v_rate*Su - omega*Sv
    dEv1 <- pi*((Iu+Iv)/N)*Sv + p*pi*((Iu+Iv)/N)*Ev2 - nu*Ev1 - theta*ve*Ev1 - mu*Ev1 + v_rate*Eu1 - omega*Ev1
    dEv2 <- nu*Ev1 - p*pi*((Iu+Iv)/N)*Ev2 - sigma*ve*Ev2 - mu*Ev2 + v_rate*Eu2 - omega*Ev2
    dIv <- theta*ve*Ev1 + sigma*ve*Ev2 + p*pi*ve*((Iu+Iv)/N)*Rv - tau*Iv - mu*Iv + v_rate*Iu - omega*Iv
    dRv <- tau*Iv - p*pi*ve*((Iu+Iv)/N)*Rv  - mu*Rv + v_rate*Ru - omega*Rv
    
    # Return rates of change
    list(c(dSu, dEu1, dEu2, dIu, dRu,
           dSv, dEv1, dEv2, dIv, dRv))
  })
}

# Initial state values
init_vac <- c(
  Su = 816, Eu1 = 2, Eu2 = 140, Iu = 2, Ru = 40,
  Sv = 0, Ev1 = 0, Ev2 = 0, Iv = 0, Rv = 0
)

# Time sequence for simulation
times <- seq(0, 1000, by = 1) # for 10000 years, by year

# Solve the system of equations
out <- lsoda(y = init_vac, times = times, func = tb_model_vac, parms = parms_vac)

# Convert output to data frame and long format for plotting
out_df <- as.data.frame(out)
out_long <- melt(as.data.frame(out_df),"time")

# Plot results with ggplot2
ggplot(out_long, aes(x = time, y = value, color = variable, group = variable)) +
  geom_line(size = 1.2) +
  xlab("Time (years)") + ylab("Population") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Calculate epidemiological measures
out_df %>%
  mutate(N = Su + Eu1 + Eu2 + Iu + Ru + Sv + Ev1 + Ev2 + Iv + Rv ,
         N_vac_p = (Sv + Ev1 + Ev2 + Iv + Rv) / N , 
         Ip100k = 1e5 * (parms_vac["theta"]*Eu1 + parms_vac["ve"]*parms_vac["theta"]*Ev1 + #cases due to fast progression
                         parms_vac["sigma"]*Eu2 + parms_vac["ve"]*parms_vac["sigma"]*Ev2 + # cases due to slow progression
                         parms_vac["p"]*parms_vac["pi"]*((Iu+Iv)/N)*Ru + parms_vac["ve"]*parms_vac["p"]*parms_vac["pi"]*((Iu+Iv)/N)*Rv) / N ) %>%   # cases due to reinfection  -
  tail(5)

# Question 1. Compare three scenarios, one with no vaccination, one with only vaccination at birth (v_b = 0.8 and ve = 0.5), 
# and one with vaccination of the entire population (‘random’) vaccination (v.rate = 0.05 and ve = 0.5) How does the incidence change? 
# Which is most effective and why?

# No vax: 140 per 100k
# Vax at birth: 0 per 100k (80% coverage) 
# Random vax: 0 per 100k (attain ~83% coverage)

# These assume *no waning*! 

# Question 2.  The random vaccination currently represented in the model acts on which states? Does this vaccine prevent infection (prevention-of-infection, ‘POI’) 
# or disease (prevention-of-disease, ‘POD’)? How would you represent a POI vaccine?

# Currently (POD): E1, E2, R
# For POI: S , E2, R (states which are subject to force of infection / reinfection)

# Question 3. Protection from BCG against pulmonary tuberculosis (the infectious form of tuberculosis) wanes after approximately 10 years. (Martinez et al. Lancet ID, 2022)  
# Modify the code to incorporate waning immunity from vaccination. Start with a model that includes both vaccination at birth (v_b = 0.8) 
# and vaccination at random (ve_d = 0.5, c = 0.05). Then, add a rate of waning immunity (omega) of 0.1 per year (equivalent to 10 years of protection). 
# Then, test a rate of 0.2 per year (5 years of protection). How does this change tuberculosis incidence? 

# omega = 0.1 per year:incidence is >0, 37% of pop is protected
# omega = 0.2 per year:incidence is higher, 22% of pop is protected

# Question 4. How would you calculate vaccine impact? (In other words, the total cases averted under a vaccination scenario?)

# a 'counterfactual' comparison - compare a baseline scenario without vaccination to one with vaccination
