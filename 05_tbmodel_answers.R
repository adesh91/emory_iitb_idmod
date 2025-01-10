library(deSolve)
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

# Define parameters - based on those from Dye et al. J Roy Soc Inter 2008 (PMCID: PMC3226985) and Dowdy et al. IJTLD 2014 (PMCID: PMC4041555)
parms <- c(
  
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
  m = 0           # Rate of mortality from TB 
)

# Define the differential equations model
tb_model <- function(times, init, parms) {
  with(as.list(c(init, parms)), {
    
    # Differential equations
    
    N = S+E1+E2+I+R
      
    dS <- mu*N - pi*(I/N)*S - mu*S
    dE1 <- pi*(I/N)*S + p*pi*(I/N)*E2 - nu*E1 - theta*E1 - mu*E1
    dE2 <- nu*E1 - p*pi*(I/N)*E2 - sigma*E2 - mu*E2
    dI <- theta*E1 + sigma*E2 + p*pi*(I/N)*R - tau*I - m*I - mu*I 
    dR <- tau*I - p*pi*(I/N)*R  - mu*R
    
    # Return rates of change
    list(c(dS, dE1, dE2, dI, dR))
  })
}

# Initial state values, N = 1000
init <- c(
  S = 816, E1 = 2, E2 = 140, I = 2, R = 40
) 


# Time sequence for simulation
times <- seq(0, 10000, by = 1) # for 10000 years, by year

# Solve the system of equations
out <- lsoda(y = init, times = times, func = tb_model, parms = parms)

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
  mutate(N = S + E1 + E2 + I + R ,
         Ip100k = 1e5 * ((parms["theta"]*E1) + # cases due to fast progression
                           (parms["sigma"]*E2) + # cases due to slow progression
                           (parms["p"]*parms["pi"]*(I/N)*R)) / N , # cases due to reinfection  
         Mp100k = 1e5 * (parms["m"]*I) / N ,
         inf_prev = (E1 + E2 + I + R) / N) %>%  
  tail(5)


# Question 1. In the code above, write the expressions to calculate incidence (per 100k), mortality (per 100k), and infection prevalence.

# Question 2. What is the incidence and prevalence of Mtb infection at the final time step in the model? 
# Is this is line with incidence and prevalence in India? If not, why not? 

# Incidence is 140 per 100k
# Mtb infection prevalence in 47%

# Question 4. What proportion of TB cases are due to fast progression? Recurrent (or reinfection) TB? Reactivation?
# Is that what you would expect?

# 56% from primary fast progression, 13% from reactivation, 31% from reinfection

# Question 4. How sensitive is incidence to changes in 
# tau (which is the inverse of the duration of infectious period),
# and theta (the rate of fast progression)? What about sigma, the rate of reactivation?
# What does this tell you about the relative importance of these parameters?

# Model should be fit to local data, to estimate values for each parameter that give steady-state dynamics. We won't talk about fitting today,
# but a critical first step before introducing an intervention to the model.

