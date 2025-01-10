# ------------------------------------------------------------------------------
# 
# GlobalMix India incidence by site
# Simple, non-calibrated SIR using EpiModel
# 
# Based on Analysis by Sara Kim, October 24, 2024
#
# ------------------------------------------------------------------------------

# --- Library
library(EpiModel)
library(ggplot2)
library(patchwork)
library(ggtext)
library(tidyr)
library(dplyr)
library(EpiEstim)


# Model structure --------------------------------------------------------------
sirmod <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # 1. track total sizes of each age group
    num.1 <- s.num.1 + i.num.1 + r.num.1
    num.2 <- s.num.2 + i.num.2 + r.num.2
    num.3 <- s.num.3 + i.num.3 + r.num.3
    num.4 <- s.num.4 + i.num.4 + r.num.4
    num.5 <- s.num.5 + i.num.5 + r.num.5
    num.6 <- s.num.6 + i.num.6 + r.num.6
    num.7 <- s.num.7 + i.num.7 + r.num.7
    num.8 <- s.num.8 + i.num.8 + r.num.8
    num.9 <- s.num.9 + i.num.9 + r.num.9
    
    # 2. define age-specific forces of infection
    # force of infection for <6m0 contacts (infecting the <60mo year olds)

    lambda.1 <- q*c11*((i.num.1)/num.1) + 
      q*c12*((i.num.2)/num.2) + # c12 is contact between <6mo with 6-11mo
      q*c13*((i.num.3)/num.3) +  # c13 is contact between <6mo with 1-4 y/o, etc
      q*c14*((i.num.4)/num.4) +
      q*c15*((i.num.5)/num.5) + 
      q*c16*((i.num.6)/num.6) + 
      q*c17*((i.num.7)/num.7) + 
      q*c18*((i.num.8)/num.8) + 
      q*c19*((i.num.9)/num.9)
    
    
    # force of infection for 6-11mo contacts (infecting the 6-11 mo)
    lambda.2 <- q*c21*((i.num.1)/num.1) + 
      q*c22*((i.num.2)/num.2) + 
      q*c23*((i.num.3)/num.3) + 
      q*c24*((i.num.4)/num.4) +
      q*c25*((i.num.5)/num.5) + 
      q*c26*((i.num.6)/num.6) +
      q*c27*((i.num.7)/num.7) +
      q*c28*((i.num.8)/num.8) +
      q*c29*((i.num.9)/num.9)
  
    
    # force of infection for 1-4 y contacts (infecting the 1-4y )
    lambda.3 <- q*c31*((i.num.1)/num.1) + 
      q*c32*((i.num.2)/num.2) + 
      q*c33*((i.num.3)/num.3) + 
      q*c34*((i.num.4)/num.4) +
      q*c35*((i.num.5)/num.5) + 
      q*c36*((i.num.6)/num.6) +
      q*c37*((i.num.7)/num.7) +
      q*c38*((i.num.8)/num.8) +
      q*c39*((i.num.9)/num.9)
 
    
    # force of infection for 5-9y contacts (infecting the 5-9y)
    lambda.4 <- q*c41*((i.num.1)/num.1) + 
      q*c42*((i.num.2)/num.2) + 
      q*c43*((i.num.3)/num.3) + 
      q*c44*((i.num.4)/num.4) +
      q*c45*((i.num.5)/num.5) + 
      q*c46*((i.num.6)/num.6) +
      q*c47*((i.num.7)/num.7) +
      q*c48*((i.num.8)/num.8) +
      q*c49*((i.num.9)/num.9)
    
    
    # force of infection for 10-19 contacts (infecting the 10-19 year olds)
    lambda.5 <- q*c51*((i.num.1)/num.1) + 
      q*c52*((i.num.2)/num.2) + 
      q*c53*((i.num.3)/num.3) + 
      q*c54*((i.num.4)/num.4) +
      q*c55*((i.num.5)/num.5) + 
      q*c56*((i.num.6)/num.6) +
      q*c57*((i.num.7)/num.7) +
      q*c58*((i.num.8)/num.8) +
      q*c59*((i.num.9)/num.9)
    
    
    # force of infection for 20-29y contacts (infecting the 20-29 year olds)
    lambda.6 <- q*c61*((i.num.1)/num.1) + 
      q*c62*((i.num.2)/num.2) + 
      q*c63*((i.num.3)/num.3) + 
      q*c64*((i.num.4)/num.4) +
      q*c65*((i.num.5)/num.5) + 
      q*c66*((i.num.6)/num.6) +
      q*c67*((i.num.7)/num.7) +
      q*c68*((i.num.8)/num.8) +
      q*c69*((i.num.9)/num.9)
    
    
    # force of infection for 30-39y contacts (infecting the 30-39 year olds)
    lambda.7 <- q*c71*((i.num.1)/num.1) + 
      q*c72*((i.num.2)/num.2) + 
      q*c73*((i.num.3)/num.3) + 
      q*c74*((i.num.4)/num.4) +
      q*c75*((i.num.5)/num.5) + 
      q*c76*((i.num.6)/num.6) +
      q*c77*((i.num.7)/num.7) +
      q*c78*((i.num.8)/num.8) +
      q*c79*((i.num.9)/num.9)
  
    
    # force of infection for 40-59y contacts (infecting the 40-59 year olds)
    lambda.8 <- q*c81*((i.num.1)/num.1) + 
      q*c82*((i.num.2)/num.2) + 
      q*c83*((i.num.3)/num.3) + 
      q*c84*((i.num.4)/num.4) +
      q*c85*((i.num.5)/num.5) + 
      q*c86*((i.num.6)/num.6) +
      q*c87*((i.num.7)/num.7) +
      q*c88*((i.num.8)/num.8) +
      q*c89*((i.num.9)/num.9)
    
    
    # force of infection for 60+y contacts (infecting the 60+ year olds)
    lambda.9 <- q*c91*((i.num.1)/num.1) + 
      q*c92*((i.num.2 )/num.2) + 
      q*c93*((i.num.3)/num.3) + 
      q*c94*((i.num.4)/num.4) +
      q*c95*((i.num.5)/num.5) + 
      q*c96*((i.num.6)/num.6) +
      q*c97*((i.num.7)/num.7) +
      q*c98*((i.num.8)/num.8) +
      q*c99*((i.num.9)/num.9)
    
    
    # 3. differential equations 
    dS.1 <- -lambda.1*s.num.1
    dI.1 <- lambda.1*s.num.1 - gamma*i.num.1
    dR.1 <- gamma*i.num.1

    
    dS.2 <- -lambda.2*s.num.2
    dI.2 <- lambda.2*s.num.2 - gamma*i.num.2
    dR.2 <- gamma*i.num.2
    
    dS.3 <- -lambda.3*s.num.3
    dI.3 <- lambda.3*s.num.3 - gamma*i.num.3
    dR.3 <- gamma*i.num.3
    
    dS.4 <- -lambda.4*s.num.4
    dI.4 <- lambda.4*s.num.4 - gamma*i.num.4
    dR.4 <- gamma*i.num.4
    
    dS.5 <- -lambda.5*s.num.5
    dI.5 <- lambda.5*s.num.5 - gamma*i.num.5
    dR.5 <- gamma*i.num.5
    
    dS.6 <- -lambda.6*s.num.6
    dI.6 <- lambda.6*s.num.6 - gamma*i.num.6
    dR.6 <- gamma*i.num.6
    
    dS.7 <- -lambda.7*s.num.7
    dI.7 <- lambda.7*s.num.7 - gamma*i.num.7
    dR.7 <- gamma*i.num.7
    
    dS.8 <- -lambda.8*s.num.8
    dI.8 <- lambda.8*s.num.8 - gamma*i.num.8
    dR.8 <- gamma*i.num.8
    
    dS.9 <- -lambda.9*s.num.9
    dI.9 <- lambda.9*s.num.9 - gamma*i.num.9
    dR.9 <- gamma*i.num.9
    
    
    # 4. List outputs
    list(c(dS.1, dI.1, dR.1, 
           dS.2, dI.2, dR.2, 
           dS.3, dI.3, dR.3, 
           dS.4, dI.4, dR.4,
           dS.5, dI.5, dR.5,
           dS.6, dI.6, dR.6, 
           dS.7, dI.7, dR.7, 
           dS.8, dI.8, dR.8, 
           dS.9, dI.9, dR.9, 
           si.flow.1 = lambda.1*s.num.1,
           si.flow.2 = lambda.2*s.num.2,
           si.flow.3 = lambda.3*s.num.3,
           si.flow.4 = lambda.4*s.num.4,
           si.flow.5 = lambda.5*s.num.5,
           si.flow.6 = lambda.6*s.num.6,
           si.flow.7 = lambda.7*s.num.7,
           si.flow.8 = lambda.8*s.num.8,
           si.flow.9 = lambda.9*s.num.9,
           si.flow = lambda.1*s.num.1 + lambda.2*s.num.2 + lambda.3*s.num.3 + lambda.4*s.num.4 + 
             lambda.5*s.num.5 + lambda.6*s.num.6 + lambda.7*s.num.7 + lambda.8*s.num.8 + lambda.9*s.num.9
    ))
  })
}


# India specific parameters and conditions -------------------------------------

## Parameters ------------------------------------------------------------------
### Urban ----------------------------------------------------------------------
param.urban.ind <- param.dcm(gamma = 1/7, q = 0.015, 
                         c11=0.15, c12=0.04, c13=1.07, c14=0.51, c15=1.02, c16=3.55, c17=1.93, c18=2.60, c19=0.82,
                         c21=0.07, c22=0.07, c23=1.09, c24=0.96, c25=0.96, c26=3.76, c27=2.42, c28=2.47, c29=0.80,
                         c31=0.22, c32=0.22, c33=2.38, c34=1.37, c35=1.14, c36=3.17, c37=2.42, c38=1.98, c39=0.82,
                         c41=0.03, c42=0.03, c43=0.92, c44=6.89, c45=2.89, c46=1.57, c47=3.65, c48=2.29, c49=0.86,
                         c51=0.06, c52=0.00, c53=0.30, c54=0.77, c55=9.66, c56=1.55, c57=2.86, c58=2.62, c59=0.71,
                         c61=0.06, c62=0.16, c63=1.03, c64=0.68, c65=1.54, c66=5.52, c67=2.11, c68=2.83, c69=0.78,
                         c71=0.10, c72=0.06, c73=0.92, c74=1.16, c75=2.25, c76=2.33, c77=4.37, c78=2.97, c79=0.87,
                         c81=0.03, c82=0.00, c83=0.80, c84=0.50, c85=1.73, c86=2.69, c87=3.16, c88=5.42, c89=1.00,
                         c91=0.03, c92=0.00, c93=0.47, c94=0.58, c95=1.16, c96=2.27, c97=2.24, c98=4.23, c99=2.08)
### Rural ----------------------------------------------------------------------
param.rural.ind <- param.dcm(gamma = 1/7, q = 0.015, 
                         c11=0.31, c12=0.10, c13=0.79, c14=0.34, c15=1.10, c16=2.76, c17=2.07, c18=2.52, c19=1.09,
                         c21=0.07, c22=0.00, c23=0.35, c24=0.20, c25=0.98, c26=2.24, c27=1.35, c28=2.28, c29=1.06,
                         c31=0.26, c32=0.13, c33=2.25, c34=1.64, c35=1.08, c36=2.33, c37=3.15, c38=2.54, c39=1.82,
                         c41=0.08, c42=0.11, c43=1.56, c44=6.60, c45=2.98, c46=1.05, c47=3.51, c48=2.08, c49=1.51,
                         c51=0.03, c52=0.03, c53=0.50, c54=0.90, c55=10.06, c56=1.82, c57=2.49, c58=3.39, c59=1.56,
                         c61=0.06, c62=0.00, c63=0.65, c64=0.63, c65=1.53, c66=4.58, c67=2.39, c68=3.74, c69=1.42,
                         c71=0.03, c72=0.03, c73=1.11, c74=1.05, c75=2.05, c76=1.35, c77=3.70, c78=3.40, c79=1.98,
                         c81=0.06, c82=0.01, c83=0.51, c84=0.55, c85=1.38, c86=2.17, c87=2.32, c88=4.46, c89=1.89,
                         c91=0.00, c92=0.00, c93=0.24, c94=0.60, c95=1.13, c96=0.84, c97=1.89, c98=4.08, c99=2.37)

## Initial conditions   --------------------------------------------------------
### Urban   --------------------------------------------------------------------
init.urb.ind <- init.dcm(s.num.1 = 10900000,i.num.1 = 1, r.num.1 = 0, 
                     s.num.2 = 10900000, i.num.2 = 1,r.num.2 = 0, 
                     s.num.3 = 44700000, i.num.3 = 1, r.num.3 = 0, 
                     s.num.4 = 59600000, i.num.4 = 1, r.num.4 = 0, 
                     s.num.5 = 109300000, i.num.5 = 1, r.num.5 = 0, 
                     s.num.6 = 89500000, i.num.6 = 1, r.num.6 = 0, 
                     s.num.7 = 64600000, i.num.7 = 1, r.num.7 = 0, 
                     s.num.8 = 69600000, i.num.8 = 1, r.num.8 = 0,  
                     s.num.9 = 39800000, i.num.9 = 1, r.num.9 = 0, 
                     si.flow.1 = 0, si.flow.2 = 0, si.flow.3 = 0, si.flow.4 = 0,
                     si.flow.5 = 0, si.flow.6 = 0,  si.flow.7 = 0, si.flow.8 = 0, si.flow.9 = 0, si.flow = 0)
### Rural   ---------------------------------------------------------------------
init.rur.ind <- init.dcm(s.num.1 = 20300000, i.num.1 = 1, r.num.1 = 0, 
                     s.num.2 = 20300000, i.num.2 = 1, r.num.2 = 0, 
                     s.num.3 = 83100000, i.num.3 = 1,  r.num.3 = 0, 
                     s.num.4 = 110800000, i.num.4 = 1,  r.num.4 = 0, 
                     s.num.5 = 203100000, i.num.5 = 1,  r.num.5 = 0, 
                     s.num.6 = 166100000, i.num.6 = 1,  r.num.6 = 0, 
                     s.num.7 = 120000000,  i.num.7 = 1,  r.num.7 = 0, 
                     s.num.8 = 129200000, i.num.8 = 1,  r.num.8 = 0, 
                     s.num.9 = 73800000, i.num.9 = 1, r.num.9 = 0, 
                     si.flow.1 = 0, si.flow.2 = 0, si.flow.3 = 0, si.flow.4 = 0,
                     si.flow.5 = 0, si.flow.6 = 0,  si.flow.7 = 0, si.flow.8 = 0, si.flow.9 = 0, si.flow = 0)


## Controls    -----------------------------------------------------------------
control <- control.dcm(nstep=730, new.mod = sirmod)

## Compile model GlobalMix urban   ---------------------------------------------
sim.no.urb.ind <- dcm(param.urban.ind, init.urb.ind, control)
df.no.urb.ind <- as.data.frame(sim.no.urb.ind)

## Compile model GlobalMix rural  ----------------------------------------------
sim.no.rur.ind <- dcm(param.rural.ind, init.rur.ind, control)
df.no.rur.ind <- as.data.frame(sim.no.rur.ind)


# TODO: Plot Incidence Curves for Urban and Rural using df.no.urb.ind and df.no.rur.ind ####

# TODO: Calculate Incidence Overall for Urban and Rural ####


# Calculate risk for each age group x site   -----------------------------------

## India rural   ---------------------------------------------------------------
ind.rur.1 <- sum(df.no.rur.ind$si.flow.1) / init.rur.ind$s.num.1
ind.rur.2 <- sum(df.no.rur.ind$si.flow.2) / init.rur.ind$s.num.2
ind.rur.3 <- sum(df.no.rur.ind$si.flow.3) / init.rur.ind$s.num.3
ind.rur.4 <- sum(df.no.rur.ind$si.flow.4) / init.rur.ind$s.num.4
ind.rur.5 <- sum(df.no.rur.ind$si.flow.5) / init.rur.ind$s.num.5
ind.rur.6 <- sum(df.no.rur.ind$si.flow.6) / init.rur.ind$s.num.6
ind.rur.7 <- sum(df.no.rur.ind$si.flow.7) / init.rur.ind$s.num.7
ind.rur.8 <- sum(df.no.rur.ind$si.flow.8) / init.rur.ind$s.num.8
ind.rur.9 <- sum(df.no.rur.ind$si.flow.9) / init.rur.ind$s.num.9

ind.rur.all <- (sum(df.no.rur.ind$si.flow.1) + sum(df.no.rur.ind$si.flow.2) +
                  sum(df.no.rur.ind$si.flow.3) + sum(df.no.rur.ind$si.flow.4) +
                  sum(df.no.rur.ind$si.flow.5) + sum(df.no.rur.ind$si.flow.6) + 
                  sum(df.no.rur.ind$si.flow.7) + sum(df.no.rur.ind$si.flow.8) +
                  sum(df.no.rur.ind$si.flow.9)) / (init.rur.ind$s.num.1 + init.rur.ind$s.num.2 + init.rur.ind$s.num.3 +
                                                     init.rur.ind$s.num.4 + init.rur.ind$s.num.5 + init.rur.ind$s.num.6 +
                                                     init.rur.ind$s.num.7 + init.rur.ind$s.num.8 + init.rur.ind$s.num.9)

ind.rur.lt5 <- (sum(df.no.rur.ind$si.flow.1) + sum(df.no.rur.ind$si.flow.2) +
                 sum(df.no.rur.ind$si.flow.3)) / (sum(df.no.rur.ind$si.flow.1) + sum(df.no.rur.ind$si.flow.2) +
                                                    sum(df.no.rur.ind$si.flow.3) + sum(df.no.rur.ind$si.flow.4) +
                                                    sum(df.no.rur.ind$si.flow.5) + sum(df.no.rur.ind$si.flow.6) + 
                                                    sum(df.no.rur.ind$si.flow.7) + sum(df.no.rur.ind$si.flow.8) +
                                                    sum(df.no.rur.ind$si.flow.9))

ind.rur.lt6m <- (sum(df.no.rur.ind$si.flow.1)) / (sum(df.no.rur.ind$si.flow.1) + sum(df.no.rur.ind$si.flow.2) +
                                                  sum(df.no.rur.ind$si.flow.3) + sum(df.no.rur.ind$si.flow.4) +
                                                  sum(df.no.rur.ind$si.flow.5) + sum(df.no.rur.ind$si.flow.6) + 
                                                  sum(df.no.rur.ind$si.flow.7) + sum(df.no.rur.ind$si.flow.8) +
                                                  sum(df.no.rur.ind$si.flow.9))


ind_rur_values <- data.frame(
  Category = factor(c("<6mo", "6-11mo", "1-4y", "5-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y"),
                    levels = c("<6mo", "6-11mo", "1-4y", "5-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")),
  Value = c(
    ind.rur.1,
    ind.rur.2,
    ind.rur.3,
    ind.rur.4,
    ind.rur.5,
    ind.rur.6,
    ind.rur.7,
    ind.rur.8,
    ind.rur.9
  ),
  Site = rep("Rural", 9)
)


## India urban   ---------------------------------------------------------------
ind.urb.1 <- sum(df.no.urb.ind$si.flow.1) / init.urb.ind$s.num.1
ind.urb.2 <- sum(df.no.urb.ind$si.flow.2) / init.urb.ind$s.num.2
ind.urb.3 <- sum(df.no.urb.ind$si.flow.3) / init.urb.ind$s.num.3
ind.urb.4 <- sum(df.no.urb.ind$si.flow.4) / init.urb.ind$s.num.4
ind.urb.5 <- sum(df.no.urb.ind$si.flow.5) / init.urb.ind$s.num.5
ind.urb.6 <- sum(df.no.urb.ind$si.flow.6) / init.urb.ind$s.num.6
ind.urb.7 <- sum(df.no.urb.ind$si.flow.7) / init.urb.ind$s.num.7
ind.urb.8 <- sum(df.no.urb.ind$si.flow.8) / init.urb.ind$s.num.8
ind.urb.9 <- sum(df.no.urb.ind$si.flow.9) / init.urb.ind$s.num.9

ind.urb.all <- (sum(df.no.urb.ind$si.flow.1) + sum(df.no.urb.ind$si.flow.2) +
                  sum(df.no.urb.ind$si.flow.3) + sum(df.no.urb.ind$si.flow.4) +
                  sum(df.no.urb.ind$si.flow.5) + sum(df.no.urb.ind$si.flow.6) + 
                  sum(df.no.urb.ind$si.flow.7) + sum(df.no.urb.ind$si.flow.8) +
                  sum(df.no.urb.ind$si.flow.9)) / (init.urb.ind$s.num.1 + init.urb.ind$s.num.2 + init.urb.ind$s.num.3 +
                                                     init.urb.ind$s.num.4 + init.urb.ind$s.num.5 + init.urb.ind$s.num.6 +
                                                     init.urb.ind$s.num.7 + init.urb.ind$s.num.8 + init.urb.ind$s.num.9)

ind.urb.lt5 <- (sum(df.no.urb.ind$si.flow.1) + sum(df.no.urb.ind$si.flow.2) +
                  sum(df.no.urb.ind$si.flow.3)) / (sum(df.no.urb.ind$si.flow.1) + sum(df.no.urb.ind$si.flow.2) +
                                                     sum(df.no.urb.ind$si.flow.3) + sum(df.no.urb.ind$si.flow.4) +
                                                     sum(df.no.urb.ind$si.flow.5) + sum(df.no.urb.ind$si.flow.6) + 
                                                     sum(df.no.urb.ind$si.flow.7) + sum(df.no.urb.ind$si.flow.8) +
                                                     sum(df.no.urb.ind$si.flow.9))

ind.urb.lt6m <- (sum(df.no.urb.ind$si.flow.1)) / (sum(df.no.urb.ind$si.flow.1) + sum(df.no.urb.ind$si.flow.2) +
                                                    sum(df.no.urb.ind$si.flow.3) + sum(df.no.urb.ind$si.flow.4) +
                                                    sum(df.no.urb.ind$si.flow.5) + sum(df.no.urb.ind$si.flow.6) + 
                                                    sum(df.no.urb.ind$si.flow.7) + sum(df.no.urb.ind$si.flow.8) +
                                                    sum(df.no.urb.ind$si.flow.9))




ind_urb_values <- data.frame(
  Category = factor(c("<6mo", "6-11mo", "1-4y", "5-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y"),
                    levels = c("<6mo", "6-11mo", "1-4y", "5-9y", "10-19y", "20-29y", "30-39y", "40-59y", "60+y")),
  Value = c(
    ind.urb.1,
    ind.urb.2,
    ind.urb.3,
    ind.urb.4,
    ind.urb.5,
    ind.urb.6,
    ind.urb.7,
    ind.urb.8,
    ind.urb.9
  ),
  Site=rep("Urban",9)
)

ind.values <- rbind(ind_rur_values, ind_urb_values)

ind.plot <- ggplot(ind.values, aes(Category, y = Value, fill = Site)) +
  geom_bar(stat = "identity", position = "dodge") +  # Position bars side by side
  labs(
    title = "India - Incidence",
    x = " ",
    y = " "
  ) +
  scale_fill_manual(values = c("Rural" = "aquamarine4", "Urban" = "steelblue3")) +  # Custom colors
  theme_minimal() +
  ylim(c(0, 1)) 

ind.plot



# Calculate R effective --------------------------------------------------------

## India -----------------------------------------------------------------------
### Rural  ---------------------------------------------------------------------
ind.rur.2 <- df.no.rur.ind %>%
  rename(dates= time, I = si.flow) %>%
  dplyr::select(dates, I)
all.dates <- as.data.frame(seq(1,730,1))
names(all.dates) <- "dates"
ind.rur.epiestim <- merge(x=ind.rur.2, y=all.dates, by="dates", all="TRUE")
ind.rur.epiestim <- ind.rur.epiestim %>% mutate(I = ifelse(is.na(I), 0, I)) 

estimates.rur <- estimate_R(ind.rur.epiestim$I, 
                        method="parametric_si",
                        config = make_config(list(
                          mean_si = 2.6, 
                          std_si = 1.5))
)
median(estimates.rur$R$`Mean(R)`)
max(estimates.rur$R$`Mean(R)`)
### Urban  ---------------------------------------------------------------------
ind.urb.2 <- df.no.urb.ind %>%
  rename(dates= time, I = si.flow) %>%
  dplyr::select(dates, I)
all.dates <- as.data.frame(seq(1,730,1))
names(all.dates) <- "dates"
ind.urb.epiestim <- merge(x=ind.urb.2, y=all.dates, by="dates", all="TRUE")
ind.urb.epiestim <- ind.urb.epiestim %>% mutate(I = ifelse(is.na(I), 0, I)) 

estimates.urb <- estimate_R(ind.urb.epiestim$I, 
                        method="parametric_si",
                        config = make_config(list(
                          mean_si = 2.6, 
                          std_si = 1.5))
)


median(estimates.urb$R$`Mean(R)`)
max(estimates.urb$R$`Mean(R)`)

estimates <- data.frame(Date = c(estimates.urb$R$t_start, estimates.urb$R$t_start),
                        Rt = c(estimates.rur$R$`Mean(R)`, estimates.urb$R$`Mean(R)`),
                        Site = c(rep("Rural", nrow(estimates.rur$R)), rep("Urban", nrow(estimates.rur$R))))
  
ind.rt.plot <- ggplot(estimates, aes(Date, y = Rt, color = Site)) +
  geom_line(lwd = 2) +  # Position bars side by side
  labs(
    title = "India - Time-Varying Reproduction Number",
    x = "Date",
    y = "Rt"
  ) +
  scale_color_manual(values = c("Rural" = "aquamarine4", "Urban" = "steelblue3")) +  # Custom colors
  theme_minimal() 

ind.rt.plot

# TODO: Compare this to the incidence plots you created before. What is the relationship between Rt and incidence? ####
# TODO: Compare this to the incidence plots you created before. How does Rural differ from Urban? ####

