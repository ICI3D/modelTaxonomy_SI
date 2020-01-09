## The SIR Model Family - master script
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
## 
## Juliet R.C. Pulliam, 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)
##
## Parts of this script are based on tutorials developed by the author for the following
## program:
##
## Clinic on Meaningful Modeling of Epidemiological Data
## International Clinics on Infectious Disease Dynamics and Data Program
## African Institute for Mathematical Sciences, Muizenberg, Cape Town, RSA
##
## https://raw.githubusercontent.com/ICI3D/RTutorials/master/ICI3D_Lab1_ODEmodels.R
##
## Juliet R.C. Pulliam
## CC-BY-NC, 2012
## 
## This file was last updated by Juliet R.C. Pulliam, January 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)

##--- The SIR Model Family ---##
##
## The SIR model world has 3 state variables (though only two of these need
## to be specified to define the state of the system, since the population
## size is constant):
##
## S - the number of susceptibles in the population
## I - the number of infected/infectious individuals in the population
## R - the number of recovered/immune individuals in the population
##
## and 2 parameters:
##
## N - the total population size; N = S + I + R
## R_0 - the basic reproductive number
##
## The latter parameter is often written in terms of two distinct parameters:
## R_0 = beta / gamma, where these parameters are defined as follows:
##
## beta - the transmission coefficient; beta = the contact rate * infectivity
## gamma - 1 / the infectious period
##
## However, the behavior of the system of equations is dependent only on R_0.
##

library(deSolve)               # Load libary to be used for numerical integration of ODE system
library(msde)                  # Load libary to be used for approximation of SDE system
library(RColorBrewer)          # Load library for colors

##--- Source scripts with model definitions and plotting functions ---##

source('SIR_Model_1_ODE.R')
source('SIR_Model_2_Gillespie.R')
# source('SIR_Model_3_SDE.R')
source('SIR_Model_4_detReedFrost.R')
source('SIR_Model_5_stochReedFrost.R')
source('SIR_Model_6_chainBinomial.R')
source('SIR_plot.R')

##--- Define parameters, control variables, and initial conditions ---##

SEED <- 20200109
set.seed(SEED)

MAXTIME <- 60
MAXINF <- 100 # For plotting
NREPS <- 15   # For stochastic models

time <- 0
N0 <- 1000

init <- c(S = N0 - 1,
          I = 1)
values <- c(R0 = 1.5,        # Basic reproduction number
            N = N0)          # population size (constant)
values.cb <- c(bb = 1.5,     # transmission coefficient
               N = N0,       # population size (constant)
               gg = 1)       # 1 / duration of infection (removal rate)

##--- Model 1: Ordinary differential equation model ---##

ts_ODE <- data.frame(lsoda(
  y = init,                  # Initial conditions for population
  times = seq(0, MAXTIME, 0.01), # Timepoints for evaluation
  func = sir_1,              # Function to evaluate
  parms = values             # Vector of parameters
))

plotSIR(ts_ODE, plotLab = 'B1')

##--- Model 2: The Gillespie algorithm ---##

ts_Ga <- data.frame(runGillespie(
  y = init,                  # Initial conditions for population
  func = sir_2,              # Function to evaluate
  parms = values             # Vector of parameters
))

plotSIR(ts_Ga, plotLab = 'B2', lwd = 1, col = brewer.pal(9, 'Set1')[1])

for(ii in 2:NREPS){
  ts_Ga <- data.frame(runGillespie(
    y = init,                # Initial conditions for population
    func = sir_2,            # Function to evaluate
    parms = values           # Vector of parameters
  ))
  lines(ts_Ga$time,          # Time on the x axis
        ts_Ga$I,             # Number infected (I) on the y axis
        col = brewer.pal(9, 'Set1')[ii%%9])
}

##--- Model 3: Stochastic differential equations ---##

### To be added

##--- Model 4: Deterministic Reed-Frost ---##

ts_dRF <- sir_4(
  y = init,                  # Initial conditions for population
  parms = values             # Vector of parameters
)

plotSIR(ts_dRF, plotLab = 'B4', lwd = 3, type = 's')

##--- Model 5: Stochastic Reed-Frost ---##

ts_sRF <- sir_5(
  y = init,                  # Initial conditions for population
  parms = values             # Vector of parameters
)

plotSIR(ts_sRF, plotLab = 'B5', lwd = 2, col = brewer.pal(9, 'Set1')[1], type = 's')

for(ii in 2:floor(NREPS/2)){
  ts_sRF <- sir_5(
    y = init,                # Initial conditions for population
    parms = values           # Vector of parameters
  )
  lines(ts_sRF$time,          # Time on the x axis
        ts_sRF$I,             # Number infected (I) on the y axis
        col = brewer.pal(9, 'Set1')[ii%%9], type = 's')
}

##--- Model 6: Generic chain binomial ---##

ts_cb <- sir_6(
  y = init,                  # Initial conditions for population
  tstep = 0.5,               # Timestep
  parms = values.cb          # Vector of parameters
)

plotSIR(ts_cb, plotLab = 'B6', lwd = 1, col = brewer.pal(9, 'Set1')[1], type = 's')

for(ii in 2:NREPS){
  ts_cb <- sir_6(
    y = init,                # Initial conditions for population
    tstep = 0.5,             # Timestep
    parms = values.cb           # Vector of parameters
  )
  lines(ts_cb$time,          # Time on the x axis
        ts_cb$I,             # Number infected (I) on the y axis
        col = brewer.pal(9, 'Set1')[ii%%9], type = 's')
}
