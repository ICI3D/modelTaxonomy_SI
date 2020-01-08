## SIR Model 1: ordinary differential equations
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## This model example is modified from:
##
## Lab 1: ODE Models in R
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
##
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

library(deSolve)                # Load libary to be used for numerical integration

##-- Model 1: ODE implementation --##
##
## sir_1() -- Function for numerical analysis of model 1
sir_1 <- function(t, y, parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    dSdt <- -R0*S*I/N
    dIdt <- R0*S*I/N - I
    # Note: Population size is constant, so don't need to specify dRdt
    list(c(dSdt, dIdt))
  })
}


time <- 0
MAXTIME <- 60
MAXINF <- 100 # For plotting
N0 <- 1000

## Define a vector init that gives the initial values of our 2 state variables:

init <- c(S = N0 - 1,
            I = 1)
values <- c(R0 = 1.5,        # Basic reproductive ratio
            N = N0)            # population size (constant)

time.out <- seq(0, MAXTIME, 0.01)

ts <- data.frame(lsoda(
  y = init,               # Initial conditions for population
  times = time.out,             # Timepoints for evaluation
  func = sir_1,                   # Function to evaluate
  parms = values                # Vector of parameters
))

par(lwd = 3, mar = c(5, 5.1, 1, 1), cex.lab = 3, xaxt = 'n', yaxt = 'n')
plot(ts$time,               # Time on the x axis
     ts$I,                  # Number infected (I) on the y axis
     xlab = "Time",         # Label the x axis
     ylab = NA,             # Label the y axis
     xlim = c(0, MAXTIME),  # Set limits of x axis    
     ylim = c(0, MAXINF),   # Set limits of y axis        
     type = "l",            # Use a line plot
     bty = "n")             # Remove the box around the plot
axis(1, at = seq(0, MAXTIME, MAXTIME / 3), labels = NA, xaxt = 's', lwd = 2)
axis(2, at = seq(0, MAXINF, MAXINF / 5), labels = NA, yaxt = 's', lwd = 2)
mtext(text = 'B1', side = 3, outer = F, adj = -.35, cex = 3, padj = 0.6)
mtext(text = 'Infected', side = 2, cex = 3, line = 2) # Label the y axis
