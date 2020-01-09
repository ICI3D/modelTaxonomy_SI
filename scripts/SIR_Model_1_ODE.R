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

