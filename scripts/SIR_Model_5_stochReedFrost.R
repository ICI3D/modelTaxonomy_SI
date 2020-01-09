## SIR Model 5: Stochastic Reed Frost
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## This model example is modified from:
##
## Introduction to Infectious Disease Dynamics (Chain Binomial - Tutorial 3)
## Clinic on Meaningful Modeling of Epidemiological Data
## International Clinics on Infectious Disease Dynamics and Data Program
## African Institute for Mathematical Sciences, Muizenberg, Cape Town, RSA
##
## https://raw.githubusercontent.com/ICI3D/RTutorials/master/ICI3D_Example_chainBinom.R
##
## Juliet R.C. Pulliam
## CC-BY-NC, 2012
## 
## Last updated by Juliet R.C. Pulliam, January 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)

##-- Model 5: Stochastic Reed-Frost --##

sir_5 <- function(y, parms, mtime = MAXTIME){
  with(c(as.list(y), parms),{
    qq <- 1-R0/N  # Pairwise probability of avoiding potentially infectious contact
    
    Cases <- I
    Sus <- S
    
    for(Time in 1 : mtime){
      Cases <- c(Cases, rbinom(1, Sus[Time], (1-qq^Cases[Time])))
      Sus <- c(Sus, Sus[Time] - Cases[Time+1])
    }
    return(data.frame(time = 0:mtime, S = Sus, I = Cases, R = N-Sus-Cases))
  })
}
