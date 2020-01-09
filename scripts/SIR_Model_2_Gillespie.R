## SIR Model 2: the Gillespie approxmation
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## This model example is modified from:
##
## Gillespie algorithm benchmark question: key for part (c)
## Clinic on Meaningful Modeling of Epidemiological Data
## International Clinics on Infectious Disease Dynamics and Data Program
## African Institute for Mathematical Sciences, Muizenberg, Cape Town, RSA
##
## https://raw.githubusercontent.com/ICI3D/RTutorials/master/gillespie/ICI3D_Example_StochSIRSgillespie_PartC.R
##
## Juliet R.C. Pulliam
## CC-BY-NC, 2012
## 
## Last updated by Juliet R.C. Pulliam, January 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)

##-- Model 2: Gillespie implementation --##

sir_2 <- function(time, svars, parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(svars), parms),{
    trans.rate <- R0*S*I/N
    remov.rate <- I
    
    total.rate <- trans.rate + remov.rate # Total rate
    event.time <- time + rexp(1, total.rate) # Time of next event
    
    dd <- runif(1)
    if(dd < trans.rate/total.rate){
      # Next event is transmission
      S <- S-1
      I <- I+1
      event <- 'transmission'
    }else{
      # Next event is removal
      I <- I-1
      R <- R+1  	  
      event <- 'removal'
    }
    
    return(data.frame(time = event.time, event = event, S = S, I = I, R = R))
  })
}

runGillespie <- function(y, func, parms, tmax = MAXTIME){
  with(c(as.list(y), parms),{
    next.time <- data.frame(time = 0, event = NA, S = S, I = I, R = N - S - I)
    out <- next.time
    while(next.time$time < tmax & next.time$I > 0){
      pop <- c(S = next.time$S, I = next.time$I, R = next.time$R)
      next.time <- sir_2(next.time$time, pop, parms)
      out <- rbind(out, next.time)
    }
    return(out)
  })
}  
