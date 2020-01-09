## SIR Model 6: Generic chain binomial
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## Juliet R.C. Pulliam, 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)

##-- Model 6: Generic chain binomial --##

# Note that this impelementation, unlike the others, requires the separation of time scales.

sir_6 <- function(y, tstep, parms, mtime = MAXTIME){
  with(c(as.list(y), parms),{
    qq <- exp(-bb/N*tstep)  # Pairwise probability of avoiding potentially infectious contact
    
    Cases <- I
    Sus <- S
    
    for(Time in seq(0, mtime, tstep)){
      newCases <- rbinom(1, Sus[length(Sus)], (1-qq^Cases[length(Cases)]))
      removedCases <- rbinom(1, Cases[length(Cases)], 1 - exp(-gg*tstep))
      Sus <- c(Sus, Sus[length(Sus)] - newCases)
      Cases <- c(Cases, Cases[length(Cases)] + newCases - removedCases)
    }
    return(data.frame(time = seq(0, mtime + tstep, tstep), S = Sus, I = Cases, R = N-Sus-Cases))
  })
}
