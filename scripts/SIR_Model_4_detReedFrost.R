## SIR Model 4: Deterministic Reed Frost
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## Juliet R.C. Pulliam, 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)


##-- Model 4: Deterministic Reed-Frost --##

sir_4 <- function(y, parms, mtime = MAXTIME){
  with(c(as.list(y), parms),{
    qq <- 1-R0/N  # Pairwise probability of avoiding potentially infectious contact
    
    Cases <- I
    Sus <- S
    
    for(Time in 1 : mtime){
      Cases <- c(Cases, Sus[Time] * (1-qq^Cases[Time]))
      Sus <- c(Sus, Sus[Time] - Cases[Time+1])
    }
    return(data.frame(time = 0:mtime, S = Sus, I = Cases, R = N-Sus-Cases))
  })
}
