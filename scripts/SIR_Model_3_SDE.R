## SIR Model 3: stochastic differential equations
##
## Pulliam et al. (2020) Compartmental models for infectious disease dynamics: a taxonomy
##
## Juliet R.C. Pulliam, January 2020
## Some Rights Reserved
## CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/)

library(diffeqr)

## sir_1() -- Function for numerical analysis of model 1
sir_3 <- function(y, parms, t){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    dSdt <- -R0*S*I/N
    dIdt <- R0*S*I/N - I
    # Note: Population size is constant, so don't need to specify dRdt
    return(c(dSdt, dIdt))
  })
}

gnoise <- function(y, parms, t){
  return(c(0.3*y[1]*y[2]/parms[2],0))
}

time <- 0
MAXTIME <- 60
MAXINF <- 100 # For plotting
N0 <- 1000

init <- c(S = N0 - 1,
          I = 1)
values <- c(R0 = 1.5,        # Basic reproductive ratio
            N = N0)            # population size (constant)

time.out <- seq(0, MAXTIME)

diffeqr::sde.solve()

ts <- data.frame(diffeqr::sde.solve(sir_3, gnoise, init, time.out, p = values, saveat=0.005))

plotly::plot_ly(x = sol$t, y = sol$u, type = 'scatter', mode = 'lines')