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
##

SEED <- 20200108
set.seed(SEED)

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

N0 <- 1000
R_0 <- 1.5
MAXTIME <- 60
MAXINF <- 100

init <- c(S = N0 - 1,
          I = 1)
values <- c(R0 = R_0,        # Basic reproductive ratio
            N = N0)          # population size (constant)

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

ts <- data.frame(runGillespie(
  y = init,               # Initial conditions for population
  func = sir_2,                   # Function to evaluate
  parms = values                # Vector of parameters
))



par(lwd = 1, mar = c(5, 5.1, 1, 1), cex.lab = 3, xaxt = 'n', yaxt = 'n')
plot(ts$time,               # Time on the x axis
     ts$I,                  # Number infected (I) on the y axis
     xlab = "Time",         # Label the x axis
     ylab = NA,             # Label the y axis
     xlim = c(0, MAXTIME),  # Set limits of x axis    
     ylim = c(0, MAXINF),   # Set limits of y axis        
     type = "l",            # Use a line plot
     bty = "n",             # Remove the box around the plot
     col = 2)               # Set color of line
axis(1, at = seq(0, MAXTIME, MAXTIME / 3), labels = NA, xaxt = 's', lwd = 2)
axis(2, at = seq(0, MAXINF, MAXINF / 5), labels = NA, yaxt = 's', lwd = 2)
mtext(text = 'B2', side = 3, outer = F, adj = -.35, cex = 3, padj = 0.6)
mtext(text = 'Infected', side = 2, cex = 3, line = 2) # Label the y axis

for(ii in 2:15){
  ts <- data.frame(runGillespie(
    y = init,               # Initial conditions for population
    func = sir_2,                   # Function to evaluate
    parms = values                # Vector of parameters
  ))
  lines(ts$time,               # Time on the x axis
       ts$I,                  # Number infected (I) on the y axis
       col = ii+1)             # Remove the box around the plot
}
