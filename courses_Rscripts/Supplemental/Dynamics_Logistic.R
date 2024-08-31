# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The continuous-time logistic
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# install.packages('deSolve')
library(deSolve)

# Define the function for the continuous logistic.
# Notice that this function has three inputs:  
#    - the starting population size'N0',
#    - the total time 'T', 
#    - variable 'p' representing the two parameters of our equation: r and K.

ClogisK <- function(t, y, p) {
  N <- y[1]
  with(as.list(p), {
    dNdt <- r * N * (1 - N / K)
    return(list(dNdt))
  })
}

# Include both r and K in the params variable:
params <- c(r = 0.1, 
            K = 100)

# Specify the starting abundance and the time-points at which to determine the population sizes:
N0 <- c(N = 1)
t <- seq(1, 100)

# Integrate over time using the ode() function of deSolve
out <- ode(
  y = N0,
  times = t,
  func = ClogisK,
  parms = params
)

# Convert this output to a data.frame to ease plotting and data extraction
out <- data.frame(out)

# Determine the change in N between times t and t+1
N <- out$N
deltaN <- N[-1] - N[-length(N)]

# Plot the dynamics and growth rates side-by-side
par(mfrow = c(1, 2))
plot(out$time, out$N, 
     type = 'l', 
     lwd = 4)
plot(N[-length(N)], deltaN,
     type = 'l',
     lwd = 4,
     xlab = 'N')

################
# Try out various combinations of N0, T, r, and K.
# How do the values of r, K, and N0 affect the maximum value of deltaN (i.e. dNdt)?

#########################################################################################
#########################################################################################
#########################################################################################
