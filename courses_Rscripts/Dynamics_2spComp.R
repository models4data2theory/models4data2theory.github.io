# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Dynamics of 2-species Lotka-Volterra mutualism
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(deSolve)

model <- function(t, x, params) {
  N1 <- x[1]
  N2 <- x[2]
  with(as.list(params), {
    dN1dt <- r1 * N1 * (1 - N1 / K1 + a12 * N2 / K1)
    dN2dt <- r2 * N2 * (1 - N2 / K2 + a21 * N1 / K2)
  out <- c(dN1dt, dN2dt)
  list(out)
  })
}


# Time points at which to evaluate
t <- seq(0, 30, by = 1)

# Starting population sizes at t = 0
xstart <- c(N1 = 0.01, 
            N2 = 0.01) 

# Parameter values
parameters <- c(a12 = 0.5,
                a21 = 0.7,
                r1 = 1,
                r2 = 1,
                K1 = 1,
                K2 = 1)

# Integrate
out <- as.data.frame(ode(xstart, t, model, parameters))

# Plot dynamics over time
plot(out$time,
     out$N1,
     ylab = 'Popn size',
     xlab = 'Time',
     type = 'l',
     ylim = c(0, max(out$N1, out$N2) * 1.1)
)
lines(out$time, 
      out$N2, 
      lty = 2)
legend('bottomright',
       legend = c('Sp1', 'Sp2'),
       lty = c(1, 2))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~