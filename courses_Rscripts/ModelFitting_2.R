# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fitting models f1-f3 to a randomly-generated visitation rate dataset -------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set.seed(1)

n <- 40 # Sample size
P <- runif(n, 0, 100) # Plant densities

a <- 3.333 # True value of 'a'
h <- 0.025 # True value of 'h'

k <- NULL # Visitation counts
for(i in 1:n){
  # Generate a 'k' for each 'P' assuming the truth is Type 2 (model f2)
  k[i] <- rpois(1, a * P[i] / (1 + a * h * P[i]) )
}


plot(P, k,
     xlim = c(0, max(P)),
     ylim = c(0, max(k)),
     xlab = 'Plant density (P)',
     ylab = 'Visitation count (k)')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Type I (model f1) ---------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlL.pois.f1 <- function(par){
  a <- par['a']
  -sum(dpois(k, a * P, log = TRUE))
}

fit.f1 <- optim(par = list(a = 1),
                fn = nlL.pois.f1)
print(fit.f1)

plot.f1 <- function(x, fit){
  a <- fit$par['a']
  return(a * x)
}
curve(plot.f1(x, fit.f1), 
      add = TRUE,
      lwd = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Continue yourself ------------------------------------------------------------
# Now fit the Type II (f2) and Type III (f3) models for yourself and compare the
# 1) parameter MLEs and 2) likelihoods.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Type II (model f2) --------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlL.pois.f2 <- function(par){
  a <- par['a']
  h <- par['h']
  -sum(dpois(k, a * P / (1 + a * h * P), log = TRUE))
}

fit.f2 <- optim(par = list(a = fit.f1$par['a'] * 5, 
                           h = 1/50),
                fn = nlL.pois.f2)
print(fit.f2)

plot.f2 <- function(x, fit){
  a <- fit$par['a']
  h <- fit$par['h']
  return(a * x / (1 + a * h * x))
}
curve(plot.f2(x, fit.f2), 
      add = TRUE,
      lwd = 2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Type III (model f3) -------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nlL.pois.f3 <- function(par){
  a <- par['a']
  h <- par['h']
  theta <- par['theta']
  -sum(dpois(k, a * P^theta / (1 + a * h * P^theta), log = TRUE))
}

fit.f3 <- optim(par = list(a = fit.f2$par['a'], 
                           h = fit.f2$par['h'], 
                           theta = 1),
                fn = nlL.pois.f3)
print(fit.f3)

plot.f3 <- function(x, fit){
  a <- fit$par['a']
  h <- fit$par['h']
  theta <- fit$par['theta']
  return(a * x^theta / (1 + a * h * x^theta))
}
curve(plot.f3(x, fit.f3), 
      add = TRUE,
      lwd = 2,
      col = 'orange')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparisons -----------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameter estimates ---------
print(c(a = a, h = h, theta = 1))
fit.f1$par
fit.f2$par
fit.f3$par

# Negative log-Likelihoods -----
fit.f1$value
fit.f2$value
fit.f3$value

##############################################################################
