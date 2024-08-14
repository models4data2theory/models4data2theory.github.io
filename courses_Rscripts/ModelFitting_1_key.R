# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use maximum likelihood to estimate 'lambda' 
# of a Poisson given hypothetical count data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# True value of 'lambda'
lambda.true <- 30.333 

# Sample size 'n'
n <- 5

# Generate 'n' random visitation counts 'k'
k <- rpois(n, lambda.true )

# Visualize the counts
hist(k, breaks = 0:max(k))

abline(v = lambda.true,
       lwd = 2,
       col = 'red')

abline(v = mean(k),
       lwd = 2,
       col = 'blue')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To do:

# 1) Define negative log likelihood function
nlL.pois <- function(lambda){
  -sum(dpois(k, lambda, log = TRUE))
}

# 2) Provide initial guess for 'lambda'
init.par <- list(lambda = 2)

# 3) Fit the model
suppressWarnings(
  optim(init.par, nlL.pois)
)

# Note that you can suppress the optim() 1-D warning by using
# suppressWarnings( optim(.....) )



