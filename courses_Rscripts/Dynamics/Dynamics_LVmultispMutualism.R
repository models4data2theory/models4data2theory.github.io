# Load matrix:
interactionMatrix = matrix(c(1, 1, 1, 1, 1, 0, 1, 0, 0), nrow = 3, byrow = TRUE)

# Initialize variables
time_steps <- 12
results <- matrix(0, nrow = time_steps, ncol = 2)

# Get the number of rows (plant species) and columns (animal species) of the matrix encoding the plant-pollinator network (interactionMatrix)
num_plants <- nrow(interactionMatrix)
num_animals <- ncol(interactionMatrix)

# Number of interactions
L=num_plants*num_animals


# Set the mean and variance for uniform random distributions
mean_alpha <- 0.5
var_alpha <- 0.2

mean_r <- 0.3
var_r <- 0.1

mean_s <- 0.2
var_s <- 0.1

mean_N <- 0.5
var_N <- 0.1

# Generate the 3x3 matrix with values from a uniform random distribution
matrix_alpha <- runif(L, mean_alpha - var_alpha/2, mean_alpha + var_alpha/2)
matrix_alpha <- matrix(matrix_alpha, nrow = num_plants, byrow = TRUE)

# Generate intrinsic growth rates and self-limitations vectors for plants and animals with values from their respective uniform random distributions
rP <- runif(num_plants, mean_r - var_r/2, mean_r + var_r/2)
rA <- runif(num_animals, mean_r - var_r/2, mean_r + var_r/2)
sP <- runif(num_plants, mean_s - var_s/2, mean_s + var_s/2)
sA <- runif(num_animals, mean_s - var_s/2, mean_s + var_s/2)

# Generate initial conditions vectors for plants and animals with values from their respective uniform random distributions
N_P <- runif(num_plants, mean_N - var_N/2, mean_N + var_N/2)
N_A <- runif(num_animals, mean_N - var_N/2, mean_N + var_N/2)

# Combine initial conditions into a single vector
initial_conditions <- c(N_P, N_A)

# Define the system of ODEs
plant_pollinator_ode <- function(time, state, parameters) {
  
  N_P <- state[1:num_plants] # Plant populations
  N_A <- state[(num_plants + 1):(num_plants + num_animals)] # Pollinator populations
  # Calculate the rate of change
  dPdt <- (rP - sP * N_P + (matrix_alpha * interactionMatrix) %*% N_A) * N_P
  dAdt <- (rA - sA * N_A + t(matrix_alpha * interactionMatrix) %*% N_P) * N_A
  # Return the rate of change
  list(c(dPdt, dAdt))
}

# Define the time sequence for which we want the solution
times <- seq(0, 12, by = 0.01)

# Solve the system using ode
results <- ode(y = initial_conditions, times = times, func = plant_pollinator_ode, parms = NULL, method="rk4")

# Print the results
print(as.data.frame(results))