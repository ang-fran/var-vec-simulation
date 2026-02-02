rm(list = ls())
# load required packages
library(MASS)
library(vars) # for mvnorm

set.seed(0)
T = 500

# Covariance matrix
Sigma = matrix(c(1, 0.9,
                  0.9, 1), nrow = 2)

# Simulation function
simulate_var = function(A_list, T, Sigma) {
  
  p = length(A_list)
  k = nrow(A_list[[1]])
  
  # sanity checks
  stopifnot(all(sapply(A_list, nrow) == k))
  stopifnot(all(sapply(A_list, ncol) == k))
  stopifnot(nrow(Sigma) == k)
  
  y = matrix(0, nrow = T, ncol = k)
  
  for (t in (p + 1):T) { # Loop to carry out out vectorized summation
    y[t, ] = Reduce(
      `+`,
      lapply(1:p, function(j) A_list[[j]] %*% y[t - j, ])
    ) + MASS::mvrnorm(1, mu = rep(0, k), Sigma = Sigma)
  }
  
  colnames(y) = paste0("Y", 1:k)
  return(y)
}


# ---- VAR(1) ----
## Creating AR matrices for 2-variable scenario
# Stationary 
A1_stat = matrix(c(0.5, 0.1,
                    0.2, 0.4), 2, byrow = TRUE)

# Nonstationary
A1_nonstat = matrix(c(1.05, 0.1,
                       0.2, 0.95), 2, byrow = TRUE)

y_stat = simulate_var(list(A1_stat), T, Sigma)
y_nonstat = simulate_var(list(A1_nonstat), T, Sigma)

