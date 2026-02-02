rm(list = ls())
# load required packages
library(MASS)
library(vars) # for mvnorm

set.seed(0)
T = 500

# Covariance matrix
Sigma = matrix(c(1, 0.9,
                  0.9, 1), nrow = 2)

# ---- VAR(1) ----
## Creating AR matrices for 2-variable scenario
# Stationary 
A1_stat = matrix(c(0.5, 0.1,
                    0.2, 0.4), 2, byrow = TRUE)

# Nonstationary
A1_nonstat = matrix(c(1.05, 0.1,
                       0.2, 0.95), 2, byrow = TRUE)


simulate_var1 = function(A, T, Sigma) {
  y = matrix(0, T, nrow(A))
  for (t in 2:T) {
    y[t, ] = A %*% y[t - 1, ] +
      mvrnorm(1, rep(0, nrow(A)), Sigma)
  }
  y
}

y_stat = simulate_var1(A1_stat, T, Sigma)
y_nonstat = simulate_var1(A1_nonstat, T, Sigma)
