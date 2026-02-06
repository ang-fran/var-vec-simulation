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

y_var1_stat = simulate_var(list(A1_stat), T, Sigma)
y_var1_nonstat = simulate_var(list(A1_nonstat), T, Sigma)


# ---- VAR(2)----
# Stationary
A2_1_stat = matrix(c(0.3, 0.1,
                      0.2, 0.4), nrow = 2, byrow = TRUE)

A2_2_stat = matrix(c(0.2, 0.05,
                      0.1, 0.1), nrow = 2, byrow = TRUE)

# Nonstationary
A2_1_nonstat = matrix(c(0.8, 0.3,
                        0.2, 0.6), nrow=2, byrow=TRUE)
A2_2_nonstat = matrix(c(0.3, 0.1,
                        0.1, 0.2), nrow=2, byrow=TRUE)

y_var2_stat = simulate_var(list(A2_1_stat, A2_2_stat), T, Sigma)
y_var2_nonstat = simulate_var(list(A2_1_nonstat, A2_2_nonstat), T, Sigma)


# ---- Stationarity check ----
check_var_stationarity = function(A_list) {
  # Inputs: A_list = list comprising of all AR coefficient matrices
  
  p = length(A_list) # lag order
  k = nrow(A_list[[1]]) # number of time series variables
  
  # Forming companion matrix
  F = do.call(cbind, A_list)
  
  # When p = 1, stationarity only depends on eigenvalues of A_1
  # otherwise, stationarity depends on eigenvalues of companion matrix
  # Companion matrix is made up of AR coefficient matrices in top row,
  # and identity matrices in rest of matrix
  
  if (p > 1) { 
    F = rbind(
      F,
      cbind(diag(k * (p - 1)), matrix(0, k * (p - 1), k))
    )
  }
  
  eig = eigen(F, only.values = TRUE)$values
  
  # Outputs: 
  # a. Eigenvalues of companion matrix
  # b. Stationary: True or False?
  
  list(
    eigenvalues = eig,
    stationary = all(Mod(eig) < 1)
  )
}

# ---- Visualization ----
png("figures/var1_stat.png", width = 800, height = 600)
ts.plot(y_var1_stat, main = "Simulated Stationary VAR(1)", ylab = "Y",
        col = c("violetred", "seagreen"), lty = c(1,3), lwd = c(0.5, 1))
legend(
  "topleft",
  legend = c("y_1t", "y_2t"),
  col = c("violetred", "seagreen"),
  lty = c(1, 3),
  lwd = c(0.5, 1),
  bty = "o"
)
dev.off()

png("figures/var1_nonstat.png", width = 800, height = 600)
ts.plot(y_var1_nonstat, main = "Simulated Nonstationary VAR(1)", ylab = "Y",
        col = c("red", "cyan"), lty = c(1,3), lwd = c(0.5, 1))
legend(
  "topleft",
  legend = c("y_1t", "y_2t"),
  col = c("red", "cyan"),
  lty = c(1, 3),
  lwd = c(0.5, 1),
  bty = "o"
)
dev.off()

png("figures/var2_stat.png", width = 800, height = 600)
ts.plot(y_var2_stat, main = "Simulated Stationary VAR(2)", ylab = "Y",
        col = c("royalblue", "hotpink3"), lty = c(1,3), lwd = c(0.5, 1))
legend(
  "topleft",
  legend = c("y_1t", "y_2t"),
  col = c("royalblue", "hotpink3"),
  lty = c(1, 3),
  lwd = c(0.5, 1),
  bty = "o"
)
dev.off()

png("figures/var2_nonstat.png", width = 800, height = 600)
ts.plot(y_var2_nonstat, main = "Simulated Nonstationary VAR(2)", ylab = "Y",
        col = c("chocolate4", "darkolivegreen"), lty = c(1,3), lwd = c(0.5, 1))
legend(
  "topleft",
  legend = c("y_1t", "y_2t"),
  col = c("chocolate4", "darkolivegreen"),
  lty = c(1, 3),
  lwd = c(0.5, 1),
  bty = "o"
)
dev.off()
