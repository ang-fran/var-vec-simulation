rm(list = ls())

library(tseries)
library(vars)
library(urca)

set.seed(0)

T  = 2000
n  = 3
p  = 2            # levels VAR order
K  = 2            # ca.jo K = levels VAR order
ecdet = "none"

# ---- 1) Functions ----
var_roots_from_A = function(A_list) {
  k = nrow(A_list[[1]])
  p = length(A_list)
  
  top_block = do.call(cbind, A_list)  # k x (k*p)
  if (p == 1) {
    F = top_block
  } else {
    lower_block = cbind(diag(k * (p - 1)),
                        matrix(0, nrow = k * (p - 1), ncol = k))
    F = rbind(top_block, lower_block)
  }
  eigen(F)$values
}

count_unit_roots = function(roots, tol = 1e-3) sum(abs(Mod(roots) - 1) < tol)

extract_A_list_varest = function(v) Acoef(v)

adf_print = function(x, name) {
  cat("\nADF:", name, "\n")
  print(tseries::adf.test(as.numeric(x)))
}

eigvals_mat = function(M) as.vector(eigen(M)$values)

print_eigs = function(M, label) {
  cat("\nEigenvalues of", label, ":\n")
  print(eigvals_mat(M))
  invisible(NULL)
}

# VECM(p=2) => VAR(2) identities:
#   Δy_t = Π y_{t-1} + Γ Δy_{t-1} + ε_t
#   y_t  = A1 y_{t-1} + A2 y_{t-2} + ε_t
# where
#   A1 = I + Π + Γ
#   A2 = -Γ
# hence
#   Γ = -A2
#   Π = A1 + A2 - I

vecm_params_from_A = function(A1, A2) {
  I = diag(nrow(A1))
  Gamma_hat = -A2
  Pi_hat    = A1 + A2 - I
  list(Pi = Pi_hat, Gamma = Gamma_hat)
}

# 2) ---- DGP constructors (p=2) ----
#    r=0,1,2 are I(1) cases; r=3 is I(0) stationary

make_case_mats = function(r) {
  I = diag(n)
  
  # Default Gamma for the I(1) cases (small eigenvalues)
  Gamma_default = 0.20 * diag(n)
  
  if (r == 0) {
    alpha = matrix(0, n, 3)
    beta  = matrix(0, n, 3)
    Pi = alpha %*% t(beta)           # = 0
    Gamma = Gamma_default
    A1 = I + Pi + Gamma
    A2 = -Gamma
    return(list(r=r, Pi=Pi, alpha=alpha, beta=beta, Gamma=Gamma, A1=A1, A2=A2,
                mode="I(1) no cointegration (r=0)"))
  }
  
  if (r == 1) {
    beta  = matrix(c(1,0,0), nrow=n, ncol=1)           # z = y1
    alpha = matrix(c(-0.80,0,0), nrow=n, ncol=1)       # adjust y1 strongly
    Pi = alpha %*% t(beta)
    Gamma = Gamma_default
    A1 = I + Pi + Gamma
    A2 = -Gamma
    return(list(r=r, Pi=Pi, alpha=alpha, beta=beta, Gamma=Gamma, A1=A1, A2=A2,
                mode="I(1) cointegrated (r=1)"))
  }
  
  if (r == 2) {
    beta  = cbind(c(1,0,0), c(0,1,0))                  # z1=y1, z2=y2
    alpha = cbind(c(-0.70,0,0), c(0,-0.90,0))          # adjust y1,y2
    Pi = alpha %*% t(beta)
    Gamma = Gamma_default
    A1 = I + Pi + Gamma
    A2 = -Gamma
    return(list(r=r, Pi=Pi, alpha=alpha, beta=beta, Gamma=Gamma, A1=A1, A2=A2,
                mode="I(1) cointegrated (r=2)"))
  }
  
  if (r == 3) {
    # r=3: FULL-RANK Pi = alpha %*% t(beta) and FULL-RANK Gamma
    beta = matrix(c(
      1,  1,  0,
      0,  1,  1,
      1,  0,  1
    ), nrow = 3, byrow = TRUE)
    
    alpha = matrix(c(
      -0.70,  0.00,  0.00,
      0.00, -0.80,  0.00,
      0.00,  0.00, -0.60
    ), nrow = 3, byrow = TRUE)
    
    Pi = alpha %*% t(beta)  # full rank
    
    Gamma = matrix(c(
      0.20, 0.02, 0.00,
      0.01, 0.15, 0.01,
      0.00, 0.01, 0.10
    ), nrow = 3, byrow = TRUE)  # full rank, "small"
    
    A1 = I + Pi + Gamma
    A2 = -Gamma
    
    return(list(r=r, Pi=Pi, alpha=alpha, beta=beta, Gamma=Gamma, A1=A1, A2=A2,
                mode="Stationary VAR(2) in levels via full-rank Pi=alpha beta' (r=3; I(0))"))
  }
  
  stop("r must be 0..3")
}

# ---- 3) Simulation ----

sim_vecm_p2 = function(T, Pi, Gamma, y1=NULL, y2=NULL) {
  if (is.null(y1)) y1 = rnorm(n)
  if (is.null(y2)) y2 = y1 + rnorm(n)
  
  Y = matrix(0, T, n)
  colnames(Y) = c("X","Y","Z")
  Y[1,] = y1
  Y[2,] = y2
  
  e = matrix(rnorm(T*n), T, n)
  
  for (t in 3:T) {
    dy_lag = Y[t-1,] - Y[t-2,]
    dy_t   = as.numeric(Pi %*% Y[t-1,] + Gamma %*% dy_lag + e[t,])
    Y[t,]  = Y[t-1,] + dy_t
  }
  ts(Y)
}

# 4) ---- Run case: simulate -> fit VAR(2) on levels -> recover A1,A2,Pi,Gamma ----
#    plus Johansen summaries for r=1,2

run_case = function(r) {
  mats = make_case_mats(r)
  
  # Simulate from the VECM(p=2) DGP (works for all r here)
  Y = sim_vecm_p2(T, Pi = mats$Pi, Gamma = mats$Gamma)
  
  # Fit plain VAR(2) on levels (same order as DGP)
  v_lvl = VAR(Y, p = 2, type = "none")
  A_hat_list = extract_A_list_varest(v_lvl)
  A1_hat = A_hat_list[[1]]
  A2_hat = A_hat_list[[2]]
  
  # implied VECM params from estimated A1,A2
  implied_hat = vecm_params_from_A(A1_hat, A2_hat)
  Pi_hat = implied_hat$Pi
  Gamma_hat = implied_hat$Gamma
  
  # companion roots
  roots = var_roots_from_A(A_hat_list)
  
  # Johansen (only meaningful for I(1) cases r=1,2; r=0 has no coint; r=3 is I(0))
  jo = NULL
  if (r %in% c(1,2)) {
    jo = ca.jo(Y, type="trace", ecdet=ecdet, K=K, spec="transitory")
  }
  
  list(r=r, mats=mats, Y=Y,
       v_lvl=v_lvl,
       A_hat_list=A_hat_list,
       A1_hat=A1_hat, A2_hat=A2_hat,
       Pi_hat=Pi_hat, Gamma_hat=Gamma_hat,
       roots=roots, jo=jo)
}

# 5) ---- Reporting ----

report_case = function(res) {
  r = res$r
  mats = res$mats
  Y = res$Y
  roots = res$roots
  
  cat("\n====================================================\n")
  cat("CASE r =", r, " | ", mats$mode, "\n", sep="")
  cat("====================================================\n")
  
  # ---- True objects ----
  cat("\nTRUE A1:\n"); print(mats$A1)
  print_eigs(mats$A1, "TRUE A1")
  
  cat("\nTRUE A2:\n"); print(mats$A2)
  print_eigs(mats$A2, "TRUE A2")
  
  cat("\nTRUE Gamma (from DGP):\n"); print(mats$Gamma)
  cat("\nTRUE Pi (= A1 + A2 - I):\n"); print(mats$Pi)
  cat("rank(Pi) =", qr(mats$Pi)$rank, "\n")
  cat("rank(Gamma) =", qr(mats$Gamma)$rank, "\n")
  
  if (!is.null(mats$alpha) && !is.null(mats$beta)) {
    cat("\nTRUE alpha:\n"); print(mats$alpha)
    cat("\nTRUE beta:\n");  print(mats$beta)
  }
  
  # ---- Estimated objects ----
  cat("\nESTIMATED A1_hat (VAR(2) on levels):\n"); print(res$A1_hat)
  print_eigs(res$A1_hat, "ESTIMATED A1_hat")
  
  cat("\nESTIMATED A2_hat (VAR(2) on levels):\n"); print(res$A2_hat)
  print_eigs(res$A2_hat, "ESTIMATED A2_hat")
  
  cat("\nESTIMATED Gamma_hat (= -A2_hat):\n"); print(res$Gamma_hat)
  cat("\nESTIMATED Pi_hat (= A1_hat + A2_hat - I):\n"); print(res$Pi_hat)
  
  # ---- Errors ----
  cat("\nMax abs error in A1:", max(abs(res$A1_hat - mats$A1)), "\n")
  cat("Max abs error in A2:", max(abs(res$A2_hat - mats$A2)), "\n")
  cat("Max abs error in Gamma:", max(abs(res$Gamma_hat - mats$Gamma)), "\n")
  cat("Max abs error in Pi:", max(abs(res$Pi_hat - mats$Pi)), "\n")
  
  # ---- Roots / unit roots ----
  cat("\nCompanion roots of fitted LEVEL VAR(2):\n")
  print(roots)
  cat("Max modulus:", max(Mod(roots)), "\n")
  cat("# roots near 1 (tol=1e-3):", count_unit_roots(roots, tol=1e-3), "\n")
  if (r %in% c(0,1,2)) cat("Expected # unit roots = n - r =", (n - r), "\n")
  
  # ---- ADF on levels ----
  cat("\n--- ADF on LEVELS ---\n")
  adf_print(Y[,1], "X")
  adf_print(Y[,2], "Y")
  adf_print(Y[,3], "Z")
  
  # ---- ADF on TRUE cointegrating relations ----
  if (r %in% c(1,2)) {
    Ztrue = as.matrix(Y) %*% mats$beta
    cat("\n--- ADF on TRUE cointegrating relations z_t = beta' y_t ---\n")
    for (j in 1:ncol(Ztrue)) {
      adf_print(Ztrue[,j], paste0("z", j))
    }
    
    cat("\n--- Johansen summary (should indicate rank ~", r, ") ---\n", sep="")
    print(summary(res$jo))
  }
  
  invisible(NULL)
}

# 6) ---- Run all cases ----

results = lapply(0:3, run_case)
names(results) = paste0("r", 0:3)

invisible(lapply(results, report_case))
