# Project:   compareHD_MVN
# Objective: function to generate data according to a CFA model
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

dataGen <- function(parms, cond){

  ## Inputs ##
  # cond <- conds[1, ]

  ## Body ##
  # 1. Generate covariance matrix -----------------------------------------
    Sigma <- diag(cond$p)
  # Block 1: highly correlated variables
    Sigma[unlist(parms$vmap_miss), ] <- parms$rho_high
  # Block 2: not so highly correlated variables
    Sigma[-unlist(parms$vmap_miss), ] <- parms$rho_low
  # Fix diagonal
    diag(Sigma) <- 1
  # Make symmetric
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]

  # 2. Gen n x p data -----------------------------------------------------
  Z <- mvrnorm(n     = parms$n,
               mu    = rep(0, cond$p),
               Sigma = Sigma)
  colnames(Z) <- paste0("z", 1:ncol(Z))

  # 3. Scale according to your evs examination ----------------------------
  Z_sc <- Z * sqrt(parms$item_var) + parms$item_mean

  return( as.data.frame( Z_sc ) )
}
