# Project:   compareHD_MVN
# Objective: costum MI pooling function for saturated model
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

miPool <- function(mi_fits, m, N){
  ## Description:
  # Given a list of outputs from a lavaan::cfa or sem model, fitted on multiple
  # datasets obtained by MI, it returns all parameter estiamtes, confidence
  # intervals, and FMI
  ## Example Inputs
  # mi_fits = miFitSat(mi_data = complete(mids_out, "all"),
  #                       model = satModWrite(names(dat_miss_lv[,
  #                                                       parms$vmap_miss$ta])))
  # m = parms$mice_ndt
  # N = parms$n

  ## Pool estiamtes
  # Create Indeces vectors for the parameters
  est <- parameterEstimates(mi_fits[[1]])
  par_names <- apply(est[, 1:3], 1, paste0, collapse = "")
  indx_avg <- est$op == "~1"
  indx_var <- est$lhs == est$rhs
  indx_cor <- est$op == "~~" & est$lhs != est$rhs

  #### Means ####
  est_avg <- sapply(mi_fits, function(x) {
    est <- parameterEstimates(x)
    return(est = est[indx_avg, "est"])
  })

  # Pool estimates
  Q_bar <- rowMeans(est_avg)

  # Pool Confidence Intervals
  all_vcov <- lapply(X = mi_fits, vcov)
  U_bar <- diag(Reduce('+', all_vcov) / m)[indx_avg]
  B <- diag(1 / (m-1) * (est_avg - Q_bar) %*% t(est_avg - Q_bar))
  T_var <- U_bar + B + B/m
  CI <- poolCI(m, N, Q_bar, B, T_var)

  # FMI and RIV
  fmi_out <- fmi(m = m, b = B, t = T_var)
  riv_out <- riv(m = m, b = B, u = U_bar)

  # Put together
  res_avg <- cbind(Q_bar = Q_bar, CI,
                   fmi = fmi_out, riv = riv_out)

  #### Variances ####
  est_var <- sapply(mi_fits, function(x) {
    est <- parameterEstimates(x)
    return(est[indx_var, "est"])
  })

  # Pool estimates
  Q_bar <- rowMeans(est_var)

  # Pool Confidence Intervals
  U_bar <- diag(Reduce('+', all_vcov) / m)[indx_var]
  B <- diag(1 / (m-1) * (est_var - Q_bar) %*% t(est_var - Q_bar))
  T_var <- U_bar + B + B/m
  CI <- poolCI(m, N, Q_bar, B, T_var)

  # FMI and RIV
  fmi_out <- fmi(m = m, b = B, t = T_var)
  riv_out <- riv(m = m, b = B, u = U_bar)

  # Put together
  res_var <- cbind(Q_bar = Q_bar, CI,
                   fmi = fmi_out, riv = riv_out)

  #### Correlations ####
  est_cor <- sapply(mi_fits, function(x) {
    est <- standardizedSolution(x)
    return(est[indx_cor, "est.std"])
  })

  # Transform before pooling
  est_cor_z <- fisher_z(est_cor)

  # Pool
  z_bar <- rowMeans(est_cor_z)
  z_U <- 1/(N-3)
  B <- diag(1 / (m-1) * (est_cor_z - z_bar) %*% t(est_cor_z - z_bar))
  T_var <- z_U + B + B/m

  # Degrees of freedom and FMI
  fmi_out <- fmi(m = m, b = B, t = T_var)
  riv_out <- riv(m = m, b = B, u = z_U)
  nu_com <- N - length(z_bar) # n - k where k number of paramteres estimated
  nu <- miDf(m, b = B, t = T_var, nu_com)

  # CI computation
  t_nu <- qt(1 - (1-.95)/2, df = nu)
  CI <- data.frame(lwr = z_bar - t_nu * sqrt(T_var),
                   upr = z_bar + t_nu * sqrt(T_var))
  res_z <- cbind(Q_bar = z_bar, CI)

  # Backtransform to correlations
  res_cor <- cbind(fisher_z_inv(res_z), fmi = fmi_out, riv = riv_out)

  #### Store ####
  pooled <- cbind(par = par_names, rbind(res_avg, res_var, res_cor))
  rownames(pooled) <- NULL

  return(pooled)
}