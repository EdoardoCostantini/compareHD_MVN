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
  ests <- lapply(mi_fits, function(x) {
    est_all <- parameterEstimates(x, standardized = TRUE)
    est <- est_all[, c("est", "std.all")]
    rownames(est) <- apply(est_all[, 1:3], 1, paste0, collapse = "")
    return(est)
  })

  ests_both <- do.call(cbind, ests)
  coefs_raw <- as.matrix(ests_both[, grep("est", colnames(ests_both))])
  coefs_std <- ests_both[, grep("std.all", colnames(ests_both))]
  Q_bar <- rowMeans(coefs_raw)
  Q_bar_std <- rowMeans(coefs_std)

  ## Pool Confidence Intervals
  all_vcov <- lapply(X = mi_fits, vcov)
  U_bar <- diag(Reduce('+', all_vcov) / m)
  B <- diag(1 / (m-1) * (coefs_raw - Q_bar) %*% t(coefs_raw - Q_bar))
  T_var <- U_bar + B + B/m

  # Degrees of freedom and FMI
  fmi_out <- fmi(m = length(mi_fits), b = B, t = T_var)
  riv_out <- riv(m = length(mi_fits), b = B, u = U_bar)
  nu_com <- N - length(Q_bar) # n - k where k number of paramteres estimated
  nu <- miDf(length(mi_fits), b = B, t = T_var, nu_com)

  # CI computation
  t_nu <- qt(1 - (1-.95)/2, df = nu)
  CI <- data.frame(lwr = Q_bar - t_nu * sqrt(T_var),
                   upr = Q_bar + t_nu * sqrt(T_var))

  ## Store
  pooled <- cbind(names(Q_bar), Q_bar_std, Q_bar, CI, fmi_out, riv_out)
  colnames(pooled) <- c("par", "std.all", "est", "ci.lower", "ci.upper",
                        "fmi", "riv")

  return(pooled)
}