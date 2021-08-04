# Project:   compareHD_MVN
# Objective: Collection of functions needed for MI tasks
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

fmi <- function(m, b, t){
  # proportion of variation attributable to the missing data
  # aka fmi
  # (van Buuren, 2018, p. 46)
  fmi <- (1 + 1/m) * b/t
  return(fmi)
}

riv <- function(m, b, u){
  # relative increase in variance due to nonresponse
  # (van Buuren, 2018, p. 47)
  riv <- (1 + 1/m) * b/u
  return(riv)
}

miDf <- function(m, b, t, dfCom) {
  fmi   <- fmi(m, b, t)
  df0   <- (m - 1) * (1 / fmi^2)
  dfObs <- (dfCom + 1) / (dfCom + 3) * dfCom * (1 - fmi)

  df0 / (1 + (df0 / dfObs))
}

fisher_z <- function(x){
  # Transformation of correlation to fisher Z pre-pooling
  .5 * log((1+x) / (1-x))
}

fisher_z_inv <- function(Z_bar){
  # Back transformation of fisher Z after pooling
  (exp(2 * Z_bar) - 1) / (exp(2 * Z_bar) + 1)
}

poolCI <- function (m, N, Q_bar, B, T_var){
  # Pool confidence intervals for a normally distributed estimate
  nu_com <- N - length(Q_bar) # n - k where k number of paramteres estimated
  nu <- miDf(m, b = B, t = T_var, nu_com) # Degrees of freedom
  t_nu <- qt(1 - (1-.95)/2, df = nu) # critical value
  CI <- data.frame(lwr = Q_bar - t_nu * sqrt(T_var),
                   upr = Q_bar + t_nu * sqrt(T_var))
  return(CI)
}
