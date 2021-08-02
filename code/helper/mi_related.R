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