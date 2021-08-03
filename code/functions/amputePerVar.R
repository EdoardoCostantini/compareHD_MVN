# Project:   compareHD_MVN
# Objective: function to ampute a collection of target variables
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

amputePerVar <- function(data, targets, preds, pm = .5, type = "high"){
  ## Description
  # Uses the simMissingness() function to impose missing values on
  # targets set of variables based on preds.
  # It returns the variables in targets with missing values imposed.
  ## Example Inputs
  # cond = conds[1, ]
  # data = dataGen(parms = parms, cond = cond)
  # targets = parms$vmap_miss$ta
  # preds = parms$vmap_miss$mp
  # pm = .5
  # type = "high"

  ## Body
  preds_object <- model.matrix(
    ~ ., as.data.frame(data[, preds])
  )[, -1, drop = FALSE]

  for (i in targets) {
    nR <- simMissingness(pm    = pm,
                         data  = preds_object,
                         type  = type)

    # Fill in NAs
    data[nR, i] <- NA
  }

  # Result
  return( data )
}