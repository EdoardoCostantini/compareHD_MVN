# Project:   compareHD_MVN
# Objective: subroutine runCell a single condition for a single repetition
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-08-03

runCell <- function(rp, cond, parms, fs) {

# Example Internals -------------------------------------------------------

  # rp   = 1
  # cond = conds[6, ]

# Data Generation ---------------------------------------------------------

  ## Data according to condition
  data_OG <- dataGen(parms = parms, cond = cond)

  ## Impose Missingness according to condition
  data_miss <- amputePerVar(data = data_OG,
                            targets = parms$vmap_miss$ta,
                            preds = parms$vmap_miss$mp,
                            pm = cond$pm,
                            type = "high")

# Imputation --------------------------------------------------------------

  # Are we doing imputation?
  if(!cond$method %in% c("OG", "cca")){
    # Are we using the optimal model?
    if(cond$method == "norm.optimal"){
      predmat <- make.predictorMatrix(data_miss)
      predmat[, -unlist(parms$vmap_miss)] <- 0 # use only good predictors
      mids_out <- mice(data_miss,
                       method = "norm",
                       predictorMatrix = predmat,
                       m = parms$mice_ndt,
                       maxit = parms$mice_iters,
                       printFlag = TRUE)
    }
    # Are we using blasso?
    if(cond$method == "blasso"){
      mids_out <- imputeBlasso(data_miss,
                               m = parms$blasso_ndt,
                               maxit = parms$blasso_iters)
    }
    # Are we using any of the other HD-MI methods?
    if(cond$method %in%
      c("durr.gaus",
        "iurr.gaus",
        "norm", # bridge
        "pcr.boot",
        "cart",
        "rf")){
      mids_out <- mice(data_miss,
                       method = as.vector(cond$method),
                       m = parms$mice_ndt,
                       maxit = parms$mice_iters,
                       printFlag = TRUE,
                       # Bypass linear dependency check
                       eps = 0,
                       # Argument used only by norm method
                       ls.meth = "ridge",
                       ridge   = 1e-3)
    }
  }

# Analyze and pool --------------------------------------------------------

  # Are we doing imputation?
  if(!cond$method %in% c("OG", "cca")){
    mi_sat_fits <- miFitSat(mi_data = complete(mids_out, "all"),
                            model = satModWrite(names(data_miss[,
                                                        parms$vmap_miss$ta]))
    )
    mi_sat_pool <- miPool(mi_fits = mi_sat_fits,
                          m = parms$mice_ndt,
                          N = parms$n)
    result <- mi_sat_pool
  } else {
    # Are we doing OG or cca?
    if(cond$method %in% "OG"){
      data_noMI <- data_OG
    } else {
      data_noMI <- na.omit(data_miss)
    }
      # Fit on single complete data
      og_sat_fit <- miFitSat(mi_data = list(data_noMI),
                             model = satModWrite(names(data_noMI[,
                                                         parms$vmap_miss$ta]))
      )[[1]]
      og_est_all <- parameterEstimates(og_sat_fit, standardized = TRUE)
      og_est <- og_est_all[, c("std.all", "est", "ci.lower", "ci.upper")]
      result <- cbind(par = apply(og_est_all[, 1:3], 1, paste0, collapse = ""),
                      og_est,
                      fmi = NA,
                      riv = NA)
  }

  # Attach descriptor
  row.names(cond) <- NULL # to avoid a warning in the cbind
  result <- cbind(rp = rp, cond, result)

# Store Output ------------------------------------------------------------

  ## Store simulation results
  if(parms$goal == "simulation"){
    saveRDS(result,
            file = paste0(fs$outDir,
                          "rep_", rp, "_", cond$tag,
                          ".rds")
    )
  }
  ## or store convergence check results
  if(parms$goal == "conv_check"){
    saveRDS(mids_out,
            file = paste0(fs$outDir,
                          "rep_", rp, "_", cond$tag,
                          ".rds")
    )
  }
}